# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

"""
    plate_eqs!(s, eqs, psys, wing; twist_surfaces, twist_angle,
               R_b_to_w, aero_force_b, aero_moment_b,
               aero_force_point_b, pos, vel, com_w, wind_vec_gnd, height)

Generate symbolic flat-plate aerodynamic equations for a `PlateWing` from its
flat-plate `TwistSurface` sections.
"""
function plate_eqs!(s, eqs, psys, wing; twist_surfaces, twist_angle,
                    R_b_to_w, aero_force_b, aero_moment_b,
                    aero_force_point_b, pos, vel, com_w,
                    wind_vec_gnd, height)
    section_idxs = wing.twist_surface_idxs
    Rbw = collect(R_b_to_w[:, :, wing.idx])

    n_surf = length(section_idxs)
    @variables begin
        plate_x_w(t)[1:3, 1:n_surf]
        plate_y_w(t)[1:3, 1:n_surf]
        plate_z_w(t)[1:3, 1:n_surf]
        plate_va_w(t)[1:3, 1:n_surf]
        plate_v_tan(t)[1:n_surf]
        plate_v_norm(t)[1:n_surf]
        plate_alpha(t)[1:n_surf]
        plate_cl(t)[1:n_surf]
        plate_cd(t)[1:n_surf]
        plate_q(t)[1:n_surf]
        plate_q_drag(t)[1:n_surf]
        plate_lift(t)[1:3, 1:n_surf]
        plate_drag(t)[1:3, 1:n_surf]
        plate_force_w(t)[1:3, 1:n_surf]
    end

    for (section_idx, ts_idx) in enumerate(section_idxs)
        ts = twist_surfaces[ts_idx]
        pidx = ts.point_idxs[1]

        x_b = smooth_normalize(collect(get_twist_surface_chord(psys, ts_idx)))
        y_b = collect(get_twist_surface_y_airf(psys, ts_idx))

        cos_twist = cos(twist_angle[ts_idx])
        sin_twist = sin(twist_angle[ts_idx])
        x_twisted = cos_twist * x_b + sin_twist * (y_b × x_b)
        z_twisted = x_twisted × y_b

        eqs = [
            eqs
            plate_x_w[:, section_idx] ~ Rbw * x_twisted
            plate_y_w[:, section_idx] ~ Rbw * y_b
            plate_z_w[:, section_idx] ~ Rbw * z_twisted
        ]

        wind_at_h = calc_wind_factor(
            s.am, pos[1, pidx], pos[2, pidx],
            pos[3, pidx], psys) * wind_vec_gnd
        eqs = [
            eqs
            plate_va_w[:, section_idx] ~ wind_at_h - vel[:, pidx]
        ]

        apparent_wind = collect(plate_va_w[:, section_idx])
        x_axis_w = collect(plate_x_w[:, section_idx])
        z_axis_w = collect(plate_z_w[:, section_idx])
        eqs = [
            eqs
            plate_v_tan[section_idx] ~ apparent_wind ⋅ x_axis_w
            plate_v_norm[section_idx] ~ apparent_wind ⋅ z_axis_w
            plate_alpha[section_idx] ~
                rad2deg(atan(plate_v_norm[section_idx],
                             plate_v_tan[section_idx]))
        ]

        eqs = [
            eqs
            plate_cl[section_idx] ~ get_plate_cl(
                psys, wing.idx, plate_alpha[section_idx])
            plate_cd[section_idx] ~
                get_plate_drag_corr(psys, wing.idx) *
                get_plate_cd(psys, wing.idx, plate_alpha[section_idx])
        ]

        eqs = [
            eqs
            plate_q[section_idx] ~ 0.5 *
                calc_rho(s.am, height[pidx]) *
                (plate_v_tan[section_idx]^2 + plate_v_norm[section_idx]^2)
            plate_q_drag[section_idx] ~ 0.5 *
                calc_rho(s.am, height[pidx]) *
                (apparent_wind ⋅ apparent_wind)
        ]

        alpha_rad = atan(plate_v_norm[section_idx], plate_v_tan[section_idx])
        va_airf_dir = cos(alpha_rad) * x_axis_w + sin(alpha_rad) * z_axis_w
        y_axis_w = collect(plate_y_w[:, section_idx])
        lift_dir = smooth_normalize(va_airf_dir × y_axis_w)
        drag_dir = smooth_normalize(y_axis_w × lift_dir)

        area = get_twist_surface_area(psys, ts_idx)
        eqs = [
            eqs
            plate_lift[:, section_idx] ~
                plate_q[section_idx] * area *
                plate_cl[section_idx] * lift_dir
            plate_drag[:, section_idx] ~
                plate_q_drag[section_idx] * area *
                plate_cd[section_idx] * drag_dir
            plate_force_w[:, section_idx] ~
                plate_lift[:, section_idx] + plate_drag[:, section_idx]
        ]
    end

    if wing.dynamics_type == PARTICLE_DYNAMICS
        aero_force_point = aero_force_point_b::AbstractArray
        for (section_idx, ts_idx) in enumerate(section_idxs)
            pidx = twist_surfaces[ts_idx].point_idxs[1]
            eqs = [
                eqs
                aero_force_point[:, pidx] ~
                    Rbw' * plate_force_w[:, section_idx]
            ]
        end
        eqs = [
            eqs
            aero_force_b[:, wing.idx] ~
                sum([Rbw' * plate_force_w[:, section_idx]
                     for section_idx in 1:n_surf])
            aero_moment_b[:, wing.idx] ~ zeros(3)
        ]
    elseif wing.dynamics_type == RIGID_DYNAMICS
        force_sum = sum([
            Rbw' * plate_force_w[:, section_idx]
            for section_idx in 1:n_surf])
        moment_sum = sum([
            Rbw' * ((pos[:, twist_surfaces[section_idxs[section_idx]].point_idxs[1]] -
                     com_w[:, wing.idx]) ×
                    plate_force_w[:, section_idx])
            for section_idx in 1:n_surf])
        eqs = [
            eqs
            aero_force_b[:, wing.idx] ~ force_sum
            aero_moment_b[:, wing.idx] ~ moment_sum
        ]
    end

    return eqs
end
