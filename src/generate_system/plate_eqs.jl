# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Symbolic equation generation for PlateWing flat-plate aerodynamics.
#
# For each PlateSurface: compute twist, rotate axes, calculate AoA
# from atan of wind components, look up CL/CD, compute lift/drag.
# All equations are symbolic — evaluated every ODE timestep.

"""
    plate_eqs!(eqs, wing, psys, R_b_to_w, va_point_b,
               aero_force_b, aero_moment_b,
               aero_force_point_b, pos, com_w, height)

Generate symbolic flat-plate aerodynamic equations for a PlateWing.

For each PlateSurface:
1. Read twist angle directly from surface
2. Rotate x_airf/z_airf around y_airf by twist
3. AoA = atan(v_normal, v_tangential)
4. CL/CD from registered interpolation lookup
5. Lift perpendicular to flow in airfoil plane, drag along flow

Forces are applied per-point (PARTICLE_DYNAMICS) or summed to wing
(RIGID_DYNAMICS).
"""
function plate_eqs!(s, eqs, psys, wing;
                    R_b_to_w, aero_force_b, aero_moment_b,
                    aero_force_point_b, pos, vel, com_w,
                    wind_vec_gnd, height)
    surfaces = wing.surfaces
    Rbw = collect(R_b_to_w[:, :, wing.idx])

    # Declare per-surface symbolic variables
    n_surf = length(surfaces)
    @variables begin
        plate_twist(t)[1:n_surf]
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

    for (section_idx, surf) in enumerate(surfaces)
        pidx = surf.point_idx

        # Twist from surface
        eqs = [
            eqs
            plate_twist[section_idx] ~
                get_surface_twist(psys, wing.idx, section_idx)
        ]

        # Get surface axes in body frame
        x_b = collect(
            get_surface_x_airf(psys, wing.idx, section_idx))
        y_b = collect(
            get_surface_y_airf(psys, wing.idx, section_idx))
        # z_b = x_b × y_b (normal)

        # Rotate x and z around y by twist (Rodrigues)
        cos_twist = cos(plate_twist[section_idx])
        sin_twist = sin(plate_twist[section_idx])
        x_twisted = cos_twist * x_b + sin_twist * (y_b × x_b)
        z_twisted = x_twisted × y_b

        # Transform to world frame
        eqs = [
            eqs
            plate_x_w[:, section_idx] ~ Rbw * x_twisted
            plate_y_w[:, section_idx] ~ Rbw * y_b
            plate_z_w[:, section_idx] ~ Rbw * z_twisted
        ]

        # Apparent wind at surface point
        wind_at_h = calc_wind_factor(
            s.am,
            pos[1, pidx], pos[2, pidx],
            pos[3, pidx], psys) * wind_vec_gnd
        eqs = [
            eqs
            plate_va_w[:, section_idx] ~
                wind_at_h - vel[:, pidx]
        ]

        # AoA from wind projection (VSM convention)
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

        # CL/CD lookup
        eqs = [
            eqs
            plate_cl[section_idx] ~ get_plate_cl(
                psys, wing.idx, plate_alpha[section_idx])
            plate_cd[section_idx] ~
                get_plate_drag_corr(psys, wing.idx) *
                get_plate_cd(
                    psys, wing.idx, plate_alpha[section_idx])
        ]

        # Dynamic pressure: in-plane for lift, full for drag
        apparent_wind = collect(plate_va_w[:, section_idx])
        eqs = [
            eqs
            plate_q[section_idx] ~ 0.5 *
                calc_rho(s.am, height[pidx]) *
                (plate_v_tan[section_idx]^2 + plate_v_norm[section_idx]^2)
            plate_q_drag[section_idx] ~ 0.5 *
                calc_rho(s.am, height[pidx]) *
                (apparent_wind ⋅ apparent_wind)
        ]

        # Lift and drag directions (VSM convention):
        # Effective flow in airfoil plane, then lift
        # perpendicular and drag along flow.
        alpha_rad = atan(plate_v_norm[section_idx], plate_v_tan[section_idx])
        va_airf_dir = cos(alpha_rad) * x_axis_w +
                      sin(alpha_rad) * z_axis_w
        y_axis_w = collect(plate_y_w[:, section_idx])
        lift_dir = smooth_normalize(va_airf_dir × y_axis_w)
        drag_dir = smooth_normalize(y_axis_w × lift_dir)

        area = get_surface_area(psys, wing.idx, section_idx)
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

    # Apply forces depending on wing type
    if wing.dynamics_type == PARTICLE_DYNAMICS
        # Per-point forces in body frame
        aero_force_point = aero_force_point_b::AbstractArray
        for (section_idx, surf) in enumerate(surfaces)
            pidx = surf.point_idx
            eqs = [
                eqs
                aero_force_point[:, pidx] ~
                    Rbw' * plate_force_w[:, section_idx]
            ]
        end
        # Wing-level: sum of all surface forces (body frame)
        eqs = [
            eqs
            aero_force_b[:, wing.idx] ~
                sum([Rbw' * plate_force_w[:, section_idx]
                     for section_idx in 1:n_surf])
            aero_moment_b[:, wing.idx] ~ zeros(3)
        ]
    elseif wing.dynamics_type == RIGID_DYNAMICS
        # Sum forces and moments about COM
        force_sum = sum([
            Rbw' * plate_force_w[:, section_idx]
            for section_idx in 1:n_surf])
        moment_sum = sum([
            Rbw' * ((pos[:, surfaces[section_idx].point_idx] -
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
