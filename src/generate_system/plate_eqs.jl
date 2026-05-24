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

    for (si, surf) in enumerate(surfaces)
        pidx = surf.point_idx

        # Twist from surface
        eqs = [
            eqs
            plate_twist[si] ~
                get_surface_twist(psys, wing.idx, si)
        ]

        # Get surface axes in body frame
        x_b = collect(
            get_surface_x_airf(psys, wing.idx, si))
        y_b = collect(
            get_surface_y_airf(psys, wing.idx, si))
        # z_b = x_b × y_b (normal)

        # Rotate x and z around y by twist (Rodrigues)
        ct = cos(plate_twist[si])
        st = sin(plate_twist[si])
        x_twisted = ct * x_b + st * (y_b × x_b)
        z_twisted = x_twisted × y_b

        # Transform to world frame
        eqs = [
            eqs
            plate_x_w[:, si] ~ Rbw * x_twisted
            plate_y_w[:, si] ~ Rbw * y_b
            plate_z_w[:, si] ~ Rbw * z_twisted
        ]

        # Apparent wind at surface point
        wind_at_h = calc_wind_factor(
            s.am,
            pos[1, pidx], pos[2, pidx],
            pos[3, pidx], psys) * wind_vec_gnd
        eqs = [
            eqs
            plate_va_w[:, si] ~
                wind_at_h - vel[:, pidx]
        ]

        # AoA from wind projection (VSM convention)
        va = collect(plate_va_w[:, si])
        xw = collect(plate_x_w[:, si])
        zw = collect(plate_z_w[:, si])
        eqs = [
            eqs
            plate_v_tan[si] ~ va ⋅ xw
            plate_v_norm[si] ~ va ⋅ zw
            plate_alpha[si] ~
                rad2deg(atan(plate_v_norm[si],
                             plate_v_tan[si]))
        ]

        # CL/CD lookup
        eqs = [
            eqs
            plate_cl[si] ~ get_plate_cl(
                psys, wing.idx, plate_alpha[si])
            plate_cd[si] ~
                get_plate_drag_corr(psys, wing.idx) *
                get_plate_cd(
                    psys, wing.idx, plate_alpha[si])
        ]

        # Dynamic pressure: in-plane for lift, full for drag
        va = collect(plate_va_w[:, si])
        eqs = [
            eqs
            plate_q[si] ~ 0.5 *
                calc_rho(s.am, height[pidx]) *
                (plate_v_tan[si]^2 + plate_v_norm[si]^2)
            plate_q_drag[si] ~ 0.5 *
                calc_rho(s.am, height[pidx]) *
                (va ⋅ va)
        ]

        # Lift and drag directions (VSM convention):
        # Effective flow in airfoil plane, then lift
        # perpendicular and drag along flow.
        alpha_rad = atan(plate_v_norm[si], plate_v_tan[si])
        va_airf_dir = cos(alpha_rad) * xw +
                      sin(alpha_rad) * zw
        yw = collect(plate_y_w[:, si])
        lift_dir = smooth_normalize(va_airf_dir × yw)
        drag_dir = smooth_normalize(yw × lift_dir)

        area = get_surface_area(psys, wing.idx, si)
        eqs = [
            eqs
            plate_lift[:, si] ~
                plate_q[si] * area *
                plate_cl[si] * lift_dir
            plate_drag[:, si] ~
                plate_q_drag[si] * area *
                plate_cd[si] * drag_dir
            plate_force_w[:, si] ~
                plate_lift[:, si] + plate_drag[:, si]
        ]
    end

    # Apply forces depending on wing type
    if wing.dynamics_type == PARTICLE_DYNAMICS
        # Per-point forces in body frame
        afpb = aero_force_point_b::AbstractArray
        for (si, surf) in enumerate(surfaces)
            pidx = surf.point_idx
            eqs = [
                eqs
                afpb[:, pidx] ~
                    Rbw' * plate_force_w[:, si]
            ]
        end
        # Wing-level: sum of all surface forces (body frame)
        eqs = [
            eqs
            aero_force_b[:, wing.idx] ~
                sum([Rbw' * plate_force_w[:, si]
                     for si in 1:n_surf])
            aero_moment_b[:, wing.idx] ~ zeros(3)
        ]
    elseif wing.dynamics_type == RIGID_DYNAMICS
        # Sum forces and moments about COM
        force_sum = sum([
            Rbw' * plate_force_w[:, si]
            for si in 1:n_surf])
        moment_sum = sum([
            Rbw' * ((pos[:, surfaces[si].point_idx] -
                     com_w[:, wing.idx]) ×
                    plate_force_w[:, si])
            for si in 1:n_surf])
        eqs = [
            eqs
            aero_force_b[:, wing.idx] ~ force_sum
            aero_moment_b[:, wing.idx] ~ moment_sum
        ]
    end

    return eqs
end
