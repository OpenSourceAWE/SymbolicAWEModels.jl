# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Scalar kinematic equation generation

"""
    scalar_eqs!(s, eqs, psys; kwargs...)

Generate equations for derived scalar kinematic quantities useful for control and
analysis.

This includes elevation, azimuth, heading, course, angle of attack, and their time
derivatives, as well as apparent wind calculations.

# Arguments
- `s::SymbolicAWEModel`: The main model object.
- `eqs`, `psys`: Accumulating vectors and symbolic parameter.
- `kwargs...`: Symbolic variables for the system's state.

# Returns
- `eqs`: The updated list of system equations.
"""
function scalar_eqs!(
    s, eqs, psys;
    R_b_to_w, wind_vec_gnd, va_wing_b, wing_pos,
    wing_vel, wing_acc, twist_angle, ω_b, α_b,
    R_v_to_w, pos
)
    (; wings) = s.sys_struct
    @variables begin
        # Body frame axes and apparent wind (column-major: [1:3, wing_idx])
        e_x(t)[1:3, eachindex(wings)]
        e_y(t)[1:3, eachindex(wings)]
        e_z(t)[1:3, eachindex(wings)]
        wind_vel_wing(t)[1:3, eachindex(wings)]
        wind_disturb(t)[1:3, eachindex(wings)]
        va_wing(t)[1:3, eachindex(wings)]
    end
    eqs = [
        eqs
        wind_vec_gnd ~ get_wind_vec(psys)
    ]
    for wing in wings
        eqs = [
            eqs
            e_x[:, wing.idx] ~ R_b_to_w[:, 1, wing.idx]
            e_y[:, wing.idx] ~ R_b_to_w[:, 2, wing.idx]
            e_z[:, wing.idx] ~ R_b_to_w[:, 3, wing.idx]
            wind_vel_wing[:, wing.idx] ~
                calc_wind_factor(s.am, wing_pos[1, wing.idx], wing_pos[2, wing.idx],
                                 wing_pos[3, wing.idx], psys) * wind_vec_gnd
            wind_disturb[:, wing.idx] ~ get_wind_disturb(psys, wing.idx)
            va_wing[:, wing.idx] ~
                wind_vel_wing[:, wing.idx] - wing_vel[:, wing.idx] +
                wind_disturb[:, wing.idx]
            va_wing_b[:, wing.idx] ~ R_b_to_w[:, :, wing.idx]' * va_wing[:, wing.idx]
        ]
    end
    @variables begin
        # Kinematic quantities
        heading(t)[eachindex(wings)]
        turn_rate(t)[1:3, eachindex(wings)]
        turn_acc(t)[1:3, eachindex(wings)]
        azimuth(t)[eachindex(wings)]
        azimuth_vel(t)[eachindex(wings)]
        azimuth_acc(t)[eachindex(wings)]
        elevation(t)[eachindex(wings)]
        elevation_vel(t)[eachindex(wings)]
        elevation_acc(t)[eachindex(wings)]
        course(t)[eachindex(wings)]
        angle_of_attack(t)[eachindex(wings)]
        R_t_to_w(t)[1:3, 1:3, eachindex(wings)]
        distance(t)[eachindex(wings)]
        distance_vel(t)[eachindex(wings)]
        distance_acc(t)[eachindex(wings)]
    end

    for wing in wings
        # Compute position relative to transform base point.
        # Spherical coordinates (heading, elevation, azimuth,
        # distance, R_t_to_w) are defined on the sphere centered
        # at the base, not at the world origin.
        transforms = s.sys_struct.transforms
        if wing.transform_idx != 0 &&
                wing.transform_idx <= length(transforms)
            tf = transforms[wing.transform_idx]
            bp_idx = tf.base_point_idx
            rel_pos = wing_pos[:, wing.idx] .-
                pos[:, bp_idx]
        else
            rel_pos = wing_pos[:, wing.idx]
        end

        x, y, _ = rel_pos
        has_twist_surfaces = !isempty(wing.twist_surface_idxs)
        half_len = 0
        if has_twist_surfaces
            half_len = wing.twist_surface_idxs[1] +
                length(wing.twist_surface_idxs) ÷ 2 - 1
        end

        # Calculate heading using tangential sphere frame.
        # Projects e_x onto the tangent plane via R_t_to_w:
        # x_t = elevation dir (away from zenith),
        # y_t = azimuthal.
        # heading = 0 when e_x aligns with x_t (nose toward
        # GS).
        heading_t_1 = e_x[:, wing.idx] ⋅
            R_t_to_w[:, 1, wing.idx]
        heading_t_2 = e_x[:, wing.idx] ⋅
            R_t_to_w[:, 2, wing.idx]
        # Course: velocity direction in same tangential frame.
        course_t_1 = wing_vel[:, wing.idx] ⋅
            R_t_to_w[:, 1, wing.idx]
        course_t_2 = wing_vel[:, wing.idx] ⋅
            R_t_to_w[:, 2, wing.idx]

        # Unified equations for both RIGID_DYNAMICS and PARTICLE_DYNAMICS
        # wings. PARTICLE_DYNAMICS wings have ω_b=α_b=0 (set in
        # wing_eqs!), so turn_rate/turn_acc naturally evaluate
        # to zero.
        eqs = [
            eqs
            vec(R_v_to_w[:, :, wing.idx]) ~
                vec(calc_R_v_to_w(
                    rel_pos, e_x[:, wing.idx]))
            vec(R_t_to_w[:, :, wing.idx]) ~
                vec(sym_calc_R_t_to_w(rel_pos))
            heading[wing.idx] ~
                atan(heading_t_2, heading_t_1)
            turn_rate[:, wing.idx] ~
                R_v_to_w[:, :, wing.idx]' *
                (R_b_to_w[:, :, wing.idx] *
                    ω_b[:, wing.idx])
            turn_acc[:, wing.idx] ~
                R_v_to_w[:, :, wing.idx]' *
                (R_b_to_w[:, :, wing.idx] *
                    α_b[:, wing.idx])
            distance[wing.idx] ~ smooth_norm(rel_pos)
            distance_vel[wing.idx] ~
                wing_vel[:, wing.idx] ⋅
                    R_t_to_w[:, 3, wing.idx]
            distance_acc[wing.idx] ~
                wing_acc[:, wing.idx] ⋅
                    R_t_to_w[:, 3, wing.idx]
            elevation[wing.idx] ~
                KiteUtils.calc_elevation(rel_pos)
            elevation_vel[wing.idx] ~
                dot(wing_vel[:, wing.idx],
                    -R_t_to_w[:, 1, wing.idx]) /
                distance[wing.idx]
            elevation_acc[wing.idx] ~
                dot(wing_acc[:, wing.idx],
                    -R_t_to_w[:, 1, wing.idx]) /
                distance[wing.idx]
            azimuth[wing.idx] ~
                KiteUtils.azimuth_east(rel_pos)
            azimuth_vel[wing.idx] ~
                dot(wing_vel[:, wing.idx],
                    -R_t_to_w[:, 2, wing.idx]) /
                smooth_norm([x, y])
            azimuth_acc[wing.idx] ~
                dot(wing_acc[:, wing.idx],
                    -R_t_to_w[:, 2, wing.idx]) /
                smooth_norm([x, y])
            course[wing.idx] ~
                atan(course_t_2, course_t_1)
            angle_of_attack[wing.idx] ~
                calc_angle_of_attack(
                    va_wing_b[:, wing.idx]) +
                (has_twist_surfaces ?
                    0.5 * twist_angle[half_len] +
                    0.5 * twist_angle[half_len + 1] :
                    0)
        ]
    end
    return eqs
end
