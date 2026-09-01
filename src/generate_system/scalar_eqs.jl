# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Scalar kinematic equation generation

"""
    scalar_eqs!(s, eqs, params; kwargs...)

Generate equations for the derived scalar kinematics used in control and analysis:
elevation, azimuth, heading, course, angle of attack, their time derivatives, and the
apparent wind.

# Arguments
- `s::SymbolicAWEModel`: The main model object.
- `eqs`: Accumulating equation vector.
- `kwargs...`: Symbolic variables for the system's state.

# Returns
- `eqs`: The updated list of system equations.
"""
function scalar_eqs!(
    s, eqs, params;
    R_b_to_w, wind_vec_gnd, va_wing_b, wing_pos,
    wing_vel, wing_acc, twist_angle, ω_b, α_b,
    R_v_to_w, pos
)
    (; wings) = s.sys_struct
    wind_factor = param_computed!(params.reg, :wind_factor, WindFactorReader())
    @variables begin
        # Body frame axes and apparent wind (column-major: [1:3, wing_idx])
        e_x(t)[1:3, eachindex(wings)]
        e_y(t)[1:3, eachindex(wings)]
        e_z(t)[1:3, eachindex(wings)]
        wind_vel_wing(t)[1:3, eachindex(wings)]
        wind_disturb(t)[1:3, eachindex(wings)]
        va_wing(t)[1:3, eachindex(wings)]
    end
    eqs = [eqs; collect(wind_vec_gnd) .~ ground_wind_vec(params)]
    for wing in wings
        eqs = [
            eqs
            e_x[:, wing.idx] ~ R_b_to_w[:, 1, wing.idx]
            e_y[:, wing.idx] ~ R_b_to_w[:, 2, wing.idx]
            e_z[:, wing.idx] ~ R_b_to_w[:, 3, wing.idx]
            wind_vel_wing[:, wing.idx] ~
                wind_factor(wing_pos[3, wing.idx]) * wind_vec_gnd
            wind_disturb[:, wing.idx] ~ params.wings[wing.idx].wind_disturb
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
        # Spherical coords are centered at the transform base point, not origin.
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

        has_stations = !isempty(wing.station_idxs)
        half_len = has_stations ?
            wing.station_idxs[1] +
            length(wing.station_idxs) ÷ 2 - 1 : 0
        twist_offset = has_stations ?
            0.5 * twist_angle[half_len] + 0.5 * twist_angle[half_len + 1] : 0

        scalars = wing_scalar_kinematics(;
            rel_pos,
            e_x = e_x[:, wing.idx],
            R_t_to_w = R_t_to_w[:, :, wing.idx],
            R_v_to_w = R_v_to_w[:, :, wing.idx],
            R_b_to_w = R_b_to_w[:, :, wing.idx],
            vel = wing_vel[:, wing.idx],
            acc = wing_acc[:, wing.idx],
            omega_b = ω_b[:, wing.idx],
            alpha_b = α_b[:, wing.idx],
            va_b = va_wing_b[:, wing.idx],
            twist_offset)

        eqs = [
            eqs
            vec(R_v_to_w[:, :, wing.idx]) ~
                vec(calc_R_v_to_w(rel_pos, e_x[:, wing.idx]))
            vec(R_t_to_w[:, :, wing.idx]) ~ vec(sym_calc_R_t_to_w(rel_pos))
            heading[wing.idx] ~ scalars.heading
            turn_rate[:, wing.idx] ~ scalars.turn_rate
            turn_acc[:, wing.idx] ~ scalars.turn_acc
            distance[wing.idx] ~ scalars.distance
            distance_vel[wing.idx] ~ scalars.distance_vel
            distance_acc[wing.idx] ~ scalars.distance_acc
            elevation[wing.idx] ~ scalars.elevation
            elevation_vel[wing.idx] ~ scalars.elevation_vel
            elevation_acc[wing.idx] ~ scalars.elevation_acc
            azimuth[wing.idx] ~ scalars.azimuth
            azimuth_vel[wing.idx] ~ scalars.azimuth_vel
            azimuth_acc[wing.idx] ~ scalars.azimuth_acc
            course[wing.idx] ~ scalars.course
            angle_of_attack[wing.idx] ~ scalars.angle_of_attack
        ]
    end
    return eqs
end
