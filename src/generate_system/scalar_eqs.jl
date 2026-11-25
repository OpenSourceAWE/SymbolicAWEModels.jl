# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Scalar kinematic equation generation

"""
    scalar_eqs!(s, eqs, params; kwargs...)

Generate equations for the derived scalar kinematics used in control and analysis:
elevation, azimuth, heading, course, angle of attack, their time derivatives, and the
apparent wind.

# Arguments
- `s::SymbolicAWEModel`: The run.
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
    wings = s.definition.wings
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
    for wing_idx in eachindex(wings)
        body_idx = wing_body(s, wing_idx)
        eqs = [
            eqs
            e_x[:, wing_idx] ~ R_b_to_w[:, 1, body_idx]
            e_y[:, wing_idx] ~ R_b_to_w[:, 2, body_idx]
            e_z[:, wing_idx] ~ R_b_to_w[:, 3, body_idx]
            wind_vel_wing[:, wing_idx] ~
                wing_wind_source(params, wing_idx, wind_vec_gnd)(
                    wing_pos[3, body_idx])
            wind_disturb[:, wing_idx] ~ params.wings[wing_idx].wind_disturb
            va_wing[:, wing_idx] ~
                wind_vel_wing[:, wing_idx] - wing_vel[:, body_idx] +
                wind_disturb[:, wing_idx]
            va_wing_b[:, wing_idx] ~ R_b_to_w[:, :, body_idx]' * va_wing[:, wing_idx]
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

    for wing_idx in eachindex(wings)
        body_idx = wing_body(s, wing_idx)
        # Spherical coords are centered at the transform base point, not origin.
        base_point = transform_base_point(s.transforms, body_idx)
        rel_pos = base_point == 0 ? wing_pos[:, body_idx] :
            wing_pos[:, body_idx] .- pos[:, base_point]

        scalars = wing_scalar_kinematics(;
            rel_pos,
            e_x = e_x[:, wing_idx],
            R_t_to_w = R_t_to_w[:, :, wing_idx],
            R_v_to_w = R_v_to_w[:, :, wing_idx],
            R_b_to_w = R_b_to_w[:, :, body_idx],
            vel = wing_vel[:, body_idx],
            acc = wing_acc[:, body_idx],
            omega_b = ω_b[:, body_idx],
            alpha_b = α_b[:, body_idx],
            va_b = va_wing_b[:, wing_idx],
            twist_offset = mid_span_twist(i -> twist_angle[i],
                                          wing_station_idxs(s, wing_idx)))

        eqs = [
            eqs
            vec(R_v_to_w[:, :, wing_idx]) ~
                vec(calc_R_v_to_w(rel_pos, e_x[:, wing_idx]))
            vec(R_t_to_w[:, :, wing_idx]) ~ vec(sym_calc_R_t_to_w(rel_pos))
            heading[wing_idx] ~ scalars.heading
            turn_rate[:, wing_idx] ~ scalars.turn_rate
            turn_acc[:, wing_idx] ~ scalars.turn_acc
            distance[wing_idx] ~ scalars.distance
            distance_vel[wing_idx] ~ scalars.distance_vel
            distance_acc[wing_idx] ~ scalars.distance_acc
            elevation[wing_idx] ~ scalars.elevation
            elevation_vel[wing_idx] ~ scalars.elevation_vel
            elevation_acc[wing_idx] ~ scalars.elevation_acc
            azimuth[wing_idx] ~ scalars.azimuth
            azimuth_vel[wing_idx] ~ scalars.azimuth_vel
            azimuth_acc[wing_idx] ~ scalars.azimuth_acc
            course[wing_idx] ~ scalars.course
            angle_of_attack[wing_idx] ~ scalars.angle_of_attack
        ]
    end
    return eqs
end
