# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

# Implementation of the ram air wing model using ModelingToolkit.jl

"""
    calc_speed_acc(winch::AsyncMachine, tether_vel, norm_, set_speed)

Calculate winch acceleration for a speed-controlled asynchronous machine model.
"""
function calc_speed_acc(winch::AsyncMachine, tether_vel, norm_, set_speed)
    calc_acceleration(winch, tether_vel, norm_; set_speed, set_torque=nothing, use_brake=false) # TODO: add brake setting
end

"""
    calc_moment_acc(winch::TorqueControlledMachine, tether_vel, norm_, set_torque)

Calculate winch acceleration for a torque-controlled machine model.
"""
function calc_moment_acc(winch::TorqueControlledMachine, tether_vel, norm_, set_torque)
    calc_acceleration(winch, tether_vel, norm_; set_speed=nothing, set_torque, use_brake=false)
end

"""
    calc_angle_of_attack(va_wing_b)

Calculate the angle of attack [rad] from the apparent wind vector in the body frame.
"""
function calc_angle_of_attack(va_wing_b)
    return atan(va_wing_b[3], va_wing_b[1])
end

"""
    sym_normalize(vec)

Symbolic-safe normalization of a vector.
"""
function sym_normalize(vec)
    return vec / norm(vec)
end

"""
    quaternion_to_rotation_matrix(q)

Convert a quaternion `q` (scalar-first format [w, x, y, z]) to a 3x3 rotation matrix.
"""
function quaternion_to_rotation_matrix(q)
    w, x, y, z = q[1], q[2], q[3], q[4]
    
    return [
        1 - 2*(y*y + z*z)  2*(x*y - z*w)      2*(x*z + y*w);
        2*(x*y + z*w)      1 - 2*(x*x + z*z)  2*(y*z - x*w);
        2*(x*z - y*w)      2*(y*z + x*w)      1 - 2*(x*x + y*y)
    ]
end

"""
    rotation_matrix_to_quaternion(R)

Convert a 3x3 rotation matrix `R` to a quaternion (scalar-first format [w, x, y, z]).
"""
function rotation_matrix_to_quaternion(R)
    tr_ = tr(R)
    
    if tr_ > 0
        S = sqrt(tr_ + 1.0) * 2
        w = 0.25 * S
        x = (R[3,2] - R[2,3]) / S
        y = (R[1,3] - R[3,1]) / S
        z = (R[2,1] - R[1,2]) / S
    elseif (R[1,1] > R[2,2]) && (R[1,1] > R[3,3])
        S = sqrt(1.0 + R[1,1] - R[2,2] - R[3,3]) * 2
        w = (R[3,2] - R[2,3]) / S
        x = 0.25 * S
        y = (R[1,2] + R[2,1]) / S
        z = (R[1,3] + R[3,1]) / S
    elseif R[2,2] > R[3,3]
        S = sqrt(1.0 + R[2,2] - R[1,1] - R[3,3]) * 2
        w = (R[1,3] - R[3,1]) / S
        x = (R[1,2] + R[2,1]) / S
        y = 0.25 * S
        z = (R[2,3] + R[3,2]) / S
    else
        S = sqrt(1.0 + R[3,3] - R[1,1] - R[2,2]) * 2
        w = (R[2,1] - R[1,2]) / S
        x = (R[1,3] + R[3,1]) / S
        y = (R[2,3] + R[3,2]) / S
        z = 0.25 * S
    end
    
    return [w, x, y, z]
end

# The following functions are registered for use within ModelingToolkit.jl's symbolic context.
# They act as symbolic placeholders for accessing fields from the SystemStructure (`psys`)
# and Settings (`pset`) parameter objects during equation generation.
get_pos_w(sys::SystemStructure, idx::Int16) = sys.points[idx].pos_w
@register_array_symbolic get_pos_w(sys::SystemStructure, idx::Int16) begin
    size=(3,)
    eltype=SimFloat
end
get_vel_w(sys::SystemStructure, idx::Int16) = sys.points[idx].vel_w
@register_array_symbolic get_vel_w(sys::SystemStructure, idx::Int16) begin
    size=(3,)
    eltype=SimFloat
end
get_pos_b(sys::SystemStructure, idx::Int16) = sys.points[idx].pos_b
@register_array_symbolic get_pos_b(sys::SystemStructure, idx::Int16) begin
    size=(3,)
    eltype=SimFloat
end
get_wing_pos_w(sys::SystemStructure, idx::Int16) = sys.wings[idx].pos_w
@register_array_symbolic get_wing_pos_w(sys::SystemStructure, idx::Int16) begin
    size=(3,)
    eltype=SimFloat
end
get_wing_vel_w(sys::SystemStructure, idx::Int16) = sys.wings[idx].vel_w
@register_array_symbolic get_wing_vel_w(sys::SystemStructure, idx::Int16) begin
    size=(3,)
    eltype=SimFloat
end
get_Q_b_w(sys::SystemStructure, idx::Int16) = sys.wings[idx].Q_b_w
@register_array_symbolic get_Q_b_w(sys::SystemStructure, idx::Int16) begin
    size=(4,)
    eltype=SimFloat
end
get_ω_b(sys::SystemStructure, idx::Int16) = sys.wings[idx].ω_b
@register_array_symbolic get_ω_b(sys::SystemStructure, idx::Int16) begin
    size=(3,)
    eltype=SimFloat
end
get_wind_disturb(sys::SystemStructure, idx::Int16) = sys.wings[idx].wind_disturb
@register_array_symbolic get_wind_disturb(sys::SystemStructure, idx::Int16) begin
    size=(3,)
    eltype=SimFloat
end
get_disturb(sys::SystemStructure, idx::Int16) = sys.points[idx].disturb
@register_array_symbolic get_disturb(sys::SystemStructure, idx::Int16) begin
    size=(3,)
    eltype=SimFloat
end
get_vsm_y(sys::SystemStructure, idx::Int16, iy::Int) = sys.wings[idx].vsm_y[iy]
@register_symbolic get_vsm_y(sys::SystemStructure, idx::Int16, iy::Int)
get_vsm_x(sys::SystemStructure, idx::Int16, ix::Int) = sys.wings[idx].vsm_x[ix]
@register_symbolic get_vsm_x(sys::SystemStructure, idx::Int16, ix::Int)
get_vsm_jac(sys::SystemStructure, idx::Int16, ix::Int, iy::Int) = sys.wings[idx].vsm_jac[ix,iy]
@register_symbolic get_vsm_jac(sys::SystemStructure, idx::Int16, ix::Int, iy::Int)
get_pulley_len(sys::SystemStructure, idx::Int16) = sys.pulleys[idx].len
@register_symbolic get_pulley_len(sys::SystemStructure, idx::Int16)
get_pulley_vel(sys::SystemStructure, idx::Int16) = sys.pulleys[idx].vel
@register_symbolic get_pulley_vel(sys::SystemStructure, idx::Int16)
get_set_value(sys::SystemStructure, idx::Int16) = sys.winches[idx].set_value
@register_symbolic get_set_value(sys::SystemStructure, idx::Int16)
get_twist(sys::SystemStructure, idx::Int16) = sys.groups[idx].twist
@register_symbolic get_twist(sys::SystemStructure, idx::Int16)
get_twist_ω(sys::SystemStructure, idx::Int16) = sys.groups[idx].twist_ω
@register_symbolic get_twist_ω(sys::SystemStructure, idx::Int16)
get_mass(sys::SystemStructure, idx::Int16) = sys.points[idx].mass
@register_symbolic get_mass(sys::SystemStructure, idx::Int16)
get_l0(sys::SystemStructure, idx::Int16) = sys.segments[idx].l0
@register_symbolic get_l0(sys::SystemStructure, idx::Int16)
get_diameter(sys::SystemStructure, idx::Int16) = sys.segments[idx].diameter
@register_symbolic get_diameter(sys::SystemStructure, idx::Int16)
get_compression_frac(sys::SystemStructure, idx::Int16) = sys.segments[idx].compression_frac
@register_symbolic get_compression_frac(sys::SystemStructure, idx::Int16)
get_moment_frac(sys::SystemStructure, idx::Int16) = sys.groups[idx].moment_frac
@register_symbolic get_moment_frac(sys::SystemStructure, idx::Int16)
get_sum_len(sys::SystemStructure, idx::Int16) = sys.pulleys[idx].sum_len
@register_symbolic get_sum_len(sys::SystemStructure, idx::Int16)
get_tether_len(sys::SystemStructure, idx::Int16) = sys.winches[idx].tether_len
@register_symbolic get_tether_len(sys::SystemStructure, idx::Int16)
get_tether_vel(sys::SystemStructure, idx::Int16) = sys.winches[idx].tether_vel
@register_symbolic get_tether_vel(sys::SystemStructure, idx::Int16)
get_axial_stiffness(sys::SystemStructure, idx::Int16) = sys.segments[idx].axial_stiffness
@register_symbolic get_axial_stiffness(sys::SystemStructure, idx::Int16)
get_axial_damping(sys::SystemStructure, idx::Int16) = sys.segments[idx].axial_damping
@register_symbolic get_axial_damping(sys::SystemStructure, idx::Int16)
get_bridle_damping(sys::SystemStructure, idx::Int16) = sys.points[idx].bridle_damping
@register_symbolic get_bridle_damping(sys::SystemStructure, idx::Int16)
get_brake(sys::SystemStructure, idx::Int16) = sys.winches[idx].brake
@register_symbolic get_brake(sys::SystemStructure, idx::Int16)
get_fix_point_sphere(sys::SystemStructure, idx::Int16) = sys.points[idx].fix_sphere
@register_symbolic get_fix_point_sphere(sys::SystemStructure, idx::Int16)
get_fix_wing_sphere(sys::SystemStructure, idx::Int16) = sys.wings[idx].fix_sphere
@register_symbolic get_fix_wing_sphere(sys::SystemStructure, idx::Int16)

get_set_mass(set::Settings) = set.mass
@register_symbolic get_set_mass(set::Settings)
get_rho_tether(set::Settings) = set.rho_tether
@register_symbolic get_rho_tether(set::Settings)
get_e_tether(set::Settings) = set.e_tether
@register_symbolic get_e_tether(set::Settings)
get_damping(set::Settings) = set.damping
@register_symbolic get_damping(set::Settings)
get_c_spring(set::Settings) = set.c_spring
@register_symbolic get_c_spring(set::Settings)
get_cd_tether(set::Settings) = set.cd_tether
@register_symbolic get_cd_tether(set::Settings)
get_v_wind(set::Settings) = set.v_wind
@register_symbolic get_v_wind(set::Settings)
get_upwind_dir(set::Settings) = set.upwind_dir
@register_symbolic get_upwind_dir(set::Settings)
get_g_earth(set::Settings) = set.g_earth
@register_symbolic get_g_earth(set::Settings)

"""
    force_eqs!(s, system, psys, pset, eqs, defaults, guesses; R_b_w, wing_pos, ...)

Generate the force and constraint equations for the mass-spring-damper components.

This function builds the core equations for:
- **Points**: Newton's second law for dynamic points, force balance for quasi-static
  points, and kinematic constraints for points attached to the wing.
- **Segments**: Spring-damper forces (Hooke's law) and aerodynamic drag forces.
- **Pulleys**: Length redistribution dynamics or constraints.
- **Winches**: Tether length and velocity dynamics based on the winch model.
- **Groups**: Rotational dynamics for wing twist deformation.

# Arguments
- `s::SymbolicAWEModel`: The main model object.
- `system::SystemStructure`: The physical structure definition.
- `psys`, `pset`: Symbolic parameters representing `system` and `s.set`.
- `eqs`, `defaults`, `guesses`: The accumulating vectors for the MTK system.
- `R_b_w`, `wing_pos`, etc.: Symbolic variables for the wing's state.

# Returns
- `(eqs, defaults, guesses, tether_wing_force, tether_wing_moment)`: A tuple containing
  the updated equation lists and the calculated aggregate forces and moments exerted
  by the tethers on the wing.
"""
function force_eqs!(s, system, psys, pset, eqs, defaults, guesses; 
                    R_b_w, wing_pos, wing_vel, wind_vec_gnd, group_aero_moment, 
                    twist_angle, twist_ω, set_values, fix_wing)

    @unpack points, groups, segments, pulleys, tethers, winches, wings = system
    
    # ==================== POINTS ==================== #
    tether_wing_force = zeros(Num, length(wings), 3)
    tether_wing_moment = zeros(Num, length(wings), 3)
    @variables begin
        pos(t)[1:3, eachindex(points)]
        vel(t)[1:3, eachindex(points)]
        acc(t)[1:3, eachindex(points)]
        point_force(t)[1:3, eachindex(points)]
        disturb_force(t)[1:3, eachindex(points)]
        tether_r(t)[1:3, eachindex(points)]
        point_mass(t)[eachindex(points)]
        chord_b(t)[1:3, eachindex(points)]
        normal(t)[1:3, eachindex(points)]
        pos_b(t)[1:3, eachindex(points)]
        fix_point_sphere(t)[eachindex(points)]

        spring_force_vec(t)[1:3, eachindex(segments)]
        drag_force(t)[1:3, eachindex(segments)]
        l0(t)[eachindex(segments)]
    end
    for point in points
        F::Vector{Num} = zeros(Num, 3)
        mass = get_mass(psys, point.idx)
        for segment in segments
            if point.idx in segment.point_idxs
                mass_per_meter = get_rho_tether(pset) * π * (get_diameter(psys, segment.idx)/2)^2
                inverted = segment.point_idxs[2] == point.idx
                if inverted
                    F .-= spring_force_vec[:, segment.idx]
                else
                    F .+= spring_force_vec[:, segment.idx]
                end
                mass += mass_per_meter * l0[segment.idx] / 2
                F .+= 0.5drag_force[:, segment.idx]
            end
        end

        eqs = [
            eqs
            point_mass[point.idx] ~ mass
            disturb_force[:, point.idx] ~ get_disturb(psys, point.idx)
            point_force[:, point.idx]  ~ F + [0, 0, -get_g_earth(pset) * mass] + 
                                            disturb_force[:, point.idx]
        ]

        if point.type == WING
            found = 0
            group = nothing
            for group_ in groups
                if point.idx in group_.point_idxs
                    group = group_
                    found += 1
                end
            end
            !(found in [0,1]) && error("Kite point number $(point.idx) is part of $found groups, 
                  and should be part of exactly 0 or 1 groups.")

            if found == 1
                found = 0
                wing = nothing
                for wing_ in wings
                    if group.idx in wing_.group_idxs
                        wing = wing_
                        found += 1
                    end
                end
                !(found == 1) && error("Kite group number $(group.idx) is part of $found wings, 
                      and should be part of exactly 1 wing.")

                fixed_pos = group.le_pos
                eqs = [
                    eqs
                    chord_b[:, point.idx]   ~ get_pos_b(psys, point.idx) .- fixed_pos
                    normal[:, point.idx]   ~ chord_b[:, point.idx] × group.y_airf
                    pos_b[:, point.idx]     ~ fixed_pos .+ cos(twist_angle[group.idx]) * chord_b[:, point.idx] - sin(twist_angle[group.idx]) * normal[:, point.idx]
                ]
            elseif found == 0
                eqs = [
                    eqs
                    pos_b[:, point.idx]     ~ get_pos_b(psys, point.idx)
                #     tether_r[:, point.idx]  ~ pos[:, point.idx] - wing_pos[point.wing_idx, :]
                ]
                # tether_wing_moment[point.wing_idx, :] .+= tether_r[:, point.idx] × point_force[:, point.idx]
            end
            eqs = [
                eqs
                # pos_b[:, point.idx]     ~ get_pos_b(psys, point.idx)
                tether_r[:, point.idx]  ~ pos[:, point.idx] - wing_pos[point.wing_idx, :]
            ]
            tether_wing_moment[point.wing_idx, :] .+= tether_r[:, point.idx] × point_force[:, point.idx]
            tether_wing_force[point.wing_idx, :] .+= point_force[:, point.idx]
            
            eqs = [
                eqs
                pos[:, point.idx]   ~ wing_pos[point.wing_idx, :] + R_b_w[point.wing_idx, :, :] * pos_b[:, point.idx]
                vel[:, point.idx]   ~ zeros(3)
                acc[:, point.idx]   ~ zeros(3)
            ]
        elseif point.type == STATIC
            eqs = [
                eqs
                pos[:, point.idx]   ~ get_pos_w(psys, point.idx)
                vel[:, point.idx]   ~ zeros(3)
                acc[:, point.idx]   ~ zeros(3)
            ]
        elseif point.type == DYNAMIC
            # p = pos[:, point.idx]
            # n = sym_normalize(wing_pos)
            # n = n * (p ⋅ n)
            # r = (p - n) # https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line#Vector_formulation
            if length(wings) > 0
                bridle_damp_vec = get_bridle_damping(psys, point.idx) * (vel[:, point.idx] - wing_vel[point.wing_idx, :])
            else
                bridle_damp_vec = zeros(3)
            end
            axis = sym_normalize(pos[:, point.idx])
            eqs = [
                eqs
                fix_point_sphere[point.idx] ~ get_fix_point_sphere(psys, point.idx)
                D(pos[:, point.idx]) ~  ifelse.(fix_point_sphere[point.idx]==true,
                                                vel[:, point.idx] ⋅ axis * axis,
                                                vel[:, point.idx]
                                        )
                D(vel[:, point.idx]) ~  ifelse.(fix_point_sphere[point.idx]==true,
                                                acc[:, point.idx] ⋅ axis * axis,
                                                acc[:, point.idx]
                                        )
                acc[:, point.idx]    ~ point_force[:, point.idx] ./ mass - bridle_damp_vec
            ]
            defaults = [
                defaults
                [pos[j, point.idx] => get_pos_w(psys, point.idx)[j] for j in 1:3]
                [vel[j, point.idx] => get_vel_w(psys, point.idx)[j] for j in 1:3]
            ]
        elseif point.type == QUASI_STATIC
            eqs = [
                eqs
                vel[:, point.idx]   ~ zeros(3)
                acc[:, point.idx]   ~ zeros(3)
                acc[:, point.idx]   ~ point_force[:, point.idx] / mass
            ]
            guesses = [
                guesses
                [acc[j, point.idx] => 0 for j in 1:3]
                [pos[j, point.idx] => get_pos_w(psys, point.idx)[j] for j in 1:3]
                [point_force[j, point.idx] => 0 for j in 1:3]
            ]
        else
            error("Unknown point type: $(typeof(point))")
        end
    end

    # ==================== GROUPS ==================== #
    if length(groups) > 0
        @variables begin
            trailing_edge_angle(t)[eachindex(groups)] # angle left / right
            trailing_edge_ω(t)[eachindex(groups)] # angular rate
            trailing_edge_α(t)[eachindex(groups)] # angular acc
            free_twist_angle(t)[eachindex(groups)]
            twist_α(t)[eachindex(groups)] # angular acc
            group_tether_force(t)[eachindex(groups)]
            group_tether_moment(t)[eachindex(groups)]
            tether_force(t)[eachindex(groups), eachindex(groups[1].point_idxs)]
            tether_moment(t)[eachindex(groups), eachindex(groups[1].point_idxs)]
            r_group(t)[eachindex(groups), eachindex(groups[1].point_idxs)]
            r_vec(t)[eachindex(groups), eachindex(groups[1].point_idxs), 1:3]
        end
    end
    
    for group in groups
        found = 0
        wing = nothing
        for wing_ in wings
            if group.idx in wing_.group_idxs
                wing = wing_
                found += 1
            end
        end
        !(found == 1) && error("Kite group number $(group.idx) is part of $found wings, 
              and should be part of exactly 1 wing.")

        all(iszero.(tether_wing_moment[wing.idx, :])) && 
            error("Tether wing moment is zero. At least one of the wing connection points should not be part of a deforming group.")

        x_airf = normalize(group.chord)
        init_z_airf = x_airf × group.y_airf
        z_airf = x_airf * sin(twist_angle[group.idx]) + init_z_airf * cos(twist_angle[group.idx])
        for (i, point_idx) in enumerate(group.point_idxs)
            eqs = [
                eqs
                r_vec[group.idx, i, :]      ~ (get_pos_b(psys, point_idx) .- (group.le_pos + get_moment_frac(psys, group.idx)*group.chord))
                r_group[group.idx, i]       ~ r_vec[group.idx, i, :] ⋅ normalize(group.chord)
                tether_force[group.idx, i]  ~ (point_force[:, point_idx] ⋅ (R_b_w[wing.idx, :, :] * -z_airf))
                tether_moment[group.idx, i] ~ r_group[group.idx, i] * tether_force[group.idx, i]
            ]
        end
        
        inertia = 1/3 * (get_set_mass(pset)/length(groups)) * (norm(group.chord))^2 # plate inertia around leading edge
        @parameters twist_damp = 50
        @parameters max_twist = deg2rad(90)

        eqs = [
            eqs
            group_tether_force[group.idx] ~ sum(tether_force[group.idx, :])
            group_tether_moment[group.idx] ~ sum(tether_moment[group.idx, :])
            twist_α[group.idx] ~ (group_aero_moment[group.idx] + group_tether_moment[group.idx]) / inertia
            twist_angle[group.idx] ~ clamp(free_twist_angle[group.idx], -max_twist, max_twist)
        ]
        if group.type == DYNAMIC
            eqs = [
                eqs
                D(free_twist_angle[group.idx]) ~ ifelse(fix_wing==true, 0, twist_ω[group.idx])
                D(twist_ω[group.idx]) ~ ifelse(fix_wing==true, 0, twist_α[group.idx] - twist_damp * twist_ω[group.idx])
            ]
            defaults = [
                defaults
                free_twist_angle[group.idx] => get_twist(psys, group.idx)
                twist_ω[group.idx] => get_twist_ω(psys, group.idx)
            ]
        elseif group.type == QUASI_STATIC
            eqs = [
                eqs
                twist_ω[group.idx] ~ 0
                twist_α[group.idx] ~ 0
            ]
            guesses = [
                guesses
                free_twist_angle[group.idx] => 0
                twist_angle[group.idx] => 0
            ]
        else
            error("Wrong group type.")
        end
    end

    # ==================== SEGMENTS ==================== #
    @variables begin
        segment_vec(t)[1:3, eachindex(segments)]
        unit_vec(t)[1:3, eachindex(segments)]
        len(t)[eachindex(segments)]
        rel_vel(t)[1:3, eachindex(segments)]
        spring_vel(t)[eachindex(segments)]
        spring_force(t)[eachindex(segments)]
        stiffness(t)[eachindex(segments)]
        damping(t)[eachindex(segments)]

        height(t)[eachindex(segments)]
        segment_vel(t)[1:3, eachindex(segments)]
        segment_rho(t)[eachindex(segments)]
        wind_vel(t)[1:3, eachindex(segments)]
        va(t)[1:3, eachindex(segments)]
        area(t)[eachindex(segments)]
        app_perp_vel(t)[1:3, eachindex(segments)]
        drag_force(t)[1:3, eachindex(segments)]

        pulley_len(t)[eachindex(pulleys)]

        tether_len(t)[eachindex(winches)]
    end
    for segment in segments
        p1, p2 = segment.point_idxs[1], segment.point_idxs[2]
        guesses = [
            guesses
            [segment_vec[i, segment.idx] => get_pos_w(psys, p2)[i] - get_pos_w(psys, p1)[i] for i in 1:3]
        ]

        in_pulley = 0
        for pulley in pulleys
            if segment.idx == pulley.segment_idxs[1] # each bridle segment has to be part of no pulley or one pulley
                eqs = [
                    eqs
                    l0[segment.idx] ~ pulley_len[pulley.idx]
                ]
                in_pulley += 1
            end
            if segment.idx == pulley.segment_idxs[2]
                eqs = [
                    eqs
                    l0[segment.idx] ~ get_sum_len(psys, pulley.idx) - pulley_len[pulley.idx]
                ]
                in_pulley += 1
            end
        end
        (in_pulley > 1) && error("Bridle segment number $(segment.idx) is part of
              $in_pulley pulleys, and should be part of either 0 or 1 pulleys.")

        if in_pulley == 0
            in_tether = 0
            for tether in tethers
                if segment.idx in tether.segment_idxs # each tether segment has to be part of exactly one tether
                    in_winch = 0
                    winch_idx = 0
                    for winch in winches
                        if tether.idx in winch.tether_idxs
                            winch_idx = winch.idx
                            in_winch += 1
                        end
                    end
                    (in_winch != 1) && error("Tether number $(tether.idx) is connected to
                          $(in_winch) winches, and should have 1 winch connected.")

                    eqs = [
                        eqs
                        l0[segment.idx] ~ tether_len[winch_idx] / length(tether.segment_idxs)
                    ]
                    in_tether += 1
                end
            end
            !(in_tether in [0,1]) && error("Segment number $(segment.idx) is part of 
                  $in_tether tethers, and should be part of exactly 0 or 1 tether.")
            if in_tether == 0
                eqs = [
                    eqs
                    l0[segment.idx] ~ get_l0(psys, segment.idx)
                ]
            end
        end

        eqs = [
            eqs
            # spring force equations
            segment_vec[:, segment.idx]  ~ pos[:, p2] - pos[:, p1]
            len[segment.idx]             ~ norm(segment_vec[:, segment.idx])
            unit_vec[:, segment.idx]   ~ segment_vec[:, segment.idx]/len[segment.idx]
            rel_vel[:, segment.idx]      ~ vel[:, p1] - vel[:, p2]
            spring_vel[segment.idx]      ~ rel_vel[:, segment.idx] ⋅ unit_vec[:, segment.idx]
            damping[segment.idx]         ~ get_axial_damping(psys, segment.idx) / len[segment.idx]
            stiffness[segment.idx]       ~ ifelse(len[segment.idx] > l0[segment.idx],
                                            get_axial_stiffness(psys, segment.idx) / len[segment.idx],
                                            get_compression_frac(psys, segment.idx) * get_axial_stiffness(psys, segment.idx) / len[segment.idx])
            spring_force[segment.idx] ~  (stiffness[segment.idx] * (len[segment.idx] - l0[segment.idx]) - 
                                     damping[segment.idx] * spring_vel[segment.idx])
            spring_force_vec[:, segment.idx]  ~ spring_force[segment.idx] * unit_vec[:, segment.idx]
            
            # drag force equations
            height[segment.idx]          ~ max(0.0, 0.5(pos[:, p1][3] + pos[:, p2][3]))
            segment_vel[:, segment.idx]  ~ 0.5(vel[:, p1] + vel[:, p2])
            segment_rho[segment.idx]     ~ calc_rho(s.am, height[segment.idx])
            wind_vel[:, segment.idx]     ~ AtmosphericModels.calc_wind_factor(s.am, 
                                           max(height[segment.idx], 1.0), s.set.profile_law) * wind_vec_gnd
            va[:, segment.idx]           ~ wind_vel[:, segment.idx] - segment_vel[:, segment.idx]
            area[segment.idx]            ~ len[segment.idx] * get_diameter(psys, segment.idx)
            app_perp_vel[:, segment.idx] ~ va[:, segment.idx] - 
                                           (va[:, segment.idx] ⋅ unit_vec[:, segment.idx]) * unit_vec[:, segment.idx]
            drag_force[:, segment.idx]   ~ (0.5 * segment_rho[segment.idx] * get_cd_tether(pset) * norm(va[:, segment.idx]) * area[segment.idx]) * app_perp_vel[:, segment.idx]
        ]
    end

    # ==================== PULLEYS ==================== #
    @variables begin
        pulley_len(t)[eachindex(pulleys)]
        pulley_vel(t)[eachindex(pulleys)]
        pulley_force(t)[eachindex(pulleys)]
        pulley_acc(t)[eachindex(pulleys)]
    end
    @parameters pulley_damp = 5.0
    for pulley in pulleys
        segment = segments[pulley.segment_idxs[1]]
        mass_per_meter = get_rho_tether(pset) * π * (get_diameter(psys, segment.idx)/2)^2
        mass = get_sum_len(psys, pulley.idx) * mass_per_meter
        eqs = [
            eqs
            pulley_force[pulley.idx]    ~ spring_force[pulley.segment_idxs[1]] - spring_force[pulley.segment_idxs[2]]
            pulley_acc[pulley.idx]      ~ pulley_force[pulley.idx] / mass
        ]
        if pulley.type == DYNAMIC
            eqs = [
                eqs
                D(pulley_len[pulley.idx])  ~ pulley_vel[pulley.idx]
                D(pulley_vel[pulley.idx]) ~ pulley_acc[pulley.idx] - pulley_damp * pulley_vel[pulley.idx]
            ]
            defaults = [
                defaults
                pulley_len[pulley.idx] => get_pulley_len(psys, pulley.idx)
                pulley_vel[pulley.idx] => get_pulley_vel(psys, pulley.idx)
            ]
        elseif pulley.type == QUASI_STATIC
            eqs = [
                eqs 
                pulley_vel[pulley.idx] ~ 0
                pulley_acc[pulley.idx] ~ 0
            ]
            guesses = [
                guesses
                pulley_len[pulley.idx] => get_l0(psys, pulley.segment_idxs[1])
            ]
        else
            error("Wrong pulley type")
        end
    end

    # ==================== WINCHES ==================== #
    @variables begin
        tether_vel(t)[eachindex(winches)]
        tether_acc(t)[eachindex(winches)]
        winch_force(t)[eachindex(winches)]
        winch_force_vec(t)[1:3, eachindex(winches)]
        brake(t)[eachindex(winches)]
    end
    for winch in winches
        F = zeros(Num, 3)
        for tether_idx in winch.tether_idxs
            winch_idx = tethers[tether_idx].winch_idx
            (winch_idx > length(points)) && 
                error("Point number $winch_idx does not exist.")
            F .+= point_force[:, winch_idx]
        end
        eqs = [
            eqs
            brake[winch.idx] ~ get_brake(psys, winch.idx)
            D(tether_len[winch.idx]) ~ ifelse(brake[winch.idx]==true,
                                              0, tether_vel[winch.idx])
            D(tether_vel[winch.idx]) ~ ifelse(brake[winch.idx]==true,
                                              0, tether_acc[winch.idx])

            tether_acc[winch.idx] ~ calc_moment_acc( # TODO: moment and speed control
                winch.model, tether_vel[winch.idx], 
                winch_force[winch.idx], 
                set_values[winch.idx]
            )
            winch_force_vec[:, winch.idx] ~ F
            winch_force[winch.idx] ~ norm(winch_force_vec[:, winch.idx])
        ]
        defaults = [
            defaults
            tether_len[winch.idx] => get_tether_len(psys, winch.idx)
            tether_vel[winch.idx] => get_tether_vel(psys, winch.idx)
        ]
    end

    # ==================== TETHERS ==================== #
    @variables begin
        stretched_len(t)[eachindex(tethers)]
        tether_spring_force(t)[eachindex(tethers)]
    end
    for tether in tethers
        slen = zero(Num)
        tforce = zero(Num)
        for segment_idx in tether.segment_idxs
            slen += len[segment_idx]
            tforce += spring_force[segment_idx]
        end
        tforce /= length(tether.segment_idxs)
        eqs = [
            eqs
            stretched_len[tether.idx] ~ slen
            tether_spring_force[tether.idx] ~ tforce
        ]
    end

    return eqs, defaults, guesses, tether_wing_force, tether_wing_moment
end

"""
    wing_eqs!(s, eqs, psys, pset, defaults; tether_wing_force, ...)

Generate the differential equations for the wing's rigid body dynamics.

This function builds the equations for:
- Quaternion kinematics for the wing's orientation.
- Euler's rotation equations for the angular acceleration, including gyroscopic effects.
- Newton's second law for the translational motion of the wing's center of mass.

# Arguments
- `s::SymbolicAWEModel`: The main model object.
- `eqs`, `psys`, `pset`, `defaults`: Accumulating vectors and symbolic parameters.
- `tether_wing_force`, `tether_wing_moment`: Aggregate forces/moments from tethers.
- `aero_force_b`, `aero_moment_b`: Aerodynamic forces/moments.
- `ω_b`, `α_b`, `R_b_w`, `wing_pos`, `wing_vel`, `wing_acc`: Symbolic state variables.

# Returns
- `(eqs, defaults)`: A tuple containing the updated equation and default value lists.
"""
function wing_eqs!(s, eqs, psys, pset, defaults; tether_wing_force, tether_wing_moment, aero_force_b, 
    aero_moment_b, ω_b, α_b, R_b_w, wing_pos, wing_vel, wing_acc, fix_wing
)
    wings = s.sys_struct.wings
    @variables begin
        # potential differential variables
        wing_acc_b(t)[eachindex(wings), 1:3]
        α_b_damped(t)[eachindex(wings), 1:3]
        ω_b_stable(t)[eachindex(wings), 1:3]

        # rotations and frames
        Q_b_w(t)[eachindex(wings), 1:4] # quaternion orientation of the wing body frame relative to the world frame
        Q_vel(t)[eachindex(wings), 1:4] # quaternion rate of change

        # rest: forces, moments, vectors and scalar values
        moment_b(t)[eachindex(wings), 1:3] # moment in principal frame
        wing_mass(t)[eachindex(wings)]
        fix_wing_sphere(t)[eachindex(wings)]
    end
    @parameters ω_damp = 150

    Ω(ω) = [0      -ω[1]  -ω[2]  -ω[3];
            ω[1]    0      ω[3]  -ω[2];
            ω[2]   -ω[3]   0      ω[1];
            ω[3]    ω[2]  -ω[1]   0]

    for wing in wings
        vsm_wing = wing.vsm_wing
        I_b = [vsm_wing.inertia_tensor[1,1], vsm_wing.inertia_tensor[2,2], vsm_wing.inertia_tensor[3,3]]
        axis = sym_normalize(wing_pos[wing.idx, :])
        axis_b = R_b_w[wing.idx, :, :]' * axis
        eqs = [
            eqs
            fix_wing_sphere[wing.idx] ~ get_fix_wing_sphere(psys, wing.idx)
            [D(Q_b_w[wing.idx, i]) ~ Q_vel[wing.idx, i] for i in 1:4]
            [Q_vel[wing.idx, i] ~ 0.5 * sum(Ω(ω_b_stable[wing.idx, :])[i, j] * Q_b_w[wing.idx, j] for j in 1:4) for i in 1:4]
            ω_b_stable[wing.idx, :] ~ ifelse.(fix_wing==true,
                zeros(3),
                ifelse.(fix_wing_sphere[wing.idx]==true,
                    ω_b[wing.idx, :] - ω_b[wing.idx, :] ⋅ axis_b * axis_b,
                    ω_b[wing.idx, :]
                )
            )
            D(ω_b[wing.idx, :]) ~ ifelse.(fix_wing==true,
                zeros(3),
                ifelse.(fix_wing_sphere[wing.idx]==true,
                    α_b_damped[wing.idx, :] - α_b_damped[wing.idx, :] ⋅ axis_b * axis_b,
                    α_b_damped[wing.idx, :]
                )
            )
            α_b_damped[wing.idx, :] ~ [α_b[wing.idx, 1], α_b[wing.idx, 2] - ω_damp*ω_b[wing.idx, 2], α_b[wing.idx, 3]]
    
            [R_b_w[wing.idx, :, i] ~ quaternion_to_rotation_matrix(Q_b_w[wing.idx, :])[:, i] for i in 1:3]
            
            α_b[wing.idx, 1] ~ (moment_b[wing.idx, 1] + (I_b[2] - I_b[3]) * ω_b[wing.idx, 2] * ω_b[wing.idx, 3]) / I_b[1]
            α_b[wing.idx, 2] ~ (moment_b[wing.idx, 2] + (I_b[3] - I_b[1]) * ω_b[wing.idx, 3] * ω_b[wing.idx, 1]) / I_b[2]
            α_b[wing.idx, 3] ~ (moment_b[wing.idx, 3] + (I_b[1] - I_b[2]) * ω_b[wing.idx, 1] * ω_b[wing.idx, 2]) / I_b[3]
            moment_b[wing.idx, :] ~ aero_moment_b[wing.idx, :] + R_b_w[wing.idx, :, :]' * tether_wing_moment[wing.idx, :]
            
            D(wing_pos[wing.idx, :]) ~ ifelse.(fix_wing==true,
                zeros(3),
                ifelse.(fix_wing_sphere[wing.idx]==true,
                    wing_vel[wing.idx, :] ⋅ axis * axis,
                    wing_vel[wing.idx, :]
                )
            )
            D(wing_vel[wing.idx, :]) ~ ifelse.(fix_wing==true,
                zeros(3),
                ifelse.(fix_wing_sphere[wing.idx]==true,
                    wing_acc[wing.idx, :] ⋅ axis * axis,
                    wing_acc[wing.idx, :]
                )
            )
            wing_mass[wing.idx] ~ get_set_mass(pset)
            wing_acc[wing.idx, :] ~ (tether_wing_force[wing.idx, :] + R_b_w[wing.idx, :, :] * aero_force_b[wing.idx, :]) / wing_mass[wing.idx]
        ]
        defaults = [
            defaults
            [Q_b_w[wing.idx, i] => get_Q_b_w(psys, wing.idx)[i] for i in 1:4]
            [ω_b[wing.idx, i] => get_ω_b(psys, wing.idx)[i] for i in 1:3]
            [wing_pos[wing.idx, i] => get_wing_pos_w(psys, wing.idx)[i] for i in 1:3]
            [wing_vel[wing.idx, i] => get_wing_vel_w(psys, wing.idx)[i] for i in 1:3]
        ]
    end
    
    return eqs, defaults
end

"""
    rotate_v_around_k(v, k, θ)

Rotate vector `v` around axis `k` by angle `θ` using Rodrigues' rotation formula.
"""
function rotate_v_around_k(v, k, θ)
    k = sym_normalize(k)
    v_rot = v * cos(θ) + (k × v) * sin(θ)  + k * (k ⋅ v) * (1 - cos(θ))
    return v_rot
end

"""
    calc_R_v_w(wing_pos, e_x)

Calculate the rotation matrix from the view frame (`_v`) to the world frame (`_w`).

The view frame is defined with its z-axis pointing from the origin to the wing,
and its x-axis aligned with the wing's x-axis projected onto the view plane.
"""
function calc_R_v_w(wing_pos, e_x)
    z = sym_normalize(wing_pos)
    y = sym_normalize(z × e_x)
    x = y × z
    return [x y z]
end

"""
    calc_R_t_w(elevation, azimuth)

Calculate the rotation matrix from the tether frame (`_t`) to the world frame (`_w`).

The tether frame is a spherical coordinate system defined by elevation and azimuth angles.
"""
function calc_R_t_w(elevation, azimuth)
    x = rotate_around_z(rotate_around_y([0, 0, -1], -elevation), azimuth)
    z = rotate_around_z(rotate_around_y([1, 0, 0], -elevation), azimuth)
    y = z × x
    return [x y z]
end

"""
    calc_heading(R_t_w, R_v_w)

Calculate the heading angle [rad] of the wing.

Heading is defined as the rotation angle between the tether frame and the view frame.
"""
function calc_heading(R_t_w::AbstractMatrix, R_v_w::AbstractMatrix)
    heading_vec = R_t_w' * R_v_w[:, 1]
    heading = atan(heading_vec[2], heading_vec[1])
    return heading
end

"""
    scalar_eqs!(s, eqs, psys, pset; R_b_w, wind_vec_gnd, ...)

Generate equations for derived scalar kinematic quantities.

This function calculates variables that are useful for control and analysis but are not
fundamental states of the system, such as elevation, azimuth, heading, and course angles,
along with their time derivatives.

# Arguments
- `s::SymbolicAWEModel`: The main model object.
- `eqs`, `psys`, `pset`: Accumulating vectors and symbolic parameters.
- `R_b_w`, `wind_vec_gnd`, etc.: Symbolic variables for the system's state.

# Returns
- `eqs`: The updated list of system equations.
"""
function scalar_eqs!(s, eqs, psys, pset; R_b_w, wind_vec_gnd, va_wing_b, wing_pos, wing_vel, wing_acc, twist_angle, twist_ω, ω_b, α_b)
    @unpack wings = s.sys_struct
    wind_scale_gnd = get_v_wind(pset)
    @variables begin
        e_x(t)[eachindex(wings), 1:3]
        e_y(t)[eachindex(wings), 1:3]
        e_z(t)[eachindex(wings), 1:3]
        wind_vel_wing(t)[eachindex(wings), 1:3]
        wind_disturb(t)[eachindex(wings), 1:3]
        va_wing(t)[eachindex(wings), 1:3]
        upwind_dir(t)
    end
    eqs = [
        eqs
        upwind_dir ~ deg2rad(get_upwind_dir(pset))
        wind_vec_gnd ~ max(wind_scale_gnd, 1e-6) * rotate_around_z([0, -1, 0], -upwind_dir)
    ]
    for wing in wings
        eqs = [
            eqs
            e_x[wing.idx, :]       ~ R_b_w[wing.idx, :,1]
            e_y[wing.idx, :]       ~ R_b_w[wing.idx, :,2]
            e_z[wing.idx, :]       ~ R_b_w[wing.idx, :,3]
            wind_vel_wing[wing.idx, :] ~ AtmosphericModels.calc_wind_factor(s.am, 
                    max(wing_pos[wing.idx, 3], 1.0), s.set.profile_law) * wind_vec_gnd
            wind_disturb[wing.idx, :] ~ get_wind_disturb(psys, wing.idx)
            va_wing[wing.idx, :] ~ wind_vel_wing[wing.idx, :] - wing_vel[wing.idx, :] + 
                                       wind_disturb[wing.idx, :]
            va_wing_b[wing.idx, :] ~ R_b_w[wing.idx, :, :]' * va_wing[wing.idx, :]
        ]
    end
    @variables begin
        heading(t)[eachindex(wings)]
        turn_rate(t)[eachindex(wings), 1:3]
        turn_acc(t)[eachindex(wings), 1:3]
        azimuth(t)[eachindex(wings)]
        azimuth_vel(t)[eachindex(wings)]
        azimuth_acc(t)[eachindex(wings)]
        elevation(t)[eachindex(wings)]
        elevation_vel(t)[eachindex(wings)]
        elevation_acc(t)[eachindex(wings)]
        course(t)[eachindex(wings)]
        x_acc(t)[eachindex(wings)]
        y_acc(t)[eachindex(wings)]
        sphere_pos(t)[eachindex(wings), 1:2, 1:2]
        sphere_vel(t)[eachindex(wings), 1:2, 1:2]
        sphere_acc(t)[eachindex(wings), 1:2, 1:2]
        angle_of_attack(t)[eachindex(wings)]
        R_v_w(t)[eachindex(wings), 1:3, 1:3]
        R_t_w(t)[eachindex(wings), 1:3, 1:3]
        distance(t)[eachindex(wings)]
        distance_vel(t)[eachindex(wings)]
        distance_acc(t)[eachindex(wings)]
    end

    for wing in wings
        x, y, z = wing_pos[wing.idx, :]
        x´, y´, z´ = wing_vel[wing.idx, :]
        x´´, y´´, z´´ = wing_acc[wing.idx, :]

        half_len = wing.group_idxs[1] + length(wing.group_idxs)÷2 - 1

        eqs = [
            eqs
            vec(R_v_w[wing.idx, :, :])    .~ vec(calc_R_v_w(wing_pos[wing.idx, :], e_x[wing.idx, :]))
            vec(R_t_w[wing.idx, :, :])    .~ vec(calc_R_t_w(elevation[wing.idx], azimuth[wing.idx]))
            heading[wing.idx]         ~ calc_heading(R_t_w[wing.idx, :, :], R_v_w[wing.idx, :, :])
            turn_rate[wing.idx, :]      ~ R_v_w[wing.idx, :, :]' * (R_b_w[wing.idx, :, :] * ω_b[wing.idx, :]) # Project angular velocity onto view frame
            turn_acc[wing.idx, :]       ~ R_v_w[wing.idx, :, :]' * (R_b_w[wing.idx, :, :] * α_b[wing.idx, :])
            distance[wing.idx]        ~ norm(wing_pos[wing.idx, :])
            distance_vel[wing.idx]    ~ wing_vel[wing.idx, :] ⋅ R_v_w[wing.idx, :, 3]
            distance_acc[wing.idx]    ~ wing_acc[wing.idx, :] ⋅ R_v_w[wing.idx, :, 3]

            elevation[wing.idx]           ~ KiteUtils.calc_elevation(wing_pos[wing.idx, :])
            elevation_vel[wing.idx]       ~ (x*z´ - z*x´) / 
                                            (x^2 + z^2)
            elevation_acc[wing.idx]       ~ ((x^2 + z^2)*(x*z´´ - z*x´´) + 2(z*x´ - x*z´)*(x*x´ + z*z´))/(x^2 + z^2)^2
            azimuth[wing.idx]             ~ -KiteUtils.azimuth_east(wing_pos[wing.idx, :])
            azimuth_vel[wing.idx]         ~ (-y*x´ + x*y´) / 
                                            (x^2 + y^2)
            azimuth_acc[wing.idx]         ~ ((x^2 + y^2)*(-y*x´´ + x*y´´) + 2(y*x´ - x*y´)*(x*x´ + y*y´))/(x^2 + y^2)^2
            course[wing.idx]              ~ atan(-azimuth_vel[wing.idx], elevation_vel[wing.idx])
            x_acc[wing.idx]               ~ wing_acc ⋅ e_x
            y_acc[wing.idx]               ~ wing_acc ⋅ e_y

            angle_of_attack[wing.idx]     ~ calc_angle_of_attack(va_wing_b[wing.idx, :]) + 
                                            0.5twist_angle[half_len] + 0.5twist_angle[half_len+1]
        ]
    end
    return eqs
end

"""
    Base.getindex(x::ModelingToolkit.Symbolics.SymArray, idxs::Vector{Int16})

Extend `getindex` to allow indexing a symbolic array with a vector of indices.
"""
function Base.getindex(x::ModelingToolkit.Symbolics.SymArray, idxs::Vector{Int16})
    Num[Base.getindex(x, idx) for idx in idxs]
end

"""
    linear_vsm_eqs!(s, eqs, guesses, psys; aero_force_b, ...)

Generate linearized aerodynamic equations using the Vortex Step Method (VSM).

This function approximates the complex, nonlinear aerodynamic forces and moments
by using a first-order Taylor expansion around the current operating point. The
Jacobian of the aerodynamic forces with respect to the state variables (`va_wing_b`,
`ω_b`, `twist_angle`) is pre-calculated and provided as a parameter (`vsm_jac`).

# Arguments
- `s::SymbolicAWEModel`: The main model object.
- `eqs`, `guesses`, `psys`: Accumulating vectors and symbolic parameters.
- `aero_force_b`, `aero_moment_b`, etc.: Symbolic variables for aerodynamic and state quantities.

# Returns
- `(eqs, guesses)`: A tuple containing the updated equation and guess lists.
"""
function linear_vsm_eqs!(s, eqs, guesses, psys; aero_force_b, aero_moment_b, group_aero_moment, twist_angle, va_wing_b, ω_b)
    @unpack groups, wings = s.sys_struct
    if length(wings) == 0
        return eqs, guesses
    end

    ny = 3+length(wings[1].group_idxs)+3
    nx = 3+3+length(wings[1].group_idxs)

    @variables begin
        y(t)[eachindex(wings), 1:ny]
        dy(t)[eachindex(wings), 1:ny]
        last_y(t)[eachindex(wings), 1:ny] 
        last_x(t)[eachindex(wings), 1:nx]
        vsm_jac(t)[eachindex(wings), 1:nx, 1:ny]
    end

    for wing in wings
        eqs = [
            eqs
            [last_y[wing.idx, iy] ~ get_vsm_y(psys, wing.idx, iy) for iy in 1:ny]
            [last_x[wing.idx, ix] ~ get_vsm_x(psys, wing.idx, ix) for ix in 1:nx]
            [vsm_jac[wing.idx, ix, iy] ~ get_vsm_jac(psys, wing.idx, ix, iy) for ix in 1:nx for iy in 1:ny]
            y[wing.idx, :] ~ [va_wing_b[wing.idx, :]; ω_b[wing.idx, :]; twist_angle[wing.group_idxs]]
            dy[wing.idx, :] ~ y[wing.idx, :] - last_y[wing.idx, :]
            [aero_force_b[wing.idx, :]; aero_moment_b[wing.idx, :]; group_aero_moment[wing.group_idxs]] ~ 
                last_x[wing.idx, :] + vsm_jac[wing.idx, :, :] * dy[wing.idx, :]
        ]
    
        if s.set.quasi_static
            guesses = [guesses; [y[wing.idx, iy] => get_vsm_y(psys, wing.idx, iy) for iy in 1:ny]]
        end
    end
    return eqs, guesses
end

"""
    create_sys!(s::SymbolicAWEModel, system::SystemStructure; prn=true)

Create the full `ModelingToolkit.System` for the AWE model.

This is the main top-level function that orchestrates the generation of the entire
set of differential-algebraic equations (DAEs). It calls specialized sub-functions
to build the equations for each part of the system (forces, wing dynamics, scalar
kinematics, linearized aerodynamics) and assembles them into a single `System` object.

# Arguments
- `s::SymbolicAWEModel`: The main model object to be populated.
- `system::SystemStructure`: The physical structure definition.
- `prn::Bool=true`: If true, print progress information.

# Returns
- `set_values`: The symbolic variable representing the control inputs.
"""
function create_sys!(s::SymbolicAWEModel, system::SystemStructure; prn=true)
    eqs = []
    defaults = Pair{Num, Any}[]
    guesses = Pair{Num, Any}[]

    @unpack wings, groups, winches = system

    @parameters begin
        psys::SystemStructure = system
        pset::Settings = s.set
        fix_wing = false
    end
    @variables begin
        # potential differential variables
        set_values(t)[eachindex(winches)] = zeros(length(winches))
        wing_pos(t)[eachindex(wings), 1:3] # xyz pos of wing in world frame
        wing_vel(t)[eachindex(wings), 1:3]
        wing_acc(t)[eachindex(wings), 1:3]
        ω_b(t)[eachindex(wings), 1:3] # turn rate in principal frame
        α_b(t)[eachindex(wings), 1:3]

        # rotations and frames
        R_b_w(t)[eachindex(wings), 1:3, 1:3] # rotation of the wing body frame relative to the world frame

        # rest: forces, moments, vectors and scalar values
        aero_force_b(t)[eachindex(wings), 1:3]
        aero_moment_b(t)[eachindex(wings), 1:3]
        twist_angle(t)[eachindex(groups)]
        twist_ω(t)[eachindex(groups)]
        group_aero_moment(t)[eachindex(groups)]
        wind_vec_gnd(t)[1:3]
        va_wing_b(t)[eachindex(wings), 1:3]
    end

    eqs, defaults, guesses, tether_wing_force, tether_wing_moment = 
        force_eqs!(s, system, psys, pset, eqs, defaults, guesses; 
            R_b_w, wing_pos, wing_vel, wind_vec_gnd, group_aero_moment, 
            twist_angle, twist_ω, set_values, fix_wing)
    eqs, guesses = linear_vsm_eqs!(s, eqs, guesses, psys; aero_force_b, 
            aero_moment_b, group_aero_moment, twist_angle, va_wing_b, ω_b)
    eqs, defaults = wing_eqs!(s, eqs, psys, pset, defaults; 
            tether_wing_force, tether_wing_moment, aero_force_b, aero_moment_b, 
            ω_b, α_b, R_b_w, wing_pos, wing_vel, wing_acc, fix_wing)
    eqs = scalar_eqs!(s, eqs, psys, pset; 
            R_b_w, wind_vec_gnd, va_wing_b, wing_pos, wing_vel, wing_acc, 
            twist_angle, twist_ω, ω_b, α_b)
    
    # te_I = (1/3 * (get_set_mass(pset)/8) * te_len^2)
    # # -damping / I * ω = α_damping
    # # solve for c: (c * (k*m/s^2) / (k*m^2)) * (m/s)=m/s^2 in wolframalpha
    # # damping should be in N*m*s
    # rot_damping = 0.1s.damping * te_len

    # eqs = [
    #     eqs
    #     trailing_edge_α[1] ~ (force[:, s.i_A]) ⋅ e_te_A * te_len / te_I - (rot_damping[1] / te_I) * trailing_edge_ω[1] # TODO: add trailing edge
    #     trailing_edge_α[2] ~ (force[:, s.i_B]) ⋅ e_te_B * te_len / te_I - (rot_damping[2] / te_I) * trailing_edge_ω[2]
    # ]
    
    eqs = Symbolics.scalarize.(reduce(vcat, Symbolics.scalarize.(eqs)))

    # discrete_events = [
    #     true => [
    #         [Q_b_w[i] ~ normalize(Q_b_w)[i] for i in 1:4]
    #         [twist_angle[i] ~ clamp(twist_angle[i], -π/2, π/2) for i in eachindex(s.point_groups)]
    #         ]
    #     ]

    # @named sys = System(eqs, t; discrete_events)
    time = @elapsed @named sys = System(eqs, t)
    prn && println("\tCreated System in $time seconds.")

    defaults = [
        defaults
        [set_values[winch.idx] => get_set_value(psys, winch.idx) for winch in winches]
    ]

    s.defaults = defaults
    s.guesses = guesses
    s.full_sys = sys
    return set_values
end
