# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Point dynamics equation generation

"""
    point_damping_accel(point, params, R_b_to_w, wing_idx, vel_w, vel_diff_w)

Per-mass damping acceleration for a DYNAMIC point. Each frame's term is built
only when its coefficient is set; a `nothing` coefficient keeps that term out of
the equation entirely (no zero-valued parameter to prune). The body-frame term
also needs a wing frame — pass `vel_diff_w = nothing` (point velocity relative to
its wing) to skip it when no wing is available.
"""
function point_damping_accel(point, params, R_b_to_w, wing_idx, vel_w, vel_diff_w)
    accel = zeros(Num, 3)
    if !isnothing(point.body_frame_damping) && !isnothing(vel_diff_w)
        R = R_b_to_w[:, :, wing_idx]
        coeff = params.points[point.idx].body_frame_damping
        accel = accel + R * (coeff .* (R' * vel_diff_w))
    end
    if !isnothing(point.world_frame_damping)
        coeff = params.points[point.idx].world_frame_damping
        accel = accel + coeff .* vel_w
    end
    return accel
end

"""
    body_ride_eqs(point, body, force_on_body, params;
                  pos, vel, acc, body_pos_w, body_R_b_to_w, body_com_w,
                  body_com_vel, body_ω_b, body_force, body_moment)

Kinematic pose equations for a point rigidly anchored to `body` at
`anchor_b`, and the load it feeds back to that body. The point tracks the
body's rigid motion (`vel = com_vel + ω×arm`), and its `force_on_body` is
applied at the anchor with the moment about the body COM. Shared by the
`BODY_STATIC` single-body ride and PARTICLE-wing beam-node coupling; mutates
`body_force`/`body_moment` in place and returns the pos/vel/acc equations.
"""
function body_ride_eqs(point, body, force_on_body, params;
                       pos, vel, acc, body_pos_w, body_R_b_to_w, body_com_w,
                       body_com_vel, body_ω_b, body_force, body_moment)
    anchor = collect(params.points[point.idx].anchor_b)
    R_body = collect(body_R_b_to_w[:, :, body])
    anchor_w = collect(body_pos_w[:, body]) .+ R_body * anchor
    ω_w = R_body * collect(body_ω_b[:, body])
    arm = anchor_w .- collect(body_com_w[:, body])
    body_force[:, body] .+= force_on_body
    body_moment[:, body] .+= arm × force_on_body
    return [
        pos[:, point.idx] ~ anchor_w
        vel[:, point.idx] ~ collect(body_com_vel[:, body]) .+ (ω_w × arm)
        acc[:, point.idx] ~ zeros(3)
    ]
end

"""
    point_eqs!(s, eqs, defaults, points, segments, twist_surfaces, wings, params, initial;
               R_b_to_w, wing_vel, wind_vec_gnd, twist_angle,
               pos, vel, acc, point_force, point_mass, spring_force_vec, drag_force, l0,
               spring_sum_force, point_drag_force, total_drag,
               disturb_force, tether_r, chord_b, fixed_pos, normal, pos_b,
               fix_point_sphere, fix_static,
               va_point_b, va_point_w, wind_at_point, height,
               aero_force_point_b,
               twist_surface_y_airf)

Generate equations for all point types (STATIC, DYNAMIC, WING).

# Arguments
- `s::SymbolicAWEModel`: The main model object (for atmospheric model).
- `eqs`, `defaults`: Accumulating vectors for the MTK system.
- `points`, `segments`, `twist_surfaces`, `wings`: System components.
- `R_b_to_w`: Symbolic rotation matrix (body to world).
- `wing_vel`: Symbolic wing center of mass velocity.
- `wind_vec_gnd`: Symbolic ground-level wind vector.
- `twist_angle`: Symbolic twist_surface twist angle.
- `pos`, `vel`, `acc`: Pre-declared point state variables.
- `point_force`, `point_mass`: Pre-declared point force and mass variables.
- `spring_force_vec`, `drag_force`, `l0`: Pre-declared segment force variables.
- `spring_sum_force`: Pre-declared accumulated spring/drag forces variable.
- Other variables: Various point-specific symbolic variables.
- `body_force`, `body_moment`: Mutable arrays to accumulate WING-point loads onto bodies.

# Returns
- Tuple `(eqs, defaults)` with updated equation vectors.
  Note: `body_force` and `body_moment` are modified in-place.
"""
function point_eqs!(s, eqs, defaults, points, segments, twist_surfaces, wings, params, initial;
                    R_b_to_w, com_w,
                    wing_vel, wind_vec_gnd, twist_angle,
                    pos, vel, acc, point_force, point_mass, spring_force_vec, drag_force, l0,
                    spring_sum_force, point_drag_force, total_drag,
                    disturb_force, tether_r, chord_b, fixed_pos, normal, pos_b,
                    fix_point_sphere, fix_static,
                    va_point_b, va_point_w, wind_at_point, height,
                    aero_force_point_b,
                    twist_surface_y_airf,
                    body_force, body_moment, body_pos_w, body_com_w, body_R_b_to_w,
                    body_com_vel, body_ω_b)

    wind_factor = param_computed!(params.reg, :wind_factor, WindFactorReader())
    for point in points
        F::Vector{Num} = zeros(Num, 3)
        seg_drag::Vector{Num} = zeros(Num, 3)
        mass = params.points[point.idx].extra_mass
        for segment in segments
            if point.idx in segment.point_idxs
                mass_per_meter =
                    params.segments[segment.idx].density * π *
                    (params.segments[segment.idx].diameter / 2)^2
                inverted = segment.point_idxs[2] == point.idx
                if inverted
                    F .-= spring_force_vec[:, segment.idx]
                else
                    F .+= spring_force_vec[:, segment.idx]
                end
                mass += mass_per_meter * l0[segment.idx] / 2
                half_seg_drag = 0.5 * drag_force[:, segment.idx]
                F .+= half_seg_drag
                seg_drag .+= half_seg_drag
            end
        end

        eqs = [
            eqs
            spring_sum_force[:, point.idx] ~ F
            point_mass[point.idx] ~ mass
            disturb_force[:, point.idx] ~ params.points[point.idx].ext_force_w
        ]

        # Apparent velocity for ALL points (PARTICLE_DYNAMICS wings need body frame).
        wing_idx_for_transform = if point.is_wing_node
            point.wing_idx
        elseif length(wings) > 0
            # Use first wing for non-wing points
            Int64(1)
        else
            nothing
        end

        drag_rhs = 0.5 * calc_rho(s.am, height[point.idx]) *
            params.points[point.idx].drag_coeff *
            smooth_norm(va_point_w[:, point.idx]) *
            params.points[point.idx].area *
            va_point_w[:, point.idx]
        va_point_b_rhs = isnothing(wing_idx_for_transform) ? zeros(3) :
            R_b_to_w[:, :, wing_idx_for_transform]' * va_point_w[:, point.idx]
        eqs = [
            eqs
            height[point.idx] ~ max(0.0, pos[3, point.idx])
            wind_at_point[:, point.idx] ~
                wind_factor(pos[3, point.idx]) * wind_vec_gnd
            va_point_w[:, point.idx] ~
                wind_at_point[:, point.idx] - vel[:, point.idx]
            va_point_b[:, point.idx] ~ va_point_b_rhs
            point_drag_force[:, point.idx] ~ drag_rhs
        ]

        # Total drag: point aero drag + share of segment drag
        eqs = [
            eqs
            total_drag[:, point.idx] ~
                point_drag_force[:, point.idx] + seg_drag
        ]

        if point.type == BODY_STATIC && !point.is_wing_node
            eqs = [
                eqs
                point_force[:, point.idx] ~
                    spring_sum_force[:, point.idx] +
                    Num[0, 0, -params.set.g_earth * mass] +
                    disturb_force[:, point.idx] + point_drag_force[:, point.idx]
            ]
            force_on_point = collect(point_force[:, point.idx])
            if point.joint_idx > 0
                # Rides the joint's corotational-Hermite centerline; load split by axial fraction to both bodies.
                joint = s.sys_struct.timoshenko_joints[point.joint_idx]
                a = joint.body_a_idx
                b = joint.body_b_idx
                R_a = collect(body_R_b_to_w[:, :, a])
                R_b = collect(body_R_b_to_w[:, :, b])
                jp = params.timoshenko_joints[joint.idx]
                x_a = collect(body_pos_w[:, a]) .+ R_a * collect(jp.anchor_a_b)
                x_b = collect(body_pos_w[:, b]) .+ R_b * collect(jp.anchor_b_b)
                e1, e2, e3, beam_len = timoshenko_element_frame(x_a, x_b, R_a)
                element_frame = [e1[1] e2[1] e3[1];
                                 e1[2] e2[2] e3[2];
                                 e1[3] e2[3] e3[3]]
                Da = (element_frame' * R_a) * collect(jp.R_a_rel0)'
                Db = (element_frame' * R_b) * collect(jp.R_b_rel0)'
                θ_a = [0.5 * (Da[3, 2] - Da[2, 3]),
                       0.5 * (Da[1, 3] - Da[3, 1]),
                       0.5 * (Da[2, 1] - Da[1, 2])]
                θ_b = [0.5 * (Db[3, 2] - Db[2, 3]),
                       0.5 * (Db[1, 3] - Db[3, 1]),
                       0.5 * (Db[2, 1] - Db[1, 2])]
                sfrac = params.points[point.idx].beam_frac
                # Cubic-Hermite transverse deflection from the two end slopes.
                N2 = beam_len * (sfrac - 2sfrac^2 + sfrac^3)
                N4 = beam_len * (-sfrac^2 + sfrac^3)
                v_defl = N2 * θ_a[3] + N4 * θ_b[3]
                w_defl = -(N2 * θ_a[2] + N4 * θ_b[2])
                x_center = x_a .+ (sfrac * beam_len) .* e1 .+
                    v_defl .* e2 .+ w_defl .* e3
                offset = collect(params.points[point.idx].beam_offset_b)
                pos_point = x_center .+ element_frame * offset
                ω_a_w = R_a * collect(body_ω_b[:, a])
                ω_b_w = R_b * collect(body_ω_b[:, b])
                vel_a = collect(body_com_vel[:, a]) .+
                    (ω_a_w × (pos_point .- collect(body_com_w[:, a])))
                vel_b = collect(body_com_vel[:, b]) .+
                    (ω_b_w × (pos_point .- collect(body_com_w[:, b])))
                eqs = [
                    eqs
                    pos[:, point.idx] ~ pos_point
                    vel[:, point.idx] ~ (1 - sfrac) .* vel_a .+ sfrac .* vel_b
                    acc[:, point.idx] ~ zeros(3)
                ]
                force_a = (1 - sfrac) .* force_on_point
                force_b = sfrac .* force_on_point
                body_force[:, a] .+= force_a
                body_force[:, b] .+= force_b
                body_moment[:, a] .+=
                    (pos_point .- collect(body_com_w[:, a])) × force_a
                body_moment[:, b] .+=
                    (pos_point .- collect(body_com_w[:, b])) × force_b
            else
                # Rides a Body kinematically; feeds its force and COM moment to the body.
                eqs = [eqs; body_ride_eqs(point, point.body_idx, force_on_point,
                    params; pos, vel, acc, body_pos_w, body_R_b_to_w, body_com_w,
                    body_com_vel, body_ω_b, body_force, body_moment)]
            end
        elseif point.is_wing_node
            # The wing is a body (looked up in the full bodies, incl. AeroNone wings).
            wing = s.sys_struct.bodies[point.wing_idx]

            if wing.dynamics_type == PARTICLE_DYNAMICS
                # AeroNone particle wings have no per-point aero array, so zero force.
                aero_force_w = is_wing(wing) ?
                    R_b_to_w[:, :, wing.idx] * aero_force_point_b[:, point.idx] :
                    zeros(Num, 3)

                eqs = [
                    eqs
                    point_force[:, point.idx] ~
                        spring_sum_force[:, point.idx] + aero_force_w + Num[0, 0, -params.set.g_earth * mass] + disturb_force[:, point.idx] + point_drag_force[:, point.idx]
                ]

                if point.body_idx > 0
                    # Rides a Body (beam-node coupling): kinematic pose, loads to body.
                    eqs = [eqs; body_ride_eqs(point, point.body_idx,
                        collect(point_force[:, point.idx]), params;
                        pos, vel, acc, body_pos_w, body_R_b_to_w, body_com_w,
                        body_com_vel, body_ω_b, body_force, body_moment)]
                    continue
                end

                damp_accel = point_damping_accel(
                    point, params, R_b_to_w, wing.idx, vel[:, point.idx],
                    vel[:, point.idx] - wing_vel[:, point.wing_idx])

                # DYNAMIC point equations
                axis = smooth_normalize(pos[:, point.idx])
                eqs = [
                    eqs
                    fix_point_sphere[point.idx] ~ params.points[point.idx].fix_sphere
                    fix_static[point.idx] ~ params.points[point.idx].fix_static
                    D(pos[:, point.idx]) ~ ifelse.(
                        fix_static[point.idx] == true,
                        zeros(3),
                        ifelse.(fix_point_sphere[point.idx]==true,
                                vel[:, point.idx] ⋅ axis * axis,
                                vel[:, point.idx]
                        )
                    )
                    D(vel[:, point.idx]) ~ ifelse.(
                        fix_static[point.idx] == true,
                        zeros(3),
                        ifelse.(fix_point_sphere[point.idx]==true,
                                acc[:, point.idx] ⋅ axis * axis,
                                acc[:, point.idx]
                        )
                    )
                    acc[:, point.idx] ~ point_force[:, point.idx] ./ mass - damp_accel
                ]
                defaults = [
                    defaults
                    bind_initial!(initial.points[point.idx].pos_w, collect(pos[:, point.idx]))
                    bind_initial!(initial.points[point.idx].vel_w, collect(vel[:, point.idx]))
                ]

            elseif wing.dynamics_type == RIGID_DYNAMICS
                # RIGID_DYNAMICS point feeds force to body accumulator; gravity at COM.
                eqs = [
                    eqs
                    point_force[:, point.idx] ~
                        spring_sum_force[:, point.idx] + disturb_force[:, point.idx] + point_drag_force[:, point.idx]
                ]

                found = 0
                twist_surface = nothing
                for twist_surface_ in twist_surfaces
                    if point.idx in twist_surface_.point_idxs
                        twist_surface = twist_surface_
                        found += 1
                    end
                end
                in_group = found == 1
                !(found in [0, 1]) && error(
                    "Kite point number $(point.idx) is part of $found twist_surfaces, " *
                    "and should be part of exactly 0 or 1 twist_surfaces.",
                )

                if found == 1
                    found = 0
                    for wing_ in s.sys_struct.bodies
                        if twist_surface.idx in wing_.twist_surface_idxs
                            found += 1
                        end
                    end
                    !(found == 1) && error(
                        "Kite twist_surface number $(twist_surface.idx) is part of $found bodies, " *
                        "and should be part of exactly 1 body.",
                    )

                    eqs = [
                        eqs
                        fixed_pos[:, point.idx] ~ params.twist_surfaces[twist_surface.idx].le_pos
                        chord_b[:, point.idx] ~
                            params.points[point.idx].pos_b .- fixed_pos[:, point.idx]
                        normal[:, point.idx] ~ chord_b[:, point.idx] × twist_surface_y_airf[:, twist_surface.idx]
                        pos_b[:, point.idx] ~
                            fixed_pos[:, point.idx] .+
                            cos(twist_angle[twist_surface.idx]) * chord_b[:, point.idx] -
                            sin(twist_angle[twist_surface.idx]) * normal[:, point.idx]
                    ]
                elseif found == 0
                    eqs = [eqs; pos_b[:, point.idx] ~ params.points[point.idx].pos_b]
                end
                # Moment arm about COM (world frame)
                eqs = [
                    eqs
                    tether_r[:, point.idx] ~
                        pos[:, point.idx] -
                        com_w[:, point.wing_idx]
                ]
                # group_points_moment may exclude in-group points from the wing moment.
                point_moment = tether_r[:, point.idx] ×
                    point_force[:, point.idx]
                if in_group
                    point_moment = ifelse.(
                        params.bodies[point.wing_idx].group_points_moment == true,
                        point_moment, zeros(3))
                end
                # Feed the body load accumulator (read by body_eqs!).
                body_force[:, point.wing_idx] .+= point_force[:, point.idx]
                body_moment[:, point.wing_idx] .+= point_moment

                # pos_b is the offset from COM in the body frame.
                eqs = [
                    eqs
                    pos[:, point.idx] ~
                        com_w[:, point.wing_idx] +
                        R_b_to_w[:, :, point.wing_idx] *
                        pos_b[:, point.idx]
                    vel[:, point.idx] ~ zeros(3)
                    acc[:, point.idx] ~ zeros(3)
                ]
            else
                error("Unsupported dynamics_type $(wing.dynamics_type) " *
                      "for WING point $(point.idx)")
            end
        elseif point.type == STATIC
            # Define point_force for STATIC points
            eqs = [
                eqs
                point_force[:, point.idx] ~
                    spring_sum_force[:, point.idx] + Num[0, 0, -params.set.g_earth * mass] + disturb_force[:, point.idx] + point_drag_force[:, point.idx]
                pos[:, point.idx] ~ params.points[point.idx].pos_w
                vel[:, point.idx] ~ zeros(3)
                acc[:, point.idx] ~ zeros(3)
            ]
        elseif point.type == DYNAMIC
            # Define point_force for DYNAMIC points
            eqs = [
                eqs
                point_force[:, point.idx] ~
                    spring_sum_force[:, point.idx] + Num[0, 0, -params.set.g_earth * mass] + disturb_force[:, point.idx] + point_drag_force[:, point.idx]
            ]

            wing_idx = length(wings) > 0 ? point.wing_idx : 0
            vel_diff_w = length(wings) > 0 ?
                vel[:, point.idx] - wing_vel[:, point.wing_idx] : nothing
            damp_accel = point_damping_accel(
                point, params, R_b_to_w, wing_idx, vel[:, point.idx], vel_diff_w)

            axis = smooth_normalize(pos[:, point.idx])
            eqs = [
                eqs
                fix_point_sphere[point.idx] ~ params.points[point.idx].fix_sphere
                fix_static[point.idx] ~ params.points[point.idx].fix_static
                D(pos[:, point.idx]) ~ ifelse.(
                    fix_static[point.idx] == true,
                    zeros(3),
                    ifelse.(fix_point_sphere[point.idx]==true,
                            vel[:, point.idx] ⋅ axis * axis,
                            vel[:, point.idx]
                    )
                )
                D(vel[:, point.idx]) ~ ifelse.(
                    fix_static[point.idx] == true,
                    zeros(3),
                    ifelse.(fix_point_sphere[point.idx]==true,
                            acc[:, point.idx] ⋅ axis * axis,
                            acc[:, point.idx]
                    )
                )
                acc[:, point.idx]    ~ point_force[:, point.idx] ./ mass - damp_accel
            ]
            defaults = [
                defaults
                bind_initial!(initial.points[point.idx].pos_w, collect(pos[:, point.idx]))
                bind_initial!(initial.points[point.idx].vel_w, collect(vel[:, point.idx]))
            ]
        else
            error("Unknown point type: $(typeof(point))")
        end
    end

    return eqs, defaults
end
