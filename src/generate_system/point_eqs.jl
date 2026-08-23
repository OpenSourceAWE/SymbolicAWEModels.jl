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
    beam_hermite_ride_eqs(point, force_on_point, s, params; kwargs...)

Kinematics of a `point` that rides its `TimoshenkoJoint`'s corotational cubic-Hermite
centerline at `beam_frac` (transverse deflection from the two end slopes, plus a
frame-carried `beam_offset_b`), and the load it feeds to the two end bodies split by
axial fraction. Mutates `body_force`/`body_moment` in place; returns the point's
pos/vel/acc equations. Shared by beam-anchored bridle points and beam-anchored
wing-node aero receivers, so both track the deformed beam identically.
"""
function beam_hermite_ride_eqs(point, force_on_point, s, params;
                               pos, vel, acc, body_pos_w, body_R_b_to_w, body_com_w,
                               body_com_vel, body_ω_b, body_force, body_moment)
    joint = s.sys_struct.timoshenko_joints[point.joint_idx]
    a = joint.body_a_idx
    b = joint.body_b_idx
    R_a = collect(body_R_b_to_w[:, :, a])
    R_b = collect(body_R_b_to_w[:, :, b])
    ex = beam_hermite_ride_expressions(joint, params, point.idx;
        pos_a = collect(body_pos_w[:, a]), R_a, com_a = collect(body_com_w[:, a]),
        com_vel_a = collect(body_com_vel[:, a]), omega_a_w = R_a * collect(body_ω_b[:, a]),
        pos_b = collect(body_pos_w[:, b]), R_b, com_b = collect(body_com_w[:, b]),
        com_vel_b = collect(body_com_vel[:, b]), omega_b_w = R_b * collect(body_ω_b[:, b]))
    pos_point = ex.pos_point
    force_a = (1 - ex.sfrac) .* force_on_point
    force_b = ex.sfrac .* force_on_point
    body_force[:, a] .+= force_a
    body_force[:, b] .+= force_b
    body_moment[:, a] .+= (pos_point .- collect(body_com_w[:, a])) × force_a
    body_moment[:, b] .+= (pos_point .- collect(body_com_w[:, b])) × force_b
    return [
        pos[:, point.idx] ~ pos_point
        vel[:, point.idx] ~ ex.vel_point
        acc[:, point.idx] ~ zeros(3)
    ]
end

"""
    point_eqs!(s, eqs, defaults, points, segments, twist_surfaces, wings, params, initial;
               R_b_to_w, wing_vel, wind_vec_gnd, twist_angle,
               pos, vel, acc, point_force, point_mass, spring_force_vec, drag_force, l0,
               spring_sum_force, point_aero_drag, total_drag,
               disturb_force, tether_r, chord_b, fixed_pos, normal, pos_b,
               fix_point_sphere, fix_static,
               va_point_b, va_point_w, wind_at_point, height,
               aero_force_point_b,
               twist_surface_y_airf)

Generate equations for all point types (STATIC, DYNAMIC, BODY_STATIC).

Each point's net force is the shared [`point_net_force`](@ref): the structural load
gathered from its incident segments (their spring force with the endpoint sign and
half their drag), plus per-node aero, its own drag and gravity. A point that rides a
rigid body has zero gravitational mass here, since that mass is carried at the body's
COM. Free particles integrate through [`confined_derivatives`](@ref), which applies
`fix_static` and `fix_sphere`; a rigid wing node is placed instead by
[`twist_deformed_offset`](@ref) from its body's COM.

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
- `body_force`, `body_moment`: Mutable arrays to accumulate wing-node loads onto bodies.

# Returns
- Tuple `(eqs, defaults)` with updated equation vectors.
  Note: `body_force` and `body_moment` are modified in-place.
"""
function point_eqs!(s, eqs, defaults, points, segments, twist_surfaces, wings, params, initial;
                    R_b_to_w, com_w,
                    wing_vel, wind_vec_gnd, twist_angle,
                    pos, vel, acc, point_force, point_mass, spring_force_vec, drag_force, l0,
                    spring_sum_force, point_aero_drag, total_drag,
                    disturb_force, tether_r, chord_b, fixed_pos, normal, pos_b,
                    fix_point_sphere, fix_static,
                    va_point_b, va_point_w, wind_at_point, height,
                    aero_force_point_b,
                    twist_surface_y_airf,
                    body_force, body_moment, body_pos_w, body_com_w, body_R_b_to_w,
                    body_com_vel, body_ω_b)

    wind_gnd = collect(wind_vec_gnd)
    for point in points
        F::Vector{Num} = zeros(Num, 3)
        seg_drag::Vector{Num} = zeros(Num, 3)
        mass = params.points[point.idx].extra_mass
        for segment in segments
            if point.idx in segment.point_idxs
                inverted = segment.point_idxs[2] == point.idx
                if inverted
                    F .-= spring_force_vec[:, segment.idx]
                else
                    F .+= spring_force_vec[:, segment.idx]
                end
                mass += segment_half_mass(l0[segment.idx],
                    params.segments[segment.idx].diameter,
                    params.segments[segment.idx].density)
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

        drag_coeff = params.points[point.idx].drag_coeff
        area = params.points[point.idx].area
        wind_source = point_wind_source(params, point.idx, wind_gnd)
        drag_rhs = point_drag_force(collect(va_point_w[:, point.idx]),
            air_density(s.am, height[point.idx]), drag_coeff, area)
        va_point_b_rhs = isnothing(wing_idx_for_transform) ? zeros(3) :
            R_b_to_w[:, :, wing_idx_for_transform]' * va_point_w[:, point.idx]
        eqs = [
            eqs
            height[point.idx] ~ max(0.0, pos[3, point.idx])
            wind_at_point[:, point.idx] ~ wind_source(pos[3, point.idx])
            va_point_w[:, point.idx] ~
                wind_at_point[:, point.idx] - vel[:, point.idx]
            va_point_b[:, point.idx] ~ va_point_b_rhs
            point_aero_drag[:, point.idx] ~ drag_rhs
        ]

        # Total drag: point aero drag + share of segment drag
        eqs = [
            eqs
            total_drag[:, point.idx] ~
                point_aero_drag[:, point.idx] + seg_drag
        ]

        # A rigid wing's own VSM wrench replaces the per-node aero force.
        rigid_wing_node = point.is_wing_node &&
            s.sys_struct.bodies[point.wing_idx].dynamics_type == RIGID_DYNAMICS
        aero_force_w = (!rigid_wing_node && point.is_wing_node &&
                        is_wing(s.sys_struct.bodies[point.wing_idx])) ?
            collect(R_b_to_w[:, :, point.wing_idx] *
                    aero_force_point_b[:, point.idx]) : zeros(Num, 3)
        rides_own_wing = point.body_idx > 0 &&
            wing_frame_member(point, point.body_idx)
        carries_gravity = !(rigid_wing_node || rides_own_wing)
        eqs = [
            eqs
            point_force[:, point.idx] ~ point_net_force(s,
                collect(pos[:, point.idx]), collect(vel[:, point.idx]),
                collect(spring_sum_force[:, point.idx]) .+ aero_force_w .+
                    collect(disturb_force[:, point.idx]),
                carries_gravity ? mass : 0.0, drag_coeff, area, wind_source,
                carries_gravity ? params.set.g_earth : 0.0)
        ]

        # EXCEPTION to the anchor rule: a wing node on a RIGID_DYNAMICS wing is
        # placed by a twist-deformed offset in the wing body frame — the wing is
        # itself the rigid body, so there is no node body or beam to ride, and the
        # section-twist DOF deforms the offset. Every other point follows the
        # anchor rule below (joint → beam, body → rigid ride, STATIC → fixed,
        # else → free particle).
        if rigid_wing_node
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
            if in_group
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
                    fixed_pos[:, point.idx] ~
                        params.twist_surfaces[twist_surface.idx].le_pos
                    chord_b[:, point.idx] ~
                        params.points[point.idx].pos_undeformed_b .-
                            fixed_pos[:, point.idx]
                    normal[:, point.idx] ~ chord_b[:, point.idx] ×
                        twist_surface_y_airf[:, twist_surface.idx]
                ]
            end
            surface_idx = in_group ? twist_surface.idx : 0
            eqs = [eqs; pos_b[:, point.idx] ~ twist_deformed_offset(
                params, point.idx, surface_idx,
                in_group ? twist_angle[twist_surface.idx] : 0.0)]
            eqs = [
                eqs
                tether_r[:, point.idx] ~
                    pos[:, point.idx] - com_w[:, point.wing_idx]
            ]
            point_moment = tether_r[:, point.idx] × point_force[:, point.idx]
            if in_group
                point_moment = ifelse.(
                    params.bodies[point.wing_idx].group_points_moment == true,
                    point_moment, zeros(3))
            end
            body_force[:, point.wing_idx] .+= point_force[:, point.idx]
            body_moment[:, point.wing_idx] .+= point_moment
            eqs = [
                eqs
                pos[:, point.idx] ~
                    com_w[:, point.wing_idx] +
                    R_b_to_w[:, :, point.wing_idx] * pos_b[:, point.idx]
                vel[:, point.idx] ~ zeros(3)
                acc[:, point.idx] ~ zeros(3)
            ]
            continue
        end

        # Placement by anchor — the single rule for every non-rigid-wing point.
        force_on_point = collect(point_force[:, point.idx])
        if point.joint_idx > 0
            eqs = [eqs; beam_hermite_ride_eqs(point, force_on_point, s, params;
                pos, vel, acc, body_pos_w, body_R_b_to_w, body_com_w,
                body_com_vel, body_ω_b, body_force, body_moment)]
        elseif point.body_idx > 0
            eqs = [eqs; body_ride_eqs(point, point.body_idx, force_on_point,
                params; pos, vel, acc, body_pos_w, body_R_b_to_w, body_com_w,
                body_com_vel, body_ω_b, body_force, body_moment)]
        elseif point.type == STATIC
            eqs = [
                eqs
                pos[:, point.idx] ~ params.points[point.idx].pos_w
                vel[:, point.idx] ~ zeros(3)
                acc[:, point.idx] ~ zeros(3)
            ]
        else
            # Free particle: integrated position/velocity (DYNAMIC point or an
            # unanchored surface node).
            pars = point_particle_params(params, point.idx)
            wing_idx_damp = length(wings) > 0 ? point.wing_idx : 0
            vel_diff_w = length(wings) > 0 ?
                vel[:, point.idx] - wing_vel[:, point.wing_idx] : nothing
            damp_accel = point_damping_accel(
                point, params, R_b_to_w, wing_idx_damp, vel[:, point.idx], vel_diff_w)
            velocity, acceleration = confined_derivatives(
                pos[:, point.idx], vel[:, point.idx], collect(acc[:, point.idx]),
                (; fix_sphere = fix_point_sphere[point.idx],
                   fix_static = fix_static[point.idx]))
            eqs = [
                eqs
                fix_point_sphere[point.idx] ~ pars.fix_sphere
                fix_static[point.idx] ~ pars.fix_static
                D(pos[:, point.idx]) ~ velocity
                D(vel[:, point.idx]) ~ acceleration
                acc[:, point.idx] ~ point_force[:, point.idx] ./ mass - damp_accel
            ]
            defaults = [
                defaults
                bind_initial!(initial.points[point.idx].pos_w, collect(pos[:, point.idx]))
                bind_initial!(initial.points[point.idx].vel_w, collect(vel[:, point.idx]))
            ]
        end
    end

    return eqs, defaults
end
