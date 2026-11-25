# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Point dynamics equation generation

"""
    point_damping_accel(sam, i, params, R_b_to_w, wing_vel, vel_w)

Per-mass damping acceleration for DYNAMIC point `i` of `sam`: a world-frame term against
its world velocity `vel_w`, plus, for a point that belongs to a wing, a body-frame term
in the frame of that wing's body against the point's velocity relative to it.
"""
function point_damping_accel(sam, i, params, R_b_to_w, wing_vel, vel_w)
    accel = collect(params.points[i].world_frame_damping .* vel_w)
    wing_idx = point_wing(sam, i)
    isnothing(wing_idx) && return accel
    body_idx = wing_body(sam, wing_idx)
    R = R_b_to_w[:, :, body_idx]
    coeff = params.points[i].body_frame_damping
    vel_diff_w = vel_w - wing_vel[:, body_idx]
    return accel + R * (coeff .* (R' * vel_diff_w))
end

"""
    body_ride_eqs(i, body, force_on_body, params;
                  pos, vel, acc, body_pos_w, body_R_b_to_w, body_com_w,
                  body_com_vel, body_ω_b, body_force, body_moment)

Kinematic pose equations for point `i` rigidly anchored to `body` at
`anchor_b`, and the load it feeds back to that body. The point tracks the
body's rigid motion (`vel = com_vel + ω×arm`), and its `force_on_body` is
applied at the anchor with the moment about the body COM. Shared by the
`BODY_STATIC` single-body ride and PARTICLE-wing beam-node coupling; mutates
`body_force`/`body_moment` in place and returns the pos/vel/acc equations.
"""
function body_ride_eqs(i, body, force_on_body, params;
                       pos, vel, acc, body_pos_w, body_R_b_to_w, body_com_w,
                       body_com_vel, body_ω_b, body_force, body_moment)
    anchor = collect(params.points[i].anchor_b)
    R_body = collect(body_R_b_to_w[:, :, body])
    anchor_w = collect(body_pos_w[:, body]) .+ R_body * anchor
    ω_w = R_body * collect(body_ω_b[:, body])
    arm = anchor_w .- collect(body_com_w[:, body])
    body_force[:, body] .+= force_on_body
    body_moment[:, body] .+= arm × force_on_body
    return [
        pos[:, i] ~ anchor_w
        vel[:, i] ~ collect(body_com_vel[:, body]) .+ (ω_w × arm)
        acc[:, i] ~ zeros(3)
    ]
end

"""
    beam_hermite_ride_eqs(sam, i, force_on_point, params; kwargs...)

Kinematics of point `i` of `sam` that rides its Timoshenko tube's corotational cubic-Hermite
centerline at `beam_frac` (transverse deflection from the two end slopes, plus a
frame-carried `beam_offset_b`), and the load it feeds to the two end bodies split by
axial fraction. Mutates `body_force`/`body_moment` in place; returns the point's
pos/vel/acc equations. Shared by beam-anchored bridle points and beam-anchored
wing-node aero receivers, so both track the deformed beam identically.
"""
function beam_hermite_ride_eqs(sam, i, force_on_point, params;
                               pos, vel, acc, body_pos_w, body_R_b_to_w, body_com_w,
                               body_com_vel, body_ω_b, body_force, body_moment)
    tube_idx = sam.definition.points[i].tube
    a, b = sam.definition.tubes[tube_idx].bodies
    R_a = collect(body_R_b_to_w[:, :, a])
    R_b = collect(body_R_b_to_w[:, :, b])
    ex = beam_hermite_ride_expressions(tube_idx, params, i;
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
        pos[:, i] ~ pos_point
        vel[:, i] ~ ex.vel_point
        acc[:, i] ~ zeros(3)
    ]
end

"""
    point_eqs!(s, eqs, defaults, params, initial;
               R_b_to_w, com_w, wing_vel, wind_vec_gnd, twist_angle,
               pos, vel, acc, point_force, point_mass, spring_force_vec, drag_force, l0,
               spring_sum_force, point_aero_drag, total_drag,
               disturb_force, tether_r, chord_b, fixed_pos, normal, pos_b,
               fix_point_sphere, fix_static,
               va_point_b, va_point_w, wind_at_point, height,
               aero_force_point_b, station_y_airf,
               body_force, body_moment, body_pos_w, body_com_w, body_R_b_to_w,
               body_com_vel, body_ω_b)

Generate equations for all point types (STATIC, DYNAMIC, BODY_STATIC) of the run `s`.

Each point's net force is the shared [`point_net_force`](@ref): the structural load
gathered from its incident segments (their spring force with the endpoint sign and
half their drag), plus per-node aero, its own drag and gravity. A point a rigid body
carries ([`carrier_body_idx`](@ref)) has zero gravitational mass here, since its mass
is part of the body's `total_mass` and weighs at the body's COM. Free particles
integrate through [`confined_derivatives`](@ref), which applies `fix_static` and
`fix_sphere`; a rigid wing node is placed instead by [`twist_deformed_offset`](@ref)
from its body's COM.

# Arguments
- `s::SymbolicAWEModel`: The run.
- `eqs`, `defaults`: Accumulating vectors for the MTK system.
- `R_b_to_w`, `com_w`, `wing_vel`: Symbolic body rotation, centre of mass and velocity,
  by body.
- `wind_vec_gnd`: Symbolic ground-level wind vector.
- `twist_angle`: Symbolic station twist angle.
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
function point_eqs!(s, eqs, defaults, params, initial;
                    R_b_to_w, com_w,
                    wing_vel, wind_vec_gnd, twist_angle,
                    pos, vel, acc, point_force, point_mass, spring_force_vec, drag_force, l0,
                    spring_sum_force, point_aero_drag, total_drag,
                    disturb_force, tether_r, chord_b, fixed_pos, normal, pos_b,
                    fix_point_sphere, fix_static,
                    va_point_b, va_point_w, wind_at_point, height,
                    aero_force_point_b,
                    station_y_airf,
                    body_force, body_moment, body_pos_w, body_com_w, body_R_b_to_w,
                    body_com_vel, body_ω_b)
    sys = s.definition
    wind_gnd = collect(wind_vec_gnd)
    for (i, point) in enumerate(sys.points)
        F::Vector{Num} = zeros(Num, 3)
        seg_drag::Vector{Num} = zeros(Num, 3)
        mass = params.points[i].extra_mass
        for k in s.topology.point_segments[i]
            if sys.segments[k].points[2] == i
                F .-= spring_force_vec[:, k]
            else
                F .+= spring_force_vec[:, k]
            end
            mass += segment_half_mass(l0[k], params.segments[k].diameter,
                                      params.segments[k].density)
            half_seg_drag = 0.5 * drag_force[:, k]
            F .+= half_seg_drag
            seg_drag .+= half_seg_drag
        end

        eqs = [
            eqs
            spring_sum_force[:, i] ~ F
            point_mass[i] ~ mass
            disturb_force[:, i] ~ params.points[i].ext_force_w
        ]

        wing_idx = point_wing(s, i)
        body_idx = isnothing(wing_idx) ? 0 : wing_body(s, wing_idx)
        drag_coeff = params.points[i].drag_coeff
        area = params.points[i].area
        wind_source = point_wind_source(params, i, wind_gnd)
        drag_rhs = point_drag_force(collect(va_point_w[:, i]),
            air_density(s.am, height[i]), drag_coeff, area)
        va_point_b_rhs = body_idx != 0 ?
            R_b_to_w[:, :, body_idx]' * va_point_w[:, i] : zeros(3)
        eqs = [
            eqs
            height[i] ~ max(0.0, pos[3, i])
            wind_at_point[:, i] ~ wind_source(pos[3, i])
            va_point_w[:, i] ~ wind_at_point[:, i] - vel[:, i]
            va_point_b[:, i] ~ va_point_b_rhs
            point_aero_drag[:, i] ~ drag_rhs
        ]

        # Total drag: point aero drag + share of segment drag
        eqs = [eqs; total_drag[:, i] ~ point_aero_drag[:, i] + seg_drag]

        # A rigid wing's own VSM wrench replaces the per-node aero force.
        is_wing_node = wing_node(sys, i) && body_idx != 0
        rigid_wing_node = is_wing_node &&
            body_dynamics_type(s, body_idx) == RIGID_DYNAMICS
        aero_force_w = (!rigid_wing_node && is_wing_node && has_aero(s.aero[wing_idx])) ?
            collect(R_b_to_w[:, :, body_idx] * aero_force_point_b[:, i]) : zeros(Num, 3)
        carries_gravity = carrier_body_idx(s, i) == 0
        eqs = [
            eqs
            point_force[:, i] ~ point_net_force(s,
                collect(pos[:, i]), collect(vel[:, i]),
                collect(spring_sum_force[:, i]) .+ aero_force_w .+
                    collect(disturb_force[:, i]),
                carries_gravity ? mass : 0.0, drag_coeff, area, wind_source,
                carries_gravity ? params.set.g_earth : 0.0)
        ]

        # EXCEPTION to the anchor rule: a wing node on a RIGID_DYNAMICS wing is
        # placed by a twist-deformed offset in the wing body frame — the wing is
        # itself the rigid body, so there is no node body or beam to ride, and the
        # section-twist DOF deforms the offset. Every other point follows the
        # anchor rule below (tube → beam, body → rigid ride, STATIC → fixed,
        # else → free particle).
        if rigid_wing_node
            stations = [station_idx
                        for (station_idx, station) in enumerate(sys.spanwise_stations)
                        if i in station.points]
            length(stations) in (0, 1) || error(
                "Kite point number $i is part of $(length(stations)) stations, " *
                "and should be part of exactly 0 or 1 stations.")
            station_idx = isempty(stations) ? 0 : only(stations)
            if station_idx != 0
                owners = count(wing -> station_idx in wing_station_idxs(s, wing),
                               eachindex(sys.wings))
                owners == 1 || error("Kite station number $station_idx is part of " *
                    "$owners wings, and should be part of exactly 1 wing.")
                eqs = [
                    eqs
                    fixed_pos[:, i] ~ params.spanwise_stations[station_idx].le_pos
                    chord_b[:, i] ~ params.points[i].pos_undeformed_b .- fixed_pos[:, i]
                    normal[:, i] ~ chord_b[:, i] × station_y_airf[:, station_idx]
                ]
            end
            eqs = [eqs; pos_b[:, i] ~ twist_deformed_offset(params, i, station_idx,
                station_idx != 0 ? twist_angle[station_idx] : 0.0)]
            eqs = [eqs; tether_r[:, i] ~ pos[:, i] - com_w[:, body_idx]]
            point_moment = tether_r[:, i] × point_force[:, i]
            if station_idx != 0
                point_moment = ifelse.(
                    params.bodies[body_idx].group_points_moment == true,
                    point_moment, zeros(3))
            end
            body_force[:, body_idx] .+= point_force[:, i]
            body_moment[:, body_idx] .+= point_moment
            eqs = [
                eqs
                pos[:, i] ~ com_w[:, body_idx] + R_b_to_w[:, :, body_idx] * pos_b[:, i]
                vel[:, i] ~ zeros(3)
                acc[:, i] ~ zeros(3)
            ]
            continue
        end

        # Placement by anchor — the single rule for every non-rigid-wing point.
        force_on_point = collect(point_force[:, i])
        if !isnothing(point.tube)
            eqs = [eqs; beam_hermite_ride_eqs(s, i, force_on_point, params;
                pos, vel, acc, body_pos_w, body_R_b_to_w, body_com_w,
                body_com_vel, body_ω_b, body_force, body_moment)]
        elseif !isnothing(point.body)
            eqs = [eqs; body_ride_eqs(i, point.body, force_on_point,
                params; pos, vel, acc, body_pos_w, body_R_b_to_w, body_com_w,
                body_com_vel, body_ω_b, body_force, body_moment)]
        elseif point_type(point) == STATIC
            eqs = [
                eqs
                pos[:, i] ~ params.points[i].pos_w
                vel[:, i] ~ zeros(3)
                acc[:, i] ~ zeros(3)
            ]
        else
            # Free particle: integrated position/velocity (DYNAMIC point or an
            # unanchored surface node).
            pars = point_particle_params(params, i)
            damp_accel = point_damping_accel(s, i, params, R_b_to_w, wing_vel, vel[:, i])
            velocity, acceleration = confined_derivatives(
                pos[:, i], vel[:, i], collect(acc[:, i]),
                (; fix_sphere = fix_point_sphere[i], fix_static = fix_static[i]))
            eqs = [
                eqs
                fix_point_sphere[i] ~ pars.fix_sphere
                fix_static[i] ~ pars.fix_static
                D(pos[:, i]) ~ velocity
                D(vel[:, i]) ~ acceleration
                acc[:, i] ~ point_force[:, i] ./ (mass + pars.apparent_mass) - damp_accel
            ]
            defaults = [
                defaults
                bind_initial!(initial.points[i].pos_w, collect(pos[:, i]))
                bind_initial!(initial.points[i].vel_w, collect(vel[:, i]))
            ]
        end
    end

    return eqs, defaults
end
