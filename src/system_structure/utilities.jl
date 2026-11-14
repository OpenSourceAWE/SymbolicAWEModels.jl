# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# ==================== SHARED HELPERS ==================== #

"""
    segment_initial_length(segment, points)

Compute segment length from endpoint `pos_ENU` positions.
"""
function segment_initial_length(segment, points)
    point1 = points[segment.point_idxs[1]]
    point2 = points[segment.point_idxs[2]]
    return norm(point1.pos_ENU - point2.pos_ENU)
end

"""
    segment_world_length(segment, points)

Compute segment length from endpoint `pos_w` positions.
"""
function segment_world_length(segment, points)
    point1 = points[segment.point_idxs[1]]
    point2 = points[segment.point_idxs[2]]
    return norm(point1.pos_w - point2.pos_w)
end




# ==================== VALIDATION ==================== #

"""
    validate_sys_struct(sys_struct::SystemStructure)

Error on values that break a run and warn on suspicious ones: NaN or non-positive point and
wing mass, wing and winch inertia, NaN positions, rigid wings without stations, non-positive
`drum_radius` or `gear_ratio`, segment `diameter` outside (0, 1) m, non-positive `l0` or
stiffness, negative damping, zero pulley `sum_len`, unequal station `moment_frac`.
"""
function validate_sys_struct(sys_struct)
    (; points, stations, segments, pulleys, wings, winches) = sys_struct

    # ==================== POINT VALIDATIONS ==================== #
    for point in points
        # Check for NaN extra_mass
        if isnan(point.extra_mass)
            error("Point $(point.name) has NaN extra_mass")
        end

        # Warn about negative extra_mass (physically nonsensical but still works)
        if point.extra_mass < 0
            @warn "Point $(point.name) has negative extra_mass $(point.extra_mass) kg. " *
                  "This is physically nonsensical."
        end

        # Check mass before NaN position: NaN pos is often caused by zero mass.
        if point.type == DYNAMIC && point.total_mass <= 0
            error("Point $(point.name) has non-positive total_mass ($(point.total_mass)). " *
                  "This will cause division by zero in acceleration calculations.")
        end

        # Check for NaN position (often a symptom of zero mass)
        if any(isnan.(point.pos_w))
            error("Point $(point.name) has NaN position: pos_w = $(point.pos_w)")
        end
    end

    # ==================== WING VALIDATIONS ==================== #
    for wing in wings
        # Check mass/inertia before NaN position: NaN pos is often caused by zero mass.
        if wing.total_mass <= 0
            error("Wing $(wing.name) has non-positive total_mass ($(wing.total_mass)). " *
                  "This will cause division by zero in acceleration calculations.")
        end

        if wing.dynamics_type == RIGID_DYNAMICS
            I_b = wing.inertia_principal

            # Check for zero or suspiciously small inertia (before NaN checks)
            for i in 1:3
                if I_b[i] ≈ 0.0
                    error("Wing $(wing.name) has zero inertia component " *
                          "I_b[$i] = $(I_b[i]). " *
                          "All principal inertia components must be non-zero.")
                elseif I_b[i] < 1e-6
                    @warn "Wing $(wing.name) has very small inertia component " *
                          "I_b[$i] = $(I_b[i]) kg⋅m², may cause numerical issues"
                end
            end

            # Check for NaN inertia
            if any(isnan.(I_b))
                error("Wing $(wing.name) has NaN inertia: I_b = $I_b")
            end

            # AeroNone does not couple to sections, so missing stations is fine.
            if isempty(wing.station_idxs) &&
               couples_to_sections(wing.aero)
                @warn "Wing $(wing.name) (RIGID_DYNAMICS)" *
                    " has no stations"
            end
        end

        # Check for NaN position (often a symptom of zero mass/inertia)
        if any(isnan.(wing.pos_w))
            error("Wing $(wing.name) has NaN position: pos_w = $(wing.pos_w)")
        end
    end

    # ==================== WINCH VALIDATIONS ==================== #
    for winch in winches
        # Check for NaN inertia
        if isnan(winch.inertia_total)
            error("Winch $(winch.name) has NaN inertia_total")
        end

        # Check for zero or negative inertia
        if winch.inertia_total ≈ 0.0
            error("Winch $(winch.name) has zero inertia_total. " *
                  "All winches must have non-zero inertia.")
        elseif winch.inertia_total < 0
            error("Winch $(winch.name) has negative inertia_total " *
                  "$(winch.inertia_total) kg⋅m². Inertia must be positive.")
        elseif winch.inertia_total < 1e-6
            @warn "Winch $(winch.name) has very small inertia_total " *
                  "$(winch.inertia_total) kg⋅m², may cause numerical issues"
        end

        # Check for non-positive drum_radius
        if winch.drum_radius <= 0
            error("Winch $(winch.name) has non-positive drum_radius ($(winch.drum_radius)). " *
                  "This will cause division by zero in torque conversions.")
        end

        # Check for non-positive gear_ratio
        if winch.gear_ratio <= 0
            error("Winch $(winch.name) has non-positive gear_ratio ($(winch.gear_ratio)). " *
                  "This will cause division by zero in speed/torque calculations.")
        end
    end

    # ==================== SEGMENT VALIDATIONS ==================== #
    for segment in segments
        # Wing structural segments don't use diameter (stiffness explicit, drag from VSM)
        wing_structural = all(points[i].is_wing_node for i in segment.point_idxs)
        if !wing_structural && !(0 < segment.diameter < 1)
            @warn "Segment $(segment.name) has unusual diameter " *
                  "$(segment.diameter) m (expected range: 0 to 1 m)"
        end

        # Rest length must be positive
        if segment.l0 <= 0
            error("Segment $(segment.name) has non-positive rest length " *
                  "l0 = $(segment.l0) m. This will cause division by zero.")
        end

        # Warn about zero or negative stiffness/damping (callable laws self-validate).
        if segment.unit_stiffness isa Real
            if segment.unit_stiffness ≈ 0.0
                @warn "Segment $(segment.name) has zero stiffness"
            elseif segment.unit_stiffness < 0
                @warn "Segment $(segment.name) has negative stiffness " *
                      "$(segment.unit_stiffness) N"
            end
        end

        if segment.unit_damping < 0
            @warn "Segment $(segment.name) has negative damping " *
                  "$(segment.unit_damping) N⋅s"
        end
    end

    # ==================== PULLEY VALIDATIONS ==================== #
    for pulley in pulleys
        if pulley.sum_len ≈ 0
            error("Pulley $(pulley.name) has zero total length constraint " *
                  "(sum_len = $(pulley.sum_len) m). " *
                  "Pulley constraints must have non-zero total length.")
        end
    end

    # ==================== STATION VALIDATIONS ==================== #
    if length(stations) > 0
        first_moment_frac = stations[1].moment_frac
        for station in stations
            if !(station.moment_frac ≈ first_moment_frac)
                error("Station $(station.name) has moment_frac = " *
                      "$(station.moment_frac), but all stations must have the " *
                      "same moment_frac (first station has $(first_moment_frac))")
            end
        end
    end

    return nothing
end

"""
    set_unstretched_length!(sys_struct::SystemState, tether, len)

Set `tether`'s unstretched length [m] and share it equally over its segments'
`l0`. Point positions, body poses, tube rest geometry and station flap
references are left as they are.
"""
function set_unstretched_length!(sys_struct::SystemState, tether, len)
    len > 0 || error("Tether $(tether.name): unstretched length $len m " *
        "must be positive")
    isempty(tether.segment_idxs) && error("Tether $(tether.name): " *
        "has no segments to carry an unstretched length")
    tether.len = len
    l0 = len / length(tether.segment_idxs)
    for seg_idx in tether.segment_idxs
        sys_struct.segments[seg_idx].l0 = l0
    end
    return nothing
end

"""
    seed_per_point_wind!(state::SystemState)

Give every point and every wing the ground wind `set.wind_vec` as its own wind, so a
[`PerPointWind`](@ref) model that is never written to flies in the uniform wind the
settings describe. Called from [`init_wind!`](@ref); from then on these winds belong to the
caller, who writes them between steps.
"""
function seed_per_point_wind!(state)
    for point in state.points
        point.wind_vec .= state.set.wind_vec
    end
    for wing in state.wings
        wing.wind_vec .= state.set.wind_vec
    end
    return nothing
end

# ==================== PLACEMENT ==================== #

"""
    reset_to_initial_pose!(state::SystemState; reset_vel=true)

Put every point and body back at its initial pose (`pos_ENU`, `origin_ENU`,
`Q_KA_to_ENU`), zero the twist of every non-`STATIC` station and set every winch's
reel-out speed to `init_vel`. `reset_vel` zeroes point and body velocities as well.
"""
function reset_to_initial_pose!(state::SystemState; reset_vel::Bool=true)
    for winch in state.winches
        winch.vel = winch.init_vel
    end
    for point in state.points
        point.pos_w .= point.pos_ENU
        reset_vel && (point.vel_w .= 0.0)
    end
    for body in state.bodies
        body.pos_w .= body.origin_ENU
        body.Q_b_to_w .= body.Q_KA_to_ENU
        reset_vel && (body.vel_w .= 0.0; body.ω_b .= 0.0)
    end
    for station in state.stations
        station.type == STATIC && continue
        station.twist = 0.0
        station.twist_ω = 0.0
    end
    return nothing
end

"""
    remake_wing_aero!(state::SystemState, set::Settings)

Rebuild every wing's aero engine from `set` and the structure's `vsm_set`, e.g.
after editing `aero_geometry.yaml`. A no-op for aero modes without an engine.
"""
function remake_wing_aero!(state, set::Settings)
    for wing in state.wings
        remake_aero!(wing.aero, wing, set, state.vsm_set, state)
    end
    return nothing
end

"""
    init_wind!(state::SystemState, set::Settings)

Set each wing's wind from the ground wind `set.wind_vec` at the wing's height, or
seed every point and wing with it under [`PerPointWind`](@ref), and initialise each
wing's aerodynamic operating point from that wind.
"""
function init_wind!(state, set::Settings)
    (; wings) = state
    if per_point_wind(state)
        seed_per_point_wind!(state)
    else
        wind_factor = WindFactor(state.am, set.profile_law)
        for wing in wings
            wing.wind_vec .= wind_factor(wing.pos_w[3]) * set.wind_vec
        end
    end
    for wing in wings
        wing.dynamics_type == PARTICLE_DYNAMICS && continue
        R_b_to_w = wing.R_b_to_w::Matrix{SimFloat}
        init_aero_state!(wing.aero, wing, R_b_to_w' *
            (per_point_wind(state) ? wing.wind_vec : set.wind_vec))
    end
    return nothing
end

"""
    init_sys_struct!(state::SystemState, set::Settings; remake_vsm=false, prn=true)

Bring a placed structure to the start of a run: derive its mass properties and rest
geometry ([`init_derived_properties!`](@ref)), rebuild every wing's aero engine when
`remake_vsm` ([`remake_wing_aero!`](@ref)), and set the wind and aero operating point
([`init_wind!`](@ref)). [`init!`](@ref) runs it.
"""
function init_sys_struct!(state, set::Settings;
                          remake_vsm::Bool=false, prn::Bool=true)
    init_derived_properties!(state; prn)
    remake_vsm && remake_wing_aero!(state, set)
    init_wind!(state, set)
    return nothing
end

"""
    init_rest_geometry!(sys_struct::SystemStructure)

Capture every tube's rest configuration and every flap station's rest deflection
from the initial pose, so the structure is unstrained there.
"""
function init_rest_geometry!(sys_struct)
    init_tube_rest!.(sys_struct.tubes, Ref(sys_struct.bodies))
    init_station_flap!.(sys_struct.stations, Ref(sys_struct))
    return nothing
end

"""
    init_derived_properties!(sys_struct::SystemStructure; prn=true)

Derive from the initial pose every body's mass properties
([`update_mass_properties!`](@ref)) and the rest geometry of its tubes and flap
stations ([`init_rest_geometry!`](@ref)).
"""
function init_derived_properties!(sys_struct; prn::Bool=true)
    update_mass_properties!(sys_struct; prn)
    init_rest_geometry!(sys_struct)
    return nothing
end

"""
    start_run!(state::SystemState; reset_vel=true)

Put `state` at the start of a run from the pose its structure is placed in: every
point and body at its initial pose ([`reset_to_initial_pose!`](@ref)), the twist of a
`STATIC` station, masses, segment, tether and pulley lengths from its segments' rest
lengths, and
points and bodies moving as the rates of its transforms say
([`set_transform_velocities!`](@ref)). `reset_vel=false` keeps the velocities `state`
holds.
"""
function start_run!(state; reset_vel::Bool=true)
    foreach(station -> station.twist = station.row.twist, state.stations)
    reset_to_initial_pose!(state; reset_vel)
    for point in state.points
        point.apparent_mass = point.row.apparent_mass
    end
    foreach(body -> body.apparent_mass = body.row.apparent_mass, state.bodies)
    foreach(segment -> segment.l0 = segment.row.l0, state.segments)
    for tether in state.tethers
        tether.len = sum(state.segments[idx].l0 for idx in tether.segment_idxs; init=0.0)
    end
    for pulley in state.pulleys
        segment1, segment2 = (state.segments[idx] for idx in pulley.segment_idxs)
        pulley.sum_len = segment1.l0 + segment2.l0
        pulley.len = segment1.l0
        pulley.vel = 0.0
    end
    fit_body_frames!(state.points, state.bodies)
    reset_vel && set_transform_velocities!(state, state.transforms)
    return state
end

"""
    Pose(state::SystemState)

The KiteGeometry `Pose` `state` stands in: the position of each point, and the origin
and frame of each body.
"""
KiteGeometry.Pose(state::SystemState) =
    Pose([KiteGeometry.Vec3(point.pos_w) for point in state.points],
         [KiteGeometry.Vec3(body.pos_w) for body in state.bodies],
         [KiteGeometry.Mat3(body.R_b_to_w) for body in state.bodies])

"""Put the points and bodies of `state` where `pose` has them."""
function set_pose!(state::SystemState, pose)
    foreach((point, pos) -> point.pos_w .= pos, state.points, pose.pos)
    for (body, pos, R) in zip(state.bodies, pose.body_pos, pose.body_R)
        body.pos_w .= pos
        body.R_b_to_w = R
    end
    return state
end

"""
    reposition!(state::SystemState, transforms=state.transforms; update_vel=false)

Turn the points and bodies of a running `state` to the current elevation, azimuth and
heading of each of `transforms` about its base where it stands, with KiteGeometry's
`reposition!`, keeping their velocities. `update_vel` sets them to the transforms' rates
instead ([`set_transform_velocities!`](@ref)).
"""
function KiteGeometry.reposition!(state::SystemState, transforms=state.transforms;
                                  update_vel::Bool=false)
    fit_body_frames!(state.points, state.bodies)
    set_pose!(state, reposition!(Pose(state), transforms))
    update_vel && set_transform_velocities!(state, transforms)
    fit_body_frames!(state.points, state.bodies)
    return state
end

"""
    set_transform_velocities!(state::SystemState, transforms)

Set every point and body of `state` to the velocity KiteGeometry's
`transform_velocities` gives it, zero for one no transform moves; a `RIGID_DYNAMICS`
body takes the matching `ω_b` as well.
"""
function set_transform_velocities!(state::SystemState, transforms)
    point_vel, body_vel, body_spin = transform_velocities(Pose(state), transforms)
    foreach((point, vel) -> point.vel_w .= vel, state.points, point_vel)
    for (body, vel, spin) in zip(state.bodies, body_vel, body_spin)
        body.vel_w .= vel
        body.ω_b .= body.dynamics_type == RIGID_DYNAMICS ? body.R_b_to_w' * spin : 0.0
    end
    return state
end

"""
    fit_body_frames!(points, bodies)

Fit each `PARTICLE_DYNAMICS` body's frame and origin to its reference points and set its
nodes' undeformed offsets from them, and the rigid wings' nodes' offsets from the initial
pose ([`init_undeformed_offsets!`](@ref)).
"""
function fit_body_frames!(points, bodies)
    for body in bodies
        body.dynamics_type == PARTICLE_DYNAMICS || continue
        (isnothing(body.z_ref_points) || isnothing(body.y_ref_points) ||
         isnothing(body.origin)) && continue
        R_b_to_w, origin = calc_particle_dynamics_wing_frame(
            points, body.z_ref_points, body.y_ref_points, body.origin)
        body.R_b_to_w = R_b_to_w
        body.pos_w .= origin
        for point in points
            if point.is_wing_node && point.wing_idx == body.idx
                point.pos_undeformed_b .= R_b_to_w' * (point.pos_w - origin)
            end
        end
    end
    init_undeformed_offsets!(bodies, points)
end

"""
    transform_base_point(transforms, body_idx) -> Int

The base point of the transform among `transforms` moving body `body_idx`, from which its
elevation and azimuth are measured; 0 for none.
"""
function transform_base_point(transforms, body_idx)
    for transform in transforms
        body_idx in transform.bodies && return something(transform.base_point, 0)
    end
    return 0
end

# ==================== SYSSTATE INTEROP ==================== #

"""Put the world-frame aero force of each point back from `sys_state`; a log without
those channels leaves them as they are."""
function restore_point_aero_forces!(sys, sys_state)
    length(sys_state.aero_force_x) == length(sys_state.X) || return nothing
    columns = sys.sysstate
    copyto!(columns.aero_force_x, sys_state.aero_force_x)
    copyto!(columns.aero_force_y, sys_state.aero_force_y)
    copyto!(columns.aero_force_z, sys_state.aero_force_z)
    return nothing
end

"""
    update_from_sysstate!(sys::SystemState, sys_state::SysState)

Copy the state a `SysState` carries (point positions, wing orientations, winch
lengths, twist angles, segment spring forces) into a run's state, e.g. to plot one
snapshot of a `SysLog` with the Makie extension. Fields `SysState`
cannot supply (net point forces, winch forces, station moments) are set to `NaN` so
they are not plotted.

`sys` must have been created with the same model configuration as the simulation that
produced the log; its point count must match the parametric type `P` of `SysState{P}`.

Every field the log can reach is written, whether or not the compiler integrates it,
because a derived field can still seed `u0` as a torn variable.
"""
function update_from_sysstate!(sys::SystemState, sys_state::SysState{P}) where P
    (; points, segments, stations, pulleys, tethers, winches, wings, bodies) = sys

    # Position slot layout (points, panel corners, wing origins, body origins).
    slots = position_slots(sys)
    n_points = length(points)
    n_panel_corners = count_aero_log_points(wings)
    n_wings = length(wings)
    total_without_wings = n_points + n_panel_corners
    has_wing_slots = P == slots.total

    if !has_wing_slots && P != total_without_wings
        error("SystemStructure expects $(slots.total) points " *
              "($n_points regular + $n_panel_corners corners + " *
              "$n_wings wings + $(length(bodies)) bodies) or " *
              "$total_without_wings without wing slots, but SysState has $P points")
    end

    # Logs written before the velocity columns existed restart from rest.
    has_velocity = length(sys_state.VX) == length(sys_state.X)
    for point in points
        point.pos_w[1] = sys_state.X[point.idx]
        point.pos_w[2] = sys_state.Y[point.idx]
        point.pos_w[3] = sys_state.Z[point.idx]
        point.vel_w .= has_velocity ?
            (sys_state.VX[point.idx], sys_state.VY[point.idx],
             sys_state.VZ[point.idx]) : (0.0, 0.0, 0.0)
    end

    # Update wing state if wings exist
    if length(wings) > 0 && length(wings) == 1  # Currently only support single-wing systems
        wing = wings[1]

        # Copy orientation quaternion
        wing.Q_b_to_w .= sys_state.orient

        # Copy spherical coordinates
        wing.elevation = Float64(sys_state.elevation)
        wing.azimuth = Float64(sys_state.azimuth)
        wing.heading = Float64(sys_state.heading)

        if has_wing_slots
            wing_slot = slots.wings[wing.idx]
            wing.pos_w[1] = sys_state.X[wing_slot]
            wing.pos_w[2] = sys_state.Y[wing_slot]
            wing.pos_w[3] = sys_state.Z[wing_slot]
        else
            wing.pos_w .= [mean(sys_state.X), mean(sys_state.Y), mean(sys_state.Z)]
        end

        # Copy velocity if available in vel_kite
        wing.vel_w .= sys_state.vel_kite

        # Set angular velocity to NaN (turn_rates in SysState, but need conversion)
        wing.ω_b .= sys_state.turn_rates

        wing.aero_force_b .= sys_state.aero_force_KA
        wing.aero_moment_b .= sys_state.aero_moment_KA
        wing.wind_vec .= sys_state.v_wind_kite
        wing.aoa = Float64(sys_state.AoA)
        wing.course = Float64(sys_state.course)
    end

    restore_point_aero_forces!(sys, sys_state)

    # Standalone bodies: origins after the wing slots, orientations after the wings.
    if has_wing_slots
        for rigid_body in bodies
            slot = slots.bodies[rigid_body.idx]
            frame = n_wings + rigid_body.idx
            rigid_body.pos_w[1] = sys_state.X[slot]
            rigid_body.pos_w[2] = sys_state.Y[slot]
            rigid_body.pos_w[3] = sys_state.Z[slot]
            has_velocity && (rigid_body.vel_w .= (sys_state.VX[slot],
                sys_state.VY[slot], sys_state.VZ[slot]))
            rigid_body.Q_b_to_w .= sys_state.orients[frame]
            frame <= length(sys_state.turn_rate_x) &&
                (rigid_body.ω_b .= (sys_state.turn_rate_x[frame],
                    sys_state.turn_rate_y[frame], sys_state.turn_rate_z[frame]))
        end
    end

    # Logs written before twist_vel existed restart the surfaces at rest.
    has_twist_vel = length(sys_state.twist_vel) == length(stations)
    for station in stations
        i = station.idx
        if i <= length(sys_state.twist_angles)
            station.twist = Float64(sys_state.twist_angles[i])
            station.twist_ω = has_twist_vel ?
                Float64(sys_state.twist_vel[i]) : 0.0
        end
        station.aero_moment = NaN
    end
    for pulley in pulleys
        pulley.idx > length(sys_state.pulley_len) && break
        pulley.len = Float64(sys_state.pulley_len[pulley.idx])
        pulley.vel = Float64(sys_state.pulley_vel[pulley.idx])
    end

    for wing in wings
        restore_aero_twist!(wing.aero, wing, stations)
    end

    # Update tether lengths from SysState (per-tether)
    for (tether_idx, tether) in enumerate(tethers)
        tether_idx > length(sys_state.l_tether) && break
        tether.len = Float64(sys_state.l_tether[tether_idx])
    end

    # Update winch state from SysState (per-winch)
    for i in eachindex(winches)
        i > length(sys_state.v_reelout) && continue
        winches[i].vel = Float64(sys_state.v_reelout[i])
        winches[i].set_value = Float64(sys_state.set_torque[i])
    end

    corner_idx = n_points
    for wing in wings
        corner_idx = read_aero_log_points!(wing.aero, wing, sys,
                                           sys_state, corner_idx)
        restore_flap_delta!(wing.aero, wing, sys_state)
        restore_live_shape!(wing.aero, wing, sys.points)
    end

    set = sys.set
    if !isnothing(set) && set.use_wind_vec
        set.wind_vec = MVec3(sys_state.v_wind_gnd)
    end

    # Only the force is logged; `segment.len` still comes from the symbolic getters.
    has_spring_force = length(sys_state.spring_force) == length(segments)
    for segment in segments
        segment.force = has_spring_force ?
            Float64(sys_state.spring_force[segment.idx]) : NaN
    end

    return nothing
end

# ==================== DAMPING SETTERS ==================== #

"""
    assign_damping!(components, field, damping, idxs)

Write `damping` (scalar or 3-vector) into `field` of `components[idx]` for every
`idx` in `idxs`. A fresh vector reaches each component, so later edits to one
never leak into the others.
"""
function assign_damping!(components, field::Symbol,
                         damping::Union{Real, AbstractVector}, idxs)
    damp_vec = damping isa Real ?
        SVector{3,SimFloat}(damping, damping, damping) : SVector{3,SimFloat}(damping)
    for idx in idxs
        setproperty!(components[idx], field, damp_vec)
    end
    return nothing
end

"""
    set_world_frame_damping(sys::SystemStructure, damping[, point_idxs])
    set_world_frame_damping(components::AbstractVector, damping[, idxs])

Set the world-frame damping coefficient, on `sys.points` by default or on any
component vector given directly — pass `sys.bodies` to damp the rigid bodies of a
beam wing.

World-frame damping applies a velocity-dependent drag in the global reference
frame: ``\\mathbf{a}_{damp} = -c_{damp} \\odot \\mathbf{v}``, where ``c_{damp}``
is the per-mass damping vector [1/s] and ``\\odot`` is element-wise
multiplication. A scalar `damping` applies to all three axes.
"""
set_world_frame_damping(sys::SystemStructure, damping::Union{Real, AbstractVector},
                        point_idxs) =
    assign_damping!(sys.points, :world_frame_damping, damping, point_idxs)

set_world_frame_damping(sys::SystemStructure, damping::Union{Real, AbstractVector}) =
    set_world_frame_damping(sys, damping, eachindex(sys.points))

set_world_frame_damping(components::AbstractVector, damping::Union{Real, AbstractVector},
                        idxs=eachindex(components)) =
    assign_damping!(components, :world_frame_damping, damping, idxs)

"""
    set_body_frame_damping(sys::SystemStructure, damping[, point_idxs])
    set_body_frame_damping(components::AbstractVector, damping[, idxs])

Set the body-frame damping coefficient, on `sys.points` by default or on any
component vector given directly — pass `sys.bodies` to damp the rigid bodies of a
beam wing.

The velocity is resolved on the body axes before damping, so a per-axis
coefficient damps selectively: `[0, 0, 20]` on a wing body resists motion normal
to the wing without slowing flight along it. Coefficients are per-mass [1/s]; a
scalar applies to all three axes.
"""
set_body_frame_damping(sys::SystemStructure, damping::Union{Real, AbstractVector},
                       point_idxs) =
    assign_damping!(sys.points, :body_frame_damping, damping, point_idxs)

set_body_frame_damping(sys::SystemStructure, damping::Union{Real, AbstractVector}) =
    set_body_frame_damping(sys, damping, eachindex(sys.points))

set_body_frame_damping(components::AbstractVector, damping::Union{Real, AbstractVector},
                       idxs=eachindex(components)) =
    assign_damping!(components, :body_frame_damping, damping, idxs)

"""
    set_angular_damping(bodies::AbstractVector, damping[, idxs])

Set the angular damping coefficient of rigid bodies: the spin is resolved on the
body's own axes and damped per axis, `dω/dt -= c .* ω`, so `[0, 20, 0]` resists
rotation about the body `y` axis alone. Coefficients are [1/s]; a scalar applies
to all three axes.

Unlike the tube Rayleigh damping this resists the body's *absolute* spin, so it
brakes rigid rotation of whatever the body belongs to as well as deformation.
"""
set_angular_damping(bodies::AbstractVector, damping::Union{Real, AbstractVector},
                    idxs=eachindex(bodies)) =
    assign_damping!(bodies, :angular_damping, damping, idxs)

# ==================== SEGMENT STATISTICS ==================== #

"""
    segment_stretch_stats(sys::SystemStructure)

Maximum and mean relative stretch `(len - l0) / l0` over the segments in tension
(`len > l0`), plus the index of the most stretched one.

For pulley segments, the combined length of both legs is used against the pulley's
`sum_l0`, since the pulley constraint distributes length between them.

# Arguments
- `sys::SystemStructure`: System structure with current segment states

# Returns
- `(max_stretch, mean_stretch, max_idx)`: Tuple of maximum stretch, mean stretch,
  and index of the segment with maximum stretch (or first pulley segment index)
"""
function segment_stretch_stats(sys::SystemState)
    if isempty(sys.segments)
        return (0.0, 0.0, 0)
    end

    # Build set of segment indices that belong to pulleys
    pulley_seg_idxs = Set{Int64}()
    for pulley in sys.pulleys
        push!(pulley_seg_idxs, pulley.segment_idxs[1])
        push!(pulley_seg_idxs, pulley.segment_idxs[2])
    end

    stretch_data = Tuple{Int64, Float64}[]

    # Add pulley stretches (combined length of both segments)
    for pulley in sys.pulleys
        seg1 = sys.segments[pulley.segment_idxs[1]]
        seg2 = sys.segments[pulley.segment_idxs[2]]
        combined_len = seg1.len + seg2.len
        sum_l0 = pulley.sum_len  # sum_len stores seg1.l0 + seg2.l0
        if combined_len > sum_l0
            push!(stretch_data, (pulley.segment_idxs[1],
                                 (combined_len - sum_l0) / sum_l0))
        end
    end

    # Add non-pulley segment stretches
    for seg in sys.segments
        if seg.idx ∉ pulley_seg_idxs && seg.len > seg.l0
            push!(stretch_data, (seg.idx, (seg.len - seg.l0) / seg.l0))
        end
    end

    if isempty(stretch_data)
        return (0.0, 0.0, 0)
    end

    stretches = [entry[2] for entry in stretch_data]
    max_stretch = maximum(stretches)
    mean_stretch = sum(stretches) / length(stretches)
    max_idx = stretch_data[argmax(stretches)][1]

    return (max_stretch, mean_stretch, max_idx)
end

"""
    write_total_mass!(sys_struct)

Write each point's effective translational mass — `extra_mass` plus the half-mass
of every incident segment — into `point.total_mass`.
"""
write_total_mass!(::SystemState) = nothing
function write_total_mass!(sys_struct)
    for point in sys_struct.points
        point.total_mass = point.extra_mass
    end
    for segment in sys_struct.segments
        material = sys_struct.materials[segment.material]
        half = segment_half_mass(segment.l0, segment.diameter,
                                 material.parameters["density"])
        sys_struct.points[segment.point_idxs[1]].total_mass += half
        sys_struct.points[segment.point_idxs[2]].total_mass += half
    end
    return nothing
end

"""
    frame_quaternion(x_dir, ref) -> Vector{Float64}

Quaternion of the right-handed frame with `x̂` along `x_dir` and `ẑ` along the
component of `ref` orthogonal to `x̂` (`ref` must not be parallel to `x_dir`).
"""
function frame_quaternion(x_dir, ref)
    x_axis = normalize(Vector{Float64}(x_dir))
    ref_vec = Vector{Float64}(ref)
    z_axis = normalize(ref_vec - (ref_vec ⋅ x_axis) * x_axis)
    y_axis = z_axis × x_axis
    return rotation_matrix_to_quaternion(hcat(x_axis, y_axis, z_axis))
end

"""
    frame_quaternion_xy(x_dir, y_ref) -> Vector{Float64}

Quaternion of the right-handed frame with `x̂` along `x_dir` (kept exact) and `ŷ`
along the component of `y_ref` orthogonal to `x̂`, so `ẑ = x̂ × ŷ`. Use when the
primary axis is `x` — e.g. a node body whose `x̂` is the strut/chord direction and
`ŷ` the in-plane span, matching `flap_axis=[0,1,0]`. `y_ref` must not be parallel
to `x_dir`.
"""
function frame_quaternion_xy(x_dir, y_ref)
    x_axis = normalize(Vector{Float64}(x_dir))
    y_vec = Vector{Float64}(y_ref)
    y_axis = normalize(y_vec - (y_vec ⋅ x_axis) * x_axis)
    z_axis = x_axis × y_axis
    return rotation_matrix_to_quaternion(hcat(x_axis, y_axis, z_axis))
end
