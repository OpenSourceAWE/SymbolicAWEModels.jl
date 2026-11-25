# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# ==================== SHARED HELPERS ==================== #

"""
    segment_initial_length(segment, points)

Compute segment length from endpoint `pos_ENU` positions.
"""
function segment_initial_length(segment, points)
    point1 = points[segment.points[1]]
    point2 = points[segment.points[2]]
    return norm(point1.pos_ENU - point2.pos_ENU)
end

# ==================== VALIDATION ==================== #

"""
    validate_run(sam)

Error on values of the run `sam` that break it and warn on suspicious ones: NaN or
non-positive point and wing mass, wing and winch inertia, NaN positions, rigid wings
without stations, non-positive `drum_radius` or `gear_ratio`, segment `diameter` outside
(0, 1) m, non-positive `rest_length` or stiffness, negative damping, zero pulley
`sum_len`, unequal station `moment_frac`.
"""
function validate_run(sam)
    sys = sam.definition
    value(collection, i, quantity) = run_value(sam, Val(collection), i, Val(quantity))

    for (i, point) in enumerate(sys.points)
        extra_mass = value(:points, i, :extra_mass)
        isnan(extra_mass) && error("Point $(point.name) has NaN extra_mass")
        extra_mass < 0 && @warn "Point $(point.name) has negative extra_mass " *
            "$extra_mass kg. This is physically nonsensical."
        # Check mass before NaN position: NaN pos is often caused by zero mass.
        total_mass = point_total_mass(sam, i)
        if point_type(point) == DYNAMIC && total_mass <= 0
            error("Point $(point.name) has non-positive total_mass ($total_mass). " *
                  "This will cause division by zero in acceleration calculations.")
        end
        pos_w = value(:points, i, :pos_w)
        any(isnan, pos_w) && error("Point $(point.name) has NaN position: pos_w = $pos_w")
    end

    for (w, wing) in enumerate(sys.wings)
        b = wing_body(sam, w)
        mass = body_mass_properties(sam, b)
        if mass.total_mass <= 0
            error("Wing $(wing.name) has non-positive total_mass ($(mass.total_mass)). " *
                  "This will cause division by zero in acceleration calculations.")
        end
        if body_dynamics_type(sam, b) == RIGID_DYNAMICS
            I_b = mass.inertia_principal
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
            any(isnan, I_b) && error("Wing $(wing.name) has NaN inertia: I_b = $I_b")
            # AeroNone does not couple to sections, so missing stations is fine.
            isempty(wing_station_idxs(sam, w)) && couples_to_sections(sam.aero[w]) &&
                @warn "Wing $(wing.name) (RIGID_DYNAMICS) has no stations"
        end
        pos_w = value(:bodies, b, :pos_w)
        any(isnan, pos_w) && error("Wing $(wing.name) has NaN position: pos_w = $pos_w")
    end

    for (i, winch) in enumerate(sys.winches)
        inertia_total = value(:winches, i, :inertia_total)
        isnan(inertia_total) && error("Winch $(winch.name) has NaN inertia_total")
        if inertia_total ≈ 0.0
            error("Winch $(winch.name) has zero inertia_total. " *
                  "All winches must have non-zero inertia.")
        elseif inertia_total < 0
            error("Winch $(winch.name) has negative inertia_total " *
                  "$inertia_total kg⋅m². Inertia must be positive.")
        elseif inertia_total < 1e-6
            @warn "Winch $(winch.name) has very small inertia_total " *
                  "$inertia_total kg⋅m², may cause numerical issues"
        end
        drum_radius = value(:winches, i, :drum_radius)
        drum_radius <= 0 && error("Winch $(winch.name) has non-positive drum_radius " *
            "($drum_radius). This will cause division by zero in torque conversions.")
        gear_ratio = value(:winches, i, :gear_ratio)
        gear_ratio <= 0 && error("Winch $(winch.name) has non-positive gear_ratio " *
            "($gear_ratio). This will cause division by zero in speed/torque " *
            "calculations.")
    end

    for (i, segment) in enumerate(sys.segments)
        # Wing structural segments don't use diameter (stiffness explicit, drag from VSM)
        diameter = value(:segments, i, :diameter)
        if !wing_structural_segment(sys, i) && !(0 < diameter < 1)
            @warn "Segment $(segment.name) has unusual diameter $diameter m " *
                  "(expected range: 0 to 1 m)"
        end
        rest_length = value(:segments, i, :rest_length)
        rest_length <= 0 && error("Segment $(segment.name) has non-positive rest " *
            "length rest_length = $rest_length m. This will cause division by zero.")
        # Callable laws self-validate.
        stiffness = value(:segments, i, :unit_stiffness)
        if stiffness isa Real
            if stiffness ≈ 0.0
                @warn "Segment $(segment.name) has zero stiffness"
            elseif stiffness < 0
                @warn "Segment $(segment.name) has negative stiffness $stiffness N"
            end
        end
        unit_damping = value(:segments, i, :unit_damping)
        unit_damping < 0 && @warn "Segment $(segment.name) has negative damping " *
            "$unit_damping N⋅s"
    end

    for (i, pulley) in enumerate(sys.pulleys)
        sum_len = pulley_sum_len(sam, i)
        sum_len ≈ 0 && error("Pulley $(pulley.name) has zero total length constraint " *
            "(sum_len = $sum_len m). Pulley constraints must have non-zero total length.")
    end

    stations = sys.spanwise_stations
    if length(stations) > 0
        first_moment_frac = value(:spanwise_stations, 1, :moment_frac)
        for (s, station) in enumerate(stations)
            moment_frac = value(:spanwise_stations, s, :moment_frac)
            moment_frac ≈ first_moment_frac || error("Station $(station.name) has " *
                "moment_frac = $moment_frac, but all stations must have the same " *
                "moment_frac (first station has $first_moment_frac)")
        end
    end
    return nothing
end

"""
    set_unstretched_length!(sam, tether, len)

Set the unstretched length [m] of tether `tether` of the run `sam`, an index or a name,
which its segments share equally as their rest length. Point positions, body poses and
rest geometry are left as they are.
"""
function set_unstretched_length!(sam, tether, len)
    i = item_index(sam, :tethers, tether)
    name = sam.definition.tethers[i].name
    len > 0 || error("Tether $name: unstretched length $len m must be positive")
    isempty(sam.definition.tethers[i].segments) &&
        error("Tether $name: has no segments to carry an unstretched length")
    set_run_value!(sam, :tethers, i, :len, len)
    return nothing
end

"""
    seed_per_point_wind!(sam)

Give every point and every body of the run `sam` the ground wind `set.wind_vec` as its
own wind, so a [`PerPointWind`](@ref) model that is never written to flies in the uniform
wind the settings describe. Called from [`init_wind!`](@ref); from then on these winds
belong to the caller, who writes them between steps.
"""
function seed_per_point_wind!(sam)
    for i in eachindex(sam.definition.points)
        set_run_value!(sam, :points, i, :wind_vec, sam.set.wind_vec)
    end
    for w in eachindex(sam.definition.wings)
        set_run_value!(sam, :bodies, wing_body(sam, w), :wind_vec, sam.set.wind_vec)
    end
    return nothing
end

# ==================== PLACEMENT ==================== #

"""
    reset_to_initial_pose!(sam; reset_vel=true)

Put every point and body of the run `sam` back at its initial pose (`pos_ENU`,
`origin_ENU`, `Q_KA2ENU`), zero the twist of every non-`STATIC` station and set every
winch's reel-out speed to `init_vel`. `reset_vel` zeroes point and body velocities as
well.
"""
function reset_to_initial_pose!(sam; reset_vel::Bool=true)
    sys = sam.definition
    for (i, winch) in enumerate(sys.winches)
        set_run_value!(sam, :winches, i, :vel, something(winch.init_vel, 0.0))
    end
    for (i, point) in enumerate(sys.points)
        set_run_value!(sam, :points, i, :pos_w, point.pos_ENU)
        reset_vel && set_run_value!(sam, :points, i, :vel_w, zero(SVec3))
    end
    for (b, body) in enumerate(sys.bodies)
        set_run_value!(sam, :bodies, b, :pos_w, origin_ENU(sam, b))
        set_run_value!(sam, :bodies, b, :Q_b_to_w, body.Q_KA2ENU)
        if reset_vel
            set_run_value!(sam, :bodies, b, :vel_w, zero(SVec3))
            set_run_value!(sam, :bodies, b, :ω_b, zero(SVec3))
        end
    end
    for (s, station) in enumerate(sys.spanwise_stations)
        station_type(station) == STATIC && continue
        set_run_value!(sam, :spanwise_stations, s, :twist, 0.0)
        set_run_value!(sam, :spanwise_stations, s, :twist_ω, 0.0)
    end
    return nothing
end

"""
    remake_wing_aero!(sam)

Rebuild every wing's aero engine of the run `sam` from its settings and its `vsm_set`,
e.g. after editing `aero_geometry.yaml`. A no-op for aero modes without an engine.
"""
function remake_wing_aero!(sam)
    for (w, mode) in enumerate(sam.aero)
        remake_aero!(mode, sam, w)
    end
    return nothing
end

"""
    init_wind!(sam)

Set each wing's wind of the run `sam` from the ground wind `set.wind_vec` at the wing's
height, or seed every point and wing with it under [`PerPointWind`](@ref), and
initialise each wing's aerodynamic operating point from that wind.
"""
function init_wind!(sam)
    set = sam.set
    if per_point_wind(sam)
        seed_per_point_wind!(sam)
    else
        wind_factor = WindFactor(sam.am, set.profile_law)
        for w in eachindex(sam.definition.wings)
            b = wing_body(sam, w)
            height = run_value(sam, :bodies, b, :pos_w)[3]
            set_run_value!(sam, :bodies, b, :wind_vec, wind_factor(height) * set.wind_vec)
        end
    end
    for (w, mode) in enumerate(sam.aero)
        b = wing_body(sam, w)
        body_dynamics_type(sam, b) == PARTICLE_DYNAMICS && continue
        wind = per_point_wind(sam) ? run_value(sam, :bodies, b, :wind_vec) : set.wind_vec
        init_aero_state!(mode, sam, w, body_rotation(sam, b)' * wind)
    end
    return nothing
end

"""
    start_lengths!(sam)

Start each tether of the run `sam` at its segments' rest lengths, and each pulley's rope
split as its two segments' rest lengths give it, at rest.
"""
function start_lengths!(sam)
    sys = sam.definition
    rest_length(k) = something(sys.segments[k].rest_length, 0.0)
    for (i, tether) in enumerate(sys.tethers)
        set_run_value!(sam, :tethers, i, :len, sum(rest_length, tether.segments; init=0.0))
    end
    for (i, pulley) in enumerate(sys.pulleys)
        set_run_value!(sam, :pulleys, i, :len, rest_length(pulley.segments[1]))
        set_run_value!(sam, :pulleys, i, :vel, 0.0)
    end
    return sam
end

"""
    start_run!(sam; reset_vel=true)

Put the run `sam` at its start from the pose its definition is placed in: every point
and body at its initial pose ([`reset_to_initial_pose!`](@ref)), the twist of a `STATIC`
station, tether and pulley lengths from its segments' rest lengths, and points and bodies
moving as the rates of its transforms say ([`set_transform_velocities!`](@ref)).
`reset_vel=false` keeps the velocities the run holds.
"""
function start_run!(sam; reset_vel::Bool=true)
    for (s, station) in enumerate(sam.definition.spanwise_stations)
        set_run_value!(sam, :spanwise_stations, s, :twist, something(station.twist, 0.0))
    end
    reset_to_initial_pose!(sam; reset_vel)
    start_lengths!(sam)
    fit_body_frames!(sam)
    reset_vel && set_transform_velocities!(sam, sam.transforms)
    return sam
end

"""
    Pose(sam::SymbolicAWEModel)

The KiteGeometry `Pose` the run `sam` stands in: the position of each point, and the
origin and frame of each body.
"""
function KiteGeometry.Pose(sam::SymbolicAWEModel)
    sys = sam.definition
    return Pose([SVec3(run_value(sam, :points, i, :pos_w)) for i in eachindex(sys.points)],
                [SVec3(run_value(sam, :bodies, b, :pos_w)) for b in eachindex(sys.bodies)],
                [SMatrix{3, 3, Float64}(body_rotation(sam, b))
                 for b in eachindex(sys.bodies)])
end

"""Put the points and bodies of the run `sam` where `pose` has them."""
function set_pose!(sam, pose)
    for (i, pos) in enumerate(pose.pos)
        set_run_value!(sam, :points, i, :pos_w, pos)
    end
    for (b, (pos, R)) in enumerate(zip(pose.body_pos, pose.body_R))
        set_run_value!(sam, :bodies, b, :pos_w, pos)
        set_run_value!(sam, :bodies, b, :Q_b_to_w, rotation_matrix_to_quaternion(R))
    end
    return sam
end

"""
    reposition!(sam::SymbolicAWEModel, transforms=sam.transforms; update_vel=false)

Turn the points and bodies of the running `sam` to the current elevation, azimuth and
heading of each of `transforms` about its base where it stands, with KiteGeometry's
`reposition!`, keeping their velocities. `update_vel` sets them to the transforms' rates
instead ([`set_transform_velocities!`](@ref)).
"""
function KiteGeometry.reposition!(sam::SymbolicAWEModel, transforms=sam.transforms;
                                  update_vel::Bool=false)
    fit_body_frames!(sam)
    set_pose!(sam, reposition!(Pose(sam), transforms))
    update_vel && set_transform_velocities!(sam, transforms)
    fit_body_frames!(sam)
    return sam
end

"""
    set_transform_velocities!(sam, transforms)

Set every point and body of the run `sam` to the velocity KiteGeometry's
`transform_velocities` gives it, zero for one no transform moves; a `RIGID_DYNAMICS`
body takes the matching `ω_b` as well.
"""
function set_transform_velocities!(sam, transforms)
    point_vel, body_vel, body_spin = transform_velocities(Pose(sam), transforms)
    for (i, vel) in enumerate(point_vel)
        set_run_value!(sam, :points, i, :vel_w, vel)
    end
    for (b, (vel, spin)) in enumerate(zip(body_vel, body_spin))
        set_run_value!(sam, :bodies, b, :vel_w, vel)
        ω_b = body_dynamics_type(sam, b) == RIGID_DYNAMICS ?
            body_rotation(sam, b)' * spin : zero(SVec3)
        set_run_value!(sam, :bodies, b, :ω_b, ω_b)
    end
    return sam
end

"""
    fit_body_frames!(sam)

Fit the frame and origin of the body of each `PARTICLE_DYNAMICS` wing of the run `sam`
to its reference points.
"""
function fit_body_frames!(sam)
    for w in eachindex(sam.definition.wings)
        b = wing_body(sam, w)
        frame = wing_frame(sam, w)
        (body_dynamics_type(sam, b) == PARTICLE_DYNAMICS && !isnothing(frame)) || continue
        R_b_to_w, origin = KiteGeometry.wing_frame(run_positions(sam), frame)
        set_run_value!(sam, :bodies, b, :Q_b_to_w, rotation_matrix_to_quaternion(R_b_to_w))
        set_run_value!(sam, :bodies, b, :pos_w, origin)
    end
    return sam
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

"""Put the world-frame aero force of each point of the run `sam` back from `sys_state`;
a log without those channels leaves them as they are."""
function restore_point_aero_forces!(sam, sys_state)
    length(sys_state.aero_force_x) == length(sys_state.X) || return nothing
    columns = sam.sysstate
    copyto!(columns.aero_force_x, sys_state.aero_force_x)
    copyto!(columns.aero_force_y, sys_state.aero_force_y)
    copyto!(columns.aero_force_z, sys_state.aero_force_z)
    return nothing
end

"""
    update_from_sysstate!(sam::SymbolicAWEModel, sys_state::SysState)

Copy the state a `SysState` carries (point positions, wing and body poses, winch
lengths, twist angles, segment spring forces, the wings' circulation) into the run `sam`,
e.g. to plot one snapshot of a `SysLog` with the Makie extension or to resume from it,
reading each wing's pose at its wing slot and frame and every other body's at its body
slot and frame, as [`update_sys_state!`](@ref) logs them. Fields `SysState` cannot supply
(winch forces) are left as they are, a log without spring forces leaves them `NaN`.

`sam` must have been created with the same model configuration as the simulation that
produced the log; its point count must match the parametric type `P` of `SysState{P}`.
"""
function update_from_sysstate!(sam::SymbolicAWEModel, sys_state::SysState{P}) where P
    sys = sam.definition
    slots = sam.topology.slots
    n_wings = length(sys.wings)
    P == slots.total || error("$(sys.metadata.name) logs $(slots.total) positions " *
        "($(length(slots.points)) points, $(length(slots.panel_corners)) aero log " *
        "points, $n_wings wings and $(length(sys.bodies)) bodies), but the SysState " *
        "has $P")

    # Logs written before the velocity columns existed restart from rest.
    has_velocity = length(sys_state.VX) == length(sys_state.X)
    logged(column, slot) = SVec3(getfield(sys_state, column[1])[slot],
                                 getfield(sys_state, column[2])[slot],
                                 getfield(sys_state, column[3])[slot])
    for i in eachindex(sys.points)
        set_run_value!(sam, :points, i, :pos_w, logged((:X, :Y, :Z), i))
        set_run_value!(sam, :points, i, :vel_w,
                       has_velocity ? logged((:VX, :VY, :VZ), i) : zero(SVec3))
    end

    for b in eachindex(sys.bodies)
        w = body_wing(sam, b)
        slot, frame = w == 0 ? (slots.bodies[b], n_wings + b) : (slots.wings[w], w)
        set_run_value!(sam, :bodies, b, :pos_w, logged((:X, :Y, :Z), slot))
        has_velocity && set_run_value!(sam, :bodies, b, :vel_w,
                                       logged((:VX, :VY, :VZ), slot))
        set_run_value!(sam, :bodies, b, :Q_b_to_w, sys_state.orients[frame])
        frame <= length(sys_state.turn_rate_x) && set_run_value!(sam, :bodies, b, :ω_b,
            logged((:turn_rate_x, :turn_rate_y, :turn_rate_z), frame))
    end
    if n_wings == 1
        set_run_value!(sam, :wings, 1, :elevation, Float64(sys_state.elevation))
        set_run_value!(sam, :wings, 1, :azimuth, Float64(sys_state.azimuth))
        set_run_value!(sam, :wings, 1, :heading, Float64(sys_state.heading))
        set_run_value!(sam, :wings, 1, :aero_force_b, sys_state.aero_force_KA)
        set_run_value!(sam, :wings, 1, :aero_moment_b, sys_state.aero_moment_KA)
        set_run_value!(sam, :bodies, wing_body(sam, 1), :wind_vec, sys_state.v_wind_kite)
        set_run_value!(sam, :wings, 1, :aoa, Float64(sys_state.AoA))
        set_run_value!(sam, :wings, 1, :course, Float64(sys_state.course))
    end

    restore_point_aero_forces!(sam, sys_state)

    gamma = sam.sysstate.gamma_distribution
    length(sys_state.gamma_distribution) == length(gamma) &&
        copyto!(gamma, sys_state.gamma_distribution)

    # Logs written before twist_vel existed restart the surfaces at rest.
    has_twist_vel = length(sys_state.twist_vel) == length(sys.spanwise_stations)
    for s in eachindex(sys.spanwise_stations)
        s <= length(sys_state.twist_angles) || continue
        set_run_value!(sam, :spanwise_stations, s, :twist,
                       Float64(sys_state.twist_angles[s]))
        set_run_value!(sam, :spanwise_stations, s, :twist_ω,
                       has_twist_vel ? Float64(sys_state.twist_vel[s]) : 0.0)
    end
    for i in eachindex(sys.pulleys)
        i > length(sys_state.pulley_len) && break
        set_run_value!(sam, :pulleys, i, :len, Float64(sys_state.pulley_len[i]))
        set_run_value!(sam, :pulleys, i, :vel, Float64(sys_state.pulley_vel[i]))
    end

    for (w, mode) in enumerate(sam.aero)
        restore_aero_twist!(mode, sam, w)
    end

    for i in eachindex(sys.tethers)
        i > length(sys_state.l_tether) && break
        set_run_value!(sam, :tethers, i, :len, Float64(sys_state.l_tether[i]))
    end
    for i in eachindex(sys.winches)
        i > length(sys_state.v_reelout) && continue
        set_run_value!(sam, :winches, i, :vel, Float64(sys_state.v_reelout[i]))
        set_run_value!(sam, :winches, i, :set_value, Float64(sys_state.set_torque[i]))
    end

    corner_idx = length(sys.points)
    for (w, mode) in enumerate(sam.aero)
        corner_idx = read_aero_log_points!(mode, sam, w, sys_state, corner_idx)
        restore_flap_delta!(mode, sam, w, sys_state)
        restore_live_shape!(mode, sam, w)
    end

    set = sam.set
    if !isnothing(set) && set.use_wind_vec
        set.wind_vec = MVec3(sys_state.v_wind_gnd)
    end

    # Only the force is logged; a segment's length follows from its points.
    has_spring_force = length(sys_state.spring_force) == length(sys.segments)
    for i in eachindex(sys.segments)
        set_run_value!(sam, :segments, i, :force,
                       has_spring_force ? Float64(sys_state.spring_force[i]) : NaN)
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
    for idx in idxs
        setproperty!(components[idx], field, axis_column(damping))
    end
    return nothing
end

"""
    set_world_frame_damping(sys::SystemDefinition, damping[, point_idxs])
    set_world_frame_damping(components::AbstractVector, damping[, idxs])

Set the world-frame damping coefficient, on `sys.points` by default or on any
component vector given directly — pass `sys.bodies` to damp the rigid bodies of a
beam wing.

World-frame damping applies a velocity-dependent drag in the global reference
frame: ``\\mathbf{a}_{damp} = -c_{damp} \\odot \\mathbf{v}``, where ``c_{damp}``
is the per-mass damping vector [1/s] and ``\\odot`` is element-wise
multiplication. A scalar `damping` applies to all three axes.
"""
set_world_frame_damping(sys::SystemDefinition, damping::Union{Real, AbstractVector},
                        point_idxs) =
    assign_damping!(sys.points, :world_frame_damping, damping, point_idxs)

set_world_frame_damping(sys::SystemDefinition, damping::Union{Real, AbstractVector}) =
    set_world_frame_damping(sys, damping, eachindex(sys.points))

set_world_frame_damping(components::AbstractVector, damping::Union{Real, AbstractVector},
                        idxs=eachindex(components)) =
    assign_damping!(components, :world_frame_damping, damping, idxs)

"""
    set_body_frame_damping(sys::SystemDefinition, damping[, point_idxs])
    set_body_frame_damping(components::AbstractVector, damping[, idxs])

Set the body-frame damping coefficient, on `sys.points` by default or on any
component vector given directly — pass `sys.bodies` to damp the rigid bodies of a
beam wing.

The velocity is resolved on the body axes before damping, so a per-axis
coefficient damps selectively: `[0, 0, 20]` on a wing body resists motion normal
to the wing without slowing flight along it. Coefficients are per-mass [1/s]; a
scalar applies to all three axes.
"""
set_body_frame_damping(sys::SystemDefinition, damping::Union{Real, AbstractVector},
                       point_idxs) =
    assign_damping!(sys.points, :body_frame_damping, damping, point_idxs)

set_body_frame_damping(sys::SystemDefinition, damping::Union{Real, AbstractVector}) =
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
    segment_stretch_stats(sam)

Maximum and mean relative stretch `(len - rest_length) / rest_length` over the segments of
the run `sam` in tension (`len > rest_length`), plus the index of the most stretched one.
For pulley segments, the combined length of both legs is used against the pulley's
`sum_len`, since the pulley constraint distributes length between them.

# Returns
- `(max_stretch, mean_stretch, max_idx)`: Tuple of maximum stretch, mean stretch,
  and index of the segment with maximum stretch (or first pulley segment index)
"""
function segment_stretch_stats(sam)
    sys = sam.definition
    isempty(sys.segments) && return (0.0, 0.0, 0)
    pulley_seg_idxs = Set{Int64}()
    stretch_data = Tuple{Int64, Float64}[]
    for (i, pulley) in enumerate(sys.pulleys)
        union!(pulley_seg_idxs, pulley.segments)
        combined_len = sum(segment_world_length(sam, k) for k in pulley.segments)
        sum_len = pulley_sum_len(sam, i)
        combined_len > sum_len &&
            push!(stretch_data, (pulley.segments[1], (combined_len - sum_len) / sum_len))
    end
    for k in eachindex(sys.segments)
        k in pulley_seg_idxs && continue
        len = segment_world_length(sam, k)
        rest_length = run_value(sam, :segments, k, :rest_length)
        len > rest_length && push!(stretch_data, (k, (len - rest_length) / rest_length))
    end
    isempty(stretch_data) && return (0.0, 0.0, 0)
    stretches = [entry[2] for entry in stretch_data]
    return (maximum(stretches), sum(stretches) / length(stretches),
            stretch_data[argmax(stretches)][1])
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
