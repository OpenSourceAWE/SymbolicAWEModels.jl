# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# The assembler: the one layer that knows both a run's `SystemDefinition` and the
# runtime. It classifies the topology, compiles one kernel per component type,
# adds an instance per component, wires the outputs into the inputs, and binds
# each instance's parameters and initial state. Everything below it — kernel,
# instance, wiring, schedule — sees only integers and buffers.

const PARTICLE_INPUTS = [:force_in, :mass_in, :drag_in]
const PARTICLE_OUTPUTS = [:pos, :vel]
const PULLEY_POINT_INPUTS = [PARTICLE_INPUTS; :tension_in; :line_in]
const PULLEY_POINT_OUTPUTS = [:pos, :vel, :pulley_len_out, :pulley_vel_out]
const WINCH_INPUTS = [PARTICLE_INPUTS; :tension_in]
const SEGMENT_INPUTS = [:src_pos, :src_vel, :dst_pos, :dst_vel]
const SEGMENT_OUTPUTS = [:src_force, :dst_force, :half_mass, :half_drag]
const BODY_INPUTS = [:force_in, :moment_in]
const PARENTED_BODY_INPUTS = [BODY_INPUTS; :wing_frame; :wing_velocity]
const BODY_OUTPUTS = [:pos, :vel, :frame, :com, :com_velocity, :omega_w]
const RIDE_INPUTS = [:pose_pos, :pose_frame, :pose_com, :pose_com_velocity,
                     :pose_omega]
const RIDE_OUTPUTS = [:pos, :vel, :arm, :height]
const WRENCH_INPUTS = [:height, :vel, :arm, :force_in, :mass_in, :drag_in]
const WRENCH_OUTPUTS = [:force_out, :moment_out]
const JOINT_INPUTS = [:a_pos, :a_frame, :a_com, :a_com_velocity, :a_omega,
                      :b_pos, :b_frame, :b_com, :b_com_velocity, :b_omega]
const JOINT_OUTPUTS = [:force_a, :moment_a, :force_b, :moment_b]
const HERMITE_RIDE_OUTPUTS = [:pos, :vel, :arm_a, :arm_b, :height]
const HERMITE_WRENCH_INPUTS = [:height, :vel, :arm_a, :arm_b, :force_in,
                               :mass_in, :drag_in]
const KINEMATIC_INPUTS = [:z1_pos, :z2_pos, :y1_pos, :y2_pos, :origin_pos,
                          :origin_vel]
const WING_NODE_INPUTS = [:force_in, :mass_in, :drag_in, :wing_frame,
                          :wing_velocity]
"""A flap's two frames, nine scalars each rather than two vectors: a hinge axis
leaves some entries unread, and `mtkcompile` will not take part of a declared
array."""
const FLAP_INPUTS = [[Symbol(:main_frame_, k) for k in 1:9];
                     [Symbol(:flap_frame_, k) for k in 1:9]]
"""A point flap's three chord positions and the wing frame its hinge axis is in."""
const POINT_FLAP_INPUTS = [[Symbol(:point_pos_, k) for k in 1:3];
                           [Symbol(:wing_frame_, k) for k in 1:9]]
const WING_AERO_POSE_INPUTS = [:wing_pos, :wing_frame]
const AERO_INFLOW_POINT_INPUTS = [:pos, :vel, :wing_pos, :wing_frame]
const AERO_INFLOW_POINT_OUTPUTS = [:pos_b, :va_b, :rho]
const AERO_INFLOW_INPUTS = [:va_in, :rho_in]
const AERO_INFLOW_OUTPUTS = [:va, :rho]
const AERO_PANEL_INPUTS = [:le_a, :te_a, :le_b, :te_b,
                           :va_a, :va_b, :rho_a, :rho_b]
const AERO_PANEL_PITCH_INPUTS = [:dva_a, :dva_b]
const AERO_PANEL_WAGNER_INPUTS = [:wagner_deficiency]
const WAGNER_LAG_INPUTS = [:va_in]
const WAGNER_LAG_OUTPUTS = [:deficiency]
const AERO_PANEL_OUTPUTS = [:force_out, :couple_out]
const AERO_POINT_FORCE_INPUTS = [:pos_b, :wing_frame, :force_b_in]
const AERO_POINT_FORCE_OUTPUTS = [:force, :force_b, :moment_b]
const WING_AERO_SUM_INPUTS = [:force_b_in, :moment_b_in]
const TWIST_NODE_INPUTS = [RIDE_INPUTS; :twist_angle]
const TWIST_WRENCH_INPUTS = [:height, :vel, :arm, :frame, :twist_angle, :force_in,
                             :mass_in, :drag_in]
const TWIST_OUTPUTS = [:twist_angle, :twist_vel]
const STATION_INPUTS = [:aero_moment_in, :node_moment_in, :node_force_in,
                              :node_mass_in]

"""
    PointRole(kind, pulley_idx, segment_idx, winch_idx, body_idx, element_idx)

How one point is realised: `kind` is `:particle`, `:anchor`, `:pulley`, `:winch`,
`:wing_node`, `:ride` or `:hermite`. A `:pulley` point carries the pulley it splits
and one of that pulley's segments (whose material gives the rope mass); a `:winch`
point carries its winch; a `:ride` point carries the body it is anchored to, and a
`:wing_node` or `:pulley` point the body of the fitted wing it is damped against.
`element_idx` is the Timoshenko tube a `:hermite` point rides, or the station of a
`:twist_node`.
"""
struct PointRole
    kind::Symbol
    pulley_idx::Int
    segment_idx::Int
    winch_idx::Int
    body_idx::Int
    element_idx::Int
end

"""
    classify_points(sam) -> Vector{PointRole}

Decide each point's component type, following `point_eqs!`'s anchor rule: a point
anchored to a beam rides its Timoshenko element, one anchored to a body rides that
body, a point that splits a pulley is a pulley particle, a point carrying a winch is
a winch anchor, and the rest follow their `DynamicsType`. A wing node's role holds
its wing's body.
"""
function classify_points(sam)
    sys = sam.definition
    roles = [PointRole(:none, 0, 0, 0, 0, 0) for _ in sys.points]
    for (i, point) in enumerate(sys.points)
        if rigid_wing_node(sam, i)
            roles[i] = PointRole(:twist_node, 0, 0, 0,
                                 wing_body(sam, point_wing(sam, i)), station_of(sam, i))
            continue
        end
        if !isnothing(point.tube)
            roles[i] = PointRole(:hermite, 0, 0, 0, 0, point.tube)
            continue
        end
        if point_type(point) == BODY_STATIC
            !isnothing(point.body) || error(
                "KernelBackend: BODY_STATIC point $(point.name) is anchored to " *
                "neither a body nor a beam.")
            roles[i] = PointRole(:ride, 0, 0, 0, point.body, 0)
            continue
        end
        point_type(point) in (STATIC, DYNAMIC) || error(
            "KernelBackend: point $(point.name) has type $(point_type(point)); only " *
            "STATIC, DYNAMIC and BODY_STATIC are supported so far.")
        body_idx = kinematic_wing_of(sam, i)
        kind = point_type(point) == STATIC ? :anchor :
            body_idx == 0 ? :particle : :wing_node
        roles[i] = PointRole(kind, 0, 0, 0, body_idx, 0)
    end
    for (i, pulley) in enumerate(sys.pulleys)
        pulley_type(pulley) == DYNAMIC || error(
            "KernelBackend: pulley $(pulley.name) has type $(pulley_type(pulley)); " *
            "only DYNAMIC pulleys exist.")
        point = pulley_point_index(sys, i)
        roles[point] = PointRole(:pulley, i, pulley.segments[1], 0,
                                 kinematic_wing_of(sam, point), 0)
    end
    for (i, winch) in enumerate(sys.winches)
        point_type(sys.points[winch.winch_point]) == STATIC || error(
            "KernelBackend: winch $(winch.name) is at a non-STATIC point.")
        roles[winch.winch_point] = PointRole(:winch, 0, 0, i, 0, 0)
    end
    return roles
end

"""
    rigid_wing_node(sam, i) -> Bool

Whether point `i` of `sam` is a structural node of a `RIGID_DYNAMICS` wing.
`point_eqs!` takes such a node out of the anchor rule entirely — the wing *is* the
rigid body, so the node is placed by a twist-deformed body-frame offset instead of
riding anything — and so does [`classify_points`](@ref).
"""
function rigid_wing_node(sam, i)
    wing_idx = point_wing(sam, i)
    return !isnothing(wing_idx) && wing_node(sam.definition, i) &&
        body_dynamics_type(sam, wing_body(sam, wing_idx)) == RIGID_DYNAMICS
end

"""
    station_of(sam, i) -> Int

The one station point `i` of `sam` belongs to, or `0`. Two would make its section
twist ambiguous, which `point_eqs!` also rejects.
"""
function station_of(sam, i)
    stations = sam.definition.spanwise_stations
    found = [station_idx for (station_idx, station) in enumerate(stations)
             if i in station.points]
    length(found) <= 1 || error(
        "KernelBackend: point $(sam.definition.points[i].name) is in " *
        "$(length(found)) stations; expected 0 or 1.")
    return isempty(found) ? 0 : only(found)
end

"""
    kinematic_wing_of(sam, i) -> Int

The body of the fitted (`KINEMATIC`) wing whose frame point `i` of `sam` needs, or
`0`. `point_eqs!` damps *any* DYNAMIC point against its wing's frame, not only
aerodynamic surface nodes, so a steering or pulley point with `body_frame_damping`
needs it too.
"""
function kinematic_wing_of(sam, i)
    wing_idx = point_wing(sam, i)
    isnothing(wing_idx) && return 0
    body_idx = wing_body(sam, wing_idx)
    return body_type(sam, body_idx) == KINEMATIC ? body_idx : 0
end

"""
    reads_wing_frame(sys, i, role) -> Bool

Whether point `i` of definition `sys`, realised as `role`, reads its fitted wing's
frame: every wing node does, and a pulley point on such a wing when it damps in it.
"""
reads_wing_frame(sys, i, role) = role.kind === :wing_node ||
    (role.kind === :pulley && role.body_idx != 0 &&
     sys.points[i].body_frame_damping !== nothing)

"""
    pulley_point_index(sys, i) -> Int

The point pulley `i` of definition `sys` splits its rope over: the single point its
two segments share.
"""
function pulley_point_index(sys, i)
    pulley = sys.pulleys[i]
    first_points = sys.segments[pulley.segments[1]].points
    second_points = sys.segments[pulley.segments[2]].points
    shared = intersect(first_points, second_points)
    length(shared) == 1 || error(
        "KernelBackend: pulley $(pulley.name) segments share $(length(shared)) " *
        "points; expected exactly 1.")
    return only(shared)
end

"""
    SegmentRole(kind, pulley_idx, pulley_side, tether_idx, segment_count, winch_point)

How one segment is realised: `kind` is `:spring` (fixed rest length), `:structural`
(drag-free wing link), `:pulley` (rest length from a pulley split) or `:tether`
(rest length from a winch). `pulley_side` is `+1` for a pulley's first segment and
`−1` for its second; `segment_count` is its tether's segment count and
`winch_point` the winch point it touches, or `0`.
"""
struct SegmentRole
    kind::Symbol
    pulley_idx::Int
    pulley_side::Float64
    tether_idx::Int
    segment_count::Int
    winch_point::Int
end

"""
    classify_segments(sam) -> Vector{SegmentRole}

Decide each segment's component type. A pulley's two segments take their rest
length from the pulley split; a winched tether's segments take theirs from the
winch. An *unwinched* tether's rest length never changes, so its segments stay
plain springs whose `rest_length` parameter the assembler reads from the tether
instead.
"""
function classify_segments(sam)
    sys = sam.definition
    winch_point_of_tether = Dict{Int, Int}()
    for winch in sys.winches, tether_idx in winch.tethers
        winch_point_of_tether[tether_idx] = winch.winch_point
    end
    roles = SegmentRole[]
    for k in eachindex(sys.segments)
        tether_idx = sam.topology.segment_tethers[k]
        count = tether_idx == 0 ? 0 : length(sys.tethers[tether_idx].segments)
        pulley_idx, side = sam.topology.segment_pulleys[k]
        if pulley_idx > 0
            push!(roles, SegmentRole(:pulley, pulley_idx, side == 1 ? 1.0 : -1.0,
                                     0, 0, 0))
        elseif haskey(winch_point_of_tether, tether_idx)
            push!(roles, SegmentRole(:tether, 0, 0.0, tether_idx, count,
                                     winch_point_of_tether[tether_idx]))
        else
            kind = wing_structural_segment(sys, k) ? :structural : :spring
            push!(roles, SegmentRole(kind, 0, 0.0, tether_idx, count, 0))
        end
    end
    return roles
end

"""
    KernelEntry(index, registry, source)

A compiled component type: its position in the runtime's kernel table, the
`ParamRegistry` its build filled, and the component index it was built from — the
index every instance's parameter paths are remapped away from.
"""
struct KernelEntry
    index::Int
    registry::ParamRegistry
    source::Int
end

"""
    kernel!(builder, table, sam, key, source, make, inputs, outputs) -> KernelEntry

Return the kernel registered under `key`, compiling it on first use.
`make(params)` builds the component `System` from a fresh parameter view, so every
field the component reads is recorded against `source` and can be remapped to each
instance.
"""
function kernel!(builder, table, sam, key, source, make, inputs, outputs)
    entry = get(table, key, nothing)
    entry === nothing || return entry
    registry = ParamRegistry(sam)
    seconds = @elapsed kernel = compile_kernel(make(ParamView(registry)), inputs,
                                               outputs; name = key)
    fresh = KernelEntry(add_kernel!(builder, kernel, seconds), registry, source)
    builder.verbose && @printf("\t  compiled %-28s %7.2f s\n", key, seconds)
    table[key] = fresh
    return fresh
end

"""
    callable_field_key(key, component, fields) -> Symbol

`key` narrowed by which of `fields` the component holds a callable in rather than a
`Real`. A component's equations branch on that when they are built — a `Real`
rigidity becomes a numeric parameter and a callable one a callable slot — so two
components that differ there cannot share a kernel even though they are the same
kind, and writing one's law into the other's slot fails on type.
"""
function callable_field_key(key::Symbol, component, fields)
    callable = filter(field -> !(getproperty(component, field) isa Real), fields)
    isempty(callable) && return key
    return Symbol(key, :_, join(callable, '_'))
end

"""
    KernelModel

An assembled model: the runtime `system`, its initial state and parameters, the
parameter sync, and the instance index of every point, segment and body, which is
all the state getter and the control setter need to find their values. A
`BODY_STATIC` point has two: `point_instances` holds its kinematics and
`wrench_instances` the load it feeds back (`0` for every other point).
`aero_instances` holds each wing's aero instance and `twist_instances` each twist
surface's, `0` where there is none. `aero_force_outputs` holds each point's
[`AeroPointForce`](@ref) and `inflow_instances` its [`AeroInflowPoint`](@ref), whose
`va_b` output is the apparent wind the aero is solved on; `0` for a point with
neither.
"""
struct KernelModel{S, P}
    system::S
    u0::Vector{SimFloat}
    params::P
    param_sync::KernelParamSync
    point_instances::Vector{Int}
    wrench_instances::Vector{Int}
    segment_instances::Vector{Int}
    body_instances::Vector{Int}
    point_roles::Vector{PointRole}
    segment_roles::Vector{SegmentRole}
    aero_instances::Vector{Int}
    twist_instances::Vector{Int}
    inflow_instances::Vector{Int}
    aero_force_outputs::Vector{Tuple{Int, Symbol}}
end

"""
    assemble(sam) -> KernelModel

Translate the run `sam` into a [`KernelModel`](@ref). Kernels are compiled once per
component type, one instance is added per component, the wiring is recorded, and
the parameters and initial state are bound from the run afterwards, when every
instance's buffer slice is final.
"""
function assemble(sam; verbose = false)
    sys = sam.definition
    point_roles = classify_points(sam)
    segment_roles = classify_segments(sam)
    builder = SystemBuilder(; verbose)
    table = Dict{Symbol, KernelEntry}()
    bindings = Tuple{Int, KernelEntry, Dict{Symbol, Int}}[]

    body_instances = [add_body!(builder, table, bindings, sam, body_idx)
                      for body_idx in eachindex(sys.bodies)]
    twist_instances = add_stations!(builder, table, bindings, sam)
    wrench_instances = zeros(Int, length(sys.points))
    point_instances = [add_point!(builder, table, bindings, sam, i, point_roles[i],
                                  body_instances, wrench_instances, twist_instances)
                       for i in eachindex(sys.points)]
    segment_instances = [add_segment!(builder, table, bindings, sam, k,
                                      segment_roles[k])
                         for k in eachindex(sys.segments)]
    for k in eachindex(sys.segments)
        wire_segment!(builder, sys, k, segment_roles[k],
                      point_instances, segment_instances, wrench_instances)
    end
    for body_idx in eachindex(sys.bodies)
        body_type(sam, body_idx) == KINEMATIC || continue
        wire_kinematic_body!(builder, sam, body_idx, body_instances[body_idx],
                             point_instances)
    end
    for (i, role) in enumerate(point_roles)
        reads_wing_frame(sys, i, role) || continue
        connect!(builder, body_instances[role.body_idx], :frame,
                 point_instances[i], :wing_frame)
        connect!(builder, body_instances[role.body_idx], :vel,
                 point_instances[i], :wing_velocity)
    end
    for body_idx in eachindex(sys.bodies)
        parent = parent_body(sam, body_idx)
        (body_type(sam, body_idx) == DYNAMIC && parent != 0) || continue
        connect!(builder, body_instances[parent], :frame,
                 body_instances[body_idx], :wing_frame)
        connect!(builder, body_instances[parent], :vel,
                 body_instances[body_idx], :wing_velocity)
    end
    for tube_idx in eachindex(sys.tubes)
        add_tube!(builder, table, bindings, sam, tube_idx, body_instances)
    end
    flap_instances = add_flap_deltas!(builder, table, bindings, sam, body_instances,
                                      point_instances)
    inflow_instances = zeros(Int, length(sys.points))
    aero_force_outputs = fill((0, :force), length(sys.points))
    aero_instances = [add_wing_aero!(builder, table, bindings, sam, wing_idx,
                                     body_instances, point_instances, flap_instances,
                                     twist_instances, wrench_instances,
                                     inflow_instances, aero_force_outputs)
                      for wing_idx in eachindex(sys.wings)]

    system = build_system(builder)
    params, sync = bind_params(system, sam, bindings, segment_roles,
                               segment_instances)
    apply_constants!(params, system, segment_roles, segment_instances)
    u0 = zeros(SimFloat, system.n_states)
    initial_state!(u0, system, sam, point_roles, point_instances,
                   body_instances, twist_instances)
    return KernelModel(system, u0, params, sync, point_instances,
                          wrench_instances, segment_instances, body_instances,
                          point_roles, segment_roles, aero_instances,
                          twist_instances, inflow_instances, aero_force_outputs)
end

"""
    parent_body(sam, body_idx) -> Int

The body of the wing whose frame the damping of body `body_idx` of `sam` acts in, its
`wing` column, or `0` where it names none.
"""
function parent_body(sam, body_idx)
    parent_wing = sam.definition.bodies[body_idx].wing
    return isnothing(parent_wing) ? 0 :
        wing_body(sam, item_index(sam, :wings, parent_wing))
end

"""
    add_flap_deltas!(builder, table, bindings, sam, bodies, points) -> Dict{Int, Int}

Add one flap-deflection kernel per flapped station and wire what it reads: a point
flap's three chord positions and its wing's frame ([`PointFlapDelta`](@ref)), or a
body flap's two orientations ([`FlapDelta`](@ref)). Returns the instance of each such
surface; a surface with no flap is absent, and the aero input it would feed stays
unconnected and so reads the zero `station_delta_eqs!` binds it to.
"""
function add_flap_deltas!(builder, table, bindings, sam, bodies, points)
    instances = Dict{Int, Int}()
    for (station_idx, station) in enumerate(sam.definition.spanwise_stations)
        has_flap(station) || continue
        point_flap = has_point_flap(station)
        key = point_flap ? :point_flap_delta : :flap_delta
        entry = kernel!(builder, table, sam, key, station_idx,
                        params -> point_flap ?
                            PointFlapDelta(sam, params, station_idx; name = key) :
                            FlapDelta(sam, params, station_idx; name = key),
                        point_flap ? POINT_FLAP_INPUTS : FLAP_INPUTS, [:delta])
        instance = add_instance!(builder, entry.index)
        push!(bindings, (instance, entry, Dict(:spanwise_stations => station_idx)))
        if point_flap
            for (k, i) in enumerate(station.flap_points)
                connect!(builder, points[i], :pos, instance, Symbol(:point_pos_, k))
            end
            connect!(builder, bodies[wing_body(sam, station.wing)], :frame, instance,
                     :wing_frame)
        else
            main, flap = station.flap_bodies
            connect!(builder, bodies[main], :frame, instance, :main_frame)
            connect!(builder, bodies[flap], :frame, instance, :flap_frame)
        end
        instances[station_idx] = instance
    end
    return instances
end

"""
    add_wing_aero!(builder, table, bindings, sam, wing_idx, bodies, points, flaps,
                   twists, wrenches, inflow_instances, aero_force_outputs) -> Int

Add the aero instance of wing `wing_idx` and wire it, dispatching on how the wing
carries its loads: a `PARTICLE_DYNAMICS` wing's aero delivers a force to each
structural point ([`ParticleWingAero`](@ref)), a `RIGID_DYNAMICS` wing's delivers one
wrench to its body ([`add_rigid_wing_aero!`](@ref)).

Either way the aero is a component of its own, not a term inside the points or the
body: nothing about it is a cycle, so the schedule runs the structure, then the wing
frame, then the aero, then the derivatives.
"""
function add_wing_aero!(builder, table, bindings, sam, wing_idx, bodies, points, flaps,
                        twists, wrenches, inflow_instances, aero_force_outputs)
    wing_dynamics_type(sam.definition.wings[wing_idx]) == PARTICLE_DYNAMICS ||
        return add_rigid_wing_aero!(builder, table, bindings, sam, wing_idx, bodies,
                                    twists)
    supports_panel_decomposition(sam.aero[wing_idx]) && return add_panel_wing_aero!(
        builder, table, bindings, sam, wing_idx, bodies, points, flaps, wrenches,
        inflow_instances, aero_force_outputs)
    nodes = wing_points(sam, wing_idx)
    surfaces = wing_flap_surfaces(sam, wing_idx)
    body = bodies[wing_body(sam, wing_idx)]
    key = Symbol(:wing_aero_, wing_idx)
    inputs = [WING_AERO_POSE_INPUTS
              [Symbol(:point_pos_, k) for k in eachindex(nodes)]
              [Symbol(:point_vel_, k) for k in eachindex(nodes)]
              [Symbol(:flap_delta_, surface) for surface in surfaces]]
    outputs = [Symbol(:point_force_, k) for k in eachindex(nodes)]
    entry = kernel!(builder, table, sam, key, wing_idx,
                    params -> ParticleWingAero(sam, params, wing_idx; name = key),
                    inputs, outputs)
    instance = add_instance!(builder, entry.index)
    push!(bindings, (instance, entry, Dict{Symbol, Int}()))
    connect!(builder, body, :pos, instance, :wing_pos)
    connect!(builder, body, :frame, instance, :wing_frame)
    for (k, node) in enumerate(nodes)
        connect!(builder, points[node], :pos, instance, Symbol(:point_pos_, k))
        connect!(builder, points[node], :vel, instance, Symbol(:point_vel_, k))
        connect!(builder, instance, Symbol(:point_force_, k),
                 load_target(points, wrenches, node), :force_in)
        aero_force_outputs[node] = (instance, Symbol(:point_force_, k))
    end
    for surface in surfaces
        haskey(flaps, surface) || continue
        connect!(builder, flaps[surface], :delta, instance,
                 Symbol(:flap_delta_, surface))
    end
    return instance
end

"""
    add_panel_wing_aero!(builder, table, bindings, sam, wing_idx, bodies, points,
                         flaps, wrenches, inflow_instances, aero_force_outputs) -> Int

Add the aerodynamics of `PARTICLE_DYNAMICS` wing `wing_idx` as small repeated
components: one [`AeroInflowPoint`](@ref) per structural point, one
[`AeroInflow`](@ref) per inflow group, one [`AeroPanel`](@ref) per refined panel, one
[`AeroPointForce`](@ref) per point again, and one [`WingAeroSum`](@ref) for the
readouts. Everything between them — the strut interpolation, the inflow average, the
load scatter — is a constant weight, so it is wiring rather than equations. Returns
the sum's instance, which is where the wing's readouts hang.

The one-component equivalent, [`ParticleWingAero`](@ref), is superlinear in the wing's
size; this decomposition is not, which is what makes a large wing buildable.
"""
function add_panel_wing_aero!(builder, table, bindings, sam, wing_idx, bodies, points,
                              flaps, wrenches, inflow_instances,
                              aero_force_outputs)
    mode = sam.aero[wing_idx]
    nodes = wing_points(sam, wing_idx)
    body_idx = wing_body(sam, wing_idx)
    body = bodies[body_idx]
    inflow_points = add_aero_inflow_points!(builder, table, bindings, sam, nodes,
                                            points, body_idx, body)
    for (node, instance) in zip(nodes, inflow_points)
        inflow_instances[node] = instance
    end
    groups, section_group = aero_inflow_groups(mode, sam, wing_idx)
    inflows = add_aero_inflows!(builder, table, bindings, sam, groups, inflow_points)
    pitches = if flow_curvature_enabled(mode)
        add_aero_inflows!(builder, table, bindings, sam,
                          aero_pitch_groups(mode, sam, wing_idx)[1], inflow_points)
    else
        nothing
    end
    wagner = wagner_enabled(mode) ?
        add_wagner_lag!(builder, table, bindings, sam, wing_idx, inflow_points) :
        nothing
    panels = add_aero_panels!(builder, table, bindings, sam, wing_idx, inflow_points,
                              inflows, pitches, section_group, flaps, wagner)
    forces = add_aero_point_forces!(builder, table, bindings, sam, wing_idx, nodes,
                                    inflow_points, body, points, wrenches)
    for (node, instance) in zip(nodes, forces)
        aero_force_outputs[node] = (instance, :force)
    end
    for (panel, node, force_weight, couple_weight) in
            aero_scatter_entries(mode, sam, wing_idx)
        force_weight == 0 || connect!(builder, panels[panel], :force_out,
                                      forces[node], :force_b_in; weight = force_weight)
        couple_weight == 0 || connect!(builder, panels[panel], :couple_out,
                                       forces[node], :force_b_in;
                                       weight = couple_weight)
    end
    return add_wing_aero_sum!(builder, table, bindings, sam, forces)
end

"""
    add_aero_inflow_points!(builder, table, bindings, sam, nodes, points, body_idx,
                            body)

One [`AeroInflowPoint`](@ref) per structural point of a wing, wired from that point's
kinematics and the pose of the wing's body `body_idx`, instance `body`. The only
per-component field the kernel reads is the point's prescribed wind, remapped per
instance, so every point of every wing shares it. Returns the instances, indexed as
`nodes` is.
"""
function add_aero_inflow_points!(builder, table, bindings, sam, nodes, points, body_idx,
                                 body)
    source = first(nodes)
    entry = kernel!(builder, table, sam, :aero_inflow_point, source,
                    params -> AeroInflowPoint(sam, params, source;
                                              name = :aero_inflow_point),
                    AERO_INFLOW_POINT_INPUTS, AERO_INFLOW_POINT_OUTPUTS)
    instances = Int[]
    for node in nodes
        instance = add_instance!(builder, entry.index)
        push!(bindings, (instance, entry, Dict(:points => node, :bodies => body_idx)))
        connect!(builder, points[node], :pos, instance, :pos)
        connect!(builder, points[node], :vel, instance, :vel)
        connect!(builder, body, :pos, instance, :wing_pos)
        connect!(builder, body, :frame, instance, :wing_frame)
        push!(instances, instance)
    end
    return instances
end

"""
    add_aero_inflows!(builder, table, bindings, sam, groups, inflow_points) -> Vector{Int}

One [`AeroInflow`](@ref) per group of [`aero_inflow_groups`](@ref), gathering its
points' apparent wind and density at the group's weights.
"""
function add_aero_inflows!(builder, table, bindings, sam, groups, inflow_points)
    entry = kernel!(builder, table, sam, :aero_inflow, 0,
                    params -> AeroInflow(; name = :aero_inflow),
                    AERO_INFLOW_INPUTS, AERO_INFLOW_OUTPUTS)
    instances = Int[]
    for group in groups
        instance = add_instance!(builder, entry.index)
        push!(bindings, (instance, entry, Dict{Symbol, Int}()))
        for (node, weight) in group
            connect!(builder, inflow_points[node], :va_b, instance, :va_in; weight)
            connect!(builder, inflow_points[node], :rho, instance, :rho_in; weight)
        end
        push!(instances, instance)
    end
    return instances
end

"""
    add_aero_panels!(builder, table, bindings, sam, wing_idx, inflow_points, inflows,
                     pitches, section_group, flaps, wagner=nothing) -> Vector{Int}

One [`AeroPanel`](@ref) per refined panel of wing `wing_idx`, reading its two
sections' corners from the strut interpolation ([`aero_geometry_entries`](@ref)),
their inflow from the group they belong to, and its flap deflection from the station
it deflects with. Every panel of a wing shares one kernel.

`pitches` are the [`aero_pitch_groups`](@ref) gathers, or `nothing` on a wing without
[`flow_curvature_enabled`](@ref), which then has no such inputs to connect. `wagner`
is the wing's [`add_wagner_lag!`](@ref) instance, or `nothing` in the same way.
"""
function add_aero_panels!(builder, table, bindings, sam, wing_idx, inflow_points,
                          inflows, pitches, section_group, flaps, wagner=nothing)
    mode = sam.aero[wing_idx]
    with_flap = !isempty(wing_flap_surfaces(sam, wing_idx))
    inputs = isnothing(pitches) ? AERO_PANEL_INPUTS :
        [AERO_PANEL_INPUTS; AERO_PANEL_PITCH_INPUTS]
    isnothing(wagner) || (inputs = [inputs; AERO_PANEL_WAGNER_INPUTS])
    with_flap && (inputs = [inputs; :flap_delta])
    instances = Int[]
    key = Symbol(:aero_panel_, wing_idx)
    for panel_idx in 1:Int(mode.vsm_wing.n_panels)
        entry = kernel!(builder, table, sam, key, panel_idx,
                        params -> AeroPanel(sam, params, wing_idx, panel_idx;
                                            name = key, with_flap),
                        inputs, AERO_PANEL_OUTPUTS)
        instance = add_instance!(builder, entry.index)
        push!(bindings, (instance, entry, Dict(:panels => panel_idx)))
        for (side, section) in ((:a, panel_idx), (:b, panel_idx + 1))
            inflow = inflows[section_group[section]]
            connect!(builder, inflow, :va, instance, Symbol(:va_, side))
            connect!(builder, inflow, :rho, instance, Symbol(:rho_, side))
            isnothing(pitches) || connect!(builder, pitches[section_group[section]],
                                           :va, instance, Symbol(:dva_, side))
        end
        isnothing(wagner) || connect!(builder, wagner, :deficiency, instance,
                                      :wagner_deficiency)
        push!(instances, instance)
    end
    for (panel, corner, node, weight) in aero_geometry_entries(mode, sam, wing_idx)
        connect!(builder, inflow_points[node], :pos_b, instances[panel], corner; weight)
    end
    with_flap && wire_panel_flaps!(builder, mode, instances, flaps)
    return instances
end

"""
    add_wagner_lag!(builder, table, bindings, sam, wing_idx, inflow_points) -> Int

The one [`WagnerLag`](@ref) instance of wing `wing_idx`, its `va_in` gathered at equal
weight over every node of the wing so the lag rides the wing's mean apparent wind.
Returns the instance the panels read their deficiency from.
"""
function add_wagner_lag!(builder, table, bindings, sam, wing_idx, inflow_points)
    key = Symbol(:wagner_lag_, wing_idx)
    entry = kernel!(builder, table, sam, key, 0,
                    params -> WagnerLag(sam, params, wing_idx; name = key),
                    WAGNER_LAG_INPUTS, WAGNER_LAG_OUTPUTS)
    instance = add_instance!(builder, entry.index)
    push!(bindings, (instance, entry, Dict{Symbol, Int}()))
    weight = 1.0 / length(inflow_points)
    for point in inflow_points
        connect!(builder, point, :va_b, instance, :va_in; weight)
    end
    return instance
end

"""
    wire_panel_flaps!(builder, mode, panels, flaps)

Connect each panel's flap deflection to the [`FlapDelta`](@ref) of the station it
deflects with. A panel mapped to no surface keeps its input unconnected, which reads
zero — the same deflection `flap_delta_inputs` gives it.
"""
function wire_panel_flaps!(builder, mode, panels, flaps)
    for (panel, surface) in enumerate(mode.panel_station)
        haskey(flaps, surface) || continue
        connect!(builder, flaps[surface], :delta, panels[panel], :flap_delta)
    end
    return nothing
end

"""
    add_aero_point_forces!(builder, table, bindings, sam, wing_idx, nodes,
                           inflow_points, body, points, wrenches) -> Vector{Int}

One [`AeroPointForce`](@ref) per structural point of wing `wing_idx`, gathering the
panels' scattered body-frame load and delivering the world force to the point (or to
its wrench half, when the point rides a body).
"""
function add_aero_point_forces!(builder, table, bindings, sam, wing_idx, nodes,
                                inflow_points, body, points, wrenches)
    key = Symbol(:aero_point_force_, wing_idx)
    instances = Int[]
    for (k, node) in enumerate(nodes)
        entry = kernel!(builder, table, sam, key, node,
                        params -> AeroPointForce(sam, params, wing_idx, node;
                                                 name = key),
                        AERO_POINT_FORCE_INPUTS, AERO_POINT_FORCE_OUTPUTS)
        instance = add_instance!(builder, entry.index)
        push!(bindings, (instance, entry, Dict(:points => node)))
        connect!(builder, inflow_points[k], :pos_b, instance, :pos_b)
        connect!(builder, body, :frame, instance, :wing_frame)
        connect!(builder, instance, :force, load_target(points, wrenches, node),
                 :force_in)
        push!(instances, instance)
    end
    return instances
end

"""
    add_wing_aero_sum!(builder, table, bindings, sam, forces) -> Int

The [`WingAeroSum`](@ref) gathering every point's body-frame force and moment, whose
observables are the wing's `aero_force_b` and `aero_moment_b` readouts.
"""
function add_wing_aero_sum!(builder, table, bindings, sam, forces)
    entry = kernel!(builder, table, sam, :wing_aero_sum, 0,
                    params -> WingAeroSum(; name = :wing_aero_sum),
                    WING_AERO_SUM_INPUTS, Symbol[])
    instance = add_instance!(builder, entry.index)
    push!(bindings, (instance, entry, Dict{Symbol, Int}()))
    for force in forces
        connect!(builder, force, :force_b, instance, :force_b_in)
        connect!(builder, force, :moment_b, instance, :moment_b_in)
    end
    return instance
end

"""
    add_rigid_wing_aero!(builder, table, bindings, sam, wing_idx, bodies, twists) -> Int

Add the [`WingAero`](@ref) of `RIGID_DYNAMICS` wing `wing_idx` and wire it: its
body's pose and each of its stations' angle and rate in, the world wrench about the
body COM back into the body, and each surface's aerodynamic hinge moment on to that
surface. Returns the instance. `aero_eqs!` drives every surface's hinge moment
except a prescribed one with no aero sections, which has nothing to drive it with.
"""
function add_rigid_wing_aero!(builder, table, bindings, sam, wing_idx, bodies, twists)
    surfaces = wing_station_idxs(sam, wing_idx)
    body = bodies[wing_body(sam, wing_idx)]
    key = Symbol(:wing_aero_, wing_idx)
    inputs = [RIDE_INPUTS
              [Symbol(:twist_angle_, surface) for surface in surfaces]
              [Symbol(:twist_vel_, surface) for surface in surfaces]]
    outputs = [WRENCH_OUTPUTS
               [Symbol(:twist_moment_, surface) for surface in surfaces]]
    entry = kernel!(builder, table, sam, key, wing_idx,
                    params -> WingAero(sam, params, wing_idx; name = key),
                    inputs, outputs)
    instance = add_instance!(builder, entry.index)
    push!(bindings, (instance, entry, Dict{Symbol, Int}()))
    for (source, target) in ((:pos, :pose_pos), (:frame, :pose_frame),
                             (:com, :pose_com), (:com_velocity, :pose_com_velocity),
                             (:omega_w, :pose_omega))
        connect!(builder, body, source, instance, target)
    end
    connect!(builder, instance, :force_out, body, :force_in)
    connect!(builder, instance, :moment_out, body, :moment_in)
    for surface in surfaces
        twists[surface] == 0 && continue
        connect!(builder, twists[surface], :twist_angle, instance,
                 Symbol(:twist_angle_, surface))
        connect!(builder, twists[surface], :twist_vel, instance,
                 Symbol(:twist_vel_, surface))
        station_aero_driven(sam, surface) &&
            connect!(builder, instance, Symbol(:twist_moment_, surface),
                     twists[surface], :aero_moment_in)
    end
    return instance
end

"""
    add_body!(builder, table, bindings, sam, body_idx) -> Int

Add the instance for body `body_idx`. Each of the three kinds is its own compiled
type, because which one a body is, is topology: a `DYNAMIC` body integrates, a
`STATIC` one is clamped, and a `KINEMATIC` one is fitted from reference points.
"""
function add_body!(builder, table, bindings, sam, body_idx)
    kind = body_type(sam, body_idx)
    kind == KINEMATIC && return add_kinematic_body!(builder, table, bindings, sam,
                                                    body_idx)
    kind in (DYNAMIC, STATIC) || error(
        "KernelBackend: body $(sam.definition.bodies[body_idx].name) has type $kind; " *
        "only DYNAMIC, STATIC and KINEMATIC are supported so far.")
    parented = kind != STATIC && parent_body(sam, body_idx) != 0
    key = kind == STATIC ? :static_body :
        parented ? :parented_rigid_body : :rigid_body
    make = kind == STATIC ?
        (params -> StaticBody(sam, params, body_idx; name = key)) :
        (params -> RigidBody(sam, params, body_idx; name = key, parented))
    entry = kernel!(builder, table, sam, key, body_idx, make,
                    parented ? PARENTED_BODY_INPUTS : BODY_INPUTS, BODY_OUTPUTS)
    instance = add_instance!(builder, entry.index)
    push!(bindings, (instance, entry, Dict(:bodies => body_idx)))
    return instance
end

"""
    add_kinematic_body!(builder, table, bindings, sam, body_idx) -> Int

Add a `KINEMATIC` wing body. It has no state: its frame is fitted from four
structural reference points and its origin pose read from a fifth, so the whole
component is wiring. Each reference is a weighted blend of real points, delivered by
the weights in the [`Wiring`](@ref).
"""
function add_kinematic_body!(builder, table, bindings, sam, body_idx)
    entry = kernel!(builder, table, sam, :kinematic_body, body_idx,
                    params -> KinematicBody(sam, params, body_idx;
                                            name = :kinematic_body),
                    KINEMATIC_INPUTS, BODY_OUTPUTS)
    instance = add_instance!(builder, entry.index)
    push!(bindings, (instance, entry, Dict(:bodies => body_idx)))
    return instance
end

"""
    wire_kinematic_body!(builder, sam, body_idx, body_instance, point_instances)

Wire the reference points of the fitted wing whose body is `body_idx` into it: the
four frame references and the origin's position and velocity, each a weighted blend
of point outputs.
"""
function wire_kinematic_body!(builder, sam, body_idx, body_instance, point_instances)
    frame = wing_frame(sam, body_wing(sam, body_idx))
    references = ((frame.z_ref_points[1], :z1_pos, :pos),
                  (frame.z_ref_points[2], :z2_pos, :pos),
                  (frame.y_ref_points[1], :y1_pos, :pos),
                  (frame.y_ref_points[2], :y2_pos, :pos),
                  (frame.origin, :origin_pos, :pos),
                  (frame.origin, :origin_vel, :vel))
    for (reference, target, source) in references, (point, weight) in reference
        connect!(builder, point_instances[point], source, body_instance, target; weight)
    end
    return nothing
end

"""
    add_point!(builder, table, bindings, sam, i, role, bodies, wrenches, twists) -> Int

Add the instance for point `i` and record which container indices its parameters
must be remapped to. A `:ride` point becomes two instances — the kinematics, whose
index is returned so segments wire to it, and the wrench, recorded in `wrenches` and
wired to its body here.
"""
function add_point!(builder, table, bindings, sam, i, role, bodies, wrenches, twists)
    sys = sam.definition
    index_map = Dict(:points => i)
    role.kind === :ride && return add_ride_point!(builder, table, bindings, sam, i,
                                                  role, bodies, wrenches)
    role.kind === :twist_node && return add_twist_node!(builder, table, bindings, sam,
                                                        i, role, bodies, wrenches,
                                                        twists)
    role.kind === :hermite && return add_hermite_ride_point!(builder, table, bindings,
                                                             sam, i, role, bodies,
                                                             wrenches)
    with_damping = sys.points[i].body_frame_damping !== nothing
    entry = if role.kind === :wing_node
        key = with_damping ? :wing_node : :wing_node_undamped
        index_map[:bodies] = role.body_idx
        kernel!(builder, table, sam, key, i,
                params -> WingNodePoint(sam, params, i; name = key, with_damping),
                WING_NODE_INPUTS, PARTICLE_OUTPUTS)
    elseif role.kind === :particle
        kernel!(builder, table, sam, :particle, i,
                params -> Particle(sam, params, i; name = :particle),
                PARTICLE_INPUTS, PARTICLE_OUTPUTS)
    elseif role.kind === :anchor
        kernel!(builder, table, sam, :anchor, i,
                params -> Anchor(sam, params, i; name = :anchor),
                PARTICLE_INPUTS, PARTICLE_OUTPUTS)
    elseif role.kind === :pulley
        index_map[:pulleys] = role.pulley_idx
        index_map[:segments] = role.segment_idx
        damped = reads_wing_frame(sys, i, role)
        key = damped ? :pulley_point_wing : :pulley_point
        inputs = damped ? [PULLEY_POINT_INPUTS; :wing_frame; :wing_velocity] :
            PULLEY_POINT_INPUTS
        kernel!(builder, table, sam, key, i,
                params -> PulleyParticle(sam, params, i, role.pulley_idx,
                                         role.segment_idx; name = key,
                                         with_damping = damped),
                inputs, PULLEY_POINT_OUTPUTS)
    elseif role.kind === :winch
        tethers = sys.winches[role.winch_idx].tethers
        index_map[:winches] = role.winch_idx
        key = Symbol(:winch_, role.winch_idx)
        outputs = [PARTICLE_OUTPUTS;
                   [Symbol(:tether_len_, k) for k in eachindex(tethers)];
                   [Symbol(:tether_vel_, k) for k in eachindex(tethers)]]
        kernel!(builder, table, sam, key, i,
                params -> WinchAnchor(sam, params, role.winch_idx, i; name = key),
                WINCH_INPUTS, outputs)
    else
        error("KernelBackend: point $(sys.points[i].name) has no component for role " *
              "$(role.kind)")
    end
    instance = add_instance!(builder, entry.index)
    push!(bindings, (instance, entry, index_map))
    return instance
end

"""
    add_ride_point!(builder, table, bindings, sam, i, role, bodies, wrenches)

Add the two instances `BODY_STATIC` point `i` needs and wire them to its body: the
[`RidePoint`](@ref) reads the body's pose, the [`RideWrench`](@ref) reads the ride
point's `pos`/`vel`/`arm` and delivers force and moment back. Returns the ride
point, which is what the incident segments connect to.
"""
function add_ride_point!(builder, table, bindings, sam, i, role, bodies, wrenches)
    body = bodies[role.body_idx]
    kinematics = kernel!(builder, table, sam, :ride_point, i,
                         params -> RidePoint(sam, params, i; name = :ride_point),
                         RIDE_INPUTS, RIDE_OUTPUTS)
    gravity = carrier_body_idx(sam, i) == 0
    key = gravity ? :ride_wrench : :riding_wrench
    statics = kernel!(builder, table, sam, key, i,
                      params -> RideWrench(sam, params, i; name = key,
                                           with_gravity = gravity),
                      WRENCH_INPUTS, WRENCH_OUTPUTS)
    ride = add_instance!(builder, kinematics.index)
    wrench = add_instance!(builder, statics.index)
    push!(bindings, (ride, kinematics, Dict(:points => i)))
    push!(bindings, (wrench, statics, Dict(:points => i)))
    wrenches[i] = wrench
    for (source, target) in ((:pos, :pose_pos), (:vel, :pose_pos),
                             (:frame, :pose_frame), (:com, :pose_com),
                             (:com_velocity, :pose_com_velocity),
                             (:omega_w, :pose_omega))
        source === :vel && continue
        connect!(builder, body, source, ride, target)
    end
    connect!(builder, ride, :height, wrench, :height)
    connect!(builder, ride, :vel, wrench, :vel)
    connect!(builder, ride, :arm, wrench, :arm)
    connect!(builder, wrench, :force_out, body, :force_in)
    connect!(builder, wrench, :moment_out, body, :moment_in)
    return ride
end

"""
    add_stations!(builder, table, bindings, sam) -> Vector{Int}

Add one twist instance per station that has a section twist to report: a
[`StationDOF`](@ref) for a `DYNAMIC` surface, whose twist is a state, and a
[`PrescribedTwist`](@ref) for a `STATIC` one, whose twist is a parameter. A
`KINEMATIC` surface has neither and gets `0`; its deflection is a
[`FlapDelta`](@ref) instead.
"""
function add_stations!(builder, table, bindings, sam)
    stations = sam.definition.spanwise_stations
    instances = zeros(Int, length(stations))
    for (station_idx, station) in enumerate(stations)
        station_type(station) in (DYNAMIC, STATIC) || continue
        dynamic = station_type(station) == DYNAMIC
        key = dynamic ? :station : :prescribed_twist
        make = dynamic ? StationDOF : PrescribedTwist
        entry = kernel!(builder, table, sam, key, station_idx,
                        params -> make(sam, params, station_idx; name = key),
                        dynamic ? STATION_INPUTS : [:aero_moment_in],
                        TWIST_OUTPUTS)
        instance = add_instance!(builder, entry.index)
        push!(bindings, (instance, entry, Dict(:spanwise_stations => station_idx)))
        instances[station_idx] = instance
    end
    return instances
end

"""
    add_twist_node!(builder, table, bindings, sam, i, role, bodies, wrenches, twists)

Add the two instances structural node `i` of a `RIGID_DYNAMICS` wing needs and wire
them: the [`TwistNodePoint`](@ref) places it from the wing body's pose and its
surface's twist, the [`TwistNodeWrench`](@ref) sends its load and moment back to the
body and its bridle couple and mass on to the surface. Returns the node, which is
what the incident segments connect to. Only a `DYNAMIC` surface takes the couple; a
prescribed one holds its twist whatever the nodes pull.
"""
function add_twist_node!(builder, table, bindings, sam, i, role, bodies, wrenches,
                         twists)
    body = bodies[role.body_idx]
    surface = role.element_idx
    gated = surface > 0 &&
        !run_value(sam, :bodies, role.body_idx, :group_points_moment)
    index_map = Dict(:points => i, :bodies => role.body_idx)
    surface > 0 && (index_map[:spanwise_stations] = surface)
    suffix = surface > 0 ? "" : "_free"
    kinematics = kernel!(builder, table, sam, Symbol(:twist_node, suffix), i,
                         params -> TwistNodePoint(sam, params, i;
                                                  name = :twist_node,
                                                  surface_idx = surface),
                         TWIST_NODE_INPUTS, RIDE_OUTPUTS)
    outputs = surface > 0 ?
              [WRENCH_OUTPUTS; :node_force; :node_moment; :node_mass] :
              WRENCH_OUTPUTS
    key = Symbol(:twist_wrench, suffix, gated ? "_ungated" : "")
    statics = kernel!(builder, table, sam, key, i,
                      params -> TwistNodeWrench(sam, params, i; name = key,
                                                surface_idx = surface, gated),
                      TWIST_WRENCH_INPUTS, outputs)
    node = add_instance!(builder, kinematics.index)
    wrench = add_instance!(builder, statics.index)
    push!(bindings, (node, kinematics, index_map))
    push!(bindings, (wrench, statics, index_map))
    wrenches[i] = wrench
    for (source, target) in ((:pos, :pose_pos), (:frame, :pose_frame),
                             (:com, :pose_com), (:com_velocity, :pose_com_velocity),
                             (:omega_w, :pose_omega))
        connect!(builder, body, source, node, target)
    end
    connect!(builder, body, :frame, wrench, :frame)
    connect!(builder, node, :height, wrench, :height)
    connect!(builder, node, :vel, wrench, :vel)
    connect!(builder, node, :arm, wrench, :arm)
    connect!(builder, wrench, :force_out, body, :force_in)
    connect!(builder, wrench, :moment_out, body, :moment_in)
    if surface > 0 && twists[surface] != 0
        connect!(builder, twists[surface], :twist_angle, node, :twist_angle)
        connect!(builder, twists[surface], :twist_angle, wrench, :twist_angle)
        if station_type(sam.definition.spanwise_stations[surface]) == DYNAMIC
            connect!(builder, wrench, :node_force, twists[surface], :node_force_in)
            connect!(builder, wrench, :node_moment, twists[surface], :node_moment_in)
            connect!(builder, wrench, :node_mass, twists[surface], :node_mass_in)
        end
    end
    return node
end

"""
    add_hermite_ride_point!(builder, table, bindings, sam, i, role, bodies, wrenches)

Add the two instances beam-anchored point `i` needs and wire them to the two end
bodies of the Timoshenko element it rides: the [`HermiteRidePoint`](@ref) reads both
poses, the [`HermiteRideWrench`](@ref) reads the ride point's `pos`/`vel` and two
moment arms and delivers a share of its load to each body. Returns the ride point,
which is what the incident segments connect to.
"""
function add_hermite_ride_point!(builder, table, bindings, sam, i, role, bodies,
                                 wrenches)
    tube_idx = role.element_idx
    index_map = Dict(:points => i, :tubes => tube_idx)
    kinematics = kernel!(builder, table, sam, :hermite_ride_point, i,
                         params -> HermiteRidePoint(sam, params, i;
                                                    name = :hermite_ride_point),
                         JOINT_INPUTS, HERMITE_RIDE_OUTPUTS)
    statics = kernel!(builder, table, sam, :hermite_ride_wrench, i,
                      params -> HermiteRideWrench(sam, params, i;
                                                  name = :hermite_ride_wrench),
                      HERMITE_WRENCH_INPUTS, JOINT_OUTPUTS)
    ride = add_instance!(builder, kinematics.index)
    wrench = add_instance!(builder, statics.index)
    push!(bindings, (ride, kinematics, index_map))
    push!(bindings, (wrench, statics, index_map))
    wrenches[i] = wrench
    body_a, body_b = sam.definition.tubes[tube_idx].bodies
    for (prefix, body) in ((:a, body_a), (:b, body_b))
        for source in (:pos, :frame, :com, :com_velocity)
            connect!(builder, bodies[body], source, ride, Symbol(prefix, :_, source))
        end
        connect!(builder, bodies[body], :omega_w, ride, Symbol(prefix, :_omega))
        connect!(builder, ride, Symbol(:arm_, prefix), wrench, Symbol(:arm_, prefix))
        connect!(builder, wrench, Symbol(:force_, prefix), bodies[body], :force_in)
        connect!(builder, wrench, Symbol(:moment_, prefix), bodies[body], :moment_in)
    end
    connect!(builder, ride, :height, wrench, :height)
    connect!(builder, ride, :vel, wrench, :vel)
    return ride
end

"""
    add_tube!(builder, table, bindings, sam, tube_idx, bodies)

Add tube `tube_idx` and wire it: both end bodies' poses in, the restoring wrench on
each back out. Its row model picks the component ([`tube_kernel`](@ref)); its parameters
are remapped over `tubes`.
"""
function add_tube!(builder, table, bindings, sam, tube_idx, bodies)
    columns = tube_columns(sam.definition.tubes[tube_idx])
    kind, make = tube_kernel(columns)
    fields = rigidity_fields(columns)
    rigidities = NamedTuple{fields}(map(field -> run_value(sam, :tubes, tube_idx, field),
                                        fields))
    key = callable_field_key(kind, rigidities, fields)
    entry = kernel!(builder, table, sam, key, tube_idx,
                    params -> make(sam, params, tube_idx; name = kind),
                    JOINT_INPUTS, JOINT_OUTPUTS)
    instance = add_instance!(builder, entry.index)
    push!(bindings, (instance, entry, Dict(:tubes => tube_idx)))
    body_a, body_b = sam.definition.tubes[tube_idx].bodies
    for (prefix, body) in ((:a, body_a), (:b, body_b))
        for (source, target) in ((:pos, :pos), (:frame, :frame), (:com, :com),
                                 (:com_velocity, :com_velocity),
                                 (:omega_w, :omega))
            connect!(builder, bodies[body], source, instance,
                     Symbol(prefix, :_, target))
        end
        connect!(builder, instance, Symbol(:force_, prefix), bodies[body], :force_in)
        connect!(builder, instance, Symbol(:moment_, prefix), bodies[body],
                 :moment_in)
    end
    return instance
end

"""
    tube_kernel(columns) -> (kind, make)

The kernel kind a tube with the row model `columns` is built as and the component
constructor building it.
"""
tube_kernel(::ElasticTubeModel) = (:elastic_tube, ElasticTubeComponent)
tube_kernel(::TimoshenkoTubeModel) = (:timoshenko_tube, TimoshenkoTubeComponent)

"""
    add_segment!(builder, table, bindings, sam, k, role) -> Int

Add the instance for segment `k` and record its parameter remapping.
"""
function add_segment!(builder, table, bindings, sam, k, role)
    index_map = Dict(:segments => k)
    stiffness = (:unit_stiffness,)
    segment = (; unit_stiffness = run_value(sam, :segments, k, :unit_stiffness))
    entry = if role.kind === :pulley
        index_map[:pulleys] = role.pulley_idx
        kernel!(builder, table, sam,
                callable_field_key(:pulley_segment, segment, stiffness), k,
                params -> PulleySegment(sam, params, k, role.pulley_idx;
                                        name = :pulley_segment),
                [SEGMENT_INPUTS; :rest_len; :rest_vel],
                [SEGMENT_OUTPUTS; :tension; :line])
    elseif role.kind === :tether
        kernel!(builder, table, sam,
                callable_field_key(:tether_segment, segment, stiffness), k,
                params -> TetherSegment(sam, params, k; name = :tether_segment),
                [SEGMENT_INPUTS; :rest_len; :rest_vel],
                [SEGMENT_OUTPUTS; :src_tension; :dst_tension])
    else
        with_drag = role.kind === :spring
        key = with_drag ? :spring_segment : :structural_segment
        kernel!(builder, table, sam, callable_field_key(key, segment, stiffness),
                k, params -> SpringSegment(sam, params, k; name = key, with_drag),
                SEGMENT_INPUTS, SEGMENT_OUTPUTS)
    end
    instance = add_instance!(builder, entry.index)
    push!(bindings, (instance, entry, index_map))
    return instance
end

"""
    wire_segment!(builder, sys, k, role, point_instances, segment_instances,
                  wrench_instances)

Connect segment `k` of definition `sys` to its two endpoints: their pose into its
inputs, its endpoint force, mass share and drag share back into theirs — a clamped
endpoint included, which does not move in response but does report the load it
carries. A pulley segment also exchanges the rope split and its tension with the
pulley point; a tether segment takes its rest length from the winch and delivers its
tension there.
"""
function wire_segment!(builder, sys, k, role, point_instances, segment_instances,
                       wrench_instances)
    segment = segment_instances[k]
    source, target = sys.segments[k].points
    source == target && error(
        "KernelBackend: segment $(sys.segments[k].name) has both endpoints at point " *
        "$source.")
    for (endpoint, prefix) in ((source, :src), (target, :dst))
        point = point_instances[endpoint]
        connect!(builder, point, :pos, segment, Symbol(prefix, :_pos))
        connect!(builder, point, :vel, segment, Symbol(prefix, :_vel))
        loaded = load_target(point_instances, wrench_instances, endpoint)
        connect!(builder, segment, :half_drag, loaded, :drag_in)
        connect!(builder, segment, Symbol(prefix, :_force), loaded, :force_in)
        connect!(builder, segment, :half_mass, loaded, :mass_in)
    end
    if role.kind === :pulley
        pulley_point = point_instances[pulley_point_index(sys, role.pulley_idx)]
        connect!(builder, pulley_point, :pulley_len_out, segment, :rest_len)
        connect!(builder, pulley_point, :pulley_vel_out, segment, :rest_vel)
        connect!(builder, segment, :tension, pulley_point, :tension_in)
        connect!(builder, segment, :line, pulley_point, :line_in)
    elseif role.kind === :tether
        winch_point = point_instances[role.winch_point]
        winch = sys.winches[find_winch(sys, role.winch_point)]
        position = findfirst(==(role.tether_idx), winch.tethers)
        connect!(builder, winch_point, Symbol(:tether_len_, position), segment,
                 :rest_len)
        connect!(builder, winch_point, Symbol(:tether_vel_, position), segment,
                 :rest_vel)
        source == role.winch_point &&
            connect!(builder, segment, :src_tension, winch_point, :tension_in)
        target == role.winch_point &&
            connect!(builder, segment, :dst_tension, winch_point, :tension_in)
    end
    return nothing
end

"""
    load_target(points, wrenches, idx) -> Int

The instance a point's loads are delivered to: its wrench half when it is anchored
to a body or a beam, otherwise the point itself. An anchored point is two
components, and only the second one accepts force, mass and drag.
"""
load_target(points, wrenches, idx) = wrenches[idx] == 0 ? points[idx] : wrenches[idx]

"""The index of the winch of definition `sys` sitting at point `i`."""
find_winch(sys, i) = findfirst(winch -> winch.winch_point == i, sys.winches)

"""
    buffer_slots(system, instance, buffer, name) -> Vector{Int}

The global indices of `name` in `instance`'s slice of `buffer` (`:states`,
`:inputs`, `:outputs`, `:observables` or `:params`).
"""
function buffer_slots(system, instance::Int, buffer::Symbol, name::Symbol)
    inst = system.instances[instance]
    kernel = system.kernels[inst.kernel]
    return first(getfield(inst, buffer)) .- 1 .+ slots(getfield(kernel, buffer), name)
end

"""
    winch_state_slots(system, instance, name, alias) -> Vector{Int}

The state slots of `name` on a winch anchor, taking `alias` when `name` is not the
unknown that survived. [`WinchAnchor`](@ref) binds `motor.vel` to its own
`winch_vel` and, for a single tether, `motor.len` to `tether_len_1`; either member
of such a pair can be the one `mtkcompile` keeps, so the reader may not assume.
`alias` is `nothing` where no pair exists.
"""
function winch_state_slots(system, instance::Int, name::Symbol, alias)
    kernel = system.kernels[system.instances[instance].kernel]
    key = alias === nothing || has_slot(kernel.states, name) ? name : alias
    return buffer_slots(system, instance, :states, key)
end

"""The [`sync_class`](@ref) of the parameter `reader` fills: `:run` for a computed one."""
function reader_class(reader)
    inner = path_reader(reader)
    return isnothing(inner) ? :run : sync_class(inner.path)
end

"""
    bind_params(system, sam, bindings, segment_roles, segment_instances)

Fill the parameter buffer with every kernel's build-time defaults, then record the
live reader for each parameter that reads a value of the run `sam`, remapped to its
own instance. Returns `(KernelParams, KernelParamSync)`.
"""
function bind_params(system, sam, bindings, segment_roles, segment_instances)
    numeric = zeros(SimFloat, system.n_params)
    callables = callable_store(system)
    slots = Int[]
    readers = Any[]
    callable_targets = Any[]
    callable_readers = Any[]
    for (instance, entry, index_map) in bindings
        inst = system.instances[instance]
        kernel = system.kernels[inst.kernel]
        numeric[inst.params] .= kernel.param_defaults
        numeric_bound, callable_bound = instance_readers(kernel, entry.registry,
                                                         index_map)
        for (slot, reader) in numeric_bound
            push!(slots, first(inst.params) - 1 + slot)
            push!(readers, reader)
        end
        for (slot, reader) in callable_bound
            push!(callable_targets,
                  CallableSlot(callables[inst.kernel], inst.position, slot))
            push!(callable_readers, reader)
        end
    end
    retarget_tether_rest_lengths!(readers, segment_roles)
    bind_segment_winds!(slots, readers, system, sam, segment_instances)
    classes = map(reader_class, readers)
    of(class) = group_readers(slots[classes .=== class], readers[classes .=== class])
    sync = KernelParamSync(of(:live), of(:aero), of(:run),
                           group_callables(callable_targets, callable_readers))
    return KernelParams(numeric, callables), sync
end

"""
    callable_store(system) -> Tuple

The callable parameters, seeded from each kernel's build-time defaults: one entry
per kernel, holding one *tuple* per instance of that kernel. The tuple is what
keeps a polar call inside a generated kernel statically dispatched — a vector would
widen to the join of the slots' types, and a wing panel's `cl`, `cd` and `cm` are
three different types, so every coefficient lookup would become a dynamic dispatch.
"""
function callable_store(system)
    counts = zeros(Int, length(system.kernels))
    for inst in system.instances
        counts[inst.kernel] += 1
    end
    return ntuple(length(system.kernels)) do k
        defaults = Tuple(system.kernels[k].callable_defaults)
        fill(defaults, counts[k])
    end
end

"""
    bind_segment_winds!(slots, readers, system, sam, segment_instances)

Point each segment's `src_wind`/`dst_wind` parameters at its two endpoints'
`wind_vec`, which is how a [`PerPointWind`](@ref) segment reads the wind of the
points it spans: a segment kernel is instanced over `:segments`, so its endpoints
cannot both be reached through the registry's per-instance index remapping.
[`segment_wind_params`](@ref) mints these parameters, so nothing else binds them; a
segment without tether drag has none.
"""
function bind_segment_winds!(slots, readers, system, sam, segment_instances)
    per_point_wind(sam) || return nothing
    for (k, instance) in enumerate(segment_instances)
        kernel = system.kernels[system.instances[instance].kernel]
        has_slot(kernel.params, :src_wind) || continue
        endpoints = sam.definition.segments[k].points
        for (endpoint, name) in zip(endpoints, (:src_wind, :dst_wind))
            for (k, slot) in enumerate(buffer_slots(system, instance, :params, name))
                push!(slots, slot)
                push!(readers, PathReader((:points, endpoint, :wind_vec, k)))
            end
        end
    end
    return nothing
end

"""
    retarget_tether_rest_lengths!(readers, segment_roles)

Point an unwinched tether segment's rest-length reader at its tether's length
instead of the segment's own `l0`, matching `segment_eqs!`'s `l0 = tether_len /
n_segs` for a tether with no winch to integrate that length. The reader is found by
the path it reads, which names the segment and so is unique.
"""
function retarget_tether_rest_lengths!(readers, segment_roles)
    for (idx, role) in enumerate(segment_roles)
        (role.tether_idx > 0 && role.kind !== :tether) || continue
        path = (:segments, idx, :rest_length)
        position = findfirst(reader -> reader isa PathReader && reader.path == path,
                             readers)
        position === nothing && error(
            "KernelBackend: segment $idx has no rest-length parameter to retarget")
        readers[position] = ScaledReader(
            PathReader((:tethers, role.tether_idx, :len)), 1 / role.segment_count)
    end
    return nothing
end

"""
    apply_constants!(params, system, segment_roles, segment_instances)

Write the parameters fixed by the topology straight into the buffer: a pulley
segment's side and a tether segment's segment count. They never change, so no
reader syncs them each step.
"""
function apply_constants!(params, system, segment_roles, segment_instances)
    for (idx, role) in enumerate(segment_roles)
        instance = segment_instances[idx]
        if role.kind === :pulley
            slot = only(buffer_slots(system, instance, :params, :pulley_side))
            params.numeric[slot] = role.pulley_side
        elseif role.kind === :tether
            slot = only(buffer_slots(system, instance, :params, :segment_count))
            params.numeric[slot] = role.segment_count
        end
    end
    return nothing
end

"""
    initial_state!(u0, system, sam, point_roles, point_instances, body_instances,
                   twist_instances)

Fill `u0` with the initial state of the run `sam`: each `DYNAMIC` station's twist and
rate, each body's principal pose, each particle's `pos`/`vel`, each pulley's split
`pulley_len`/`pulley_vel` and each winch's `winch_vel` and per-tether lengths. Called
again by [`KernelInitialSync`](@ref) whenever a problem is reused, so a slot the run
does not reach has to be zeroed here.
"""
function initial_state!(u0, system, sam, point_roles, point_instances, body_instances,
                        twist_instances)
    sys = sam.definition
    fill!(u0, zero(SimFloat))
    states(instance, name) = buffer_slots(system, instance, :states, name)
    for (station_idx, instance) in enumerate(twist_instances)
        instance == 0 && continue
        station_type(sys.spanwise_stations[station_idx]) == DYNAMIC || continue
        u0[only(states(instance, :free_twist_angle))] =
            run_value(sam, :spanwise_stations, station_idx, :twist)
        u0[only(states(instance, :twist_omega))] =
            run_value(sam, :spanwise_stations, station_idx, :twist_ω)
    end
    for (body_idx, instance) in enumerate(body_instances)
        # Only a DYNAMIC body integrates; a clamped or fitted one has no state.
        body_type(sam, body_idx) == DYNAMIC || continue
        u0[states(instance, :com_w)] .= run_value(sam, :bodies, body_idx, :com_w)
        u0[states(instance, :com_vel)] .= run_value(sam, :bodies, body_idx, :com_vel)
        u0[states(instance, :Q)] .= run_value(sam, :bodies, body_idx, :Q_p_to_w)
        u0[states(instance, :omega_p)] .= run_value(sam, :bodies, body_idx, :ω_p)
    end
    for (i, role) in enumerate(point_roles)
        instance = point_instances[i]
        if role.kind in (:particle, :pulley, :wing_node)
            u0[states(instance, :pos)] .= run_value(sam, :points, i, :pos_w)
            u0[states(instance, :vel)] .= run_value(sam, :points, i, :vel_w)
        end
        if role.kind === :pulley
            u0[only(states(instance, :pulley_len))] =
                run_value(sam, :pulleys, role.pulley_idx, :len)
            u0[only(states(instance, :pulley_vel))] =
                run_value(sam, :pulleys, role.pulley_idx, :vel)
        elseif role.kind === :winch
            tethers = sys.winches[role.winch_idx].tethers
            single = length(tethers) == 1
            u0[only(winch_state_slots(system, instance, :winch_vel,
                                      Symbol("motor₊vel")))] =
                run_value(sam, :winches, role.winch_idx, :vel)
            for (k, tether_idx) in enumerate(tethers)
                alias = single ? Symbol("motor₊len") : nothing
                slot = only(winch_state_slots(system, instance,
                                              Symbol(:tether_len_, k), alias))
                u0[slot] = run_value(sam, :tethers, tether_idx, :len)
            end
        end
    end
    return u0
end
