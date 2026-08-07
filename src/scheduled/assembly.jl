# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# The assembler: the one layer that knows both the `SystemStructure` and the
# runtime. It classifies the topology, compiles one kernel per component type,
# adds an instance per component, wires the outputs into the inputs, and binds
# each instance's parameters and initial state. Everything below it — kernel,
# instance, wiring, schedule — sees only integers and buffers.

const PARTICLE_INPUTS = [:force_in, :mass_in, :drag_in]
const PARTICLE_OUTPUTS = [:pos, :vel]
const ANCHOR_INPUTS = [:drag_in]
const PULLEY_POINT_INPUTS = [:force_in, :mass_in, :drag_in, :tension_in]
const PULLEY_POINT_OUTPUTS = [:pos, :vel, :pulley_len_out]
const WINCH_INPUTS = [:drag_in, :tension_in]
const SEGMENT_INPUTS = [:src_pos, :src_vel, :dst_pos, :dst_vel]
const SEGMENT_OUTPUTS = [:src_force, :dst_force, :half_mass, :half_drag]
const BODY_INPUTS = [:force_in, :moment_in]
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
const FITTED_INPUTS = [:z1_pos, :z2_pos, :y1_pos, :y2_pos, :origin_pos,
                       :origin_vel]
const WING_NODE_INPUTS = [:force_in, :mass_in, :drag_in, :wing_frame,
                          :wing_velocity]

"""
    PointRole(kind, pulley_idx, segment_idx, winch_idx, body_idx, joint_idx)

How one point is realised: `kind` is `:particle`, `:anchor`, `:pulley`, `:winch`,
`:wing_node`, `:ride` or `:hermite`. A `:pulley` point carries the pulley it splits
and one of that pulley's segments (whose material gives the rope mass); a `:winch`
point carries its winch; a `:ride` point carries the body it is anchored to and a
`:hermite` point the Timoshenko joint whose beam it rides.
"""
struct PointRole
    kind::Symbol
    pulley_idx::Int
    segment_idx::Int
    winch_idx::Int
    body_idx::Int
    joint_idx::Int
end

"""
    classify_points(sys_struct) -> Vector{PointRole}

Decide each point's component type, following `point_eqs!`'s anchor rule: a point
anchored to a beam rides its Timoshenko element, one anchored to a body rides that
body, a point that splits a pulley is a pulley particle, a point carrying a winch is
a winch anchor, and the rest follow their `DynamicsType`.
"""
function classify_points(sys_struct)
    roles = [PointRole(:none, 0, 0, 0, 0, 0) for _ in sys_struct.points]
    for (i, point) in enumerate(sys_struct.points)
        if point.joint_idx > 0
            roles[i] = PointRole(:hermite, 0, 0, 0, 0, point.joint_idx)
            continue
        end
        if point.type == BODY_STATIC
            point.body_idx > 0 || error(
                "ScheduledBackend: BODY_STATIC point $(point.name) is anchored to " *
                "neither a body nor a beam.")
            roles[i] = PointRole(:ride, 0, 0, 0, point.body_idx, 0)
            continue
        end
        point.type in (STATIC, DYNAMIC) || error(
            "ScheduledBackend: point $(point.name) has type $(point.type); only " *
            "STATIC, DYNAMIC and BODY_STATIC are supported so far.")
        wing = fitted_wing_of(sys_struct, point)
        kind = point.type == STATIC ? :anchor : wing == 0 ? :particle : :wing_node
        roles[i] = PointRole(kind, 0, 0, 0, wing, 0)
    end
    for pulley in sys_struct.pulleys
        pulley.type == DYNAMIC || error(
            "ScheduledBackend: pulley $(pulley.name) has type $(pulley.type); " *
            "only DYNAMIC pulleys exist.")
        index = pulley_point_index(sys_struct, pulley)
        roles[index] = PointRole(:pulley, pulley.idx, pulley.segment_idxs[1],
                                 0, 0, 0)
    end
    for winch in sys_struct.winches
        index = winch.winch_point_idx
        sys_struct.points[index].type == STATIC || error(
            "ScheduledBackend: winch $(winch.name) is at a non-STATIC point.")
        roles[index] = PointRole(:winch, 0, 0, 0, winch.idx, 0)
    end
    return roles
end

"""
    fitted_wing_of(sys_struct, point) -> Int

The fitted (`KINEMATIC`) wing whose frame `point` needs, or `0`. `point_eqs!` damps
*any* DYNAMIC point against its wing's frame, not only aerodynamic surface nodes, so
a steering or pulley point with `body_frame_damping` needs it too.
"""
function fitted_wing_of(sys_struct, point)
    idx = point.wing_idx
    (idx > 0 && idx <= length(sys_struct.bodies)) || return 0
    sys_struct.bodies[idx].type == KINEMATIC || return 0
    needs = point.is_wing_node || point.body_frame_damping !== nothing
    return needs ? idx : 0
end

"""
    pulley_point_index(sys_struct, pulley) -> Int

The point a pulley splits its rope over: the single point its two segments share.
"""
function pulley_point_index(sys_struct, pulley)
    first_points = sys_struct.segments[pulley.segment_idxs[1]].point_idxs
    second_points = sys_struct.segments[pulley.segment_idxs[2]].point_idxs
    shared = intersect(first_points, second_points)
    length(shared) == 1 || error(
        "ScheduledBackend: pulley $(pulley.name) segments share $(length(shared)) " *
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
    classify_segments(sys_struct) -> Vector{SegmentRole}

Decide each segment's component type. A pulley's two segments take their rest
length from the pulley split; a winched tether's segments take theirs from the
winch. An *unwinched* tether's rest length never changes, so its segments stay
plain springs whose `l0` parameter the assembler reads from the tether instead.
"""
function classify_segments(sys_struct)
    tether_of_segment = Dict{Int, Int}()
    for tether in sys_struct.tethers, index in tether.segment_idxs
        tether_of_segment[index] = tether.idx
    end
    winch_point_of_tether = Dict{Int, Int}()
    for winch in sys_struct.winches, index in winch.tether_idxs
        winch_point_of_tether[index] = winch.winch_point_idx
    end
    side_of_segment = Dict{Int, Tuple{Int, Float64}}()
    for pulley in sys_struct.pulleys
        side_of_segment[pulley.segment_idxs[1]] = (pulley.idx, 1.0)
        side_of_segment[pulley.segment_idxs[2]] = (pulley.idx, -1.0)
    end
    roles = SegmentRole[]
    for segment in sys_struct.segments
        tether_idx = get(tether_of_segment, segment.idx, 0)
        count = tether_idx == 0 ? 0 :
            length(sys_struct.tethers[tether_idx].segment_idxs)
        if haskey(side_of_segment, segment.idx)
            pulley_idx, side = side_of_segment[segment.idx]
            push!(roles, SegmentRole(:pulley, pulley_idx, side, 0, 0, 0))
        elseif haskey(winch_point_of_tether, tether_idx)
            push!(roles, SegmentRole(:tether, 0, 0.0, tether_idx, count,
                                     winch_point_of_tether[tether_idx]))
        else
            kind = wing_structural_segment(sys_struct, segment.idx) ?
                :structural : :spring
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
    registry = ParamRegistry(sam.sys_struct)
    kernel = compile_kernel(make(ParamView(registry)), inputs, outputs; name = key)
    fresh = KernelEntry(add_kernel!(builder, kernel), registry, source)
    table[key] = fresh
    return fresh
end

"""
    ScheduledModel

An assembled model: the runtime `system`, its initial state and parameters, the
parameter sync, and the instance index of every point, segment and body, which is
all the state getter and the control setter need to find their values. A
`BODY_STATIC` point has two: `point_instances` holds its kinematics and
`wrench_instances` the load it feeds back (`0` for every other point).
"""
struct ScheduledModel{S, P}
    system::S
    u0::Vector{SimFloat}
    params::P
    param_sync::ScheduledParamSync
    point_instances::Vector{Int}
    wrench_instances::Vector{Int}
    segment_instances::Vector{Int}
    body_instances::Vector{Int}
    point_roles::Vector{PointRole}
    segment_roles::Vector{SegmentRole}
end

"""
    assemble(sam) -> ScheduledModel

Translate `sam.sys_struct` into a scheduled model. Kernels are compiled once per
component type, one instance is added per component, the wiring is recorded, and
the parameters and initial state are bound from the struct afterwards, when every
instance's buffer slice is final.
"""
function assemble(sam)
    sys_struct = sam.sys_struct
    point_roles = classify_points(sys_struct)
    segment_roles = classify_segments(sys_struct)
    builder = SystemBuilder()
    table = Dict{Symbol, KernelEntry}()
    bindings = Tuple{Int, KernelEntry, Dict{Symbol, Int}}[]

    body_instances = [add_body!(builder, table, bindings, sam, i)
                      for i in eachindex(sys_struct.bodies)]
    wrench_instances = zeros(Int, length(sys_struct.points))
    point_instances = [add_point!(builder, table, bindings, sam, i, point_roles[i],
                                  body_instances, wrench_instances)
                       for i in eachindex(sys_struct.points)]
    segment_instances = [add_segment!(builder, table, bindings, sam, i,
                                      segment_roles[i])
                         for i in eachindex(sys_struct.segments)]
    for i in eachindex(sys_struct.segments)
        wire_segment!(builder, sys_struct, i, segment_roles[i], point_roles,
                      point_instances, segment_instances, wrench_instances)
    end
    for (idx, body) in enumerate(sys_struct.bodies)
        body.type == KINEMATIC || continue
        wire_fitted_body!(builder, sys_struct, idx, body_instances[idx],
                          point_instances)
    end
    for (idx, role) in enumerate(point_roles)
        role.kind === :wing_node || continue
        connect!(builder, body_instances[role.body_idx], :frame,
                 point_instances[idx], :wing_frame)
        connect!(builder, body_instances[role.body_idx], :vel,
                 point_instances[idx], :wing_velocity)
    end
    for joint in sys_struct.elastic_joints
        add_joint!(builder, table, bindings, sam, joint, body_instances,
                   :elastic_joints, ElasticJointComponent)
    end
    for joint in sys_struct.timoshenko_joints
        add_joint!(builder, table, bindings, sam, joint, body_instances,
                   :timoshenko_joints, TimoshenkoJointComponent)
    end

    system = build_system(builder)
    params, sync = bind_params(system, sys_struct, bindings, segment_roles,
                               segment_instances)
    apply_constants!(params, system, segment_roles, segment_instances)
    u0 = initial_state(system, sys_struct, point_roles, point_instances,
                       body_instances)
    return ScheduledModel(system, u0, params, sync, point_instances,
                          wrench_instances, segment_instances, body_instances,
                          point_roles, segment_roles)
end

"""
    add_body!(builder, table, bindings, sam, idx) -> Int

Add the instance for body `idx`. Each of the three kinds is its own compiled type,
because which one a body is, is topology: a `DYNAMIC` body integrates, a `STATIC`
one is clamped, and a `KINEMATIC` one is fitted from reference points.
"""
function add_body!(builder, table, bindings, sam, idx)
    body = sam.sys_struct.bodies[idx]
    body.type == KINEMATIC && return add_fitted_body!(builder, table, bindings,
                                                      sam, idx)
    body.type in (DYNAMIC, STATIC) || error(
        "ScheduledBackend: body $(body.name) has type $(body.type); only DYNAMIC, " *
        "STATIC and KINEMATIC are supported so far.")
    make = body.type == STATIC ? StaticBody : RigidBody
    key = body.type == STATIC ? :static_body : :rigid_body
    entry = kernel!(builder, table, sam, key, idx,
                    params -> make(sam, params, idx; name = key),
                    BODY_INPUTS, BODY_OUTPUTS)
    instance = add_instance!(builder, entry.index)
    push!(bindings, (instance, entry, Dict(:bodies => idx)))
    return instance
end

"""
    add_fitted_body!(builder, table, bindings, sam, idx) -> Int

Add a `KINEMATIC` wing body. It has no state: its frame is fitted from four
structural reference points and its origin pose read from a fifth, so the whole
component is wiring. Each reference is a weighted blend of real points, delivered by
the weights in the [`Wiring`](@ref).
"""
function add_fitted_body!(builder, table, bindings, sam, idx)
    entry = kernel!(builder, table, sam, :fitted_body, idx,
                    params -> KinematicBody(sam, params, idx; name = :fitted_body),
                    FITTED_INPUTS, BODY_OUTPUTS)
    instance = add_instance!(builder, entry.index)
    push!(bindings, (instance, entry, Dict(:bodies => idx)))
    return instance
end

"""
    wire_fitted_body!(builder, sys_struct, idx, body_instance, point_instances)

Wire a fitted wing's reference points into it: the four frame references and the
origin's position and velocity, each a weighted blend of point outputs.
"""
function wire_fitted_body!(builder, sys_struct, idx, body_instance, point_instances)
    wing = sys_struct.bodies[idx]
    references = ((wing.z_ref_points[1], :z1_pos, :pos),
                  (wing.z_ref_points[2], :z2_pos, :pos),
                  (wing.y_ref_points[1], :y1_pos, :pos),
                  (wing.y_ref_points[2], :y2_pos, :pos),
                  (wing.origin, :origin_pos, :pos),
                  (wing.origin, :origin_vel, :vel))
    for (reference, target, source) in references
        for (position, point) in enumerate(reference.ids)
            weight = length(reference.ids) == 1 ? 1.0 : reference.weights[position]
            connect!(builder, point_instances[point], source, body_instance, target;
                     weight)
        end
    end
    return nothing
end

"""
    add_point!(builder, table, bindings, sam, idx, role, bodies, wrenches) -> Int

Add the instance for point `idx` and record which container indices its parameters
must be remapped to. A `:ride` point becomes two instances — the kinematics, whose
index is returned so segments wire to it, and the wrench, recorded in `wrenches` and
wired to its body here.
"""
function add_point!(builder, table, bindings, sam, idx, role, bodies, wrenches)
    sys_struct = sam.sys_struct
    index_map = Dict(:points => idx)
    role.kind === :ride && return add_ride_point!(builder, table, bindings, sam, idx,
                                                  role, bodies, wrenches)
    role.kind === :hermite && return add_hermite_ride_point!(builder, table, bindings,
                                                             sam, idx, role, bodies,
                                                             wrenches)
    entry = if role.kind === :wing_node
        point = sys_struct.points[idx]
        with_aero = point.is_wing_node && is_wing(sys_struct.bodies[point.wing_idx])
        with_damping = point.body_frame_damping !== nothing
        key = Symbol(:wing_node, with_aero ? "" : "_free",
                     with_damping ? "" : "_undamped")
        kernel!(builder, table, sam, key, idx,
                params -> WingNodePoint(sam, params, idx; name = key, with_aero,
                                        with_damping),
                WING_NODE_INPUTS, PARTICLE_OUTPUTS)
    elseif role.kind === :particle
        kernel!(builder, table, sam, :particle, idx,
                params -> Particle(sam, params, idx; name = :particle),
                PARTICLE_INPUTS, PARTICLE_OUTPUTS)
    elseif role.kind === :anchor
        kernel!(builder, table, sam, :anchor, idx,
                params -> Anchor(sam, params, idx; name = :anchor),
                ANCHOR_INPUTS, PARTICLE_OUTPUTS)
    elseif role.kind === :pulley
        index_map[:pulleys] = role.pulley_idx
        index_map[:segments] = role.segment_idx
        kernel!(builder, table, sam, :pulley_point, idx,
                params -> PulleyParticle(sam, params, idx, role.pulley_idx,
                                         role.segment_idx; name = :pulley_point),
                PULLEY_POINT_INPUTS, PULLEY_POINT_OUTPUTS)
    else
        winch = sys_struct.winches[role.winch_idx]
        index_map[:winches] = role.winch_idx
        key = Symbol(:winch_, role.winch_idx)
        outputs = [PARTICLE_OUTPUTS;
                   [Symbol(:tether_len_, k) for k in eachindex(winch.tether_idxs)]]
        kernel!(builder, table, sam, key, idx,
                params -> WinchAnchor(sam, params, winch, idx; name = key),
                WINCH_INPUTS, outputs)
    end
    instance = add_instance!(builder, entry.index)
    push!(bindings, (instance, entry, index_map))
    return instance
end

"""
    add_ride_point!(builder, table, bindings, sam, idx, role, bodies, wrenches)

Add the two instances a `BODY_STATIC` point needs and wire them to its body: the
[`RidePoint`](@ref) reads the body's pose, the [`RideWrench`](@ref) reads the ride
point's `pos`/`vel`/`arm` and delivers force and moment back. Returns the ride
point, which is what the incident segments connect to.
"""
function add_ride_point!(builder, table, bindings, sam, idx, role, bodies, wrenches)
    body = bodies[role.body_idx]
    point = sam.sys_struct.points[idx]
    kinematics = kernel!(builder, table, sam, :ride_point, idx,
                         params -> RidePoint(sam, params, idx; name = :ride_point),
                         RIDE_INPUTS, RIDE_OUTPUTS)
    gravity = !wing_frame_member(point, point.body_idx)
    key = gravity ? :ride_wrench : :riding_wrench
    statics = kernel!(builder, table, sam, key, idx,
                      params -> RideWrench(sam, params, idx; name = key,
                                           with_gravity = gravity),
                      WRENCH_INPUTS, WRENCH_OUTPUTS)
    ride = add_instance!(builder, kinematics.index)
    wrench = add_instance!(builder, statics.index)
    push!(bindings, (ride, kinematics, Dict(:points => idx)))
    push!(bindings, (wrench, statics, Dict(:points => idx)))
    wrenches[idx] = wrench
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
    add_hermite_ride_point!(builder, table, bindings, sam, idx, role, bodies, wrenches)

Add the two instances a beam-anchored point needs and wire them to the two end
bodies of the Timoshenko element it rides: the [`HermiteRidePoint`](@ref) reads both
poses, the [`HermiteRideWrench`](@ref) reads the ride point's `pos`/`vel` and two
moment arms and delivers a share of its load to each body. Returns the ride point,
which is what the incident segments connect to.
"""
function add_hermite_ride_point!(builder, table, bindings, sam, idx, role, bodies,
                                 wrenches)
    joint = sam.sys_struct.timoshenko_joints[role.joint_idx]
    index_map = Dict(:points => idx, :timoshenko_joints => role.joint_idx)
    kinematics = kernel!(builder, table, sam, :hermite_ride_point, idx,
                         params -> HermiteRidePoint(sam, params, idx;
                                                    name = :hermite_ride_point),
                         JOINT_INPUTS, HERMITE_RIDE_OUTPUTS)
    statics = kernel!(builder, table, sam, :hermite_ride_wrench, idx,
                      params -> HermiteRideWrench(sam, params, idx;
                                                  name = :hermite_ride_wrench),
                      HERMITE_WRENCH_INPUTS, JOINT_OUTPUTS)
    ride = add_instance!(builder, kinematics.index)
    wrench = add_instance!(builder, statics.index)
    push!(bindings, (ride, kinematics, index_map))
    push!(bindings, (wrench, statics, index_map))
    wrenches[idx] = wrench
    for (prefix, body) in ((:a, joint.body_a_idx), (:b, joint.body_b_idx))
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
    add_joint!(builder, table, bindings, sam, joint, bodies, container, make)

Add one body-to-body joint and wire it: both end bodies' poses in, the restoring
wrench on each back out. `container` names the joint collection its parameters are
remapped over and `make` builds the component.
"""
function add_joint!(builder, table, bindings, sam, joint, bodies, container, make)
    entry = kernel!(builder, table, sam, container, joint.idx,
                    params -> make(sam, params, joint.idx; name = container),
                    JOINT_INPUTS, JOINT_OUTPUTS)
    instance = add_instance!(builder, entry.index)
    push!(bindings, (instance, entry, Dict(container => joint.idx)))
    for (prefix, body) in ((:a, joint.body_a_idx), (:b, joint.body_b_idx))
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
    add_segment!(builder, table, bindings, sam, idx, role) -> Int

Add the instance for segment `idx` and record its parameter remapping.
"""
function add_segment!(builder, table, bindings, sam, idx, role)
    index_map = Dict(:segments => idx)
    entry = if role.kind === :pulley
        index_map[:pulleys] = role.pulley_idx
        kernel!(builder, table, sam, :pulley_segment, idx,
                params -> PulleySegment(sam, params, idx, role.pulley_idx;
                                        name = :pulley_segment),
                [SEGMENT_INPUTS; :rest_len], [SEGMENT_OUTPUTS; :tension])
    elseif role.kind === :tether
        kernel!(builder, table, sam, :tether_segment, idx,
                params -> TetherSegment(sam, params, idx; name = :tether_segment),
                [SEGMENT_INPUTS; :rest_len],
                [SEGMENT_OUTPUTS; :src_tension; :dst_tension])
    else
        with_drag = role.kind === :spring
        key = with_drag ? :spring_segment : :structural_segment
        kernel!(builder, table, sam, key, idx,
                params -> SpringSegment(sam, params, idx; name = key, with_drag),
                SEGMENT_INPUTS, SEGMENT_OUTPUTS)
    end
    instance = add_instance!(builder, entry.index)
    push!(bindings, (instance, entry, index_map))
    return instance
end

"""
    wire_segment!(builder, sys_struct, idx, role, point_roles, points, segments)

Connect segment `idx` to its two endpoints: their pose into its inputs, its
endpoint force, mass share and drag share back into theirs. A pulley segment also
exchanges the rope split and its tension with the pulley point; a tether segment
takes its rest length from the winch and delivers its tension there.
"""
function wire_segment!(builder, sys_struct, idx, role, point_roles, point_instances,
                       segment_instances, wrench_instances)
    segment = segment_instances[idx]
    source, target = sys_struct.segments[idx].point_idxs
    source == target && error(
        "ScheduledBackend: segment $(sys_struct.segments[idx].name) has both " *
        "endpoints at point $source.")
    for (endpoint, prefix) in ((source, :src), (target, :dst))
        point = point_instances[endpoint]
        connect!(builder, point, :pos, segment, Symbol(prefix, :_pos))
        connect!(builder, point, :vel, segment, Symbol(prefix, :_vel))
        loaded = point_roles[endpoint].kind === :ride ?
            wrench_instances[endpoint] : point
        connect!(builder, segment, :half_drag, loaded, :drag_in)
        if point_roles[endpoint].kind in (:particle, :pulley, :ride, :hermite,
                                          :wing_node)
            connect!(builder, segment, Symbol(prefix, :_force), loaded, :force_in)
            connect!(builder, segment, :half_mass, loaded, :mass_in)
        end
    end
    if role.kind === :pulley
        pulley = sys_struct.pulleys[role.pulley_idx]
        pulley_point = point_instances[pulley_point_index(sys_struct, pulley)]
        connect!(builder, pulley_point, :pulley_len_out, segment, :rest_len)
        connect!(builder, segment, :tension, pulley_point, :tension_in)
    elseif role.kind === :tether
        winch_point = point_instances[role.winch_point]
        winch = sys_struct.winches[find_winch(sys_struct, role.winch_point)]
        position = findfirst(==(role.tether_idx), winch.tether_idxs)
        connect!(builder, winch_point, Symbol(:tether_len_, position), segment,
                 :rest_len)
        source == role.winch_point &&
            connect!(builder, segment, :src_tension, winch_point, :tension_in)
        target == role.winch_point &&
            connect!(builder, segment, :dst_tension, winch_point, :tension_in)
    end
    return nothing
end

"""The index of the winch sitting at point `point_idx`."""
find_winch(sys_struct, point_idx) =
    findfirst(winch -> winch.winch_point_idx == point_idx, sys_struct.winches)

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
    bind_params(system, sys_struct, bindings, segment_roles, segment_instances)

Fill the parameter buffer with every kernel's build-time defaults, then record the
live reader for each parameter that is a struct field read, remapped to its own
instance. Returns `(ScheduledParams, ScheduledParamSync)`.
"""
function bind_params(system, sys_struct, bindings, segment_roles, segment_instances)
    numeric = zeros(SimFloat, system.n_params)
    callables = callable_store(system)
    slots = Int[]
    readers = Any[]
    callable_targets = Any[]
    callable_slots = Int[]
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
            push!(callable_targets, callables[inst.kernel][inst.position])
            push!(callable_slots, slot)
            push!(callable_readers, reader)
        end
    end
    retarget_tether_rest_lengths!(readers, segment_roles)
    sync = ScheduledParamSync(slots, readers, callable_targets, callable_slots,
                              callable_readers)
    return ScheduledParams(numeric, callables), sync
end

"""
    callable_store(system) -> Tuple

The callable parameters, seeded from each kernel's build-time defaults: one entry
per kernel, holding one vector per instance of that kernel. Narrowing each vector
to the defaults' own element type is what keeps a polar call inside a generated
kernel statically dispatched.
"""
function callable_store(system)
    counts = zeros(Int, length(system.kernels))
    for inst in system.instances
        counts[inst.kernel] += 1
    end
    return ntuple(length(system.kernels)) do k
        defaults = identity.(system.kernels[k].callable_defaults)
        [copy(defaults) for _ in 1:counts[k]]
    end
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
        path = (:segments, idx, :l0)
        position = findfirst(reader -> reader isa PathReader && reader.path == path,
                             readers)
        position === nothing && error(
            "ScheduledBackend: segment $idx has no rest-length parameter to retarget")
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
    initial_state(system, sys_struct, point_roles, point_instances)

The initial state vector: each body's principal pose, each particle's
`pos`/`vel`, each pulley's split `pulley_len`/`pulley_vel` and each winch's
`winch_vel` and per-tether lengths, read from the struct.
"""
function initial_state(system, sys_struct, point_roles, point_instances,
                       body_instances)
    u0 = zeros(SimFloat, system.n_states)
    for (idx, instance) in enumerate(body_instances)
        body = sys_struct.bodies[idx]
        # Only a DYNAMIC body integrates; a clamped or fitted one has no state.
        body.type == DYNAMIC || continue
        u0[buffer_slots(system, instance, :states, :com_w)] .= body.com_w
        u0[buffer_slots(system, instance, :states, :com_vel)] .= body.com_vel
        u0[buffer_slots(system, instance, :states, :Q)] .= body.Q_p_to_w
        u0[buffer_slots(system, instance, :states, :omega_p)] .= body.ω_p
    end
    for (idx, role) in enumerate(point_roles)
        instance = point_instances[idx]
        point = sys_struct.points[idx]
        if role.kind in (:particle, :pulley, :wing_node)
            u0[buffer_slots(system, instance, :states, :pos)] .= point.pos_w
            u0[buffer_slots(system, instance, :states, :vel)] .= point.vel_w
        end
        if role.kind === :pulley
            pulley = sys_struct.pulleys[role.pulley_idx]
            u0[only(buffer_slots(system, instance, :states, :pulley_len))] = pulley.len
            u0[only(buffer_slots(system, instance, :states, :pulley_vel))] = pulley.vel
        elseif role.kind === :winch
            winch = sys_struct.winches[role.winch_idx]
            u0[only(buffer_slots(system, instance, :states, :winch_vel))] = winch.vel
            for (k, tether_idx) in enumerate(winch.tether_idxs)
                slot = only(buffer_slots(system, instance, :states,
                                         Symbol(:tether_len_, k)))
                u0[slot] = sys_struct.tethers[tether_idx].len
            end
        end
    end
    return u0
end
