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

"""
    PointRole(kind, pulley_idx, segment_idx, winch_idx)

How one point is realised: `kind` is `:particle`, `:anchor`, `:pulley` or `:winch`.
A `:pulley` point carries the pulley it splits and one of that pulley's segments
(whose material gives the rope mass); a `:winch` point carries its winch.
"""
struct PointRole
    kind::Symbol
    pulley_idx::Int
    segment_idx::Int
    winch_idx::Int
end

"""
    classify_points(sys_struct) -> Vector{PointRole}

Decide each point's component type. A point that splits a pulley is a pulley
particle, a point carrying a winch is a winch anchor, and the rest follow their
`DynamicsType`. `BODY_STATIC` and `KINEMATIC` points belong to the structural and
aero layers and are rejected here.
"""
function classify_points(sys_struct)
    roles = [PointRole(:none, 0, 0, 0) for _ in sys_struct.points]
    for (i, point) in enumerate(sys_struct.points)
        point.type in (STATIC, DYNAMIC) || error(
            "ScheduledBackend: point $(point.name) has type $(point.type); only " *
            "STATIC and DYNAMIC are supported so far.")
        roles[i] = PointRole(point.type == STATIC ? :anchor : :particle, 0, 0, 0)
    end
    for pulley in sys_struct.pulleys
        pulley.type == DYNAMIC || error(
            "ScheduledBackend: pulley $(pulley.name) has type $(pulley.type); " *
            "only DYNAMIC pulleys exist.")
        index = pulley_point_index(sys_struct, pulley)
        roles[index] = PointRole(:pulley, pulley.idx, pulley.segment_idxs[1], 0)
    end
    for winch in sys_struct.winches
        index = winch.winch_point_idx
        sys_struct.points[index].type == STATIC || error(
            "ScheduledBackend: winch $(winch.name) is at a non-STATIC point.")
        roles[index] = PointRole(:winch, 0, 0, winch.idx)
    end
    return roles
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
all the state getter and the control setter need to find their values.
"""
struct ScheduledModel{S, P}
    system::S
    u0::Vector{SimFloat}
    params::P
    param_sync::ScheduledParamSync
    point_instances::Vector{Int}
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
    point_instances = [add_point!(builder, table, bindings, sam, i, point_roles[i])
                       for i in eachindex(sys_struct.points)]
    segment_instances = [add_segment!(builder, table, bindings, sam, i,
                                      segment_roles[i])
                         for i in eachindex(sys_struct.segments)]
    for i in eachindex(sys_struct.segments)
        wire_segment!(builder, sys_struct, i, segment_roles[i], point_roles,
                      point_instances, segment_instances)
    end

    system = build_system(builder)
    params, sync = bind_params(system, sys_struct, bindings, segment_roles,
                               segment_instances)
    apply_constants!(params, system, segment_roles, segment_instances)
    u0 = initial_state(system, sys_struct, point_roles, point_instances,
                       body_instances)
    return ScheduledModel(system, u0, params, sync, point_instances,
                          segment_instances, body_instances, point_roles,
                          segment_roles)
end

"""
    add_body!(builder, table, bindings, sam, idx) -> Int

Add the instance for body `idx`. A `STATIC` body is clamped, which is topology, so
it gets its own compiled type rather than a runtime branch.
"""
function add_body!(builder, table, bindings, sam, idx)
    body = sam.sys_struct.bodies[idx]
    body.type in (DYNAMIC, STATIC) || error(
        "ScheduledBackend: body $(body.name) has type $(body.type); only DYNAMIC " *
        "and STATIC are supported so far.")
    frozen = body.type == STATIC
    key = frozen ? :static_body : :rigid_body
    entry = kernel!(builder, table, sam, key, idx,
                    params -> RigidBody(sam, params, idx; name = key, frozen),
                    BODY_INPUTS, BODY_OUTPUTS)
    instance = add_instance!(builder, entry.index)
    push!(bindings, (instance, entry, Dict(:bodies => idx)))
    return instance
end

"""
    add_point!(builder, table, bindings, sam, idx, role) -> Int

Add the instance for point `idx` and record which container indices its parameters
must be remapped to.
"""
function add_point!(builder, table, bindings, sam, idx, role)
    sys_struct = sam.sys_struct
    index_map = Dict(:points => idx)
    entry = if role.kind === :particle
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
                       segment_instances)
    segment = segment_instances[idx]
    source, target = sys_struct.segments[idx].point_idxs
    source == target && error(
        "ScheduledBackend: segment $(sys_struct.segments[idx].name) has both " *
        "endpoints at point $source.")
    for (endpoint, prefix) in ((source, :src), (target, :dst))
        point = point_instances[endpoint]
        connect!(builder, point, :pos, segment, Symbol(prefix, :_pos))
        connect!(builder, point, :vel, segment, Symbol(prefix, :_vel))
        connect!(builder, segment, :half_drag, point, :drag_in)
        if point_roles[endpoint].kind in (:particle, :pulley)
            connect!(builder, segment, Symbol(prefix, :_force), point, :force_in)
            connect!(builder, segment, :half_mass, point, :mass_in)
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
        u0[buffer_slots(system, instance, :states, :com_w)] .= body.com_w
        u0[buffer_slots(system, instance, :states, :com_vel)] .= body.com_vel
        u0[buffer_slots(system, instance, :states, :Q)] .= body.Q_p_to_w
        u0[buffer_slots(system, instance, :states, :omega_p)] .= body.ω_p
    end
    for (idx, role) in enumerate(point_roles)
        instance = point_instances[idx]
        point = sys_struct.points[idx]
        if role.kind in (:particle, :pulley)
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
