# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# The assembler: the one layer that knows both the `SystemStructure` and the
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
"""The joint fields whose equations differ between a `Real` and a callable, so a
kernel is shared only by joints that agree on them ([`callable_field_key`](@ref))."""
const ELASTIC_RIGIDITIES = [:stiffness_axial, :stiffness_shear,
                            :stiffness_torsion, :stiffness_bending]
const TIMOSHENKO_RIGIDITIES = [:EA, :GA, :GJ, :EIy, :EIz]
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
const TWIST_SURFACE_INPUTS = [:aero_moment_in, :node_moment_in, :node_force_in,
                              :node_mass_in]

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
        if rigid_wing_node(sys_struct, point)
            roles[i] = PointRole(:twist_node, 0, 0, 0, point.wing_idx,
                                 twist_surface_of(sys_struct, point.idx))
            continue
        end
        if point.joint_idx > 0
            roles[i] = PointRole(:hermite, 0, 0, 0, 0, point.joint_idx)
            continue
        end
        if point.type == BODY_STATIC
            point.body_idx > 0 || error(
                "KernelBackend: BODY_STATIC point $(point.name) is anchored to " *
                "neither a body nor a beam.")
            roles[i] = PointRole(:ride, 0, 0, 0, point.body_idx, 0)
            continue
        end
        point.type in (STATIC, DYNAMIC) || error(
            "KernelBackend: point $(point.name) has type $(point.type); only " *
            "STATIC, DYNAMIC and BODY_STATIC are supported so far.")
        wing = kinematic_wing_of(sys_struct, point)
        kind = point.type == STATIC ? :anchor : wing == 0 ? :particle : :wing_node
        roles[i] = PointRole(kind, 0, 0, 0, wing, 0)
    end
    for pulley in sys_struct.pulleys
        pulley.type == DYNAMIC || error(
            "KernelBackend: pulley $(pulley.name) has type $(pulley.type); " *
            "only DYNAMIC pulleys exist.")
        index = pulley_point_index(sys_struct, pulley)
        roles[index] = PointRole(:pulley, pulley.idx, pulley.segment_idxs[1],
                                 0, 0, 0)
    end
    for winch in sys_struct.winches
        index = winch.winch_point_idx
        sys_struct.points[index].type == STATIC || error(
            "KernelBackend: winch $(winch.name) is at a non-STATIC point.")
        roles[index] = PointRole(:winch, 0, 0, winch.idx, 0, 0)
    end
    return roles
end

"""
    rigid_wing_node(sys_struct, point) -> Bool

Whether `point` is a structural node of a `RIGID_DYNAMICS` wing. `point_eqs!` takes
such a node out of the anchor rule entirely — the wing *is* the rigid body, so the
node is placed by a twist-deformed body-frame offset instead of riding anything —
and so does [`classify_points`](@ref).
"""
rigid_wing_node(sys_struct, point) =
    point.is_wing_node && point.wing_idx > 0 &&
    sys_struct.bodies[point.wing_idx].dynamics_type == RIGID_DYNAMICS

"""
    twist_surface_of(sys_struct, idx) -> Int

The one twist surface point `idx` belongs to, or `0`. Two would make its section
twist ambiguous, which `point_eqs!` also rejects.
"""
function twist_surface_of(sys_struct, idx)
    found = [surface.idx for surface in sys_struct.twist_surfaces
             if idx in surface.point_idxs]
    length(found) <= 1 || error(
        "KernelBackend: point $(sys_struct.points[idx].name) is in " *
        "$(length(found)) twist surfaces; expected 0 or 1.")
    return isempty(found) ? 0 : only(found)
end

"""
    kinematic_wing_of(sys_struct, point) -> Int

The fitted (`KINEMATIC`) wing whose frame `point` needs, or `0`. `point_eqs!` damps
*any* DYNAMIC point against its wing's frame, not only aerodynamic surface nodes, so
a steering or pulley point with `body_frame_damping` needs it too.
"""
function kinematic_wing_of(sys_struct, point)
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
    callable = filter(field -> !(getfield(component, field) isa Real), fields)
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
surface's, `0` where there is none. `aero_force_instances` holds each point's
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
    aero_force_instances::Vector{Int}
end

"""
    assemble(sam) -> KernelModel

Translate `sam.sys_struct` into a [`KernelModel`](@ref). Kernels are compiled once per
component type, one instance is added per component, the wiring is recorded, and
the parameters and initial state are bound from the struct afterwards, when every
instance's buffer slice is final.
"""
function assemble(sam; verbose = false)
    sys_struct = sam.sys_struct
    point_roles = classify_points(sys_struct)
    segment_roles = classify_segments(sys_struct)
    builder = SystemBuilder(; verbose)
    table = Dict{Symbol, KernelEntry}()
    bindings = Tuple{Int, KernelEntry, Dict{Symbol, Int}}[]

    body_instances = [add_body!(builder, table, bindings, sam, i)
                      for i in eachindex(sys_struct.bodies)]
    twist_instances = add_twist_surfaces!(builder, table, bindings, sam)
    wrench_instances = zeros(Int, length(sys_struct.points))
    point_instances = [add_point!(builder, table, bindings, sam, i, point_roles[i],
                                  body_instances, wrench_instances, twist_instances)
                       for i in eachindex(sys_struct.points)]
    segment_instances = [add_segment!(builder, table, bindings, sam, i,
                                      segment_roles[i])
                         for i in eachindex(sys_struct.segments)]
    for i in eachindex(sys_struct.segments)
        wire_segment!(builder, sys_struct, i, segment_roles[i],
                      point_instances, segment_instances, wrench_instances)
    end
    for (idx, body) in enumerate(sys_struct.bodies)
        body.type == KINEMATIC || continue
        wire_kinematic_body!(builder, sys_struct, idx, body_instances[idx],
                          point_instances)
    end
    for (idx, role) in enumerate(point_roles)
        role.kind === :wing_node || continue
        connect!(builder, body_instances[role.body_idx], :frame,
                 point_instances[idx], :wing_frame)
        connect!(builder, body_instances[role.body_idx], :vel,
                 point_instances[idx], :wing_velocity)
    end
    for (idx, body) in enumerate(sys_struct.bodies)
        (body.type == DYNAMIC && body.wing_idx != 0) || continue
        connect!(builder, body_instances[body.wing_idx], :frame,
                 body_instances[idx], :wing_frame)
        connect!(builder, body_instances[body.wing_idx], :vel,
                 body_instances[idx], :wing_velocity)
    end
    for joint in sys_struct.elastic_joints
        add_joint!(builder, table, bindings, sam, joint, body_instances,
                   :elastic_joints, ElasticJointComponent, ELASTIC_RIGIDITIES)
    end
    for joint in sys_struct.timoshenko_joints
        add_joint!(builder, table, bindings, sam, joint, body_instances,
                   :timoshenko_joints, TimoshenkoJointComponent,
                   TIMOSHENKO_RIGIDITIES)
    end
    flap_instances = add_flap_deltas!(builder, table, bindings, sam, body_instances)
    inflow_instances = zeros(Int, length(sys_struct.points))
    aero_force_instances = zeros(Int, length(sys_struct.points))
    aero_instances = [add_wing_aero!(builder, table, bindings, sam, wing,
                                     body_instances, point_instances, flap_instances,
                                     twist_instances, wrench_instances,
                                     inflow_instances, aero_force_instances)
                      for wing in sys_struct.wings]

    system = build_system(builder)
    params, sync = bind_params(system, sys_struct, bindings, segment_roles,
                               segment_instances)
    apply_constants!(params, system, segment_roles, segment_instances)
    u0 = zeros(SimFloat, system.n_states)
    initial_state!(u0, system, sys_struct, point_roles, point_instances,
                   body_instances, twist_instances)
    return KernelModel(system, u0, params, sync, point_instances,
                          wrench_instances, segment_instances, body_instances,
                          point_roles, segment_roles, aero_instances,
                          twist_instances, inflow_instances, aero_force_instances)
end

"""
    add_flap_deltas!(builder, table, bindings, sam, bodies) -> Dict{Int, Int}

Add one [`FlapDelta`](@ref) per flapped twist surface and wire its two flap bodies'
orientations in. Returns the instance of each such surface; a surface with no flap
is absent, and the aero input it would feed stays unconnected and so reads the zero
`twist_surface_delta_eqs!` binds it to.
"""
function add_flap_deltas!(builder, table, bindings, sam, bodies)
    instances = Dict{Int, Int}()
    for surface in sam.sys_struct.twist_surfaces
        has_flap(surface) || continue
        entry = kernel!(builder, table, sam, :flap_delta, surface.idx,
                        params -> FlapDelta(sam, params, surface.idx;
                                            name = :flap_delta),
                        FLAP_INPUTS, [:delta])
        instance = add_instance!(builder, entry.index)
        push!(bindings, (instance, entry, Dict(:twist_surfaces => surface.idx)))
        main, flap = surface.flap_body_idxs
        connect!(builder, bodies[main], :frame, instance, :main_frame)
        connect!(builder, bodies[flap], :frame, instance, :flap_frame)
        instances[surface.idx] = instance
    end
    return instances
end

"""
    add_wing_aero!(builder, table, bindings, sam, wing, bodies, points, flaps,
                   twists, wrenches) -> Int

Add a wing's aero instance and wire it, dispatching on how the wing carries its
loads: a `PARTICLE_DYNAMICS` wing's aero delivers a force to each structural point
([`ParticleWingAero`](@ref)), a `RIGID_DYNAMICS` wing's delivers one wrench to its
body ([`add_rigid_wing_aero!`](@ref)).

Either way the aero is a component of its own, not a term inside the points or the
body, because nothing about it is a cycle: neither a point's position nor a body's
pose depends on the force it receives, so the schedule simply runs the structure,
then the wing frame, then the aero, then the derivatives.
"""
function add_wing_aero!(builder, table, bindings, sam, wing, bodies, points, flaps,
                        twists, wrenches, inflow_instances, aero_force_instances)
    wing.dynamics_type == PARTICLE_DYNAMICS || return add_rigid_wing_aero!(
        builder, table, bindings, sam, wing, bodies, twists)
    supports_panel_decomposition(wing.aero) && return add_panel_wing_aero!(
        builder, table, bindings, sam, wing, bodies, points, flaps, wrenches,
        inflow_instances, aero_force_instances)
    sys_struct = sam.sys_struct
    nodes = wing_points(sys_struct, wing)
    surfaces = wing_flap_surfaces(wing)
    key = Symbol(:wing_aero_, wing.idx)
    inputs = [WING_AERO_POSE_INPUTS
              [Symbol(:point_pos_, k) for k in eachindex(nodes)]
              [Symbol(:point_vel_, k) for k in eachindex(nodes)]
              [Symbol(:flap_delta_, surface) for surface in surfaces]]
    outputs = [Symbol(:point_force_, k) for k in eachindex(nodes)]
    entry = kernel!(builder, table, sam, key, wing.idx,
                    params -> ParticleWingAero(sam, params, wing.idx; name = key),
                    inputs, outputs)
    instance = add_instance!(builder, entry.index)
    push!(bindings, (instance, entry, Dict{Symbol, Int}()))
    connect!(builder, bodies[wing.idx], :pos, instance, :wing_pos)
    connect!(builder, bodies[wing.idx], :frame, instance, :wing_frame)
    for (k, node) in enumerate(nodes)
        connect!(builder, points[node.idx], :pos, instance, Symbol(:point_pos_, k))
        connect!(builder, points[node.idx], :vel, instance, Symbol(:point_vel_, k))
        connect!(builder, instance, Symbol(:point_force_, k),
                 load_target(points, wrenches, node.idx), :force_in)
    end
    for surface in surfaces
        haskey(flaps, surface) || continue
        connect!(builder, flaps[surface], :delta, instance,
                 Symbol(:flap_delta_, surface))
    end
    return instance
end

"""
    add_panel_wing_aero!(builder, table, bindings, sam, wing, bodies, points, flaps,
                         wrenches) -> Int

Add a `PARTICLE_DYNAMICS` wing's aerodynamics as small repeated components: one
[`AeroInflowPoint`](@ref) per structural point, one [`AeroInflow`](@ref) per inflow
group, one [`AeroPanel`](@ref) per refined panel, one [`AeroPointForce`](@ref) per
point again, and one [`WingAeroSum`](@ref) for the readouts. Everything between them —
the strut interpolation, the inflow average, the load scatter — is a constant weight,
so it is wiring rather than equations. Returns the sum's instance, which is where the
wing's readouts hang.

This is what [`ParticleWingAero`](@ref) does in one component. One is superlinear in
the wing's size and the other is not, which on a wing of any real size is the whole
difference between a build that finishes and one that does not.
"""
function add_panel_wing_aero!(builder, table, bindings, sam, wing, bodies, points,
                              flaps, wrenches, inflow_instances,
                              aero_force_instances)
    nodes = wing_points(sam.sys_struct, wing)
    body = bodies[wing.idx]
    inflow_points = add_aero_inflow_points!(builder, table, bindings, sam, nodes,
                                            points, body)
    for (node, instance) in zip(nodes, inflow_points)
        inflow_instances[node.idx] = instance
    end
    groups, section_group = aero_inflow_groups(wing.aero, wing, nodes)
    inflows = add_aero_inflows!(builder, table, bindings, sam, groups, inflow_points)
    pitches = if flow_curvature_enabled(wing)
        add_aero_inflows!(builder, table, bindings, sam,
                          aero_pitch_groups(wing.aero, wing, nodes)[1], inflow_points)
    else
        nothing
    end
    wagner = wagner_enabled(wing) ?
        add_wagner_lag!(builder, table, bindings, sam, wing, inflow_points) : nothing
    panels = add_aero_panels!(builder, table, bindings, sam, wing, nodes,
                              inflow_points, inflows, pitches, section_group, flaps,
                              wagner)
    forces = add_aero_point_forces!(builder, table, bindings, sam, wing, nodes,
                                    inflow_points, body, points, wrenches)
    for (node, instance) in zip(nodes, forces)
        aero_force_instances[node.idx] = instance
    end
    for (panel, node, force_weight, couple_weight) in
            aero_scatter_entries(wing.aero, wing, nodes)
        force_weight == 0 || connect!(builder, panels[panel], :force_out,
                                      forces[node], :force_b_in; weight = force_weight)
        couple_weight == 0 || connect!(builder, panels[panel], :couple_out,
                                       forces[node], :force_b_in;
                                       weight = couple_weight)
    end
    return add_wing_aero_sum!(builder, table, bindings, sam, forces)
end

"""
    add_aero_inflow_points!(builder, table, bindings, sam, nodes, points, body)

One [`AeroInflowPoint`](@ref) per structural point of a wing, wired from that point's
kinematics and its wing body's pose. The kernel reads no per-component field, so every
point of every wing shares it. Returns the instances, indexed as `nodes` is.
"""
function add_aero_inflow_points!(builder, table, bindings, sam, nodes, points, body)
    entry = kernel!(builder, table, sam, :aero_inflow_point, 0,
                    params -> AeroInflowPoint(sam, params; name = :aero_inflow_point),
                    AERO_INFLOW_POINT_INPUTS, AERO_INFLOW_POINT_OUTPUTS)
    instances = Int[]
    for node in nodes
        instance = add_instance!(builder, entry.index)
        push!(bindings, (instance, entry, Dict{Symbol, Int}()))
        connect!(builder, points[node.idx], :pos, instance, :pos)
        connect!(builder, points[node.idx], :vel, instance, :vel)
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
    add_aero_panels!(builder, table, bindings, sam, wing, nodes, inflow_points,
                     inflows, pitches, section_group, flaps) -> Vector{Int}

One [`AeroPanel`](@ref) per refined panel, reading its two sections' corners from the
strut interpolation ([`aero_geometry_entries`](@ref)), their inflow from the group they
belong to, and its flap deflection from the twist surface it deflects with. Panels
differ only in their `±1` span sign, so a wing needs at most two kernels.

`pitches` are the [`aero_pitch_groups`](@ref) gathers, or `nothing` on a wing without
[`flow_curvature_enabled`](@ref), which then has no such inputs to connect. `wagner`
is the wing's [`add_wagner_lag!`](@ref) instance, or `nothing` in the same way.
"""
function add_aero_panels!(builder, table, bindings, sam, wing, nodes, inflow_points,
                          inflows, pitches, section_group, flaps, wagner=nothing)
    spanwise = collect(SimFloat, wing.vsm_wing.spanwise_direction)
    with_flap = !isempty(wing_flap_surfaces(wing))
    inputs = isnothing(pitches) ? AERO_PANEL_INPUTS :
        [AERO_PANEL_INPUTS; AERO_PANEL_PITCH_INPUTS]
    isnothing(wagner) || (inputs = [inputs; AERO_PANEL_WAGNER_INPUTS])
    with_flap && (inputs = [inputs; :flap_delta])
    instances = Int[]
    for (panel_idx, orient) in enumerate(panel_span_signs(wing, spanwise))
        key = Symbol(:aero_panel_, wing.idx, orient > 0 ? :_up : :_down)
        entry = kernel!(builder, table, sam, key, panel_idx,
                        params -> AeroPanel(sam, params, wing.idx, panel_idx, orient;
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
    for (panel, corner, node, weight) in aero_geometry_entries(wing.aero, wing, nodes)
        connect!(builder, inflow_points[node], :pos_b, instances[panel], corner; weight)
    end
    with_flap && wire_panel_flaps!(builder, wing.aero, instances, flaps)
    return instances
end

"""
    add_wagner_lag!(builder, table, bindings, sam, wing, inflow_points) -> Int

The wing's one [`WagnerLag`](@ref) instance, its `va_in` gathered at equal weight
over every node of the wing so the lag rides the wing's mean apparent wind. Returns
the instance the panels read their deficiency from.
"""
function add_wagner_lag!(builder, table, bindings, sam, wing, inflow_points)
    key = Symbol(:wagner_lag_, wing.idx)
    entry = kernel!(builder, table, sam, key, 0,
                    params -> WagnerLag(sam, params, wing.idx; name = key),
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

Connect each panel's flap deflection to the [`FlapDelta`](@ref) of the twist surface it
deflects with. A panel mapped to no surface keeps its input unconnected, which reads
zero — the same deflection `flap_delta_inputs` gives it.
"""
function wire_panel_flaps!(builder, mode, panels, flaps)
    for (panel, surface) in enumerate(mode.panel_twist_surface)
        haskey(flaps, surface) || continue
        connect!(builder, flaps[surface], :delta, panels[panel], :flap_delta)
    end
    return nothing
end

"""
    add_aero_point_forces!(builder, table, bindings, sam, wing, nodes, inflow_points,
                           body, points, wrenches) -> Vector{Int}

One [`AeroPointForce`](@ref) per structural point, gathering the panels' scattered
body-frame load and delivering the world force to the point (or to its wrench half,
when the point rides a body).
"""
function add_aero_point_forces!(builder, table, bindings, sam, wing, nodes,
                                inflow_points, body, points, wrenches)
    key = Symbol(:aero_point_force_, wing.idx)
    instances = Int[]
    for (k, node) in enumerate(nodes)
        entry = kernel!(builder, table, sam, key, node.idx,
                        params -> AeroPointForce(sam, params, wing.idx, node.idx;
                                                 name = key),
                        AERO_POINT_FORCE_INPUTS, AERO_POINT_FORCE_OUTPUTS)
        instance = add_instance!(builder, entry.index)
        push!(bindings, (instance, entry, Dict(:points => node.idx)))
        connect!(builder, inflow_points[k], :pos_b, instance, :pos_b)
        connect!(builder, body, :frame, instance, :wing_frame)
        connect!(builder, instance, :force, load_target(points, wrenches, node.idx),
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
    add_rigid_wing_aero!(builder, table, bindings, sam, wing, bodies, twists) -> Int

Add a `RIGID_DYNAMICS` wing's [`WingAero`](@ref) and wire it: the wing body's pose
and each of its twist surfaces' angle and rate in, the world wrench about the body
COM back into the body, and each surface's aerodynamic hinge moment on to that
surface. Returns the instance. `aero_eqs!` drives every surface's hinge moment
except a prescribed one with no aero sections, which has nothing to drive it with.
"""
function add_rigid_wing_aero!(builder, table, bindings, sam, wing, bodies, twists)
    surfaces = wing.twist_surface_idxs
    key = Symbol(:wing_aero_, wing.idx)
    inputs = [RIDE_INPUTS
              [Symbol(:twist_angle_, surface) for surface in surfaces]
              [Symbol(:twist_vel_, surface) for surface in surfaces]]
    outputs = [WRENCH_OUTPUTS
               [Symbol(:twist_moment_, surface) for surface in surfaces]]
    entry = kernel!(builder, table, sam, key, wing.idx,
                    params -> WingAero(sam, params, wing.idx; name = key),
                    inputs, outputs)
    instance = add_instance!(builder, entry.index)
    push!(bindings, (instance, entry, Dict{Symbol, Int}()))
    for (source, target) in ((:pos, :pose_pos), (:frame, :pose_frame),
                             (:com, :pose_com), (:com_velocity, :pose_com_velocity),
                             (:omega_w, :pose_omega))
        connect!(builder, bodies[wing.idx], source, instance, target)
    end
    connect!(builder, instance, :force_out, bodies[wing.idx], :force_in)
    connect!(builder, instance, :moment_out, bodies[wing.idx], :moment_in)
    for surface in surfaces
        twists[surface] == 0 && continue
        connect!(builder, twists[surface], :twist_angle, instance,
                 Symbol(:twist_angle_, surface))
        connect!(builder, twists[surface], :twist_vel, instance,
                 Symbol(:twist_vel_, surface))
        twist_surface_aero_driven(sam.sys_struct.twist_surfaces[surface]) &&
            connect!(builder, instance, Symbol(:twist_moment_, surface),
                     twists[surface], :aero_moment_in)
    end
    return instance
end

"""
    add_body!(builder, table, bindings, sam, idx) -> Int

Add the instance for body `idx`. Each of the three kinds is its own compiled type,
because which one a body is, is topology: a `DYNAMIC` body integrates, a `STATIC`
one is clamped, and a `KINEMATIC` one is fitted from reference points.
"""
function add_body!(builder, table, bindings, sam, idx)
    body = sam.sys_struct.bodies[idx]
    body.type == KINEMATIC && return add_kinematic_body!(builder, table, bindings,
                                                      sam, idx)
    body.type in (DYNAMIC, STATIC) || error(
        "KernelBackend: body $(body.name) has type $(body.type); only DYNAMIC, " *
        "STATIC and KINEMATIC are supported so far.")
    parented = body.type != STATIC && body.wing_idx != 0
    key = body.type == STATIC ? :static_body :
        parented ? :parented_rigid_body : :rigid_body
    make = body.type == STATIC ?
        (params -> StaticBody(sam, params, idx; name = key)) :
        (params -> RigidBody(sam, params, idx; name = key, parented))
    entry = kernel!(builder, table, sam, key, idx, make,
                    parented ? PARENTED_BODY_INPUTS : BODY_INPUTS, BODY_OUTPUTS)
    instance = add_instance!(builder, entry.index)
    push!(bindings, (instance, entry, Dict(:bodies => idx)))
    return instance
end

"""
    add_kinematic_body!(builder, table, bindings, sam, idx) -> Int

Add a `KINEMATIC` wing body. It has no state: its frame is fitted from four
structural reference points and its origin pose read from a fifth, so the whole
component is wiring. Each reference is a weighted blend of real points, delivered by
the weights in the [`Wiring`](@ref).
"""
function add_kinematic_body!(builder, table, bindings, sam, idx)
    entry = kernel!(builder, table, sam, :kinematic_body, idx,
                    params -> KinematicBody(sam, params, idx; name = :kinematic_body),
                    KINEMATIC_INPUTS, BODY_OUTPUTS)
    instance = add_instance!(builder, entry.index)
    push!(bindings, (instance, entry, Dict(:bodies => idx)))
    return instance
end

"""
    wire_kinematic_body!(builder, sys_struct, idx, body_instance, point_instances)

Wire a fitted wing's reference points into it: the four frame references and the
origin's position and velocity, each a weighted blend of point outputs.
"""
function wire_kinematic_body!(builder, sys_struct, idx, body_instance, point_instances)
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
    add_point!(builder, table, bindings, sam, idx, role, bodies, wrenches, twists)
        -> Int

Add the instance for point `idx` and record which container indices its parameters
must be remapped to. A `:ride` point becomes two instances — the kinematics, whose
index is returned so segments wire to it, and the wrench, recorded in `wrenches` and
wired to its body here.
"""
function add_point!(builder, table, bindings, sam, idx, role, bodies, wrenches,
                    twists)
    sys_struct = sam.sys_struct
    index_map = Dict(:points => idx)
    role.kind === :ride && return add_ride_point!(builder, table, bindings, sam, idx,
                                                  role, bodies, wrenches)
    role.kind === :twist_node && return add_twist_node!(builder, table, bindings, sam,
                                                        idx, role, bodies, wrenches,
                                                        twists)
    role.kind === :hermite && return add_hermite_ride_point!(builder, table, bindings,
                                                             sam, idx, role, bodies,
                                                             wrenches)
    entry = if role.kind === :wing_node
        with_damping = sys_struct.points[idx].body_frame_damping !== nothing
        key = with_damping ? :wing_node : :wing_node_undamped
        kernel!(builder, table, sam, key, idx,
                params -> WingNodePoint(sam, params, idx; name = key, with_damping),
                WING_NODE_INPUTS, PARTICLE_OUTPUTS)
    elseif role.kind === :particle
        kernel!(builder, table, sam, :particle, idx,
                params -> Particle(sam, params, idx; name = :particle),
                PARTICLE_INPUTS, PARTICLE_OUTPUTS)
    elseif role.kind === :anchor
        kernel!(builder, table, sam, :anchor, idx,
                params -> Anchor(sam, params, idx; name = :anchor),
                PARTICLE_INPUTS, PARTICLE_OUTPUTS)
    elseif role.kind === :pulley
        index_map[:pulleys] = role.pulley_idx
        index_map[:segments] = role.segment_idx
        kernel!(builder, table, sam, :pulley_point, idx,
                params -> PulleyParticle(sam, params, idx, role.pulley_idx,
                                         role.segment_idx; name = :pulley_point),
                PULLEY_POINT_INPUTS, PULLEY_POINT_OUTPUTS)
    elseif role.kind === :winch
        winch = sys_struct.winches[role.winch_idx]
        index_map[:winches] = role.winch_idx
        key = Symbol(:winch_, role.winch_idx)
        outputs = [PARTICLE_OUTPUTS;
                   [Symbol(:tether_len_, k) for k in eachindex(winch.tether_idxs)];
                   [Symbol(:tether_vel_, k) for k in eachindex(winch.tether_idxs)]]
        kernel!(builder, table, sam, key, idx,
                params -> WinchAnchor(sam, params, winch, idx; name = key),
                WINCH_INPUTS, outputs)
    else
        error("KernelBackend: point $(sys_struct.points[idx].name) has no " *
              "component for role $(role.kind)")
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
    add_twist_surfaces!(builder, table, bindings, sam) -> Vector{Int}

Add one twist instance per twist surface that has a section twist to report: a
[`TwistSurfaceDOF`](@ref) for a `DYNAMIC` surface, whose twist is a state, and a
[`PrescribedTwist`](@ref) for a `STATIC` one, whose twist is a parameter. A
`KINEMATIC` surface has neither and gets `0`; its deflection is a
[`FlapDelta`](@ref) instead.
"""
function add_twist_surfaces!(builder, table, bindings, sam)
    instances = zeros(Int, length(sam.sys_struct.twist_surfaces))
    for surface in sam.sys_struct.twist_surfaces
        surface.type in (DYNAMIC, STATIC) || continue
        dynamic = surface.type == DYNAMIC
        key = dynamic ? :twist_surface : :prescribed_twist
        make = dynamic ? TwistSurfaceDOF : PrescribedTwist
        entry = kernel!(builder, table, sam, key, surface.idx,
                        params -> make(sam, params, surface.idx; name = key),
                        dynamic ? TWIST_SURFACE_INPUTS : [:aero_moment_in],
                        TWIST_OUTPUTS)
        instance = add_instance!(builder, entry.index)
        push!(bindings, (instance, entry, Dict(:twist_surfaces => surface.idx)))
        instances[surface.idx] = instance
    end
    return instances
end

"""
    add_twist_node!(builder, table, bindings, sam, idx, role, bodies, wrenches, twists)

Add the two instances a `RIGID_DYNAMICS` wing's structural node needs and wire them:
the [`TwistNodePoint`](@ref) places it from the wing body's pose and its surface's
twist, the [`TwistNodeWrench`](@ref) sends its load and moment back to the body and
its bridle couple and mass on to the surface. Returns the node, which is what the
incident segments connect to. Only a `DYNAMIC` surface takes the couple; a
prescribed one holds its twist whatever the nodes pull.
"""
function add_twist_node!(builder, table, bindings, sam, idx, role, bodies, wrenches,
                         twists)
    body = bodies[role.body_idx]
    surface = role.joint_idx
    gated = surface > 0 && !sam.sys_struct.bodies[role.body_idx].group_points_moment
    index_map = Dict(:points => idx)
    surface > 0 && (index_map[:twist_surfaces] = surface)
    suffix = surface > 0 ? "" : "_free"
    kinematics = kernel!(builder, table, sam, Symbol(:twist_node, suffix), idx,
                         params -> TwistNodePoint(sam, params, idx;
                                                  name = :twist_node,
                                                  surface_idx = surface),
                         TWIST_NODE_INPUTS, RIDE_OUTPUTS)
    outputs = surface > 0 ?
              [WRENCH_OUTPUTS; :node_force; :node_moment; :node_mass] :
              WRENCH_OUTPUTS
    key = Symbol(:twist_wrench, suffix, gated ? "_ungated" : "")
    statics = kernel!(builder, table, sam, key, idx,
                      params -> TwistNodeWrench(sam, params, idx; name = key,
                                                surface_idx = surface, gated),
                      TWIST_WRENCH_INPUTS, outputs)
    node = add_instance!(builder, kinematics.index)
    wrench = add_instance!(builder, statics.index)
    push!(bindings, (node, kinematics, index_map))
    push!(bindings, (wrench, statics, index_map))
    wrenches[idx] = wrench
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
        if sam.sys_struct.twist_surfaces[surface].type == DYNAMIC
            connect!(builder, wrench, :node_force, twists[surface], :node_force_in)
            connect!(builder, wrench, :node_moment, twists[surface], :node_moment_in)
            connect!(builder, wrench, :node_mass, twists[surface], :node_mass_in)
        end
    end
    return node
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
function add_joint!(builder, table, bindings, sam, joint, bodies, container, make,
                    rigidity_fields)
    key = callable_field_key(container, joint, rigidity_fields)
    entry = kernel!(builder, table, sam, key, joint.idx,
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
    stiffness = (:unit_stiffness,)
    segment = sam.sys_struct.segments[idx]
    entry = if role.kind === :pulley
        index_map[:pulleys] = role.pulley_idx
        kernel!(builder, table, sam,
                callable_field_key(:pulley_segment, segment, stiffness), idx,
                params -> PulleySegment(sam, params, idx, role.pulley_idx;
                                        name = :pulley_segment),
                [SEGMENT_INPUTS; :rest_len; :rest_vel],
                [SEGMENT_OUTPUTS; :tension; :line])
    elseif role.kind === :tether
        kernel!(builder, table, sam,
                callable_field_key(:tether_segment, segment, stiffness), idx,
                params -> TetherSegment(sam, params, idx; name = :tether_segment),
                [SEGMENT_INPUTS; :rest_len; :rest_vel],
                [SEGMENT_OUTPUTS; :src_tension; :dst_tension])
    else
        with_drag = role.kind === :spring
        key = with_drag ? :spring_segment : :structural_segment
        kernel!(builder, table, sam, callable_field_key(key, segment, stiffness),
                idx, params -> SpringSegment(sam, params, idx; name = key, with_drag),
                SEGMENT_INPUTS, SEGMENT_OUTPUTS)
    end
    instance = add_instance!(builder, entry.index)
    push!(bindings, (instance, entry, index_map))
    return instance
end

"""
    wire_segment!(builder, sys_struct, idx, role, points, segments)

Connect segment `idx` to its two endpoints: their pose into its inputs, its
endpoint force, mass share and drag share back into theirs — a clamped endpoint
included, which does not move in response but does report the load it carries. A
pulley segment also exchanges the rope split and its tension with the pulley point;
a tether segment takes its rest length from the winch and delivers its tension
there.
"""
function wire_segment!(builder, sys_struct, idx, role, point_instances,
                       segment_instances, wrench_instances)
    segment = segment_instances[idx]
    source, target = sys_struct.segments[idx].point_idxs
    source == target && error(
        "KernelBackend: segment $(sys_struct.segments[idx].name) has both " *
        "endpoints at point $source.")
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
        pulley = sys_struct.pulleys[role.pulley_idx]
        pulley_point = point_instances[pulley_point_index(sys_struct, pulley)]
        connect!(builder, pulley_point, :pulley_len_out, segment, :rest_len)
        connect!(builder, pulley_point, :pulley_vel_out, segment, :rest_vel)
        connect!(builder, segment, :tension, pulley_point, :tension_in)
        connect!(builder, segment, :line, pulley_point, :line_in)
    elseif role.kind === :tether
        winch_point = point_instances[role.winch_point]
        winch = sys_struct.winches[find_winch(sys_struct, role.winch_point)]
        position = findfirst(==(role.tether_idx), winch.tether_idxs)
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

"""
    bind_params(system, sys_struct, bindings, segment_roles, segment_instances)

Fill the parameter buffer with every kernel's build-time defaults, then record the
live reader for each parameter that is a struct field read, remapped to its own
instance. Returns `(KernelParams, KernelParamSync)`.
"""
function bind_params(system, sys_struct, bindings, segment_roles, segment_instances)
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
    sync = KernelParamSync(group_readers(slots, readers),
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
    initial_state!(u0, system, sys_struct, point_roles, point_instances,
                   body_instances, twist_instances)

Fill `u0` with the initial state: each `DYNAMIC` twist surface's twist and rate,
each body's principal pose, each particle's `pos`/`vel`, each pulley's split
`pulley_len`/`pulley_vel` and each winch's `winch_vel` and per-tether lengths, read
from the struct. Called again by [`KernelInitialSync`](@ref) whenever a problem is
reused, so a slot the struct does not reach has to be zeroed here.
"""
function initial_state!(u0, system, sys_struct, point_roles, point_instances,
                        body_instances, twist_instances)
    fill!(u0, zero(SimFloat))
    for (idx, instance) in enumerate(twist_instances)
        instance == 0 && continue
        surface = sys_struct.twist_surfaces[idx]
        surface.type == DYNAMIC || continue
        u0[only(buffer_slots(system, instance, :states, :free_twist_angle))] =
            surface.twist
        u0[only(buffer_slots(system, instance, :states, :twist_omega))] =
            surface.twist_ω
    end
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
            single = length(winch.tether_idxs) == 1
            u0[only(winch_state_slots(system, instance, :winch_vel,
                                      Symbol("motor₊vel")))] = winch.vel
            for (k, tether_idx) in enumerate(winch.tether_idxs)
                alias = single ? Symbol("motor₊len") : nothing
                slot = only(winch_state_slots(system, instance,
                                              Symbol(:tether_len_, k), alias))
                u0[slot] = sys_struct.tethers[tether_idx].len
            end
        end
    end
    return u0
end
