# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

"""
    SymbolicAWEModelsNetworkDynamicsExt

Network backend: assemble a `SymbolicAWEModel` as a NetworkDynamics `Network`,
compiling one kernel per component *type* (flat in node count). Loaded when both
`NetworkDynamics` and `Graphs` are available. Implements
`build_prob!(::NetworkBackend, sam)` and `init_backend!(::NetworkBackend, …)`.

The per-type `System`s here are the network's assembly of the same physics the
monolith connector components use: they call the shared `point_acceleration` /
`segment_endpoint_loads` helpers (D2). The only difference is the aggregation
sign — ND sums edge outputs *into* the vertex input, so edges emit the **positive**
endpoint loads (the monolith's `Flow` convention negates them).

Kernels use array-valued I/O variables (`pos(t)[1:3]`), which the forked ND MTK
integration scalarizes for indexing to `:pos_1, :pos_2, :pos_3`. Parameters are
declared through the flat-params constructors (`SAM.make_param` /
`SAM.make_array_param`), never raw `@parameters`; each per-instance parameter is
re-read from the `SystemStructure` every step by the same `SAM.ParamGroup`
`sync_params!` the monolith uses, only with a `VIndex`/`EIndex` setter.

## Uniform interface

NetworkDynamics fixes one vertex-output width and one edge-output width for the
whole network. Every vertex outputs `pos[1:3], vel[1:3], pulley_len_out` (plain
points emit `pulley_len_out = 0`) and every edge outputs, per endpoint,
`force[1:3], mass, tension` — `tension` being a role-signed scalar spring force a
pulley reads as `spring[seg1] − spring[seg2]` and a winch reads as its incident
tether tension.

## Coverage

Points (`STATIC`/`DYNAMIC`), spring-damper segments, dynamic pulleys, and reeling
winches are network components. A pulley is a local vertex owning `pulley_len`
whose two incident segments read it as their `l0`. A winch is a local vertex
owning `winch_vel` and one `tether_len` per connected tether; every segment of a
winched tether reads that `tether_len` through a NetworkDynamics **external
input** (`extin`, an exact within-step non-local read) and sets `l0 =
tether_len / n_segments`. `BODY_STATIC`/`KINEMATIC` points, rigid bodies, joints
and VSM aero are still monolith-only.
"""
module SymbolicAWEModelsNetworkDynamicsExt

using SymbolicAWEModels
using NetworkDynamics
using Graphs
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using SymbolicIndexingInterface: setp, getu
using LinearAlgebra: cross, ×

const SAM = SymbolicAWEModels

# ======================= live struct-backed parameters ======================= #
# Per-type ND kernels carry generic parameters (one symbol per field, one value
# per instance in `NWParameter`). To keep the `SystemStructure` the single source
# of truth — so mutating a field takes effect live, as on the monolith — the flat
# parameter buffer is re-read from the struct every step by the monolith's own
# `SAM.ParamGroup`/`sync_group!`, using the same `PathReader` (getproperty) primitive
# and `setp` (setproperty) setter — only the setter indexes `VIndex`/`EIndex` onto the
# `Network` instead of bare symbols. Every synced parameter is a plain field read
# (`PathReader`); wind, tether drag and pulley rope mass are computed *in-equation*
# from such reads (`ground_wind_vec`, `set.cd_tether`, `pulley_rope_mass`). The
# assembly-fixed structural constants (pulley side/sign, segment count, tension sign,
# body-damp fallback) are topology, not struct data, so they are written into the
# parameter vector once by [`set_const_params!`](@ref) and never re-read.

"""
    network_view(ss)

A build-time [`SAM.ParamView`](@ref) tagged with [`SAM.NetworkBackend`](@ref):
component equations read `params.points[idx].field` (etc.), minting one *generic*
symbol per struct field (one kernel per component type) and recording each as a
`SAM.ParamEntry` whose `path` drives per-instance sync ([`replay_fields!`](@ref)).
"""
network_view(ss) = SAM.ParamView{SAM.NetworkBackend}(SAM.ParamRegistry(ss))

"""Alias for [`SAM.param_unknowns`](@ref): a kernel's recorded `System` param list."""
param_unknowns(params) = SAM.param_unknowns(params)

"""Runtime address of a scalarized parameter `sym`: `EIndex` on an edge, else
`VIndex` on a vertex, at network index `index`."""
param_addr(edge::Bool, index, sym) = edge ? EIndex(index, sym) : VIndex(index, sym)

"""
    replay_fields!(builder, reg, container, addr_index, cont_index, ss; edge=false, skip=())

Bind, for one real instance, every plain struct-field parameter a kernel's `reg`
recorded under `container` (`:points`/`:segments`). Each entry's `path` is
re-pointed to `cont_index` for the live reader and addressed at `addr_index`
([`param_addr`](@ref)); vector fields expand to one scalar `field_k` per component and
matrix fields to `field_row_col` (column-major), matching NetworkDynamics' parameter
scalarization. Computed entries (non-`PathReader`,
e.g. `wind_gnd`) and any field in `skip` are left for the caller to bind explicitly.
"""
function replay_fields!(builder, reg, container::Symbol, addr_index, cont_index, ss;
                        edge=false, skip=())
    for entry in reg.entries
        entry.read isa SAM.PathReader || continue
        path = entry.read.path
        (length(path) >= 3 && path[1] === container) || continue
        field = path[3]
        field in skip && continue
        if entry.kind === :scalar
            add_param!(builder, param_addr(edge, addr_index, field),
                       SAM.PathReader((container, cont_index, field)))
        elseif entry.kind === :array
            value = entry.read(ss)
            if ndims(value) == 2
                rows, cols = size(value)
                for col in 1:cols, row in 1:rows
                    add_param!(builder,
                        param_addr(edge, addr_index, Symbol(field, :_, row, :_, col)),
                        SAM.PathReader((container, cont_index, field,
                                        (col - 1) * rows + row)))
                end
            else
                for k in 1:length(value)
                    add_param!(builder,
                        param_addr(edge, addr_index, Symbol(field, :_, k)),
                        SAM.PathReader((container, cont_index, field, k)))
                end
            end
        end
    end
    return nothing
end

"""
    record_edge_callables!(callables, reg, container, addr_index, cont_index, ss)

Record every callable (nonnumeric) struct-field parameter a joint-edge kernel's `reg`
recorded under `container` (`:timoshenko_joints`/`:elastic_joints`) — the nonlinear
rigidity/stiffness laws (`EIy(κ)` etc.). Each is bound at `EIndex(addr_index, field)` with
a live reader on `(container, cont_index, field)`, routed through the fork's nonnumeric
SII path ([`add_callable!`](@ref)). Numeric fields are handled by [`replay_fields!`](@ref).
"""
function record_edge_callables!(callables, reg, container::Symbol, addr_index,
                                cont_index, ss)
    for entry in reg.entries
        entry.read isa SAM.PathReader || continue
        entry.kind === :callable || continue
        path = entry.read.path
        (length(path) >= 3 && path[1] === container) || continue
        field = path[3]
        add_callable!(callables, EIndex(addr_index, field),
                      SAM.PathReader((container, cont_index, field)))
    end
    return nothing
end

"""
    ParamBuilder

Accumulates `(index, reader)` pairs while the network is assembled, then builds a
`SAM.ParamGroup`. `index` is a `VIndex`/`EIndex` onto a scalar parameter symbol;
`reader(sys_struct)` returns its live value.
"""
struct ParamBuilder
    indices::Vector{Any}
    readers::Vector{Any}
end
ParamBuilder() = ParamBuilder(Any[], Any[])

"""
    add_param!(builder, index, reader)

Record one live parameter: written to `index`, read as `reader(sys_struct)`.
"""
function add_param!(builder::ParamBuilder, index, reader)
    push!(builder.indices, index)
    push!(builder.readers, reader)
    return nothing
end

"""
    CallableBuilder

Accumulates callable (nonnumeric) parameters — the `ContinuousPolar` cl/cd/cm the live
aero modes carry. Each is set through the ND fork's uniform SII callable route
(`setp(nw, VIndex(i, sym))`), one setter per parameter because a batched `setp` over a
vector does not route element-wise through the nonnumeric path. Distinct from
[`ParamBuilder`](@ref), whose values live in the flat `Float64` buffer.
"""
struct CallableBuilder
    indices::Vector{Any}
    readers::Vector{Any}
end
CallableBuilder() = CallableBuilder(Any[], Any[])

"""
    add_callable!(builder, index, reader)

Record one callable parameter: written to `index` (a `VIndex` onto a nonnumeric
symbol), read as `reader(sys_struct)` (the live polar object).
"""
function add_callable!(builder::CallableBuilder, index, reader)
    push!(builder.indices, index)
    push!(builder.readers, reader)
    return nothing
end

"""
    MultiCallableSetter(setters)

A [`SAM.ParamGroup`](@ref) setter over several callable parameters: applies each
single-index `setp(nw, VIndex(...))` (which routes through the ND fork's nonnumeric
store) to the matching buffer entry. Mirrors the batched-`setp` interface
`(target, buffer)` the shared `sync_group!` calls, but element-wise so each callable
takes the nonnumeric path.
"""
struct MultiCallableSetter{S}
    setters::S
end
function (m::MultiCallableSetter)(target, buffer)
    @inbounds for k in eachindex(m.setters)
        m.setters[k](target, buffer[k])
    end
    return nothing
end

"""
    build_network_param_sync(nw, builder, callables)

Turn the recorded parameters into a `SAM.ParamSync` — a `Float64` scalar group (all
scalar + scalarized-array-element params, via a batched `setp`) and a callable group
(the polar callables, via a [`MultiCallableSetter`](@ref)). Either may be `nothing`;
returns `nothing` when both are empty. It is the same `ParamSync` the monolith syncs,
only with `VIndex`/`EIndex` setters, applied through the shared `ProbWithAttributes`
machinery every step (and after every VSM refresh).
"""
function build_network_param_sync(nw, builder::ParamBuilder, callables::CallableBuilder)
    scalar = nothing
    if !isempty(builder.indices)
        setter = setp(nw, builder.indices)
        buffer = Vector{SAM.SimFloat}(undef, length(builder.indices))
        scalar = SAM.ParamGroup(setter, builder.readers, buffer)
    end
    callable = nothing
    if !isempty(callables.indices)
        setters = [setp(nw, idx) for idx in callables.indices]
        buffer = Vector{Any}(undef, length(callables.indices))
        callable = SAM.ParamGroup(MultiCallableSetter(setters), callables.readers, buffer)
    end
    (scalar === nothing && callable === nothing) && return nothing
    return SAM.ParamSync(scalar, nothing, callable)
end

# ======================= vertex Systems ======================= #

"""
    network_dynamic_point(s, params, idx; name)

Particle vertex: the shared [`SAM.DynamicPoint`](@ref) component wrapped as a network
vertex. Both backends assemble the identical component; the network wraps it as a
`VertexModel`, the monolith `@named`-instantiates it.
"""
network_dynamic_point(s, params, idx; name, wide = false) =
    SAM.DynamicPoint(s, params, idx; name, wide)

"""
    network_static_point(s, params, idx; name)

Ground-anchored vertex: the shared [`SAM.StaticPoint`](@ref) component wrapped as a
network vertex.
"""
network_static_point(s, params, idx; name, wide = false) =
    SAM.StaticPoint(s, params, idx; name, wide)

"""
    pulley_rope_mass(params, idx)

The rope mass `sum_len · ρ · π (d/2)²` of the pulley at vertex `idx`, built
*in-equation* from the struct-field parameters `params.pulleys[…].sum_len` and its
first segment's `density`/`diameter` (a `pulley_mass` reader is thus unneeded). The
mass is a per-pulley constant but varies between pulleys, so it is read live rather
than baked into the shared kernel default.
"""
function pulley_rope_mass(params, idx)
    ss = params.reg.sys_struct
    pulley_idx = findfirst(p -> pulley_point_idx(ss, p) == idx, ss.pulleys)
    seg = params.segments[ss.pulleys[pulley_idx].segment_idxs[1]]
    return params.pulleys[pulley_idx].sum_len * seg.density *
           π * (seg.diameter / 2)^2
end

"""
    network_pulley_point(s, params, idx; name)

Dynamic pulley vertex: the shared [`SAM.PulleyPoint`](@ref) wrapped as a network
vertex, with the rope mass supplied by [`pulley_rope_mass`](@ref). The aggregated
`tension_in` is the imbalance `spring[seg1] − spring[seg2]` between its two incident
segments; `pulley_len_out` exposes the split so the incident segments read it as `l0`.
"""
network_pulley_point(s, params, idx; name, wide = false) =
    SAM.PulleyPoint(s, params, idx, pulley_rope_mass(params, idx); name, wide)

"""
    BODY_DAMP_EXTIN

The 15 external-input symbols a wing node reads: the four ref-point positions
`zp1/zp2/yp1/yp2` (12) and the wing origin velocity `ovel` (3). Matches the input
variable names of [`SAM.wing_node_inputs`](@ref); used to wire the NetworkDynamics
`extin` for each wing node.
"""
const BODY_DAMP_EXTIN = [
    :zp1x, :zp1y, :zp1z, :zp2x, :zp2y, :zp2z,
    :yp1x, :yp1y, :yp1z, :yp2x, :yp2y, :yp2z, :ovx, :ovy, :ovz]

"""
    network_wing_node_point(s, params, idx; name)

A KINEMATIC wing node: the shared [`SAM.WingNodePoint`](@ref) wrapped as a network
vertex. The wing frame and wing velocity it reads through [`SAM.wing_node_inputs`](@ref)
are supplied via NetworkDynamics `extin` from the wing's ref points ([`BODY_DAMP_EXTIN`]).
"""
network_wing_node_point(s, params, idx; name) =
    SAM.WingNodePoint(s, params, idx; name)

"""
    network_wing_node_pulley_point(s, params, idx; name)

A pulley vertex that is also a wing node: the shared [`SAM.WingNodePulleyPoint`](@ref)
wrapped as a network vertex, with the rope mass from [`pulley_rope_mass`](@ref).
"""
network_wing_node_pulley_point(s, params, idx; name) =
    SAM.WingNodePulleyPoint(s, params, idx, pulley_rope_mass(params, idx); name)

"""
    network_winch_point(s, winch, winch_point; name)

Reeling winch vertex: the shared [`SAM.WinchPoint`](@ref) wrapped as a network vertex.
Each `tether_len_k` state is read by that tether's segments through an `extin`.
"""
network_winch_point(s, winch, winch_point; name) =
    SAM.WinchPoint(s, winch, winch_point; name)

# ======================= edge Systems ======================= #

"""
    network_segment(s, params, idx; name)

Plain/structural spring-damper edge: the shared [`SAM.SpringDamperSegment`](@ref)
component wrapped as a network edge. Emits zero tension (neither endpoint is a
pulley/winch reader); a [`SAM.wing_structural_segment`](@ref) representative drops the
drag term, giving the `:structural` edge kind its own drag-free compiled kernel.
"""
network_segment(s, params, idx; name) = SAM.SpringDamperSegment(s, params, idx; name)

"""
    network_pulley_segment(s, params, idx; name)

Pulley spring-damper edge. `l0` is driven by the pulley vertex's `pulley_len`
output (read from whichever endpoint is the pulley, `pulley_at_src`): the first
pulley segment gets `l0 = pulley_len`, the second `l0 = sum_len − pulley_len`
(`pulley_side = ±1`). It emits its role-signed spring tension `pulley_side·spring`
to the pulley endpoint so the pulley aggregates `spring[seg1] − spring[seg2]`.
"""
function network_pulley_segment(s, params, idx; name, wide = false)
    io, extras = wide ? SAM.segment_io_wide() : (SAM.segment_io(), nothing)
    vars, _, _, src_pulley_len, _, _, dst_pulley_len = io
    spring, wind = SAM.segment_spring_params(params, idx)
    pulley_sum_len = SAM.make_param(:pulley_sum_len, 1.0)
    pulley_side = SAM.make_param(:pulley_side, 1.0)
    pulley_at_src = SAM.make_param(:pulley_at_src, 0.0)
    endpoint_len = ifelse(pulley_at_src > 0.5, src_pulley_len, dst_pulley_len)
    l0 = ifelse(pulley_side > 0.0, endpoint_len, pulley_sum_len - endpoint_len)
    fsrc, fdst, half, spring_scalar = SAM.segment_loads(s, io, l0, spring, wind)
    src_tension_val = ifelse(pulley_at_src > 0.5, pulley_side * spring_scalar, 0.0)
    dst_tension_val = ifelse(pulley_at_src > 0.5, 0.0, pulley_side * spring_scalar)
    eqs = SAM.endpoint_load_eqs(io, fsrc, fdst, half, src_tension_val, dst_tension_val)
    if wide
        eqs = [eqs; collect(extras.src_moment) .~ 0; collect(extras.dst_moment) .~ 0]
    end
    return System(eqs, t, vars,
        [param_unknowns(params); pulley_sum_len; pulley_side; pulley_at_src]; name)
end

"""
    network_wrench_segment(s, params, idx, ride_idx, body_at_src; name)

Wide spring-damper edge between a free point and a rigid body's `BODY_STATIC` ride point
(absorbed into the body vertex). It reads the body endpoint's pose (`pose_R`/`pose_com`/
`pose_com_vel`/`pose_omega`), reconstructs the ride position `body_pos + R·anchor_b` and
its rigid-motion velocity, computes the spring-damper load between the free point and
that ride position (shared [`SAM.segment_endpoint_loads`](@ref)), and delivers the
**wrench** to the body — the force at the anchor plus the moment `arm × force` about the
body COM — while the free point receives the equal-and-opposite force. `anchor_b` is read
live from the ride point's struct field; `body_at_src` (which endpoint is the body) is a
build-time `Bool` **baked into the kernel** — no runtime `ifelse`, so the symbolic
expressions stay small (an `ifelse`-selected pose blows up `mtkcompile`).
"""
function network_wrench_segment(s, params, idx, ride_idx, body_at_src::Bool; name)
    io, extras = SAM.segment_io_wide()
    vars = io[1]
    (_, src_pos, src_vel, _, dst_pos, dst_vel, _,
     src_force, src_mass, src_tension, dst_force, dst_mass, dst_tension) = io
    anchor_b = collect(params.points[ride_idx].anchor_b)
    body_pos, body_R, body_com, body_com_vel, body_omega, free_pos, free_vel =
        body_at_src ?
            (src_pos, extras.src_pose_R, extras.src_pose_com, extras.src_pose_com_vel,
             extras.src_pose_omega, dst_pos, dst_vel) :
            (dst_pos, extras.dst_pose_R, extras.dst_pose_com, extras.dst_pose_com_vel,
             extras.dst_pose_omega, src_pos, src_vel)
    R_mat = reshape(collect(body_R), 3, 3)
    anchor_w = collect(body_pos) .+ R_mat * anchor_b
    arm = anchor_w .- collect(body_com)
    ride_vel = collect(body_com_vel) .+ (collect(body_omega) × arm)
    spring, wind = SAM.segment_spring_params(params, idx; with_drag = true)
    l0 = params.segments[idx].l0
    force_on_free, force_on_ride, half_mass, _ = SAM.segment_endpoint_loads(
        s, collect(free_pos), collect(free_vel), anchor_w, ride_vel,
        spring.unit_stiffness, spring.unit_damping, spring.compression_frac, l0,
        spring.diameter, spring.density, spring.cd_tether, collect(wind);
        with_drag = true)
    moment_on_body = arm × force_on_ride
    body_force, free_force = force_on_ride, force_on_free
    src_f, dst_f = body_at_src ? (body_force, free_force) : (free_force, body_force)
    src_m = body_at_src ? moment_on_body : zeros(3)
    dst_m = body_at_src ? zeros(3) : moment_on_body
    eqs = [
        collect(src_force) .~ src_f
        src_mass ~ (body_at_src ? 0.0 : half_mass); src_tension ~ 0.0
        collect(dst_force) .~ dst_f
        dst_mass ~ (body_at_src ? half_mass : 0.0); dst_tension ~ 0.0
        collect(extras.src_moment) .~ src_m
        collect(extras.dst_moment) .~ dst_m
    ]
    return System(eqs, t, vars, param_unknowns(params); name)
end

"""
    joint_edge_ab_poses(io, extras, a_at_src)

The `(a_pose, b_pose)` tuples a body↔body joint edge reads, each `(pos, pose_R, pose_com,
pose_com_vel, pose_omega)`, mapping joint body A/B to the edge `src`/`dst` endpoints by the
build-time `a_at_src` flag.
"""
function joint_edge_ab_poses(io, extras, a_at_src)
    (_, src_pos, _, _, dst_pos, _, _) = io
    src = (src_pos, extras.src_pose_R, extras.src_pose_com, extras.src_pose_com_vel,
           extras.src_pose_omega)
    dst = (dst_pos, extras.dst_pose_R, extras.dst_pose_com, extras.dst_pose_com_vel,
           extras.dst_pose_omega)
    return a_at_src ? (src, dst) : (dst, src)
end

"""
    joint_edge_emit_eqs(io, extras, a_at_src, ex)

The output equations a body↔body joint edge writes: the shared-helper wrench `ex`
(`force_on_a`/`moment_on_a`/`force_on_b`/`moment_on_b`) mapped to the edge's `src`/`dst`
force+moment outputs (`mass`/`tension` zero), following the `a_at_src` orientation.
"""
function joint_edge_emit_eqs(io, extras, a_at_src, ex)
    (_, _, _, _, _, _, _, src_force, src_mass, src_tension,
     dst_force, dst_mass, dst_tension) = io
    src_f, src_m, dst_f, dst_m = a_at_src ?
        (ex.force_on_a, ex.moment_on_a, ex.force_on_b, ex.moment_on_b) :
        (ex.force_on_b, ex.moment_on_b, ex.force_on_a, ex.moment_on_a)
    return [
        collect(src_force) .~ src_f
        src_mass ~ 0.0; src_tension ~ 0.0
        collect(dst_force) .~ dst_f
        dst_mass ~ 0.0; dst_tension ~ 0.0
        collect(extras.src_moment) .~ src_m
        collect(extras.dst_moment) .~ dst_m
    ]
end

"""
    network_timoshenko_joint_edge(s, params, jidx, a_at_src; name)

Wide body↔body edge for a `TimoshenkoJoint` (a corotational beam element). It reads both
bodies' poses, evaluates the shared [`SAM.timoshenko_element_wrench`](@ref) (frame,
per-node deformations, consistent Timoshenko stiffness + damping) and delivers the equal-
and-opposite restoring wrench to each body. `a_at_src` (whether joint body A is the edge's
`src`) is baked in at build time. The `frame`/`theta`/`force`/`moment` torn variables are
edge-internal (mirroring the monolith's tearing). Only constant (`Real`) rigidities are
supported so far.
"""
function network_timoshenko_joint_edge(s, params, jidx, a_at_src::Bool; name)
    io, extras = SAM.segment_io_wide()
    joint = params.reg.sys_struct.timoshenko_joints[jidx]
    torn = @variables begin
        tj_frame(t)[1:3, 1:3]
        tj_theta_a(t)[1:3]
        tj_theta_b(t)[1:3]
        tj_force_a(t)[1:3]
        tj_force_b(t)[1:3]
        tj_moment_a(t)[1:3]
        tj_moment_b(t)[1:3]
    end
    (apos, aR, acom, acomvel, aomega), (bpos, bR, bcom, bcomvel, bomega) =
        joint_edge_ab_poses(io, extras, a_at_src)
    ex = SAM.timoshenko_element_wrench(joint, params;
        frame = torn[1], theta_a = torn[2], theta_b = torn[3],
        force_a = torn[4], force_b = torn[5], moment_a = torn[6], moment_b = torn[7],
        pos_a = apos, R_a = reshape(collect(aR), 3, 3), com_a = acom,
        com_vel_a = acomvel, omega_a_w = aomega,
        pos_b = bpos, R_b = reshape(collect(bR), 3, 3), com_b = bcom,
        com_vel_b = bcomvel, omega_b_w = bomega)
    eqs = [ex.tear_eqs; joint_edge_emit_eqs(io, extras, a_at_src, ex)]
    return System(eqs, t, [io[1]; torn], param_unknowns(params); name)
end

"""
    network_elastic_joint_edge(s, params, jidx, a_at_src; name)

Wide body↔body edge for a lumped 6-DOF `ElasticJoint`. It reads both bodies' poses,
evaluates the shared [`SAM.elastic_joint_wrench`](@ref) (per-DOF axial/shear/torsion/
bending stiffness + damping from the relative anchor pose) and delivers the equal-and-
opposite restoring wrench to each body. `a_at_src` is baked in at build time; the
world-frame force/torque torn variables are edge-internal. Only constant (`Real`)
stiffnesses are supported so far.
"""
function network_elastic_joint_edge(s, params, jidx, a_at_src::Bool; name)
    io, extras = SAM.segment_io_wide()
    joint = params.reg.sys_struct.elastic_joints[jidx]
    torn = @variables ej_force(t)[1:3] ej_torque(t)[1:3]
    (apos, aR, acom, acomvel, aomega), (bpos, bR, bcom, bcomvel, bomega) =
        joint_edge_ab_poses(io, extras, a_at_src)
    ex = SAM.elastic_joint_wrench(joint, params;
        force_w = torn[1], torque_w = torn[2],
        pos_a = apos, R_a = reshape(collect(aR), 3, 3), com_a = acom,
        com_vel_a = acomvel, omega_a_w = aomega,
        pos_b = bpos, R_b = reshape(collect(bR), 3, 3), com_b = bcom,
        com_vel_b = bcomvel, omega_b_w = bomega)
    eqs = [ex.tear_eqs; joint_edge_emit_eqs(io, extras, a_at_src, ex)]
    return System(eqs, t, [io[1]; torn], param_unknowns(params); name)
end

"""
    network_tether_segment(s, params, idx; name)

Winched-tether spring-damper edge. Its rest length is `l0 = tether_len / n_segs`,
where `tether_len` arrives as a NetworkDynamics external input (`tether_len_ext`)
read from the winch vertex. The tether segment incident to the winch point emits
`+spring` there (`tension_sign_src`/`tension_sign_dst = 1`) so the winch reads the
tension; every other tether segment emits zero.
"""
function network_tether_segment(s, params, idx; name, wide = false)
    io, extras = wide ? SAM.segment_io_wide() : (SAM.segment_io(), nothing)
    vars = io[1]
    tether_len_ext = only(@variables tether_len_ext(t) [input = true])
    push!(vars, tether_len_ext)
    spring, wind = SAM.segment_spring_params(params, idx)
    n_segs = SAM.make_param(:n_segs, 1.0)
    tension_sign_src = SAM.make_param(:tension_sign_src, 0.0)
    tension_sign_dst = SAM.make_param(:tension_sign_dst, 0.0)
    l0 = tether_len_ext / n_segs
    fsrc, fdst, half, spring_scalar = SAM.segment_loads(s, io, l0, spring, wind)
    eqs = SAM.endpoint_load_eqs(io, fsrc, fdst, half,
        tension_sign_src * spring_scalar, tension_sign_dst * spring_scalar)
    if wide
        eqs = [eqs; collect(extras.src_moment) .~ 0; collect(extras.dst_moment) .~ 0]
    end
    return System(eqs, t, vars,
        [param_unknowns(params); n_segs; tension_sign_src; tension_sign_dst]; name)
end

# ======================= I/O symbol lists ======================= #
# Array-valued I/O symbols; the forked ND MTK integration scalarizes each to its
# `_1/_2/_3` components while preserving the vertex-output / edge-output ordering.

const VERTEX_INPUTS = [:force_in, :mass_in, :tension_in]
const VERTEX_OUTPUTS = [:pos, :vel, :pulley_len_out]
# Wide superset interface (body-containing networks, §8.5): every vertex shares this
# input/output width; a point zero-fills the pose slots, a body zero-fills pulley_len.
const WIDE_VERTEX_INPUTS = [:force_in, :mass_in, :tension_in, :moment_in]
const WIDE_VERTEX_OUTPUTS =
    [:pos, :vel, :pulley_len_out, :pose_R, :pose_com, :pose_com_vel, :pose_omega]
const WIDE_EDGE_SRC_IN = [:src_pos, :src_vel, :src_pulley_len, :src_pose_R,
    :src_pose_com, :src_pose_com_vel, :src_pose_omega]
const WIDE_EDGE_DST_IN = [:dst_pos, :dst_vel, :dst_pulley_len, :dst_pose_R,
    :dst_pose_com, :dst_pose_com_vel, :dst_pose_omega]
const WIDE_EDGE_SRC_OUT = [:src_force, :src_mass, :src_tension, :src_moment]
const WIDE_EDGE_DST_OUT = [:dst_force, :dst_mass, :dst_tension, :dst_moment]
const EDGE_SRC_IN = [:src_pos, :src_vel, :src_pulley_len]
const EDGE_DST_IN = [:dst_pos, :dst_vel, :dst_pulley_len]
const EDGE_SRC_OUT = [:src_force, :src_mass, :src_tension]
const EDGE_DST_OUT = [:dst_force, :dst_mass, :dst_tension]

# ======================= helpers ======================= #

"""
    segment_stiffness(seg)

Frozen linear stiffness of `seg`. A callable (nonlinear) `unit_stiffness` force
law cannot be passed as a scalar ND parameter, so it is rejected until the
network segment supports it.
"""
function segment_stiffness(seg)
    seg.unit_stiffness isa Real && return Float64(seg.unit_stiffness)
    error("NetworkBackend: segment $(seg.name) has a callable (nonlinear) " *
          "unit_stiffness; only linear (Real) stiffness is supported so far.")
end

"""
    SegmentRoles

Per-segment classification the network build derives once: `kind` is `:plain`,
`:structural` (drag-free wing link), `:pulley` or `:tether`; the remaining fields
carry the pulley split / winched-tether data an edge needs (`0`/`nothing` when not
applicable).
"""
struct SegmentRoles
    kind::Symbol
    pulley_idx::Int
    pulley_side::Float64
    pulley_at_src::Bool
    tether_idx::Int
    n_segs::Int
    tension_sign_src::Float64
    tension_sign_dst::Float64
end

"""
    classify_segments(ss)

Return a `SegmentRoles` per segment. A segment is `:pulley` if it is one of a
pulley's two segments, `:tether` if it belongs to a winched tether, `:structural`
if both endpoints are wing nodes (drag-free), else `:plain`. The tether segment
incident to the winch point is marked to emit `+spring` there.
"""
function classify_segments(ss)
    winch_of_tether = Dict{Int, Int}()
    for winch in ss.winches, tether_idx in winch.tether_idxs
        winch_of_tether[tether_idx] = winch.winch_point_idx
    end
    tether_of_segment = Dict{Int, Int}()
    for tether in ss.tethers, seg_idx in tether.segment_idxs
        tether_of_segment[seg_idx] = tether.idx
    end
    pulley_of_segment = Dict{Int, Tuple{Int, Float64}}()
    for pulley in ss.pulleys
        pulley_of_segment[pulley.segment_idxs[1]] = (pulley.idx, 1.0)
        pulley_of_segment[pulley.segment_idxs[2]] = (pulley.idx, -1.0)
    end

    roles = SegmentRoles[]
    for seg in ss.segments
        # NetworkDynamics orients each SimpleGraph edge min→max, which may differ
        # from `seg.point_idxs`; the pulley/winch signs must follow the edge.
        src_point, dst_point = minmax(seg.point_idxs[1], seg.point_idxs[2])
        if haskey(pulley_of_segment, seg.idx)
            pidx, side = pulley_of_segment[seg.idx]
            pulley_point = pulley_point_idx(ss, ss.pulleys[pidx])
            at_src = (src_point == pulley_point)
            push!(roles, SegmentRoles(:pulley, pidx, side, at_src, 0, 0, 0.0, 0.0))
        elseif haskey(tether_of_segment, seg.idx) &&
               haskey(winch_of_tether, tether_of_segment[seg.idx])
            tidx = tether_of_segment[seg.idx]
            n_segs = length(ss.tethers[tidx].segment_idxs)
            winch_point = winch_of_tether[tidx]
            sign_src = (src_point == winch_point) ? 1.0 : 0.0
            sign_dst = (dst_point == winch_point) ? 1.0 : 0.0
            push!(roles, SegmentRoles(:tether, 0, 0.0, false, tidx, n_segs,
                                      sign_src, sign_dst))
        else
            kind = SAM.wing_structural_segment(ss, seg.idx) ? :structural : :plain
            push!(roles, SegmentRoles(kind, 0, 0.0, false, 0, 0, 0.0, 0.0))
        end
    end
    return roles
end

"""
    pulley_point_idx(ss, pulley)

The pulley's vertex point: the single point shared by its two segments.
"""
function pulley_point_idx(ss, pulley)
    a = ss.segments[pulley.segment_idxs[1]].point_idxs
    b = ss.segments[pulley.segment_idxs[2]].point_idxs
    shared = intersect(a, b)
    length(shared) == 1 || error(
        "NetworkBackend: pulley $(pulley.name) segments share $(length(shared)) " *
        "points; expected exactly 1.")
    return only(shared)
end

"""
    ref_single_id(ref_points)

The single point id of a `WeightedRefPoints`; errors on multi-point weighted refs
(the network wing frame supports single-point refs so far).
"""
function ref_single_id(ref_points)
    length(ref_points.ids) == 1 || error(
        "NetworkBackend: wing ref points must be single-point (got " *
        "$(length(ref_points.ids))-point weighted ref); not supported yet.")
    return ref_points.ids[1]
end

"""
    network_wing_node(ss, point)

The point's `KINEMATIC` wing body if the point needs that wing's frame through
`extin`, else `nothing`. A DYNAMIC point needs the wing frame when it is a wing node
(`is_wing_node` — carries the frozen aero force) **or** it has a nonzero
`body_frame_damping` (`point_eqs.jl`'s `point_damping_accel` damps *any* DYNAMIC
point relative to its wing frame, not only wing nodes — e.g. steering/pulley points).
Selects the wing-node vertex kernel.
"""
function network_wing_node(ss, point)
    (point.type == SAM.DYNAMIC && 0 < point.wing_idx <= length(ss.bodies)) ||
        return nothing
    ss.bodies[point.wing_idx].type == SAM.KINEMATIC || return nothing
    bd = point.body_frame_damping
    needs_damp = bd !== nothing && !all(iszero, bd)
    (point.is_wing_node || needs_damp) || return nothing
    return ss.bodies[point.wing_idx]
end

"""
    point_body_damp(ss, point)

The point's `body_frame_damping` coefficients if it is a wing node
([`network_wing_node`](@ref)) with a nonzero coefficient, else `nothing`. Used to
decide whether the `body_damp` parameter reads the struct field or is held at zero.
"""
function point_body_damp(ss, point)
    network_wing_node(ss, point) === nothing && return nothing
    bd = point.body_frame_damping
    (bd === nothing || all(iszero, bd)) && return nothing
    return bd
end

"""
    body_damp_extin(ss, wing)

The `extin` pair list binding a body-damped point's ref-point inputs
([`BODY_DAMP_EXTIN`]) to the wing's structural ref points (z/y frame points'
positions and the origin velocity), read as the neighbours' scalarized `pos_k`/
`vel_k` outputs.
"""
function body_damp_extin(ss, wing)
    zp1 = ref_single_id(wing.z_ref_points[1]); zp2 = ref_single_id(wing.z_ref_points[2])
    yp1 = ref_single_id(wing.y_ref_points[1]); yp2 = ref_single_id(wing.y_ref_points[2])
    origin = ref_single_id(wing.origin)
    return [
        :zp1x => VIndex(zp1, :pos_1), :zp1y => VIndex(zp1, :pos_2), :zp1z => VIndex(zp1, :pos_3),
        :zp2x => VIndex(zp2, :pos_1), :zp2y => VIndex(zp2, :pos_2), :zp2z => VIndex(zp2, :pos_3),
        :yp1x => VIndex(yp1, :pos_1), :yp1y => VIndex(yp1, :pos_2), :yp1z => VIndex(yp1, :pos_3),
        :yp2x => VIndex(yp2, :pos_1), :yp2y => VIndex(yp2, :pos_2), :yp2z => VIndex(yp2, :pos_3),
        :ovx => VIndex(origin, :vel_1), :ovy => VIndex(origin, :vel_2),
        :ovz => VIndex(origin, :vel_3),
    ]
end

"""
    needs_live_aero(wing)

Whether a PARTICLE wing computes its per-point aero force **live** in the RHS (so its
wing nodes need [`SAM.LiveAeroWingNodePoint`](@ref)) rather than reading a frozen
per-point force. True for every real-aero mode except the frozen-override
[`SAM.AeroDirect`](@ref); false for [`SAM.AeroNone`](@ref) (no force, `is_wing` false).
"""
needs_live_aero(wing) = SAM.is_wing(wing) && !SAM.provides_aero_override(wing.aero)

"""
    wing_node_points(ss, wing)

The wing's aero points in the order [`SAM.wing_points`](@ref) uses — the order the aero
component's `point_pos`/`point_force` connectors are indexed by, so a live wing node's
`aero_slot` is its position in this list.
"""
wing_node_points(ss, wing) =
    [p.idx for p in ss.points if p.is_wing_node && p.wing_idx == wing.idx]

"""
    live_aero_extin(ss, wing, point_idxs)

The `extin` pair list a live-aero wing node reads: the ref-point frame and origin
velocity ([`body_damp_extin`](@ref)), the wing origin position (`opx/opy/opz`), and
every wing point's position/velocity (`wpos_k_c`/`wvel_k_c`, matching
[`SAM.live_aero_node_inputs`](@ref)). Every live-aero node of the wing reads the same
list, so they share one compiled kernel (only the `aero_slot` parameter differs).
"""
function live_aero_extin(ss, wing, point_idxs)
    origin = ref_single_id(wing.origin)
    pairs = body_damp_extin(ss, wing)
    push!(pairs, :opx => VIndex(origin, :pos_1))
    push!(pairs, :opy => VIndex(origin, :pos_2))
    push!(pairs, :opz => VIndex(origin, :pos_3))
    for (k, pidx) in enumerate(point_idxs)
        for c in 1:3
            push!(pairs, Symbol(:wpos_, k, :_, c) => VIndex(pidx, Symbol(:pos_, c)))
            push!(pairs, Symbol(:wvel_, k, :_, c) => VIndex(pidx, Symbol(:vel_, c)))
        end
    end
    return pairs
end

"""
    build_live_aero_vmodels(sam, ss)

One `VertexModel` per PARTICLE wing that computes live aero ([`needs_live_aero`](@ref)),
shared by all that wing's aero nodes (they differ only in the `aero_slot` parameter).
Returns `(vmodel_of, info_of, regs_of)`: `vmodel_of[point_idx]` the shared model,
`info_of[point_idx] = (wing, slot)`, and `regs_of[wing_idx] = (point_reg, aero_reg,
wing, point_idxs)` for per-instance sync. Empty when no wing computes live aero.
"""
function build_live_aero_vmodels(sam, ss)
    vmodel_of = Dict{Int, Any}()
    info_of = Dict{Int, Any}()
    regs_of = Dict{Int, Any}()
    for wing in ss.wings
        (wing.dynamics_type == SAM.PARTICLE_DYNAMICS && needs_live_aero(wing)) ||
            continue
        point_idxs = wing_node_points(ss, wing)
        npts = length(point_idxs)
        repr = point_idxs[1]
        preg = network_view(ss)
        # Full-path (monolith-style) names for the aero subsystem: a live mode may read
        # several twist-surfaces in one kernel, whose leaf names would collide under the
        # network's per-field memoization (AeroPlate's `chord`/`area`), so name by path.
        areg = SAM.ParamView(SAM.ParamRegistry(ss))
        extin = live_aero_extin(ss, wing, point_idxs)
        base = VertexModel(
            SAM.LiveAeroWingNodePoint(sam, preg, areg, repr, wing, 1, npts; name = :live),
            VERTEX_INPUTS, VERTEX_OUTPUTS; extin, mtkcompile = true, name = :live)
        regs_of[wing.idx] = (preg.reg, areg.reg, wing, point_idxs)
        for (slot, pidx) in enumerate(point_idxs)
            vmodel_of[pidx] = base
            info_of[pidx] = (wing, slot)
        end
    end
    return vmodel_of, info_of, regs_of
end

# ======================= network assembly ======================= #

"""
    build_network(sam)

Translate `sam.sys_struct` into a NetworkDynamics `Network` with per-instance
parameters and initial state. Returns `(nw, u0, p0, meta)` with flat state/parameter
vectors and a `meta` named tuple carrying the live `param_sync` and the index maps
the getter/control setters need.
"""
function build_network(sam)
    ss = sam.sys_struct
    points = ss.points
    segments = ss.segments
    n = length(points)

    body_idxs = [b.idx for b in ss.bodies if b.type == SAM.DYNAMIC]
    if !isempty(body_idxs)
        return build_body_mixed_network(sam, ss, body_idxs)
    end

    for point in points
        (point.type == SAM.STATIC || point.type == SAM.DYNAMIC) && continue
        error("NetworkBackend: point $(point.name) has unsupported type " *
              "$(point.type); only STATIC and DYNAMIC are supported so far.")
    end

    pulley_of_point = Dict{Int, Int}()
    for pulley in ss.pulleys
        pulley_of_point[pulley_point_idx(ss, pulley)] = pulley.idx
    end
    winch_of_point = Dict{Int, Int}()
    for winch in ss.winches
        winch_of_point[winch.winch_point_idx] = winch.idx
    end

    graph = SimpleGraph(n)
    seg_of = Dict{Tuple{Int, Int}, Any}()
    for seg in segments
        a, b = seg.point_idxs
        a == b && error("NetworkBackend: segment $(seg.name) is a self-loop " *
                        "(both endpoints are point $a).")
        key = minmax(a, b)
        haskey(seg_of, key) && error("NetworkBackend: points $key are joined " *
            "by more than one segment; parallel edges are not supported.")
        add_edge!(graph, a, b)
        seg_of[key] = seg
    end

    wing_node_vmodel_of, wing_kind_of, kregs =
        build_wing_node_vmodels(sam, ss, pulley_of_point)
    live_vmodel_of, live_info_of, live_regs = build_live_aero_vmodels(sam, ss)

    kind_of = Vector{Symbol}(undef, n)
    for i in 1:n
        if haskey(live_vmodel_of, i)
            kind_of[i] = :live
        elseif haskey(wing_node_vmodel_of, i)
            kind_of[i] = wing_kind_of[i]
        elseif haskey(winch_of_point, i)
            kind_of[i] = :winch
        elseif haskey(pulley_of_point, i)
            kind_of[i] = :pul
        else
            kind_of[i] = points[i].type == SAM.STATIC ? :stat : :dyn
        end
    end

    dyn = build_vertex_kernel!(kregs, sam, ss, :dyn, kind_of, network_dynamic_point)
    stat = build_vertex_kernel!(kregs, sam, ss, :stat, kind_of, network_static_point)
    pulley_v = build_vertex_kernel!(kregs, sam, ss, :pul, kind_of, network_pulley_point)
    winch_v = Dict{Int, Any}()
    for winch in ss.winches
        winch_v[winch.winch_point_idx] = VertexModel(
            network_winch_point(sam, winch, points[winch.winch_point_idx];
                                name = :wch),
            VERTEX_INPUTS, VERTEX_OUTPUTS; mtkcompile = true, name = :wch)
    end

    vmodels = Vector{VertexModel}(undef, n)
    for i in 1:n
        vmodels[i] = kind_of[i] === :live ? live_vmodel_of[i] :
            kind_of[i] === :winch ? winch_v[i] :
            haskey(wing_node_vmodel_of, i) ? wing_node_vmodel_of[i] :
            kind_of[i] === :pul ? pulley_v :
            kind_of[i] === :stat ? stat : dyn
    end

    roles = classify_segments(ss)
    role_of_seg = Dict(segments[k].idx => roles[k] for k in eachindex(segments))
    plain_edge = build_edge_kernel!(kregs, sam, ss, :seg, :plain, segments, role_of_seg,
        network_segment)
    structural_edge = build_edge_kernel!(kregs, sam, ss, :sseg, :structural, segments,
        role_of_seg, network_segment)
    pulley_edge = build_edge_kernel!(kregs, sam, ss, :pseg, :pulley, segments,
        role_of_seg, network_pulley_segment)
    tether_edge_of_tether = build_tether_edges(sam, ss, winch_of_point, kregs,
        segments, role_of_seg)

    edgelist = collect(edges(graph))
    emodels = Vector{EdgeModel}(undef, length(edgelist))
    for (j, e) in enumerate(edgelist)
        seg = seg_of[minmax(src(e), dst(e))]
        role = role_of_seg[seg.idx]
        if role.kind == :pulley
            emodels[j] = pulley_edge
        elseif role.kind == :tether
            emodels[j] = tether_edge_of_tether[role.tether_idx]
        elseif role.kind == :structural
            emodels[j] = structural_edge
        else
            emodels[j] = plain_edge
        end
    end
    nw = Network(graph, vmodels, emodels)

    write_total_mass!(ss)
    param, state = NWParameter(nw), NWState(nw)
    set_states!(ss, state, winch_of_point, pulley_of_point)
    builder = ParamBuilder()
    callables = CallableBuilder()
    record_vertex_params!(builder, ss, winch_of_point, pulley_of_point, kind_of, kregs)
    record_live_aero_params!(builder, callables, ss, live_info_of, live_regs)
    record_edge_params!(builder, ss, edgelist, seg_of, role_of_seg, kregs)
    set_const_params!(nw, param, ss, edgelist, seg_of, role_of_seg, kind_of)
    set_live_aero_slots!(nw, param, live_info_of)
    param_sync = build_network_param_sync(nw, builder, callables)

    meta = (; param_sync, winch_of_point, pulley_of_point,
            winch_tethers = Dict(w.winch_point_idx => collect(w.tether_idxs)
                                 for w in ss.winches),
            body_idxs = Int[])
    return nw, uflat(state), pflat(param), meta
end

"""
    MixedEdgeInfo

Per-edge data the mixed body/point assembly derives once. `kind`: `:plain` (point↔point
spring), `:wrench` (point↔body spring), `:dual_wrench` (body↔body spring, both ends ride
bodies), `:timo_joint`/`:elastic_joint` (body↔body joint). `seg` is the segment (`nothing`
for a joint). `ride_src`/`ride_dst` are the `BODY_STATIC` points absorbed at the edge's src
/dst endpoints (`0` where the endpoint is a free point). `joint_idx` is the joint (`0` for a
segment). `a_at_src` marks joint body A at `src` (baked into the joint kernel).
"""
struct MixedEdgeInfo
    kind::Symbol
    seg::Any
    ride_src::Int
    ride_dst::Int
    joint_idx::Int
    a_at_src::Bool
end

"""
    build_wing_body_vertices(sam, ss, wing_body_idxs)

One wide `VertexModel` per rigid-WING body (`is_wing`): the shared
[`SAM.WingBodyVertex`](@ref) (rigid body + VSM aero) with its own body param registry and
a separate monolith-style aero registry (full-path names, so a multi-twist-surface aero's
leaf names do not collide). Returns `(vm_of, pv_of, areg_of)` keyed by body index, all
empty when no wing body exists.
"""
function build_wing_body_vertices(sam, ss, wing_body_idxs)
    vm_of = Dict{Int, Any}()
    pv_of = Dict{Int, Any}()
    areg_of = Dict{Int, Any}()
    for bidx in wing_body_idxs
        bpv = network_view(ss)
        apv = SAM.ParamView(SAM.ParamRegistry(ss))
        nm = Symbol(:wbody_, bidx)
        vm_of[bidx] = VertexModel(
            SAM.WingBodyVertex(sam, bpv, apv, bidx; name = nm),
            WIDE_VERTEX_INPUTS, WIDE_VERTEX_OUTPUTS; mtkcompile = true, name = nm)
        pv_of[bidx] = bpv
        areg_of[bidx] = apv.reg
    end
    return vm_of, pv_of, areg_of
end

"""
    record_wing_twist_params!(builder, ss, vertex, wing)

Bind a wing-body vertex's per-twist-surface `twist_k`/`twist_vel_k` params (distinct names
per surface, so multiple surfaces do not alias) to live reads of each mapped twist
surface's `twist`/`twist_ω`.
"""
function record_wing_twist_params!(builder, ss, vertex, wing)
    for (k, gidx) in enumerate(wing.twist_surface_idxs)
        add_param!(builder, VIndex(vertex, Symbol(:twist_, k)),
                   SAM.PathReader((:twist_surfaces, gidx, :twist)))
        add_param!(builder, VIndex(vertex, Symbol(:twist_vel_, k)),
                   SAM.PathReader((:twist_surfaces, gidx, :twist_ω)))
    end
    return nothing
end

"""
    record_wing_ride_drag!(builder, ss, vertex, body_idx)

Bind a wing-body vertex's per-ride-point drag params (`drag_anchor_k_c`/`drag_area_k`/
`drag_cd_k`) to live reads of each `BODY_STATIC` ride point (`area > 0`) that
[`SAM.body_ride_drag_wrench`](@ref) drags, in the same filter order.
"""
function record_wing_ride_drag!(builder, ss, vertex, body_idx)
    ride = [p for p in ss.points
            if p.type == SAM.BODY_STATIC && p.body_idx == body_idx && p.area > 0]
    for (k, point) in enumerate(ride)
        for c in 1:3
            add_param!(builder, VIndex(vertex, Symbol(:drag_anchor_, k, :_, c)),
                       SAM.PathReader((:points, point.idx, :anchor_b, c)))
        end
        add_param!(builder, VIndex(vertex, Symbol(:drag_area_, k)),
                   SAM.PathReader((:points, point.idx, :area)))
        add_param!(builder, VIndex(vertex, Symbol(:drag_cd_, k)),
                   SAM.PathReader((:points, point.idx, :drag_coeff)))
    end
    return nothing
end

"""
    build_body_mixed_network(sam, ss, body_idxs)

Assemble a `Network` of integrated rigid bodies (`type == DYNAMIC`) together with free
points, using the wide vertex/edge superset (§8.5). Free `STATIC`/`DYNAMIC` points are
vertices; `BODY_STATIC` points are **absorbed** into their body vertex (their motion is
a body-pose function). A segment between two free points is a `:plain` wide
[`SAM.SpringDamperSegment`](@ref); a segment touching a ride point is a `:wrench`
[`network_wrench_segment`](@ref) delivering force+moment to the body. The bare-body case
(no free points, no segments) is the empty-edge subcase. Pulleys, winches, wing nodes
and rigid-wing aero are not yet supported here. Returns `(nw, u0, p0, meta)`.
"""
function build_body_mixed_network(sam, ss, body_idxs)
    points = ss.points
    pulley_of_point = Dict{Int, Int}()
    for pulley in ss.pulleys
        pulley_of_point[pulley_point_idx(ss, pulley)] = pulley.idx
    end
    winch_of_point = Dict{Int, Int}()
    for winch in ss.winches
        winch_of_point[winch.winch_point_idx] = winch.idx
    end
    free_idxs = Int[]
    for (i, p) in enumerate(points)
        p.type == SAM.BODY_STATIC && continue
        (p.type == SAM.STATIC || p.type == SAM.DYNAMIC) || error(
            "NetworkBackend(body): point $(p.name) type $(p.type) unsupported yet.")
        (p.type == SAM.DYNAMIC && p.is_wing_node) && error(
            "NetworkBackend(body): wing-node point $(p.name) not supported yet.")
        push!(free_idxs, i)
    end
    roles = classify_segments(ss)
    role_of_seg = Dict(ss.segments[k].idx => roles[k] for k in eachindex(ss.segments))
    n_free = length(free_idxs)
    static_body_idxs = [b.idx for b in ss.bodies if b.type == SAM.STATIC]
    all_body_idxs = [body_idxs; static_body_idxs]
    vertex_of_point = Dict(i => v for (v, i) in enumerate(free_idxs))
    vertex_of_body = Dict(bidx => n_free + k for (k, bidx) in enumerate(all_body_idxs))
    nv = n_free + length(all_body_idxs)

    endpoint(pidx) = points[pidx].type == SAM.BODY_STATIC ?
        (vertex_of_body[points[pidx].body_idx], pidx) : (vertex_of_point[pidx], 0)

    graph = SimpleGraph(nv)
    edge_info = Dict{Tuple{Int, Int}, MixedEdgeInfo}()
    add_mixed_edge!(graph, edge_info, key, info) = begin
        haskey(edge_info, key) && error("NetworkBackend(body): parallel edges " *
            "between vertices $key are not supported.")
        add_edge!(graph, key[1], key[2])
        edge_info[key] = info
    end
    is_hermite(pidx) = points[pidx].type == SAM.BODY_STATIC && points[pidx].joint_idx > 0
    for seg in ss.segments
        p1, p2 = seg.point_idxs
        if is_hermite(p1) || is_hermite(p2)
            (is_hermite(p1) && is_hermite(p2)) && error("NetworkBackend(body): " *
                "segment $(seg.name) between two Hermite ride points not supported.")
            hpt = is_hermite(p1) ? p1 : p2
            freept = is_hermite(p1) ? p2 : p1
            points[freept].type == SAM.BODY_STATIC && error("NetworkBackend(body): " *
                "Hermite segment $(seg.name)'s other end must be a free point.")
            vfree = vertex_of_point[freept]
            joint = ss.timoshenko_joints[points[hpt].joint_idx]
            for (to_a, bidx) in ((true, joint.body_a_idx), (false, joint.body_b_idx))
                add_mixed_edge!(graph, edge_info, minmax(vfree, vertex_of_body[bidx]),
                    MixedEdgeInfo(:hermite, seg, hpt, 0, joint.idx, to_a))
            end
            continue
        end
        va, ride_a = endpoint(p1)
        vb, ride_b = endpoint(p2)
        va == vb && error("NetworkBackend(body): segment $(seg.name) is a self-loop " *
            "(both endpoints on vertex $va).")
        ride_src, ride_dst = va < vb ? (ride_a, ride_b) : (ride_b, ride_a)
        role = role_of_seg[seg.idx]
        if ride_src == 0 && ride_dst == 0
            kind = role.kind
        else
            (role.kind == :plain || role.kind == :structural) || error(
                "NetworkBackend(body): a $(role.kind) segment ($(seg.name)) touching a " *
                "body ride point is not supported yet.")
            kind = ride_src > 0 && ride_dst > 0 ? :dual_wrench : :wrench
        end
        add_mixed_edge!(graph, edge_info, minmax(va, vb),
            MixedEdgeInfo(kind, seg, ride_src, ride_dst, 0, false))
    end
    add_joint_edge!(kind, joint) = begin
        haskey(vertex_of_body, joint.body_a_idx) &&
            haskey(vertex_of_body, joint.body_b_idx) || error(
            "NetworkBackend(body): joint $(joint.name) connects a body that is not " *
            "an integrated DYNAMIC body.")
        va = vertex_of_body[joint.body_a_idx]
        vb = vertex_of_body[joint.body_b_idx]
        va == vb && error("NetworkBackend(body): joint $(joint.name) connects a " *
            "body to itself.")
        add_mixed_edge!(graph, edge_info, minmax(va, vb),
            MixedEdgeInfo(kind, nothing, 0, 0, joint.idx, va < vb))
    end
    for joint in ss.timoshenko_joints
        add_joint_edge!(:timo_joint, joint)
    end
    for joint in ss.elastic_joints
        add_joint_edge!(:elastic_joint, joint)
    end

    free_special = union(keys(winch_of_point), keys(pulley_of_point))
    dyn_vm, dyn_reg = build_wide_vertex(sam, ss, free_idxs, SAM.DYNAMIC,
        network_dynamic_point; exclude = free_special)
    stat_vm, stat_reg = build_wide_vertex(sam, ss, free_idxs, SAM.STATIC,
        network_static_point; exclude = free_special)
    pulley_vm, pulley_reg = build_wide_pulley_vertex(sam, ss, free_idxs, pulley_of_point)
    winch_vm_of = build_wide_winch_vertices(sam, ss, free_idxs, winch_of_point)
    wing_body_idxs = [b for b in body_idxs if SAM.is_wing(ss.bodies[b])]
    plain_body_idxs = [b for b in body_idxs if !SAM.is_wing(ss.bodies[b])]
    body_pv = nothing
    body_vm = nothing
    if !isempty(plain_body_idxs)
        body_pv = network_view(ss)
        body_vm = VertexModel(
            SAM.BodyVertex(sam, body_pv, plain_body_idxs[1]; name = :body),
            WIDE_VERTEX_INPUTS, WIDE_VERTEX_OUTPUTS; mtkcompile = true, name = :body)
    end
    wing_vm_of, wing_pv_of, wing_areg_of = build_wing_body_vertices(sam, ss, wing_body_idxs)
    static_body_pv = nothing
    static_body_vm = nothing
    if !isempty(static_body_idxs)
        static_body_pv = network_view(ss)
        static_body_vm = VertexModel(
            SAM.StaticBody(sam, static_body_pv, static_body_idxs[1]; name = :sbody),
            WIDE_VERTEX_INPUTS, WIDE_VERTEX_OUTPUTS; mtkcompile = true, name = :sbody)
    end

    vmodels = Vector{VertexModel}(undef, nv)
    for (v, i) in enumerate(free_idxs)
        vmodels[v] = haskey(winch_of_point, i) ? winch_vm_of[i] :
            haskey(pulley_of_point, i) ? pulley_vm :
            points[i].type == SAM.STATIC ? stat_vm : dyn_vm
    end
    for (k, bidx) in enumerate(all_body_idxs)
        vmodels[n_free + k] = ss.bodies[bidx].type == SAM.STATIC ? static_body_vm :
            haskey(wing_vm_of, bidx) ? wing_vm_of[bidx] : body_vm
    end

    plain_em, plain_reg = build_wide_plain_edge(sam, ss, edge_info)
    wrench_of, wrench_reg_of = build_wrench_edges(sam, ss, edge_info)
    dual_em, dual_reg = build_dual_wrench_edges(sam, ss, edge_info)
    joint_of, joint_reg_of = build_joint_edges(sam, ss, edge_info)
    hermite_of, hermite_reg_of = build_hermite_edges(sam, ss, edge_info, vertex_of_body)
    tether_of, tether_reg = build_wide_tether_edges(sam, ss, winch_of_point, edge_info,
        vertex_of_point)
    pulley_em, pulley_ereg = build_wide_pulley_edge(sam, ss, edge_info)

    edgelist = collect(edges(graph))
    emodels = Vector{EdgeModel}(undef, length(edgelist))
    for (j, e) in enumerate(edgelist)
        key = minmax(src(e), dst(e))
        info = edge_info[key]
        emodels[j] = (info.kind == :timo_joint || info.kind == :elastic_joint) ?
            joint_of[(info.kind, info.a_at_src)] :
            info.kind == :hermite ? hermite_of[key] :
            info.kind == :dual_wrench ? dual_em :
            info.kind == :wrench ? wrench_of[wrench_body_at_src(info)] :
            info.kind == :tether ? tether_of[role_of_seg[info.seg.idx].tether_idx] :
            info.kind == :pulley ? pulley_em : plain_em
    end
    nw = Network(graph, vmodels, emodels)

    body_vertices = [(vertex_of_body[bidx], bidx) for bidx in all_body_idxs]
    dyn_body_vertices = [(vertex_of_body[bidx], bidx) for bidx in body_idxs]
    write_total_mass!(ss)
    param, state = NWParameter(nw), NWState(nw)
    set_body_states!(ss, state, dyn_body_vertices)
    for (v, i) in enumerate(free_idxs)
        if haskey(winch_of_point, i)
            winch = ss.winches[winch_of_point[i]]
            state.v[v, :winch_vel] = winch.vel
            for (pos, tidx) in enumerate(winch.tether_idxs)
                state.v[v, Symbol(:tether_len_, pos)] = ss.tethers[tidx].len
            end
        elseif haskey(pulley_of_point, i)
            pulley = ss.pulleys[pulley_of_point[i]]
            set_particle_state!(state, v, points[i])
            state.v[v, :pulley_len] = pulley.len
            state.v[v, :pulley_vel] = pulley.vel
        elseif points[i].type == SAM.DYNAMIC
            set_particle_state!(state, v, points[i])
        end
    end

    builder = ParamBuilder()
    for (v, i) in enumerate(free_idxs)
        if haskey(winch_of_point, i)
            record_winch_params!(builder, v, i, winch_of_point[i])
        elseif haskey(pulley_of_point, i)
            replay_fields!(builder, pulley_reg, :points, v, i, ss;
                           skip = (:body_frame_damping,))
            record_wind_params!(builder, v)
            record_pulley_mass_params!(builder, ss, v, pulley_of_point[i])
        elseif points[i].type == SAM.STATIC
            replay_fields!(builder, stat_reg, :points, v, i, ss)
        else
            replay_fields!(builder, dyn_reg, :points, v, i, ss;
                           skip = (:body_frame_damping,))
            record_wind_params!(builder, v)
        end
    end
    callables = CallableBuilder()
    for (vertex, bidx) in body_vertices
        if ss.bodies[bidx].type == SAM.STATIC
            record_body_params!(builder, static_body_pv.reg, ss, vertex, bidx;
                                gravity = false)
        elseif haskey(wing_vm_of, bidx)
            record_body_params!(builder, wing_pv_of[bidx].reg, ss, vertex, bidx;
                                gravity = true)
            record_wind_params!(builder, vertex)
            record_aero_params!(builder, callables, ss, vertex, ss.bodies[bidx],
                                wing_areg_of[bidx])
            record_wing_twist_params!(builder, ss, vertex, ss.bodies[bidx])
            record_wing_ride_drag!(builder, ss, vertex, bidx)
        else
            record_body_params!(builder, body_pv.reg, ss, vertex, bidx; gravity = true)
        end
    end
    record_mixed_edge_params!(builder, callables, ss, edgelist, edge_info,
                              plain_reg, wrench_reg_of, dual_reg, joint_reg_of,
                              hermite_reg_of, tether_reg, pulley_ereg, role_of_seg)
    set_mixed_const_params!(nw, param, ss, edgelist, edge_info, role_of_seg)
    param_sync = build_network_param_sync(nw, builder, callables)

    body_static = [(i, points[i].body_idx) for i in eachindex(points)
                   if points[i].type == SAM.BODY_STATIC && points[i].joint_idx == 0]
    winch_tethers = Dict(vertex_of_point[w.winch_point_idx] => collect(w.tether_idxs)
                         for w in ss.winches)
    meta = (; param_sync, winch_of_point, pulley_of_point,
            winch_tethers, body_idxs,
            body_vertices = dyn_body_vertices, body_static, vertex_of_point)
    return nw, uflat(state), pflat(param), meta
end

"""
    build_wide_vertex(sam, ss, free_idxs, ptype, kernelfn)

Compile the wide `VertexModel` for one free-point `ptype` (`STATIC`/`DYNAMIC`) using the
first such point as representative, or `(nothing, nothing)` when none. Returns
`(vmodel, reg)` for per-instance sync.
"""
function build_wide_vertex(sam, ss, free_idxs, ptype, kernelfn; exclude = ())
    repr = findfirst(i -> ss.points[i].type == ptype && !(i in exclude), free_idxs)
    repr === nothing && return nothing, nothing
    pv = network_view(ss)
    vm = VertexModel(kernelfn(sam, pv, free_idxs[repr]; name = :pt, wide = true),
        WIDE_VERTEX_INPUTS, WIDE_VERTEX_OUTPUTS; mtkcompile = true, name = :pt)
    return vm, pv.reg
end

"""
    build_wide_pulley_vertex(sam, ss, free_idxs, pulley_of_point)

Compile the wide pulley `VertexModel` (shared kernel, first pulley free point as
representative) so pulleys coexist with rigid bodies in the mixed network, or
`(nothing, nothing)` when the model has no pulley. Returns `(vmodel, reg)`.
"""
function build_wide_pulley_vertex(sam, ss, free_idxs, pulley_of_point)
    repr = findfirst(i -> haskey(pulley_of_point, i), free_idxs)
    repr === nothing && return nothing, nothing
    pv = network_view(ss)
    vm = VertexModel(
        network_pulley_point(sam, pv, free_idxs[repr]; name = :pul, wide = true),
        WIDE_VERTEX_INPUTS, WIDE_VERTEX_OUTPUTS; mtkcompile = true, name = :pul)
    return vm, pv.reg
end

"""
    build_wide_winch_vertices(sam, ss, free_idxs, winch_of_point)

Compile one wide winch `VertexModel` per winch (each carries its own motor subsystem
and tether-length states), keyed by the winch point index, so winches coexist with
rigid bodies in the mixed network. Returns a `Dict{Int, VertexModel}` (empty if none).
"""
function build_wide_winch_vertices(sam, ss, free_idxs, winch_of_point)
    vm_of = Dict{Int, Any}()
    for winch in ss.winches
        wp = winch.winch_point_idx
        vm_of[wp] = VertexModel(
            SAM.WinchPoint(sam, winch, ss.points[wp]; name = :wch, wide = true),
            WIDE_VERTEX_INPUTS, WIDE_VERTEX_OUTPUTS; mtkcompile = true, name = :wch)
    end
    return vm_of
end

"""
    build_wide_plain_edge(sam, ss, edge_info)

Compile the wide plain point↔point `EdgeModel` (first `:plain` edge as representative),
or `(nothing, nothing)` when none. Returns `(emodel, reg)`.
"""
function build_wide_plain_edge(sam, ss, edge_info)
    repr = nothing
    for info in values(edge_info)
        info.kind == :plain && (repr = info.seg.idx; break)
    end
    repr === nothing && return nothing, nothing
    pv = network_view(ss)
    em = EdgeModel(SAM.SpringDamperSegment(sam, pv, repr; name = :seg, wide = true),
        WIDE_EDGE_SRC_IN, WIDE_EDGE_DST_IN, WIDE_EDGE_SRC_OUT, WIDE_EDGE_DST_OUT;
        mtkcompile = true, name = :seg)
    return em, pv.reg
end

"""
    build_wide_pulley_edge(sam, ss, edge_info)

Compile the wide pulley `EdgeModel` (first `:pulley` edge as representative; the shared
kernel handles both pulley-at-src orientations by a scalar `pulley_at_src` param), or
`(nothing, nothing)` when the model has no pulley. Returns `(emodel, reg)`.
"""
function build_wide_pulley_edge(sam, ss, edge_info)
    repr = nothing
    for info in values(edge_info)
        info.kind == :pulley && (repr = info.seg.idx; break)
    end
    repr === nothing && return nothing, nothing
    pv = network_view(ss)
    em = EdgeModel(network_pulley_segment(sam, pv, repr; name = :pseg, wide = true),
        WIDE_EDGE_SRC_IN, WIDE_EDGE_DST_IN, WIDE_EDGE_SRC_OUT, WIDE_EDGE_DST_OUT;
        mtkcompile = true, name = :pseg)
    return em, pv.reg
end

"""
    build_wide_tether_edges(sam, ss, winch_of_point, edge_info, vertex_of_point)

Compile one wide winched-tether `EdgeModel` per tether, rest length wired to the winch
vertex's `tether_len_k` state through `extin` (the winch **vertex** index, mapped via
`vertex_of_point`). The kernel is compiled once (first `:tether` edge as representative)
and rebound per tether with `EdgeModel(base; extin=…)`. Returns `(edges_of, reg)` keyed
by tether index (empty `Dict`, `nothing` reg when the model has no winch).
"""
function build_wide_tether_edges(sam, ss, winch_of_point, edge_info, vertex_of_point)
    edges_of = Dict{Int, Any}()
    isempty(ss.winches) && return edges_of, nothing
    repr = nothing
    for info in values(edge_info)
        info.kind == :tether && (repr = info.seg.idx; break)
    end
    repr === nothing && return edges_of, nothing
    winch_tether_pos = Dict{Tuple{Int, Int}, Int}()
    for winch in ss.winches, (pos, tidx) in enumerate(winch.tether_idxs)
        winch_tether_pos[(winch.winch_point_idx, tidx)] = pos
    end
    base = nothing
    reg = nothing
    for winch in ss.winches, tidx in winch.tether_idxs
        pos = winch_tether_pos[(winch.winch_point_idx, tidx)]
        sym = Symbol(:tether_len_, pos)
        wv = vertex_of_point[winch.winch_point_idx]
        if base === nothing
            pv = network_view(ss)
            base = EdgeModel(
                network_tether_segment(sam, pv, repr; name = :tseg, wide = true),
                WIDE_EDGE_SRC_IN, WIDE_EDGE_DST_IN, WIDE_EDGE_SRC_OUT, WIDE_EDGE_DST_OUT;
                extin = [:tether_len_ext => VIndex(wv, sym)], mtkcompile = true,
                name = :tseg)
            reg = pv.reg
            edges_of[tidx] = base
        else
            edges_of[tidx] = EdgeModel(base; extin = [VIndex(wv, sym)])
        end
    end
    return edges_of, reg
end

"""
    build_wrench_edges(sam, ss, edge_info)

Compile one wide point↔body wrench `EdgeModel` per distinct `body_at_src` orientation
present (the flag is baked into the kernel, so the two orientations are separate
kernels). Returns `(em_of, reg_of)` — `Dict{Bool,EdgeModel}` and `Dict{Bool,reg}` keyed
by `body_at_src`, using the first wrench edge of each orientation as representative.
"""
wrench_body_at_src(info) = info.ride_src > 0
wrench_ride_idx(info) = max(info.ride_src, info.ride_dst)

function build_wrench_edges(sam, ss, edge_info)
    em_of = Dict{Bool, EdgeModel}()
    reg_of = Dict{Bool, Any}()
    for info in values(edge_info)
        info.kind == :wrench || continue
        bas = wrench_body_at_src(info)
        haskey(em_of, bas) && continue
        pv = network_view(ss)
        em_of[bas] = EdgeModel(
            network_wrench_segment(sam, pv, info.seg.idx, wrench_ride_idx(info), bas;
                                   name = :wseg),
            WIDE_EDGE_SRC_IN, WIDE_EDGE_DST_IN, WIDE_EDGE_SRC_OUT, WIDE_EDGE_DST_OUT;
            mtkcompile = true, name = :wseg)
        reg_of[bas] = pv.reg
    end
    return em_of, reg_of
end

"""
    network_dual_wrench_segment(s, params, idx, ride_src, ride_dst; name)

Wide body↔body spring-damper edge: a segment whose **both** ends are `BODY_STATIC` ride
points, on the `src` and `dst` body vertices. It reconstructs each ride position from that
body's pose + its `anchor_b`, computes the spring-damper load between them (shared
[`SAM.segment_endpoint_loads`](@ref)), and delivers the force + moment-about-COM to each
body. The two anchors use explicitly-named params (`anchor_b_src`/`anchor_b_dst`) since the
network's per-field param naming would otherwise alias two `anchor_b` reads.
"""
function network_dual_wrench_segment(s, params, idx, ride_src, ride_dst; name)
    io, extras = SAM.segment_io_wide()
    (_, src_pos, _, _, dst_pos, _, _,
     src_force, src_mass, src_tension, dst_force, dst_mass, dst_tension) = io
    anchor_src = [SAM.make_param(Symbol(:anchor_b_src_, k), 0.0) for k in 1:3]
    anchor_dst = [SAM.make_param(Symbol(:anchor_b_dst_, k), 0.0) for k in 1:3]
    torn = @variables dw_xsrc(t)[1:3] dw_xdst(t)[1:3] dw_fsrc(t)[1:3] dw_fdst(t)[1:3]
    xsrc, xdst, fsrc, fdst = torn
    R_src = reshape(collect(extras.src_pose_R), 3, 3)
    R_dst = reshape(collect(extras.dst_pose_R), 3, 3)
    com_src = collect(extras.src_pose_com)
    com_dst = collect(extras.dst_pose_com)
    xsrc_expr = collect(src_pos) .+ R_src * collect(anchor_src)
    xdst_expr = collect(dst_pos) .+ R_dst * collect(anchor_dst)
    vel_src = collect(extras.src_pose_com_vel) .+
        (collect(extras.src_pose_omega) × (collect(xsrc) .- com_src))
    vel_dst = collect(extras.dst_pose_com_vel) .+
        (collect(extras.dst_pose_omega) × (collect(xdst) .- com_dst))
    spring, wind = SAM.segment_spring_params(params, idx; with_drag = true)
    l0 = params.segments[idx].l0
    force_src, force_dst, half_mass, _ = SAM.segment_endpoint_loads(
        s, collect(xsrc), vel_src, collect(xdst), vel_dst,
        spring.unit_stiffness, spring.unit_damping, spring.compression_frac, l0,
        spring.diameter, spring.density, spring.cd_tether, collect(wind);
        with_drag = true)
    eqs = [
        collect(xsrc) .~ xsrc_expr
        collect(xdst) .~ xdst_expr
        collect(fsrc) .~ force_src
        collect(fdst) .~ force_dst
        collect(src_force) .~ collect(fsrc)
        src_mass ~ 0.0; src_tension ~ 0.0
        collect(dst_force) .~ collect(fdst)
        dst_mass ~ 0.0; dst_tension ~ 0.0
        collect(extras.src_moment) .~ (collect(xsrc) .- com_src) × collect(fsrc)
        collect(extras.dst_moment) .~ (collect(xdst) .- com_dst) × collect(fdst)
    ]
    return System(eqs, t, [io[1]; torn],
        [param_unknowns(params); anchor_src; anchor_dst]; name)
end

"""
    build_dual_wrench_edges(sam, ss, edge_info)

Compile the wide body↔body `:dual_wrench` `EdgeModel` (one kernel, first such edge
representative), or `(nothing, nothing)` when none. Returns `(emodel, reg)`.
"""
function build_dual_wrench_edges(sam, ss, edge_info)
    repr = nothing
    for info in values(edge_info)
        info.kind == :dual_wrench && (repr = info; break)
    end
    repr === nothing && return nothing, nothing
    pv = network_view(ss)
    em = EdgeModel(
        network_dual_wrench_segment(sam, pv, repr.seg.idx, repr.ride_src, repr.ride_dst;
                                    name = :dwseg),
        WIDE_EDGE_SRC_IN, WIDE_EDGE_DST_IN, WIDE_EDGE_SRC_OUT, WIDE_EDGE_DST_OUT;
        mtkcompile = true, name = :dwseg)
    return em, pv.reg
end

"""
    hermite_extin(ss, other_body_idx, vertex_of_body)

The `extin` pair list binding a Hermite-ride edge's `oth_*` inputs to the *other* beam
body's wide pose outputs (`pos`, `pose_R`, `pose_com`, `pose_com_vel`, `pose_omega`) at
that body's vertex.
"""
function hermite_extin(ss, other_body_idx, vertex_of_body)
    v = vertex_of_body[other_body_idx]
    xyz = ('x', 'y', 'z')
    pairs = Pair{Symbol, Any}[]
    for c in 1:3
        push!(pairs, Symbol(:oth_p, xyz[c]) => VIndex(v, Symbol(:pos_, c)))
        push!(pairs, Symbol(:oth_c, xyz[c]) => VIndex(v, Symbol(:pose_com_, c)))
        push!(pairs, Symbol(:oth_v, xyz[c]) => VIndex(v, Symbol(:pose_com_vel_, c)))
        push!(pairs, Symbol(:oth_w, xyz[c]) => VIndex(v, Symbol(:pose_omega_, c)))
    end
    for k in 1:9
        push!(pairs, Symbol(:oth_R, k) => VIndex(v, Symbol(:pose_R_, k)))
    end
    return pairs
end

"""
    network_hermite_wrench_segment(s, params, seg_idx, joint_idx, ride_idx, to_body_a; name)

One of the **two** wide edges a segment touching a Hermite `BODY_STATIC` ride point spawns
(the point rides `joint`'s deformed centerline, feedthrough on both end bodies). This edge
runs from the free endpoint `K` (`src`) to one beam body (`dst`, body A if `to_body_a`);
it reads the **other** body's pose via `extin` (`oth_*`), reconstructs the identical ride
position via the shared [`SAM.beam_hermite_ride_expressions`](@ref), computes the spring
load between `K` and the ride point, and delivers the axial-fraction share (`1−sfrac` to A,
`sfrac` to B) of both the reaction on `K` and the force+moment on this body. The two edges
together give `K` the full reaction and each body its beam-fraction split. `pos`/force
torn internal; `to_body_a` baked at build time.
"""
function network_hermite_wrench_segment(s, params, seg_idx, joint_idx, ride_idx,
                                        to_body_a::Bool; name)
    io, extras = SAM.segment_io_wide()
    (_, src_pos, src_vel, _, dst_pos, dst_vel, _,
     src_force, src_mass, src_tension, dst_force, dst_mass, dst_tension) = io
    ext = @variables begin
        oth_px(t), [input = true]; oth_py(t), [input = true]; oth_pz(t), [input = true]
        oth_cx(t), [input = true]; oth_cy(t), [input = true]; oth_cz(t), [input = true]
        oth_vx(t), [input = true]; oth_vy(t), [input = true]; oth_vz(t), [input = true]
        oth_wx(t), [input = true]; oth_wy(t), [input = true]; oth_wz(t), [input = true]
        oth_R1(t), [input = true]; oth_R2(t), [input = true]; oth_R3(t), [input = true]
        oth_R4(t), [input = true]; oth_R5(t), [input = true]; oth_R6(t), [input = true]
        oth_R7(t), [input = true]; oth_R8(t), [input = true]; oth_R9(t), [input = true]
    end
    (oth_px, oth_py, oth_pz, oth_cx, oth_cy, oth_cz, oth_vx, oth_vy, oth_vz,
     oth_wx, oth_wy, oth_wz, oth_R1, oth_R2, oth_R3, oth_R4, oth_R5, oth_R6,
     oth_R7, oth_R8, oth_R9) = ext
    oth_pos = [oth_px, oth_py, oth_pz]
    oth_com = [oth_cx, oth_cy, oth_cz]
    oth_comvel = [oth_vx, oth_vy, oth_vz]
    oth_omega = [oth_wx, oth_wy, oth_wz]
    oth_R = [oth_R1, oth_R2, oth_R3, oth_R4, oth_R5, oth_R6, oth_R7, oth_R8, oth_R9]
    joint = params.reg.sys_struct.timoshenko_joints[joint_idx]
    this_body = (dst_pos, extras.dst_pose_R, extras.dst_pose_com,
                 extras.dst_pose_com_vel, extras.dst_pose_omega)
    other_body = (oth_pos, oth_R, oth_com, oth_comvel, oth_omega)
    (apos, aR, acom, acomvel, aomega), (bpos, bR, bcom, bcomvel, bomega) =
        to_body_a ? (this_body, other_body) : (other_body, this_body)
    torn = @variables hm_pos(t)[1:3] hm_fk(t)[1:3] hm_fp(t)[1:3] hm_frame(t)[1:3, 1:3] hm_ta(t)[1:3] hm_tb(t)[1:3]
    hpos, hfk, hfp, hframe, hta, htb = torn
    hk = SAM.beam_hermite_ride_expressions(joint, params, ride_idx;
        pos_a = apos, R_a = reshape(collect(aR), 3, 3), com_a = acom,
        com_vel_a = acomvel, omega_a_w = aomega,
        pos_b = bpos, R_b = reshape(collect(bR), 3, 3), com_b = bcom,
        com_vel_b = bcomvel, omega_b_w = bomega,
        frame = hframe, theta_a = hta, theta_b = htb)
    spring, wind = SAM.segment_spring_params(params, seg_idx; with_drag = true)
    l0 = params.segments[seg_idx].l0
    fk, fp, half_mass, _ = SAM.segment_endpoint_loads(
        s, collect(src_pos), collect(src_vel), collect(hpos),
        hk.ride_velocity(collect(hpos)),
        spring.unit_stiffness, spring.unit_damping, spring.compression_frac, l0,
        spring.diameter, spring.density, spring.cd_tether, collect(wind);
        with_drag = true)
    frac = to_body_a ? (1 - hk.sfrac) : hk.sfrac
    knot_mass = params.points[ride_idx].extra_mass + half_mass
    gravity_knot = [0.0, 0.0, -params.set.g_earth * knot_mass]
    force_on_point = collect(hfp) .+ gravity_knot
    force_body = frac .* force_on_point
    moment_body = (collect(hpos) .- collect(extras.dst_pose_com)) × force_body
    eqs = [
        hk.tear_eqs
        collect(hpos) .~ hk.pos_point
        collect(hfk) .~ fk
        collect(hfp) .~ fp
        collect(src_force) .~ frac .* collect(hfk)
        src_mass ~ frac * half_mass; src_tension ~ 0.0
        collect(dst_force) .~ force_body
        dst_mass ~ 0.0; dst_tension ~ 0.0
        collect(extras.src_moment) .~ zeros(3)
        collect(extras.dst_moment) .~ moment_body
    ]
    return System(eqs, t, [io[1]; ext; torn], param_unknowns(params); name)
end

"""
    build_hermite_edges(sam, ss, edge_info, vertex_of_body)

Compile the wide Hermite-ride `EdgeModel`s. Each `:hermite` edge is a free-point→beam-body
edge that reads the *other* beam body via `extin`; the kernel is compiled once per
`to_body_a` orientation and rebound per edge with `EdgeModel(base; extin=…)` (like the
winched-tether edges), since only the other body's vertex differs. Returns
`(em_of, reg_of)` — `em_of` keyed by the graph vertex-pair, `reg_of` by `to_body_a`.
"""
function build_hermite_edges(sam, ss, edge_info, vertex_of_body)
    em_of = Dict{Tuple{Int, Int}, EdgeModel}()
    reg_of = Dict{Bool, Any}()
    base_of = Dict{Bool, Any}()
    for (key, info) in edge_info
        info.kind == :hermite || continue
        to_a = info.a_at_src
        joint = ss.timoshenko_joints[info.joint_idx]
        other = to_a ? joint.body_b_idx : joint.body_a_idx
        extin = hermite_extin(ss, other, vertex_of_body)
        if haskey(base_of, to_a)
            em_of[key] = EdgeModel(base_of[to_a]; extin = last.(extin))
        else
            pv = network_view(ss)
            base = EdgeModel(
                network_hermite_wrench_segment(sam, pv, info.seg.idx, info.joint_idx,
                                               info.ride_src, to_a; name = :hmseg),
                WIDE_EDGE_SRC_IN, WIDE_EDGE_DST_IN, WIDE_EDGE_SRC_OUT, WIDE_EDGE_DST_OUT;
                extin, mtkcompile = true, name = :hmseg)
            base_of[to_a] = base
            reg_of[to_a] = pv.reg
            em_of[key] = base
        end
    end
    return em_of, reg_of
end

"""
    build_joint_edges(sam, ss, edge_info)

Compile one wide body↔body joint `EdgeModel` per distinct `(kind, body_at_src)` present —
`:timo_joint` → [`network_timoshenko_joint_edge`](@ref), `:elastic_joint` →
[`network_elastic_joint_edge`](@ref) — using the first edge of each key as representative.
Returns `(em_of, reg_of)` keyed by `(kind, body_at_src)`.
"""
function build_joint_edges(sam, ss, edge_info)
    em_of = Dict{Tuple{Symbol, Bool}, EdgeModel}()
    reg_of = Dict{Tuple{Symbol, Bool}, Any}()
    kernelfn = Dict(:timo_joint => network_timoshenko_joint_edge,
                    :elastic_joint => network_elastic_joint_edge)
    for info in values(edge_info)
        (info.kind == :timo_joint || info.kind == :elastic_joint) || continue
        key = (info.kind, info.a_at_src)
        haskey(em_of, key) && continue
        pv = network_view(ss)
        em_of[key] = EdgeModel(
            kernelfn[info.kind](sam, pv, info.joint_idx, info.a_at_src;
                                name = :joint),
            WIDE_EDGE_SRC_IN, WIDE_EDGE_DST_IN, WIDE_EDGE_SRC_OUT, WIDE_EDGE_DST_OUT;
            mtkcompile = true, name = :joint)
        reg_of[key] = pv.reg
    end
    return em_of, reg_of
end

"""
    set_body_states!(ss, state, body_vertices)

Fill each body vertex's 13 principal-frame states (`com_w`, `com_vel`, `Q_p_to_w`,
`omega_p`) from the struct fields `reinit!`/`init_principal_frame!` already populated,
so the network starts from the same pose as the monolith. `body_vertices` is a list of
`(graph_vertex, body_idx)` pairs (a body vertex sits after the free-point vertices).
"""
function set_body_states!(ss, state, body_vertices)
    for (vertex, bidx) in body_vertices
        body = ss.bodies[bidx]
        for k in 1:3
            state.v[vertex, Symbol(:com_w_, k)] = body.com_w[k]
            state.v[vertex, Symbol(:com_vel_, k)] = body.com_vel[k]
            state.v[vertex, Symbol(:omega_p_, k)] = body.ω_p[k]
        end
        for k in 1:4
            state.v[vertex, Symbol(:Q_, k)] = body.Q_p_to_w[k]
        end
    end
    return nothing
end

"""
    record_body_params!(builder, reg, ss, vertex, bidx)

Record body vertex `vertex`'s per-instance parameters: every plain struct field the
kernel read from `bodies[bidx]` (replayed generically, matrix `R_b_to_p` included) plus
the world gravity `set.g_earth`, each bound as a live read so mutating the body struct
takes effect next step, as on the monolith.
"""
function record_body_params!(builder, reg, ss, vertex, bidx; gravity = true)
    replay_fields!(builder, reg, :bodies, vertex, bidx, ss)
    gravity && add_param!(builder, VIndex(vertex, :g_earth),
                          SAM.PathReader((:set, :g_earth)))
    return nothing
end

"""
    record_mixed_edge_params!(builder, callables, ss, edgelist, edge_info, plain_reg,
                              wrench_reg_of, dual_reg, joint_reg_of, hermite_reg_of,
                              tether_reg, pulley_reg, role_of_seg)

Record each mixed body/point edge's live parameters: the segment's spring/drag fields
(replayed from the matching kernel registry — `wrench_reg_of[body_at_src]` for a wrench
edge) and `cd_tether`/wind. A `:wrench` edge also binds its ride point's `anchor_b`; a
`:tether`/`:pulley` edge replays from its own registry (`pulley` also binds the pulley
`sum_len`). `body_at_src` and the pulley/tether topology constants are baked/const, so
nothing is set here for them.
"""
function record_mixed_edge_params!(builder, callables, ss, edgelist, edge_info,
                                   plain_reg, wrench_reg_of, dual_reg, joint_reg_of,
                                   hermite_reg_of, tether_reg, pulley_reg, role_of_seg)
    for (j, e) in enumerate(edgelist)
        info = edge_info[minmax(src(e), dst(e))]
        if info.kind == :tether
            replay_fields!(builder, tether_reg, :segments, j, info.seg.idx, ss;
                           edge = true)
            add_param!(builder, EIndex(j, :cd_tether),
                       SAM.PathReader((:set, :cd_tether)))
            record_wind_params!(builder, j; edge = true)
            continue
        end
        if info.kind == :pulley
            replay_fields!(builder, pulley_reg, :segments, j, info.seg.idx, ss;
                           edge = true)
            add_param!(builder, EIndex(j, :cd_tether),
                       SAM.PathReader((:set, :cd_tether)))
            record_wind_params!(builder, j; edge = true)
            add_param!(builder, EIndex(j, :pulley_sum_len),
                SAM.PathReader((:pulleys, role_of_seg[info.seg.idx].pulley_idx,
                               :sum_len)))
            continue
        end
        if info.kind == :timo_joint || info.kind == :elastic_joint
            container = info.kind == :timo_joint ? :timoshenko_joints : :elastic_joints
            reg = joint_reg_of[(info.kind, info.a_at_src)]
            replay_fields!(builder, reg, container, j, info.joint_idx, ss; edge = true)
            record_edge_callables!(callables, reg, container, j, info.joint_idx, ss)
            continue
        end
        if info.kind == :hermite
            reg = hermite_reg_of[info.a_at_src]
            replay_fields!(builder, reg, :segments, j, info.seg.idx, ss; edge = true)
            replay_fields!(builder, reg, :points, j, info.ride_src, ss; edge = true)
            replay_fields!(builder, reg, :timoshenko_joints, j, info.joint_idx, ss;
                           edge = true)
            add_param!(builder, EIndex(j, :cd_tether),
                       SAM.PathReader((:set, :cd_tether)))
            add_param!(builder, EIndex(j, :g_earth),
                       SAM.PathReader((:set, :g_earth)))
            record_wind_params!(builder, j; edge = true)
            continue
        end
        seg = info.seg
        segment_stiffness(seg)
        reg = info.kind == :dual_wrench ? dual_reg :
            info.kind == :wrench ? wrench_reg_of[wrench_body_at_src(info)] : plain_reg
        replay_fields!(builder, reg, :segments, j, seg.idx, ss; edge = true)
        add_param!(builder, EIndex(j, :cd_tether), SAM.PathReader((:set, :cd_tether)))
        record_wind_params!(builder, j; edge = true)
        if info.kind == :wrench
            for k in 1:3
                add_param!(builder, EIndex(j, Symbol(:anchor_b_, k)),
                    SAM.PathReader((:points, wrench_ride_idx(info), :anchor_b, k)))
            end
        elseif info.kind == :dual_wrench
            for k in 1:3
                add_param!(builder, EIndex(j, Symbol(:anchor_b_src_, k)),
                           SAM.PathReader((:points, info.ride_src, :anchor_b, k)))
                add_param!(builder, EIndex(j, Symbol(:anchor_b_dst_, k)),
                           SAM.PathReader((:points, info.ride_dst, :anchor_b, k)))
            end
        end
    end
    return nothing
end

"""
    build_tether_edges(sam, ss, winch_of_point, kregs, segments, role_of_seg)

One `network_tether_segment` `EdgeModel` per winched tether, its rest length wired
to the winch vertex's matching `tether_len_k` state via `extin`. The kernel is
compiled once (with the first tether segment as representative, its `params`
registry stashed in `kregs[:tseg]`) and rebound per tether with
`EdgeModel(base; extin=…)` so only the external-input index differs.
"""
function build_tether_edges(sam, ss, winch_of_point, kregs, segments, role_of_seg)
    edges_of = Dict{Int, Any}()
    isempty(ss.winches) && return edges_of
    repr = role_representative(segments, role_of_seg, :tether)
    winch_tether_pos = Dict{Tuple{Int, Int}, Int}()
    for winch in ss.winches, (pos, tidx) in enumerate(winch.tether_idxs)
        winch_tether_pos[(winch.winch_point_idx, tidx)] = pos
    end
    base = nothing
    for winch in ss.winches, tidx in winch.tether_idxs
        pos = winch_tether_pos[(winch.winch_point_idx, tidx)]
        sym = Symbol(:tether_len_, pos)
        extin = [:tether_len_ext => VIndex(winch.winch_point_idx, sym)]
        if base === nothing
            pv = network_view(ss)
            base = EdgeModel(network_tether_segment(sam, pv, repr; name = :tseg),
                EDGE_SRC_IN, EDGE_DST_IN, EDGE_SRC_OUT, EDGE_DST_OUT;
                extin, mtkcompile = true, name = :tseg)
            kregs[:tseg] = pv.reg
            edges_of[tidx] = base
        else
            edges_of[tidx] = EdgeModel(base;
                extin = [VIndex(winch.winch_point_idx, sym)])
        end
    end
    return edges_of
end

"""
    build_vertex_kernel!(kregs, sam, ss, kind, kind_of, kernelfn)

Compile the vertex `System` for one non-winch `kind` (`:dyn`/`:stat`/`:pul`) into a
`VertexModel`, using the first vertex of that kind as the representative index for
`kernelfn(sam, params, idx)` and stashing the fresh `params` registry in `kregs` for
per-instance sync. Returns `nothing` when no vertex has that kind.
"""
function build_vertex_kernel!(kregs, sam, ss, kind::Symbol, kind_of, kernelfn)
    repr = findfirst(==(kind), kind_of)
    repr === nothing && return nothing
    pv = network_view(ss)
    vm = VertexModel(kernelfn(sam, pv, repr; name = kind),
        VERTEX_INPUTS, VERTEX_OUTPUTS; mtkcompile = true, name = kind)
    kregs[kind] = pv.reg
    return vm
end

"""
    role_representative(segments, role_of_seg, role_kind)

The `idx` of the first segment whose role is `role_kind` (`:plain`/`:pulley`/
`:tether`), used as the representative index for that edge kernel; `nothing` when no
segment has that role.
"""
function role_representative(segments, role_of_seg, role_kind::Symbol)
    for seg in segments
        role_of_seg[seg.idx].kind === role_kind && return seg.idx
    end
    return nothing
end

"""
    build_edge_kernel!(kregs, sam, ss, name, role_kind, segments, role_of_seg, kernelfn)

Compile the edge `System` for one `role_kind` into an `EdgeModel`, using the first
segment of that role as the representative index for `kernelfn(sam, params, idx)` and
stashing the fresh `params` registry in `kregs[name]` for per-instance sync. Returns
`nothing` when no segment has that role.
"""
function build_edge_kernel!(kregs, sam, ss, name::Symbol, role_kind::Symbol,
                            segments, role_of_seg, kernelfn)
    repr = role_representative(segments, role_of_seg, role_kind)
    repr === nothing && return nothing
    pv = network_view(ss)
    em = EdgeModel(kernelfn(sam, pv, repr; name),
        EDGE_SRC_IN, EDGE_DST_IN, EDGE_SRC_OUT, EDGE_DST_OUT;
        mtkcompile = true, name)
    kregs[name] = pv.reg
    return em
end

"""
    build_wing_node_vmodels(sam, ss, pulley_of_point)

One `VertexModel` per DYNAMIC wing node on a KINEMATIC wing (selected by
[`network_wing_node`](@ref)), reading its wing's ref-point frame through `extin` for
the frozen aero force and body-frame damping. The dynamic and pulley kernels are
each compiled once and rebound per point with `VertexModel(base; extin=…)`, so only
the wing's ref-point indices differ. Returns a `point_idx => VertexModel` map (empty
when no point is a wing node).
"""
function build_wing_node_vmodels(sam, ss, pulley_of_point)
    vmodel_of = Dict{Int, Any}()
    kind_of = Dict{Int, Symbol}()
    kregs = Dict{Symbol, Any}()
    dyn_base = nothing
    pulley_base = nothing
    for (i, point) in enumerate(ss.points)
        wing = network_wing_node(ss, point)
        wing === nothing && continue
        # Live-aero wing nodes get the live kernel (build_live_aero_vmodels), not the
        # frozen wing-node kernel; other damped nodes on a live wing stay frozen.
        (point.is_wing_node && needs_live_aero(wing)) && continue
        extin = body_damp_extin(ss, wing)
        if haskey(pulley_of_point, i)
            kind_of[i] = :wnpul
            if pulley_base === nothing
                pv = network_view(ss)
                pulley_base = VertexModel(
                    network_wing_node_pulley_point(sam, pv, i; name = :wnpul),
                    VERTEX_INPUTS, VERTEX_OUTPUTS; extin, mtkcompile = true,
                    name = :wnpul)
                kregs[:wnpul] = pv.reg
                vmodel_of[i] = pulley_base
            else
                vmodel_of[i] = VertexModel(pulley_base; extin = last.(extin))
            end
        else
            kind_of[i] = :wnode
            if dyn_base === nothing
                pv = network_view(ss)
                dyn_base = VertexModel(network_wing_node_point(sam, pv, i; name = :wnode),
                    VERTEX_INPUTS, VERTEX_OUTPUTS; extin, mtkcompile = true,
                    name = :wnode)
                kregs[:wnode] = pv.reg
                vmodel_of[i] = dyn_base
            else
                vmodel_of[i] = VertexModel(dyn_base; extin = last.(extin))
            end
        end
    end
    return vmodel_of, kind_of, kregs
end

"""
    record_aero_params!(builder, callables, ss, vertex, wing, areg)

Record a live-aero node's aero-component parameters (namespaced `aero₊…` after nesting)
for per-instance sync: the frozen induced velocity `v_ind` (a `3×n_panels` matrix,
scalarized column-major to `aero₊v_ind_c_i` with a linear reader index), any scalar/
vector aero params, and the `cl/cd/cm` polar callables (routed to `callables`, the
nonnumeric SII path). Every entry reads live from `ss.wings[wing.idx].aero`.
"""
function record_aero_params!(builder, callables, ss, vertex, wing, areg)
    for entry in areg.entries
        entry.read isa SAM.PathReader || continue
        path = entry.read.path
        field = SAM.param_name(path)  # full-path name (matches the monolith-style areg)
        if entry.kind === :callable
            add_callable!(callables, VIndex(vertex, Symbol("aero₊", field)),
                          SAM.PathReader(path))
        elseif entry.kind === :scalar
            add_param!(builder, VIndex(vertex, Symbol("aero₊", field)),
                       SAM.PathReader(path))
        elseif entry.kind === :array
            value = entry.read(ss)
            if ndims(value) == 2
                rows, cols = size(value)
                for col in 1:cols, row in 1:rows
                    add_param!(builder,
                        VIndex(vertex, Symbol("aero₊", field, "_", row, "_", col)),
                        SAM.PathReader((path..., (col - 1) * rows + row)))
                end
            else
                for k in 1:length(value)
                    add_param!(builder,
                        VIndex(vertex, Symbol("aero₊", field, "_", k)),
                        SAM.PathReader((path..., k)))
                end
            end
        end
    end
    return nothing
end

"""
    record_live_aero_params!(builder, callables, ss, live_info_of, live_regs)

Record every live-aero wing node's per-instance parameters: its point mass/drag/damping
(replayed from the point registry, plus the ground wind and body-frame damping like a
frozen wing node) and its aero-component parameters ([`record_aero_params!`](@ref)). The
`aero_slot` selecting this node's own force is a topology constant set by
[`set_live_aero_slots!`](@ref).
"""
function record_live_aero_params!(builder, callables, ss, live_info_of, live_regs)
    for (i, (wing, _)) in live_info_of
        preg, areg, _, _ = live_regs[wing.idx]
        replay_fields!(builder, preg, :points, i, i, ss;
                       skip = (:body_frame_damping,))
        record_wind_params!(builder, i)
        bind_body_damp!(builder, ss, i, ss.points[i])
        record_aero_params!(builder, callables, ss, i, wing, areg)
    end
    return nothing
end

"""
    set_live_aero_slots!(nw, param, live_info_of)

Write each live-aero node's `aero_slot` — its position in the wing's aero-point list, so
one shared kernel selects the right per-point force — straight into `param` once (a fixed
assembly constant, not re-synced).
"""
function set_live_aero_slots!(nw, param, live_info_of)
    isempty(live_info_of) && return nothing
    indices = Any[]
    values = SAM.SimFloat[]
    for (i, (_, slot)) in live_info_of
        push!(indices, VIndex(i, :aero_slot))
        push!(values, SAM.SimFloat(slot))
    end
    setp(nw, indices)(param, values)
    return nothing
end

"""
    write_total_mass!(ss)

Effective translational mass per point (`extra_mass` + incident half-masses),
written to the struct so `validate_sys_struct` and diagnostics see it.
"""
function write_total_mass!(ss)
    for point in ss.points
        point.total_mass = point.extra_mass
    end
    for seg in ss.segments
        half = SAM.segment_half_mass(seg.l0, seg.diameter, seg.density)
        ss.points[seg.point_idxs[1]].total_mass += half
        ss.points[seg.point_idxs[2]].total_mass += half
    end
    return nothing
end

# ======================= initial state ======================= #

"""
    set_states!(ss, state, winch_of_point, pulley_of_point)

Fill each vertex's initial *state* from the struct: particle `pos`/`vel`, the
pulley split `pulley_len`/`pulley_vel`, and the winch `winch_vel` + `tether_len_k`.
Static/winch positions are algebraic (a `pos_w` parameter), so carry no state.
"""
function set_states!(ss, state, winch_of_point, pulley_of_point)
    for (i, point) in enumerate(ss.points)
        if haskey(winch_of_point, i)
            winch = ss.winches[winch_of_point[i]]
            state.v[i, :winch_vel] = winch.vel
            for (pos, tidx) in enumerate(winch.tether_idxs)
                state.v[i, Symbol(:tether_len_, pos)] = ss.tethers[tidx].len
            end
        elseif haskey(pulley_of_point, i)
            pulley = ss.pulleys[pulley_of_point[i]]
            set_particle_state!(state, i, point)
            state.v[i, :pulley_len] = pulley.len
            state.v[i, :pulley_vel] = pulley.vel
        elseif point.type == SAM.STATIC
            # algebraic position; nothing to initialise
        else
            set_particle_state!(state, i, point)
        end
    end
    return nothing
end

"""
    set_particle_state!(state, i, point)

Set the initial `pos`/`vel` (scalarized `pos_k`/`vel_k`) for particle vertex `i`.
"""
function set_particle_state!(state, i, point)
    for k in 1:3
        state.v[i, Symbol(:pos_, k)] = point.pos_w[k]
        state.v[i, Symbol(:vel_, k)] = point.vel_w[k]
    end
    return nothing
end

# ======================= parameter recording ======================= #

"""
    record_vertex_params!(builder, ss, winch_of_point, pulley_of_point, kind_of, kregs)

Record every vertex's per-instance parameters (index + live reader) into `builder`,
dispatched on `kind_of[i]`. Struct fields are replayed generically from the kernel's
registry (`kregs`); the `set.wind_vec` ground wind, the pulley rope-mass fields, the
winch controls and the `body_frame_damping` `nothing`-fallback are bound explicitly.
"""
function record_vertex_params!(builder, ss, winch_of_point, pulley_of_point,
                               kind_of, kregs)
    for (i, point) in enumerate(ss.points)
        kind = kind_of[i]
        kind === :live && continue  # recorded by record_live_aero_params!
        if kind === :winch
            record_winch_params!(builder, i, i, winch_of_point[i])
            continue
        end
        if kind === :stat
            replay_fields!(builder, kregs[:stat], :points, i, i, ss)
            continue
        end
        replay_fields!(builder, kregs[kind], :points, i, i, ss;
                       skip = (:body_frame_damping,))
        record_wind_params!(builder, i)
        (kind === :wnode || kind === :wnpul) && bind_body_damp!(builder, ss, i, point)
        (kind === :pul || kind === :wnpul) &&
            record_pulley_mass_params!(builder, ss, i, pulley_of_point[i])
    end
    return nothing
end

"""
    record_wind_params!(builder, addr; edge=false)

Bind the ground wind `set.wind_vec` (scalarized `wind_vec_1..3`) that
[`SAM.ground_wind_vec`](@ref) reads, on vertex (`edge=false`) or edge (`edge=true`)
`addr`.
"""
function record_wind_params!(builder, addr; edge = false)
    for k in 1:3
        add_param!(builder, param_addr(edge, addr, Symbol(:wind_vec_, k)),
                   SAM.PathReader((:set, :wind_vec, k)))
    end
    return nothing
end

"""
    record_pulley_mass_params!(builder, ss, i, pulley_idx)

Bind vertex `i`'s rope-mass parameters — the pulley `sum_len` and its first segment's
`density`/`diameter` — that [`pulley_rope_mass`](@ref) multiplies in-equation.
"""
function record_pulley_mass_params!(builder, ss, i, pulley_idx)
    seg1 = ss.pulleys[pulley_idx].segment_idxs[1]
    add_param!(builder, VIndex(i, :sum_len),
               SAM.PathReader((:pulleys, pulley_idx, :sum_len)))
    add_param!(builder, VIndex(i, :density),
               SAM.PathReader((:segments, seg1, :density)))
    add_param!(builder, VIndex(i, :diameter),
               SAM.PathReader((:segments, seg1, :diameter)))
    return nothing
end

"""
    bind_body_damp!(builder, ss, i, point)

Bind a wing-node vertex `i`'s body-frame damping `body_frame_damping_k` as a live
struct read from the point's `body_frame_damping`. An aero-only wing node (damping
`nothing`) is skipped here and zeroed once by [`set_const_params!`](@ref).
"""
function bind_body_damp!(builder, ss, i, point)
    point_body_damp(ss, point) === nothing && return nothing
    for k in 1:3
        add_param!(builder, VIndex(i, Symbol(:body_frame_damping_, k)),
                   SAM.PathReader((:points, i, :body_frame_damping, k)))
    end
    return nothing
end

"""
    record_pos_w_params!(builder, vertex, point_idx)

Record the anchored position `pos_w_k` for a static/winch `vertex`, reading the struct
field of point `point_idx` (which may differ from the vertex index in a body network).
"""
function record_pos_w_params!(builder, vertex, point_idx)
    for k in 1:3
        add_param!(builder, VIndex(vertex, Symbol(:pos_w_, k)),
                   SAM.PathReader((:points, point_idx, :pos_w, k)))
    end
    return nothing
end

"""
    record_winch_params!(builder, vertex, point_idx, widx)

Record the winch vertex parameters for `vertex` (winch `widx` at `point_idx`): anchored
position plus `brake` and `speed_controlled`. `set_value` is deliberately excluded — it
is a control input owned by `set_set_values`/`get_set_values`, and syncing it here would
clobber the setpoint the control setter writes each step.
"""
function record_winch_params!(builder, vertex, point_idx, widx)
    record_pos_w_params!(builder, vertex, point_idx)
    add_param!(builder, VIndex(vertex, :brake),
               SAM.PathReader((:winches, widx, :brake)))
    add_param!(builder, VIndex(vertex, :speed_controlled),
               SAM.PathReader((:winches, widx, :speed_controlled)))
    return nothing
end

"""
    record_edge_params!(builder, ss, edgelist, seg_of, role_of_seg)

Record every edge's live spring/drag parameters (struct-field reads) into `builder`.
The pulley `sum_len` is a struct read; the assembly-fixed pulley/tether constants
(side, at-source flag, segment count, tension signs) are written once by
[`set_const_params!`](@ref), not synced.
"""
function record_edge_params!(builder, ss, edgelist, seg_of, role_of_seg, kregs)
    reg_of_role = Dict(:plain => :seg, :structural => :sseg,
                       :pulley => :pseg, :tether => :tseg)
    for (j, e) in enumerate(edgelist)
        seg = seg_of[minmax(src(e), dst(e))]
        role = role_of_seg[seg.idx]
        segment_stiffness(seg)  # validate linear stiffness up front
        replay_fields!(builder, kregs[reg_of_role[role.kind]], :segments, j, seg.idx, ss;
                       edge = true)
        if role.kind !== :structural
            add_param!(builder, EIndex(j, :cd_tether), SAM.PathReader((:set, :cd_tether)))
        end
        record_wind_params!(builder, j; edge = true)
        role.kind == :pulley && add_param!(builder, EIndex(j, :pulley_sum_len),
            SAM.PathReader((:pulleys, role.pulley_idx, :sum_len)))
    end
    return nothing
end

"""
    set_const_params!(nw, param, ss, edgelist, seg_of, role_of_seg, kind_of)

Write the assembly-fixed structural constants straight into `param` once — no reader,
no per-step sync. Each pulley edge gets its `pulley_side`/`pulley_at_src`, each
winched-tether edge its `n_segs`/`tension_sign_src`/`tension_sign_dst`, and every
aero-only wing node (body damping `nothing`) a zero `body_frame_damping`. These are
topology, fixed at assembly, so they are set here rather than re-read every step.
"""
function set_const_params!(nw, param, ss, edgelist, seg_of, role_of_seg, kind_of)
    indices = Any[]
    values = SAM.SimFloat[]
    for (j, e) in enumerate(edgelist)
        role = role_of_seg[seg_of[minmax(src(e), dst(e))].idx]
        if role.kind == :pulley
            append!(indices, [EIndex(j, :pulley_side), EIndex(j, :pulley_at_src)])
            append!(values, [role.pulley_side, role.pulley_at_src ? 1.0 : 0.0])
        elseif role.kind == :tether
            append!(indices, [EIndex(j, :n_segs), EIndex(j, :tension_sign_src),
                              EIndex(j, :tension_sign_dst)])
            append!(values, SAM.SimFloat[role.n_segs, role.tension_sign_src,
                                         role.tension_sign_dst])
        end
    end
    for (i, point) in enumerate(ss.points)
        (kind_of[i] === :wnode || kind_of[i] === :wnpul || kind_of[i] === :live) ||
            continue
        point_body_damp(ss, point) === nothing || continue
        for k in 1:3
            push!(indices, VIndex(i, Symbol(:body_frame_damping_, k)))
            push!(values, 0.0)
        end
    end
    isempty(indices) || setp(nw, indices)(param, values)
    return nothing
end

"""
    set_mixed_const_params!(nw, param, ss, edgelist, edge_info, role_of_seg)

Write the assembly-fixed pulley/tether constants into a mixed body network's `param`
once (no reader): each `:pulley` edge its `pulley_side`/`pulley_at_src`, each `:tether`
edge its `n_segs`/`tension_sign_src`/`tension_sign_dst`. Mirrors the edge portion of
[`set_const_params!`](@ref), keyed by the mixed-edge `edge_info`.
"""
function set_mixed_const_params!(nw, param, ss, edgelist, edge_info, role_of_seg)
    indices = Any[]
    values = SAM.SimFloat[]
    for (j, e) in enumerate(edgelist)
        info = edge_info[minmax(src(e), dst(e))]
        (info.kind == :pulley || info.kind == :tether) || continue
        role = role_of_seg[info.seg.idx]
        if info.kind == :pulley
            append!(indices, [EIndex(j, :pulley_side), EIndex(j, :pulley_at_src)])
            append!(values, [role.pulley_side, role.pulley_at_src ? 1.0 : 0.0])
        else
            append!(indices, [EIndex(j, :n_segs), EIndex(j, :tension_sign_src),
                              EIndex(j, :tension_sign_dst)])
            append!(values, SAM.SimFloat[role.n_segs, role.tension_sign_src,
                                         role.tension_sign_dst])
        end
    end
    isempty(indices) || setp(nw, indices)(param, values)
    return nothing
end

# ======================= state scatter ======================= #

"""
    NetworkWingGeom

The ref-point indices a PARTICLE_DYNAMICS wing needs so the state getter can
reconstruct its kinematic state from the struct points (a KINEMATIC wing is fitted,
not integrated): the four frame ref points (`zp1`/`zp2`/`yp1`/`yp2`), the origin, and
the aero points that need per-point `va_b`.
"""
struct NetworkWingGeom
    wing::Any
    zp1::Int; zp2::Int; yp1::Int; yp2::Int; origin::Int
    aero_points::Vector{Int}
    live_vertex::Int
end

"""
    network_wing_geoms(ss)

One [`NetworkWingGeom`](@ref) per PARTICLE_DYNAMICS wing, resolving its ref-point
indices ([`ref_single_id`](@ref)) and aero points (its `is_wing_node` points) once at
build so the state getter can call [`SAM.wing_kinematics_from_points!`](@ref).
"""
function network_wing_geoms(ss)
    geoms = NetworkWingGeom[]
    for wing in ss.wings
        wing.dynamics_type == SAM.PARTICLE_DYNAMICS || continue
        aero_pts = [i for (i, p) in enumerate(ss.points)
                    if p.is_wing_node && p.wing_idx == wing.idx]
        live_vertex = needs_live_aero(wing) && !isempty(aero_pts) ? aero_pts[1] : 0
        push!(geoms, NetworkWingGeom(wing,
            ref_single_id(wing.z_ref_points[1]), ref_single_id(wing.z_ref_points[2]),
            ref_single_id(wing.y_ref_points[1]), ref_single_id(wing.y_ref_points[2]),
            ref_single_id(wing.origin), aero_pts, live_vertex))
    end
    return geoms
end

"""
    NetworkStateGetter(nw, sys_struct, meta)

Callable that scatters the ND integrator state back into the `SystemStructure`,
mirroring the monolith's `get_all_state`: each `DYNAMIC` point's `pos_w`/`vel_w`,
each pulley's `len`/`vel`, each winch's `vel`/`force`/`set_value` and its tethers'
`len`, and each PARTICLE_DYNAMICS wing's reconstructed kinematics
([`SAM.wing_kinematics_from_points!`](@ref)) so `refresh_aero!` sees fresh apparent
wind — the network's stand-in for the monolith's symbolic `va_point_b`.
"""
struct NetworkStateGetter{NW}
    nw::NW
    dyn_idxs::Vector{Tuple{Int, Int}}
    pulley_idxs::Vector{Tuple{Int, Int}}
    winch_idxs::Vector{Tuple{Int, Int}}
    winch_tethers::Dict{Int, Vector{Int}}
    wing_geoms::Vector{NetworkWingGeom}
    wing_aero_readers::Vector{Any}
    body_readers::Vector{Any}
    body_static::Vector{Tuple{Int, Int}}
end

"""
    body_output_reader(nw, body, vertex)

Build the getter tuple `(body, vertex, pos, vel, R, omega, com, com_vel)` that reads a
body vertex's wide pose outputs (`pos`, `vel`, `pose_R`, `pose_omega`, `pose_com`,
`pose_com_vel`) off `vertex`, so [`NetworkStateGetter`](@ref) can scatter them back into
the struct.
"""
function body_output_reader(nw, body, vertex)
    getvec(sym) = getu(nw, [VIndex(vertex, Symbol(sym, :_, c)) for c in 1:3])
    R = getu(nw, [VIndex(vertex, Symbol(:pose_R_, k)) for k in 1:9])
    return (body, vertex, getvec(:pos), getvec(:vel), R, getvec(:pose_omega),
            getvec(:pose_com), getvec(:pose_com_vel))
end

"""
    wing_aero_reader(nw, wing, vertex)

Build the observed-variable getter tuple `(wing, force_getu, moment_getu)` that reads
a live wing's aggregate body-frame aero force/moment (`wing_aero_force_b_1..3`,
`wing_aero_moment_b_1..3`, exposed by [`SAM.LiveAeroWingNodePoint`](@ref)) off `vertex`.
"""
function wing_aero_reader(nw, wing, vertex)
    fsyms = [VIndex(vertex, Symbol(:wing_aero_force_b_, c)) for c in 1:3]
    msyms = [VIndex(vertex, Symbol(:wing_aero_moment_b_, c)) for c in 1:3]
    return (wing, getu(nw, fsyms), getu(nw, msyms))
end

function NetworkStateGetter(nw, ss, meta)
    vop = hasproperty(meta, :vertex_of_point) ? meta.vertex_of_point : nothing
    vmap = i -> vop === nothing ? i : vop[i]
    dyn_idxs = [(vmap(i), i) for (i, p) in enumerate(ss.points) if p.type == SAM.DYNAMIC]
    pulley_idxs = [(vmap(pulley_point_idx(ss, p)), p.idx) for p in ss.pulleys]
    winch_idxs = [(vmap(w.winch_point_idx), w.idx) for w in ss.winches]
    geoms = network_wing_geoms(ss)
    readers = Any[wing_aero_reader(nw, wg.wing, wg.live_vertex)
                  for wg in geoms if wg.live_vertex != 0]
    body_vertices = hasproperty(meta, :body_vertices) ? meta.body_vertices :
        [(v, b) for (v, b) in enumerate(meta.body_idxs)]
    body_readers = Any[body_output_reader(nw, ss.bodies[bidx], vertex)
                       for (vertex, bidx) in body_vertices]
    body_static = hasproperty(meta, :body_static) ? meta.body_static :
        Tuple{Int, Int}[]
    return NetworkStateGetter(nw, dyn_idxs, pulley_idxs, winch_idxs,
                              meta.winch_tethers, geoms, readers, body_readers,
                              body_static)
end

function (g::NetworkStateGetter)(integ, ss)
    s = NWState(g.nw, integ.u, integ.p)
    points = ss.points
    for (vi, i) in g.dyn_idxs
        point = points[i]
        for k in 1:3
            point.pos_w[k] = s.v[vi, Symbol(:pos_, k)]
            point.vel_w[k] = s.v[vi, Symbol(:vel_, k)]
        end
    end
    for wg in g.wing_geoms
        SAM.wing_kinematics_from_points!(wg.wing, points, ss.set, ss.am;
            zp1 = wg.zp1, zp2 = wg.zp2, yp1 = wg.yp1, yp2 = wg.yp2,
            origin = wg.origin, aero_points = wg.aero_points)
    end
    for (wing, force_getu, moment_getu) in g.wing_aero_readers
        wing.aero_force_b .= force_getu(integ)
        wing.aero_moment_b .= moment_getu(integ)
    end
    wind_factor = SAM.WindFactor(ss.am, ss.set.profile_law)
    for (body, vertex, pos_g, vel_g, R_g, omega_g, com_g, comvel_g) in g.body_readers
        body.pos_w .= pos_g(integ)
        body.vel_w .= vel_g(integ)
        body.com_w .= com_g(integ)
        body.com_vel .= comvel_g(integ)
        R_b_to_w = reshape(collect(R_g(integ)), 3, 3)
        body.Q_b_to_w .= SAM.rotation_matrix_to_quaternion(R_b_to_w)
        body.ω_b .= R_b_to_w' * omega_g(integ)
        for k in 1:4
            body.Q_p_to_w[k] = s.v[vertex, Symbol(:Q_, k)]
        end
        for k in 1:3
            body.ω_p[k] = s.v[vertex, Symbol(:omega_p_, k)]
        end
        if SAM.is_wing(body)
            body.R_b_to_w .= R_b_to_w
            body.v_wind .= wind_factor(body.pos_w[3]) .* ss.set.wind_vec
            body.va_b .= R_b_to_w' * (body.v_wind .- body.vel_w .+ body.wind_disturb)
        end
    end
    for (ride_idx, body_idx) in g.body_static
        point = points[ride_idx]
        body = ss.bodies[body_idx]
        R_b_to_w = SAM.quaternion_to_rotation_matrix(body.Q_b_to_w)
        anchor_w = body.pos_w .+ R_b_to_w * point.anchor_b
        ω_w = R_b_to_w * body.ω_b
        point.pos_w .= anchor_w
        point.vel_w .= body.com_vel .+ (ω_w × (anchor_w .- body.com_w))
    end
    for (vi, pidx) in g.pulley_idxs
        pulley = ss.pulleys[pidx]
        pulley.len = s.v[vi, :pulley_len]
        pulley.vel = s.v[vi, :pulley_vel]
    end
    for (vi, widx) in g.winch_idxs
        winch = ss.winches[widx]
        winch.vel = s.v[vi, :winch_vel]
        fill!(winch.force, 0.0)
        winch.force[1] = s.v[vi, :winch_force]
        winch.set_value = s.p.v[vi, :set_value]
        for (pos, tidx) in enumerate(g.winch_tethers[vi])
            ss.tethers[tidx].len = s.v[vi, Symbol(:tether_len_, pos)]
        end
    end
    return nothing
end

# ======================= control setter ======================= #

"""
    NetworkControlSetter(nw, ss, meta)

Callable `(integ, values)` writing each winch's control `set_value` into the ND
parameter vector, mirroring the monolith `set_set_values`. `values[winch.idx]` is
the setpoint for the winch; the winch vertex index is mapped through
`meta.vertex_of_point` (identity in a point-only network).
"""
struct NetworkControlSetter{S}
    setter::S
    order::Vector{Int}
end

function NetworkControlSetter(nw, ss, meta)
    isempty(ss.winches) && return nothing
    vop = hasproperty(meta, :vertex_of_point) ? meta.vertex_of_point : nothing
    vmap = i -> vop === nothing ? i : vop[i]
    idxs = [VIndex(vmap(w.winch_point_idx), :set_value) for w in ss.winches]
    order = [w.idx for w in ss.winches]
    return NetworkControlSetter(setp(nw, idxs), order)
end

function (c::NetworkControlSetter)(integ, values)
    c.setter(integ, [values[i] for i in c.order])
    return nothing
end

# ======================= backend hooks ======================= #

function SAM.build_prob!(::SAM.NetworkBackend, sam; prn = true)
    nw, u0, p0, meta = build_network(sam)
    dt = SAM.SimFloat(1 / sam.set.sample_freq)
    prob = ODEProblem(nw, u0, (0.0, dt), p0)
    SAM.sync_params!(meta.param_sync, prob, sam.sys_struct)
    getter = NetworkStateGetter(nw, sam.sys_struct, meta)
    setter = NetworkControlSetter(nw, sam.sys_struct, meta)
    sam.prob = SAM.ProbWithAttributes(; prob,
        param_sync = meta.param_sync, initial_sync = nothing,
        set_set_values = setter, get_set_values = nothing,
        get_aero_input = nothing, get_all_state = getter)
    return true
end

function SAM.init_backend!(::SAM.NetworkBackend, sam, solver;
        adaptive = true, prn = true, reinit_sys = true, reset_vel = true,
        ignore_l0 = false, apply_tether_lengths = true, remake_vsm = true,
        reset_integrator = true, vsm_min_wind = 0.5, lin_vsm = true)
    if reinit_sys
        SAM.reinit!(sam.sys_struct, sam.set;
            ignore_l0, remake_vsm, reset_vel, apply_tether_lengths, prn)
    end
    SAM.build_prob!(SAM.NetworkBackend(), sam; prn)
    integrator, _ = SAM.reinit!(sam, sam.prob, solver;
        adaptive, reset_integrator, lin_vsm, vsm_min_wind, prn)
    return integrator
end

end # module
