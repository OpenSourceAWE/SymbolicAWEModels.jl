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
([`param_addr`](@ref)); array fields expand to one scalar `field_k` per component,
matching NetworkDynamics' parameter scalarization. Computed entries (non-`PathReader`,
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
            for k in 1:length(entry.read(ss))
                add_param!(builder, param_addr(edge, addr_index, Symbol(field, :_, k)),
                           SAM.PathReader((container, cont_index, field, k)))
            end
        end
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
network_dynamic_point(s, params, idx; name) = SAM.DynamicPoint(s, params, idx; name)

"""
    network_static_point(s, params, idx; name)

Ground-anchored vertex: the shared [`SAM.StaticPoint`](@ref) component wrapped as a
network vertex.
"""
network_static_point(s, params, idx; name) = SAM.StaticPoint(s, params, idx; name)

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
network_pulley_point(s, params, idx; name) =
    SAM.PulleyPoint(s, params, idx, pulley_rope_mass(params, idx); name)

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
function network_pulley_segment(s, params, idx; name)
    io = SAM.segment_io()
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
    return System(eqs, t, vars,
        [param_unknowns(params); pulley_sum_len; pulley_side; pulley_at_src]; name)
end

"""
    network_tether_segment(s, params, idx; name)

Winched-tether spring-damper edge. Its rest length is `l0 = tether_len / n_segs`,
where `tether_len` arrives as a NetworkDynamics external input (`tether_len_ext`)
read from the winch vertex. The tether segment incident to the winch point emits
`+spring` there (`tension_sign_src`/`tension_sign_dst = 1`) so the winch reads the
tension; every other tether segment emits zero.
"""
function network_tether_segment(s, params, idx; name)
    io = SAM.segment_io()
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
    return System(eqs, t, vars,
        [param_unknowns(params); n_segs; tension_sign_src; tension_sign_dst]; name)
end

# ======================= I/O symbol lists ======================= #
# Array-valued I/O symbols; the forked ND MTK integration scalarizes each to its
# `_1/_2/_3` components while preserving the vertex-output / edge-output ordering.

const VERTEX_INPUTS = [:force_in, :mass_in, :tension_in]
const VERTEX_OUTPUTS = [:pos, :vel, :pulley_len_out]
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
                                 for w in ss.winches))
    return nw, uflat(state), pflat(param), meta
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
            record_winch_params!(builder, i, winch_of_point[i])
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
    record_pos_w_params!(builder, i)

Record the anchored position `pos_w_k` for a static/winch vertex `i`.
"""
function record_pos_w_params!(builder, i)
    for k in 1:3
        add_param!(builder, VIndex(i, Symbol(:pos_w_, k)),
                   SAM.PathReader((:points, i, :pos_w, k)))
    end
    return nothing
end

"""
    record_winch_params!(builder, i, widx)

Record the winch vertex parameters for vertex `i` (winch `widx`): anchored position
plus `brake` and `speed_controlled`. `set_value` is deliberately excluded — it is a
control input owned by `set_set_values`/`get_set_values`, and syncing it here would
clobber the setpoint the control setter writes each step.
"""
function record_winch_params!(builder, i, widx)
    record_pos_w_params!(builder, i)
    add_param!(builder, VIndex(i, :brake), SAM.PathReader((:winches, widx, :brake)))
    add_param!(builder, VIndex(i, :speed_controlled),
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
    dyn_idxs::Vector{Int}
    pulley_idxs::Vector{Tuple{Int, Int}}
    winch_idxs::Vector{Tuple{Int, Int}}
    winch_tethers::Dict{Int, Vector{Int}}
    wing_geoms::Vector{NetworkWingGeom}
    wing_aero_readers::Vector{Any}
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
    dyn_idxs = [i for (i, p) in enumerate(ss.points) if p.type == SAM.DYNAMIC]
    pulley_idxs = [(pulley_point_idx(ss, p), p.idx) for p in ss.pulleys]
    winch_idxs = [(w.winch_point_idx, w.idx) for w in ss.winches]
    geoms = network_wing_geoms(ss)
    readers = Any[wing_aero_reader(nw, wg.wing, wg.live_vertex)
                  for wg in geoms if wg.live_vertex != 0]
    return NetworkStateGetter(nw, dyn_idxs, pulley_idxs, winch_idxs,
                              meta.winch_tethers, geoms, readers)
end

function (g::NetworkStateGetter)(integ, ss)
    s = NWState(g.nw, integ.u, integ.p)
    points = ss.points
    for i in g.dyn_idxs
        point = points[i]
        for k in 1:3
            point.pos_w[k] = s.v[i, Symbol(:pos_, k)]
            point.vel_w[k] = s.v[i, Symbol(:vel_, k)]
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
    NetworkControlSetter(nw, ss)

Callable `(integ, values)` writing each winch's control `set_value` into the ND
parameter vector, mirroring the monolith `set_set_values`. `values[winch.idx]` is
the setpoint for the winch at vertex `winch_idxs[k]`.
"""
struct NetworkControlSetter{S}
    setter::S
    order::Vector{Int}
end

function NetworkControlSetter(nw, ss)
    isempty(ss.winches) && return nothing
    idxs = [VIndex(w.winch_point_idx, :set_value) for w in ss.winches]
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
    setter = NetworkControlSetter(nw, sam.sys_struct)
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
