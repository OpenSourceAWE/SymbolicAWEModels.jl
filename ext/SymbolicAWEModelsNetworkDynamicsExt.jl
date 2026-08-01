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
endpoint loads (the monolith's `Flow` convention negates them). Variables are
scalarized to match the ND MTK-integration convention (bare `[input]`/`[output]`
scalar names, not `@connector` variables).

## Uniform interface

NetworkDynamics fixes one vertex-output width and one edge-output width for the
whole network. Every vertex therefore outputs `pos1..3, vel1..3, pulley_len_out`
(plain points emit `pulley_len_out = 0`) and every edge outputs, per endpoint,
`force1..3, mass, tension` — `tension` being a role-signed scalar spring force a
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
using SymbolicIndexingInterface: setp

const SAM = SymbolicAWEModels

# ======================= vertex Systems ======================= #

"""
    dynamic_point_dynamics(s, pos, vel, force, mass, pars)

Shared body of the DYNAMIC point/pulley vertices: `D(pos)=vel` and
`D(vel)=point_acceleration(...)` from the shared kernel, reading the vertex's
drag/damping/wind parameters `pars` (a named tuple of the scalar parameters).
"""
function dynamic_point_dynamics(s, pos, vel, force, mass, pars)
    accel = SAM.point_acceleration(s, pos, vel, force, mass, pars.drag_coeff,
        pars.area, [pars.world_damping1, pars.world_damping2, pars.world_damping3],
        [pars.wind_gnd1, pars.wind_gnd2, pars.wind_gnd3])
    return [
        D(pos[1]) ~ vel[1], D(pos[2]) ~ vel[2], D(pos[3]) ~ vel[3],
        D(vel[1]) ~ accel[1], D(vel[2]) ~ accel[2], D(vel[3]) ~ accel[3],
    ]
end

"""
    vertex_io()

The uniform vertex input/output variables (states `pos`/`vel`, aggregated edge
inputs `force_in`/`mass_in`/`tension_in`, and the `pulley_len_out` output). Returns
`(vars, pos, vel, force, mass_in, tension_in)`.
"""
function vertex_io()
    vars = @variables begin
        pos1(t); pos2(t); pos3(t); vel1(t); vel2(t); vel3(t)
        pulley_len_out(t)
        force_in1(t), [input = true]
        force_in2(t), [input = true]
        force_in3(t), [input = true]
        mass_in(t), [input = true]
        tension_in(t), [input = true]
    end
    return vars, [pos1, pos2, pos3], [vel1, vel2, vel3],
           [force_in1, force_in2, force_in3], mass_in, tension_in, pulley_len_out
end

"""
    network_dynamic_point(s; name)

Particle vertex `System`: states `pos1..3`/`vel1..3`, aggregated edge inputs, and
`pulley_len_out = 0`. `D(vel)` comes from the shared `point_acceleration`.
"""
function network_dynamic_point(s; name)
    vars, pos, vel, force, mass_in, _, pulley_len_out = vertex_io()
    pars = @parameters begin
        extra_mass = 0.0
        drag_coeff = 0.0
        area = 0.0
        world_damping1 = 0.0
        world_damping2 = 0.0
        world_damping3 = 0.0
        wind_gnd1 = 0.0
        wind_gnd2 = 0.0
        wind_gnd3 = 0.0
    end
    pnt = (; drag_coeff, area, world_damping1, world_damping2, world_damping3,
           wind_gnd1, wind_gnd2, wind_gnd3)
    eqs = [dynamic_point_dynamics(s, pos, vel, force, extra_mass + mass_in, pnt)
           pulley_len_out ~ 0.0]
    return System(eqs, t, vars, pars; name)
end

"""
    network_static_point(s; name)

Ground-anchored vertex `System`: `pos = pos_w`, `vel = 0`, `pulley_len_out = 0`.
Declares (and ignores) the aggregated edge inputs so edges may deliver to it.
"""
function network_static_point(s; name)
    vars, pos, vel, _, _, _, pulley_len_out = vertex_io()
    pars = @parameters pos_w1 = 0.0 pos_w2 = 0.0 pos_w3 = 0.0
    eqs = [
        pos[1] ~ pos_w1, pos[2] ~ pos_w2, pos[3] ~ pos_w3,
        vel[1] ~ 0.0, vel[2] ~ 0.0, vel[3] ~ 0.0,
        pulley_len_out ~ 0.0,
    ]
    return System(eqs, t, vars, pars; name)
end

"""
    network_pulley_point(s; name)

Dynamic pulley vertex: a particle (`pos`/`vel`) that additionally owns the pulley
rope split `pulley_len`/`pulley_vel`. The aggregated `tension_in` is the imbalance
`spring[seg1] − spring[seg2]` between its two incident segments (each emits its
role-signed spring force), driving `D(pulley_vel) = tension_in/mass − damp·vel`.
`pulley_len_out` exposes the split so the incident segments read it as their `l0`.
"""
function network_pulley_point(s; name)
    vars, pos, vel, force, mass_in, tension_in, pulley_len_out = vertex_io()
    extra = @variables pulley_len(t) pulley_vel(t)
    append!(vars, extra)
    pars = @parameters begin
        extra_mass = 0.0
        drag_coeff = 0.0
        area = 0.0
        world_damping1 = 0.0
        world_damping2 = 0.0
        world_damping3 = 0.0
        wind_gnd1 = 0.0
        wind_gnd2 = 0.0
        wind_gnd3 = 0.0
        pulley_mass = 1.0
        pulley_damp = 5.0
    end
    pnt = (; drag_coeff, area, world_damping1, world_damping2, world_damping3,
           wind_gnd1, wind_gnd2, wind_gnd3)
    eqs = [
        dynamic_point_dynamics(s, pos, vel, force, extra_mass + mass_in, pnt)
        D(pulley_len) ~ pulley_vel
        D(pulley_vel) ~ tension_in / pulley_mass - pulley_damp * pulley_vel
        pulley_len_out ~ pulley_len
    ]
    return System(eqs, t, vars, pars; name)
end

"""
    network_winch_point(s, winch, winch_point; name)

Reeling winch vertex at a `STATIC` winch point. Owns the motor speed `winch_vel`
and one `tether_len` state per connected tether (`tether_len_1`, …). The winch
motor law is `winch_component(winch.model, …)` reused verbatim with a fresh
`ParamView` (drum parameters baked as defaults); it reads the summed tether
tension `smooth_norm(tension_in)`, the mean tether length, the control `set_value`
and the `brake`, and returns the drum acceleration `acc`. Integrates
`D(winch_vel) = brake·0 + acc` and `D(tether_len_k) = brake·0 + winch_vel`; each
`tether_len_k` is read by that tether's segments through an `extin`.
"""
function network_winch_point(s, winch, winch_point; name)
    winch_point.type == SAM.STATIC || error(
        "NetworkBackend: winch $(winch.name) is at a non-STATIC point; only " *
        "STATIC winch points are supported so far.")
    n_tethers = length(winch.tether_idxs)
    vars, pos, vel, _, _, tension_in, pulley_len_out = vertex_io()
    winch_vel, winch_force = @variables winch_vel(t) winch_force(t)
    tether_lens = map(1:n_tethers) do k
        nm = Symbol(:tether_len_, k)
        only(@variables $nm(t))
    end
    append!(vars, [winch_vel, winch_force])
    append!(vars, tether_lens)
    pars = @parameters begin
        pos_w1 = 0.0
        pos_w2 = 0.0
        pos_w3 = 0.0
        set_value = 0.0
        brake = 0.0
        speed_controlled = 0.0
    end

    view = SAM.ParamView(SAM.ParamRegistry(s.sys_struct))
    motor = SAM.winch_component(winch.model, s.sys_struct, winch.idx;
                                name = :motor, params = view)
    SAM.validate_winch_component(motor, winch)
    winch_acc = ifelse(speed_controlled > 0.5, 0.0, motor.acc)

    eqs = [
        pos[1] ~ pos_w1, pos[2] ~ pos_w2, pos[3] ~ pos_w3,
        vel[1] ~ 0.0, vel[2] ~ 0.0, vel[3] ~ 0.0,
        pulley_len_out ~ 0.0,
        winch_force ~ SAM.smooth_norm(tension_in),
        motor.vel ~ winch_vel,
        motor.len ~ sum(tether_lens) / n_tethers,
        motor.force ~ winch_force,
        motor.set_value ~ set_value,
        motor.brake ~ brake,
        D(winch_vel) ~ ifelse(brake > 0.5, 0.0, winch_acc),
    ]
    for tl in tether_lens
        eqs = [eqs; D(tl) ~ ifelse(brake > 0.5, 0.0, winch_vel)]
    end
    return System(eqs, t, vars, pars; name, systems = [motor])
end

# ======================= edge Systems ======================= #

"""
    edge_io()

The uniform edge input/output variables: the src/dst vertex outputs read as inputs
(`src_pos`/`src_vel`/`src_pulley_len`, likewise `dst_…`) and the endpoint loads
written as outputs (`src_force`/`src_mass`/`src_tension`, likewise `dst_…`).
Returns `(vars, src_pos, src_vel, src_pulley_len, dst_pos, dst_vel,
dst_pulley_len, src_force, src_mass, src_tension, dst_force, dst_mass,
dst_tension)`.
"""
function edge_io()
    vars = @variables begin
        src_pos1(t), [input = true]
        src_pos2(t), [input = true]
        src_pos3(t), [input = true]
        src_vel1(t), [input = true]
        src_vel2(t), [input = true]
        src_vel3(t), [input = true]
        src_pulley_len(t), [input = true]
        dst_pos1(t), [input = true]
        dst_pos2(t), [input = true]
        dst_pos3(t), [input = true]
        dst_vel1(t), [input = true]
        dst_vel2(t), [input = true]
        dst_vel3(t), [input = true]
        dst_pulley_len(t), [input = true]
        src_force1(t), [output = true]
        src_force2(t), [output = true]
        src_force3(t), [output = true]
        src_mass(t), [output = true]
        src_tension(t), [output = true]
        dst_force1(t), [output = true]
        dst_force2(t), [output = true]
        dst_force3(t), [output = true]
        dst_mass(t), [output = true]
        dst_tension(t), [output = true]
    end
    return vars,
        [src_pos1, src_pos2, src_pos3], [src_vel1, src_vel2, src_vel3], src_pulley_len,
        [dst_pos1, dst_pos2, dst_pos3], [dst_vel1, dst_vel2, dst_vel3], dst_pulley_len,
        [src_force1, src_force2, src_force3], src_mass, src_tension,
        [dst_force1, dst_force2, dst_force3], dst_mass, dst_tension
end

"""
    edge_load_eqs(s, l0, spring_pars, wind, io; src_tension_val, dst_tension_val)

Assemble the endpoint-load equations shared by every edge type from the scalarized
`io` tuple: compute the positive endpoint forces, half-masses and scalar spring
tension from `segment_endpoint_loads`, and bind the outputs. `src_tension_val` /
`dst_tension_val` are the role-signed tensions to emit to each endpoint (`0` unless
that endpoint is a pulley/winch reader).
"""
function edge_load_eqs(s, l0, spring_pars, wind, io, src_tension_val, dst_tension_val)
    (_, src_pos, src_vel, _, dst_pos, dst_vel, _,
     src_force, src_mass, src_tension, dst_force, dst_mass, dst_tension) = io
    fsrc, fdst, half, _ = SAM.segment_endpoint_loads(
        s, src_pos, src_vel, dst_pos, dst_vel, spring_pars.unit_stiffness,
        spring_pars.unit_damping, spring_pars.compression_frac, l0,
        spring_pars.diameter, spring_pars.density, spring_pars.cd_tether, wind)
    return [
        src_force[1] ~ fsrc[1], src_force[2] ~ fsrc[2], src_force[3] ~ fsrc[3],
        src_mass ~ half, src_tension ~ src_tension_val,
        dst_force[1] ~ fdst[1], dst_force[2] ~ fdst[2], dst_force[3] ~ fdst[3],
        dst_mass ~ half, dst_tension ~ dst_tension_val,
    ]
end

"""
    spring_parameters()

The spring-damper + tether-drag parameters common to all edge types, plus the
ground wind. Returns `(pars, spring_pars_named, wind_vec)`.
"""
function spring_parameters()
    pars = @parameters begin
        unit_stiffness = 0.0
        unit_damping = 0.0
        compression_frac = 0.1
        diameter = 0.0
        density = 0.0
        cd_tether = 1.0
        wind_gnd1 = 0.0
        wind_gnd2 = 0.0
        wind_gnd3 = 0.0
    end
    spring = (; unit_stiffness, unit_damping, compression_frac, diameter,
              density, cd_tether)
    return pars, spring, [wind_gnd1, wind_gnd2, wind_gnd3]
end

"""
    network_segment(s; name)

Plain spring-damper edge with a frozen rest length `l0` (parameter). Emits zero
tension (neither endpoint is a pulley/winch reader).
"""
function network_segment(s; name)
    io = edge_io()
    vars = io[1]
    spars, spring, wind = spring_parameters()
    l0 = only(@parameters l0 = 1.0)
    eqs = edge_load_eqs(s, l0, spring, wind, io, 0.0, 0.0)
    return System(eqs, t, vars, [spars; l0]; name)
end

"""
    network_pulley_segment(s; name)

Pulley spring-damper edge. `l0` is driven by the pulley vertex's `pulley_len`
output (read from whichever endpoint is the pulley, `pulley_at_src`): the first
pulley segment gets `l0 = pulley_len`, the second `l0 = sum_len − pulley_len`
(`pulley_side = ±1`). It emits its role-signed spring tension `pulley_side·spring`
to the pulley endpoint so the pulley aggregates `spring[seg1] − spring[seg2]`.
"""
function network_pulley_segment(s; name)
    io = edge_io()
    vars, _, _, src_pulley_len, _, _, dst_pulley_len = io
    spars, spring, wind = spring_parameters()
    extra = @parameters pulley_sum_len = 1.0 pulley_side = 1.0 pulley_at_src = 0.0
    endpoint_len = ifelse(pulley_at_src > 0.5, src_pulley_len, dst_pulley_len)
    l0 = ifelse(pulley_side > 0.0, endpoint_len, pulley_sum_len - endpoint_len)
    fsrc, fdst, half, spring_scalar = SAM.segment_endpoint_loads(
        s, io[2], io[3], io[5], io[6], spring.unit_stiffness, spring.unit_damping,
        spring.compression_frac, l0, spring.diameter, spring.density,
        spring.cd_tether, wind)
    src_tension_val = ifelse(pulley_at_src > 0.5, pulley_side * spring_scalar, 0.0)
    dst_tension_val = ifelse(pulley_at_src > 0.5, 0.0, pulley_side * spring_scalar)
    (_, _, _, _, _, _, _, src_force, src_mass, src_tension,
     dst_force, dst_mass, dst_tension) = io
    eqs = [
        src_force[1] ~ fsrc[1], src_force[2] ~ fsrc[2], src_force[3] ~ fsrc[3],
        src_mass ~ half, src_tension ~ src_tension_val,
        dst_force[1] ~ fdst[1], dst_force[2] ~ fdst[2], dst_force[3] ~ fdst[3],
        dst_mass ~ half, dst_tension ~ dst_tension_val,
    ]
    return System(eqs, t, vars, [spars; extra]; name)
end

"""
    network_tether_segment(s; name)

Winched-tether spring-damper edge. Its rest length is `l0 = tether_len / n_segs`,
where `tether_len` arrives as a NetworkDynamics external input (`tether_len_ext`)
read from the winch vertex. The tether segment incident to the winch point emits
`+spring` there (`tension_sign_src`/`tension_sign_dst = 1`) so the winch reads the
tension; every other tether segment emits zero.
"""
function network_tether_segment(s; name)
    io = edge_io()
    vars = io[1]
    ext = @variables tether_len_ext(t) [input = true]
    tether_len_ext = only(ext)
    push!(vars, tether_len_ext)
    spars, spring, wind = spring_parameters()
    extra = @parameters n_segs = 1.0 tension_sign_src = 0.0 tension_sign_dst = 0.0
    l0 = tether_len_ext / n_segs
    (_, _, _, _, _, _, _, src_force, src_mass, src_tension,
     dst_force, dst_mass, dst_tension) = io
    fsrc, fdst, half, spring_scalar = SAM.segment_endpoint_loads(
        s, io[2], io[3], io[5], io[6], spring.unit_stiffness, spring.unit_damping,
        spring.compression_frac, l0, spring.diameter, spring.density,
        spring.cd_tether, wind)
    eqs = [
        src_force[1] ~ fsrc[1], src_force[2] ~ fsrc[2], src_force[3] ~ fsrc[3],
        src_mass ~ half, src_tension ~ tension_sign_src * spring_scalar,
        dst_force[1] ~ fdst[1], dst_force[2] ~ fdst[2], dst_force[3] ~ fdst[3],
        dst_mass ~ half, dst_tension ~ tension_sign_dst * spring_scalar,
    ]
    return System(eqs, t, vars, [spars; extra]; name)
end

# ======================= I/O symbol lists ======================= #

const VERTEX_INPUTS = [:force_in1, :force_in2, :force_in3, :mass_in, :tension_in]
const VERTEX_OUTPUTS = [:pos1, :pos2, :pos3, :vel1, :vel2, :vel3, :pulley_len_out]
const EDGE_SRC_IN = [:src_pos1, :src_pos2, :src_pos3, :src_vel1, :src_vel2,
                     :src_vel3, :src_pulley_len]
const EDGE_DST_IN = [:dst_pos1, :dst_pos2, :dst_pos3, :dst_vel1, :dst_vel2,
                     :dst_vel3, :dst_pulley_len]
const EDGE_SRC_OUT = [:src_force1, :src_force2, :src_force3, :src_mass, :src_tension]
const EDGE_DST_OUT = [:dst_force1, :dst_force2, :dst_force3, :dst_mass, :dst_tension]

# ======================= helpers ======================= #

"""
    ground_wind(sam)

Ground-level wind vector, matching the monolith's `wind_vec_gnd` (`scalar_eqs!`):
`set.wind_vec` verbatim, with a tiny +x fallback when it is exactly zero (avoids
a normalize-by-zero in the apparent-wind direction).
"""
function ground_wind(sam)
    wind_vec = collect(Float64, sam.set.wind_vec)
    sum(abs2, wind_vec) < 1e-20 && return [1e-10, 0.0, 0.0]
    return wind_vec
end

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
`:pulley` or `:tether`; the remaining fields carry the pulley split / winched-tether
data an edge needs (`0`/`nothing` when not applicable).
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
pulley's two segments, `:tether` if it belongs to a winched tether, else `:plain`.
The tether segment incident to the winch point is marked to emit `+spring` there.
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
            push!(roles, SegmentRoles(:plain, 0, 0.0, false, 0, 0, 0.0, 0.0))
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
    pulley_line_mass(ss, pulley)

Total rope mass carried by a pulley: `sum_len · ρ · π (d/2)²` of its first
segment, matching `pulley_eqs!`.
"""
function pulley_line_mass(ss, pulley)
    seg = ss.segments[pulley.segment_idxs[1]]
    return pulley.sum_len * seg.density * π * (seg.diameter / 2)^2
end

# ======================= network assembly ======================= #

"""
    build_network(sam)

Translate `sam.sys_struct` into a NetworkDynamics `Network` with per-instance
parameters and initial state. Returns `(nw, u0, p0, meta)` with flat state/parameter
vectors and a `meta` named tuple of index maps the getter/control setters need.
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

    dyn = VertexModel(network_dynamic_point(sam; name = :dyn),
        VERTEX_INPUTS, VERTEX_OUTPUTS; mtkcompile = true, name = :dyn)
    stat = VertexModel(network_static_point(sam; name = :stat),
        VERTEX_INPUTS, VERTEX_OUTPUTS; mtkcompile = true, name = :stat)
    pulley_v = isempty(ss.pulleys) ? nothing :
        VertexModel(network_pulley_point(sam; name = :pul),
            VERTEX_INPUTS, VERTEX_OUTPUTS; mtkcompile = true, name = :pul)
    winch_v = Dict{Int, Any}()
    for winch in ss.winches
        winch_v[winch.winch_point_idx] = VertexModel(
            network_winch_point(sam, winch, points[winch.winch_point_idx];
                                name = :wch),
            VERTEX_INPUTS, VERTEX_OUTPUTS; mtkcompile = true, name = :wch)
    end

    vmodels = Vector{VertexModel}(undef, n)
    for i in 1:n
        if haskey(winch_of_point, i)
            vmodels[i] = winch_v[i]
        elseif haskey(pulley_of_point, i)
            vmodels[i] = pulley_v
        elseif points[i].type == SAM.STATIC
            vmodels[i] = stat
        else
            vmodels[i] = dyn
        end
    end

    roles = classify_segments(ss)
    role_of_seg = Dict(segments[k].idx => roles[k] for k in eachindex(segments))
    plain_edge = EdgeModel(network_segment(sam; name = :seg),
        EDGE_SRC_IN, EDGE_DST_IN, EDGE_SRC_OUT, EDGE_DST_OUT;
        mtkcompile = true, name = :seg)
    pulley_edge = isempty(ss.pulleys) ? nothing :
        EdgeModel(network_pulley_segment(sam; name = :pseg),
            EDGE_SRC_IN, EDGE_DST_IN, EDGE_SRC_OUT, EDGE_DST_OUT;
            mtkcompile = true, name = :pseg)
    tether_edge_of_tether = build_tether_edges(sam, ss, winch_of_point)

    edgelist = collect(edges(graph))
    emodels = Vector{EdgeModel}(undef, length(edgelist))
    for (j, e) in enumerate(edgelist)
        seg = seg_of[minmax(src(e), dst(e))]
        role = role_of_seg[seg.idx]
        if role.kind == :pulley
            emodels[j] = pulley_edge
        elseif role.kind == :tether
            emodels[j] = tether_edge_of_tether[role.tether_idx]
        else
            emodels[j] = plain_edge
        end
    end
    nw = Network(graph, vmodels, emodels)

    write_total_mass!(ss)
    param, state = NWParameter(nw), NWState(nw)
    wind = ground_wind(sam)
    set_vertex_data!(sam, ss, param, state, wind, winch_of_point, pulley_of_point)
    set_edge_data!(sam, ss, param, edgelist, seg_of, role_of_seg, wind)

    meta = (; winch_of_point, pulley_of_point,
            winch_tethers = Dict(w.winch_point_idx => collect(w.tether_idxs)
                                 for w in ss.winches))
    return nw, uflat(state), pflat(param), meta
end

"""
    build_tether_edges(sam, ss, winch_of_point)

One `network_tether_segment` `EdgeModel` per winched tether, its rest length wired
to the winch vertex's matching `tether_len_k` state via `extin`. The kernel is
compiled once and rebound per tether with `EdgeModel(base; extin=…)` so only the
external-input index differs.
"""
function build_tether_edges(sam, ss, winch_of_point)
    edges_of = Dict{Int, Any}()
    isempty(ss.winches) && return edges_of
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
            base = EdgeModel(network_tether_segment(sam; name = :tseg),
                EDGE_SRC_IN, EDGE_DST_IN, EDGE_SRC_OUT, EDGE_DST_OUT;
                extin, mtkcompile = true, name = :tseg)
            edges_of[tidx] = base
        else
            edges_of[tidx] = EdgeModel(base;
                extin = [VIndex(winch.winch_point_idx, sym)])
        end
    end
    return edges_of
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

"""
    set_vertex_data!(sam, ss, param, state, wind, winch_of_point, pulley_of_point)

Fill each vertex's per-instance parameters and initial state: drag/damping/wind
and `pos`/`vel` for particles; the pulley split `pulley_len`/`pulley_vel` and rope
mass for pulleys; `winch_vel`, `tether_len_k`, `set_value`, `brake` for winches.
"""
function set_vertex_data!(sam, ss, param, state, wind, winch_of_point, pulley_of_point)
    for (i, point) in enumerate(ss.points)
        if haskey(winch_of_point, i)
            winch = ss.winches[winch_of_point[i]]
            param.v[i, :pos_w1] = point.pos_w[1]
            param.v[i, :pos_w2] = point.pos_w[2]
            param.v[i, :pos_w3] = point.pos_w[3]
            param.v[i, :set_value] = winch.set_value
            param.v[i, :brake] = winch.brake
            param.v[i, :speed_controlled] = winch.speed_controlled ? 1.0 : 0.0
            state.v[i, :winch_vel] = winch.vel
            for (pos, tidx) in enumerate(winch.tether_idxs)
                state.v[i, Symbol(:tether_len_, pos)] = ss.tethers[tidx].len
            end
        elseif haskey(pulley_of_point, i)
            pulley = ss.pulleys[pulley_of_point[i]]
            set_particle_params!(param, state, i, point, wind)
            param.v[i, :pulley_mass] = pulley_line_mass(ss, pulley)
            state.v[i, :pulley_len] = pulley.len
            state.v[i, :pulley_vel] = pulley.vel
        elseif point.type == SAM.STATIC
            param.v[i, :pos_w1] = point.pos_w[1]
            param.v[i, :pos_w2] = point.pos_w[2]
            param.v[i, :pos_w3] = point.pos_w[3]
        else
            set_particle_params!(param, state, i, point, wind)
        end
    end
    return nothing
end

"""
    set_particle_params!(param, state, i, point, wind)

Set the shared DYNAMIC-particle parameters (mass/drag/damping/wind) and initial
`pos`/`vel` for vertex `i`.
"""
function set_particle_params!(param, state, i, point, wind)
    param.v[i, :extra_mass] = point.extra_mass
    param.v[i, :drag_coeff] = point.drag_coeff
    param.v[i, :area] = point.area
    damping = point.world_frame_damping
    param.v[i, :world_damping1] = damping[1]
    param.v[i, :world_damping2] = damping[2]
    param.v[i, :world_damping3] = damping[3]
    param.v[i, :wind_gnd1] = wind[1]
    param.v[i, :wind_gnd2] = wind[2]
    param.v[i, :wind_gnd3] = wind[3]
    state.v[i, :pos1] = point.pos_w[1]
    state.v[i, :pos2] = point.pos_w[2]
    state.v[i, :pos3] = point.pos_w[3]
    state.v[i, :vel1] = point.vel_w[1]
    state.v[i, :vel2] = point.vel_w[2]
    state.v[i, :vel3] = point.vel_w[3]
    return nothing
end

"""
    set_edge_data!(sam, ss, param, edgelist, seg_of, role_of_seg, wind)

Fill each edge's spring/drag parameters and role-specific data (frozen `l0`,
pulley split, or winched-tether `n_segs`/tension signs).
"""
function set_edge_data!(sam, ss, param, edgelist, seg_of, role_of_seg, wind)
    points = ss.points
    for (j, e) in enumerate(edgelist)
        seg = seg_of[minmax(src(e), dst(e))]
        role = role_of_seg[seg.idx]
        wing_structural = points[seg.point_idxs[1]].is_wing_node &&
                          points[seg.point_idxs[2]].is_wing_node
        param.e[j, :unit_stiffness] = segment_stiffness(seg)
        param.e[j, :unit_damping] = seg.unit_damping
        param.e[j, :compression_frac] = seg.compression_frac
        param.e[j, :diameter] = seg.diameter
        param.e[j, :density] = seg.density
        param.e[j, :cd_tether] = wing_structural ? 0.0 : sam.set.cd_tether
        param.e[j, :wind_gnd1] = wind[1]
        param.e[j, :wind_gnd2] = wind[2]
        param.e[j, :wind_gnd3] = wind[3]
        if role.kind == :pulley
            param.e[j, :pulley_sum_len] = ss.pulleys[role.pulley_idx].sum_len
            param.e[j, :pulley_side] = role.pulley_side
            param.e[j, :pulley_at_src] = role.pulley_at_src ? 1.0 : 0.0
        elseif role.kind == :tether
            param.e[j, :n_segs] = role.n_segs
            param.e[j, :tension_sign_src] = role.tension_sign_src
            param.e[j, :tension_sign_dst] = role.tension_sign_dst
        else
            param.e[j, :l0] = seg.l0
        end
    end
    return nothing
end

# ======================= state scatter ======================= #

"""
    NetworkStateGetter(nw, sys_struct, meta)

Callable that scatters the ND integrator state back into the `SystemStructure`,
mirroring the monolith's `get_all_state`: each `DYNAMIC` point's `pos_w`/`vel_w`,
each pulley's `len`/`vel`, and each winch's `vel`/`force`/`set_value` and its
tethers' `len`.
"""
struct NetworkStateGetter{NW}
    nw::NW
    dyn_idxs::Vector{Int}
    pulley_idxs::Vector{Tuple{Int, Int}}
    winch_idxs::Vector{Tuple{Int, Int}}
    winch_tethers::Dict{Int, Vector{Int}}
end

function NetworkStateGetter(nw, ss, meta)
    dyn_idxs = [i for (i, p) in enumerate(ss.points) if p.type == SAM.DYNAMIC]
    pulley_idxs = [(pulley_point_idx(ss, p), p.idx) for p in ss.pulleys]
    winch_idxs = [(w.winch_point_idx, w.idx) for w in ss.winches]
    return NetworkStateGetter(nw, dyn_idxs, pulley_idxs, winch_idxs,
                              meta.winch_tethers)
end

function (g::NetworkStateGetter)(integ, ss)
    s = NWState(g.nw, integ.u, integ.p)
    points = ss.points
    for i in g.dyn_idxs
        point = points[i]
        point.pos_w[1] = s.v[i, :pos1]
        point.pos_w[2] = s.v[i, :pos2]
        point.pos_w[3] = s.v[i, :pos3]
        point.vel_w[1] = s.v[i, :vel1]
        point.vel_w[2] = s.v[i, :vel2]
        point.vel_w[3] = s.v[i, :vel3]
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
    NetworkControlSetter(nw, winch_idxs)

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
    getter = NetworkStateGetter(nw, sam.sys_struct, meta)
    setter = NetworkControlSetter(nw, sam.sys_struct)
    sam.prob = SAM.ProbWithAttributes(; prob,
        param_sync = nothing, initial_sync = nothing,
        set_set_values = setter, get_set_values = nothing,
        get_aero_input = nothing, get_all_state = getter)
    return true
end

function SAM.init_backend!(::SAM.NetworkBackend, sam, solver;
        adaptive = true, prn = true, reinit_sys = true, reset_vel = true,
        ignore_l0 = false, apply_tether_lengths = true, remake_vsm = true,
        reset_integrator = true, vsm_min_wind = 0.5)
    if reinit_sys
        SAM.reinit!(sam.sys_struct, sam.set;
            ignore_l0, remake_vsm, reset_vel, apply_tether_lengths, prn)
    end
    SAM.build_prob!(SAM.NetworkBackend(), sam; prn)
    integrator, _ = SAM.reinit!(sam, sam.prob, solver;
        adaptive, reset_integrator, lin_vsm = false, vsm_min_wind, prn)
    return integrator
end

end # module
