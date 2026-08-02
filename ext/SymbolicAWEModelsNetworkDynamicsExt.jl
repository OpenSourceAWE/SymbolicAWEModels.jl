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
using SymbolicIndexingInterface: setp
using LinearAlgebra: cross

const SAM = SymbolicAWEModels

# ======================= live struct-backed parameters ======================= #
# Per-type ND kernels carry generic parameters (one symbol per field, one value
# per instance in `NWParameter`). To keep the `SystemStructure` the single source
# of truth — so mutating a field takes effect live, as on the monolith — the flat
# parameter buffer is re-read from the struct every step by the monolith's own
# `SAM.ParamGroup`/`sync_group!`, using the same `PathReader` (getproperty) primitive
# and `setp` (setproperty) setter — only the setter indexes `VIndex`/`EIndex` onto the
# `Network` instead of bare symbols. Values that are not a plain field read
# (ground wind, tether drag, pulley rope mass, structural constants) get a small
# named reader below; each reader maps a `SystemStructure` to one scalar.

"""
    ground_wind(ss)

Ground-level wind vector, matching the monolith's `wind_vec_gnd` (`scalar_eqs!`):
`set.wind_vec` verbatim, with a tiny +x fallback when it is exactly zero (avoids
a normalize-by-zero in the apparent-wind direction).
"""
function ground_wind(ss)
    wind_vec = collect(Float64, ss.set.wind_vec)
    sum(abs2, wind_vec) < 1e-20 && return [1e-10, 0.0, 0.0]
    return wind_vec
end

"""
    GroundWindReader(index)

Reader returning component `index` of [`ground_wind`](@ref).
"""
struct GroundWindReader
    index::Int
end
(reader::GroundWindReader)(ss) = ground_wind(ss)[reader.index]

"""
    TetherDragReader(wing_structural)

Reader for a segment's tether drag coefficient: `0` for a wing-structural segment
(the wing's aero owns those loads), else `set.cd_tether`.
"""
struct TetherDragReader
    wing_structural::Bool
end
(reader::TetherDragReader)(ss) = reader.wing_structural ? 0.0 : ss.set.cd_tether

"""
    PulleyLineMassReader(pulley_idx)

Reader for a pulley's rope mass `sum_len · ρ · π (d/2)²` (its first segment),
matching `pulley_eqs!`.
"""
struct PulleyLineMassReader
    pulley_idx::Int
end
function (reader::PulleyLineMassReader)(ss)
    pulley = ss.pulleys[reader.pulley_idx]
    seg = ss.segments[pulley.segment_idxs[1]]
    return pulley.sum_len * seg.density * π * (seg.diameter / 2)^2
end

"""
    ConstReader(value)

Reader returning a fixed structural constant (pulley side, segment count, tension
sign) that never changes after assembly.
"""
struct ConstReader
    value::SAM.SimFloat
end
(reader::ConstReader)(ss) = reader.value

SAM.sync_params!(g::SAM.ParamGroup, target, ss::SAM.SystemStructure) =
    SAM.sync_group!(g, target, ss)

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
    build_network_param_sync(nw, builder)

Turn the recorded parameters into a `SAM.ParamGroup` (a `setp` over all recorded
indices plus the readers), or `nothing` when empty. It is the same group type the
monolith syncs, only with a `VIndex`/`EIndex` setter; applied to the problem in
`build_prob!` and re-applied every step through the shared `ProbWithAttributes`
machinery.
"""
function build_network_param_sync(nw, builder::ParamBuilder)
    isempty(builder.indices) && return nothing
    setter = setp(nw, builder.indices)
    buffer = Vector{SAM.SimFloat}(undef, length(builder.indices))
    return SAM.ParamGroup(setter, builder.readers, buffer)
end

# ======================= vertex Systems ======================= #

"""
    particle_params()

The shared DYNAMIC-particle parameters (mass/drag/damping/wind). Returns the flat
parameter vector (with `world_damping`/`wind_gnd` as `[1:3]` array parameters) and
a named tuple of the same symbols.
"""
function particle_params()
    extra_mass = SAM.make_param(:extra_mass, 0.0)
    drag_coeff = SAM.make_param(:drag_coeff, 0.0)
    area = SAM.make_param(:area, 0.0)
    world_damping = SAM.make_array_param(:world_damping, zeros(3))
    wind_gnd = SAM.make_array_param(:wind_gnd, zeros(3))
    pars = [extra_mass, drag_coeff, area, world_damping, wind_gnd]
    named = (; extra_mass, drag_coeff, area, world_damping, wind_gnd)
    return pars, named
end

"""
    dynamic_point_dynamics(s, pos, vel, force, mass, pars)

Shared body of the DYNAMIC point/pulley vertices: `D(pos)=vel` and
`D(vel)=point_acceleration(...)` from the shared kernel, reading the vertex's
drag/damping/wind parameters `pars` (a [`particle_params`](@ref) named tuple).
"""
function dynamic_point_dynamics(s, pos, vel, force, mass, pars)
    accel = SAM.point_acceleration(s, collect(pos), collect(vel), collect(force),
        mass, pars.drag_coeff, pars.area, collect(pars.world_damping),
        collect(pars.wind_gnd))
    return [D.(collect(pos)) .~ collect(vel);
            D.(collect(vel)) .~ accel]
end

"""
    vertex_io()

The uniform vertex input/output variables (states `pos`/`vel`, aggregated edge
inputs `force_in`/`mass_in`/`tension_in`, and the `pulley_len_out` output), all as
array-valued variables where they are 3-vectors. Returns `(vars, pos, vel,
force_in, mass_in, tension_in, pulley_len_out)`.
"""
function vertex_io()
    vars = @variables begin
        pos(t)[1:3]
        vel(t)[1:3]
        pulley_len_out(t)
        force_in(t)[1:3], [input = true]
        mass_in(t), [input = true]
        tension_in(t), [input = true]
    end
    return vars, pos, vel, force_in, mass_in, tension_in, pulley_len_out
end

"""
    network_dynamic_point(s; name)

Particle vertex `System`: states `pos`/`vel`, aggregated edge inputs, and
`pulley_len_out = 0`. `D(vel)` comes from the shared `point_acceleration`.
"""
function network_dynamic_point(s; name)
    vars, pos, vel, force, mass_in, _, pulley_len_out = vertex_io()
    pars, named = particle_params()
    eqs = [dynamic_point_dynamics(s, pos, vel, force, named.extra_mass + mass_in, named);
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
    pos_w = SAM.make_array_param(:pos_w, zeros(3))
    eqs = [
        collect(pos) .~ collect(pos_w);
        collect(vel) .~ zeros(3);
        pulley_len_out ~ 0.0;
    ]
    return System(eqs, t, vars, [pos_w]; name)
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
    pars, named = particle_params()
    pulley_mass = SAM.make_param(:pulley_mass, 1.0)
    pulley_damp = SAM.make_param(:pulley_damp, 5.0)
    eqs = [
        dynamic_point_dynamics(s, pos, vel, force, named.extra_mass + mass_in, named);
        D(pulley_len) ~ pulley_vel;
        D(pulley_vel) ~ tension_in / pulley_mass - pulley_damp * pulley_vel;
        pulley_len_out ~ pulley_len;
    ]
    return System(eqs, t, vars, [pars; pulley_mass; pulley_damp]; name)
end

"""
    wing_frame_columns(zp1, zp2, yp1, yp2)

The particle-wing body→world rotation columns fitted from the four structural ref
points, matching `wing_eqs.jl`: `z = normalize(zp2−zp1)`,
`x = normalize(normalize(yp2−yp1) × z)`, `y = z × x`. Returns `(xaxis, yaxis,
zaxis)`, each a 3-vector (the columns of `R_b_to_w`).
"""
function wing_frame_columns(zp1, zp2, yp1, yp2)
    zaxis = SAM.smooth_normalize(zp2 .- zp1)
    xaxis = SAM.smooth_normalize(cross(SAM.smooth_normalize(yp2 .- yp1), zaxis))
    yaxis = cross(zaxis, xaxis)
    return xaxis, yaxis, zaxis
end

"""
    BODY_DAMP_EXTIN

The 15 external-input symbols a body-damped point reads: the four ref-point
positions `zp1/zp2/yp1/yp2` (12) and the wing origin velocity `ovel` (3).
"""
const BODY_DAMP_EXTIN = [
    :zp1x, :zp1y, :zp1z, :zp2x, :zp2y, :zp2z,
    :yp1x, :yp1y, :yp1z, :yp2x, :yp2y, :yp2z, :ovx, :ovy, :ovz]

"""
    body_damp_inputs()

The 15 scalar external-input variables ([`BODY_DAMP_EXTIN`]) a body-damped point
reads from its wing's ref points. Returns `(ext, zp1, zp2, yp1, yp2, ovel)` with
each ref point as a 3-vector of the input variables.
"""
function body_damp_inputs()
    ext = @variables begin
        zp1x(t), [input = true]; zp1y(t), [input = true]; zp1z(t), [input = true]
        zp2x(t), [input = true]; zp2y(t), [input = true]; zp2z(t), [input = true]
        yp1x(t), [input = true]; yp1y(t), [input = true]; yp1z(t), [input = true]
        yp2x(t), [input = true]; yp2y(t), [input = true]; yp2z(t), [input = true]
        ovx(t), [input = true]; ovy(t), [input = true]; ovz(t), [input = true]
    end
    (zp1x, zp1y, zp1z, zp2x, zp2y, zp2z, yp1x, yp1y, yp1z,
     yp2x, yp2y, yp2z, ovx, ovy, ovz) = ext
    return ext, [zp1x, zp1y, zp1z], [zp2x, zp2y, zp2z], [yp1x, yp1y, yp1z],
           [yp2x, yp2y, yp2z], [ovx, ovy, ovz]
end

"""
    body_frame_damp_accel(vel, body_damp, zp1, zp2, yp1, yp2, ovel)

The body-frame damping acceleration `R·(coeff ⊙ (Rᵀ·(vel − wing_vel)))`, with the
wing frame `R` fitted from the ref points ([`wing_frame_columns`](@ref)) and the
wing velocity taken as the origin velocity `ovel`.
"""
function body_frame_damp_accel(vel, body_damp, zp1, zp2, yp1, yp2, ovel)
    xaxis, yaxis, zaxis = wing_frame_columns(zp1, zp2, yp1, yp2)
    rot = [xaxis[1] yaxis[1] zaxis[1];
           xaxis[2] yaxis[2] zaxis[2];
           xaxis[3] yaxis[3] zaxis[3]]
    return rot * (collect(body_damp) .* (rot' * (collect(vel) .- ovel)))
end

"""
    network_body_damped_point(s; name)

A DYNAMIC particle that additionally carries `body_frame_damping` — a damping of
its velocity *relative to its wing*, expressed in the fitted wing frame
(`point_eqs.jl`'s `point_damping_accel`). The wing frame and wing velocity are read
through NetworkDynamics external inputs from the wing's ref points
([`BODY_DAMP_EXTIN`]); the damping acceleration is subtracted from the shared
`point_acceleration`.
"""
function network_body_damped_point(s; name)
    vars, pos, vel, force, mass_in, _, pulley_len_out = vertex_io()
    ext, zp1, zp2, yp1, yp2, ovel = body_damp_inputs()
    append!(vars, ext)
    pars, named = particle_params()
    body_damp = SAM.make_array_param(:body_damp, zeros(3))
    accel = SAM.point_acceleration(s, collect(pos), collect(vel), collect(force),
        named.extra_mass + mass_in, named.drag_coeff, named.area,
        collect(named.world_damping), collect(named.wind_gnd))
    damp = body_frame_damp_accel(vel, body_damp, zp1, zp2, yp1, yp2, ovel)
    eqs = [
        D.(collect(pos)) .~ collect(vel);
        D.(collect(vel)) .~ accel .- damp;
        pulley_len_out ~ 0.0;
    ]
    return System(eqs, t, vars, [pars; body_damp]; name)
end

"""
    network_body_damped_pulley_point(s; name)

A dynamic pulley vertex (as [`network_pulley_point`](@ref)) whose particle motion
also carries `body_frame_damping`, read from the wing ref points through
[`BODY_DAMP_EXTIN`] like [`network_body_damped_point`](@ref). Used for pulley
points that belong to a wing (e.g. steering points).
"""
function network_body_damped_pulley_point(s; name)
    vars, pos, vel, force, mass_in, tension_in, pulley_len_out = vertex_io()
    extra = @variables pulley_len(t) pulley_vel(t)
    ext, zp1, zp2, yp1, yp2, ovel = body_damp_inputs()
    append!(vars, extra); append!(vars, ext)
    pars, named = particle_params()
    pulley_mass = SAM.make_param(:pulley_mass, 1.0)
    pulley_damp = SAM.make_param(:pulley_damp, 5.0)
    body_damp = SAM.make_array_param(:body_damp, zeros(3))
    accel = SAM.point_acceleration(s, collect(pos), collect(vel), collect(force),
        named.extra_mass + mass_in, named.drag_coeff, named.area,
        collect(named.world_damping), collect(named.wind_gnd))
    damp = body_frame_damp_accel(vel, body_damp, zp1, zp2, yp1, yp2, ovel)
    eqs = [
        D.(collect(pos)) .~ collect(vel);
        D.(collect(vel)) .~ accel .- damp;
        D(pulley_len) ~ pulley_vel;
        D(pulley_vel) ~ tension_in / pulley_mass - pulley_damp * pulley_vel;
        pulley_len_out ~ pulley_len;
    ]
    return System(eqs, t, vars, [pars; pulley_mass; pulley_damp; body_damp]; name)
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
    pos_w = SAM.make_array_param(:pos_w, zeros(3))
    set_value = SAM.make_param(:set_value, 0.0)
    brake = SAM.make_param(:brake, 0.0)
    speed_controlled = SAM.make_param(:speed_controlled, 0.0)
    pars = [pos_w, set_value, brake, speed_controlled]

    view = SAM.ParamView(SAM.ParamRegistry(s.sys_struct))
    motor = SAM.winch_component(winch.model, s.sys_struct, winch.idx;
                                name = :motor, params = view)
    SAM.validate_winch_component(motor, winch)
    winch_acc = ifelse(speed_controlled > 0.5, 0.0, motor.acc)

    eqs = [
        collect(pos) .~ collect(pos_w);
        collect(vel) .~ zeros(3);
        pulley_len_out ~ 0.0;
        winch_force ~ SAM.smooth_norm(tension_in);
        motor.vel ~ winch_vel;
        motor.len ~ sum(tether_lens) / n_tethers;
        motor.force ~ winch_force;
        motor.set_value ~ set_value;
        motor.brake ~ brake;
        D(winch_vel) ~ ifelse(brake > 0.5, 0.0, winch_acc);
    ]
    for tl in tether_lens
        push!(eqs, D(tl) ~ ifelse(brake > 0.5, 0.0, winch_vel))
    end
    return System(eqs, t, vars, pars; name, systems = [motor])
end

# ======================= edge Systems ======================= #

"""
    edge_io()

The uniform edge input/output variables: the src/dst vertex outputs read as inputs
(`src_pos`/`src_vel`/`src_pulley_len`, likewise `dst_…`) and the endpoint loads
written as outputs (`src_force`/`src_mass`/`src_tension`, likewise `dst_…`), all as
array-valued variables where they are 3-vectors. Returns `(vars, src_pos, src_vel,
src_pulley_len, dst_pos, dst_vel, dst_pulley_len, src_force, src_mass, src_tension,
dst_force, dst_mass, dst_tension)`.
"""
function edge_io()
    vars = @variables begin
        src_pos(t)[1:3], [input = true]
        src_vel(t)[1:3], [input = true]
        src_pulley_len(t), [input = true]
        dst_pos(t)[1:3], [input = true]
        dst_vel(t)[1:3], [input = true]
        dst_pulley_len(t), [input = true]
        src_force(t)[1:3], [output = true]
        src_mass(t), [output = true]
        src_tension(t), [output = true]
        dst_force(t)[1:3], [output = true]
        dst_mass(t), [output = true]
        dst_tension(t), [output = true]
    end
    return vars, src_pos, src_vel, src_pulley_len, dst_pos, dst_vel, dst_pulley_len,
           src_force, src_mass, src_tension, dst_force, dst_mass, dst_tension
end

"""
    spring_parameters()

The spring-damper + tether-drag parameters common to all edge types, plus the
ground wind (a `[1:3]` array parameter). Returns `(pars, spring_named, wind_gnd)`.
"""
function spring_parameters()
    unit_stiffness = SAM.make_param(:unit_stiffness, 0.0)
    unit_damping = SAM.make_param(:unit_damping, 0.0)
    compression_frac = SAM.make_param(:compression_frac, 0.1)
    diameter = SAM.make_param(:diameter, 0.0)
    density = SAM.make_param(:density, 0.0)
    cd_tether = SAM.make_param(:cd_tether, 1.0)
    wind_gnd = SAM.make_array_param(:wind_gnd, zeros(3))
    pars = [unit_stiffness, unit_damping, compression_frac, diameter, density,
            cd_tether, wind_gnd]
    spring = (; unit_stiffness, unit_damping, compression_frac, diameter, density,
              cd_tether)
    return pars, spring, wind_gnd
end

"""
    segment_loads(s, io, l0, spring, wind)

Compute the positive endpoint forces, half-mass and scalar spring tension from the
shared `segment_endpoint_loads`, reading the array-valued endpoint states out of
`io`. Returns `(fsrc, fdst, half, spring_scalar)`.
"""
function segment_loads(s, io, l0, spring, wind)
    (_, src_pos, src_vel, _, dst_pos, dst_vel, _) = io
    return SAM.segment_endpoint_loads(
        s, collect(src_pos), collect(src_vel), collect(dst_pos), collect(dst_vel),
        spring.unit_stiffness, spring.unit_damping, spring.compression_frac, l0,
        spring.diameter, spring.density, spring.cd_tether, collect(wind))
end

"""
    endpoint_load_eqs(io, fsrc, fdst, half, src_tension_val, dst_tension_val)

Bind the edge output variables (`src_force`/`src_mass`/`src_tension`, `dst_…`) from
the computed endpoint loads and the role-signed tensions to emit to each endpoint.
"""
function endpoint_load_eqs(io, fsrc, fdst, half, src_tension_val, dst_tension_val)
    (_, _, _, _, _, _, _, src_force, src_mass, src_tension,
     dst_force, dst_mass, dst_tension) = io
    return [
        collect(src_force) .~ fsrc;
        src_mass ~ half; src_tension ~ src_tension_val;
        collect(dst_force) .~ fdst;
        dst_mass ~ half; dst_tension ~ dst_tension_val;
    ]
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
    l0 = SAM.make_param(:l0, 1.0)
    fsrc, fdst, half, _ = segment_loads(s, io, l0, spring, wind)
    eqs = endpoint_load_eqs(io, fsrc, fdst, half, 0.0, 0.0)
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
    pulley_sum_len = SAM.make_param(:pulley_sum_len, 1.0)
    pulley_side = SAM.make_param(:pulley_side, 1.0)
    pulley_at_src = SAM.make_param(:pulley_at_src, 0.0)
    endpoint_len = ifelse(pulley_at_src > 0.5, src_pulley_len, dst_pulley_len)
    l0 = ifelse(pulley_side > 0.0, endpoint_len, pulley_sum_len - endpoint_len)
    fsrc, fdst, half, spring_scalar = segment_loads(s, io, l0, spring, wind)
    src_tension_val = ifelse(pulley_at_src > 0.5, pulley_side * spring_scalar, 0.0)
    dst_tension_val = ifelse(pulley_at_src > 0.5, 0.0, pulley_side * spring_scalar)
    eqs = endpoint_load_eqs(io, fsrc, fdst, half, src_tension_val, dst_tension_val)
    return System(eqs, t, vars, [spars; pulley_sum_len; pulley_side; pulley_at_src];
                  name)
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
    tether_len_ext = only(@variables tether_len_ext(t) [input = true])
    push!(vars, tether_len_ext)
    spars, spring, wind = spring_parameters()
    n_segs = SAM.make_param(:n_segs, 1.0)
    tension_sign_src = SAM.make_param(:tension_sign_src, 0.0)
    tension_sign_dst = SAM.make_param(:tension_sign_dst, 0.0)
    l0 = tether_len_ext / n_segs
    fsrc, fdst, half, spring_scalar = segment_loads(s, io, l0, spring, wind)
    eqs = endpoint_load_eqs(io, fsrc, fdst, half,
        tension_sign_src * spring_scalar, tension_sign_dst * spring_scalar)
    return System(eqs, t, vars, [spars; n_segs; tension_sign_src; tension_sign_dst];
                  name)
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
    point_body_damp(ss, point)

The point's `body_frame_damping` coefficients if it is a DYNAMIC wing node on a
`KINEMATIC` wing with a nonzero coefficient, else `nothing`. Gated on
`is_wing_node` to match the monolith, which only damps twist-surface members;
steering/pulley points carry a `body_frame_damping` field but are not wing nodes.
"""
function point_body_damp(ss, point)
    (point.type == SAM.DYNAMIC && point.is_wing_node &&
        0 < point.wing_idx <= length(ss.bodies)) || return nothing
    ss.bodies[point.wing_idx].type == SAM.KINEMATIC || return nothing
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

    bd_vmodel_of = build_body_damped_vmodels(sam, ss, pulley_of_point)

    vmodels = Vector{VertexModel}(undef, n)
    for i in 1:n
        if haskey(bd_vmodel_of, i)
            vmodels[i] = bd_vmodel_of[i]
        elseif haskey(winch_of_point, i)
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
    set_states!(ss, state, winch_of_point, pulley_of_point)
    builder = ParamBuilder()
    record_vertex_params!(builder, ss, winch_of_point, pulley_of_point)
    record_edge_params!(builder, ss, edgelist, seg_of, role_of_seg)
    param_sync = build_network_param_sync(nw, builder)

    meta = (; param_sync, winch_of_point, pulley_of_point,
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
    build_body_damped_vmodels(sam, ss, pulley_of_point)

One `VertexModel` per DYNAMIC point carrying `body_frame_damping` on a KINEMATIC
wing (selected by [`point_body_damp`](@ref)), reading its wing's ref-point frame
through `extin`. The dynamic and pulley kernels are each compiled once and rebound
per point with `VertexModel(base; extin=…)`, so only the wing's ref-point indices
differ. Returns a `point_idx => VertexModel` map (empty when no point is damped).
"""
function build_body_damped_vmodels(sam, ss, pulley_of_point)
    vmodel_of = Dict{Int, Any}()
    dyn_base = nothing
    pulley_base = nothing
    for (i, point) in enumerate(ss.points)
        point_body_damp(ss, point) === nothing && continue
        extin = body_damp_extin(ss, ss.bodies[point.wing_idx])
        if haskey(pulley_of_point, i)
            if pulley_base === nothing
                pulley_base = VertexModel(
                    network_body_damped_pulley_point(sam; name = :bdpul),
                    VERTEX_INPUTS, VERTEX_OUTPUTS; extin, mtkcompile = true,
                    name = :bdpul)
                vmodel_of[i] = pulley_base
            else
                vmodel_of[i] = VertexModel(pulley_base; extin = last.(extin))
            end
        elseif dyn_base === nothing
            dyn_base = VertexModel(network_body_damped_point(sam; name = :bdyn),
                VERTEX_INPUTS, VERTEX_OUTPUTS; extin, mtkcompile = true,
                name = :bdyn)
            vmodel_of[i] = dyn_base
        else
            vmodel_of[i] = VertexModel(dyn_base; extin = last.(extin))
        end
    end
    return vmodel_of
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
    record_vertex_params!(builder, ss, winch_of_point, pulley_of_point)

Record every vertex's per-instance parameters (index + live reader) into `builder`,
dispatched on the vertex type (winch, pulley, static, particle).
"""
function record_vertex_params!(builder, ss, winch_of_point, pulley_of_point)
    for (i, point) in enumerate(ss.points)
        if haskey(winch_of_point, i)
            record_winch_params!(builder, i, winch_of_point[i])
        elseif haskey(pulley_of_point, i)
            record_particle_params!(builder, i)
            add_param!(builder, VIndex(i, :pulley_mass),
                       PulleyLineMassReader(pulley_of_point[i]))
            point_body_damp(ss, point) === nothing || record_body_damp_params!(builder, i)
        elseif point.type == SAM.STATIC
            record_pos_w_params!(builder, i)
        else
            record_particle_params!(builder, i)
            point_body_damp(ss, point) === nothing || record_body_damp_params!(builder, i)
        end
    end
    return nothing
end

"""
    record_body_damp_params!(builder, i)

Record the `body_damp_k` coefficients for a body-damped vertex `i`, read live from
the point's `body_frame_damping`.
"""
function record_body_damp_params!(builder, i)
    for k in 1:3
        add_param!(builder, VIndex(i, Symbol(:body_damp_, k)),
                   SAM.PathReader((:points, i, :body_frame_damping, k)))
    end
    return nothing
end

"""
    record_particle_params!(builder, i)

Record the shared DYNAMIC-particle parameters for vertex `i`: mass, drag, area, the
world-frame damping (`world_damping_k`) and the ground wind (`wind_gnd_k`).
"""
function record_particle_params!(builder, i)
    add_param!(builder, VIndex(i, :extra_mass), SAM.PathReader((:points, i, :extra_mass)))
    add_param!(builder, VIndex(i, :drag_coeff), SAM.PathReader((:points, i, :drag_coeff)))
    add_param!(builder, VIndex(i, :area), SAM.PathReader((:points, i, :area)))
    for k in 1:3
        add_param!(builder, VIndex(i, Symbol(:world_damping_, k)),
                   SAM.PathReader((:points, i, :world_frame_damping, k)))
        add_param!(builder, VIndex(i, Symbol(:wind_gnd_, k)), GroundWindReader(k))
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

Record every edge's spring/drag parameters and role-specific data (frozen `l0`,
pulley split, or winched-tether `n_segs`/tension signs) into `builder`.
"""
function record_edge_params!(builder, ss, edgelist, seg_of, role_of_seg)
    points = ss.points
    for (j, e) in enumerate(edgelist)
        seg = seg_of[minmax(src(e), dst(e))]
        role = role_of_seg[seg.idx]
        wing_structural = points[seg.point_idxs[1]].is_wing_node &&
                          points[seg.point_idxs[2]].is_wing_node
        segment_stiffness(seg)  # validate linear stiffness up front
        add_param!(builder, EIndex(j, :unit_stiffness),
                   SAM.PathReader((:segments, seg.idx, :unit_stiffness)))
        add_param!(builder, EIndex(j, :unit_damping),
                   SAM.PathReader((:segments, seg.idx, :unit_damping)))
        add_param!(builder, EIndex(j, :compression_frac),
                   SAM.PathReader((:segments, seg.idx, :compression_frac)))
        add_param!(builder, EIndex(j, :diameter),
                   SAM.PathReader((:segments, seg.idx, :diameter)))
        add_param!(builder, EIndex(j, :density),
                   SAM.PathReader((:segments, seg.idx, :density)))
        add_param!(builder, EIndex(j, :cd_tether), TetherDragReader(wing_structural))
        for k in 1:3
            add_param!(builder, EIndex(j, Symbol(:wind_gnd_, k)), GroundWindReader(k))
        end
        if role.kind == :pulley
            add_param!(builder, EIndex(j, :pulley_sum_len),
                       SAM.PathReader((:pulleys, role.pulley_idx, :sum_len)))
            add_param!(builder, EIndex(j, :pulley_side), ConstReader(role.pulley_side))
            add_param!(builder, EIndex(j, :pulley_at_src),
                       ConstReader(role.pulley_at_src ? 1.0 : 0.0))
        elseif role.kind == :tether
            add_param!(builder, EIndex(j, :n_segs), ConstReader(Float64(role.n_segs)))
            add_param!(builder, EIndex(j, :tension_sign_src),
                       ConstReader(role.tension_sign_src))
            add_param!(builder, EIndex(j, :tension_sign_dst),
                       ConstReader(role.tension_sign_dst))
        else
            add_param!(builder, EIndex(j, :l0), SAM.PathReader((:segments, seg.idx, :l0)))
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
        for k in 1:3
            point.pos_w[k] = s.v[i, Symbol(:pos_, k)]
            point.vel_w[k] = s.v[i, Symbol(:vel_, k)]
        end
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
