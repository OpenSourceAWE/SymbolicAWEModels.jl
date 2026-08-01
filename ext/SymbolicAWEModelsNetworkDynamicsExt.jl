# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

"""
    SymbolicAWEModelsNetworkDynamicsExt

Network backend: assemble a `SymbolicAWEModel` as a NetworkDynamics `Network`,
compiling one kernel per component *type* (flat in node count). Loaded when both
`NetworkDynamics` and `Graphs` are available. Implements
`build_prob!(::NetworkBackend, sam)`.

The per-type `System`s here are the network's assembly of the same physics the
monolith connector components use: they call the shared `point_acceleration` /
`segment_endpoint_loads` helpers (D2). The only difference is the aggregation
sign — ND sums edge outputs *into* the vertex input, so edges emit the **positive**
endpoint loads (the monolith's `Flow` convention negates them). Variables are
scalarized to match the ND MTK-integration convention (bare `[input]`/`[output]`
scalar names, not `@connector` variables).
"""
module SymbolicAWEModelsNetworkDynamicsExt

using SymbolicAWEModels
using NetworkDynamics
using Graphs
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

const SAM = SymbolicAWEModels

"""
    network_dynamic_point(s; name)

Particle vertex `System` for the network backend: states `pos1..3`/`vel1..3`,
inputs `force_in1..3`/`mass_in` (aggregated edge loads), outputs `pos`/`vel`.
`D(vel)` comes from the shared `point_acceleration`.
"""
function network_dynamic_point(s; name)
    vars = @variables begin
        pos1(t); pos2(t); pos3(t); vel1(t); vel2(t); vel3(t)
        force_in1(t), [input = true]
        force_in2(t), [input = true]
        force_in3(t), [input = true]
        mass_in(t), [input = true]
    end
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
    pos = [pos1, pos2, pos3]
    vel = [vel1, vel2, vel3]
    force = [force_in1, force_in2, force_in3]
    mass = extra_mass + mass_in
    accel = SAM.point_acceleration(s, pos, vel, force, mass, drag_coeff, area,
        [world_damping1, world_damping2, world_damping3],
        [wind_gnd1, wind_gnd2, wind_gnd3])
    eqs = [
        D(pos1) ~ vel1, D(pos2) ~ vel2, D(pos3) ~ vel3,
        D(vel1) ~ accel[1], D(vel2) ~ accel[2], D(vel3) ~ accel[3],
    ]
    return System(eqs, t, vars, pars; name)
end

"""
    network_static_point(s; name)

Ground-anchored vertex `System`: outputs `pos = pos_w`, `vel = 0`. Declares the
same `force_in`/`mass_in` inputs as the dynamic vertex (edges deliver to it) but
ignores them.
"""
function network_static_point(s; name)
    vars = @variables begin
        pos1(t); pos2(t); pos3(t); vel1(t); vel2(t); vel3(t)
        force_in1(t), [input = true]
        force_in2(t), [input = true]
        force_in3(t), [input = true]
        mass_in(t), [input = true]
    end
    pars = @parameters pos_w1 = 0.0 pos_w2 = 0.0 pos_w3 = 0.0
    eqs = [
        pos1 ~ pos_w1, pos2 ~ pos_w2, pos3 ~ pos_w3,
        vel1 ~ 0.0, vel2 ~ 0.0, vel3 ~ 0.0,
    ]
    return System(eqs, t, vars, pars; name)
end

"""
    network_segment(s; name)

Spring-damper edge `System`: inputs are the src/dst vertices' `pos`/`vel`, outputs
are the `+force`/`+half_mass` delivered to each endpoint from the shared
`segment_endpoint_loads` (positive sign — ND aggregates into the vertex input).
"""
function network_segment(s; name)
    vars = @variables begin
        src_pos1(t), [input = true]
        src_pos2(t), [input = true]
        src_pos3(t), [input = true]
        src_vel1(t), [input = true]
        src_vel2(t), [input = true]
        src_vel3(t), [input = true]
        dst_pos1(t), [input = true]
        dst_pos2(t), [input = true]
        dst_pos3(t), [input = true]
        dst_vel1(t), [input = true]
        dst_vel2(t), [input = true]
        dst_vel3(t), [input = true]
        src_force1(t), [output = true]
        src_force2(t), [output = true]
        src_force3(t), [output = true]
        src_mass(t), [output = true]
        dst_force1(t), [output = true]
        dst_force2(t), [output = true]
        dst_force3(t), [output = true]
        dst_mass(t), [output = true]
    end
    pars = @parameters begin
        unit_stiffness = 0.0
        unit_damping = 0.0
        compression_frac = 0.1
        l0 = 1.0
        diameter = 0.0
        density = 0.0
        cd_tether = 1.0
        wind_gnd1 = 0.0
        wind_gnd2 = 0.0
        wind_gnd3 = 0.0
    end
    src_pos = [src_pos1, src_pos2, src_pos3]
    src_vel = [src_vel1, src_vel2, src_vel3]
    dst_pos = [dst_pos1, dst_pos2, dst_pos3]
    dst_vel = [dst_vel1, dst_vel2, dst_vel3]
    fsrc, fdst, half_mass = SAM.segment_endpoint_loads(
        s, src_pos, src_vel, dst_pos, dst_vel, unit_stiffness, unit_damping,
        compression_frac, l0, diameter, density, cd_tether,
        [wind_gnd1, wind_gnd2, wind_gnd3])
    eqs = [
        src_force1 ~ fsrc[1], src_force2 ~ fsrc[2], src_force3 ~ fsrc[3],
        src_mass ~ half_mass,
        dst_force1 ~ fdst[1], dst_force2 ~ fdst[2], dst_force3 ~ fdst[3],
        dst_mass ~ half_mass,
    ]
    return System(eqs, t, vars, pars; name)
end

"""
    ground_wind(sam)

Ground-level wind vector `[vx, vy, vz]` from settings, matching the monolith's
`wind_vec_gnd`: magnitude `set.v_wind` along the `set.upwind_dir` heading.
"""
function ground_wind(sam)
    set = sam.set
    dir = deg2rad(set.upwind_dir)
    return [set.v_wind * cos(dir), set.v_wind * sin(dir), 0.0]
end

"""
    build_network(sam)

Translate `sam.sys_struct` (points + segments only) into a NetworkDynamics
`Network` with per-instance parameters and initial state. Returns
`(nw, u0, p0)` where `u0`/`p0` are the flat state/parameter vectors. STATIC points
become static vertices, DYNAMIC points dynamic vertices, each segment an edge.
"""
function build_network(sam)
    ss = sam.sys_struct
    points = ss.points
    segments = ss.segments
    n = length(points)
    graph = SimpleGraph(n)
    seg_of = Dict{Tuple{Int, Int}, Any}()
    for seg in segments
        a, b = seg.point_idxs
        add_edge!(graph, a, b)
        seg_of[minmax(a, b)] = seg
    end

    dyn = VertexModel(network_dynamic_point(sam; name = :dyn),
        [:force_in1, :force_in2, :force_in3, :mass_in],
        [:pos1, :pos2, :pos3, :vel1, :vel2, :vel3]; mtkcompile = true, name = :dyn)
    stat = VertexModel(network_static_point(sam; name = :stat),
        [:force_in1, :force_in2, :force_in3, :mass_in],
        [:pos1, :pos2, :pos3, :vel1, :vel2, :vel3]; mtkcompile = true, name = :stat)
    edge = EdgeModel(network_segment(sam; name = :seg),
        [:src_pos1, :src_pos2, :src_pos3, :src_vel1, :src_vel2, :src_vel3],
        [:dst_pos1, :dst_pos2, :dst_pos3, :dst_vel1, :dst_vel2, :dst_vel3],
        [:src_force1, :src_force2, :src_force3, :src_mass],
        [:dst_force1, :dst_force2, :dst_force3, :dst_mass];
        mtkcompile = true, name = :seg)

    vmodels = [points[i].type == SAM.STATIC ? stat : dyn for i in 1:n]
    edgelist = collect(edges(graph))
    emodels = [edge for _ in edgelist]
    nw = Network(graph, vmodels, emodels)

    wind = ground_wind(sam)
    param = NWParameter(nw)
    state = NWState(nw)
    for (i, point) in enumerate(points)
        if point.type == SAM.STATIC
            param.v[i, :pos_w1] = point.pos_w[1]
            param.v[i, :pos_w2] = point.pos_w[2]
            param.v[i, :pos_w3] = point.pos_w[3]
        else
            param.v[i, :extra_mass] = point.extra_mass
            param.v[i, :drag_coeff] = point.drag_coeff
            param.v[i, :area] = point.area
            damping = something(point.world_frame_damping, zeros(3))
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
        end
    end
    for (j, e) in enumerate(edgelist)
        seg = seg_of[minmax(src(e), dst(e))]
        stiff = seg.unit_stiffness isa Real ? Float64(seg.unit_stiffness) : 0.0
        param.e[j, :unit_stiffness] = stiff
        param.e[j, :unit_damping] = seg.unit_damping
        param.e[j, :compression_frac] = seg.compression_frac
        param.e[j, :l0] = seg.l0
        param.e[j, :diameter] = seg.diameter
        param.e[j, :density] = seg.density
        param.e[j, :cd_tether] = sam.set.cd_tether
        param.e[j, :wind_gnd1] = wind[1]
        param.e[j, :wind_gnd2] = wind[2]
        param.e[j, :wind_gnd3] = wind[3]
    end
    return nw, uflat(state), pflat(param)
end

function SAM.build_prob!(::SAM.NetworkBackend, sam; prn = true)
    nw, u0, p0 = build_network(sam)
    dt = SAM.SimFloat(1 / sam.set.sample_freq)
    prob = ODEProblem(nw, u0, (0.0, dt), p0)
    sam.prob = SAM.ProbWithAttributes(; prob)
    return true
end

end # module
