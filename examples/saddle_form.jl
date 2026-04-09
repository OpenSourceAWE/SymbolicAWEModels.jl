# SPDX-FileCopyrightText: 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

"""
Saddle form relaxation: diamond mesh with fixed boundary nodes
and saddle z-profile, relaxed to equilibrium via dynamic simulation.
"""

using Pkg
if Base.active_project() != joinpath(@__DIR__, "Project.toml")
    Pkg.activate(@__DIR__)
end

using GLMakie
using KiteUtils: azimuth_east, calc_elevation, init!, next_step!, update_sys_state!
using SymbolicAWEModels
import SymbolicAWEModels: Point  # resolve ambiguity with GLMakie
using LinearAlgebra
using YAML

function load_saddle(yaml_file;
                     unit_stiffness, unit_damping, diameter,
                     rest_length, world_frame_damping,
                     compression_frac)
    data = YAML.load_file(yaml_file)
    nodes = data["nodes"]
    one_based = get(data, "one_based", true)

    points = Point[]
    for (i, n) in enumerate(nodes)
        pos = [Float64(n["x"]), Float64(n["y"]),
               Float64(n["z"])]
        is_fixed = get(n, "fixed", false)
        mass = Float64(get(n, "mass", 1))
        type = is_fixed ? STATIC : DYNAMIC
        kw = is_fixed ? (; extra_mass=mass) :
            (; extra_mass=mass,
               world_frame_damping=world_frame_damping)
        push!(points, Point(i, pos, type; kw...))
    end

    segments = Segment[]
    for (sid, c) in enumerate(data["connections"])
        pi_idx = Int(c["i"]); pj_idx = Int(c["j"])
        if !one_based
            pi_idx += 1; pj_idx += 1
        end
        push!(segments, Segment(sid, pi_idx, pj_idx,
            unit_stiffness, unit_damping, diameter;
            l0=rest_length, compression_frac))
    end

    return points, segments
end

function run_saddle(; yaml_file)
    project_dir = dirname(@__DIR__)
    set_data_path(joinpath(project_dir, "data",
                           "saddle_form"))
    set = Settings("system.yaml")
    set.v_wind = 0.0
    set.g_earth = 0.0
    set.l_tether = 10
    n_steps = 100

    points, segments = load_saddle(yaml_file;
        unit_stiffness=24.0, unit_damping=1.0,
        diameter=0.002, rest_length=0.01,
        world_frame_damping=1.0,
        compression_frac=1.0)

    sys = SystemStructure("saddle_form", set;
                          points, segments)
    sam = SymbolicAWEModel(set, sys)
    init!(sam; remake=false)

    logger = Logger(sam, n_steps)
    sys_state = SysState(sam)

    for i in 1:n_steps
        next_step!(sam)
        update_sys_state!(sys_state, sam)
        sys_state.time = i / set.sample_freq
        log!(logger, sys_state)
    end

    save_log(logger, "saddle_form")
    syslog = load_log("saddle_form")
    scene = replay(syslog, sam.sys_struct)
    display(scene)

    fixed = Set(i for (i, p) in enumerate(points)
                if p.type == STATIC)
    disp = [norm(sam.sys_struct.points[i].pos_w .-
                 points[i].pos_cad)
            for i in eachindex(points) if i ∉ fixed]
    avg = sum(disp) / length(disp)
    @info "Relaxation" max_disp=round(maximum(disp);
        digits=4) avg_disp=round(avg; digits=4)
    return sam
end

project_dir = dirname(@__DIR__)
sam = run_saddle(yaml_file=joinpath(project_dir, "data",
      "saddle_form", "saddle_gridsize4.yaml"))
@info "Type 'sam' to inspect the final model."
