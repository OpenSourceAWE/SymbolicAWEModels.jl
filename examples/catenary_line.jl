# Copyright (c) 2025 Bart van de Lint, Jelle Poland
# SPDX-License-Identifier: LGPL-3.0-only

"""
Catenary line: tether fixed at both ends under gravity, relaxing
into its equilibrium shape.
"""

using Pkg
if Base.active_project() != joinpath(@__DIR__, "Project.toml")
    Pkg.activate(@__DIR__)
end

using GLMakie
using MakieControlPlots

using KiteUtils: init!, next_step!, update_sys_state!
using SymbolicAWEModels
import SymbolicAWEModels: Point  # resolve ambiguity with GLMakie

set_data_path(joinpath(dirname(@__DIR__), "data"))
set = Settings("base/system.yaml")
set.v_wind = 0.0

horizontal_span = 8.0
n_segments = 10
total_length = 10.0

# Segment properties (soft rope, 4mm diameter)
seg_stiffness = 500.0  # EA [N]
seg_damping   = 50.0   # [N·s]
seg_diameter  = 0.004  # [m]

points = [Point(; name=:left, pos_ENU=[0, 0, 5], fixed=true)]
for i in 1:n_segments - 1
    x = i * horizontal_span / n_segments
    push!(points, Point(; name=Symbol("p$i"), pos_ENU=[x, 0, 5],
                        model=PointModel(; world_frame_damping=1.0)))
end
push!(points, Point(; name=:right, pos_ENU=[horizontal_span, 0, 5], fixed=true))

rest_length_seg = total_length / n_segments
materials = KiteGeometry.MaterialTable(KiteGeometry.Material[], Dict{String, Any}())
segments = KiteGeometry.Segment[]
for i in 1:n_segments
    push!(segments, KiteGeometry.line_segment(
        Symbol("seg$i"), (points[i].name, points[i + 1].name), set, materials;
        rest_length=rest_length_seg, diameter=seg_diameter, unit_stiffness=seg_stiffness,
        model=SegmentSpring(; unit_damping=seg_damping, compression_frac=0.01)))
end

sys = system_definition("catenary", set; points, segments, materials)
sam = SymbolicAWEModel(set, sys)
init!(sam; remake=false)

n_steps = 200
logger = Logger(sam, n_steps)
sys_state = SysState(sam)

for i in 1:n_steps
    next_step!(sam)
    update_sys_state!(sys_state, sam)
    sys_state.time = i / set.sample_freq
    log!(logger, sys_state)
end

save_log(logger, "catenary_line")
syslog = load_log("catenary_line")
scene = replay(syslog, sam)
display(scene)
