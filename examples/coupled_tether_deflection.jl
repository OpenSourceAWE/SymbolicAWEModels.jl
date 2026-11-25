# Copyright (c) 2025 Bart van de Lint, Jelle Poland
# SPDX-License-Identifier: LGPL-3.0-only

"""
Vertical tether under pure horizontal wind drag, with gravity set
to zero to isolate aerodynamic effects.
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
set.abs_tol = 1e-8
set.rel_tol = 1e-8
set.upwind_dir = 0.0
set.g_earth = 0.0
set.v_wind = 5.0

vertical_span = 8.0
n_segments = 10
total_length = 8.2

# Segment properties (Dyneema-like, 4mm diameter)
seg_stiffness = 500_000.0  # EA [N]
seg_damping   = 1000.0     # [N·s]
seg_diameter  = 0.004      # [m]

# Points laid out horizontally; transform rotates to vertical
points = [
    Point(; name=:top, pos_ENU=[vertical_span, 0, 0], fixed=true),
]
for i in 1:n_segments - 1
    x = vertical_span - i * vertical_span / n_segments
    push!(points, Point(; name=Symbol("p$i"), pos_ENU=[x, 0, 0],
                        model=PointModel(; world_frame_damping=1.0)))
end
push!(points, Point(; name=:bottom, pos_ENU=[0, 0, 0], fixed=true))

rest_length_seg = total_length / n_segments
materials = KiteGeometry.MaterialTable(KiteGeometry.Material[], Dict{String, Any}())
segments = KiteGeometry.Segment[]
for i in 1:n_segments
    push!(segments, KiteGeometry.line_segment(
        Symbol("seg$i"), (points[i].name, points[i + 1].name), set, materials;
        rest_length=rest_length_seg, diameter=seg_diameter, unit_stiffness=seg_stiffness,
        model=SegmentSpring(; unit_damping=seg_damping, compression_frac=0.01)))
end

transforms = [
    Transform(; name=:tf, elevation=-deg2rad(90), azimuth=0, heading=0,
              base_pos=[0, 0, vertical_span], base_point=1,
              rot_point=[length(points) => 1.0], points=collect(eachindex(points))),
]

sys = system_definition("wind_drag", set; points, segments, materials)
place!(sys, transforms)
sam = SymbolicAWEModel(set, sys; transforms)
init!(sam; remake=false)

n_steps = 800
logger = Logger(sam, n_steps)
sys_state = SysState(sam)

for i in 1:n_steps
    next_step!(sam)
    update_sys_state!(sys_state, sam)
    sys_state.time = i / set.sample_freq
    log!(logger, sys_state)
end

save_log(logger, "tether_deflection")
syslog = load_log("tether_deflection")
scene = replay(syslog, sam)
display(scene)
@info "Wind deflection simulation complete"
