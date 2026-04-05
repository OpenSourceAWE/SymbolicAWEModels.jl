# Copyright (c) 2025 Bart van de Lint, Jelle Poland
# SPDX-License-Identifier: MPL-2.0

"""
Vertical tether under pure horizontal wind drag, with gravity set
to zero to isolate aerodynamic effects.
"""

using Pkg
Pkg.activate(@__DIR__)

using GLMakie
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
    Point(:top, [vertical_span, 0, 0], STATIC),
]
for i in 1:n_segments - 1
    x = vertical_span - i * vertical_span / n_segments
    push!(points, Point(Symbol("p$i"), [x, 0, 0], DYNAMIC;
                        world_frame_damping=1.0))
end
push!(points, Point(:bottom, [0, 0, 0], STATIC))

l0_seg = total_length / n_segments
segments = Segment[]
for i in 1:n_segments
    push!(segments, Segment(
        Symbol("seg$i"), points[i].name,
        points[i + 1].name,
        seg_stiffness, seg_damping, seg_diameter;
        l0=l0_seg, compression_frac=0.01))
end

transforms = [
    Transform(:tf, -deg2rad(90), 0, 0;
              base_pos=[0, 0, vertical_span],
              base_point=:top, rot_point=:bottom),
]

sys = SystemStructure("wind_drag", set;
                      points, segments, transforms)
sam = SymbolicAWEModel(set, sys)
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
scene = replay(syslog, sam.sys_struct)
display(scene)
@info "Wind deflection simulation complete"
