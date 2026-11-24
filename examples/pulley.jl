# SPDX-FileCopyrightText: 2025 Bart van de Lint
#
# SPDX-License-Identifier: LGPL-3.0-only

"""
Pulley demo: a dynamic point connected through a pulley constraint
to two anchors, with a hanging mass below.
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
set.v_wind = 10.0
set.abs_tol = 1e-3
set.rel_tol = 1e-3

materials = KiteGeometry.MaterialTable(KiteGeometry.Material[], Dict{String, Any}())
points = [
    Point(; name=:anchor1, pos_ENU=[0, 0, 2], fixed=true),
    Point(; name=:anchor2, pos_ENU=[2, 0, 2], fixed=true),
    Point(; name=:pulley_pt, pos_ENU=[0.1, 0, 1], extra_mass=0.1),
    Point(; name=:mass, pos_ENU=[0.1, 0, 0], extra_mass=0.1),
]
# Rope segments: EA = 500 N, 50 N·s damping, 4 mm diameter
rope(name, ends) = KiteGeometry.line_segment(name, ends, set, materials;
    diameter=0.004, unit_stiffness=500.0, model=SegmentSpring(; unit_damping=50.0))
segments = [
    rope(:seg1, (:pulley_pt, :anchor1)),
    rope(:seg2, (:pulley_pt, :anchor2)),
    rope(:seg3, (:pulley_pt, :mass)),
]
pulleys = [Pulley(; name=:pulley, segments=(:seg1, :seg2),
                  model=PulleyModel(; type="DYNAMIC"))]
sys = system_definition("pulley", set;
                      points, segments, materials, pulleys)
sam = SymbolicAWEModel(set, sys)
init!(sam; remake=false)

n_steps = 100
logger = Logger(sam, n_steps)
sys_state = SysState(sam)

for i in 1:n_steps
    next_step!(sam)
    update_sys_state!(sys_state, sam)
    sys_state.time = i / set.sample_freq
    log!(logger, sys_state)
end

save_log(logger, "pulley")
syslog = load_log("pulley")
scene = replay(syslog, sam.sys_struct)
display(scene)
