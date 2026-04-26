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
using KiteUtils: init!, next_step!, update_sys_state!
using SymbolicAWEModels
import SymbolicAWEModels: Point  # resolve ambiguity with GLMakie

set_data_path(joinpath(dirname(@__DIR__), "data"))
set = Settings("base/system.yaml")
set.v_wind = 10.0
set.abs_tol = 1e-3
set.rel_tol = 1e-3

points = [
    Point(:anchor1, [0, 0, 2], STATIC),
    Point(:anchor2, [2, 0, 2], STATIC),
    Point(:pulley_pt, [0.1, 0, 1], DYNAMIC;
        extra_mass=0.1),
    Point(:mass, [0.1, 0, 0], DYNAMIC;
        extra_mass=0.1),
]
segments = [
    Segment(:seg1, :pulley_pt, :anchor1,
            500.0, 50.0, 0.004),
    Segment(:seg2, :pulley_pt, :anchor2,
            500.0, 50.0, 0.004),
    Segment(:seg3, :pulley_pt, :mass,
            500.0, 50.0, 0.004),
]
pulleys = [Pulley(:pulley, :seg1, :seg2, DYNAMIC)]
sys = SystemStructure("pulley", set;
                      points, segments, pulleys)
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
