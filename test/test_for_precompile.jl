# Copyright (c) 2025 Bart van de Lint, Jelle Poland, Uwe Fechner
# SPDX-License-Identifier: LGPL-3.0-only

using GLMakie
using KiteUtils: init!, next_step!, update_sys_state!
using SymbolicAWEModels
import SymbolicAWEModels: Point  # resolve ambiguity with GLMakie

set_data_path(joinpath(dirname(@__DIR__), "data"))
set = Settings("base/system.yaml")
set.v_wind = 0

points = [
    Point(:anchor, [2, 0, 5], STATIC),
    Point(:mass, [2, 0, 2], DYNAMIC; extra_mass=1.0),
]
segments = [
    Segment(:spring, :anchor, :mass,
            500.0, 50.0, 0.005; l0=4.0),
]
transforms = [
    Transform(:tf, -deg2rad(90), 0, 0;
              base_pos=[2, 0, 5], base_point=:anchor,
              rot_point=:mass),
]

sys = SystemStructure("hanging_mass", set;
                      points, segments, transforms)
sam = SymbolicAWEModel(set, sys)
init!(sam; remake=false)

n_steps = 30
logger = Logger(sam, n_steps)
sys_state = SysState(sam)

for i in 1:n_steps
    next_step!(sam)
    update_sys_state!(sys_state, sam)
    sys_state.time = i / set.sample_freq
    log!(logger, sys_state)
end

prev_data_path = get_data_path()
tmpdir = mktempdir()
syslog = try
    set_data_path(tmpdir)
    save_log(logger, "_hanging_mass")
    load_log("_hanging_mass")
finally
    rm(tmpdir; recursive=true)
    set_data_path(prev_data_path)
end

scene = replay(syslog, sam.sys_struct)
if isinteractive()
    display(scene)
end
nothing
