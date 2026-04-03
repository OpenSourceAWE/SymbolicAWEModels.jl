# Copyright (c) 2025 2-Plate Kite Coupled Simulation
# SPDX-License-Identifier: MPL-2.0

"""
2-Plate Kite: coupled aerodynamic-structural simulation with
ramped steering inputs and interactive replay.
"""

using GLMakie
using SymbolicAWEModels, VortexStepMethod, KiteUtils
using SymbolicAWEModels: init!, next_step!, update_sys_state!

MODEL_NAME = "2plate_kite"
SIM_TIME = 2.0
N_STEPS = 600
RAMP_TIME = 2.0
STEERING_MAGNITUDE = 0.1
DEPOWER_MAGNITUDE = 0.0

pkg_root = dirname(@__DIR__)
set_data_path(joinpath(pkg_root, "data", MODEL_NAME))

struc_yaml = joinpath(get_data_path(),
                      "refine_struc_geometry.yaml")
aero_yaml = joinpath(get_data_path(), "aero_geometry.yaml")
update_aero_yaml_from_struc_yaml!(struc_yaml, aero_yaml)

set = Settings("system.yaml")
set.g_earth = 0.0
vsm_set = VortexStepMethod.VSMSettings(
    joinpath(get_data_path(), "vsm_settings.yaml");
    data_prefix=false)

sys = load_sys_struct_from_yaml(struc_yaml;
    system_name=MODEL_NAME, set, vsm_set)
sys.winches[:main_winch].brake = true
sam = SymbolicAWEModel(set, sys)

l0_left = sam.sys_struct.segments[:kcu_steering_left].l0
l0_right = sam.sys_struct.segments[:kcu_steering_right].l0

init!(sam; remake=false, lin_vsm=false)

dt = SIM_TIME / N_STEPS
logger = Logger(sam, N_STEPS + 1)
sys_state = SysState(sam)
sys_state.time = 0.0
log!(logger, sys_state)

for step in 1:N_STEPS
    t = step * dt
    ramp = clamp(t / RAMP_TIME, 0.0, 1.0)
    steer = STEERING_MAGNITUDE * ramp
    depower = DEPOWER_MAGNITUDE * ramp
    sam.sys_struct.segments[:kcu_steering_left].l0 =
        l0_left - steer + depower
    sam.sys_struct.segments[:kcu_steering_right].l0 =
        l0_right + steer + depower

    next_step!(sam; dt, vsm_interval=1)

    update_sys_state!(sys_state, sam)
    sys_state.time = t
    log!(logger, sys_state)

    if step % max(1, div(N_STEPS, 10)) == 0
        @info "Step $step/$N_STEPS (t=$(round(t; digits=2))s)"
    end
end

save_log(logger, "tmp_run")
syslog = load_log("tmp_run")
scene = replay(syslog, sam.sys_struct;
               autoplay=false, loop=true)
display(scene)
