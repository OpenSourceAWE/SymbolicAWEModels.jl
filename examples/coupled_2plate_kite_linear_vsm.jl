# Copyright (c) 2025
# SPDX-License-Identifier: MPL-2.0

"""
2-Plate kite: coupled simulation using linearized VSM updates
(re-linearize every few steps instead of solving full VSM).
"""

using Pkg
Pkg.activate(@__DIR__)

using GLMakie
using KiteUtils: init!, next_step!, update_sys_state!
using SymbolicAWEModels, VortexStepMethod

MODEL_NAME = "2plate_kite"
SIM_TIME = 10.0
VSM_INTERVAL = 3

pkg_root = dirname(@__DIR__)
set_data_path(joinpath(pkg_root, "data", MODEL_NAME))

struc_yaml = joinpath(get_data_path(),
                      "quat_struc_geometry.yaml")
aero_yaml = joinpath(get_data_path(), "aero_geometry.yaml")
update_aero_yaml_from_struc_yaml!(struc_yaml, aero_yaml)

set = Settings("system.yaml")
vsm_set = VortexStepMethod.VSMSettings(
    joinpath(get_data_path(), "vsm_settings.yaml");
    data_prefix=false)

sys = load_sys_struct_from_yaml(struc_yaml;
    system_name=MODEL_NAME, set, vsm_set)
sam::SymbolicAWEModel = SymbolicAWEModel(set, sys)
init!(sam; remake=false)

dt = 1.0 / set.sample_freq
n_steps = max(1, round(Int, SIM_TIME / dt))

logger = Logger(sam, n_steps)
sys_state = SysState(sam)

for step in 1:n_steps
    next_step!(sam; dt, vsm_interval=VSM_INTERVAL)
    update_sys_state!(sys_state, sam)
    sys_state.time = step * dt
    log!(logger, sys_state)
    if step % 100 == 0 || step == n_steps
        @info "Step $step/$n_steps" t=round(
            step * dt; digits=2)
    end
end

save_log(logger, "linear_vsm")
syslog = load_log("linear_vsm")
scene = replay(syslog, sam.sys_struct)
display(scene)
@info "Done (linearized VSM, interval=$VSM_INTERVAL)"
