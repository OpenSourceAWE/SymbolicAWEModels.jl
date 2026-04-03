# Copyright (c) 2025 Bart van de Lint, Uwe Fechner
# SPDX-License-Identifier: MPL-2.0

"""
Simple linearized model: stabilize, linearize with
`simple_linearize!`, then simulate nonlinear with steering
and compare to the linear prediction.
"""

using GLMakie
using KiteUtils, SymbolicAWEModels, VortexStepMethod
using SymbolicAWEModels: calc_steady_torque, init!, next_step!, simple_linearize!, update_sys_state!
using ControlSystemsBase

# Simulation parameters
dt = 0.05
total_time = 1.0
vsm_interval = 3
steps = Int(round(total_time / dt))
steering_magnitude = 10.0

# Initialize model
set_data_path("data/2plate_kite")
struc_yaml = joinpath(get_data_path(),
                      "quat_struc_geometry.yaml")
aero_yaml = joinpath(get_data_path(), "aero_geometry.yaml")
update_aero_yaml_from_struc_yaml!(struc_yaml, aero_yaml)

set = Settings("system.yaml")
vsm_set = VortexStepMethod.VSMSettings(
    joinpath(get_data_path(), "vsm_settings.yaml");
    data_prefix=false)

sys = load_sys_struct_from_yaml(struc_yaml;
    system_name="2plate_kite", set, vsm_set)
sam = SymbolicAWEModel(set, sys)
sam.set.abs_tol = 1e-3
sam.set.rel_tol = 1e-3

init!(sam; remake=false)
find_steady_state!(sam; t=10.0, dt=1/sam.set.sample_freq)

# Simple linearization
simple_linearize!(sam; tstab=10.0)
lin = sam.serialized_model.simple_lin_model
lin_ss = ss(lin.A, lin.B, lin.C, lin.D)

# Nonlinear simulation with steering
logger = Logger(sam, steps)
sys_state = SysState(sam)
t = 0.0
steady_torque =
    calc_steady_torque(sam)
torque_damp = 0.9
u0 = copy(steady_torque)
set_values_mat = zeros(3, steps)

for i in 1:steps
    t = i * dt
    prev_steady_torque = steady_torque
    steady_torque = torque_damp * prev_steady_torque +
        (1 - torque_damp) *
        calc_steady_torque(sam)
    sign_val = t > 0.5 ? -1 : 1
    sv = steady_torque .+ sign_val .*
        [10.0, steering_magnitude, -steering_magnitude]
    set_values_mat[:, i] = sv

    next_step!(sam; set_values=sv, dt, vsm_interval)
    update_sys_state!(sys_state, sam)

    sys_state.time = t
    log!(logger, sys_state)
end

save_log(logger, "tmp")
lg = load_log("tmp")

# Linear simulation for comparison
lin_res = lsim(lin_ss, set_values_mat .- u0,
               lg.syslog.time)

@info "Simulation completed" steps=length(
    lg.syslog.time) final_heading=round(
    rad2deg(lg.syslog.heading[end]); digits=2)
