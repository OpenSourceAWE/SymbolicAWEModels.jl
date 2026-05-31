# Copyright (c) 2025 2-Plate Kite Coupled Simulation
# SPDX-License-Identifier: LGPL-3.0-only

"""
2-Plate kite: coupled simulation with full nonlinear VSM solve
each step (PARTICLE_DYNAMICS wing). Counterpart to
coupled_2plate_kite_linear_vsm.jl for performance comparison.
"""

using Pkg
if Base.active_project() != joinpath(@__DIR__, "Project.toml")
    Pkg.activate(@__DIR__)
end

using GLMakie
using LinearAlgebra
using KiteUtils: init!, next_step!, update_sys_state!
using SymbolicAWEModels, VortexStepMethod

MODEL_NAME = "2plate_kite"
SIM_TIME = 10.0
VSM_INTERVAL = 1
RAMP_TIME = 2.0
STEERING_MAGNITUDE = 0.1

pkg_root = dirname(@__DIR__)
set_data_path(joinpath(pkg_root, "data", MODEL_NAME))

struc_yaml = joinpath(get_data_path(),
                      "particle_structural_geometry.yaml")

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
init!(sam; remake=false, remake_vsm=false)
find_steady_state!(sam; dt=1/300)

dt = 1 / 300
n_steps = max(1, round(Int, SIM_TIME / dt))
info_every = max(1, div(n_steps, 10))

logger = Logger(sam, n_steps)
sys_state = SysState(sam)

elapsed = @elapsed for step in 1:n_steps
    t = step * dt
    ramp = clamp(t / RAMP_TIME, 0.0, 1.0)
    steer = STEERING_MAGNITUDE * ramp
    sam.sys_struct.segments[:kcu_steering_left].l0 =
        l0_left - steer
    sam.sys_struct.segments[:kcu_steering_right].l0 =
        l0_right + steer

    next_step!(sam; dt, vsm_interval=VSM_INTERVAL)
    update_sys_state!(sys_state, sam)
    sys_state.time = t
    log!(logger, sys_state)
    if step % info_every == 0 || step == n_steps
        @info "Step $step/$n_steps (t=$(round(t; digits=2))s)"
    end
end
sim_time = n_steps * dt
@info "Realtime factor: $(round(sim_time / elapsed; digits=2))x" *
      " ($(round(elapsed; digits=2))s wall, " *
      "$(round(sim_time; digits=2))s sim, " *
      "$(round(1e3 * elapsed / n_steps; digits=2)) ms/step)"

save_log(logger, "nonlin_vsm")
syslog = load_log("nonlin_vsm")
scene = replay(syslog, sam.sys_struct)
# display(scene)
@info "Done (nonlinear VSM, interval=$VSM_INTERVAL)"
