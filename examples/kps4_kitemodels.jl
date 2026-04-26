# Copyright (c) 2022, 2024 Uwe Fechner
# SPDX-License-Identifier: MIT

using Printf
# simple parking test without changing the control input
# shows how to log, plot, and print the simulation results

using KiteModels, LinearAlgebra
using KiteUtils: Settings

pkg_root = dirname(@__DIR__)
set_data_path(joinpath(pkg_root, "data", "kps4"))

set::Settings = Settings("system.yaml")

set.abs_tol=0.0006
set.rel_tol=0.00001

# the following values can be changed to match your interest
dt = 0.05
set.solver="DFBDF" # IDA or DFBDF
STEPS = 600
const PLOT = false
FRONT_VIEW = false
ZOOM = true
PRINT = false
STATISTIC = false
UPWIND_DIR = -pi/2 +deg2rad(10)
# end of user parameter section #

kcu = KCU(set)
kps4 = KPS4(kcu)

if PLOT
    using Pkg
    if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
        Pkg.activate("examples")
    end
    using ControlPlots
end

logger = Logger(set.segments + 5, STEPS)

function simulate(integrator, steps, plot=false)
    iter = 0
    sim_time = steps * dt
    println("Simulating $(sim_time)s parking flight...")
    total_elapsed = 0.0
    for i in 1:steps
        if PRINT
            lift, drag = KiteModels.lift_drag(kps4)
            @printf "%.2f: " round(integrator.t, digits=2)
            println("lift, drag  [N]: $(round(lift, digits=2)), $(round(drag, digits=2))")
        end

        total_elapsed += @elapsed next_step!(
            kps4, integrator;
            set_speed=0, upwind_dir=UPWIND_DIR, dt=dt
        )
        iter += kps4.iter
        if plot
            reltime = i*dt-dt
            if mod(i, 5) == 1
                plot2d(kps4.pos, reltime;
                    zoom=ZOOM, xlim=(40,60), front=FRONT_VIEW,
                    segments=set.segments,
                    fig="upwind_dir = $(rad2deg(UPWIND_DIR)) °")
            end
        end
        sys_state = SysState(kps4)
        sys_state.var_01 = kps4.pitch
        sys_state.var_02 = kps4.pitch_rate
        log!(logger, sys_state)
    end
    avg_rt = sim_time / total_elapsed
    println("Avg realtime factor: $(round(avg_rt, digits=2))")
    iter / steps
end
kps4.set.upwind_dir = rad2deg(UPWIND_DIR)
integrator = KiteModels.init!(kps4; delta=0.001,
    stiffness_factor=0.1, prn=STATISTIC,
    steady_state=false)

Rkw = hcat(kps4.x, kps4.y, kps4.z)
va_w = kps4.v_wind - kps4.vel_kite
va_b = round.(Rkw' * va_w, digits=2)
println("Initial: aoa=[$(round(kps4.alpha_2,
    digits=2)), $(round(kps4.alpha_3,
    digits=2)), $(round(kps4.alpha_4,
    digits=2))]°, elevation=$(round(rad2deg(
    calc_elevation(kps4)), digits=2))°, " *
    "va_b=$(va_b), " *
    "v_wind=$(round.(kps4.v_wind, digits=2))")
height = calc_height(kps4)
wf = calc_wind_factor(kps4.am, height)
println("height=$(round(height, digits=2))m, " *
    "wind_factor=$(round(wf, digits=4)), " *
    "v_wind_gnd=$(round.(kps4.v_wind_gnd, digits=2))")
println("R_b_to_w:\n  x=$(round.(kps4.x, digits=3))" *
    "\n  y=$(round.(kps4.y, digits=3))" *
    "\n  z=$(round.(kps4.z, digits=3))")

av_steps = if PLOT
    local av_steps_local = simulate(integrator, STEPS, true)
    local flight_log = KiteUtils.sys_log(logger)
    local p = plotx(flight_log.syslog.time, flight_log.z,
              rad2deg.(flight_log.syslog.var_01), rad2deg.(flight_log.syslog.var_02);
              xlabel="time [s]", ylabels=["z [m]", "pitch [°]", "pitch_rate [°/s]"],
              fig="plot_height_pitch")
    display(p)
    av_steps_local
else
    simulate(integrator, STEPS)
end
lift, drag = KiteModels.lift_drag(kps4)
println("Ground wind speed: $(round(set.v_wind, digits=2)) m/s")
println("lift, drag  [N]: $(round(lift, digits=2)), $(round(drag, digits=2))")
println("Average number of callbacks per time step: $(round(av_steps, digits=2))")

points = KiteUtils.get_particles(set.height_k, set.h_bridle, set.width, set.m_k, [0, 0, 0], [0, 0, -1], [10, 0, 0])
pos_A = points[3]
pos_C = points[5]
pos_D = points[6]
Pc = 0.5*(pos_C + pos_D)
distance = norm(pos_A-Pc)
println("Distance between A and Pc: $(round(distance, digits=2)) m")
