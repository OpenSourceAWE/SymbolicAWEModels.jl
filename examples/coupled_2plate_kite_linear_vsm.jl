# Copyright (c) 2025
# SPDX-License-Identifier: MPL-2.0

"""
2-Plate kite: coupled simulation using linearized VSM updates
(re-linearize every few steps instead of solving full VSM).
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
VSM_INTERVAL = 100   # re-linearise every step while debugging

pkg_root = dirname(@__DIR__)
set_data_path(joinpath(pkg_root, "data", MODEL_NAME))

struc_yaml = joinpath(get_data_path(),
                      "quat_struc_geometry.yaml")

set = Settings("system.yaml")
set.g_earth = 0.0
vsm_set = VortexStepMethod.VSMSettings(
    joinpath(get_data_path(), "vsm_settings.yaml");
    data_prefix=false)

sys = load_sys_struct_from_yaml(struc_yaml;
    system_name=MODEL_NAME, set, vsm_set)
sys.winches[:main_winch].brake = true
for wing in sys.wings
    wing.aero_mode = AERO_LINEARIZED
end
sam = SymbolicAWEModel(set, sys)
init!(sam; remake=false, remake_vsm=false)
find_steady_state!(sam)

dt = 0.001
n_steps = max(1, round(Int, SIM_TIME / dt))

logger = Logger(sam, n_steps)
sys_state = SysState(sam)

# Ring buffer of last N step diagnostics for crash forensics.
DIAG_KEEP = 8
diag_history = Vector{NamedTuple}()

function snapshot_wing(wing, t, step)
    (
        step = step,
        t = t,
        va_b = copy(wing.va_b),
        va_norm = norm(wing.va_b),
        omega_b = copy(wing.ω_b),
        aero_y = copy(wing.aero_y),
        aero_x = copy(wing.aero_x),
        jac = copy(wing.aero_jac),
        jac_max = maximum(abs, wing.aero_jac),
        jac_norm = sqrt(sum(abs2, wing.aero_jac)),
        force_b = copy(wing.aero_force_b),
        moment_b = copy(wing.aero_moment_b),
        pos_w = copy(wing.pos_w),
        com_vel = copy(wing.com_vel),
    )
end

function push_diag!(hist, snap)
    push!(hist, snap)
    length(hist) > DIAG_KEEP && popfirst!(hist)
end

function print_diag(snap)
    @info "step $(snap.step) t=$(round(snap.t; digits=4))" *
          " |va|=$(round(snap.va_norm; digits=3))" *
          " jac_max=$(round(snap.jac_max; digits=3))"
    println("  aero_y     = ", round.(snap.aero_y; digits=4))
    println("  aero_x     = ", round.(snap.aero_x; digits=4))
    println("  va_b       = ", round.(snap.va_b; digits=4))
    println("  omega_b    = ", round.(snap.omega_b; digits=4))
    println("  force_b    = ", round.(snap.force_b; digits=2))
    println("  moment_b   = ", round.(snap.moment_b; digits=2))
    println("  pos_w      = ", round.(snap.pos_w; digits=3))
    println("  com_vel    = ", round.(snap.com_vel; digits=3))
end

last_completed_step = 0
for step in 1:n_steps
    try
        next_step!(sam; dt, vsm_interval=VSM_INTERVAL)
    catch err
        @warn "next_step! failed at step $step (t=$(round(
            step * dt; digits=4))s)" exception=err
        @info "── last $(length(diag_history)) successful steps ──"
        for snap in diag_history
            print_diag(snap)
        end
        break
    end
    global last_completed_step = step
    snap = snapshot_wing(sam.sys_struct.wings[1], step * dt, step)
    push_diag!(diag_history, snap)
    if snap.jac_max > 1e6
        @warn "Jacobian explosion at step $step (t=$(round(
            snap.t; digits=4))s) — jac_max=$(snap.jac_max)"
        println("aero_y at this step: ", snap.aero_y)
        println("aero_x at this step: ", snap.aero_x)
        println("Full Jacobian (rows = ∂x/∂y, " *
                "rows: CL,CD,CS,CMx,CMy,CMz,cm_g1,cm_g2,cm_g3 ; " *
                "cols: α,β,ω1,ω2,ω3,θ_g1,θ_g2,θ_g3):")
        for (ix, row) in enumerate(eachrow(snap.jac))
            println("  row ", ix, ": ", row)
        end
        break
    end
    update_sys_state!(sys_state, sam)
    sys_state.time = step * dt
    log!(logger, sys_state)
    if step % 100 == 0 || step == n_steps
        @info "Step $step/$n_steps" t=round(
            step * dt; digits=2) jac_max=round(
            snap.jac_max; digits=2) force=round.(
            snap.force_b; digits=1)
    end
end

save_log(logger, "linear_vsm")
syslog = load_log("linear_vsm")
scene = replay(syslog, sam.sys_struct)
display(scene)
@info "Done (linearized VSM, interval=$VSM_INTERVAL)"
