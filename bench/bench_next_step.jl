# Per-phase benchmark of `next_step!` on the 2-plate kite model.
#
# Replays the exact `next_step!` body (set_values, pre-step param sync,
# integrator step!, update_sys_struct!, vsm refresh) one phase at a time so
# every phase gets its own median, then appends a markdown block to a results
# file. Builds the model from the *currently loaded* SymbolicAWEModels source,
# so checking out a tag in a worktree and running this benchmarks that tag's
# code. The 2-plate kite runs a full nonlinear VSM solve each step, so the vsm
# refresh phase is exercised (vsm_interval=1).
#
# Run directly (current source, fast — uses the cached model):
#   juliaserver run examples bench/bench_next_step.jl -o
# Run against a git ref (own worktree + precompile):
#   bench/run_on_ref.sh v0.12.0
#
# Override label / output / date / step count via env:
#   BENCH_REF, BENCH_RESULTS, BENCH_DATE, BENCH_STEPS.

using Pkg
if Base.active_project() != joinpath(dirname(@__DIR__), "examples", "Project.toml")
    Pkg.activate(joinpath(dirname(@__DIR__), "examples"))
end

using SymbolicAWEModels, VortexStepMethod
using SymbolicAWEModels: sync_params!, refresh_aero!, has_vsm_wing,
    update_sys_struct!
using KiteUtils: init!, next_step!
using Statistics
using Printf

step_integrator!(integ, dt) =
    SymbolicAWEModels.OrdinaryDiffEqCore.step!(integ, dt, true)

model_name = "2plate_kite"
set_data_path(joinpath(dirname(@__DIR__), "data", model_name))

dt = 0.05
n_samples = parse(Int, get(ENV, "BENCH_STEPS", "400"))
vsm_interval = 1

struc_yaml = joinpath(get_data_path(), "particle_structural_geometry.yaml")
set = Settings("system.yaml")
set.g_earth = 0.0
vsm_set = VortexStepMethod.VSMSettings(
    joinpath(get_data_path(), "vsm_settings.yaml"); data_prefix=false)

# ContinuousAero needs the BILLOWING distribution; 0% billowing keeps it flat.
for wing_settings in vsm_set.wings
    wing_settings.spanwise_panel_distribution = VortexStepMethod.BILLOWING
    wing_settings.billowing_percentage = 0.0
end

sys = load_sys_struct_from_yaml(struc_yaml; system_name=model_name, set, vsm_set,
    aero_mode=ContinuousAero())
sys.winches[:main_winch].brake = true
sam = SymbolicAWEModel(set, sys)
init!(sam; remake=false, remake_vsm=false, prn=false)
find_steady_state!(sam; dt)

prob = sam.prob
integrator = sam.integrator
sys_struct = sam.sys_struct
set_values = [winch.set_value for winch in sys_struct.winches]
vsm_min_wind = 0.5

"Run one full `next_step!`, returning per-phase elapsed times (seconds)."
function timed_step!(step_idx)
    set_t = @elapsed prob.set_set_values(integrator, set_values)
    presync_t = @elapsed sync_params!(prob.param_sync, integrator, sys_struct)
    step_t = @elapsed step_integrator!(integrator, dt)
    update_t = @elapsed update_sys_struct!(prob, integrator, sys_struct)
    run_vsm = vsm_interval != 0 && step_idx % vsm_interval == 0 &&
        has_vsm_wing(sys_struct)
    vsm_t = run_vsm ? @elapsed(begin
        refresh_aero!(sam; vsm_min_wind)
        sync_params!(prob.param_sync, integrator, sys_struct)
    end) : NaN
    return (; set_t, presync_t, step_t, update_t, vsm_t)
end

for warmup_idx in 1:20
    timed_step!(warmup_idx)
end

set_s = Float64[]; presync_s = Float64[]; step_s = Float64[]
update_s = Float64[]; vsm_s = Float64[]; total_s = Float64[]
for sample_idx in 1:n_samples
    phase = timed_step!(sample_idx)
    push!(set_s, phase.set_t); push!(presync_s, phase.presync_t)
    push!(step_s, phase.step_t); push!(update_s, phase.update_t)
    isnan(phase.vsm_t) || push!(vsm_s, phase.vsm_t)
    push!(total_s, phase.set_t + phase.presync_t + phase.step_t +
        phase.update_t + (isnan(phase.vsm_t) ? 0.0 : phase.vsm_t))
end

for _ in 1:3
    next_step!(sam; dt, set_values, vsm_interval)
end
allocs = @allocated next_step!(sam; dt, set_values, vsm_interval)

us(samples) = isempty(samples) ? NaN : median(samples) * 1e6
total_us = us(total_s)
realtime = dt / (total_us * 1e-6)

phase_rows = [
    ("set_values (set_set_values)", us(set_s)),
    ("param sync (pre-step)",       us(presync_s)),
    ("integrator step!",            us(step_s)),
    ("update_sys_struct! (getter)", us(update_s)),
    ("vsm refresh (when run)",      us(vsm_s))]

ref = get(ENV, "BENCH_REF", strip(read(`git describe --tags --always --dirty`, String)))
commit = strip(read(`git rev-parse --short HEAD`, String))
date = get(ENV, "BENCH_DATE", strip(read(`git log -1 --format=%cs`, String)))
results_path = get(ENV, "BENCH_RESULTS", joinpath(@__DIR__, "results.md"))

fmt(x) = isnan(x) ? "n/a" : @sprintf("%.1f", x)
pct(x) = isnan(x) ? "—" : @sprintf("%.0f%%", 100 * x / total_us)

io = IOBuffer()
println(io, "## $ref — $commit — $date")
println(io)
println(io, "model: $model_name (particle / ContinuousAero, 0% billowing), ",
    "dt=", @sprintf("%.5f", dt), ", samples=$n_samples, ",
    "vsm_interval=$vsm_interval")
println(io)
println(io, "| phase | median (µs) | % of next_step! |")
println(io, "|---|---:|---:|")
for (name, value) in phase_rows
    println(io, "| $name | $(fmt(value)) | $(pct(value)) |")
end
println(io, "| **next_step! total** | **$(fmt(total_us))** | **100%** |")
println(io)
println(io, "allocations: $allocs B/step · realtime factor: ",
    @sprintf("%.1f", realtime), "×")
println(io)
block = String(take!(io))

print(block)

open(results_path, "a") do file
    write(file, block)
end
println("appended to $results_path")
