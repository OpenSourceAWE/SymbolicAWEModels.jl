# SPDX-FileCopyrightText: 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

using SymbolicAWEModels, KiteUtils, ModelingToolkit, BenchmarkTools, Plots, Printf
using OrdinaryDiffEqCore, OrdinaryDiffEqBDF
using Suppressor, UnPack

# --- Setup ---

# Ensure data path is set correctly, copy settings if not present
set_data_path("data")
if !ispath("data/system.yaml")
    SymbolicAWEModels.copy_model_settings()
end

"""
    build_and_time(segments; solver=FBDF())

Build, compile, and solve a SymbolicAWEModel, returning timings for each major stage.
"""
function build_and_time(segments; solver=FBDF())
    set = Settings("system.yaml")
    set.segments = segments
    set.physical_model = "ram"
    sam = SymbolicAWEModel(set)
    @unpack wings, winches = sam.sys_struct
    [winch.brake=true for winch in winches]
    [wing.fix_sphere=true for wing in wings]

    suffix = segments > 1 ? "s" : ""
    println("Creating sys with $segments segment$suffix...")
    t_creation = @elapsed begin
        inputs = SymbolicAWEModels.create_sys!(sam, sam.sys_struct; prn=false)
    end
    
    println("Simplifying sys with $segments segments...")
    t_compilation = @elapsed begin
        sys = @suppress_err mtkcompile(sam.full_sys; inputs)
    end
    
    n_states = length(unknowns(sys))
    println("System has $n_states unknowns.")

    println("Benchmarking operations for $segments segments...")
    prob = ODEProblem(sys, sam.defaults, (0.0, 10.0); sam.guesses)
    
    # Benchmark init
    b_init = @benchmark OrdinaryDiffEqCore.init($prob, $solver; 
                                                reltol=1e-4, abstol=1e-4, 
                                                save_on=false) samples=100
    t_init = median(b_init).time * 1e-9
    integ = OrdinaryDiffEqCore.init(prob, solver;
                                    reltol=1e-4, abstol=1e-4, 
                                    save_on=false)

    # Benchmark reinit!
    dt = 1/sam.set.sample_freq
    b_reinit = @benchmark(OrdinaryDiffEqCore.reinit!($integ), 
                          setup=(OrdinaryDiffEqCore.step!($integ, $dt, true)), 
                          samples=100)
    t_reinit = median(b_reinit).time * 1e-9

    # Benchmark 10dt step!
    function step!(integ, dt)
        for _ in 1:10
            OrdinaryDiffEqCore.step!(integ, dt, true)
        end
    end
    b_step = @benchmark(step!($integ, $dt), 
                        setup=(OrdinaryDiffEqCore.reinit!($integ)),
                        samples=50)
    t_step = median(b_step).time * 1e-9

    # Benchmark 10dt solve
    b_solve = @benchmark solve($prob, $solver; 
                               reltol=1e-4, abstol=1e-4, 
                               saveat=$dt, tspan=(0.0, $(10dt))) samples=50
    t_solve = median(b_solve).time * 1e-9

    println("-"^20)

    return (n_states, t_creation, t_compilation, t_init, t_reinit, t_step, t_solve)
end

# --- Benchmarking Loop ---

# Define the problem sizes (number of tether segments)
N_segments = [1,2,3,4,5,6];
n_runs = length(N_segments)

# Arrays to store timings
n_states_vec = zeros(Int, n_runs)
creation_times = zeros(n_runs)
compilation_times = zeros(n_runs)
init_times = zeros(n_runs)
reinit_times = zeros(n_runs)
step_times = zeros(n_runs)
solve_times = zeros(n_runs)
total_times = zeros(n_runs)

println("--- Starting SymbolicAWEModels.jl Benchmark ---")

# Run once for precompilation
println("Precompilation run...")
build_and_time(1) 
println("Precompilation complete.")
println("-"^120)

# Print table header
@printf "%-10s | %-10s | %-10s | %-12s | %-10s | %-10s | %-10s | %-12s | %-10s\n" "Segments" "Unknowns" "Create (s)" "Compile (s)" "Init (s)" "Reinit (s)" "Step (s)" "Solve (s)" "Total (s)"
println(repeat("-", 120))

# Main benchmark loop
for (i, n) in enumerate(N_segments)
    n_states, t_creation, t_compilation, t_init, t_reinit, t_step, t_solve = build_and_time(n)
    n_states_vec[i] = n_states
    creation_times[i] = t_creation
    compilation_times[i] = t_compilation
    init_times[i] = t_init
    reinit_times[i] = t_reinit
    step_times[i] = t_step
    solve_times[i] = t_solve
    total_times[i] = t_creation + t_compilation + t_solve

    # Print results row for the current run
    @printf "%-10d | %-10d | %-10.3f | %-12.3f | %-10.3f | %-10.3e | %-10.3e | %-12.3f | %-10.3f\n" n n_states t_creation t_compilation t_init t_reinit t_step t_solve total_times[i]
end

println(repeat("-", 120))
println("--- Benchmark Complete ---")

# --- Generate Final Plots ---

p = plot(n_states_vec, creation_times, label="Model Creation (create_sys!)",
         xaxis=:linear, yaxis=:log,
         xlabel="Number of Unknowns",
         ylabel="Time (s)",
         title="SymbolicAWEModels.jl Performance",
         legend=:topleft,
         size=(900,600),
         marker=:circle)

plot!(p, n_states_vec, compilation_times, label="Symbolic Compilation (mtkcompile)", marker=:circle)
plot!(p, n_states_vec, init_times, label="Integrator Creation (init)", marker=:circle)
plot!(p, n_states_vec, reinit_times, label="Reinitialize (reinit!)", marker=:circle)
plot!(p, n_states_vec, step_times, label="Single Step (step!)", marker=:circle)
plot!(p, n_states_vec, solve_times, label="Solve Time (1s sim)", marker=:circle)
plot!(p, n_states_vec, total_times, label="Total Initial Time", linewidth=3, linestyle=:dash)

display(p)
savefig(joinpath(get_data_path(), "symbolic_awe_benchmark.png"))

# --- Appendix: Computer Information ---
println("\n\n--- Appendix ---")
println("Computer Information:\n")
versioninfo()
println("\nPackage Information:\n")
using Pkg
Pkg.status(["SymbolicAWEModels", "ModelingToolkit", "OrdinaryDiffEqCore", "Plots", "BenchmarkTools"], mode=Pkg.PKGMODE_MANIFEST)
