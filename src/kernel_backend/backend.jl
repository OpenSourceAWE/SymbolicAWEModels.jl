# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# The two hooks the package's `init!`/`next_step!` path needs from a backend.

"""
    build_prob!(::KernelBackend, sam; prn=true)

Assemble `sam.sys_struct` into a [`KernelModel`](@ref) and wrap its right-hand side as an
`ODEProblem`. The problem's parameter object is our own flat buffer plus the
callable store, so `sync_params!` writes struct fields straight into it.
"""
function build_prob!(::KernelBackend, sam; prn = true)
    time = @elapsed model = assemble(sam; verbose = prn)
    if prn
        println("\tAssembled $(length(model.system.kernels)) kernels over " *
                "$(length(model.system.instances)) instances in $time seconds.")
        print_kernel_times(model.system)
    end
    write_total_mass!(sam.sys_struct)
    rhs = KernelRHS(model.system)
    step = SimFloat(1 / sam.set.sample_freq)
    # FullSpecialize: the right-hand side is one concrete type, so wrapping it in
    # SciMLBase's function wrappers would only add indirection (and allocate).
    problem = ODEProblem{true, SciMLBase.FullSpecialize}(
        ODEFunction(rhs; mass_matrix = model.system.mass_matrix),
        model.u0, (0.0, step), model.params)
    sync_params!(model.param_sync, problem, sam.sys_struct)
    sam.prob = ProbWithAttributes(; prob = problem,
        param_sync = model.param_sync, initial_sync = nothing,
        set_set_values = KernelControlSetter(model, sam.sys_struct),
        get_set_values = nothing, get_aero_input = nothing,
        get_all_state = KernelStateGetter(model, rhs, sam.sys_struct))
    return true
end

"""
    print_kernel_times(system)

List what each kernel cost to compile, dearest first, with how many instances share
it. Compilation is per kernel, so this is where a slow build has to be read: a large
model with few kernel types is cheap, a small one with many types is not.
"""
function print_kernel_times(system::KernelSystem)
    counts = zeros(Int, length(system.kernels))
    for instance in system.instances
        counts[instance.kernel] += 1
    end
    order = sortperm(system.compile_seconds; rev = true)
    for k in order
        @printf("\t  %-28s %7.2f s  ×%d\n", system.kernels[k].name,
                system.compile_seconds[k], counts[k])
    end
    return nothing
end

"""
    init_backend!(::KernelBackend, sam, solver; kwargs...)

Full `init!` path for the [`KernelBackend`](@ref): refresh the `SystemStructure`,
assemble the problem from it, and build the integrator. The problem is rebuilt on
every `init!`, so the initial state always comes from the struct.
"""
function init_backend!(::KernelBackend, sam, solver;
        adaptive = true, prn = true, reinit_sys = true, reset_vel = true,
        ignore_l0 = false, apply_tether_lengths = true, remake_vsm = true,
        reset_integrator = true, vsm_min_wind = 0.5, lin_vsm = true)
    if reinit_sys
        reinit!(sam.sys_struct, sam.set;
                ignore_l0, remake_vsm, reset_vel, apply_tether_lengths, prn)
    end
    build_prob!(KernelBackend(), sam; prn)
    integrator, _ = reinit!(sam, sam.prob, solver;
                            adaptive, reset_integrator, lin_vsm, vsm_min_wind, prn)
    return integrator
end
