# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# The two hooks the package's `init!`/`next_step!` path needs from a backend.

"""
    build_prob!(::KernelBackend, sam; sparse=false, prn=true)

Assemble `sam.sys_struct` into a [`KernelModel`](@ref) and wrap its right-hand side as an
`ODEProblem`. The problem's parameter object is our own flat buffer plus the
callable store, so `sync_params!` writes struct fields straight into it. `sparse`
hands the solver [`state_sparsity`](@ref) as the Jacobian prototype; without it the
Jacobian is dense, as the monolith's is. `FullSpecialize` because the right-hand
side is one concrete type, so `SciMLBase`'s function wrappers would only add
indirection and allocate — it goes on the `ODEFunction` as well as the problem,
since a bare `ODEFunction` is `AutoSpecialize` and the problem's parameter does
not reach a function that is already built.
"""
function build_prob!(::KernelBackend, sam; sparse = false, prn = true)
    time = @elapsed model = assemble(sam; verbose = prn)
    if prn
        println("\tAssembled $(length(model.system.kernels)) kernels over " *
                "$(length(model.system.instances)) instances in $time seconds.")
        print_kernel_times(model.system)
    end
    rhs = KernelRHS(model.system)
    step = SimFloat(1 / sam.set.sample_freq)
    prototype = sparse ? SimFloat.(model.system.sparsity) : nothing
    problem = ODEProblem{true, SciMLBase.FullSpecialize}(
        ODEFunction{true, SciMLBase.FullSpecialize}(
            rhs; mass_matrix = model.system.mass_matrix,
            jac_prototype = prototype),
        model.u0, (0.0, step), model.params)
    sync_params!(model.param_sync, problem, sam.sys_struct)
    sam.prob = ProbWithAttributes(; prob = problem,
        param_sync = model.param_sync, initial_sync = KernelInitialSync(model),
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
assemble the problem from it, and build the integrator. The assembled problem is
serialized to `get_model_name`'s bin and read back on the next `init!` of the same
structure, exactly as the monolith's is — the kernels' `mtkcompile` is the slow
part and it is what the bin saves. `remake` forces a rebuild, `reload` re-reads the
bin over an in-memory build. The struct's positions and rest lengths are pushed
onto a reused problem by [`KernelInitialSync`](@ref), so a bin stays valid across
them; a parameter no registry reader syncs keeps the value the struct had when the
bin was built, as on the monolith.
"""
function init_backend!(::KernelBackend, sam, solver;
        adaptive = true, prn = true, reinit_sys = true, reset_vel = true,
        ignore_l0 = false, apply_tether_lengths = true, remake_vsm = true,
        reset_integrator = true, vsm_min_wind = 0.5, lin_vsm = true,
        sparse = false, remake = false, reload = false)
    model_name = get_model_name(sam.set, sam.sys_struct; sparse,
                                backend = KernelBackend())
    model_path = joinpath(KiteUtils.get_data_path(), model_name)
    prn && @info "Model bin name: $model_name"
    loaded = load_serialized_model!(sam, model_path; remake, reload, prn)
    if reinit_sys
        reinit!(sam.sys_struct, sam.set;
                ignore_l0, remake_vsm, reset_vel, apply_tether_lengths, prn)
    end
    write_total_mass!(sam.sys_struct)
    if loaded && !isnothing(sam.prob)
        prn && @info "Reusing the assembled model, so no kernel is compiled here."
    else
        prn && @info "No model to reuse; assembling and compiling the kernels."
        build_prob!(KernelBackend(), sam; sparse, prn)
        write_time = @elapsed serialize(model_path, sam.serialized_model)
        size = @sprintf("%.2f s, %.1f MB", write_time,
                        filesize(model_path) / 2^20)
        prn && @info "Wrote the assembled model ($size) to: \n\t$model_path"
    end
    integrator, _ = reinit!(sam, sam.prob, solver;
                            adaptive, reset_integrator, lin_vsm, vsm_min_wind, prn)
    return integrator
end
