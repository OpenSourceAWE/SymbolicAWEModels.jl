# Copyright (c) 2025 Bart van de Lint and Uwe Fechner
# SPDX-License-Identifier: LGPL-3.0-only

"""
    generate_prob_getters(sys_struct, sys)

Generate getter and setter functions for the state variables of the full system model.

These functions provide a convenient way to access and modify the state and parameters
of the compiled `ODESystem` (`sys`).

# Arguments
- `sys_struct::SystemStructure`: The structure defining the system topology.
- `sys::ODESystem`: The compiled ModelingToolkit system.

# Returns
- A `NamedTuple` containing various getter and setter functions for different parts of the system state.
"""
function generate_prob_getters(sys_struct, sys, param_registry=nothing,
                               initial_registry=nothing)
    (; points, wings, twist_surfaces, pulleys, winches, tethers, segments, bodies) = sys_struct
    get_aero_input, set_set_values, get_set_values = nothing, nothing, nothing

    specs = NamedTuple[]
    if length(points) > 0
        push!(specs, scatter_spec(ss -> ss.points,
            sys.pos         => (c, v) -> copy_vec!(c.pos_w, v, c.idx),
            sys.vel         => (c, v) -> copy_vec!(c.vel_w, v, c.idx),
            sys.point_force => (c, v) -> copy_vec!(c.force, v, c.idx),
            sys.va_point_b  => (c, v) -> copy_vec!(c.va_b, v, c.idx),
            sys.point_mass  => (c, v) -> (c.total_mass = v[c.idx]; nothing),
            sys.total_drag  => (c, v) -> copy_vec!(c.drag_force, v, c.idx)))
    end
    if length(pulleys) > 0
        push!(specs, scatter_spec(ss -> ss.pulleys,
            sys.pulley_len => (c, v) -> (c.len = v[c.idx]; nothing),
            sys.pulley_vel => (c, v) -> (c.vel = v[c.idx]; nothing)))
    end
    if length(segments) > 0
        push!(specs, scatter_spec(ss -> ss.segments,
            sys.spring_force => (c, v) -> (c.force = v[c.idx]; nothing),
            sys.len          => (c, v) -> (c.len = v[c.idx]; nothing),
            sys.l0           => (c, v) -> (c.l0 = v[c.idx]; nothing)))
    end
    if length(twist_surfaces) > 0
        push!(specs, scatter_spec(ss -> ss.twist_surfaces,
            sys.twist_angle                 => (c, v) -> (c.twist = v[c.idx]; nothing),
            sys.twist_ω                     => (c, v) -> (c.twist_ω = v[c.idx]; nothing),
            sys.twist_surface_tether_force  => (c, v) -> (c.tether_force = v[c.idx]; nothing),
            sys.twist_surface_tether_moment => (c, v) -> (c.tether_moment = v[c.idx]; nothing),
            sys.twist_surface_aero_moment   => (c, v) -> (c.aero_moment = v[c.idx]; nothing)))
    end
    if length(winches) > 0
        push!(specs, scatter_spec(ss -> ss.winches,
            sys.winch_acc       => (c, v) -> (c.acc = v[c.idx]; nothing),
            sys.winch_vel       => (c, v) -> (c.vel = v[c.idx]; nothing),
            sys.set_values      => (c, v) -> (c.set_value = v[c.idx]; nothing),
            sys.winch_force_vec => (c, v) -> copy_vec!(c.force, v, c.idx),
            sys.winch_friction  => (c, v) -> (c.friction = v[c.idx]; nothing)))
        set_set_values = setp(sys, sys.set_values)
        get_set_values = getp(sys, sys.set_values)
    end
    if length(tethers) > 0
        push!(specs, scatter_spec(ss -> ss.tethers,
            sys.tether_len    => (c, v) -> (c.len = v[c.idx]; nothing),
            sys.stretched_len => (c, v) -> (c.stretched_len = v[c.idx]; nothing)))
    end
    if length(wings) > 0
        # Wing rigid-body state is synced via the embedded body in the bodies spec below.
        push!(specs, scatter_spec(ss -> ss.wings,
            sys.va_wing_b        => (c, v) -> copy_vec!(c.va_b, v, c.idx),
            sys.wind_vel_wing    => (c, v) -> copy_vec!(c.v_wind, v, c.idx),
            sys.aero_force_b     => (c, v) -> copy_vec!(c.aero_force_b, v, c.idx),
            sys.aero_moment_b    => (c, v) -> copy_vec!(c.aero_moment_b, v, c.idx),
            sys.elevation        => (c, v) -> (c.elevation = v[c.idx]; nothing),
            sys.elevation_vel    => (c, v) -> (c.elevation_vel = v[c.idx]; nothing),
            sys.elevation_acc    => (c, v) -> (c.elevation_acc = v[c.idx]; nothing),
            sys.azimuth          => (c, v) -> (c.azimuth = v[c.idx]; nothing),
            sys.azimuth_vel      => (c, v) -> (c.azimuth_vel = v[c.idx]; nothing),
            sys.azimuth_acc      => (c, v) -> (c.azimuth_acc = v[c.idx]; nothing),
            sys.heading          => (c, v) -> (c.heading = v[c.idx]; nothing),
            sys.turn_rate        => (c, v) -> copy_vec!(c.turn_rate, v, c.idx),
            sys.turn_acc         => (c, v) -> copy_vec!(c.turn_acc, v, c.idx),
            sys.course           => (c, v) -> (c.course = v[c.idx]; nothing),
            sys.angle_of_attack  => (c, v) -> (c.aoa = v[c.idx]; nothing)))

        # aero_input exists only for wings whose component exposes the connector.
        aero_inputs = [
            getproperty(sys, Symbol("aero_$(wing.idx)")).aero_input
            for wing in sys_struct.wings
            if wing.dynamics_type === RIGID_DYNAMICS && hasproperty(
                getproperty(sys, Symbol("aero_$(wing.idx)")), :aero_input)]
        get_aero_input = isempty(aero_inputs) ? nothing :
            getu(sys, collect.(aero_inputs))
    end
    if length(bodies) > 0
        push!(specs, scatter_spec(ss -> ss.bodies,
            sys.body_Q_b_to_w => (c, v) -> copy_vec!(c.Q_b_to_w, v, c.idx),
            sys.body_ω_b      => (c, v) -> copy_vec!(c.ω_b, v, c.idx),
            sys.body_pos_w    => (c, v) -> copy_vec!(c.pos_w, v, c.idx),
            sys.body_vel_w    => (c, v) -> copy_vec!(c.vel_w, v, c.idx),
            sys.body_acc_w    => (c, v) -> copy_vec!(c.acc_w, v, c.idx),
            sys.body_com_w    => (c, v) -> copy_vec!(c.com_w, v, c.idx),
            sys.body_com_vel  => (c, v) -> copy_vec!(c.com_vel, v, c.idx),
            sys.body_Q_p_to_w => (c, v) -> copy_vec!(c.Q_p_to_w, v, c.idx),
            sys.body_ω_p      => (c, v) -> copy_vec!(c.ω_p, v, c.idx)))
    end

    get_all_state = build_inplace_getter(sys, specs)

    param_sync = isnothing(param_registry) ? nothing :
        build_param_sync(sys, param_registry)
    initial_sync = isnothing(initial_registry) ? nothing :
        build_initial_sync(sys, initial_registry)

    return (; get_aero_input, set_set_values, get_set_values, get_all_state,
            param_sync, initial_sync)
end

"""
    scatter_spec(selector, pairs...)

Describe one component group: `selector(sys_struct)` yields its component vector
and each `sys_array => copyfn` pair maps a symbolic output array to a
`(component, view) -> _` closure that writes that array's slice into the
component's struct field. This is the single source of truth — both the buffer
layout and the scatter derive from the same ordered list.
"""
function scatter_spec(selector, pairs::Pair...)
    return (; selector,
            arrays  = [pair.first for pair in pairs],
            copyfns = Tuple(pair.second for pair in pairs))
end

"""
    build_grouped_views(buf, group_shapes)

Build a tuple (per group) of tuples (per output array) of zero-copy reshaped
views into `buf`. Flat order is groups-in-order, arrays-in-order, column-major
within each array — shared by [`build_inplace_getter`](@ref) and deserialization
so layouts always match.
"""
function build_grouped_views(buf, group_shapes::Tuple)
    offset = 0
    return map(group_shapes) do shapes
        Tuple(map(shapes) do shp
            n = prod(shp)
            rng = (offset + 1):(offset + n)
            offset += n
            reshape(view(buf, rng), shp)
        end)
    end
end

"""
    build_inplace_getter(sys, specs)

Build one [`InplaceGetter`](@ref) from the per-group `specs` (see
[`scatter_spec`](@ref)). All component output arrays are concatenated into a
single MTK in-place observed function so shared subexpressions (e.g. the
spring/force network) are computed once.
"""
function build_inplace_getter(sys, specs)
    arrays_per_group = [collect.(spec.arrays) for spec in specs]
    group_shapes = Tuple(Tuple(size.(arrs)) for arrs in arrays_per_group)
    flat = reduce(vcat, [vec(a) for arrs in arrays_per_group for a in arrs];
                  init = Num[])
    _, iip = ModelingToolkit.build_explicit_observed_function(
        sys, flat; return_inplace=true)
    buf = Vector{SimFloat}(undef, length(flat))
    grouped_views = build_grouped_views(buf, group_shapes)
    groups = Tuple(ScatterGroup(spec.selector, spec.copyfns, views)
                   for (spec, views) in zip(specs, grouped_views))
    return InplaceGetter(iip, buf, groups)
end

"""
    serialize(s, g::InplaceGetter)

Julia's serializer does not preserve `SubArray.parent === buf` sharing, so
serialize the generated function plus a cheap layout and rebuild the aliased
buffer and views on load.
"""
function Serialization.serialize(s::Serialization.AbstractSerializer,
                                 g::InplaceGetter)
    Serialization.serialize_type(s, typeof(g))
    serialize(s, g.fn)
    serialize(s, length(g.buf))
    serialize(s, length(g.groups))
    for group in g.groups
        serialize(s, group.selector)
        serialize(s, group.copyfns)
        serialize(s, map(size, group.views))
    end
end

function Serialization.deserialize(s::Serialization.AbstractSerializer,
                                   ::Type{<:InplaceGetter})
    fn = deserialize(s)
    n = deserialize(s)
    n_groups = deserialize(s)
    selectors = Vector{Any}(undef, n_groups)
    copyfns = Vector{Any}(undef, n_groups)
    shapes = Vector{Any}(undef, n_groups)
    for i in 1:n_groups
        selectors[i] = deserialize(s)
        copyfns[i] = deserialize(s)
        shapes[i] = deserialize(s)
    end
    buf = Vector{SimFloat}(undef, n)
    grouped_views = build_grouped_views(buf, Tuple(shapes))
    groups = Tuple(ScatterGroup(selectors[i], copyfns[i], grouped_views[i])
                   for i in 1:n_groups)
    return InplaceGetter(fn, buf, groups)
end

"""
    generate_lin_getters(sys)

Generate setter functions for the parameters of a linearized system.

# Arguments
- `sys`: The linearized ModelingToolkit system.

# Returns
- A `NamedTuple` containing the setter function for the winch set-points
  (`set_set_values`).
"""
function generate_lin_getters(sys)
    set_set_values = nothing
    if hasproperty(sys, :set_values)
        set_set_values = setp(sys, sys.set_values)
    end
    return (; set_set_values)
end

"""
    generate_control_funcs(model, inputs, outputs)

Generate in-place and out-of-place control functions from a ModelingToolkit system.

This function wraps `ModelingToolkit.generate_control_function` and
`ModelingToolkit.build_explicit_observed_function` to create the necessary functions
for simulation and analysis.

# Arguments
- `model`: The full `ODESystem`.
- `inputs`: A vector of input variables.
- `outputs`: A vector of output variables.

# Returns
- A `NamedTuple` containing the generated functions (`f_oop`, `f_ip`, `h_oop`, `h_ip`),
  system dimensions (`nu`, `nx`, `ny`), and symbolic variables (`dvs`, `psym`, `io_sys`).
"""
function generate_control_funcs(model, inputs, outputs)
    (f_ip, f_oop), dvs, psym, io_sys = @suppress_err ModelingToolkit.generate_control_function(
        model, inputs; simplify=false
    )
    num_inputs, num_states, num_outputs = length(inputs), length(dvs), length(outputs)
    (h_oop, h_ip) = ModelingToolkit.build_explicit_observed_function(
        io_sys, outputs; inputs, return_inplace=true
    )
    return (f_oop=f_oop, f_ip=f_ip, h_oop=h_oop, h_ip=h_ip,
            nu=num_inputs, nx=num_states, ny=num_outputs,
            dvs=dvs, psym=psym, io_sys=io_sys)
end

"""
    load_serialized_model!(sam, model_path; remake=false, reload=false)

Load a serialized model from disk if it is valid.

A model is considered valid if its settings and system structure hashes match the
current ones in the `SymbolicAWEModel` object (`sam`).

# Arguments
- `sam::SymbolicAWEModel`: The main model object.
- `model_path::String`: The path to the serialized model file.
- `remake::Bool`: If true, forces the model to be considered invalid, triggering a rebuild.
- `reload::Bool`: If true, forces reloading from disk even if the model is already in memory.

# Returns
- `true` if a valid model was successfully loaded into `sam.serialized_model`, `false` otherwise.
"""
function load_serialized_model!(sam, model_path; remake=false, reload=false)
    set_hash = get_set_hash(sam.set)
    sys_struct_hash = get_sys_struct_hash(sam.sys_struct)
    
    if ispath(model_path) && !remake
        if set_hash != sam.serialized_model.set_hash || sys_struct_hash != sam.serialized_model.sys_struct_hash
            sam.serialized_model = SerializedModel(; set_hash, sys_struct_hash)
            return false
        end
        if !isnothing(sam.serialized_model.full_sys) && !reload
            return true
        end
        try
            serialized_model = deserialize(model_path)
            if set_hash == serialized_model.set_hash && sys_struct_hash == serialized_model.sys_struct_hash
                sam.serialized_model = serialized_model
                return true
            end
        catch exception
            backtrace = catch_backtrace()
            log_path = model_path * ".error.log"
            open(log_path, "w") do io
                println(io,
                    "Deserialization failed at ",
                    time())
                println(io, "Path: $model_path")
                showerror(io, exception, backtrace)
            end
            @warn "Failure to deserialize " *
                "$model_path: $(typeof(exception)). " *
                "Details in $log_path"
        end
    end
    sam.serialized_model = SerializedModel(; set_hash, sys_struct_hash)
    return false
end

"""
    maybe_create_prob!(sam; create_prob=true, prn=true)

Create and cache the `ODEProblem` if it does not already exist.

This function compiles the full system, creates the `ODEProblem`, and generates
the necessary getter/setter functions.

# Arguments
- `sam::SymbolicAWEModel`: The main model object.
- `create_prob::Bool`: A flag to enable or disable the creation of the problem.
- `prn::Bool`: A flag to enable or disable printing of progress messages.

# Returns
- `true` if a new problem was created, `false` otherwise.
"""
function maybe_create_prob!(sam; create_prob=true, prn=true)
    if create_prob && isnothing(sam.prob)
        isnothing(sam.full_sys) && return false
        full_sys = something(sam.full_sys)
        local sys
        time = @elapsed @suppress_err sys = mtkcompile(full_sys; inputs=sam.inputs)
        prn && println("\tSimplified the System for ODEProblem in $time seconds.")

        dt = SimFloat(1/sam.set.sample_freq)
        time = @elapsed prob = ODEProblem(sys, sam.defaults, (0.0, dt))
        prn && println("\tCreated the ODEProblem in $time seconds.")

        time = @elapsed getters = generate_prob_getters(sam.sys_struct, sys,
            sam.param_registry, sam.initial_registry)
        prn && println("\tCreated state getters and setters in $time seconds.")

        sam.prob = ProbWithAttributes(; prob, getters...)
        return true
    end
    return false
end

"""
    maybe_create_lin_prob!(sam, outputs; ...)

Create and cache the `LinearizationProblem` if it does not exist or if the outputs have changed.

# Arguments
- `sam::SymbolicAWEModel`: The main model object.
- `outputs`: A vector of output variables for the linearization.
- `create_lin_prob::Bool`: Flag to enable/disable creation.
- `outputs_changed::Bool`: Flag indicating if the output vector has changed.
- `prn::Bool`: Flag to enable/disable printing of progress messages.

# Returns
- `true` if a new problem was created, `false` otherwise.
"""
function maybe_create_lin_prob!(sam, outputs; create_lin_prob=true, prn=true)
    if create_lin_prob && isnothing(sam.lin_prob)
        isnothing(sam.full_sys) && return false
        isnothing(outputs) && return false
        full_sys = something(sam.full_sys)
        time = @elapsed @suppress_err begin
            lin_fun, lin_sys = linearization_function(full_sys, [sam.inputs...], outputs;
                                                    op=sam.defaults)
            prob = LinearizationProblem(lin_fun, 0.0)
            getters = generate_lin_getters(lin_sys)
            sam.lin_prob = LinProbWithAttributes(; prob,
                                                getters...)
        end
        prn && println("\tCreated the LinearizationProblem in $time seconds.")
        return true
    end
    return false
end

"""
    maybe_create_control_functions!(sam, outputs; ...)

Create and cache the control functions if they do not exist or if the outputs have changed.

# Arguments
- `sam::SymbolicAWEModel`: The main model object.
- `outputs`: A vector of output variables for the control functions.
- `create_control_func::Bool`: Flag to enable/disable creation.
- `outputs_changed::Bool`: Flag indicating if the output vector has changed.
- `prn::Bool`: Flag to enable/disable printing of progress messages.

# Returns
- `true` if new functions were created, `false` otherwise.
"""
function maybe_create_control_functions!(sam, outputs; create_control_func=false, prn=true)
    if create_control_func && isnothing(sam.control_functions)
        isnothing(sam.full_sys) && return false
        isnothing(outputs) && return false
        full_sys = something(sam.full_sys)
        inputs = [sam.inputs...]
        time = @elapsed result = generate_control_funcs(full_sys, inputs, outputs)
        sam.control_functions = ControlFuncWithAttributes(; result...)
        prn && println("\tCreated the control functions in $time seconds.")
        return true
    end
    return false
end

"""
    has_custom_component(sys_struct)

Return `true` when the system has a non-default winch model or a wing using a
custom aero model, in which case the compiled model cannot be reused from cache
and must be rebuilt.
"""
function has_custom_component(sys_struct)
    any(!is_builtin_winch(winch.model)
        for winch in sys_struct.winches) && return true
    any(!is_builtin_aero(wing.aero)
        for wing in sys_struct.wings) && return true
    return false
end

"""
    init!(sam::SymbolicAWEModel; ...)

Load or build the symbolic model, create the `ODEProblem` (and optionally the
`LinearizationProblem` / control functions), serialize new builds to disk, and
return a freshly initialized `ODEIntegrator`.

# Keyword Arguments
- `solver`, `adaptive`: ODE solver and time-stepping mode. `solver=nothing` picks
  a default from `sam.set.solver`.
- `prn`: print progress messages.
- `remake`: force a full rebuild, ignoring any cached compiled model. Defaults to
  `nothing`, which rebuilds automatically when a custom winch/aero component is
  present (see [`has_custom_component`](@ref)) and reuses the cache otherwise.
- `reload`: force reloading the serialized model from disk.
- `outputs`: vector of output variables (used by linearization / control funcs).
- `create_prob`, `create_lin_prob`, `create_control_func`: which artefacts to build.
- `lin_vsm`: linearize the VSM aerodynamics after init.
- `remake_vsm`: rebuild the VSM wing/aero from settings (after editing
  `aero_geometry.yaml` etc.).
- `reset_vel`, `ignore_l0`: forwarded to `reinit!(sys_struct, set)`.
- `reinit_sys`: run `reinit!(sys_struct, set)` to refresh positions, lengths, and
  transforms. Set to `false` to preserve manual adjustments to the
  `SystemStructure` (or after calling `reinit!(sys_struct, set; …)` yourself).
- `reset_integrator`: discard the existing integrator and build a fresh one. Use
  when stale BDF history would taint the next run.
- `vsm_min_wind=0.5`: minimum |va| [m/s] for the initial VSM solve. Below this the
  solve is skipped and the wing's aero outputs are zeroed (the solver fails to
  converge / the Jacobian blows up as 1/|va|).
"""
function init!(sam::SymbolicAWEModel;
    solver=nothing, adaptive=true, prn=true,
    remake=nothing, reload=false,
    outputs=nothing,
    create_prob::Bool=true,
    create_lin_prob::Bool=false,
    create_control_func::Bool=false,
    lin_vsm::Bool=true,
    ignore_l0::Bool=false,
    remake_vsm::Bool=true,
    reset_vel::Bool=true,
    reset_integrator::Bool=true,
    reinit_sys::Bool=true,
    apply_tether_lengths::Bool=true,
    vsm_min_wind=0.5
)
    prn && @info "Initializing $(sam.sys_struct.name) model..."
    sam.sys_struct isa SystemStructure || error(
        "Equation generation requires SystemStructure, " *
        "got $(typeof(sam.sys_struct)).")
    if isnothing(remake)
        remake = has_custom_component(sam.sys_struct)
        remake && prn && @info "Custom winch/aero model detected; " *
            "forcing remake (custom component equations are not " *
            "captured by the model hash)."
    end
    time = @elapsed begin
        if isnothing(solver)
            if sam.set.solver == "QNDF"
                @warn "This solver is not tested."
                solver = QNDF()
            else
                if sam.set.solver != "FBDF"
                    @warn "Unavailable solver for SymbolicAWEModel: $(sam.set.solver). Falling back to FBDF."
                end
                solver = FBDF()
            end
        end

        if isnothing(outputs)
            @variables begin
                heading(t)[1:1]
                angle_of_attack(t)[1:1]
                tether_len(t)[1:3]
                winch_force(t)[1:3]
            end
            outputs = Num[]
            if length(sam.sys_struct.wings) > 0
                push!(outputs, heading[1], angle_of_attack[1])
            end
            if length(sam.sys_struct.winches) > 0
                push!(outputs, tether_len[1], winch_force[1])
            end
        end

        model_name = get_model_name(sam.set, sam.sys_struct)
        model_path = joinpath(KiteUtils.get_data_path(), model_name)
        prn && @info "Model bin name: $model_name"
        loaded = load_serialized_model!(sam, model_path; remake, reload)
        changed = false
        if !loaded
            sam.inputs = create_sys!(sam, sam.sys_struct; prn)
            changed = true
        end
        outputs_changed = isnothing(sam.outputs) ||
                            length(outputs) != length(sam.outputs) ||
                            !all(string.(outputs) .== string.(sam.outputs))
        if outputs_changed
            sam.lin_prob = nothing
            sam.control_functions = nothing
        end
        sam.outputs = outputs
        
        changed |= outputs_changed
        changed |= maybe_create_prob!(sam; create_prob, prn)
        changed |= maybe_create_lin_prob!(sam, outputs; create_lin_prob, prn)
        changed |= maybe_create_control_functions!(sam, outputs;
            create_control_func, prn)

        # Sync deserialized prob's flat parameters to the current sys_struct.
        if !isnothing(sam.prob)
            sync_params!(sam.prob.param_sync, sam.prob.prob, sam.sys_struct)
        end

        if changed
            prn && @info "Serializing model to: \n\t$model_path"
            serialize(model_path, sam.serialized_model)
        end

        if reinit_sys
            reinit!(sam.sys_struct, sam.set;
                    ignore_l0, remake_vsm, reset_vel,
                    apply_tether_lengths, prn)
        end
        # reinit! below syncs the struct's ICs onto the problem; no rebuild needed.
        if create_prob && !isnothing(sam.prob)
            prob = something(sam.prob)
            reset_integrator |= reload
            reinit!(sam, prob, solver;
                adaptive, reset_integrator, lin_vsm, vsm_min_wind)
        end
    end
    prn && @info "$(sam.sys_struct.name) model initialized in $time seconds."
    return sam.integrator
end


"""
    reinit!(sam, prob, solver; kwargs...) -> (ODEIntegrator, Bool)

Reset the ODE integrator from new initial conditions without rebuilding the
symbolic model. See [`init!`](@ref) for `adaptive`, `reset_integrator`,
`lin_vsm`, and `vsm_min_wind`.
"""
function reinit!(
    sam::SymbolicAWEModel,
    prob::ProbWithAttributes,
    solver;
    adaptive=true,
    reset_integrator=true,
    lin_vsm=true,
    vsm_min_wind=0.5
)
    dt = SimFloat(1/sam.set.sample_freq)
    existing = sam.integrator
    fresh = isnothing(existing) || !successful_retcode(existing.sol) || reset_integrator
    seed_set_values!(target) = isnothing(prob.set_set_values) ? nothing :
        prob.set_set_values(target,
            SimFloat[winch.set_value for winch in sam.sys_struct.winches])
    if fresh
        # A trailing reinit! would re-solve the DAE init and discard this synced state.
        sync_params!(prob.param_sync, prob.prob, sam.sys_struct)
        sync_initial!(prob.initial_sync, prob.prob, sam.sys_struct)
        seed_set_values!(prob.prob)
        integrator = init(prob.prob, solver;
            adaptive, dt, tspan=(0.0, dt), abstol=sam.set.abs_tol, reltol=sam.set.rel_tol,
            save_on=false, save_everystep=false)
        sam.integrator = integrator
    else
        integrator = existing
        sam.integrator = integrator
        sync_params!(prob.param_sync, integrator, sam.sys_struct)
        seed_set_values!(integrator)
        OrdinaryDiffEqCore.reinit!(integrator; reinit_dae=true)
    end
    update_sys_struct!(prob, integrator, sam.sys_struct)
    if lin_vsm && has_vsm_wing(sam.sys_struct)
        refresh_aero!(sam; vsm_min_wind)
        sync_params!(prob.param_sync, integrator, sam.sys_struct)
    end
    validate_sys_struct(sam.sys_struct)  # Check for division-by-zero issues
    return integrator, true
end

"""
    get_set_hash(set::Settings; fields)

Calculates a SHA1 hash for structural fields in the `Settings` object.
This is used to check if a cached compiled model is still valid.

# Structural Fields (affect symbolic equations):
- `:segments`: Number of tether segments (affects state vector size)
- `:model`: Kite model name (affects geometry)
- `:foil_file`: Airfoil data file (affects VSM setup)
- `:physical_model`: Model type (ram, simple_ram, 4_attach_ram)
- `:winch_model`: Winch dynamics model (affects winch equations)

# Runtime Fields (don't affect compilation, excluded from hash):
- `:profile_law`: Wind profile law (evaluated at runtime via symbolic function)
- `:wind_vec`, `:elevation`: Initial conditions
- Other runtime parameters
"""
function get_set_hash(set::Settings;
        fields=[:segments, :model, :foil_file, :physical_model, :winch_model]
    )
    hash_acc = zeros(UInt8, 1)
    for field in fields
        value = getfield(set, field)
        hash_acc = sha1(string((value, hash_acc)))
    end
    return hash_acc
end

"""
    get_sys_struct_hash(sys_struct::SystemStructure)

Calculates a SHA1 hash for the topology and structure of a `SystemStructure`.
This is used to check if a cached compiled model is still valid.

Includes all structural properties that affect the symbolic equations:
- Point connectivity and types (STATIC, DYNAMIC, WING, BODY_STATIC)
- Segment connectivity
- TwistSurface structure and types (STATIC, DYNAMIC)
- Pulley constraints and types
- Tether topology
- Winch configuration
- Wing topology, connectivity, aerodynamic model type (RIGID_DYNAMICS vs PARTICLE_DYNAMICS), and aero mode
- Transform hierarchy

Excludes runtime-configurable properties like masses, lengths, stiffnesses.
"""
function get_sys_struct_hash(sys_struct::SystemStructure)
    (; points, twist_surfaces, segments, pulleys, tethers, winches, wings, transforms,
       bodies, elastic_joints, timoshenko_joints) = sys_struct
    data_parts = []
    for point in points
        push!(data_parts, ("point", point.idx, point.wing_idx, point.body_idx, Int(point.type)))
    end
    for segment in segments
        # Stiffness type selects the spring law (scalar k·Δ vs callable F(ε)).
        stiff_type = segment.unit_stiffness isa Real ? "float" :
                     string(typeof(segment.unit_stiffness))
        push!(data_parts, ("segment", segment.idx, segment.point_idxs, stiff_type))
    end
    for twist_surface in twist_surfaces
        push!(data_parts, ("twist_surface", twist_surface.idx, twist_surface.point_idxs, Int(twist_surface.type)))
    end
    for pulley in pulleys
        push!(data_parts, ("pulley", pulley.idx, pulley.segment_idxs, Int(pulley.type)))
    end
    for tether in tethers
        push!(data_parts, ("tether", tether.idx, tether.segment_idxs))
    end
    for winch in winches
        push!(data_parts, ("winch", winch.idx, winch.tether_idxs))
    end
    for wing in wings
        wing_data = ("wing", wing.idx, wing.twist_surface_idxs,
                     Int(wing.dynamics_type),
                     nameof(typeof(wing.aero)),
                     aero_hash_id(wing.aero))

        # Include wing reference points in hash
        ref_hash(ref) = (ref.ids, ref.weights)
        rp_hash(ref_points) = isnothing(ref_points) ? nothing :
            (ref_hash(ref_points[1]), ref_hash(ref_points[2]))
        origin_hash(origin) = isnothing(origin) ? nothing :
            ref_hash(origin)
        wing_data = (wing_data...,
            rp_hash(wing.z_ref_points),
            rp_hash(wing.y_ref_points),
            origin_hash(wing.origin))
        push!(data_parts, wing_data)
    end
    for transform in transforms
        push!(data_parts, ("transform", transform.idx, transform.wing_idx, transform.rot_point_idx,
                        transform.base_point_idx, transform.base_transform_idx))
    end
    for rigid_body in bodies
        push!(data_parts, ("rigid_body", rigid_body.idx, Int(rigid_body.type)))
    end
    for joint in elastic_joints
        # Stiffness type selects the generated law (scalar `k·Δ` vs callable `k(Δ)`).
        stiff_type(s) = s isa Real ? "float" : string(typeof(s))
        push!(data_parts, ("elastic_joint", joint.idx,
                           joint.body_a_idx, joint.body_b_idx,
                           stiff_type(joint.stiffness_axial),
                           stiff_type(joint.stiffness_shear),
                           stiff_type(joint.stiffness_torsion),
                           stiff_type(joint.stiffness_bending)))
    end
    for joint in timoshenko_joints
        # Rigidity type selects the generated law (scalar vs callable of strain).
        rigidity_type(r) = r isa Real ? "float" : string(typeof(r))
        push!(data_parts, ("timoshenko_joint", joint.idx,
                           joint.body_a_idx, joint.body_b_idx,
                           rigidity_type(joint.EA), rigidity_type(joint.GA),
                           rigidity_type(joint.GJ), rigidity_type(joint.EIy),
                           rigidity_type(joint.EIz)))
    end
    content = string(data_parts)
    return sha1(content)
end
