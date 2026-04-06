# Copyright (c) 2025 Bart van de Lint and Uwe Fechner
# SPDX-License-Identifier: MPL-2.0

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
function generate_prob_getters(sys_struct, sys)
    c = collect
    @unpack wings, groups, pulleys, winches, tethers, segments = sys_struct
    get_wing_state, get_vsm_y, get_segment_state, get_group_state, get_pulley_state,
    get_winch_state, get_tether_state, set_set_values, get_set_values = ntuple(_ -> nothing, 9)

    if length(wings) > 0
        wing_vars = c.([
            sys.Q_b_to_w, sys.ω_b, sys.wing_pos,
            sys.wing_vel, sys.wing_acc,
            sys.va_wing_b, sys.wind_vel_wing,
            sys.aero_force_b, sys.aero_moment_b,
            sys.moment_tether_wing,
            sys.force_tether_wing,
            sys.elevation, sys.elevation_vel,
            sys.elevation_acc,
            sys.azimuth, sys.azimuth_vel,
            sys.azimuth_acc,
            sys.heading, sys.turn_rate,
            sys.turn_acc, sys.course,
            sys.angle_of_attack,
            # Principal frame state
            sys.com_w, sys.com_vel,
            sys.Q_p_to_w, sys.ω_p,
        ])
        get_wing_state = getu(sys, wing_vars)

        # vsm_input_state only exists for QUATERNION + AERO_LINEARIZED wings
        has_linearized = any(
            w.wing_type == QUATERNION &&
            w.aero_mode == AERO_LINEARIZED
            for w in sys_struct.wings)
        if has_linearized
            get_vsm_y = getu(sys, sys.vsm_input_state)
        else
            get_vsm_y = nothing
        end
    end
    if length(segments) > 0; get_segment_state = getu(sys, c.([sys.spring_force, sys.len])); end
    if length(groups) > 0; get_group_state = getu(sys, c.([sys.twist_angle, sys.twist_ω, sys.group_tether_force, sys.group_tether_moment, sys.group_aero_moment])); end
    if length(pulleys) > 0; get_pulley_state = getu(sys, c.([sys.pulley_len, sys.pulley_vel])); end
    if length(winches) > 0
        get_winch_state = getu(sys, c.([sys.tether_len, sys.tether_vel, sys.tether_acc,
                                       sys.set_values, sys.winch_force_vec, sys.tau_friction]))
        set_set_values = setp(sys, sys.set_values)
        get_set_values = getp(sys, sys.set_values)
    end
    if length(tethers) > 0; get_tether_state = getu(sys, c(sys.stretched_len)); end
    set_sys = setp(sys, sys.psys)
    set_set = setp(sys, sys.pset)
    get_struct_state = getu(sys, sys.wind_vec_gnd)

    # Always include va_point_b and point_mass in point_state (calculated for all points now)
    get_point_state = getu(sys, c.([sys.pos, sys.vel, sys.point_force, sys.va_point_b, sys.point_mass]))

    return (; get_wing_state, get_vsm_y, get_segment_state, get_group_state,
            get_pulley_state, get_winch_state, get_tether_state, set_set_values,
            get_set_values, set_sys, set_set, get_struct_state, get_point_state)
end

"""
    generate_lin_getters(sys)

Generate setter functions for the parameters of a linearized system.

# Arguments
- `sys`: The linearized ModelingToolkit system.

# Returns
- A `NamedTuple` containing setter functions for the winch set-points (`set_set_values`),
  the system structure parameters (`set_sys`), and the settings parameters (`set_set`).
"""
function generate_lin_getters(sys)
    set_set_values = nothing
    if hasproperty(sys, :set_values)
        set_set_values = setp(sys, sys.set_values)
    end
    set_sys = setp(sys, sys.psys)
    set_set = setp(sys, sys.pset)
    return (; set_set_values, set_sys, set_set)
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
    nu, nx, ny = length(inputs), length(dvs), length(outputs)
    (h_oop, h_ip) = ModelingToolkit.build_explicit_observed_function(
        io_sys, outputs; inputs, return_inplace=true
    )
    return (f_oop=f_oop, f_ip=f_ip, h_oop=h_oop, h_ip=h_ip, nu=nu, nx=nx, ny=ny,
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
        catch e
            @warn "Failure to deserialize $model_path: $(typeof(e))"
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
        local sys
        time = @elapsed @suppress_err sys = mtkcompile(sam.full_sys; inputs=sam.inputs)
        prn && println("\tSimplified the System for ODEProblem in $time seconds.")

        dt = SimFloat(1/sam.set.sample_freq)
        time = @elapsed prob = ODEProblem(sys, sam.defaults, (0.0, dt); u0_prior=sam.guesses)
        prn && println("\tCreated the ODEProblem in $time seconds.")

        time = @elapsed getters = generate_prob_getters(sam.sys_struct, sys)
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
        time = @elapsed @suppress_err begin
            lin_fun, lin_sys = linearization_function(sam.full_sys, [sam.inputs...], outputs;
                                                    op=sam.defaults, guesses=sam.guesses)
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
        inputs = [sam.inputs...]
        time = @elapsed result = generate_control_funcs(sam.full_sys, inputs, outputs)
        sam.control_functions = ControlFuncWithAttributes(; result...)
        prn && println("\tCreated the control functions in $time seconds.")
        return true
    end
    return false
end

"""
    init!(sam::SymbolicAWEModel; ...)

Initialize the `SymbolicAWEModel`.

This is the main entry point for setting up the model. It handles:
- Loading or building the symbolic model (`full_sys`).
- Creating the `ODEProblem`, `LinearizationProblem`, and control functions as needed.
- Serializing the model to disk if it was newly built.
- Initializing the ODE integrator.

# Keyword Arguments
- `solver`: The ODE solver to use. If `nothing`, a default is chosen based on settings.
- `adaptive::Bool`: Enable adaptive time-stepping for the solver.
- `prn::Bool`: Enable printing of progress messages.
- `remake::Bool`: Force a full rebuild of the symbolic model, ignoring any cached versions.
- `reload::Bool`: Force reloading of the serialized model from disk.
- `outputs`: A vector of variables to be treated as system outputs.
- `create_prob::Bool`: Whether to create the `ODEProblem`.
- `create_lin_prob::Bool`: Whether to create the `LinearizationProblem`.
- `create_control_func::Bool`: Whether to generate the control functions.
- `lin_vsm::Bool`: Whether to linearize the aerodynamics using the
                   Vortex Step Method (VSM) after initialization.
- `remake_vsm::Bool`: Recreate VSM wing and aerodynamics from settings (useful after
                      modifying aero_geometry.yaml or other VSM settings).

# Returns
- The initialized `ODEIntegrator`.
"""
function init!(sam::SymbolicAWEModel;
    solver=nothing, adaptive=true, prn=true,
    remake=false, reload=false,
    outputs=nothing,
    create_prob::Bool=true,
    create_lin_prob::Bool=false,
    create_control_func::Bool=false,
    lin_vsm::Bool=true,
    ignore_l0::Bool=false,
    remake_vsm::Bool=true,
    reset_vel::Bool=true,
    tunable_params::Bool=false
)
    prn && @info "Initializing $(sam.sys_struct.name) model..."
    sam.sys_struct isa SystemStructure{VSMWing} || error(
        "Equation generation requires SystemStructure{VSMWing}, " *
        "got SystemStructure{$(eltype(sam.sys_struct.wings))}.")
    time = @elapsed begin
        if isnothing(solver)
            solver = if sam.set.solver == "FBDF"
                sam.set.quasi_static ? FBDF(nlsolve=OrdinaryDiffEqNonlinearSolve.NLNewton(relax=sam.set.relaxation)) : FBDF()
            elseif sam.set.solver == "QNDF"
                @warn "This solver is not tested."
                QNDF()
            else
                @warn "Unavailable solver for SymbolicAWEModel: $(sam.set.solver). Falling back to FBDF."
                sam.set.quasi_static ? FBDF(nlsolve=OrdinaryDiffEqNonlinearSolve.NLNewton(relax=sam.set.relaxation)) : FBDF()
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
            sam.inputs = create_sys!(sam, sam.sys_struct; prn, tunable_params)
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
        changed |= maybe_create_control_functions!(sam, outputs; create_control_func, prn)

        if changed
            prn && @info "Serializing model to: \n\t$model_path"
            serialize(model_path, sam.serialized_model)
        end

        reinit!(sam.sys_struct, sam.set;
                ignore_l0, remake_vsm, reset_vel)
        # When reset_vel=false, state-dependent u0 changed;
        # force ODEProblem recreation to pick up new defaults.
        if !reset_vel && !isnothing(sam.prob)
            sam.prob = nothing
            changed = true
            changed |= maybe_create_prob!(sam;
                create_prob, prn)
        end
        create_prob && !isnothing(sam.prob) && reinit!(sam, sam.prob, solver; adaptive, reload, lin_vsm)
        create_lin_prob && !isnothing(sam.lin_prob) && reinit!(sam, sam.lin_prob)
    end
    prn && @info "$(sam.sys_struct.name) model initialized in $time seconds."
    return sam.integrator
end

"""
    reinit!(sam::SymbolicAWEModel, lin_prob::ModelingToolkit.LinearizationProblem)

Reinitializes a `LinearizationProblem` with the current system and settings parameters.

This function updates the internal parameter vectors of the linearization problem
with the latest values from the `SymbolicAWEModel`'s `sys_struct` and `set` fields.
"""
function reinit!(sam::SymbolicAWEModel, prob::LinProbWithAttributes)
    prob.set_sys(prob.prob, sam.sys_struct)
    prob.set_set(prob.prob, sam.set)
    nothing
end

"""
    reinit!(s::SymbolicAWEModel, prob::ODEProblem, solver; prn, precompile, reload, outputs) -> (ODEIntegrator, Bool)

Reinitializes an existing kite power system model's ODE integrator.

This function resets the integrator's state with new values from `s.set`,
allowing for the simulation to be restarted from a new initial condition
without needing to rebuild the entire symbolic model.

# Arguments
- `s::SymbolicAWEModel`: The kite power system state object.
- `prob::ODEProblem`: The ODE problem to be solved.
- `solver`: The solver to be used.

# Keyword Arguments
- `adaptive::Bool=true`: Whether to use adaptive time-stepping.
- `reload::Bool=true`: Force reloading the model from disk.
- `lin_vsm::Bool=true`: If `true`, linearizes the VSM model after reinitialization.

# Returns
- `(ODEIntegrator, Bool)`: A tuple containing the reinitialized integrator and a success flag.
"""
function reinit!(
    sam::SymbolicAWEModel,
    prob::ProbWithAttributes,
    solver;
    adaptive=true,
    reload=true, 
    lin_vsm=true
)
    dt = SimFloat(1/sam.set.sample_freq)
    existing = sam.integrator
    integrator = if isnothing(existing) || !successful_retcode(existing.sol) || reload
        init(prob.prob, solver;
            adaptive, dt, tspan=(0.0, dt), abstol=sam.set.abs_tol, reltol=sam.set.rel_tol,
            save_on=false, save_everystep=false)
    else
        something(existing)
    end
    sam.integrator = integrator
    prob.set_sys(integrator, sam.sys_struct)
    prob.set_set(integrator, sam.set)
    OrdinaryDiffEqCore.reinit!(integrator; reinit_dae=true)
    lin_vsm && update_vsm!(sam, prob)
    update_sys_struct!(prob, integrator, sam.sys_struct)
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
- `:quasi_static`: Whether points are quasi-static (affects equations)
- `:winch_model`: Winch dynamics model (affects winch equations)

# Runtime Fields (don't affect compilation, excluded from hash):
- `:profile_law`: Wind profile law (evaluated at runtime via symbolic function)
- `:v_wind`, `:elevation`: Initial conditions
- Other runtime parameters
"""
function get_set_hash(set::Settings;
        fields=[:segments, :model, :foil_file, :physical_model, :quasi_static, :winch_model]
    )
    h = zeros(UInt8, 1)
    for field in fields
        value = getfield(set, field)
        h = sha1(string((value, h)))
    end
    return h
end

"""
    get_sys_struct_hash(sys_struct::SystemStructure)

Calculates a SHA1 hash for the topology and structure of a `SystemStructure`.
This is used to check if a cached compiled model is still valid.

Includes all structural properties that affect the symbolic equations:
- Point connectivity and types (STATIC, DYNAMIC, QUASI_STATIC, WING)
- Segment connectivity
- Group structure and types (FIXED, TWIST)
- Pulley constraints and types
- Tether topology
- Winch configuration
- Wing topology, connectivity, aerodynamic model type (QUATERNION vs REFINE), and aero mode
- Transform hierarchy

Excludes runtime-configurable properties like masses, lengths, stiffnesses.
"""
function get_sys_struct_hash(sys_struct::SystemStructure)
    @unpack points, groups, segments, pulleys, tethers, winches, wings, transforms = sys_struct
    data_parts = []
    for point in points
        push!(data_parts, ("point", point.idx, point.wing_idx, Int(point.type)))
    end
    for segment in segments
        push!(data_parts, ("segment", segment.idx, segment.point_idxs))
    end
    for group in groups
        push!(data_parts, ("group", group.idx, group.point_idxs, Int(group.type)))
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
        wing_data = ("wing", wing.idx, wing.group_idxs,
                     Int(wing.wing_type),
                     Int(wing.aero_mode))

        # Include wing reference points in hash
        if wing isa VSMWing
            _ref_hash(r) = (r.ids, r.weights)
            _rp_hash(rp) = isnothing(rp) ? nothing :
                (_ref_hash(rp[1]), _ref_hash(rp[2]))
            wing_data = (wing_data...,
                _rp_hash(wing.z_ref_points),
                _rp_hash(wing.y_ref_points),
                wing.origin_idx)
        end

        push!(data_parts, wing_data)
    end
    for transform in transforms
        push!(data_parts, ("transform", transform.idx, transform.wing_idx, transform.rot_point_idx,
                        transform.base_point_idx, transform.base_transform_idx))
    end
    content = string(data_parts)
    return sha1(content)
end
