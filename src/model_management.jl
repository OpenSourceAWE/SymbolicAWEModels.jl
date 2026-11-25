# Copyright (c) 2025 Bart van de Lint and Uwe Fechner
# SPDX-License-Identifier: LGPL-3.0-only

"""
    point_aero_force_array(sam, sys) -> Matrix{Num} or nothing

The world-frame aerodynamic force of every point of `sam` as one `3 x n_points` array to
fetch from the compiled system `sys`: a particle wing's node holds its body-frame
`aero_force_point_b` turned into the world, and any other point a literal zero, a
constant of the generated function. `nothing` when no wing has per-node forces.
"""
function point_aero_force_array(sam, sys)
    hasproperty(sys, :aero_force_point_b) || return nothing
    function particle_node(i)
        w = point_wing(sam, i)
        return !isnothing(w) && wing_node(sam.definition, i) && has_aero(sam.aero[w]) &&
            body_dynamics_type(sam, wing_body(sam, w)) == PARTICLE_DYNAMICS
    end
    points = eachindex(sam.definition.points)
    any(particle_node, points) || return nothing
    force = [particle_node(i) ?
                 collect(sys.body_R_b_to_w[:, :, wing_body(sam, point_wing(sam, i))] *
                         sys.aero_force_point_b[:, i]) :
                 zeros(Num, 3)
             for i in points]
    return reduce(hcat, force)
end

"""The `SysState` of run `sam`, read in a generated getter's scatter."""
@inline run_columns(sam) = getfield(sam, :sysstate)

"""The position slot before the first body's in the `SysState` of run `sam`."""
@inline body_slot0(sam) = getfield(sam, :topology).body_slot0

"""The number of wings of run `sam`, whose orientation frames come before the bodies'."""
@inline wing_count(sam) = getfield(sam, :topology).n_wings

"""
    generate_prob_getters(sam, sys)

Generate the functions that read and write the state and parameters of the compiled
`ODESystem` `sys`, whose topology the run `sam` defines. Returns them as a `NamedTuple`.
"""
function generate_prob_getters(sam, sys, param_registry=nothing,
                               initial_registry=nothing)
    definition = sam.definition
    (; points, wings, pulleys, winches, tethers, segments, bodies) = definition
    stations = definition.spanwise_stations
    get_aero_input, set_set_values, get_set_values = nothing, nothing, nothing

    specs = NamedTuple[]
    if length(points) > 0
        pairs = Pair[
            sys.pos           => (st, v) -> write_columns!(
                xyz(run_columns(st), :X, :Y, :Z), v),
            sys.vel           => (st, v) -> write_columns!(
                xyz(run_columns(st), :VX, :VY, :VZ), v),
            sys.point_force   => (st, v) -> write_columns!(
                xyz(run_columns(st), :force_x, :force_y, :force_z), v),
            sys.wind_at_point => (st, v) -> write_columns!(
                xyz(run_columns(st), :wind_x, :wind_y, :wind_z), v),
            sys.total_drag    => (st, v) -> write_columns!(xyz(run_columns(st),
                :drag_force_x, :drag_force_y, :drag_force_z), v)]
        aero_force = point_aero_force_array(sam, sys)
        isnothing(aero_force) || push!(pairs,
            aero_force => (st, v) -> write_columns!(
                xyz(run_columns(st), :aero_force_x, :aero_force_y, :aero_force_z), v))
        push!(specs, scatter_spec(pairs...))
    end
    if length(pulleys) > 0
        push!(specs, scatter_spec(
            sys.pulley_len => (st, v) -> write_column!(run_columns(st).pulley_len, v),
            sys.pulley_vel => (st, v) -> write_column!(run_columns(st).pulley_vel, v)))
    end
    if length(segments) > 0
        push!(specs, scatter_spec(
            sys.spring_force => (st, v) -> write_column!(run_columns(st).spring_force, v)))
    end
    if length(stations) > 0
        push!(specs, scatter_spec(
            sys.twist_angle => (st, v) -> write_column!(run_columns(st).twist_angles, v),
            sys.twist_ω     => (st, v) -> write_column!(run_columns(st).twist_vel, v)))
    end
    if length(winches) > 0
        push!(specs, scatter_spec(
            sys.winch_acc      => (st, v) -> write_column!(run_columns(st).winch_acc, v),
            sys.winch_vel      => (st, v) -> write_column!(run_columns(st).v_reelout, v),
            sys.set_values     => (st, v) -> write_column!(run_columns(st).set_torque, v),
            sys.winch_friction => (st, v) -> write_column!(
                run_columns(st).winch_friction, v)))
        set_set_values = setp(sys, sys.set_values)
        get_set_values = getp(sys, sys.set_values)
    end
    if length(tethers) > 0
        push!(specs, scatter_spec(
            sys.tether_len => (st, v) -> write_column!(run_columns(st).l_tether, v)))
    end
    if length(wings) > 0
        push!(specs, scatter_spec(
            sys.wind_vel_wing => (st, v) -> write_wing_body_columns!(
                xyz(run_columns(st), :wind_x, :wind_y, :wind_z), v, st),
            sys.aero_force_b  => (st, v) -> write_columns!(xyz(run_columns(st),
                :aero_force_KA_x, :aero_force_KA_y, :aero_force_KA_z), v),
            sys.aero_moment_b => (st, v) -> write_columns!(xyz(run_columns(st),
                :aero_moment_KA_x, :aero_moment_KA_y, :aero_moment_KA_z), v),
            sys.elevation     => (st, v) -> write_column!(run_columns(st).elevations, v),
            sys.azimuth       => (st, v) -> write_column!(run_columns(st).azimuths, v),
            sys.heading       => (st, v) -> write_column!(run_columns(st).headings, v),
            sys.turn_acc      => (st, v) -> write_columns!(xyz(run_columns(st),
                :turn_acc_x, :turn_acc_y, :turn_acc_z), v),
            sys.course        => (st, v) -> write_column!(run_columns(st).courses, v),
            sys.angle_of_attack => (st, v) -> write_column!(run_columns(st).AoAs, v)))

        # aero_input exists only for wings whose component exposes the connector.
        aero_inputs = [
            getproperty(sys, Symbol("aero_$w")).aero_input
            for w in eachindex(wings)
            if body_dynamics_type(sam, wing_body(sam, w)) === RIGID_DYNAMICS &&
                hasproperty(getproperty(sys, Symbol("aero_$w")), :aero_input)]
        get_aero_input = isempty(aero_inputs) ? nothing :
            getu(sys, collect.(aero_inputs))
    end
    if length(bodies) > 0
        push!(specs, scatter_spec(
            sys.body_Q_b_to_w => (st, v) -> write_columns!(
                xyz(run_columns(st), :Qw, :Qx, :Qy, :Qz), v, wing_count(st)),
            sys.body_ω_b      => (st, v) -> write_columns!(xyz(run_columns(st),
                :turn_rate_x, :turn_rate_y, :turn_rate_z), v, wing_count(st)),
            sys.body_pos_w    => (st, v) -> write_columns!(
                xyz(run_columns(st), :X, :Y, :Z), v, body_slot0(st)),
            sys.body_vel_w    => (st, v) -> write_columns!(
                xyz(run_columns(st), :VX, :VY, :VZ), v, body_slot0(st)),
            sys.body_acc_w    => (st, v) -> write_columns!(
                xyz(run_columns(st), :acc_x, :acc_y, :acc_z), v, body_slot0(st))))
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
    scatter_spec(pairs...)

Describe one component group: each `sys_array => copyfn` pair maps a symbolic output
array to a `(state, view) -> _` closure that writes that array into the run state's
columns. This is the single source of truth — both the buffer layout and the scatter
derive from the same ordered list.
"""
function scatter_spec(pairs::Pair...)
    return (; arrays  = [pair.first for pair in pairs],
            copyfns = Tuple(pair.second for pair in pairs))
end

"""The columns `names` of `sysstate`, as a tuple."""
@inline xyz(sysstate, a::Symbol, b::Symbol, c::Symbol) =
    (getfield(sysstate, a), getfield(sysstate, b), getfield(sysstate, c))
@inline xyz(sysstate, a::Symbol, b::Symbol, c::Symbol, d::Symbol) =
    (getfield(sysstate, a), getfield(sysstate, b), getfield(sysstate, c),
     getfield(sysstate, d))

"""
    write_columns!(columns, v, offset=0)

Write each column `i` of the matrix `v` into entry `offset + i` of `columns`, one
vector per row of `v`.
"""
@inline function write_columns!(columns::Tuple, v, offset=0)
    @inbounds for i in axes(v, 2), k in eachindex(columns)
        columns[k][offset + i] = v[k, i]
    end
    return nothing
end

"""
    write_wing_body_columns!(columns, v, sam)

Write each column `w` of the matrix `v`, one per wing of `sam`, into the position slot of
that wing's body in `columns`, one vector per row of `v`.
"""
@inline function write_wing_body_columns!(columns::Tuple, v, sam)
    topology = getfield(sam, :topology)
    @inbounds for w in axes(v, 2), k in eachindex(columns)
        columns[k][topology.body_slot0 + topology.wing_bodies[w]] = v[k, w]
    end
    return nothing
end

"""
    write_column!(column, v, offset=0)

Write each entry `i` of the vector `v` into entry `offset + i` of `column`.
"""
@inline function write_column!(column, v, offset=0)
    @inbounds for i in eachindex(v)
        column[offset + i] = v[i]
    end
    return nothing
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
    groups = Tuple(ScatterGroup(spec.copyfns, views)
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
        serialize(s, group.copyfns)
        serialize(s, map(size, group.views))
    end
end

function Serialization.deserialize(s::Serialization.AbstractSerializer,
                                   ::Type{<:InplaceGetter})
    fn = deserialize(s)
    n = deserialize(s)
    n_groups = deserialize(s)
    copyfns = Vector{Any}(undef, n_groups)
    shapes = Vector{Any}(undef, n_groups)
    for i in 1:n_groups
        copyfns[i] = deserialize(s)
        shapes[i] = deserialize(s)
    end
    buf = Vector{SimFloat}(undef, n)
    grouped_views = build_grouped_views(buf, Tuple(shapes))
    groups = Tuple(ScatterGroup(copyfns[i], grouped_views[i]) for i in 1:n_groups)
    return InplaceGetter(fn, buf, groups)
end

"""
    generate_lin_getters(sys)

Generate the parameter setters of the linearized system `sys`. Returns a `NamedTuple` with
`set_set_values`, the setter of the winch set-points.
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

Generate in-place and out-of-place control functions of the `ODESystem` `model` for the
variables `inputs` and `outputs`. Returns a `NamedTuple` of the functions (`f_oop`, `f_ip`,
`h_oop`, `h_ip`), the dimensions (`nu`, `nx`, `ny`) and the symbolic variables (`dvs`,
`psym`, `io_sys`).
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
    in_memory(serialized_model) -> Bool

Whether the model already carries a build, so a valid bin on disk need not be read
back. The [`MonolithBackend`](@ref) leaves its `full_sys` behind and the
[`KernelBackend`](@ref) only its assembled problem, so either one counts.
"""
function in_memory(serialized_model::SerializedModel)
    return !isnothing(serialized_model.full_sys) || !isnothing(serialized_model.prob)
end

"""
    load_serialized_model!(sam, model_path; remake=false, reload=false)

Load the model serialized at `model_path` into `sam.serialized_model` if its settings and
structure hashes match those of `sam`; returns `true` if a valid model is loaded, else
`false`. `remake` rejects the file and triggers a rebuild, `reload` reads it even with a
model in memory, and `prn` reports whether it was kept in memory, read back or rejected.
"""
function load_serialized_model!(sam, model_path; remake=false, reload=false,
                                prn=true)
    set_hash = get_set_hash(sam.set)
    structure_hash = SymbolicAWEModels.structure_hash(sam)

    if ispath(model_path) && !remake
        if set_hash != sam.serialized_model.set_hash ||
                structure_hash != sam.serialized_model.structure_hash
            prn && @info "The settings or the structure changed since this model " *
                         "was built; rebuilding."
            sam.serialized_model = SerializedModel(; set_hash, structure_hash)
            return false
        end
        if in_memory(sam.serialized_model) && !reload
            prn && @info "Model already in memory; not reading the bin."
            return true
        end
        try
            read_time = @elapsed serialized_model = deserialize(model_path)
            if set_hash == serialized_model.set_hash &&
                    structure_hash == serialized_model.structure_hash
                size = @sprintf("%.2f s, %.1f MB", read_time,
                                filesize(model_path) / 2^20)
                prn && @info "Read the model back ($size) from: \n\t$model_path"
                sam.serialized_model = serialized_model
                return true
            end
            prn && @info "The bin was built for other settings or another " *
                         "structure; rebuilding."
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
    sam.serialized_model = SerializedModel(; set_hash, structure_hash)
    return false
end

"""
    maybe_create_prob!(sam; create_prob=true, sparse=false,
                       analytic_jacobian=false, prn=true)

Compile the full system and create the `ODEProblem` with its getters and setters, if they
do not exist yet and `create_prob` is set. `sparse` and `analytic_jacobian` are those of
[`init!`](@ref), `prn` prints progress. Returns `true` if a new problem was created.
"""
function maybe_create_prob!(sam; create_prob=true, sparse=false,
                            analytic_jacobian=false, prn=true)
    if create_prob && isnothing(sam.prob)
        return build_prob!(sam.backend, sam; sparse, analytic_jacobian, prn)
    end
    return false
end

"""
    build_prob!(backend, sam; sparse=false, analytic_jacobian=false, prn=true)

Assemble `sam.prob`: the [`MonolithBackend`](@ref) `mtkcompile`s the flattened `full_sys`
into one `ODEProblem`, the [`KernelBackend`](@ref) schedules one kernel per component type.
`sparse` gives the solver a Jacobian sparsity pattern and `analytic_jacobian` a Jacobian
(MTK's symbolic one on the monolith). Returns `true` when a problem was built.
"""
function build_prob!(::MonolithBackend, sam; sparse=false, analytic_jacobian=false,
                     prn=true)
    isnothing(sam.full_sys) && return false
    analytic_jacobian && prn && @warn "`analytic_jacobian` on the MonolithBackend " *
        "is MTK's symbolic Jacobian, which is expensive to build and throws on its " *
        "first call for any model reaching a registered numerical leaf (the wind " *
        "profile, an aerodynamic polar), because those have no symbolic derivative."
    full_sys = something(sam.full_sys)
    local sys
    time = @elapsed @suppress_err sys = mtkcompile(full_sys; inputs=sam.inputs)
    prn && println("\tSimplified the System for ODEProblem in $time seconds.")

    dt = SimFloat(1/sam.set.sample_freq)
    # skip MTK's DAE init system; reinit! already sets consistent ICs
    fill_defaults = missing_param_defaults(sys, sam.defaults)
    time = @elapsed prob = ODEProblem{true, SciMLBase.FullSpecialize}(sys,
        [sam.defaults; fill_defaults], (0.0, dt);
        build_initializeprob=false, sparse, jac=analytic_jacobian)
    prn && println("\tCreated the ODEProblem in $time seconds.")

    time = @elapsed getters = generate_prob_getters(sam, sys,
        sam.param_registry, sam.initial_registry)
    prn && println("\tCreated state getters and setters in $time seconds.")

    sam.prob = ProbWithAttributes(; prob, getters...)
    return true
end

function build_prob!(backend::ModelBackend, sam; prn=true)
    throw(BackendUnsupportedError("build_prob!", backend))
end

"""
    init_backend!(backend, sam, solver; kwargs...)

Full `init!` path for a non-[`MonolithBackend`](@ref). The monolith uses the
`init!` body directly; other backends (currently [`KernelBackend`](@ref))
implement their own assembly + integrator build here and return the fresh
`ODEIntegrator`: assemble the problem from `sam`'s definition and run, and store
`sam.prob`/`sam.integrator`.
"""
function init_backend!(backend::ModelBackend, sam, solver; kwargs...)
    throw(BackendUnsupportedError("init!", backend))
end

"""
    maybe_create_lin_prob!(sam, outputs; ...)

Create and cache the `LinearizationProblem` for the variables `outputs`, if it does not
exist yet and `create_lin_prob` is set. `prn` prints progress. Returns `true` if a new
problem was created.
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

Create and cache the control functions for the variables `outputs`, if they do not exist
yet and `create_control_func` is set. `prn` prints progress. Returns `true` if new
functions were created.
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
    has_custom_component(sam)

Return `true` when the run `sam` has a non-default winch model or a wing using a
custom aero model, in which case the compiled model cannot be reused from cache
and must be rebuilt.
"""
has_custom_component(sam) =
    any(!is_builtin_winch(winch.model.dynamics) for winch in sam.definition.winches) ||
    any(!is_builtin_aero(mode) for mode in sam.aero)

"""
    place!(sam::SymbolicAWEModel; reset_vel=true, kwargs...)

Place the model's structure by its transforms with KiteGeometry's
`place!(sam.definition, sam.transforms; kwargs...)` and start its run there
([`start_run!`](@ref)). `reset_vel=false` keeps the velocities the run holds.
"""
function KiteGeometry.place!(sam::SymbolicAWEModel; reset_vel::Bool=true, kwargs...)
    place!(sam.definition, sam.transforms; kwargs...)
    start_run!(sam; reset_vel)
    return nothing
end

"""
    init!(sam::SymbolicAWEModel; ...)

Load or build the symbolic model, create the `ODEProblem` (and optionally the
`LinearizationProblem` / control functions), serialize new builds to disk, and
return a freshly initialized `ODEIntegrator` from the run `sam` as it stands. To
start from the design geometry again, [`place!`](@ref) the structure first.

# Keyword Arguments
- `solver`, `adaptive`: ODE solver and time-stepping mode. `solver=nothing` picks
  a default from `sam.set.solver`.
- `autodiff`: automatic-differentiation choice for the default solver's Jacobian.
  Defaults to [`default_autodiff`](@ref) of the backend. Ignored when `solver` is
  passed explicitly.
- `linsolve`: factorization the default solver uses. Defaults to
  [`default_linsolve`](@ref) of the backend. Ignored when `solver` is passed
  explicitly.
- `prn`: print progress messages.
- `remake`: force a full rebuild, ignoring any cached compiled model. `nothing`
  rebuilds when a custom winch/aero component is present (see
  [`has_custom_component`](@ref)) and reuses the cache otherwise.
- `reload`: force reloading the serialized model from disk.
- `outputs`: vector of output variables (used by linearization / control funcs).
- `create_prob`, `create_lin_prob`, `create_control_func`: which artefacts to build.
- `lin_vsm`: linearize the VSM aerodynamics after init.
- `remake_vsm`: rebuild the VSM wing/aero from settings (after editing
  `aero_geometry.yaml` etc.).
- `reset_integrator`: discard the existing integrator and build a fresh one, so no
  stale BDF history taints the next run.
- `vsm_min_wind=0.5`: minimum |va| [m/s] for the initial VSM solve. Below this the
  solve is skipped and the wing's aero outputs are zeroed.
- `sparse`: give the solver a Jacobian sparsity pattern, so the
  finite-difference Jacobian is coloured and its factorization sparse. The
  [`MonolithBackend`](@ref) takes MTK's structural pattern, the
  [`KernelBackend`](@ref) derives its own from the wiring
  ([`state_sparsity`](@ref)). Part of the serialized model's name, so sparse and
  dense builds are cached separately. `nothing` takes the backend's
  [`default_sparse`](@ref).
- `analytic_jacobian=nothing`: give the solver a Jacobian instead of letting it
  differentiate the right-hand side numerically. `nothing` takes the backend's
  [`default_analytic_jacobian`](@ref). Part of the serialized model's name, as
  `sparse` is.
"""
function init!(sam::SymbolicAWEModel;
    solver=nothing, autodiff=default_autodiff(sam.backend), adaptive=true, prn=true,
    remake=nothing, reload=false,
    outputs=nothing,
    create_prob::Bool=true,
    create_lin_prob::Bool=false,
    create_control_func::Bool=false,
    lin_vsm::Bool=true,
    remake_vsm::Bool=false,
    reset_integrator::Bool=true,
    vsm_min_wind=0.5,
    sparse::Union{Bool, Nothing}=nothing,
    linsolve=default_linsolve(sam.backend),
    analytic_jacobian::Union{Bool, Nothing}=nothing
)
    sparse = something(sparse, default_sparse(sam.backend))
    analytic_jacobian = something(analytic_jacobian,
                                  default_analytic_jacobian(sam.backend))
    prn && @info "Initializing $(sam.definition.metadata.name) model..."
    isnothing(sam.set) && error("init!: $(sam.definition.metadata.name) has no settings")
    if isnothing(remake)
        remake = has_custom_component(sam)
        remake && prn && @info "Custom winch/aero model detected; " *
            "forcing remake (custom component equations are not " *
            "captured by the model hash)."
    end
    time = @elapsed begin
        if isnothing(solver)
            if sam.set.solver == "QNDF"
                @warn "This solver is not tested."
                solver = QNDF(; autodiff, linsolve)
            else
                if sam.set.solver != "FBDF"
                    @warn "Unavailable solver for SymbolicAWEModel: $(sam.set.solver). Falling back to FBDF."
                end
                solver = FBDF(; autodiff, linsolve)
            end
        end

        if !(sam.backend isa MonolithBackend)
            integrator = init_backend!(sam.backend, sam, solver;
                adaptive, prn, remake_vsm, reset_integrator, vsm_min_wind,
                lin_vsm, sparse, analytic_jacobian, remake, reload)
            prn && @info "$(sam.definition.metadata.name) model initialized " *
                "($(nameof(typeof(sam.backend))))."
            return integrator
        end

        if isnothing(outputs)
            @variables begin
                heading(t)[1:1]
                angle_of_attack(t)[1:1]
                tether_len(t)[1:3]
                winch_force(t)[1:3]
            end
            outputs = Num[]
            if length(sam.definition.wings) > 0
                push!(outputs, heading[1], angle_of_attack[1])
            end
            if length(sam.definition.winches) > 0
                push!(outputs, tether_len[1], winch_force[1])
            end
        end

        model_name = get_model_name(sam; sparse, analytic_jacobian)
        model_path = joinpath(KiteUtils.get_data_path(), model_name)
        prn && @info "Model bin name: $model_name"
        loaded = load_serialized_model!(sam, model_path; remake, reload, prn)
        changed = false
        if !loaded
            sam.inputs = create_sys!(sam; prn)
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
        changed |= maybe_create_prob!(sam; create_prob, sparse, analytic_jacobian, prn)
        changed |= maybe_create_lin_prob!(sam, outputs; create_lin_prob, prn)
        changed |= maybe_create_control_functions!(sam, outputs;
            create_control_func, prn)

        if !isnothing(sam.prob)
            sync_run_params!(sam.prob.param_sync, sam.prob.prob, sam)
        end

        if changed
            prn && @info "Serializing model to: \n\t$model_path"
            serialize(model_path, sam.serialized_model)
        end

        remake_vsm && remake_wing_aero!(sam)
        init_wind!(sam)
        # reinit!(sam, prob) syncs the struct's ICs onto the problem; no rebuild needed.
        if create_prob && !isnothing(sam.prob)
            prob = something(sam.prob)
            reset_integrator |= reload
            reinit!(sam, prob, solver;
                adaptive, reset_integrator, lin_vsm, vsm_min_wind, prn)
        end
    end
    prn && @info "$(sam.definition.metadata.name) model initialized in $time seconds."
    return sam.integrator
end


"""
    reinit!(sam, integrator; solver=integrator.alg, kwargs...) -> (ODEIntegrator, Bool)

Reset `integrator`, which is `sam`'s own, to the initial state of `sam`'s definition.
`solver` defaults to the one `integrator` was built with, so the reset stays on the
compiled right-hand side instead of forcing a second compilation at another
Jacobian's element type. `kwargs` are those of
[`reinit!`](@ref)`(sam, prob, solver; …)`.
"""
function reinit!(sam::SymbolicAWEModel,
        integrator::OrdinaryDiffEqCore.ODEIntegrator;
        solver = integrator.alg, kwargs...)
    isnothing(sam.prob) && error("reinit!: $(sam.definition.metadata.name) has no " *
        "ODEProblem; call init! first.")
    return reinit!(sam, sam.prob, solver; kwargs...)
end

"""
    reinit!(sam, prob, solver; kwargs...) -> (ODEIntegrator, Bool)

Reset the ODE integrator from new initial conditions without rebuilding the
symbolic model. See [`init!`](@ref) for `adaptive`, `reset_integrator`,
`lin_vsm`, and `vsm_min_wind`. `prn` toggles the per-phase timing logs
(integrator build / first-call JIT, initial aero solve).
"""
function reinit!(
    sam::SymbolicAWEModel,
    prob::ProbWithAttributes,
    solver;
    adaptive=true,
    reset_integrator=true,
    lin_vsm=true,
    vsm_min_wind=0.5,
    prn=true
)
    dt = SimFloat(1/sam.set.sample_freq)
    existing = sam.integrator
    fresh = isnothing(existing) || !successful_retcode(existing.sol) || reset_integrator
    seed_set_values!(target) = isnothing(prob.set_set_values) ? nothing :
        prob.set_set_values(target, copy(sam.sysstate.set_torque))
    if fresh
        # A fresh run starts from the definition, not from an earlier integrator.
        sam.integrator = nothing
        # A trailing reinit! would re-solve the DAE init and discard this synced state.
        sync_run_params!(prob.param_sync, prob.prob, sam)
        sync_initial!(prob.initial_sync, prob.prob, sam)
        seed_set_values!(prob.prob)
        prn && @info "JIT-compiling the RHS and building the integrator (first " *
            "solve; this is the slow step on large beams)…"
        jit_time = @elapsed integrator = init(prob.prob, solver;
            adaptive, dt, tspan=(0.0, dt), abstol=sam.set.abs_tol, reltol=sam.set.rel_tol,
            save_on=false, save_everystep=false)
        prn && @info "RHS compiled & integrator built in $(round(jit_time; digits=1)) s."
        sam.integrator = integrator
    else
        integrator = existing
        sync_aero_params!(prob.param_sync, integrator, sam)
        seed_set_values!(integrator)
        OrdinaryDiffEqCore.reinit!(integrator; reinit_dae=true)
    end
    update_run_state!(prob, integrator, sam)
    if lin_vsm && has_vsm_wing(sam)
        prn && @info "Running initial VSM aero solve…"
        aero_time = @elapsed refresh_aero!(sam; vsm_min_wind, cold_start=true)
        prn && @info "Initial aero solved in $(round(aero_time; digits=1)) s."
        sync_aero_params!(prob.param_sync, integrator, sam)
    end
    validate_run(sam)
    return integrator, true
end

"""
    get_set_hash(set::Settings; fields)

SHA1 hash of the `fields` of `set` that shape the compiled equations, by default
`:segments`, `:model`, `:foil_file`, `:physical_model` and `:winch_model`, to check a cached
compiled model. Numeric settings and initial conditions (`:g_earth`, `:wind_vec`,
`:cd_tether`, `:v_wind`, `:profile_law`, `:elevation`) are read at runtime and not hashed.
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
    structure_hash(sam)

SHA1 hash of the topology of the run `sam`, to check a cached compiled model: the
connectivity and types of its points (with a `BODY_STATIC` point's beam joint), segments,
stations, pulleys, tethers, winches, bodies and tubes, each wing's dynamics type, aero mode
and transform base point, and [`PerPointWind`](@ref); no mass, length or stiffness value.
"""
function structure_hash(sam)
    sys = sam.definition
    data_parts = []
    for (idx, point) in enumerate(sys.points)
        # nothing-vs-set gates whether drag/damping equations are emitted (structure,
        # not value), so 0.0↔1.0 reuses the cached bin but adding/removing regenerates.
        push!(data_parts, ("point", idx, something(point_wing(sam, idx), 0),
                           something(point.body, 0), something(point.tube, 0),
                           Int(point_type(point))))
    end
    for (idx, segment) in enumerate(sys.segments)
        # Stiffness type selects the spring law (scalar k·Δ vs callable F(ε)).
        stiffness = run_value(sam, :segments, idx, :unit_stiffness)
        stiff_type = stiffness isa Real ? "float" : string(typeof(stiffness))
        push!(data_parts, ("segment", idx, segment.points, stiff_type))
    end
    for (idx, station) in enumerate(sys.spanwise_stations)
        push!(data_parts, ("station", idx, station.points,
                           Int(station_type(station)),
                           something(station.flap_bodies, Int[]),
                           something(station.flap_points, Int[])))
    end
    for (idx, pulley) in enumerate(sys.pulleys)
        push!(data_parts, ("pulley", idx, pulley.segments,
                           Int(pulley_type(pulley))))
    end
    for (idx, tether) in enumerate(sys.tethers)
        push!(data_parts, ("tether", idx, tether.segments))
    end
    for (idx, winch) in enumerate(sys.winches)
        push!(data_parts, ("winch", idx, winch.tethers))
    end
    for (w, mode) in enumerate(sam.aero)
        b = wing_body(sam, w)
        # Polar format sets the panel interpolant dimensionality baked into the bin.
        polar_format = mode isa AbstractVSMAero && !isempty(mode.vsm_aero.panels) ?
            mode.vsm_aero.panels[1].aero_model : nothing
        # flow_curvature adds a moment term and the lag adds states, so both are
        # structure; their constants are parameters, as is the apparent mass.
        frame = wing_frame(sam, w)
        frame_refs = isnothing(frame) ? nothing :
            (frame.origin, frame.z_ref_points, frame.y_ref_points)
        push!(data_parts, ("wing", w, b, wing_station_idxs(sam, w),
                           Int(body_dynamics_type(sam, b)), nameof(typeof(mode)),
                           aero_hash_id(mode), polar_format,
                           flow_curvature_enabled(mode), wagner_enabled(mode),
                           frame_refs, wing_origin(sam, w),
                           transform_base_point(sam.transforms, b)))
    end
    per_point_wind(sam) && push!(data_parts, ("wind_mode", :per_point))
    for (idx, body) in enumerate(sys.bodies)
        push!(data_parts, ("rigid_body", idx, Int(body_type(sam, idx)),
                           something(body.wing, 0)))
    end
    for (idx, tube) in enumerate(sys.tubes)
        columns = tube_columns(tube)
        # Rigidity type selects the generated law (scalar vs callable of strain).
        rigidity_types = [value isa Real ? "float" : string(typeof(value))
                          for value in (run_value(sam, :tubes, idx, field)
                                        for field in rigidity_fields(columns))]
        push!(data_parts, ("tube", idx, tube.bodies[1], tube.bodies[2],
                           nameof(typeof(columns)), rigidity_types...))
    end
    return sha1(string(data_parts))
end
