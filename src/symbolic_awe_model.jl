# Copyright (c) 2025 Bart van de Lint and Uwe Fechner
# SPDX-License-Identifier: LGPL-3.0-only

const LinType = @NamedTuple{A::Matrix{SimFloat}, B::Matrix{SimFloat}, C::Matrix{SimFloat}, D::Matrix{SimFloat}}

"""
    ScatterGroup{Fns, Views}

One component group of an [`InplaceGetter`](@ref). `copyfns` is a tuple of
`(state, view) -> _` closures, one per output array, each writing that array into the
run state's columns; `views` is the matching tuple of zero-copy reshaped views into
the getter's buffer.
"""
struct ScatterGroup{Fns, Views}
    copyfns::Fns
    views::Views
end

"""
    InplaceGetter{F, B, G}

A single zero-allocation getter that both reads and scatters all per-step
component state. `fn` is an in-place MTK observed function `fn(buf, u, p, t)`
over the concatenation of every component's output arrays; `buf` is a
preallocated flat buffer reused each call; `groups` is a tuple of
[`ScatterGroup`](@ref)s that write the freshly-computed buffer straight into the
run state's columns. One spec drives both the buffer layout and the scatter.
"""
struct InplaceGetter{F, B, G}
    fn::F
    buf::B
    groups::G
end

"Apply each `(state, view)` copy closure of one group (tuple recursion)."
@inline scatter_arrays(::Tuple{}, ::Tuple{}, state) = nothing
@inline function scatter_arrays(copyfns::Tuple, views::Tuple, state)
    first(copyfns)(state, first(views))
    scatter_arrays(Base.tail(copyfns), Base.tail(views), state)
    return nothing
end

"Scatter every group into the run state `state` (tuple recursion over the groups)."
@inline scatter_groups(::Tuple{}, state) = nothing
@inline function scatter_groups(groups::Tuple, state)
    group = first(groups)
    scatter_arrays(group.copyfns, group.views, state)
    scatter_groups(Base.tail(groups), state)
    return nothing
end

"""
    (g::InplaceGetter)(integ, sam)

Evaluate all component state at the integrator's current point and scatter it
into the run `sam` in place. Reads `integ.u`/`integ.p`/`integ.t` as direct fields;
this method is itself the function barrier that keeps the call allocation-free.
"""
function (g::InplaceGetter)(integ, sam)
    g.fn(g.buf, integ.u, integ.p, integ.t)
    scatter_groups(g.groups, sam)
    return nothing
end

"""
    @with_kw struct ProbWithAttributes{...}

A container for the main Ordinary Differential Equation (ODE) problem and its
associated getter and setter functions for the full, nonlinear physical state.
"""
@with_kw struct ProbWithAttributes{Prob, SetSetValues,
                                  GetSetValues, GetAeroInput,
                                  GetAllState, ParamSync, InitialSync}
    "The ODE problem for the full nonlinear model."
    prob::Prob

    # Setters for the ODE
    "Syncs the run's parameters into the flat buffer once per step."
    param_sync::ParamSync
    "Pushes the run's initial conditions onto the problem's `Initial` params."
    initial_sync::InitialSync
    "Setter for the control input values."
    set_set_values::SetSetValues

    # Getters for the ODE state
    get_set_values::GetSetValues
    get_aero_input::GetAeroInput
    "One monolithic zero-alloc getter for all per-step component state."
    get_all_state::GetAllState
end

"""
    Base.getproperty(pa::ProbWithAttributes, sym::Symbol)

Overloads `getproperty` to provide convenient access to the simplified system
(`sys`) contained within the ODE problem's function definition.
"""
function Base.getproperty(pa::ProbWithAttributes, sym::Symbol)
    if sym == :sys
        # Access the `prob` field of the `pa` struct to get to its contents.
        prob = getfield(pa, :prob)
        return prob.f.sys
    end
    return getfield(pa, sym)
end

"""
    @with_kw struct LinProbWithAttributes{SetLinSetValues, SetLinSys, SetLinSet, LinOut}

A container for the general-purpose linearization problem and the resulting full
linearized model (A,B,C,D matrices).

$(TYPEDFIELDS)
"""
@with_kw struct LinProbWithAttributes{Prob, SetSetValues}
    "Linearization problem of the mtk model."
    prob::Prob

    # Setters for the linearization
    set_set_values::SetSetValues
end

"""
    @with_kw struct ControlFuncWithAttributes{FIP, FOOP, HIP, HOOP, DVS, PSYM}

A container for callable control functions and their symbolic representations,
generated from the full system model.

$(TYPEDFIELDS)
"""
@with_kw struct ControlFuncWithAttributes{FIP, FOOP, HIP, HOOP, DVS, PSYM}
    "In-place dynamics function f(dx, x, u, p, t)."
    f_ip::FIP
    "Out-of-place dynamics function dx = f(x, u, p, t)."
    f_oop::FOOP
    "In-place observation function h(y, x, u, p, t)."
    h_ip::HIP
    "Out-of-place observation function y = h(x, u, p, t)."
    h_oop::HOOP
    "Number of inputs (u)."
    nu::Int
    "Number of states (x)."
    nx::Int
    "Number of outputs (y)."
    ny::Int
    "The symbolic state vector."
    dvs::DVS
    "The symbolic parameter vector."
    psym::PSYM
    "The generated input-output system."
    io_sys::ModelingToolkit.System
end

"""
    @with_kw mutable struct SerializedModel{...}

A type-stable container for the compiled and serialized components of a
`SymbolicAWEModel`: the products of the `ModelingToolkit.jl` compilation, grouped
into nested attribute structs (`ProbWithAttributes`, etc.).

$(TYPEDFIELDS)
"""
@with_kw mutable struct SerializedModel{D<:AbstractVector}
    set_hash::Vector{UInt8}
    "Hash of the topology the model was compiled for ([`structure_hash`](@ref))."
    structure_hash::Vector{UInt8}
    "Unsimplified system of the mtk model"
    full_sys::Union{ModelingToolkit.System, Nothing} = nothing
    defaults::D = Pair{Num, Any}[]
    "Symbolic representation of the control inputs."
    inputs::Union{Symbolics.Arr, Vector{Num}} = Num[]
    "Outputs of the linearization and control function."
    outputs::Union{Symbolics.Arr, Vector{Num}} = Num[]

    "Container for the ODE problem and its getters/setters."
    prob::Union{ProbWithAttributes, Nothing} = nothing
    "Container for the linearization problem and its components."
    lin_prob::Union{LinProbWithAttributes, Nothing} = nothing
    "Container for the control functions."
    control_functions::Union{ControlFuncWithAttributes, Nothing} = nothing
end

"""
    MonolithParameter(parameters)

A parameter of the [`MonolithBackend`](@ref)'s integrator: every compiled parameter
standing for it, all written alike.
"""
struct MonolithParameter
    parameters::Vector{Any}
end

"""
    KernelParameter(slots, elements)

A parameter of the [`KernelBackend`](@ref)'s integrator: the `slots` of its numeric
buffer standing for it, and the element of the value each holds, `0` for all of it.
"""
struct KernelParameter
    slots::Vector{Int}
    elements::Vector{Int}
end

"""Access to a parameter of a run's integrator, `nothing` for one it does not hold."""
const ParameterAccess = Union{Nothing, MonolithParameter, KernelParameter}

"""
    mutable struct SymbolicAWEModel <: AbstractKiteModel

A run of a kite power system: its structure, a KiteGeometry `SystemDefinition`, the
settings and aero modes it is flown with, its state and its integrator, built using
`ModelingToolkit.jl`. Interact with it through [`init!`](@ref), [`next_step!`](@ref),
[`run_value`](@ref) and [`set_run_value!`](@ref) rather than through its fields.

$(TYPEDFIELDS)
"""
@with_kw mutable struct SymbolicAWEModel{SS<:SystemDefinition, SM<:SerializedModel,
                                         ST<:SysState, TP<:NamedTuple} <: AbstractKiteModel
    "The structure the run simulates"
    definition::SS
    "The run's settings; `nothing` for a model that only draws or places its structure"
    set::Union{Nothing, Settings}
    "The KiteGeometry transforms placing the run"
    transforms::Rows{Transform}
    "The VSM settings the wings' engines are built from"
    vsm_set::Union{Nothing, VortexStepMethod.VSMSettings}
    "The atmosphere of `set`"
    am::Union{Nothing, AtmosphericModel}
    "Where the wind comes from: [`ProfileWind`](@ref) or [`PerPointWind`](@ref)"
    wind_mode::WindMode
    "The aero mode of each wing, with its VSM engine where it has one"
    aero::Vector{AbstractAeroModel}
    "The run's state, a KiteUtils `SysState`"
    sysstate::ST
    "The lookups of the structure the run reads its items by ([`run_topology`](@ref))"
    topology::TP
    "Container for the compiled and serialized model components"
    serialized_model::SM
    "The ODE integrator for the full nonlinear model"
    integrator::Union{OrdinaryDiffEqCore.ODEIntegrator, Nothing} = nothing
    "Access to the integrator's parameters by the path each is read from"
    parameter_access::Dict{Tuple{Symbol, Int, Symbol}, ParameterAccess} =
        Dict{Tuple{Symbol, Int, Symbol}, ParameterAccess}()
    "Relative start time of the current time interval"
    t_0::SimFloat = 0.0
    "Number of next_step! calls"
    iter::Int64 = 0
    "Time spent in the VSM linearization step"
    t_vsm::SimFloat  = zero(SimFloat)
    "Time spent in the ODE integration step"
    t_step::SimFloat = zero(SimFloat)
    "Build-time flattened-parameter registry (transient, never serialized)."
    param_registry::Any = nothing
    "Build-time initial-condition registry (transient, never serialized)."
    initial_registry::Any = nothing
    "Assembly/feature backend; defaults to [`MonolithBackend`](@ref)."
    backend::ModelBackend = default_backend()
end

"""
    SAM_FIELDS

The fields of `SymbolicAWEModel`, as opposed to the properties `getproperty` and
`setproperty!` take from its `serialized_model`.
"""
const SAM_FIELDS = fieldnames(SymbolicAWEModel)

"""
    Base.getproperty(sam::SymbolicAWEModel, sym::Symbol)

A field of `sam`, else the property `sym` of its `serialized_model`: its compiled
functions and other model components.
"""
function Base.getproperty(sam::SymbolicAWEModel, sym::Symbol)
    sym in SAM_FIELDS && return getfield(sam, sym)
    return getproperty(getfield(sam, :serialized_model), sym)
end

"""
    Base.setproperty!(sam::SymbolicAWEModel, sym::Symbol, val)

Set a field of `sam`, else the property `sym` of its `serialized_model`. A new
`integrator` drops the access to the parameters of the previous one, and a new
`serialized_model` or `prob` drops the integrator built from the previous one.
"""
function Base.setproperty!(sam::SymbolicAWEModel, sym::Symbol, val)
    sym in (:integrator, :serialized_model, :prob) &&
        empty!(getfield(sam, :parameter_access))
    sym in (:serialized_model, :prob) && setfield!(sam, :integrator, nothing)
    sym in SAM_FIELDS || return setproperty!(getfield(sam, :serialized_model), sym, val)
    return setfield!(sam, sym, val)
end

"""
    SymbolicAWEModel(set, definition::SystemDefinition; transforms=Transform[],
                     vsm_set=nothing, aero=Dict(), wind_mode=ProfileWind(), kwargs...)

A run of `definition`, placed by the KiteGeometry `transforms`, with the settings `set`,
at its start from the pose `definition` is placed in ([`start_run!`](@ref)). Building it
compiles nothing; [`init!`](@ref) does.

# Arguments
- `set`: The run's settings; `nothing` for a model that only draws or places its
  structure.
- `definition::SystemDefinition`: The physical system definition.
- `transforms`: The KiteGeometry transforms placing `definition`, as KiteGeometry's
  `load_authoring` returns them; their rates start the run moving.
- `vsm_set`: The VSM settings the model's engines are built from; the project's
  `vsm_settings` by default.
- `aero`: The aero mode of each wing by its name, in place of the one its row names
  ([`wing_aero_spec`](@ref)): a flat plate's own polars, or a mode's options.
- `wind_mode::WindMode`: [`PerPointWind`](@ref) makes every point's `wind_vec` a
  settable parameter instead of a height-profile output.
- `kwargs...`: Further fields of the `SymbolicAWEModel`, such as its `backend`.
"""
function SymbolicAWEModel(set::Union{Nothing, Settings}, definition::SystemDefinition;
                          transforms=Transform[], vsm_set=nothing,
                          aero=Dict{Symbol, AbstractAeroModel}(),
                          wind_mode::WindMode=ProfileWind(), kwargs...)
    isnothing(set) || (set.physical_model = definition.metadata.name)
    specs = [wing_aero_spec(wing, aero) for wing in definition.wings]
    vsm_set = model_vsm_set(specs, set, vsm_set)
    modes = AbstractAeroModel[model_aero(spec, wing,
                                         length(station_indices(definition, w)), set,
                                         vsm_set)
                              for (w, (spec, wing)) in enumerate(zip(specs,
                                                                     definition.wings))]
    n_corners = sum(n_aero_log_points(mode, definition, w)
                    for (w, mode) in enumerate(modes); init=0)
    P, counts = sysstate_counts(definition, n_corners, sum(n_vsm_panels, modes; init=0))
    serialized_model = SerializedModel(;
        set_hash = isnothing(set) ? UInt8[] : get_set_hash(set), structure_hash = UInt8[])
    sam = SymbolicAWEModel(; definition, set,
        transforms = convert(Rows{Transform}, transforms), vsm_set,
        am = isnothing(set) ? nothing : AtmosphericModel(set), wind_mode, aero = modes,
        sysstate = SysState(P; counts..., precision=SimFloat),
        topology = run_topology(definition, n_corners), serialized_model, kwargs...)
    start_lengths!(sam)
    foreach(w -> setup_model_aero!(sam, w), eachindex(modes))
    start_run!(sam)
    sam.serialized_model.structure_hash = structure_hash(sam)
    return sam
end

"""`vsm_set`, else the project's `vsm_settings` for a run of `set` with wings flown with
the aero mode `specs` of which one is VSM's."""
function model_vsm_set(specs, set, vsm_set)
    isnothing(vsm_set) || return vsm_set
    (isnothing(set) || !any(spec -> spec isa AeroSpec, specs)) && return nothing
    path = optional_project_file("vsm_settings")
    return isnothing(path) ? nothing : VortexStepMethod.VSMSettings(path; data_prefix=false)
end

"""
    sysstate_counts(sys, n_corners, panels) -> (P, counts)
    sysstate_counts(sam) -> (P, counts)

The position slot count `P` and the keywords `counts` of `KiteUtils.SysState(P; ...)`
for a run of definition `sys` whose wings log `n_corners` aero points and hold `panels`
VSM panels, shared by its `sysstate` and the log frames written from it: its
`gamma_distribution` holds the circulation of every VSM panel of its wings, in order.
"""
function sysstate_counts(sys, n_corners, panels)
    counts = (wings = length(sys.wings), orients = n_orient_frames(sys),
              deflections = n_flap_deflections(sys), pulleys = length(sys.pulleys),
              winches = length(sys.winches), tethers = length(sys.tethers),
              segments = length(sys.segments), panels)
    return position_slots(sys, n_corners).total, counts
end
sysstate_counts(sam::SymbolicAWEModel) = sysstate_counts(sam.definition,
    length(sam.topology.slots.panel_corners), sum(n_vsm_panels, sam.aero; init=0))

"""Copy the columns a run writes directly from its state `from` into the log frame `to`."""
function copy_run_columns!(to::SysState, from::SysState)
    copyto!(to.force_x, from.force_x)
    copyto!(to.force_y, from.force_y)
    copyto!(to.force_z, from.force_z)
    copyto!(to.wind_x, from.wind_x)
    copyto!(to.wind_y, from.wind_y)
    copyto!(to.wind_z, from.wind_z)
    copyto!(to.aero_force_x, from.aero_force_x)
    copyto!(to.aero_force_y, from.aero_force_y)
    copyto!(to.aero_force_z, from.aero_force_z)
    copyto!(to.acc_x, from.acc_x)
    copyto!(to.acc_y, from.acc_y)
    copyto!(to.acc_z, from.acc_z)
    copyto!(to.turn_acc_x, from.turn_acc_x)
    copyto!(to.turn_acc_y, from.turn_acc_y)
    copyto!(to.turn_acc_z, from.turn_acc_z)
    copyto!(to.winch_acc, from.winch_acc)
    copyto!(to.winch_friction, from.winch_friction)
    copyto!(to.elevations, from.elevations)
    copyto!(to.azimuths, from.azimuths)
    copyto!(to.headings, from.headings)
    copyto!(to.courses, from.courses)
    copyto!(to.AoAs, from.AoAs)
    copyto!(to.gamma_distribution, from.gamma_distribution)
    return to
end

"""
    update_sys_state!(ss::SysState, sam::SymbolicAWEModel, zoom=1.0)

Update a `SysState` from the model's integrator: converts units (e.g. radians to
degrees) and computes derived values like AoA. Roll, pitch and yaw are
`KiteUtils.euler_KS(ss.orient)`. Each body's pose, a wing's among them, is logged at
its body slot and frame.

# Arguments
- `ss::SysState`: The state struct to be updated.
- `sam::SymbolicAWEModel`: The source model.
- `zoom::SimFloat=1.0`: A scaling factor for the position coordinates.
"""
function update_sys_state!(ss::SysState, sam::SymbolicAWEModel, zoom=1.0)
    copy_run_columns!(ss, sam.sysstate)
    ss.time = isnothing(sam.integrator) ? 0.0 : sam.integrator.t
    sys = sam.definition
    value(collection, i, quantity) = run_value(sam, Val(collection), i, Val(quantity))
    copyto!(ss.l_tether, sam.sysstate.l_tether)
    for (i, winch) in enumerate(sys.winches)
        isempty(winch.tethers) && continue
        ss.v_reelout[i] = value(:winches, i, :vel)
        ss.winch_force[i] = norm(value(:winches, i, :force))
        ss.set_torque[i] = value(:winches, i, :set_value)
    end
    copyto!(ss.twist_angles, sam.sysstate.twist_angles)
    copyto!(ss.twist_vel, sam.sysstate.twist_vel)
    if length(sys.wings) > 0
        b = wing_body(sam, 1)
        ss.acc = norm(value(:bodies, b, :acc_w))
        ss.orient .= value(:bodies, b, :Q_b_to_w)
        ss.turn_rates .= value(:wings, 1, :turn_rate)
        ss.elevation = value(:wings, 1, :elevation)
        ss.azimuth = value(:wings, 1, :azimuth)
        ss.heading = value(:wings, 1, :heading)
        ss.course = value(:wings, 1, :course)
        ss.v_app = norm(value(:bodies, b, :va_b))
        ss.v_wind_kite .= value(:bodies, b, :wind_vec)
        if ss.v_app > 1e-6
            ss.AoA = calc_aoa(sam.aero[1], sam, 1)
            ss.side_slip = calc_side_slip(sam, 1)
        else
            ss.AoA = NaN
            ss.side_slip = NaN
        end
        ss.aero_force_KA .= value(:wings, 1, :aero_force_b)
        ss.aero_moment_KA .= value(:wings, 1, :aero_moment_b)
        ss.vel_kite .= value(:bodies, b, :vel_w)
    end
    run = sam.sysstate
    for i in eachindex(sys.points)
        ss.X[i] = run.X[i] * zoom
        ss.Y[i] = run.Y[i] * zoom
        ss.Z[i] = run.Z[i] * zoom
        ss.VX[i], ss.VY[i], ss.VZ[i] = run.VX[i], run.VY[i], run.VZ[i]
        ss.drag_force_x[i] = run.drag_force_x[i]
        ss.drag_force_y[i] = run.drag_force_y[i]
        ss.drag_force_z[i] = run.drag_force_z[i]
    end
    copyto!(ss.spring_force, run.spring_force)

    corner_idx = length(sys.points)
    for (w, mode) in enumerate(sam.aero)
        corner_idx = write_aero_log_points!(mode, sam, w, ss, corner_idx, zoom)
    end
    write_flap_deflections!(ss, sam)

    for b in eachindex(sys.bodies)
        write_body_state!(ss, sam, b, zoom)
    end
    copyto!(ss.pulley_len, run.pulley_len)
    copyto!(ss.pulley_vel, run.pulley_vel)

    ss.v_wind_gnd .= sam.set.wind_vec
    nothing
end

"""
    write_body_state!(ss, sam, b, zoom)

Log the pose and rates of body `b` of `sam` into `ss`: position and velocity in
`X/Y/Z` and `VX/VY/VZ` at its position slot, orientation and turn rate in `orients` and
`turn_rate_x/y/z` at its frame.
"""
function write_body_state!(ss, sam, b, zoom)
    slot, frame = body_position_slot(sam, b), body_frame(sam, b)
    pos_w = run_value(sam, Val(:bodies), b, Val(:pos_w))
    ss.X[slot] = pos_w[1] * zoom
    ss.Y[slot] = pos_w[2] * zoom
    ss.Z[slot] = pos_w[3] * zoom
    ss.VX[slot], ss.VY[slot], ss.VZ[slot] = run_value(sam, Val(:bodies), b, Val(:vel_w))
    ss.orients[frame] .= run_value(sam, Val(:bodies), b, Val(:Q_b_to_w))
    ω_b = run_value(sam, Val(:bodies), b, Val(:ω_b))
    ss.turn_rate_x[frame] = ω_b[1]
    ss.turn_rate_y[frame] = ω_b[2]
    ss.turn_rate_z[frame] = ω_b[3]
    return nothing
end

"""
    position_slots(sys, n_corners) -> NamedTuple

Index layout of a `SysState`'s `X/Y/Z` position arrays for a run of definition `sys`
whose wings log `n_corners` aero points: structural `points`, then the aero log
`panel_corners`, then `bodies`, a wing's body logging the wing's pose. Each field is the
`UnitRange` of slots for that group (empty if none); `total` is the position count.
Body `b` uses orientation frame (`orients`) `b`. A run holds its layout in
`sam.topology.slots`.
"""
function position_slots(sys, n_corners)
    n_points = length(sys.points)
    n_bodies = length(sys.bodies)
    base_bodies = n_points + n_corners
    return (points        = 1:n_points,
            panel_corners = (n_points + 1):base_bodies,
            bodies        = (base_bodies + 1):(base_bodies + n_bodies),
            total         = base_bodies + n_bodies)
end

"""Number of orientation frames, one per body, at least 1."""
n_orient_frames(sys) = max(1, length(sys.bodies))

"""
    SysState(s::SymbolicAWEModel, zoom=1.0)

Construct a `SysState` holding the current state of the model.

# Arguments
- `s::SymbolicAWEModel`: The source model.
- `zoom::SimFloat=1.0`: A scaling factor for the position coordinates.

# Returns
- `SysState`: A new state struct representing the current model state.
"""
function SysState(s::SymbolicAWEModel, zoom=1.0; precision=KiteUtils.MyFloat)
    P, counts = sysstate_counts(s)
    ss = SysState(P; counts..., precision)
    update_sys_state!(ss, s, zoom)
    ss
end

"""
    KiteUtils.Logger(sam::SymbolicAWEModel, steps::Int)

Constructs a `Logger` from a `SymbolicAWEModel` with the correct number of points.

This convenience constructor automatically calculates the total number of points
including VSM panel corners (4 corners per panel) and creates a Logger with the
appropriate size.

# Arguments
- `sam::SymbolicAWEModel`: The AWE model to create a logger for.
- `steps::Int`: The number of time steps to allocate for logging.

# Returns
- `Logger`: A new logger with size for all points including panel corners.

# Example
```julia
logger = Logger(sam, 1000)
```
"""
function KiteUtils.Logger(sam::SymbolicAWEModel, steps::Int;
                          precision=KiteUtils.MyFloat)
    P, counts = sysstate_counts(sam)
    return Logger(P, steps; counts..., precision)
end

"""
    next_step!(s::SymbolicAWEModel, integrator::ODEIntegrator; set_values, dt, vsm_interval)

Take a simulation step, using the provided integrator.
"""
function next_step!(
    s::SymbolicAWEModel,
    integrator::OrdinaryDiffEqCore.ODEIntegrator;
    set_values=nothing, dt=1/s.set.sample_freq,
    vsm_interval=1, vsm_min_wind=0.5, vsm_warn_on_fail=false
)
    !(s.integrator === integrator) && error(
        "The ODEIntegrator doesn't belong to " *
        "the SymbolicAWEModel")
    next_step!(s; set_values, dt, vsm_interval, vsm_min_wind, vsm_warn_on_fail)
end

"""
    next_step!(s::SymbolicAWEModel; set_values, dt,
               vsm_interval, vsm_min_wind, vsm_warn_on_fail)

Advance the simulation by one time step, optionally
updating control inputs and re-linearizing the VSM
model, then update `s.sysstate` from the ODE
integrator. Errors on an unstable solver retcode.

# Keyword Arguments
- `set_values=nothing`: Control input values.
    If `nothing`, current values are used.
- `dt=1/s.set.sample_freq`: Time step size [s].
- `vsm_interval=1`: Steps between VSM
    re-linearization. 0 disables re-linearization.
- `vsm_min_wind=0.5`: Minimum apparent wind [m/s] for
    a VSM solve. Below this the solve is skipped and
    the wing's aero outputs are zeroed.
- `vsm_warn_on_fail=false`: Warn instead of erroring
    when a VSM solve fails, keeping the circulation, the
    angles of attack and the frozen forces of the last
    converged solve. The next scheduled update solves
    again.
"""
function next_step!(sam::SymbolicAWEModel;
    set_values=nothing, dt=1/sam.set.sample_freq,
    vsm_interval=1, vsm_min_wind=0.5, vsm_warn_on_fail=false
)
    prob = sam.prob
    integrator = sam.integrator
    if isnothing(integrator)
        error("next_step! called before init!: integrator is not initialized")
    end
    if isnothing(set_values)
        set_values = copy(sam.sysstate.set_torque)
    end
    if prob isa ProbWithAttributes && !isnothing(prob.set_set_values)
        prob.set_set_values(integrator, set_values)
    end
    if prob isa ProbWithAttributes
        sync_params!(prob.param_sync, integrator, sam)
    end

    sam.t_0 = integrator.t
    # Setting parameters flags a derivative discontinuity; FBDF restarts on it.
    OrdinaryDiffEqCore.derivative_discontinuity!(integrator, false)
    sam.t_step = @elapsed OrdinaryDiffEqCore.step!(integrator, dt, true)
    if !successful_retcode(integrator.sol)
        throw(AssertionError("Solver unstable at t=" *
            "$(round(integrator.t; digits=4))" *
            ": $(integrator.sol.retcode)"))
    end
    sam.iter += 1
    if prob isa ProbWithAttributes
        update_run_state!(prob, integrator, sam)
        if vsm_interval != 0 && sam.iter % vsm_interval == 0 && has_vsm_wing(sam)
            sam.t_vsm = @elapsed begin
                refresh_aero!(sam; vsm_min_wind, vsm_warn_on_fail)
                sync_aero_params!(prob.param_sync, integrator, sam)
            end
        end
    end
    return nothing
end

"""
    update_run_state!(prob::ProbWithAttributes, integ::ODEIntegrator, sam)

Write the state of `integ` into the run `sam`'s `SysState`, through the generated
getter of `prob`.
"""
function update_run_state!(prob::ProbWithAttributes,
                           integ::OrdinaryDiffEqCore.ODEIntegrator, sam)
    prob.get_all_state(integ, sam)
    return nothing
end

"""
    get_model_name(sam::SymbolicAWEModel; precompile=false, sparse=false,
                   analytic_jacobian=false)

Constructs a unique filename for the serialized model based on its configuration.
The filename includes the SymbolicAWEModels version, Julia version, physical model,
wing type, dynamics type, and component counts to ensure that the correct cached
model is loaded. `sparse`, `analytic_jacobian` and the model's backend are part of
the name because the cached `ODEProblem` carries its Jacobian prototype, its Jacobian
and its backend's assembly, so those builds are different artefacts; naming them apart
keeps all of them on disk rather than invalidating one with the other.
"""
function get_model_name(sam::SymbolicAWEModel; precompile=false, sparse=false,
                        analytic_jacobian=false)
    sys = sam.definition
    suffix = precompile ? ".default" : ""
    pkg_ver = pkgversion(SymbolicAWEModels)
    ver = "$(VERSION.major).$(VERSION.minor)"

    dynamics_types = [wing_dynamics_type(wing) for wing in sys.wings]
    dynamics_type_str = if isempty(dynamics_types)
        "no_wing"
    elseif all(wt -> wt === RIGID_DYNAMICS, dynamics_types)
        "rigid"
    elseif all(wt -> wt === PARTICLE_DYNAMICS, dynamics_types)
        "particle"
    else
        "mixed"
    end

    aero_tags = unique(aero_mode_tag(mode) for mode in sam.aero)
    aero_mode_str = if isempty(aero_tags)
        ""
    elseif length(aero_tags) == 1
        only(aero_tags)
    else
        "mixed_aero_modes"
    end

    n_bodies = length(sys.bodies)
    body_tag = n_bodies > 0 ? "_$(n_bodies)bdy" : ""
    sparse_tag = sparse ? "_sparse" : ""
    jacobian_tag = analytic_jacobian ? "_analytic" : ""

    return "model_v$(pkg_ver)_jl$(ver)_$(sam.set.physical_model)_" *
        "$(dynamics_type_str)_$(aero_mode_str)_dynamic_$(length(sys.points))pnt_" *
        "$(length(sys.segments))seg_$(length(sys.spanwise_stations))grp_" *
        "$(length(sys.wings))wng_$(length(sys.winches))wch$(body_tag)$(sparse_tag)" *
        "$(jacobian_tag)$(backend_tag(sam.backend)).bin$suffix"
end

"""
    calc_steady_torque(sam::SymbolicAWEModel)

Calculates the torque for each winch that results in zero acceleration (steady state).
"""
calc_steady_torque(sam::SymbolicAWEModel) =
    [-something(winch.drum_radius, NaN) / something(winch.gear_ratio, NaN) *
     norm(run_value(sam, :winches, i, :force)) + run_value(sam, :winches, i, :friction)
     for (i, winch) in enumerate(sam.definition.winches)]

"""
    calc_heading(sam::SymbolicAWEModel)

Heading angles [rad] of the wings of `sam`, one per wing, by the tangential sphere frame
method.
"""
calc_heading(sam::SymbolicAWEModel) =
    [calc_heading(body_rotation(sam, b), run_value(sam, :bodies, b, :pos_w))
     for b in (wing_body(sam, w) for w in eachindex(sam.definition.wings))]

"""
    calc_aoa(sam::SymbolicAWEModel)

Angle of attack [rad] of the first wing, dispatched on its aero mode
([`calc_aoa(::AbstractAeroModel, sam, w)`](@ref)). `NaN` if the mode defines no AoA.
"""
calc_aoa(sam::SymbolicAWEModel) = calc_aoa(sam.aero[1], sam, 1)

"""
    unstretched_length(sam::SymbolicAWEModel)

Returns the unstretched tether length [m] for each tether.
"""
unstretched_length(sam::SymbolicAWEModel) = copy(sam.sysstate.l_tether)

"""
    tether_length(sam::SymbolicAWEModel)

Returns the current stretched tether length [m] for each tether.
"""
tether_length(sam::SymbolicAWEModel) =
    [tether_stretched_length(sam, i) for i in eachindex(sam.definition.tethers)]

"""
    calc_height(sam::SymbolicAWEModel)

Returns the height (z-position) [m] of the wing.
"""
calc_height(sam::SymbolicAWEModel) = run_value(sam, :wings, 1, :pos_w)[3]

"""
    winch_force(sam::SymbolicAWEModel)

Returns the winch force [N] for each winch.
"""
winch_force(sam::SymbolicAWEModel) =
    [norm(winch_force_vector(sam, i)) for i in eachindex(sam.definition.winches)]

"""
    spring_forces(sam::SymbolicAWEModel)

Returns the spring force [N] for each tether segment.
"""
spring_forces(sam::SymbolicAWEModel) = copy(sam.sysstate.spring_force)

"""
    pos(sam::SymbolicAWEModel)

Returns a vector of the position vectors [m] for each point in the system.
"""
pos(sam::SymbolicAWEModel) =
    [run_value(sam, :points, i, :pos_w) for i in eachindex(sam.definition.points)]
