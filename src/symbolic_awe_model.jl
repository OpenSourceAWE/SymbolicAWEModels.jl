# Copyright (c) 2025 Bart van de Lint and Uwe Fechner
# SPDX-License-Identifier: LGPL-3.0-only

const LinType = @NamedTuple{A::Matrix{SimFloat}, B::Matrix{SimFloat}, C::Matrix{SimFloat}, D::Matrix{SimFloat}}
const GetSetNothing = Union{AbstractIndexer, Nothing}

"""
    ScatterGroup{Sel, Fns, Views}

One component group of an [`InplaceGetter`](@ref). `selector(sys_struct)` returns
the group's component vector (e.g. `sys_struct.points`); `copyfns` is a tuple of
`(component, view) -> _` closures, one per output array, each copying that
array's slice into the component's struct field; `views` is the matching tuple of
zero-copy reshaped views into the getter's buffer.
"""
struct ScatterGroup{Sel, Fns, Views}
    selector::Sel
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
`SystemStructure` fields. One spec drives both the buffer layout and the scatter.
"""
struct InplaceGetter{F, B, G}
    fn::F
    buf::B
    groups::G
end

"Copy column `idx` of the matrix view `v` into the mutable vector `field` in place."
@inline copy_vec!(field, v, idx) = (@views field .= v[:, idx]; nothing)

"Apply each `(component, view)` copy closure for one component (tuple recursion)."
@inline scatter_component(::Tuple{}, ::Tuple{}, component) = nothing
@inline function scatter_component(copyfns::Tuple, views::Tuple, component)
    first(copyfns)(component, first(views))
    scatter_component(Base.tail(copyfns), Base.tail(views), component)
    return nothing
end

"Scatter every group into `sys_struct` (tuple recursion over heterogeneous groups)."
@inline scatter_groups(::Tuple{}, sys_struct) = nothing
@inline function scatter_groups(groups::Tuple, sys_struct)
    group = first(groups)
    for component in group.selector(sys_struct)
        scatter_component(group.copyfns, group.views, component)
    end
    scatter_groups(Base.tail(groups), sys_struct)
    return nothing
end

"""
    (g::InplaceGetter)(integ, sys_struct)

Evaluate all component state at the integrator's current point and scatter it
into `sys_struct` in place. Reads `integ.u`/`integ.p`/`integ.t` as direct fields;
this method is itself the function barrier that keeps the call allocation-free.
"""
function (g::InplaceGetter)(integ, sys_struct)
    g.fn(g.buf, integ.u, integ.p, integ.t)
    scatter_groups(g.groups, sys_struct)
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
    "Syncs flattened struct-field parameters into the flat buffer once per step."
    param_sync::ParamSync
    "Pushes the struct's initial conditions onto the problem's `Initial` params."
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

A type-stable container for the compiled and serialized components of a `SymbolicAWEModel`.

This struct holds the products of the `ModelingToolkit.jl` compilation process,
now organized into nested attribute structs (`ProbWithAttributes`, etc.).
This simplifies the structure and improves serialization robustness.

$(TYPEDFIELDS)
"""
@with_kw mutable struct SerializedModel{D<:AbstractVector}
    set_hash::Vector{UInt8}
    sys_struct_hash::Vector{UInt8}
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
    mutable struct SymbolicAWEModel <: AbstractKiteModel

The main state container for a kite power system model, built using `ModelingToolkit.jl`.

This struct holds the complete state of the simulation, including the physical
structure (`SystemStructure`), the compiled model (`SerializedModel`), the atmospheric
model, and the ODE integrator.

Users typically interact with this model through high-level functions like
[`init!`](@ref) and [`next_step!`](@ref) rather than accessing its fields directly.

# Type Parameters
- `S`: Scalar type, typically `SimFloat`.
- `V`: Vector type, typically `KVec3`.
- `P`: Number of tether points in the system.

$(TYPEDFIELDS)
"""
@with_kw mutable struct SymbolicAWEModel{SS<:SystemStructure, SM<:SerializedModel} <: AbstractKiteModel
    "Reference to the point mass system with points, segments, pulleys and tethers"
    sys_struct::SS
    "Container for the compiled and serialized model components"
    serialized_model::SM # Now strongly typed
    "The ODE integrator for the full nonlinear model"
    integrator::Union{OrdinaryDiffEqCore.ODEIntegrator, Nothing} = nothing
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

Tuple of field names that are direct fields of `SymbolicAWEModel` (as opposed to fields
delegated to the nested `serialized_model`). Used by `getproperty` and `setproperty!`
to dispatch field access correctly.
"""
const SAM_FIELDS = (:sys_struct, :serialized_model, :integrator, :t_0, :iter, :t_vsm, :t_step, :param_registry, :initial_registry, :backend)

"""
    Base.getproperty(sam::SymbolicAWEModel, sym::Symbol)

Overloads `getproperty` to allow direct access to fields within the nested `serialized_model`.
This provides a convenient way to access compiled functions and other model
components without explicitly referencing `sam.serialized_model`.
"""

function Base.getproperty(sam::SymbolicAWEModel, sym::Symbol)
    if sym === :set
        getfield(sam, :sys_struct).set
    elseif sym === :am
        getfield(sam, :sys_struct).am
    elseif sym in SAM_FIELDS
        getfield(sam, sym)
    else
        getproperty(getfield(sam, :serialized_model), sym)
    end
end

"""
    Base.setproperty!(sam::SymbolicAWEModel, sym::Symbol, val)

Overloads `setproperty!` to allow direct setting of fields within the nested `serialized_model`.
This allows you to change properties of the compiled model as if they were
fields of the `SymbolicAWEModel` itself.
"""
function Base.setproperty!(sam::SymbolicAWEModel, sym::Symbol, val)
    if sym === :set
        error("Cannot replace `set`: it is owned by `sys_struct` " *
              "(const field). Mutate fields directly, " *
              "e.g. `sam.set.wind_vec = ...`.")
    elseif sym in SAM_FIELDS
        setfield!(sam, sym, val)
    else
        setproperty!(getfield(sam, :serialized_model), sym, val)
    end
end

"""
    SymbolicAWEModel(set::Settings, sys_struct::SystemStructure; kwargs...)

Constructs a `SymbolicAWEModel` from an existing `SystemStructure`.

This is the primary inner constructor. It takes a `SystemStructure` that defines the
physical layout of the kite system and prepares it for symbolic model generation.

# Arguments
- `set::Settings`: Configuration parameters.
- `sys_struct::SystemStructure`: The physical system definition.
- `kwargs...`: Further keyword arguments passed to the `SymbolicAWEModel` constructor.

# Returns
- `SymbolicAWEModel`: A model ready for symbolic equation generation via [`init!`](@ref).
"""
function SymbolicAWEModel(
    set::Settings,
    sys_struct::SystemStructure;
    kwargs...
)
    @assert set === sys_struct.set "The `set` argument must be the" *
        " same object as `sys_struct.set`"
    set_hash = get_set_hash(set)
    sys_struct_hash = get_sys_struct_hash(sys_struct)
    # Initialize with an empty, but now fully typed, SerializedModel.
    serialized_model = SerializedModel(; set_hash, sys_struct_hash)
    return SymbolicAWEModel(; sys_struct, serialized_model, kwargs...)
end

"""
    write_aero_forces!(ss, sys_struct) -> ss

Fill `ss.aero_force_x/y/z` with the aerodynamic force each point carries, world
frame. A point's own `aero_force_b` holds it only for `PARTICLE_DYNAMICS` wing
nodes, so a scattering mode is asked for its own distribution instead.
"""
function write_aero_forces!(ss::SysState, sys_struct::SystemStructure)
    fill!(ss.aero_force_x, 0.0)
    fill!(ss.aero_force_y, 0.0)
    fill!(ss.aero_force_z, 0.0)
    for wing in sys_struct.wings
        rot = quaternion_to_rotation_matrix(wing.Q_b_to_w)
        for (idx, force_b) in aero_point_forces(wing.aero, wing, sys_struct)
            force = rot * force_b
            ss.aero_force_x[idx] += force[1]
            ss.aero_force_y[idx] += force[2]
            ss.aero_force_z[idx] += force[3]
        end
    end
    return ss
end

"""
    update_sys_state!(ss::SysState, s::SymbolicAWEModel, zoom=1.0)

Updates a `SysState` object with the current state values from the `SymbolicAWEModel`.

This function takes the raw data from the model's internal integrator and populates
the fields of the user-friendly `SysState` struct, converting units (e.g., radians
to degrees) and calculating derived values like AoA and roll/pitch/yaw angles.

# Arguments
- `ss::SysState`: The state struct to be updated.
- `s::SymbolicAWEModel`: The source model.
- `zoom::SimFloat=1.0`: A scaling factor for the position coordinates.
"""
function update_sys_state!(ss::SysState, sam::SymbolicAWEModel, zoom=1.0)
    ss.time = isnothing(sam.integrator) ? 0.0 : sam.integrator.t # Use integrator time
    (; points, twist_surfaces, winches, wings, tethers,
       bodies) = sam.sys_struct

    for (ti, tether) in enumerate(tethers)
        ss.l_tether[ti] = tether.len
    end
    for winch in winches
        isempty(winch.tether_idxs) && continue
        ss.v_reelout[winch.idx] = winch.vel
        ss.winch_force[winch.idx] = norm(winch.force)
        ss.set_torque[winch.idx] = winch.set_value
    end
    # Unlike flap_angle, these are the integrated state itself.
    for twist_surface in twist_surfaces
        ss.twist_angles[twist_surface.idx] = twist_surface.twist
        ss.twist_vel[twist_surface.idx] = twist_surface.twist_ω
    end
    if length(twist_surfaces) > 0
        outer = length(twist_surfaces)
        ss.depower = rad2deg(mean(ss.twist_angles))
        ss.steering = rad2deg(ss.twist_angles[outer] - ss.twist_angles[1])
    end
    if length(wings) > 0
        wing = wings[1]
        ss.acc = norm(wing.acc_w) # Use the norm of the wing's acceleration vector
        ss.orient .= wing.Q_b_to_w
        ss.turn_rates .= wing.turn_rate
        ss.elevation = wing.elevation
        ss.azimuth = wing.azimuth
        ss.heading = wing.heading
        ss.course = wing.course
        # Apparent Wind and Aerodynamics
        ss.v_app = norm(wing.va_b)
        ss.v_wind_kite .= wing.v_wind
        # Calculate AoA and Side Slip from apparent wind in body frame
        if ss.v_app > 1e-6 # Avoid division by zero
            ss.AoA = calc_aoa(wing.aero, wing)
            ss.side_slip = calc_side_slip(wing)
        else
            ss.AoA = NaN       # Apparent wind too small to define AoA
            ss.side_slip = NaN # Side slip not defined for zero apparent wind
        end
        ss.aero_force_b .= wing.aero_force_b
        ss.aero_moment_b .= wing.aero_moment_b
        ss.tether_induced_force .= wing.tether_force
        ss.tether_induced_moment .= wing.tether_moment
        ss.vel_kite .= wing.vel_w
    end
    for point in points
        ss.X[point.idx] = point.pos_w[1] * zoom
        ss.Y[point.idx] = point.pos_w[2] * zoom
        ss.Z[point.idx] = point.pos_w[3] * zoom
        ss.VX[point.idx] = point.vel_w[1]
        ss.VY[point.idx] = point.vel_w[2]
        ss.VZ[point.idx] = point.vel_w[3]
        ss.drag_force_x[point.idx] = point.drag_force[1]
        ss.drag_force_y[point.idx] = point.drag_force[2]
        ss.drag_force_z[point.idx] = point.drag_force[3]
    end
    for segment in sam.sys_struct.segments
        ss.spring_force[segment.idx] = segment.force
    end
    write_aero_forces!(ss, sam.sys_struct)

    # Store per-mode aero log points (e.g. VSM panel corners) in world frame
    corner_idx = length(points)
    for wing in wings
        corner_idx = write_aero_log_points!(wing.aero, wing, sam.sys_struct,
                                            ss, corner_idx, zoom)
    end
    write_flap_deflections!(ss, sam.sys_struct)

    # Orientation frames are ordered wings-first, then bodies (slots after corners).
    slots = position_slots(sam.sys_struct)
    n_wings = length(wings)
    for wing in wings
        # frame 1 == legacy `orient`
        write_body_state!(ss, wing, slots.wings[wing.idx], wing.idx, zoom)
    end
    for rigid_body in bodies
        write_body_state!(ss, rigid_body, slots.bodies[rigid_body.idx],
                          n_wings + rigid_body.idx, zoom)
    end
    for pulley in sam.sys_struct.pulleys
        ss.pulley_len[pulley.idx] = pulley.len
        ss.pulley_vel[pulley.idx] = pulley.vel
    end

    ss.v_wind_gnd .= sam.set.wind_vec
    nothing
end

"""
    write_body_state!(ss, body, slot, frame, zoom)

Store a body's pose and rates: position and velocity in `X/Y/Z` and `VX/VY/VZ`
at `slot`, orientation and turn rate in `orients` and `turn_rate_x/y/z` at
`frame`. Wings occupy both a wing slot and a body slot, so both are written.
"""
function write_body_state!(ss, body, slot, frame, zoom)
    ss.X[slot] = body.pos_w[1] * zoom
    ss.Y[slot] = body.pos_w[2] * zoom
    ss.Z[slot] = body.pos_w[3] * zoom
    ss.VX[slot], ss.VY[slot], ss.VZ[slot] = body.vel_w
    ss.orients[frame] .= body.Q_b_to_w
    ss.turn_rate_x[frame] = body.ω_b[1]
    ss.turn_rate_y[frame] = body.ω_b[2]
    ss.turn_rate_z[frame] = body.ω_b[3]
    return nothing
end

"""
    position_slots(sys_struct) -> NamedTuple

Index layout of a `SysState`'s `X/Y/Z` position arrays for this model:
structural `points`, then VSM `panel_corners`, then `wings` origins, then
standalone rigid `bodies` origins. Each field is the `UnitRange` of slots for
that group (empty if none); `total` is the position count. Orientation frames
(`orients`) are laid out wings-first, then bodies — i.e. wing `w` uses frame
`w`, rigid body `b` uses frame `n_wings + b`.

```julia
slots = position_slots(sam.sys_struct)
sam.sys_struct.bodies[2]  # logged at X/Y/Z slot slots.bodies[2]
```
"""
function position_slots(sys_struct)
    n_points = length(sys_struct.points)
    n_corners = count_aero_log_points(sys_struct.wings)
    n_wings = length(sys_struct.wings)
    n_bodies = length(sys_struct.bodies)
    base_wings = n_points + n_corners
    base_bodies = base_wings + n_wings
    return (points        = 1:n_points,
            panel_corners = (n_points + 1):(n_points + n_corners),
            wings         = (base_wings + 1):(base_wings + n_wings),
            bodies        = (base_bodies + 1):(base_bodies + n_bodies),
            total         = base_bodies + n_bodies)
end

"""Number of orientation frames (wings + rigid bodies, at least 1)."""
n_orient_frames(sys_struct) = max(1,
    length(sys_struct.wings) + length(sys_struct.bodies))

"""
    SysState(s::SymbolicAWEModel, zoom=1.0)

Constructs a `SysState` object from a `SymbolicAWEModel`.

This is a convenience constructor that creates a new `SysState` object and populates it
with the current state of the provided model.

# Arguments
- `s::SymbolicAWEModel`: The source model.
- `zoom::SimFloat=1.0`: A scaling factor for the position coordinates.

# Returns
- `SysState`: A new state struct representing the current model state.
"""
function SysState(s::SymbolicAWEModel, zoom=1.0; precision=KiteUtils.MyFloat)
    ss = SysState(position_slots(s.sys_struct).total;
                  orients = n_orient_frames(s.sys_struct),
                  deflections = n_flap_deflections(s.sys_struct),
                  pulleys = length(s.sys_struct.pulleys),
                  winches = length(s.sys_struct.winches),
                  tethers = length(s.sys_struct.tethers),
                  segments = length(s.sys_struct.segments), precision)
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
logger = Logger(sam, 1000)  # Instead of Logger(length(sam.sys_struct.points), 1000)
```
"""
function KiteUtils.Logger(sam::SymbolicAWEModel, steps::Int;
                          precision=KiteUtils.MyFloat)
    slots = position_slots(sam.sys_struct)
    return Logger(slots.total, steps;
                  orients = n_orient_frames(sam.sys_struct),
                  deflections = n_flap_deflections(sam.sys_struct),
                  pulleys = length(sam.sys_struct.pulleys),
                  winches = length(sam.sys_struct.winches),
                  tethers = length(sam.sys_struct.tethers),
                  segments = length(sam.sys_struct.segments),
                  precision)
end

"""
    next_step!(s::SymbolicAWEModel, integrator::ODEIntegrator; set_values, dt, vsm_interval)

Take a simulation step, using the provided integrator.

This is a convenience method that calls the main `next_step!` function.
"""
function next_step!(
    s::SymbolicAWEModel,
    integrator::OrdinaryDiffEqCore.ODEIntegrator;
    set_values=nothing, dt=1/s.set.sample_freq,
    vsm_interval=1, vsm_min_wind=0.5
)
    !(s.integrator === integrator) && error(
        "The ODEIntegrator doesn't belong to " *
        "the SymbolicAWEModel")
    next_step!(s; set_values, dt, vsm_interval, vsm_min_wind)
end

"""
    next_step!(s::SymbolicAWEModel; set_values, dt,
               vsm_interval, vsm_min_wind)

Take a simulation step forward in time.

Advances the simulation by one time step, optionally
updating control inputs and re-linearizing the VSM
model. Then updates the `SystemStructure` with the new
state from the ODE integrator. Throws an error if the
solver returns an unstable retcode.

# Keyword Arguments
- `set_values=nothing`: Control input values.
    If `nothing`, current values are used.
- `dt=1/s.set.sample_freq`: Time step size [s].
- `vsm_interval=1`: Steps between VSM
    re-linearization. 0 disables re-linearization.
- `vsm_min_wind=0.5`: Minimum apparent wind [m/s] for
    a VSM solve. Below this the solver is skipped and
    the wing's aero outputs are zeroed, since the
    solver fails to converge or returns a Jacobian
    whose norm grows as 1/|va|.
"""
function next_step!(sam::SymbolicAWEModel;
    set_values=nothing, dt=1/sam.set.sample_freq,
    vsm_interval=1, vsm_min_wind=0.5
)
    prob = sam.prob
    integrator = sam.integrator
    if isnothing(integrator)
        error("next_step! called before init!: integrator is not initialized")
    end
    if (isnothing(set_values))
        set_values = [winch.set_value
            for winch in sam.sys_struct.winches]
    end
    if prob isa ProbWithAttributes && !isnothing(prob.set_set_values)
        prob.set_set_values(integrator, set_values)
    end
    if prob isa ProbWithAttributes
        sync_params!(prob.param_sync, integrator, sam.sys_struct)
    end

    sam.t_0 = integrator.t
    sam.t_step = @elapsed OrdinaryDiffEqCore.step!(integrator, dt, true)
    if !successful_retcode(integrator.sol)
        throw(AssertionError("Solver unstable at t=" *
            "$(round(integrator.t; digits=4))" *
            ": $(integrator.sol.retcode)"))
    end
    sam.iter += 1
    if prob isa ProbWithAttributes
        update_sys_struct!(prob, integrator, sam.sys_struct)
        if vsm_interval != 0 && sam.iter % vsm_interval == 0 &&
                has_vsm_wing(sam.sys_struct)
            sam.t_vsm = @elapsed begin
                refresh_aero!(sam; vsm_min_wind)
                sync_params!(prob.param_sync, integrator, sam.sys_struct)
            end
        end
    end
    return nothing
end

"""
    update_sys_struct!(s::SymbolicAWEModel, sys_struct::SystemStructure, integ=s.integrator)

Updates the high-level `SystemStructure` from the low-level integrator state vector.

This function reads the raw state vector from the ODE integrator and uses the generated
getter functions to populate the human-readable fields in the `SystemStructure`. This
synchronization step is crucial for making the simulation results accessible.
"""
function update_sys_struct!(prob::ProbWithAttributes,
                            integ::OrdinaryDiffEqCore.ODEIntegrator,
                            sys_struct::SystemStructure)
    prob.get_all_state(integ, sys_struct)
    refresh_deformed_positions!(sys_struct.points, sys_struct.bodies)
    return nothing
end

"""
    get_model_name(set::Settings, sys_struct::SystemStructure; precompile=false,
                   sparse=false, analytic_jacobian=false, backend=MonolithBackend())

Constructs a unique filename for the serialized model based on its configuration.
The filename includes the SymbolicAWEModels version, Julia version, physical model,
wing type, dynamics type, and component counts to ensure that the correct cached
model is loaded. `sparse`, `analytic_jacobian` and `backend` are part of the name
because the cached `ODEProblem` carries its Jacobian prototype, its Jacobian and its
backend's assembly, so those builds are different artefacts; naming them apart keeps
all of them on disk rather than invalidating one with the other.
"""
function get_model_name(set::Settings, sys_struct::SystemStructure; precompile=false,
                        sparse=false, analytic_jacobian=false,
                        backend=MonolithBackend())
    suffix = ""
    pkg_ver = pkgversion(SymbolicAWEModels)
    ver = "$(VERSION.major).$(VERSION.minor)"
    if precompile
        suffix = ".default"
    end

    # Determine wing type and aero mode
    dynamics_types = [wing.dynamics_type for wing in sys_struct.wings]
    dynamics_type_str = if isempty(dynamics_types)
        "no_wing"
    elseif all(wt -> wt === RIGID_DYNAMICS, dynamics_types)
        "rigid"
    elseif all(wt -> wt === PARTICLE_DYNAMICS, dynamics_types)
        "particle"
    else
        "mixed"
    end

    aero_tags = unique(aero_mode_tag(wing.aero) for wing in sys_struct.wings)
    aero_mode_str = if isempty(aero_tags)
        ""
    elseif length(aero_tags) == 1
        only(aero_tags)
    else
        "mixed_aero_modes"
    end

    dynamics_type = "dynamic"

    # Count components
    n_points = length(sys_struct.points)
    n_segments = length(sys_struct.segments)
    n_twist_surfaces = length(sys_struct.twist_surfaces)
    n_wings = length(sys_struct.wings)
    n_winches = length(sys_struct.winches)
    n_bodies = length(sys_struct.bodies)
    body_tag = n_bodies > 0 ? "_$(n_bodies)bdy" : ""
    sparse_tag = sparse ? "_sparse" : ""
    jacobian_tag = analytic_jacobian ? "_analytic" : ""

    return "model_v$(pkg_ver)_jl$(ver)_$(set.physical_model)_$(dynamics_type_str)_$(aero_mode_str)_$(dynamics_type)_$(n_points)pnt_$(n_segments)seg_$(n_twist_surfaces)grp_$(n_wings)wng_$(n_winches)wch$(body_tag)$(sparse_tag)$(jacobian_tag)$(backend_tag(backend)).bin$suffix"
end

"""
    calc_steady_torque(sam::SymbolicAWEModel)

Calculates the torque for each winch that results in zero acceleration (steady state).
"""
function calc_steady_torque(sam::SymbolicAWEModel)
    return calc_steady_torque(sam.sys_struct)
end
function calc_steady_torque(sys_struct::SystemStructure)
    torques = [-winch.drum_radius / winch.gear_ratio * norm(winch.force) +
               winch.friction for winch in sys_struct.winches]
    return torques
end

"""
    calc_winch_force(sys, winch_vel, winch_acc, set_values)

Calculate the tensile force on each winch tether from its motion and
motor torque, using the default-component motor dynamics inverted for
the force connector.

Reads `friction` from the live winch struct (populated each step from
the component output), so the formula matches whatever the component
reported even if the component overrides friction.
"""
function calc_winch_force(sys::SystemStructure,
        winch_vel, winch_acc, set_values)
    winches = sys.winches
    winch_force = zeros(length(winches))
    for i in eachindex(winches)
        (; gear_ratio, drum_radius, inertia_total, friction) = winches[i]
        ratio = drum_radius / gear_ratio
        α_motor = winch_acc[i] / ratio
        tau_total = α_motor * inertia_total
        winch_force[i] = (-set_values[i] + tau_total + friction) / ratio
    end
    return winch_force
end

"""
    calc_aoa(s::SymbolicAWEModel)

Angle of attack [rad] of the first wing, dispatched on its aero mode
([`calc_aoa(::AbstractAeroModel, wing)`](@ref)). `NaN` if the mode defines no AoA.
"""
function calc_aoa(sam::SymbolicAWEModel)
    wing = sam.sys_struct.wings[1]
    return calc_aoa(wing.aero, wing)
end

"""
    unstretched_length(s::SymbolicAWEModel)

Returns the unstretched tether length [m] for each tether.
"""
unstretched_length(sam::SymbolicAWEModel) = [tether.len for tether in sam.sys_struct.tethers]

"""
    tether_length(s::SymbolicAWEModel)

Returns the current stretched tether length [m] for each tether.
"""
tether_length(sam::SymbolicAWEModel) = [tether.stretched_len for tether in sam.sys_struct.tethers]

"""
    calc_height(s::SymbolicAWEModel)

Returns the height (z-position) [m] of the wing.
"""
calc_height(sam::SymbolicAWEModel) = sam.sys_struct.wings[1].pos_w[3]

"""
    winch_force(s::SymbolicAWEModel)

Returns the winch force [N] for each winch.
"""
winch_force(sam::SymbolicAWEModel) = [norm(winch.force) for winch in sam.sys_struct.winches]

"""
    spring_forces(s::SymbolicAWEModel)

Returns the spring force [N] for each tether segment.
"""
spring_forces(sam::SymbolicAWEModel) = [segment.force for segment in sam.sys_struct.segments]

"""
    pos(s::SymbolicAWEModel)

Returns a vector of the position vectors [m] for each point in the system.
"""
pos(sam::SymbolicAWEModel) = [point.pos_w for point in sam.sys_struct.points]


