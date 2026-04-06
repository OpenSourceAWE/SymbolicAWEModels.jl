# Copyright (c) 2025 Bart van de Lint and Uwe Fechner
# SPDX-License-Identifier: MPL-2.0

const LinType = @NamedTuple{A::Matrix{SimFloat}, B::Matrix{SimFloat}, C::Matrix{SimFloat}, D::Matrix{SimFloat}}
const GetSetNothing = Union{AbstractIndexer, Nothing}

"""
    @with_kw struct ProbWithAttributes{...}

A container for the main Ordinary Differential Equation (ODE) problem and its
associated getter and setter functions for the full, nonlinear physical state.
"""
@with_kw struct ProbWithAttributes{Prob, SetSys, SetSetValues, SetSet,
                                  GetSetValues, GetWingState, GetVsmY, GetSegmentState,
                                  GetWinchState, GetTetherState, GetStructState, GetPointState,
                                  GetPulleyState, GetGroupState}
    "The ODE problem for the full nonlinear model."
    prob::Prob

    # Setters for the ODE
    "Setter for the system parameters."
    set_sys::SetSys
    "Setter for the control input values."
    set_set_values::SetSetValues
    "Setter for general settings."
    set_set::SetSet

    # Getters for the ODE state
    get_set_values::GetSetValues
    get_wing_state::GetWingState
    get_vsm_y::GetVsmY
    get_segment_state::GetSegmentState
    get_winch_state::GetWinchState
    get_tether_state::GetTetherState
    get_struct_state::GetStructState
    get_point_state::GetPointState
    get_pulley_state::GetPulleyState
    get_group_state::GetGroupState
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
@with_kw struct LinProbWithAttributes{Prob, SetSetValues, SetSys, SetSet}
    "Linearization problem of the mtk model."
    prob::Prob

    # Setters for the linearization
    set_set_values::SetSetValues
    set_sys::SetSys
    set_set::SetSet
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
@with_kw mutable struct SerializedModel{D<:AbstractVector, G<:AbstractVector}
    set_hash::Vector{UInt8}
    sys_struct_hash::Vector{UInt8}
    "Unsimplified system of the mtk model"
    full_sys::Union{ModelingToolkit.System, Nothing} = nothing
    defaults::D = Pair{Num, Any}[]
    guesses::G = Pair{Num, Any}[]
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
    "Reference to the settings struct"
    set::Settings
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
    "Vector of tether length set-points"
    set_tether_len::Vector{SimFloat} = zeros(SimFloat, 3)
end

"""
    _SAM_FIELDS

Tuple of field names that are direct fields of `SymbolicAWEModel` (as opposed to fields
delegated to the nested `serialized_model`). Used by `getproperty` and `setproperty!`
to dispatch field access correctly.
"""
const _SAM_FIELDS = (:set, :sys_struct, :serialized_model, :integrator, :t_0, :iter, :t_vsm, :t_step, :set_tether_len)

"""
    Base.getproperty(sam::SymbolicAWEModel, sym::Symbol)

Overloads `getproperty` to allow direct access to fields within the nested `serialized_model`.
This provides a convenient way to access compiled functions and other model
components without explicitly referencing `sam.serialized_model`.
"""

function Base.getproperty(sam::SymbolicAWEModel, sym::Symbol)
    if sym === :am
        getfield(sam, :sys_struct).am
    elseif sym in _SAM_FIELDS
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
    if sym in _SAM_FIELDS
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
    set_hash = get_set_hash(set)
    sys_struct_hash = get_sys_struct_hash(sys_struct)
    # Initialize with an empty, but now fully typed, SerializedModel.
    serialized_model = SerializedModel(; set_hash, sys_struct_hash)
    return SymbolicAWEModel(; set, sys_struct, serialized_model, kwargs...)
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
    @unpack points, groups, segments, pulleys, winches, wings = sam.sys_struct

    # Get the state vectors from the integrator
    if length(winches) > 0
        for winch in winches
            ss.l_tether[winch.idx] = winch.tether_len
            ss.v_reelout[winch.idx] = winch.tether_vel
            ss.winch_force[winch.idx] = norm(winch.force)
            ss.set_torque[winch.idx] = winch.set_value
        end
    end
    if length(groups) > 0
        # Only fill up to the size of ss.twist_angles (typically 4)
        max_groups = min(length(groups), length(ss.twist_angles))
        for group in groups[1:max_groups]
            ss.twist_angles[group.idx] = group.twist
        end
        ss.depower = rad2deg(mean(ss.twist_angles[1:max_groups])) # Average twist for depower
        ss.steering = rad2deg(ss.twist_angles[max_groups] - ss.twist_angles[1])
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
            if wing isa VSMWing
                aoa_raw = wing.vsm_solver.sol.alpha_geometric_dist[length(wing.vsm_solver.sol.alpha_dist) ÷ 2 + 
                          (length(wing.vsm_solver.sol.alpha_dist) % 2)] # version-2, likely with induction
                ss.AoA = mod(aoa_raw + π, 2π) - π  # Wrap to [-π, π]
                ss.side_slip = asin(wing.va_b[2] / norm(wing.va_b))
            else
                ss.AoA = NaN # AoA not defined for non-VSM wings
                ss.side_slip = NaN # Side slip not defined for non-VSM wings
            end
            
        else
            ss.AoA = NaN       # Apparent wind too small to define AoA
            ss.side_slip = NaN # Side slip not defined for zero apparent wind
        end
        ss.aero_force_b .= wing.aero_force_b
        ss.aero_moment_b .= wing.aero_moment_b
        ss.tether_induced_force .= wing.tether_force
        ss.tether_induced_moment .= wing.tether_moment
        ss.vel_kite .= wing.vel_w
        # Calculate Roll, Pitch, Yaw from Quaternion
        q = wing.Q_b_to_w
        sinr_cosp = 2 * (q[1] * q[2] + q[3] * q[4])
        cosr_cosp = 1 - 2 * (q[2] * q[2] + q[3] * q[3])
        ss.roll = atan(sinr_cosp, cosr_cosp)
        sinp = 2 * (q[1] * q[3] - q[4] * q[2])
        ss.pitch = abs(sinp) >= 1 ? copysign(pi / 2, sinp) : asin(sinp)
        siny_cosp = 2 * (q[1] * q[4] + q[2] * q[3])
        cosy_cosp = 1 - 2 * (q[3] * q[3] + q[4] * q[4])
        ss.yaw = atan(siny_cosp, cosy_cosp)
    end
    for point in points
        ss.X[point.idx] = point.pos_w[1] * zoom
        ss.Y[point.idx] = point.pos_w[2] * zoom
        ss.Z[point.idx] = point.pos_w[3] * zoom
    end

    # Store VSM panel corner positions in world frame
    corner_idx = length(points)
    for wing in wings
        wing isa VSMWing || continue
        R_b_to_w = wing.R_b_to_w::Matrix{SimFloat}
        for panel in wing.vsm_aero.panels
            for j in 1:4
                corner_idx += 1
                # Transform from body frame to world frame
                corner_b = panel.corner_points[:, j]
                corner_w = wing.pos_w + R_b_to_w * corner_b
                ss.X[corner_idx] = corner_w[1] * zoom
                ss.Y[corner_idx] = corner_w[2] * zoom
                ss.Z[corner_idx] = corner_w[3] * zoom
            end
        end
    end

    ss.v_wind_gnd .= sam.sys_struct.wind_vec_gnd
    nothing
end

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
function SysState(s::SymbolicAWEModel, zoom=1.0)
    # Calculate total points: regular points + 4 corners per panel
    n_points = length(s.sys_struct.points)
    n_panel_corners = isempty(s.sys_struct.wings) ? 0 : sum(
        length(wing.vsm_aero.panels) * 4 for wing in s.sys_struct.wings if wing isa VSMWing;
        init=0
    )
    total_points = n_points + n_panel_corners
    ss = SysState{total_points}()
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
function KiteUtils.Logger(sam::SymbolicAWEModel, steps::Int)
    # Calculate total points: regular points + 4 corners per panel
    n_points = length(sam.sys_struct.points)
    n_panel_corners = isempty(sam.sys_struct.wings) ? 0 : sum(
        length(wing.vsm_aero.panels) * 4 for wing in sam.sys_struct.wings if wing isa VSMWing;
        init=0
    )
    total_points = n_points + n_panel_corners
    return Logger(total_points, steps)
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
    vsm_interval=1
)
    !(s.integrator === integrator) && error(
        "The ODEIntegrator doesn't belong to " *
        "the SymbolicAWEModel")
    next_step!(s; set_values, dt, vsm_interval)
end

"""
    next_step!(s::SymbolicAWEModel; set_values, dt,
               vsm_interval)

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
"""
function next_step!(sam::SymbolicAWEModel;
    set_values=nothing, dt=1/sam.set.sample_freq,
    vsm_interval=1
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

    sam.t_0 = integrator.t
    sam.t_step = @elapsed step!(integrator, dt, true)
    if !successful_retcode(integrator.sol)
        error("Solver unstable at t=" *
            "$(round(integrator.t; digits=4))" *
            ": $(integrator.sol.retcode)")
    end
    sam.iter += 1
    if prob isa ProbWithAttributes
        update_sys_struct!(prob, integrator, sam.sys_struct)
        if vsm_interval != 0 && sam.iter % vsm_interval == 0
            sam.t_vsm = @elapsed update_vsm!(sam, prob)
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
    @unpack points, groups, segments, pulleys, winches, tethers, wings = sys_struct
    pos, vel, force, va_b, total_mass = prob.get_point_state(integ)
    for point in points
        point.pos_w .= pos[:, point.idx]
        point.vel_w .= vel[:, point.idx]
        point.force .= force[:, point.idx]
        point.va_b .= va_b[:, point.idx]
        point.total_mass = total_mass[point.idx]
    end
    if length(pulleys) > 0
        len, vel = prob.get_pulley_state(integ)
        for pulley in pulleys
            pulley.len = len[pulley.idx]
            pulley.vel = vel[pulley.idx]
        end
    end
    if length(segments) > 0
        spring_force, len = prob.get_segment_state(integ)
        for segment in segments
            segment.force = spring_force[segment.idx]
            segment.len = len[segment.idx]
        end
    end
    if length(groups) > 0
        twist, twist_ω, tether_force, tether_moment, aero_moment = prob.get_group_state(integ)
        for group in groups
            group.twist = twist[group.idx]
            group.twist_ω = twist_ω[group.idx]
            group.tether_force = tether_force[group.idx]
            group.tether_moment = tether_moment[group.idx]
            group.aero_moment = aero_moment[group.idx]
        end
    end
    if length(winches) > 0
        tether_len, tether_vel, tether_acc, set_value, winch_force_vec, friction =
            prob.get_winch_state(integ)
        for winch in winches
            winch.tether_len = tether_len[winch.idx]
            winch.tether_vel = tether_vel[winch.idx]
            winch.tether_acc = tether_acc[winch.idx]
            winch.set_value = set_value[winch.idx]
            winch.force .= winch_force_vec[:, winch.idx]
            winch.friction = friction[winch.idx]
        end
    end
    if length(tethers) > 0
        stretched_len = prob.get_tether_state(integ)
        for tether in tethers
            tether.stretched_len = stretched_len[tether.idx]
        end
    end
    if length(wings) > 0
        wing_state = prob.get_wing_state(integ)
        Q_b_to_w, ω_b, pos_w, vel_w, acc_w,
            va_b, v_wind,
            aero_force_b, aero_moment_b,
            tether_moment, tether_force,
            elevation, elevation_vel, elevation_acc,
            azimuth, azimuth_vel, azimuth_acc,
            heading, turn_rate, turn_acc,
            course, aoa,
            com_w_v, com_vel_v,
            Q_p_to_w_v, ω_p_v = wing_state
        for wing in wings
            # Body frame output
            wing.Q_b_to_w .= Q_b_to_w[:, wing.idx]
            wing.ω_b .= ω_b[:, wing.idx]
            wing.pos_w .= pos_w[:, wing.idx]
            wing.vel_w .= vel_w[:, wing.idx]
            wing.acc_w .= acc_w[:, wing.idx]
            wing.va_b .= va_b[:, wing.idx]
            wing.v_wind .= v_wind[:, wing.idx]
            wing.aero_force_b .=
                aero_force_b[:, wing.idx]
            wing.aero_moment_b .=
                aero_moment_b[:, wing.idx]
            wing.tether_moment .=
                tether_moment[:, wing.idx]
            wing.tether_force .=
                tether_force[:, wing.idx]
            wing.elevation =
                elevation[wing.idx]
            wing.elevation_vel =
                elevation_vel[wing.idx]
            wing.elevation_acc =
                elevation_acc[wing.idx]
            wing.azimuth = azimuth[wing.idx]
            wing.azimuth_vel =
                azimuth_vel[wing.idx]
            wing.azimuth_acc =
                azimuth_acc[wing.idx]
            wing.heading = heading[wing.idx]
            wing.turn_rate .=
                turn_rate[:, wing.idx]
            wing.turn_acc .=
                turn_acc[:, wing.idx]
            wing.course = course[wing.idx]
            wing.aoa = aoa[wing.idx]
            # Principal frame state
            wing.com_w .= com_w_v[:, wing.idx]
            wing.com_vel .=
                com_vel_v[:, wing.idx]
            wing.Q_p_to_w .= Q_p_to_w_v[:, wing.idx]
            wing.ω_p .= ω_p_v[:, wing.idx]
        end
    end
    sys_struct.wind_vec_gnd .= prob.get_struct_state(integ)
    return nothing
end

"""
    get_model_name(set::Settings, sys_struct::SystemStructure; precompile=false)

Constructs a unique filename for the serialized model based on its configuration.
The filename includes the Julia version, physical model, wing type, dynamics type,
and component counts to ensure that the correct cached model is loaded.
"""
function get_model_name(set::Settings, sys_struct::SystemStructure; precompile=false)
    suffix = ""
    ver = "$(VERSION.major).$(VERSION.minor)"
    if precompile
        suffix = ".default"
    end

    # Determine wing type and aero mode
    wing_types = [wing.wing_type for wing in sys_struct.wings]
    wing_type_str = if isempty(wing_types)
        "no_wing"
    elseif all(wt -> wt == QUATERNION, wing_types)
        "quat"
    elseif all(wt -> wt == REFINE, wing_types)
        "refine"
    else
        "mixed"
    end

    aero_modes = [wing.aero_mode for wing in sys_struct.wings]
    aero_mode_str = if isempty(aero_modes)
        ""
    elseif all(m -> m == AERO_LINEARIZED, aero_modes)
        "lin"
    elseif all(m -> m == AERO_DIRECT, aero_modes)
        "dir"
    elseif all(m -> m == AERO_NONE, aero_modes)
        "none"
    else
        "mixed_aero_modes"
    end

    dynamics_type = ifelse(set.quasi_static, "static", "dynamic")

    # Count components
    n_points = length(sys_struct.points)
    n_segments = length(sys_struct.segments)
    n_groups = length(sys_struct.groups)
    n_wings = length(sys_struct.wings)
    n_winches = length(sys_struct.winches)

    return "model_$(ver)_$(set.physical_model)_$(wing_type_str)_$(aero_mode_str)_$(dynamics_type)_$(n_points)pnt_$(n_segments)seg_$(n_groups)grp_$(n_wings)wng_$(n_winches)wch.bin$suffix"
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
    calc_winch_force(tether_vel, tether_acc, motor_torque, set)

Calculate the tensile force on the winch tether based on its motion and motor torque.

This function uses a settings object to define the physical parameters of the winch.

# Arguments
- `tether_vel`: The velocity of the tether [m/s].
- `tether_acc`: The acceleration of the tether [m/s²].
- `motor_torque`: The torque applied by the motor [Nm].
- `set`: A settings struct.

# Returns
- The calculated force on the winch tether [N].
"""
function calc_winch_force(sys::SystemStructure,
        tether_vel, tether_acc, set_values)
    winches = sys.winches
    smooth_sign(x, eps) = x / sqrt(x * x + eps * eps)
    winch_force = zeros(length(winches))
    for i in eachindex(winches)
        @unpack gear_ratio, drum_radius, f_coulomb,
            c_vf, inertia_total,
            friction_epsilon = winches[i]
        ω_motor = gear_ratio / drum_radius * tether_vel[i]
        tau_friction =
            smooth_sign(ω_motor, friction_epsilon) *
            f_coulomb * drum_radius / gear_ratio +
            c_vf * ω_motor *
            drum_radius^2 / gear_ratio^2
        tau_motor = set_values[i] # set_value is the motor torque
        α_motor = tether_acc[i] / drum_radius * gear_ratio
        tau_total = α_motor * inertia_total
        winch_force[i] = (-tau_motor + tau_total + tau_friction) / drum_radius * gear_ratio
    end
    return winch_force
end

"""
    calc_aoa(s::SymbolicAWEModel)

Calculates the mean angle of attack [rad] over the wingspan from the VSM solver.
"""
function calc_aoa(sam::SymbolicAWEModel)
    wing = sam.sys_struct.wings[1]
    wing isa VSMWing || error("calc_aoa: wing[1] is not a VSMWing")
    alpha_array = wing.vsm_solver.sol.alpha_dist
    middle = length(alpha_array) ÷ 2
    return iseven(length(alpha_array)) ? (0.5 * (alpha_array[middle] + alpha_array[middle+1])) : alpha_array[middle+1]
end

"""
    unstretched_length(s::SymbolicAWEModel)

Returns the unstretched tether length [m] for each winch.
"""
unstretched_length(sam::SymbolicAWEModel) = [winch.tether_len for winch in sam.sys_struct.winches]

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

"""
    min_chord_len(s::SymbolicAWEModel)

Calculates the minimum chord length of the wing at the tip.
"""
function min_chord_len(sam::SymbolicAWEModel)
    min_len = Inf
    for wing in sam.sys_struct.wings
        wing isa VSMWing || continue
        vsm_wing = wing.vsm_wing
        if hasproperty(vsm_wing, :le_interp) && hasproperty(vsm_wing, :te_interp) && hasproperty(vsm_wing, :gamma_tip)
            le_pos = [vsm_wing.le_interp[i](vsm_wing.gamma_tip) for i in 1:3]
            te_pos = [vsm_wing.te_interp[i](vsm_wing.gamma_tip) for i in 1:3]
            min_len = min(norm(le_pos - te_pos), min_len)
        elseif hasproperty(vsm_wing, :unrefined_sections) && !isempty(vsm_wing.unrefined_sections)
            for section in vsm_wing.unrefined_sections
                chord = section.TE_point - section.LE_point
                min_len = min(norm(chord), min_len)
            end
        end
    end
    return min_len
end

"""
    set_depower_steering!(s::SymbolicAWEModel, depower, steering)

Sets the kite's depower and steering by adjusting the tether length set-points.
"""
function set_depower_steering!(sam::SymbolicAWEModel, depower, steering)
    len = sam.set_tether_len
    len .= [winch.tether_len for winch in sam.sys_struct.winches]
    depower *= min_chord_len(sam)
    steering *= min_chord_len(sam)
    len[2] = 0.5 * (2*depower + 2*len[1] + steering)
    len[3] = 0.5 * (2*depower + 2*len[1] - steering)
    return nothing
end

"""
    set_v_wind_ground!(s::SymbolicAWEModel, v_wind_gnd=s.set.v_wind, upwind_dir=-π/2)

Sets the ground wind speed [m/s] and upwind direction [rad] in the model.
"""
function set_v_wind_ground!(sam::SymbolicAWEModel, v_wind_gnd=sam.set.v_wind, upwind_dir=-pi/2)
    sam.set.v_wind = v_wind_gnd
    sam.set.upwind_dir = rad2deg(upwind_dir)
    local_prob = sam.prob
    if local_prob isa ProbWithAttributes
        local_prob.set_set(sam.integrator, sam.set)
    end
    return nothing
end
