# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

"""
    VortexStepMethod.RamAirWing(set::Settings; prn=true, kwargs...)

Create a `RamAirWing` geometry object from the settings provided.

This is a constructor helper that reads the model and foil file paths from the
`Settings` object and initializes the `RamAirWing` object from `VortexStepMethod.jl`.
"""
function VortexStepMethod.RamAirWing(set::Settings; prn=true, kwargs...)
    obj_path = joinpath(dirname(get_data_path()), set.model)
    dat_path = joinpath(dirname(get_data_path()), set.foil_file)
    if set.physical_model == "simple_ram"
        n_groups=2
    else
        n_groups=4
    end
    return RamAirWing(obj_path, dat_path; 
        mass=set.mass, crease_frac=set.crease_frac, n_groups, 
        align_to_principal=true, prn, kwargs...
    )
end

"""
    SegmentType `POWER_LINE` `STEERING_LINE` `BRIDLE`

Enumeration for the type of a tether segment.

# Elements
- `POWER_LINE`: A segment belonging to a main power line.
- `STEERING_LINE`: A segment belonging to a steering line.
- `BRIDLE`: A segment belonging to the bridle system.
"""
@enum SegmentType begin
    POWER_LINE
    STEERING_LINE
    BRIDLE
end

"""
    DynamicsType `DYNAMIC` `QUASI_STATIC` `WING` `STATIC`

Enumeration for the dynamic model governing a point's motion.

# Elements
- `DYNAMIC`: The point is a dynamic point mass, moving according to Newton's second law.
- `QUASI_STATIC`: The point's acceleration is constrained to zero, representing a force equilibrium.
- `WING`: The point is rigidly attached to a wing body and moves with it.
- `STATIC`: The point's position is fixed in the world frame.
"""
@enum DynamicsType begin
    DYNAMIC
    QUASI_STATIC
    WING
    STATIC
end

"""
    mutable struct Point

A point mass, representing a node in the mass-spring system.

$(TYPEDFIELDS)
"""
mutable struct Point
    const idx::Int16
    const transform_idx::Int16 # idx of wing used for initial orientation
    const wing_idx::Int16
    const pos_cad::KVec3
    const pos_b::KVec3 # pos relative to wing COM in body frame
    const pos_w::KVec3 # pos in world frame
    const vel_w::KVec3 # vel in world frame
    const disturb::KVec3 # disturbing force
    const force::KVec3
    const type::DynamicsType
    mass::SimFloat
    bridle_damping::SimFloat
    fix_sphere::Bool
end

"""
    Point(idx, pos_cad, type; wing_idx=1, vel_w=zeros(KVec3), transform_idx=1, mass=0.0)

Constructs a `Point` object, which can be of four different [`DynamicsType`](@ref)s:
- `STATIC`: The point does not move. ``\\ddot{\\mathbf{r}} = \\mathbf{0}``
- `DYNAMIC`: The point moves according to Newton's second law. ``\\ddot{\\mathbf{r}} = \\mathbf{F}/m``
- `QUASI_STATIC`: The acceleration is constrained to be zero by solving a nonlinear problem. ``\\mathbf{F}/m = \\mathbf{0}``
- `WING`: The point has a static position in the rigid body wing frame. ``\\mathbf{r}_w = \\mathbf{r}_{wing} + \\mathbf{R}_{b\\rightarrow w} \\mathbf{r}_b``

where:
- ``\\mathbf{r}`` is the point position vector
- ``\\mathbf{F}`` is the net force acting on the point
- ``m`` is the point mass
- ``\\mathbf{r}_w`` is the position in world frame
- ``\\mathbf{r}_{wing}`` is the wing center position
- ``\\mathbf{R}_{b\\rightarrow w}`` is the rotation matrix from body to world frame
- ``\\mathbf{r}_b`` is the position in body frame

# Arguments
- `idx::Int16`: Unique identifier for the point.
- `pos_cad::KVec3`: Position of the point in the CAD frame.
- `type::DynamicsType`: Dynamics type of the point (`STATIC`, `DYNAMIC`, etc.).

# Keyword Arguments
- `wing_idx::Int16=1`: Index of the wing this point is attached to.
- `vel_w::KVec3=zeros(KVec3)`: Initial velocity of the point in world frame.
- `transform_idx::Int16=1`: Index of the transform used for initial positioning.
- `mass::Float64=0.0`: Mass of the point [kg].
- `bridle_damping::Float64=0.0`: Damping coefficient for bridle points.
- `fix_sphere::Bool=false`: If true, constrains the point to a sphere.

# Returns
- `Point`: A new `Point` object.
"""
function Point(idx, pos_cad, type;
    wing_idx=1, vel_w=zeros(KVec3), transform_idx=1, 
    mass=0.0, bridle_damping=0.0, fix_sphere=false
)
    Point(idx, transform_idx, wing_idx, pos_cad, zeros(KVec3), zeros(KVec3),
        vel_w, zeros(KVec3), zeros(KVec3), type, mass, bridle_damping, fix_sphere)
end

"""
    mutable struct Group

A set of bridle lines that share the same twist angle and trailing edge angle.

$(TYPEDFIELDS)
"""
mutable struct Group
    const idx::Int16
    const point_idxs::Vector{Int16}
    const le_pos::KVec3 # point which the group rotates around under wing deformation
    const chord::KVec3 # chord vector in body frame which the group rotates around under wing deformation
    const y_airf::KVec3 # spanwise vector in local panel frame which the group rotates around under wing deformation
    const type::DynamicsType
    moment_frac::SimFloat
    twist::SimFloat
    twist_ω::SimFloat
    tether_force::SimFloat
    tether_moment::SimFloat
    aero_moment::SimFloat
end

"""
    Group(idx, point_idxs, vsm_wing::RamAirWing, gamma, type, moment_frac)

Constructs a `Group` object representing a collection of points on a kite body that share
a common twist deformation.

A `Group` models the local deformation of a kite wing section through twist dynamics.
All points within a group undergo the same twist rotation about the chord vector.

The governing equation is:
```math
\\begin{aligned}
\\tau = \\underbrace{\\sum_{i=1}^{4} r_{b,i} \\times (\\mathbf{F}_{b,i} \\cdot \\hat{\\mathbf{z}})}_{\\text{bridles}} + \\underbrace{r_a \\times (\\mathbf{F}_a \\cdot \\hat{\\mathbf{z}})}_{\\text{aero}}
\\end{aligned}
```

![System Overview](assets/group_slice.svg)

where:
- ``\\tau`` is the total torque about the twist axis
- ``r_{b,i}`` is the position vector of bridle point ``i`` relative to the twist center
- ``\\mathbf{F}_{b,i}`` is the force at bridle point ``i``
- ``\\hat{\\mathbf{z}}`` is the unit vector along the twist axis (chord direction)
- ``r_a`` is the position vector of the aerodynamic center relative to the twist center
- ``\\mathbf{F}_a`` is the aerodynamic force at the group's aerodynamic center

The group can have two [`DynamicsType`](@ref)s:
- `DYNAMIC`: the group rotates according to Newton's second law: ``I\\ddot{\\theta} = \\tau``
- `QUASI_STATIC`: the rotational acceleration is zero: ``\\tau = 0``

# Arguments
- `idx::Int16`: Unique identifier for the group.
- `point_idxs::Vector{Int16}`: Indices of points that move together with this group's twist.
- `vsm_wing::RamAirWing`: Wing geometry object used to extract local chord and spanwise vectors.
- `gamma`: Spanwise parameter (typically -1 to 1) defining the group's location along the wing.
- `type::DynamicsType`: Dynamics type (`DYNAMIC` for time-varying twist, `QUASI_STATIC` for equilibrium).
- `moment_frac::SimFloat`: Chordwise position (0=leading edge, 1=trailing edge) about which the group rotates.

# Returns
- `Group`: A new `Group` object with twist dynamics capability.
"""
function Group(idx, point_idxs, vsm_wing::RamAirWing, gamma, type, moment_frac)
    le_pos = [vsm_wing.le_interp[i](gamma) for i in 1:3]
    chord = [vsm_wing.te_interp[i](gamma) for i in 1:3] .- le_pos
    y_airf = normalize([vsm_wing.le_interp[i](gamma-0.01) for i in 1:3] - le_pos)
    Group(idx, point_idxs, le_pos, chord, y_airf, type, moment_frac, 0.0, 0.0, 0.0, 0.0, 0.0)
end

"""
    Group(idx, point_idxs, le_pos, chord, y_airf, type, moment_frac)

Inner constructor for a `Group` object. See [`Group`](@ref) for details.
"""
function Group(idx, point_idxs, le_pos, chord, y_airf, type, moment_frac)
    Group(idx, point_idxs, le_pos, chord, y_airf, type, moment_frac, 0.0, 0.0)
end

"""
    mutable struct Segment

A segment representing a spring-damper connection from one point to another.

$(TYPEDFIELDS)
"""
mutable struct Segment
    const idx::Int16
    const point_idxs::Tuple{Int16, Int16}
    axial_stiffness::SimFloat
    axial_damping::SimFloat
    l0::SimFloat
    compression_frac::SimFloat
    diameter::SimFloat
    len::SimFloat # current len of the segment
    force::SimFloat # current force in the segment
end

"""
    Segment(idx, point_idxs, axial_stiffness, axial_damping, diameter; l0, compression_frac)

Inner constructor for a `Segment` object. See [`Segment`](@ref) for details.
"""
function Segment(idx, point_idxs, axial_stiffness, axial_damping, diameter; 
    l0=zero(SimFloat), compression_frac=0.1
)
    Segment(idx, point_idxs, axial_stiffness, axial_damping, l0, compression_frac, 
        diameter, zero(SimFloat), zero(SimFloat))
end

"""
    Segment(idx, set, point_idxs, type; l0, compression_frac, axial_stiffness, axial_damping)

Constructs a `Segment` object representing an elastic spring-damper connection between two points.

The segment follows Hooke's law with damping and aerodynamic drag:

**Spring-Damper Force:**
```math
\\mathbf{F}_{spring} = \\left[k(l - l_0) - c\\dot{l}\\right]\\hat{\\mathbf{u}}
```

**Aerodynamic Drag:**
```math
\\mathbf{F}_{drag} = \\frac{1}{2}\\rho C_d A |\\mathbf{v}_a| \\mathbf{v}_{a,\\perp}
```

**Total Force:**
```math
\\mathbf{F}_{total} = \\mathbf{F}_{spring} + \\mathbf{F}_{drag}
```

where:
- ``k = \\frac{E \\pi d^2/4}{l}`` is the axial stiffness
- ``l`` is current length, ``l_0`` is unstretched length
- ``c = \\frac{\\xi}{c_{spring}} k`` is damping coefficient
- ``\\hat{\\mathbf{u}} = \\frac{\\mathbf{r}_2 - \\mathbf{r}_1}{l}`` is unit vector along segment
- ``\\dot{l} = (\\mathbf{v}_1 - \\mathbf{v}_2) \\cdot \\hat{\\mathbf{u}}`` is extension rate
- ``\\mathbf{v}_{a,\\perp}`` is apparent wind velocity perpendicular to segment

# Arguments
- `idx::Int16`: Unique identifier for the segment.
- `set::Settings`: The settings object containing material properties.
- `point_idxs::Tuple{Int16, Int16}`: Tuple containing the indices of the two points.
- `type::SegmentType`: Type of the segment (`POWER_LINE`, `STEERING_LINE`, `BRIDLE`).

# Keyword Arguments
- `l0::SimFloat=zero(SimFloat)`: Unstretched length [m]. Calculated from point positions if zero.
- `compression_frac::SimFloat=0.0`: Stiffness reduction factor in compression.
- `axial_stiffness::Float64=NaN`: Axial stiffness [N]. If `NaN`, it's calculated from settings.
- `axial_damping::Float64=NaN`: Axial damping [Ns]. If `NaN`, it's calculated from settings.

# Returns
- `Segment`: A new `Segment` object.
"""
function Segment(idx, set, point_idxs, type;
    l0=zero(SimFloat), compression_frac=0.0, axial_stiffness=NaN, axial_damping=NaN
)
    (type == BRIDLE) && (diameter = 0.001* set.bridle_tether_diameter)
    (type == POWER_LINE) && (diameter = 0.001* set.power_tether_diameter)
    (type == STEERING_LINE) && (diameter = 0.001* set.steering_tether_diameter)

    if isnan(axial_stiffness) || isnan(axial_damping)
        axial_stiffness = set.e_tether * (diameter/2)^2 * π
        if type == BRIDLE
            stiffness_frac = 0.01
        else
            stiffness_frac = 1.0
        end
        axial_stiffness = stiffness_frac * axial_stiffness

        axial_damping = (set.damping / set.c_spring) * axial_stiffness
    end

    Segment(idx, point_idxs, axial_stiffness, axial_damping, l0, compression_frac, 
        diameter, zero(SimFloat), zero(SimFloat))
end

"""
    mutable struct Pulley

A pulley described by two segments with the common point of the segments being the pulley.

$(TYPEDFIELDS)
"""
mutable struct Pulley
    const idx::Int16
    const segment_idxs::Tuple{Int16, Int16}
    const type::DynamicsType
    sum_len::SimFloat
    len::SimFloat
    vel::SimFloat
end

"""
    Pulley(idx, segment_idxs, type)

Constructs a `Pulley` object that enforces length redistribution between two segments.

The pulley constraint maintains constant total length while allowing force transmission:

**Constraint Equations:**
```math
l_1 + l_2 = l_{total} = \\text{constant}
```

**Force Balance:**
```math
F_{pulley} = F_1 - F_2
```

**Dynamics:**
```math
m\\ddot{l}_1 = F_{pulley} = F_1 - F_2
```

where:
- ``l_1, l_2`` are the lengths of connected segments
- ``F_1, F_2`` are the spring forces in the segments
- ``m = \\rho_{tether} \\pi (d/2)^2 l_{total}`` is the total mass of both segments
- ``\\dot{l}_1 + \\dot{l}_2 = 0`` (velocity constraint)

The pulley can have two [`DynamicsType`](@ref)s:
- `DYNAMIC`: the length redistribution follows Newton's second law: ``m\\ddot{l}_1 = F_1 - F_2``
- `QUASI_STATIC`: the forces are balanced instantaneously: ``F_1 = F_2``

# Arguments
- `idx::Int16`: Unique identifier for the pulley.
- `segment_idxs::Tuple{Int16, Int16}`: Tuple containing the indices of the two segments.
- `type::DynamicsType`: Dynamics type of the pulley (`DYNAMIC` or `QUASI_STATIC`).

# Returns
- `Pulley`: A new `Pulley` object.
"""
function Pulley(idx, segment_idxs, type)
    return Pulley(idx, segment_idxs, type, 0.0, 0.0, 0.0)
end

"""
    mutable struct Tether

A collection of segments that are controlled together by a winch.

$(TYPEDFIELDS)
"""
mutable struct Tether
    const idx::Int16
    const segment_idxs::Vector{Int16}
    const winch_idx::Int16
    stretched_len::SimFloat
end

"""
    Tether(idx, segment_idxs, winch_idx)

Constructs a `Tether` object representing a flexible line composed of multiple segments.

A tether enforces a shared unstretched length constraint across all its constituent segments:

**Length Constraint:**
```math
\\sum_{i \\in \\text{segments}} l_{0,i} = L
```

**Winch Control:**
The unstretched tether length `L` is controlled by a winch.

# Arguments
- `idx::Int16`: Unique identifier for the tether.
- `segment_idxs::Vector{Int16}`: Indices of segments that form this tether.
- `winch_idx::Int16`: Index of the winch controlling this tether.

# Returns
- `Tether`: A new `Tether` object.
"""
function Tether(idx, segment_idxs, winch_idx)
    return Tether(idx, segment_idxs, winch_idx, 0.0)
end

"""
    mutable struct Winch

A set of tethers (or a single tether) connected to a winch mechanism.

$(TYPEDFIELDS)
"""
mutable struct Winch
    const idx::Int16
    const model::AbstractWinchModel
    const tether_idxs::Vector{Int16}
    tether_len::Union{SimFloat, Nothing}
    tether_vel::SimFloat
    set_value::SimFloat
    brake::Bool
    const force::KVec3
end

"""
    Winch(idx, model, tether_idxs; tether_len=0.0, tether_vel=0.0, brake=false)

Constructs a `Winch` object that controls tether length through torque or speed regulation.

The winch acceleration function `α` depends on the winch model type:
- **Torque-controlled**: Direct torque input with motor dynamics.
- **Speed-controlled**: Velocity regulation with internal control loops.

For detailed mathematical models of winch dynamics, motor characteristics, and control algorithms,
see the [WinchModels.jl documentation](https://github.com/aenarete/WinchModels.jl/blob/main/docs/winch.md).

# Arguments
- `idx::Int16`: Unique identifier for the winch.
- `model::AbstractWinchModel`: The winch model (`TorqueControlledMachine`, etc.).
- `tether_idxs::Vector{Int16}`: Vector of indices of the tethers connected to this winch.

# Keyword Arguments
- `tether_len::SimFloat=0.0`: Initial tether length [m].
- `tether_vel::SimFloat=0.0`: Initial tether velocity (reel-out rate) [m/s].
- `brake::Bool=false`: If true, the winch brake is engaged.

# Returns
- `Winch`: A new `Winch` object.
"""
function Winch(idx, model, tether_idxs; tether_len=0.0, tether_vel=0.0, brake=false)
    return Winch(idx, model, tether_idxs, tether_len, tether_vel, 0.0, brake, zeros(KVec3))
end

"""
    mutable struct Wing

A rigid wing body that can have multiple groups of points attached to it.

The wing provides a rigid body reference frame for attached points and groups.
Points with `type == WING` move rigidly with the wing body according to the
wing's orientation matrix `R_b_w` and position `pos_w`.

# Special Properties
The wing's orientation can be accessed as a rotation matrix or a quaternion:
```julia
R_matrix = wing.R_b_w
wing.R_b_w = R_matrix

quat = wing.Q_b_w
wing.Q_b_w = quat
```

$(TYPEDFIELDS)
"""
mutable struct Wing
    const idx::Int16

    # VSM aerodynamics
    vsm_aero::VortexStepMethod.BodyAerodynamics
    vsm_wing::VortexStepMethod.RamAirWing
    vsm_solver::VortexStepMethod.Solver

    # Structural information
    const group_idxs::Vector{Int16}
    const transform_idx::Int16
    const R_b_c::Matrix{SimFloat}
    const pos_cad::KVec3

    # Differential variables in world frame, updated during simulation
    const Q_b_w::Vector{SimFloat}
    const ω_b::KVec3
    const pos_w::KVec3
    const vel_w::KVec3
    const acc_w::KVec3

    # Derived variables and parameters, updated during simulation
    const vsm_y::Vector{SimFloat}
    const vsm_x::Vector{SimFloat}
    const vsm_jac::Matrix{SimFloat}
    const wind_disturb::KVec3
    const va_b::KVec3 # apparent wind in body frame
    const v_wind::KVec3 # wind velocity in world frame at the wing
    const aero_force_b::KVec3 # aerodynamic force in body frame
    const aero_moment_b::KVec3 # aerodynamic moment in body frame
    elevation::SimFloat
    elevation_vel::SimFloat
    elevation_acc::SimFloat
    azimuth::SimFloat
    azimuth_vel::SimFloat
    azimuth_acc::SimFloat
    heading::SimFloat
    const turn_rate::KVec3
    const turn_acc::KVec3
    course::SimFloat
    aoa::SimFloat
    fix_sphere::Bool
end
function Base.getproperty(wing::Wing, sym::Symbol)
    if sym == :R_b_w
        return quaternion_to_rotation_matrix(wing.Q_b_w)
    else
        return getfield(wing, sym)
    end
end
function Base.setproperty!(wing::Wing, sym::Symbol, value)
    if sym == :R_b_w
        wing.Q_b_w .= rotation_matrix_to_quaternion(value)
    else
        setfield!(wing, sym, value)
    end
end

"""
    Wing(idx, vsm_aero, vsm_wing, vsm_solver, group_idxs, R_b_c, pos_cad; transform_idx)

Constructs a `Wing` object representing a rigid body that serves as a reference frame.

A `Wing` provides a rigid body coordinate system for kite components. Points with `type == WING`
move rigidly with the wing body. Groups attached to the wing undergo local deformation
(twist) relative to the rigid wing body frame.

**Rigid Body Dynamics:**
The wing follows standard rigid body equations of motion:
```math
\\begin{aligned}
\\frac{\\delta \\mathbf{q}_b^w}{\\delta t} &= \\frac{1}{2} \\Omega(\\boldsymbol{\\omega}_b) \\mathbf{q}_b^w \\\\
\\boldsymbol{\\tau}_w &= \\mathbf{I} \\frac{\\delta \\boldsymbol{\\omega}}{\\delta t} + \\boldsymbol{\\omega}_b \\times (\\mathbf{I}\\boldsymbol{\\omega}_b)
\\end{aligned}
```
where ``\\mathbf{q}_b^w`` is the quaternion, ``\\boldsymbol{\\omega}_b`` is the angular velocity,
``\\mathbf{I}`` is the inertia tensor, and ``\\boldsymbol{\\tau}_w`` is the total applied torque.

**Coordinate Transformations:**
Points attached to the wing transform as: ``\\mathbf{r}_w = \\mathbf{r}_{wing} + \\mathbf{R}_{b \\rightarrow w} \\mathbf{r}_b``

# Arguments
- `idx::Int16`: Unique identifier for the wing.
- `vsm_aero`, `vsm_wing`, `vsm_solver`: Vortex Step Method components.
- `group_idxs::Vector{Int16}`: Indices of groups attached to this wing.
- `R_b_c::Matrix{SimFloat}`: Rotation matrix from body frame to CAD frame.
- `pos_cad::KVec3`: Position of wing center of mass in CAD frame.

# Keyword Arguments
- `transform_idx::Int16=1`: Transform used for initial positioning and orientation.

# Returns
- `Wing`: A new `Wing` object providing a rigid body reference frame.
"""
function Wing(idx, vsm_aero, vsm_wing, vsm_solver, group_idxs, R_b_c, pos_cad; 
    transform_idx=1
)
    ny = length(group_idxs)+3+3
    nx = length(group_idxs)+3+3
    return Wing(idx, 
        vsm_aero, vsm_wing, vsm_solver,
        # Structural information
        group_idxs, transform_idx, R_b_c, pos_cad, 
        # Differential variables in world frame, updated during simulation
        zeros(4), zeros(KVec3),
        zeros(KVec3), zeros(KVec3), zeros(KVec3),
        # Derived variables and parameters, updated during simulation
        zeros(SimFloat, ny), zeros(SimFloat, nx), zeros(SimFloat, nx, ny),
        zeros(KVec3), zeros(KVec3), zeros(KVec3), zeros(KVec3), zeros(KVec3),
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
        zeros(KVec3), zeros(KVec3), 0.0, 0.0, false)
end

"""
    mutable struct Transform

Describes the spatial transformation (position and orientation) of system components
relative to a base reference point.

$(TYPEDFIELDS)
"""
mutable struct Transform
    const idx::Int16
    const wing_idx::Union{Int16, Nothing}
    const rot_point_idx::Union{Int16, Nothing}
    const base_point_idx::Union{Int16, Nothing}
    const base_transform_idx::Union{Int16, Nothing}
    elevation::SimFloat # The elevation of the rotating point or kite as seen from the base point
    azimuth::SimFloat # The azimuth of the rotating point or kite as seen from the base point
    heading::SimFloat
    base_pos::Union{KVec3, Nothing}
end

"""
    Transform(idx, elevation, azimuth, heading; base_point_idx, base_pos, base_transform_idx, wing_idx, rot_point_idx)

Constructs a `Transform` object that orients system components using spherical coordinates.

All points and wings with a matching `transform_idx` are transformed together as a rigid body:
1. **Translation**: Translate such that `base_point_idx` is at the specified `base_pos`.
2. **Rotation 1**: Rotate so the target (`wing_idx` or `rot_point_idx`) is at (`elevation`, `azimuth`) relative to the base.
3. **Rotation 2**: Rotate all components by `heading` around the base-target vector.

```math
\\mathbf{r}_{transformed} = \\mathbf{r}_{base} + \\mathbf{R}_{heading} \\circ \\mathbf{R}_{elevation,azimuth}(\\mathbf{r} - \\mathbf{r}_{base})
```

# Arguments
- `idx::Int16`: Unique identifier for the transform.
- `elevation::SimFloat`: Target elevation angle from base [rad].
- `azimuth::SimFloat`: Target azimuth angle from base [rad].
- `heading::SimFloat`: Rotation around base-target vector [rad].

# Keyword Arguments
**Base Reference (choose one method):**
- `base_pos` & `base_point_idx`: Use a fixed position and a reference point.
- `base_transform_idx`: Chain to another transform's position.

**Target Object (choose one):**
- `wing_idx`: The wing to position at (`elevation`, `azimuth`).
- `rot_point_idx`: The point to position at (`elevation`, `azimuth`).

# Returns
- `Transform`: A transform affecting all components with a matching `transform_idx`.
"""
function Transform(idx, elevation, azimuth, heading;
        base_point_idx=nothing, base_pos=nothing, base_transform_idx=nothing,
        wing_idx=nothing, rot_point_idx=nothing)
    (isnothing(wing_idx) == isnothing(rot_point_idx)) && error("Either provide a wing_idx or a rot_point_idx, not both or none.")
    (isnothing(base_pos) == isnothing(base_transform_idx)) && error("Either provide the base_pos or the base_transform_idx, not both or none.")
    (isnothing(base_pos) !== isnothing(base_point_idx)) && error("When providing a base_pos, also provide a base_point_idx.")
    Transform(idx, wing_idx, rot_point_idx, base_point_idx, base_transform_idx, elevation, azimuth, heading, base_pos)
end

"""
    Transform(idx, set, base_point_idx; kwargs...)

Constructor helper to create a `Transform` from a `Settings` object.
"""
function Transform(idx, set, base_point_idx; kwargs...)
    Transform(idx, set.elevations[idx], set.azimuths[idx], set.headings[idx], base_point_idx; kwargs...)
end

"""
    get_rot_pos(transform::Transform, wings, points)

Get the position of the rotating object (wing or point) for a given transform.
"""
function get_rot_pos(transform::Transform, wings, points)
    if !isnothing(transform.wing_idx)
        return wings[transform.wing_idx].pos_w
    elseif !isnothing(transform.rot_point_idx)
        return points[transform.rot_point_idx].pos_w
    end
end

"""
    get_base_pos(transform::Transform, wings, points)

Get the base position for a given transform, resolving chained transforms if necessary.
"""
function get_base_pos(transform::Transform, wings, points)
    curr_base_pos = points[transform.base_point_idx].pos_cad
    if !isnothing(transform.base_pos)
        return transform.base_pos, curr_base_pos
    elseif !isnothing(transform.base_transform_idx)
        base_transform = transforms[transform.base_transform_idx]
        return get_rot_pos(base_transform, wings, points), curr_base_pos
    end
end

"""
    struct SystemStructure

A discrete mass-spring-damper representation of a kite system.

This struct holds all components of the physical model, including points, segments,
winches, and wings, forming a complete description of the kite system's structure.

# Components
- [`Point`](@ref): Point masses.
- [`Group`](@ref): Collections of points for wing deformation.
- [`Segment`](@ref): Spring-damper elements.
- [`Pulley`](@ref): Elements that redistribute line lengths.
- [`Tether`](@ref): Groups of segments controlled by a winch.
- [`Winch`](@ref): Ground-based winches.
- [`Wing`](@ref): Rigid wing bodies.
- [`Transform`](@ref): Spatial transformations for initial positioning.
"""
mutable struct SystemStructure
    const name::String
    const set::Settings
    const points::Vector{Point}
    const groups::Vector{Group}
    const segments::Vector{Segment}
    const pulleys::Vector{Pulley}
    const tethers::Vector{Tether}
    const winches::Vector{Winch}
    const wings::Vector{Wing}
    const transforms::Vector{Transform}
    const y::Array{Float64, 2}
    const x::Array{Float64, 2}
    const jac::Array{Float64, 3}
    const wind_vec_gnd::KVec3
    stabilize::Bool
    fix_wing::Bool
end

"""
    SystemStructure(name, set; points, groups, segments, pulleys, tethers, winches, wings, transforms)

Constructs a `SystemStructure` object representing a complete kite system.

## Physical Models
- **"ram"**: A model with 4 deformable wing groups and a complex pulley bridle system.
- **"simple_ram"**: A model with 4 deformable wing groups and direct bridle connections.

# Arguments
- `name::String`: Model identifier ("ram", "simple_ram", or a custom name).
- `set::Settings`: Configuration parameters from `KiteUtils.jl`.

# Keyword Arguments
- `points`, `groups`, `segments`, etc.: Vectors of the system components.

# Returns
- `SystemStructure`: A complete system ready for building a `SymbolicAWEModel`.
"""
function SystemStructure(name, set; 
        points=Point[], 
        groups=Group[], 
        segments=Segment[], 
        pulleys=Pulley[], 
        tethers=Tether[], 
        winches=Winch[], 
        wings=Wing[],
        transforms=Transform[],
    )
    for (i, point) in enumerate(points)
        @assert point.idx == i
        @assert point.transform_idx <= length(transforms)
    end
    for (i, group) in enumerate(groups)
        @assert group.idx == i
    end
    for (i, segment) in enumerate(segments)
        @assert segment.idx == i
        (segment.l0 ≈ 0) && (segment.l0 = norm(points[segment.point_idxs[1]].pos_cad - points[segment.point_idxs[2]].pos_cad))
    end
    for (i, pulley) in enumerate(pulleys)
        @assert pulley.idx == i
    end
    for (i, tether) in enumerate(tethers)
        @assert tether.idx == i
    end
    for (i, winch) in enumerate(winches)
        @assert winch.idx == i
        if iszero(winch.tether_len)
            for segment_idx in tethers[winch.tether_idxs[1]].segment_idxs
                winch.tether_len += segments[segment_idx].l0
            end
        end
    end
    for (i, wing) in enumerate(wings)
        @assert wing.idx == i
    end
    for (i, transform) in enumerate(transforms)
        @assert transform.idx == i
        set.elevations[i] = rad2deg(transform.elevation)
        set.azimuths[i]   = rad2deg(transform.azimuth)
        set.headings[i]   = rad2deg(transform.heading)
    end
    if length(wings) > 0
        ny = 3+length(wings[1].group_idxs)+3
        nx = 3+3+length(wings[1].group_idxs)
    else
        ny = 0
        nx = 0
    end
    y = zeros(length(wings), ny)
    x = zeros(length(wings), nx)
    jac = zeros(length(wings), nx, ny)
    set.physical_model = name
    sys_struct = SystemStructure(name, set, points, groups, segments, pulleys, tethers,
        winches, wings, transforms, y, x, jac, zeros(KVec3), false, false)
    reinit!(sys_struct, set)
    return sys_struct
end

"""
    apply_heading(vec, R_t_w, curr_R_t_w, heading)

Apply a heading rotation to a vector.
"""
function apply_heading(vec, R_t_w, curr_R_t_w, heading)
    vec_along_z = rotate_around_z(curr_R_t_w' * vec, heading)
    return R_t_w * vec_along_z
end

"""
    reinit!(transforms::Vector{Transform}, sys_struct::SystemStructure)

Apply the initial spatial transformations to all components in a `SystemStructure`.

This function iterates through all transforms and applies the specified translation
and rotation to position and orient the kite system components correctly in the
world frame at the beginning of a simulation.
"""
function reinit!(transforms::Vector{Transform}, sys_struct::SystemStructure)
    @unpack points, wings = sys_struct
    for transform in transforms
        # ==================== TRANSLATE ==================== #
        base_pos, curr_base_pos = get_base_pos(transform, wings, points)
        T = base_pos - curr_base_pos
        for point in points
            if point.transform_idx == transform.idx
                point.pos_w .= point.pos_cad .+ T
                point.vel_w .= 0.0
            end
        end
        for wing in wings
            if wing.transform_idx == transform.idx
                wing.pos_w .= wing.pos_cad .+ T
                wing.vel_w .= 0.0
            end
        end

        # ==================== ROTATE ==================== #
        curr_rot_pos = get_rot_pos(transform, wings, points)
        curr_elevation = KiteUtils.calc_elevation(curr_rot_pos - base_pos)
        curr_azimuth = -KiteUtils.azimuth_east(curr_rot_pos - base_pos)
        curr_R_t_w = calc_R_t_w(curr_elevation, curr_azimuth)
        R_t_w = calc_R_t_w(transform.elevation, transform.azimuth)

        for point in points
            if point.transform_idx == transform.idx
                vec = point.pos_w - base_pos
                point.pos_w .= base_pos + apply_heading(vec, R_t_w, curr_R_t_w, transform.heading)
            end
            if point.type == WING
                wing = wings[point.wing_idx]
                point.pos_b .= wing.R_b_c' * (point.pos_cad - wing.pos_cad) # TODO: test this
            end
        end
        for wing in wings
            if wing.transform_idx == transform.idx
                vec = wing.pos_w - base_pos
                wing.pos_w .= base_pos + apply_heading(vec, R_t_w, curr_R_t_w, transform.heading)
                R_b_w = zeros(3,3)
                for i in 1:3
                    R_b_w[:, i] .= apply_heading(wing.R_b_c[:, i], R_t_w, curr_R_t_w, transform.heading)
                end
                wing.R_b_w = R_b_w
            end
        end
    end
end

"""
    calc_pos(wing::RamAirWing, gamma, frac)

Calculate a position on the kite based on spanwise (`gamma`) and chordwise (`frac`) parameters.
"""
function calc_pos(wing::RamAirWing, gamma, frac)
    le_pos = [wing.le_interp[i](gamma) for i in 1:3]
    chord = [wing.te_interp[i](gamma) for i in 1:3] .- le_pos
    pos = le_pos .+ chord .* frac
    return pos
end

"""
    create_tether(tether_idx, set, points, segments, tethers, attach_point, type, dynamics_type; z, axial_stiffness, axial_damping)

Procedurally create a multi-segment tether.

This function builds a tether from a specified number of segments, connecting a given
`attach_point` on the kite to a new anchor point on the ground.
"""
function create_tether(tether_idx, set, points, segments, tethers, attach_point, 
                       type, dynamics_type; z=[0,0,1], axial_stiffness=NaN, 
                       axial_damping=NaN)
    winch_pos = find_axis_point(attach_point.pos_cad, set.l_tether, z)
    dir = winch_pos - attach_point.pos_cad
    segment_idxs = Int16[]
    winch_idx = 0
    for i in 1:set.segments
        frac = i / set.segments
        pos = attach_point.pos_cad + frac * dir
        point_idx = length(points)+1 # last point idx
        segment_idx = length(segments)+1 # last segment idx
        if i == 1
            last_idx = attach_point.idx
        else
            last_idx = point_idx-1
        end
        if i == set.segments
            points = [points; Point(point_idx, pos, STATIC)]
            winch_idx = points[end].idx
        else
            points = [points; Point(point_idx, pos, dynamics_type)]
        end
        segments = [segments; Segment(segment_idx, set, (last_idx, point_idx), type;
                                      axial_stiffness, axial_damping)]
        push!(segment_idxs, segment_idx)
    end
    tethers = [tethers; Tether(tether_idx, segment_idxs, winch_idx)]
    return points, segments, tethers, tethers[end].idx
end

"""
    cad_to_body_frame(wing::RamAirWing, pos)

Transform a position from the CAD frame to the wing's body frame.
"""
function cad_to_body_frame(wing::RamAirWing, pos)
    return wing.R_cad_body * (pos + wing.T_cad_body)
end

"""
    find_axis_point(P, l, v=[0,0,1])

Calculate the coordinates of a point `Q` that lies on a line defined by vector `v`
and is at a distance `l` from a given point `P`.
"""
function find_axis_point(P, l, v=[0,0,1])
    # Compute discriminant
    D = (v ⋅ P)^2 - norm(v)^2 * (norm(P)^2 - l^2)
    D < 0 && error("No real solution: l is too small or parameters invalid")
    # Solve quadratic for t, choose solution for negative direction
    t = (v ⋅ P - √D) / norm(v)^2
    # Compute point Q = t * v
    return [t * v[1], t * v[2], t * v[3]]
end

"""
    reinit!(sys_struct::SystemStructure, set::Settings)

Re-initialize a `SystemStructure` from a `Settings` object.

This function resets various component states (e.g., winch lengths, group twists,
pulley positions) to their initial values as defined in the `Settings` object. It
is typically called before starting a new simulation run.
"""
function reinit!(sys_struct::SystemStructure, set::Settings)
    @unpack points, groups, segments, pulleys, tethers, winches, wings, transforms = sys_struct

    for segment in segments
        @assert (0 < segment.diameter < 1)
    end

    for winch in winches
        winch.tether_len = set.l_tethers[winch.idx]
        winch.tether_vel    = set.v_reel_outs[winch.idx]
    end

    (length(groups) > 0) && (first_moment_frac = groups[1].moment_frac)
    for group in groups
        group.twist = 0.0
        group.twist_ω = 0.0
        @assert group.moment_frac ≈ first_moment_frac "All group.moment_frac must be the same."
    end
    
    for transform in transforms
        transform.elevation = deg2rad(set.elevations[transform.idx])
        transform.azimuth   = deg2rad(set.azimuths[transform.idx])
        transform.heading   = deg2rad(set.headings[transform.idx])
    end

    for segment in segments
        (segment.l0 ≈ 0) && (segment.l0 = norm(points[segment.point_idxs[1]].pos_cad - points[segment.point_idxs[2]].pos_cad))
        @assert (segment.l0 > 0)
    end

    for pulley in pulleys
        segment1, segment2 = segments[pulley.segment_idxs[1]], segments[pulley.segment_idxs[2]]
        pulley.sum_len = segment1.l0 + segment2.l0
        pulley.len = segment1.l0
        pulley.vel = 0.0
        @assert !(pulley.sum_len ≈ 0)
    end

    reinit!(transforms, sys_struct)
    for wing in wings
        wing.vsm_y .= 0.0
        wing.vsm_y[1:3] .= wing.R_b_w' * [set.v_wind, 0., 0.]
    end

    return nothing
end

"""
    copy!(sys1::SystemStructure, sys2::SystemStructure)

Copy the dynamic state from one `SystemStructure` (`sys1`) to another (`sys2`).

This function is designed to transfer the state (positions, velocities, etc.) between
two system models, which can have different levels of fidelity. For example, it can
copy the state from a detailed multi-segment tether model (`sys1`) to a simplified
single-segment model (`sys2`).

The function handles several cases:
- If `sys1` and `sys2` have the same structure, it performs a direct copy of all point states.
- If `sys2` is a simplified (1-segment per tether) version of `sys1`, it copies the
  positions and velocities of the tether endpoints.
- It also copies the state of wings, groups, winches, and pulleys where applicable.
"""
function copy!(sys1::SystemStructure, sys2::SystemStructure)
    simple = false

    # copy point pos and vel
    if length(sys1.points) > 0
        if length(sys1.points) == length(sys2.points)
            for (point1, point2) in zip(sys1.points, sys2.points)
                point2.pos_w .= point1.pos_w
                point2.vel_w .= point1.vel_w
                point2.disturb .= point1.disturb
            end
        # if different number of points, copy only the tether points
        elseif length(sys1.tethers) > 1 && length(sys1.tethers) == length(sys2.tethers)
            for (tether1, tether2) in zip(sys1.tethers, sys2.tethers)
                if length(tether1.segment_idxs) == length(tether2.segment_idxs)
                    # copy the points of the segments of the tethers
                    for (segment_idx1, segment_idx2) in zip(tether1.segment_idxs, tether2.segment_idxs)
                        point_idxs1 = sys1.segments[segment_idx1].point_idxs
                        point_idxs2 = sys2.segments[segment_idx2].point_idxs
                        for (point_idx1, point_idx2) in zip(point_idxs1, point_idxs2)
                            sys2.points[point_idx2].pos_w .= sys1.points[point_idx1].pos_w
                            sys2.points[point_idx2].vel_w .= sys1.points[point_idx1].vel_w
                            sys2.points[point_idx2].disturb .= sys1.points[point_idx1].disturb
                        end
                    end
                elseif length(tether2.segment_idxs) == 1
                    # copy the first and last point of the tether
                    point_idxs1 = [sys1.segments[tether1.segment_idxs[1]].point_idxs[1],
                                   sys1.segments[tether1.segment_idxs[end]].point_idxs[2]]
                    point_idxs2 = sys2.segments[tether2.segment_idxs[1]].point_idxs
                    sys2.points[point_idxs2[1]].pos_w .= sys1.points[point_idxs1[1]].pos_w
                    sys2.points[point_idxs2[2]].pos_w .= sys1.points[point_idxs1[2]].pos_w
                    sys2.points[point_idxs2[1]].vel_w .= sys1.points[point_idxs1[1]].vel_w
                    sys2.points[point_idxs2[2]].vel_w .= sys1.points[point_idxs1[2]].vel_w
                    sys2.points[point_idxs2[1]].disturb .= sys1.points[point_idxs1[1]].disturb
                    sys2.points[point_idxs2[2]].disturb .= sys1.points[point_idxs1[2]].disturb
                    simple = true
                end
            end
        end
    end

    # copy twist and twist_ω of groups
    if length(sys1.groups) > 1 && length(sys1.groups) == length(sys2.groups)
        for (group1, group2) in zip(sys1.groups, sys2.groups)
            group2.twist = group1.twist
            group2.twist_ω = group1.twist_ω
        end
    end

    # copy winch tether lengths and velocities
    if length(sys1.winches) > 1 && length(sys1.winches) == length(sys2.winches)
        for (winch2, winch1) in zip(sys2.winches, sys1.winches)
            if !simple
                winch2.tether_len = winch1.tether_len
                winch2.tether_vel = winch1.tether_vel
            else
                winch2.tether_len = 0.0
                for tether_idx in winch1.tether_idxs
                    tether2 = sys2.tethers[tether_idx]
                    segment2 = sys2.segments[tether2.segment_idxs[1]]
                    point_idxs2 = segment2.point_idxs
                    slen = norm(sys2.points[point_idxs2[1]].pos_w .-
                                        sys2.points[point_idxs2[2]].pos_w)
                    stiffness = segment2.axial_stiffness / slen
                    nt = length(winch1.tether_idxs)
                    winch2.tether_len += (slen - norm(winch1.force)/stiffness/nt) / nt
                end
                winch2.tether_vel = winch1.tether_vel
            end
        end
    end

    # copy pulley lengths and velocities
    if length(sys1.pulleys) > 1 && length(sys1.pulleys) == length(sys2.pulleys)
        for (pulley1, pulley2) in zip(sys1.pulleys, sys2.pulleys)
            pulley2.len = pulley1.len
            pulley2.vel = pulley1.vel
        end
    end

    # copy wing positions and velocities
    if length(sys1.wings) > 1 && length(sys1.wings) == length(sys2.wings)
        for (wing1, wing2) in zip(sys1.wings, sys2.wings)
            wing2.pos_w .= wing1.pos_w
            wing2.vel_w .= wing1.vel_w
            wing2.ω_b .= wing1.ω_b
            wing2.Q_b_w .= wing1.Q_b_w
        end
    end
end
