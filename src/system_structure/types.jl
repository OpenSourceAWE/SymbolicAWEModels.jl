# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

"""
Basic type definitions for the system structure components.

This file contains enums and struct definitions for:
- SegmentType, DynamicsType, WingType enums
- Point, Group, Segment, Pulley, Tether, Winch structs
"""

# ==================== ENUMS ==================== #

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
    WingType `QUATERNION` `REFINE`

Enumeration for the aerodynamic model type of a wing.

# Elements
- `QUATERNION`: Wing uses quaternion-based rigid body dynamics with twist groups.
  Aerodynamic forces/moments are applied to the wing center of mass.
- `REFINE`: Wing uses refined per-panel forces directly applied to structural points.
  VSM panel forces are lumped to WING-type points with no rigid body constraint.
"""
@enum WingType begin
    QUATERNION
    REFINE
end

"""
    AeroMode `AERO_NONE` `AERO_DIRECT` `AERO_LINEARIZED`

Enumeration for how aerodynamic forces enter the ODE system.
Orthogonal to WingType — determines the aero computation strategy at runtime.

# Elements
- `AERO_NONE`: No aerodynamic forces (returns zeros). For debugging rigid body dynamics.
- `AERO_DIRECT`: Stored forces from nonlinear VSM solve, piecewise-constant between updates.
- `AERO_LINEARIZED`: First-order Taylor expansion using Jacobian from VSM linearization.
"""
@enum AeroMode begin
    AERO_NONE
    AERO_DIRECT
    AERO_LINEARIZED
end

# ==================== POINT ==================== #

"Reference type that can be an integer index or a symbolic name"
const NameRef = Union{Int, Symbol}

"""
    mutable struct Point

A point mass, representing a node in the mass-spring system.

$(TYPEDFIELDS)
"""
mutable struct Point
    idx::Int64  # Assigned by SystemStructure based on vector position
    const name::Union{Int, Symbol, Nothing}  # Name/identifier (Int for backwards compat)
    transform_idx::Int64 # idx of transform (resolved by SystemStructure from transform_ref)
    wing_idx::Int64      # idx of wing (resolved by SystemStructure from wing_ref)
    const transform_ref::Union{Int, Symbol}  # Raw reference to transform (name or idx), 0=no transform
    const wing_ref::Union{Int, Symbol}       # Raw reference to wing (name or idx), 0=no wing
    const pos_cad::KVec3
    const pos_b::KVec3 # pos relative to wing COM in body frame
    const pos_w::KVec3 # pos in world frame
    const vel_w::KVec3 # vel in world frame
    const disturb::KVec3 # disturbing force
    const force::KVec3
    const aero_force_b::KVec3 # aerodynamic force in body frame (for REFINE WING points)
    const va_b::KVec3 # apparent velocity in body frame (for VSM per-point va)
    const type::DynamicsType
    extra_mass::SimFloat      # User-provided mass
    total_mass::SimFloat      # extra_mass + segment weights (computed during simulation)
    body_frame_damping::KVec3
    world_frame_damping::KVec3
    area::SimFloat
    drag_coeff::SimFloat
    fix_sphere::Bool
    fix_static::Bool
end

"""
    Point(name, pos_cad, type; wing=1, transform=1, ...)

Constructs a `Point` object, which can be of four different [`DynamicsType`](@ref)s:
- `STATIC`: The point does not move. ``\\ddot{\\mathbf{r}} = \\mathbf{0}``
- `DYNAMIC`: The point moves according to Newton's second law. ``\\ddot{\\mathbf{r}} = \\mathbf{F}/m``
- `QUASI_STATIC`: The acceleration is constrained to be zero by solving a nonlinear problem. ``\\mathbf{F}/m = \\mathbf{0}``
- `WING`: The point has a static position in the rigid body wing frame. ``\\mathbf{r}_w = \\mathbf{r}_{wing} + \\mathbf{R}_{b\\rightarrow w} \\mathbf{r}_b``

# Arguments
- `name::Union{Int, Symbol}`: Name/identifier for the point (e.g., `:kcu`, `:le_1`, or `1` for legacy).
- `pos_cad::KVec3`: Position of the point in the CAD frame.
- `type::DynamicsType`: Dynamics type of the point (`STATIC`, `DYNAMIC`, etc.).

# Keyword Arguments
- `wing::Union{Int, Symbol}=1`: Reference to the wing (name or index).
- `transform::Union{Int, Symbol}=1`: Reference to the transform (name or index).
- `vel_w::KVec3=zeros(KVec3)`: Initial velocity of the point in world frame.
- `extra_mass::Float64=0.0`: User-provided mass of the point [kg].
- `body_frame_damping::Union{Float64,KVec3}=zeros(KVec3)`: Per-axis damping for body frame.
- `world_frame_damping::Union{Float64,KVec3}=zeros(KVec3)`: Per-axis damping for world frame.
- `fix_sphere::Bool=false`: If true, constrains the point to a sphere.
- `fix_static::Bool=false`: If true, dynamically freezes the point.

# Returns
- `Point`: A new `Point` object. The `idx` field is assigned later by SystemStructure.
"""
function Point(name, pos_cad, type;
    wing=nothing, transform=nothing, vel_w=nothing,
    extra_mass=0.0, body_frame_damping=nothing, world_frame_damping=nothing,
    area=0.0, drag_coeff=0.0,
    fix_sphere=false, fix_static=false
)
    # Handle nothing values - wing defaults to 1, transform 0 means no transform
    wing_ref = isnothing(wing) ? 1 : wing
    transform_ref = isnothing(transform) ? 0 : transform
    vel = isnothing(vel_w) ? zeros(KVec3) : KVec3(vel_w...)

    # Convert scalar damping to vector (broadcast to all axes), handle nothing
    bf_damp = if isnothing(body_frame_damping)
        zeros(KVec3)
    elseif body_frame_damping isa Real
        SVector{3,SimFloat}(body_frame_damping, body_frame_damping, body_frame_damping)
    else
        SVector{3,SimFloat}(body_frame_damping...)
    end
    wf_damp = if isnothing(world_frame_damping)
        zeros(KVec3)
    elseif world_frame_damping isa Real
        SVector{3,SimFloat}(world_frame_damping, world_frame_damping, world_frame_damping)
    else
        SVector{3,SimFloat}(world_frame_damping...)
    end

    # idx, transform_idx, wing_idx are placeholders - resolved by SystemStructure
    Point(0, name, 0, 0, transform_ref, wing_ref, KVec3(pos_cad...), zeros(KVec3), zeros(KVec3),
        vel, zeros(KVec3), zeros(KVec3), zeros(KVec3), zeros(KVec3), type, extra_mass, 0.0,
        bf_damp, wf_damp, area, drag_coeff,
        fix_sphere, fix_static)
end

# ==================== GROUP ==================== #

"""
    mutable struct Group

A set of bridle lines that share the same twist angle and trailing edge angle.

$(TYPEDFIELDS)
"""
mutable struct Group
    idx::Int64  # Assigned by SystemStructure
    const name::Union{Int, Symbol, Nothing}
    point_idxs::Vector{Int64}  # Resolved by SystemStructure from point_refs
    const point_refs::Vector{NameRef}  # Raw references to points (names or indices)
    le_pos::KVec3  # Leading edge position in body frame (from closest VSM panel)
    chord::KVec3   # Chord vector in body frame (from closest VSM panel)
    y_airf::KVec3  # Spanwise vector in local panel frame (from closest VSM panel)
    const type::DynamicsType
    moment_frac::SimFloat
    damping::SimFloat
    twist::SimFloat
    twist_ω::SimFloat
    tether_force::SimFloat
    tether_moment::SimFloat
    aero_moment::SimFloat
    unrefined_section_idxs::Vector{Int64}  # Indices of VSM unrefined sections in this group
end

"""
    Group(name, points, type, moment_frac; damping=50.0)

Constructs a `Group` object representing a collection of points on a
kite body that share a common twist deformation.

Group geometry (le_pos, chord, y_airf) is computed later by SystemStructure
using the closest VSM panel to the group's mean point position.

# Arguments
- `name::Union{Int, Symbol}`: Name/identifier for the group.
- `points::Vector`: References to points (names or indices).
- `type::DynamicsType`: DYNAMIC or QUASI_STATIC.
- `moment_frac::SimFloat`: Chordwise rotation point (0=LE, 1=TE).

# Keyword Arguments
- `damping::SimFloat=50.0`: Damping coefficient for twist dynamics.

# Returns
- `Group`: A new `Group` object. The `idx` and `point_idxs` are resolved by SystemStructure.
  Geometry fields (le_pos, chord, y_airf) are initialized to zero and computed during
  SystemStructure construction from the closest VSM panel.
"""
function Group(name, points, type, moment_frac; damping=50.0)
    point_refs = Vector{NameRef}([p isa Integer ? Int(p) : Symbol(p) for p in points])
    Group(0, name, Int64[], point_refs,
          zeros(KVec3), zeros(KVec3), zeros(KVec3),
          type, moment_frac, damping,
          0.0, 0.0, 0.0, 0.0, 0.0,
          Int64[])
end

# ==================== SEGMENT ==================== #

"""
    mutable struct Segment

A segment representing a spring-damper connection from one point to another.

The spring-damper model uses per-unit-length stiffness and damping:
- Effective stiffness: `k = unit_stiffness / length` [N/m]
- Effective damping: `c = unit_damping / length` [N·s/m]

$(TYPEDFIELDS)
"""
mutable struct Segment
    idx::Int64  # Assigned by SystemStructure
    const name::Union{Int, Symbol, Nothing}
    point_idxs::Tuple{Int64, Int64}  # Resolved by SystemStructure from point_refs
    const point_refs::Tuple{NameRef, NameRef}  # Raw references to endpoints (names or indices)
    "Stiffness per unit length [N]. Effective k = unit_stiffness/length [N/m]."
    unit_stiffness::SimFloat
    "Damping per unit length [N·s]. Effective c = unit_damping/length [N·s/m]."
    unit_damping::SimFloat
    "Rest (unstretched) length [m]."
    l0::SimFloat
    "Stiffness reduction factor when segment is in compression (0-1)."
    compression_frac::SimFloat
    "Segment diameter [m]."
    diameter::SimFloat
    "Current length of the segment [m]."
    len::SimFloat
    "Current force in the segment [N]."
    force::SimFloat
end

"""
    Segment(name, point_i, point_j, unit_stiffness, unit_damping, diameter; l0, compression_frac)

Basic constructor for a `Segment` object.

# Arguments
- `name::Union{Int, Symbol}`: Name/identifier for the segment.
- `point_i`, `point_j`: References to the two endpoint points (names or indices).
- `unit_stiffness`: Stiffness per unit length [N]. Effective k = unit_stiffness/length [N/m].
- `unit_damping`: Damping per unit length [N·s]. Effective c = unit_damping/length [N·s/m].
- `diameter`: Segment diameter [m].
"""
function Segment(name, point_i, point_j, unit_stiffness, unit_damping, diameter;
    l0=zero(SimFloat), compression_frac=0.1
)
    p1 = point_i isa Integer ? Int(point_i) : Symbol(point_i)
    p2 = point_j isa Integer ? Int(point_j) : Symbol(point_j)
    Segment(0, name, (0, 0), (p1, p2), unit_stiffness, unit_damping, l0, compression_frac,
        diameter, zero(SimFloat), zero(SimFloat))
end

"""
    Segment(name, set, point_i, point_j, type; l0, compression_frac, unit_stiffness, unit_damping)

Constructs a `Segment` using settings for material properties.

# Arguments
- `name::Union{Int, Symbol}`: Name/identifier for the segment.
- `set::Settings`: The settings object containing material properties.
- `point_i`, `point_j`: References to the two endpoint points (names or indices).
- `type::SegmentType`: Type of the segment (`POWER_LINE`, `STEERING_LINE`, `BRIDLE`).

# Keyword Arguments
- `l0::SimFloat=zero(SimFloat)`: Unstretched length [m]. Calculated from point positions if zero.
- `compression_frac::SimFloat=0.0`: Stiffness reduction factor in compression.
- `diameter_mm::Float64=NaN`: Tether diameter [mm]. If `NaN`, uses default from settings.
- `unit_stiffness::Float64=NaN`: Stiffness per unit length [N]. Effective k = unit_stiffness/length.
- `unit_damping::Float64=NaN`: Damping per unit length [N·s]. Effective c = unit_damping/length.
"""
function Segment(name, set, point_i, point_j, type;
    l0=zero(SimFloat), compression_frac=0.0, diameter_mm=NaN, unit_stiffness=NaN, unit_damping=NaN
)
    p1 = point_i isa Integer ? Int(point_i) : Symbol(point_i)
    p2 = point_j isa Integer ? Int(point_j) : Symbol(point_j)

    # Set default diameter from settings if not specified
    if isnan(diameter_mm)
        (type == BRIDLE) && (diameter_mm = set.bridle_tether_diameter)
        (type == POWER_LINE) && (diameter_mm = set.power_tether_diameter)
        (type == STEERING_LINE) && (diameter_mm = set.steering_tether_diameter)
    end
    # Convert diameter from mm to m
    diameter_m = 0.001 * diameter_mm

    # Compute unit_stiffness if not provided
    if isnan(unit_stiffness)
        unit_stiffness = set.e_tether * (diameter_m/2)^2 * π
        if type == BRIDLE
            stiffness_frac = 0.01
        else
            stiffness_frac = 1.0
        end
        unit_stiffness *= stiffness_frac
    end

    # Compute unit_damping if not provided
    if isnan(unit_damping)
        # Use rel_damping if available, otherwise compute from unit_damping/unit_stiffness ratio
        if hasproperty(set, :rel_damping) && set.rel_damping != 0.0
            unit_damping = set.rel_damping * unit_stiffness
        elseif hasproperty(set, :unit_damping) && hasproperty(set, :unit_stiffness) &&
                set.unit_damping != 0.0
            unit_damping = (set.unit_damping / set.unit_stiffness) * unit_stiffness
        else
            @warn "damping is zero!"
            unit_damping = 0.0  # fallback if no damping info available
        end
    end

    Segment(0, name, (0, 0), (p1, p2), unit_stiffness, unit_damping, l0, compression_frac,
        diameter_m, zero(SimFloat), zero(SimFloat))
end

# ==================== PULLEY ==================== #

"""
    mutable struct Pulley

A pulley described by two segments with the common point of the segments being the pulley.

$(TYPEDFIELDS)
"""
mutable struct Pulley
    idx::Int64  # Assigned by SystemStructure
    const name::Union{Int, Symbol, Nothing}
    segment_idxs::Tuple{Int64, Int64}  # Resolved by SystemStructure from segment_refs
    const segment_refs::Tuple{NameRef, NameRef}  # Raw references to segments
    const type::DynamicsType
    sum_len::SimFloat
    len::SimFloat
    vel::SimFloat
end

"""
    Pulley(name, segment_i, segment_j, type)

Constructs a `Pulley` object that enforces length redistribution between two segments.

# Arguments
- `name::Union{Int, Symbol}`: Name/identifier for the pulley.
- `segment_i`, `segment_j`: References to the two segments (names or indices).
- `type::DynamicsType`: Dynamics type (`DYNAMIC` or `QUASI_STATIC`).
"""
function Pulley(name, segment_i, segment_j, type)
    s1 = segment_i isa Integer ? Int(segment_i) : Symbol(segment_i)
    s2 = segment_j isa Integer ? Int(segment_j) : Symbol(segment_j)
    return Pulley(0, name, (0, 0), (s1, s2), type, 0.0, 0.0, 0.0)
end

# ==================== TETHER ==================== #

"""
    mutable struct Tether

A collection of segments that are controlled together by a winch.

$(TYPEDFIELDS)
"""
mutable struct Tether
    idx::Int64  # Assigned by SystemStructure
    const name::Union{Int, Symbol, Nothing}
    segment_idxs::Vector{Int64}  # Resolved by SystemStructure from segment_refs
    const segment_refs::Vector{NameRef}  # Raw references to segments
    winch_point_idx::Int64  # Resolved by SystemStructure from winch_point_ref
    const winch_point_ref::NameRef  # Raw reference to winch point
    stretched_len::SimFloat
end

"""
    Tether(name, segments; winch_point=nothing)

Constructs a `Tether` object representing a flexible line composed of multiple segments.

# Arguments
- `name::Union{Int, Symbol}`: Name/identifier for the tether.
- `segments::Vector`: References to segments that form this tether (names or indices).

# Keyword Arguments
- `winch_point=nothing`: Reference to the ground point where tether attaches to winch.
  Defaults to point 1 if not specified.
"""
function Tether(name, segments; winch_point=nothing)
    segment_refs = Vector{NameRef}([s isa Integer ? Int(s) : Symbol(s) for s in segments])
    # Handle nothing - default to point 1
    wp = isnothing(winch_point) ? 1 : winch_point
    wp_ref = wp isa Integer ? Int(wp) : Symbol(wp)
    return Tether(0, name, Int64[], segment_refs, 0, wp_ref, 0.0)
end

# ==================== WINCH ==================== #

"""
    mutable struct Winch

A set of tethers (or a single tether) connected to a winch mechanism.

$(TYPEDFIELDS)
"""
mutable struct Winch
    idx::Int64  # Assigned by SystemStructure
    const name::Union{Int, Symbol, Nothing}
    tether_idxs::Vector{Int64}  # Resolved by SystemStructure from tether_refs
    const tether_refs::Vector{NameRef}  # Raw references to tethers
    tether_len::Union{SimFloat, Nothing}
    tether_vel::SimFloat
    tether_acc::SimFloat
    set_value::SimFloat
    brake::Bool
    speed_controlled::Bool
    const force::KVec3
    gear_ratio::SimFloat
    drum_radius::SimFloat
    f_coulomb::SimFloat
    c_vf::SimFloat
    inertia_total::SimFloat
    friction::SimFloat
    friction_epsilon::SimFloat
end

"""
    Winch(name, set, tethers; tether_len=0.0, tether_vel=0.0, brake=false, speed_controlled=false)

Constructs a `Winch` object that controls tether length through torque or speed regulation.

# Arguments
- `name::Union{Int, Symbol}`: Name/identifier for the winch.
- `set::Settings`: The main settings object, used to retrieve winch parameters.
- `tethers::Vector`: References to tethers connected to this winch (names or indices).

# Keyword Arguments
- `tether_len::SimFloat=0.0`: Initial tether length [m].
- `tether_vel::SimFloat=0.0`: Initial tether velocity (reel-out rate) [m/s].
- `brake::Bool=false`: If true, the winch brake is engaged.
- `speed_controlled::Bool=false`: If true, tether velocity is prescribed externally
  (not integrated by the ODE). `D(tether_vel) = 0`, length still tracks velocity.
- `friction_epsilon::SimFloat=6.0`: Smoothing parameter for sign function in Coulomb friction.
"""
function Winch(name, set::Settings, tethers;
               tether_len=0.0, tether_vel=0.0, brake=false,
               speed_controlled=false, friction_epsilon=6.0)
    tether_refs = Vector{NameRef}([t isa Integer ? Int(t) : Symbol(t) for t in tethers])
    return Winch(0, name, Int64[], tether_refs,
                 tether_len, tether_vel, 0.0, 0.0,
                 brake, speed_controlled, zeros(KVec3),
                 set.gear_ratio, set.drum_radius,
                 set.f_coulomb, set.c_vf,
                 set.inertia_total, zero(SimFloat),
                 friction_epsilon)
end

"""
    Winch(name, tethers, gear_ratio, drum_radius, f_coulomb, c_vf, inertia_total; tether_len=0.0, tether_vel=0.0, brake=false, speed_controlled=false)

Constructs a `Winch` object by directly providing its physical parameters.

# Arguments
- `name::Union{Int, Symbol}`: Name/identifier for the winch.
- `tethers::Vector`: References to tethers connected to this winch (names or indices).
- `gear_ratio`, `drum_radius`, `f_coulomb`, `c_vf`, `inertia_total`: Physical parameters.
"""
function Winch(name, tethers, gear_ratio, drum_radius,
               f_coulomb, c_vf, inertia_total;
               tether_len=0.0, tether_vel=0.0, brake=false,
               speed_controlled=false, friction_epsilon=6.0)
    tether_refs = Vector{NameRef}([t isa Integer ? Int(t) : Symbol(t) for t in tethers])
    return Winch(0, name, Int64[], tether_refs,
                 tether_len, tether_vel, 0.0, 0.0,
                 brake, speed_controlled, zeros(KVec3),
                 gear_ratio, drum_radius, f_coulomb,
                 c_vf, inertia_total, zero(SimFloat),
                 friction_epsilon)
end

# ==================== TRANSFORM ==================== #

"""
    mutable struct Transform

Describes the spatial transformation (position and orientation) of system components
relative to a base reference point.

$(TYPEDFIELDS)
"""
mutable struct Transform
    idx::Int64  # Assigned by SystemStructure
    const name::Union{Int, Symbol, Nothing}
    wing_idx::Union{Int64, Nothing}  # Resolved by SystemStructure
    const wing_ref::Union{NameRef, Nothing}  # Raw reference to wing
    rot_point_idx::Union{Int64, Nothing}  # Resolved by SystemStructure
    const rot_point_ref::Union{NameRef, Nothing}  # Raw reference to rotation point
    base_point_idx::Union{Int64, Nothing}  # Resolved by SystemStructure
    const base_point_ref::Union{NameRef, Nothing}  # Raw reference to base point
    base_transform_idx::Union{Int64, Nothing}  # Resolved by SystemStructure
    const base_transform_ref::Union{NameRef, Nothing}  # Raw reference to base transform
    elevation::SimFloat  # [rad]
    azimuth::SimFloat    # [rad]
    heading::SimFloat    # [rad]
    elevation_vel::SimFloat  # [rad/s] angular velocity in elevation direction
    azimuth_vel::SimFloat    # [rad/s] angular velocity in azimuth direction
    turn_rate::SimFloat      # [rad/s] angular velocity around radial axis (not yet implemented)
    base_pos::Union{KVec3, Nothing}
end

# Helper to convert ref to NameRef or nothing
_to_ref(::Nothing) = nothing
_to_ref(x::Integer) = Int(x)
_to_ref(x) = Symbol(x)

"""
    Transform(name, elevation, azimuth, heading; base_point, base_pos, base_transform, wing, rot_point)

Constructs a `Transform` object that orients system components using spherical coordinates.

# Arguments
- `name::Union{Int, Symbol}`: Name/identifier for the transform.
- `elevation`, `azimuth`, `heading`: Spherical coordinates [rad].

# Keyword Arguments
**Base Reference (choose one method):**
- `base_pos` & `base_point`: Use a fixed position and a reference point.
- `base_transform`: Chain to another transform's position.

**Target Object (choose one):**
- `wing`: Reference to the wing to position at (elevation, azimuth).
- `rot_point`: Reference to the point to position at (elevation, azimuth).
"""
function Transform(name, elevation, azimuth, heading;
        base_point=nothing, base_pos=nothing, base_transform=nothing,
        wing=nothing, rot_point=nothing,
        elevation_vel=0.0, azimuth_vel=0.0, turn_rate=0.0)
    (isnothing(wing) == isnothing(rot_point)) && error("Either provide a wing or a rot_point, not both or none.")
    (isnothing(base_pos) == isnothing(base_transform)) && error("Either provide the base_pos or the base_transform, not both or none.")
    (!isnothing(base_pos) && isnothing(base_point)) && error("When providing a base_pos, also provide a base_point.")

    wing_ref = _to_ref(wing)
    rot_point_ref = _to_ref(rot_point)
    base_point_ref = _to_ref(base_point)
    base_transform_ref = _to_ref(base_transform)

    Transform(0, name, nothing, wing_ref, nothing, rot_point_ref,
              nothing, base_point_ref, nothing, base_transform_ref,
              elevation, azimuth, heading, elevation_vel, azimuth_vel, turn_rate,
              isnothing(base_pos) ? nothing : KVec3(base_pos...))
end

"""
    Transform(name, set, base_point; kwargs...)

Constructor helper to create a `Transform` from a `Settings` object.
Note: Uses idx=1 for settings indexing (legacy compatibility).
"""
function Transform(name, set, base_point; idx_for_set=1, kwargs...)
    elevation_vel = hasfield(typeof(set), :elevation_vels) ? set.elevation_vels[idx_for_set] : 0.0
    azimuth_vel = hasfield(typeof(set), :azimuth_vels) ? set.azimuth_vels[idx_for_set] : 0.0
    turn_rate = hasfield(typeof(set), :turn_rates) ? set.turn_rates[idx_for_set] : 0.0
    Transform(name, set.elevations[idx_for_set], set.azimuths[idx_for_set], set.headings[idx_for_set];
              base_point, elevation_vel, azimuth_vel, turn_rate, kwargs...)
end

"""
    get_rot_pos(transform::Transform, wings, points)

Get the position of the rotating object (wing or point) for a given transform.
"""
function get_rot_pos(transform::Transform, wings, points)
    wing_idx = transform.wing_idx
    if !isnothing(wing_idx)
        return wings[something(wing_idx)].pos_w
    end
    rot_point_idx = transform.rot_point_idx
    if !isnothing(rot_point_idx)
        return points[something(rot_point_idx)].pos_w
    end
    error("Transform #$(transform.idx): neither wing_idx nor rot_point_idx is set")
end

"""
    get_base_pos(transform::Transform, transforms, wings, points)

Get the base position for a given transform, resolving chained transforms if necessary.
"""
function get_base_pos(transform::Transform, transforms, wings, points)
    curr_base_pos = points[something(transform.base_point_idx)].pos_cad
    base_pos = transform.base_pos
    if !isnothing(base_pos)
        return something(base_pos), curr_base_pos
    end
    base_transform_idx = transform.base_transform_idx
    if !isnothing(base_transform_idx)
        base_transform = transforms[something(base_transform_idx)]
        return get_rot_pos(base_transform, wings, points), curr_base_pos
    end
    error("Transform #$(transform.idx): neither base_pos nor base_transform_idx is set")
end
