# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

"""
Basic type definitions for the system structure components.

This file contains enums and struct definitions for:
- DynamicsType, WingType, SegmentType (deprecated) enums; AbstractAeroModel types
- Point, TwistSurface, Segment, Pulley, Tether, Winch structs
"""

# ==================== ENUMS ==================== #

"""
    SegmentType `POWER_LINE` `STEERING_LINE` `BRIDLE`

!!! warning "Deprecated"
    `SegmentType` is no longer used as a Segment constructor
    parameter. It is kept only so that old code hits a
    deprecation error instead of `UndefVarError`.
"""
@enum SegmentType begin
    POWER_LINE
    STEERING_LINE
    BRIDLE
end

"""
    DynamicsType `DYNAMIC` `STATIC` `WING` `BODY_STATIC`

Enumeration for the dynamic model governing a point's motion, a rigid body's
motion, or a twist_surface's twist. The shared idea: `DYNAMIC` quantities carry
differential state and are solved by the dynamics; the others are *prescribed* —
no state, held constant within a step, but mutable between steps.

# Elements
- `DYNAMIC`: Solved by the dynamics. A point moves by Newton's second law; a
  rigid body integrates its 6-DOF state; a twist_surface twist solves its
  equilibrium.
- `STATIC`: Prescribed, no state. A point's position is welded in the **world**
  frame; a rigid body is clamped to the world; a twist_surface twist is a
  prescribed control input read live via the registered getter.
- `WING`: A point rides a [`Wing`](@ref) — static in the wing's body frame.
- `BODY_STATIC`: A point rides a [`Body`](@ref) — static in the rigid
  body's body frame; it feeds its net force (and the moment about the body COM)
  into the body.
"""
@enum DynamicsType begin
    DYNAMIC
    STATIC
    WING
    BODY_STATIC
    KINEMATIC
end

"""
    WingType `RIGID_DYNAMICS` `PARTICLE_DYNAMICS`

Enumeration for the aerodynamic model type of a wing.

# Elements
- `RIGID_DYNAMICS`: Wing uses quaternion-based rigid body dynamics with twist twist_surfaces.
  Aerodynamic forces/moments are applied to the wing center of mass.
- `PARTICLE_DYNAMICS`: Wing uses refined per-panel forces directly applied to structural points.
  VSM panel forces are lumped to WING-type points with no rigid body constraint.
"""
@enum WingType begin
    RIGID_DYNAMICS
    PARTICLE_DYNAMICS
end

# Backwards-compatible deprecated aliases for the previous WingType names.
Base.@deprecate_binding QUATERNION RIGID_DYNAMICS
Base.@deprecate_binding REFINE PARTICLE_DYNAMICS

"""
    abstract type WingDynamics

Phantom tag for a wing's dynamics, carried as the second type parameter of
[`Wing`](@ref). It mirrors the `dynamics_type` field but lifts the
rigid/particle distinction into the type domain so aero (and other) methods can
dispatch on [`RigidWing`](@ref)/[`ParticleWing`](@ref) instead of branching on
`dynamics_type`. Subtypes: [`RigidDynamics`](@ref), [`ParticleDynamics`](@ref).
"""
abstract type WingDynamics end

"""Dynamics tag for a `RIGID_DYNAMICS` wing (integrated 6-DOF body)."""
struct RigidDynamics <: WingDynamics end

"""Dynamics tag for a `PARTICLE_DYNAMICS` wing (frame fitted from points)."""
struct ParticleDynamics <: WingDynamics end

"""
    wing_dynamics(dynamics_type::WingType) -> Type{<:WingDynamics}

Map a [`WingType`](@ref) to its [`WingDynamics`](@ref) tag.
"""
wing_dynamics(dynamics_type::WingType) =
    dynamics_type == RIGID_DYNAMICS ? RigidDynamics : ParticleDynamics

"""
    AbstractAeroModel

Supertype for a wing's aerodynamic model. The concrete subtype selects, by
dispatch, the [`aero_component`](@ref) builder that emits the wing's aero
equations. Built-in subtypes: [`AeroNone`](@ref), [`AeroDirect`](@ref),
[`AeroLinearized`](@ref), [`AeroPlate`](@ref). Subtype it and add an
`aero_component` method to plug in custom aerodynamics; see the VSM coupling
documentation for a worked example with a live-updating field.
"""
abstract type AbstractAeroModel end

"""
    mutable struct VSMEngine{BA, W, SL}

Vortex Step Method aerodynamic engine carried by a VSM aero mode
([`AbstractVSMAero`](@ref) `engine` field). Holds the VortexStepMethod objects,
the linearization state, and the structural↔panel mapping.

# Fields
- `vsm_aero`, `vsm_wing`, `vsm_solver`: VortexStepMethod objects.
- `aero_y`: operating-point inputs `[alpha, beta, ω1, ω2, ω3, twist...]`.
- `aero_x`: baseline wind-axis coefficients `[CL, CD, CS, CM1, CM2, CM3, cm...]`.
- `aero_jac`: dense Jacobian `d(aero_x)/d(aero_y)`.
- `point_to_vsm_point`, `wing_segments`: PARTICLE_DYNAMICS structural↔panel maps.
- `aero_scale_chord`: force scale compensating chord-length error (PARTICLE).
- `aero_z_offset`: body-frame z-shift of VSM panels (RIGID).
"""
mutable struct VSMEngine{BA, W, SL}
    vsm_aero::BA
    vsm_wing::W
    vsm_solver::SL
    aero_y::Vector{SimFloat}
    aero_x::Vector{SimFloat}
    aero_jac::Matrix{SimFloat}
    point_to_vsm_point::Union{Nothing, Dict{Int64, Tuple{Int64, Symbol}}}
    wing_segments::Union{Nothing, Vector{Tuple{Int64, Int64}}}
    aero_scale_chord::SimFloat
    aero_z_offset::SimFloat
end

"""
    abstract type AbstractVSMAero <: AbstractAeroModel

Aero modes backed by a [`VSMEngine`](@ref) (stored in the `engine` field). VSM
operations dispatch on this supertype; the engine's fields (`vsm_wing`, `aero_x`,
…) are forwarded from the mode. Built-in subtypes: [`AeroDirect`](@ref) and
[`AeroLinearized`](@ref).
"""
abstract type AbstractVSMAero <: AbstractAeroModel end

# VSM engine fields forwarded from a VSM aero mode to its `engine`.
const VSM_ENGINE_FIELDS = (
    :vsm_aero, :vsm_wing, :vsm_solver, :aero_y, :aero_x, :aero_jac,
    :point_to_vsm_point, :wing_segments, :aero_scale_chord, :aero_z_offset)

function Base.getproperty(mode::AbstractVSMAero, sym::Symbol)
    sym === :engine && return getfield(mode, :engine)
    if sym in VSM_ENGINE_FIELDS
        engine = getfield(mode, :engine)
        engine === nothing && error(
            "$(typeof(mode)) has no VSM engine yet; field `$sym` unavailable.")
        return getproperty(engine, sym)
    end
    return getfield(mode, sym)
end

function Base.setproperty!(mode::AbstractVSMAero, sym::Symbol, value)
    sym === :engine && return setfield!(mode, :engine, value)
    if sym in VSM_ENGINE_FIELDS
        engine = getfield(mode, :engine)
        engine === nothing && error(
            "$(typeof(mode)) has no VSM engine yet; cannot set `$sym`.")
        return setproperty!(engine, sym, value)
    end
    return setfield!(mode, sym, value)
end

# Concrete aero modes and their dispatch methods live in src/aero_modes/.

"""
    is_builtin_aero(mode::AbstractAeroModel) -> Bool

`true` for the package's built-in aero models. Custom models return `false`,
which forces a model rebuild (the compiled cache cannot be reused for
user-supplied equations). Each built-in mode sets this in its `aero_modes/` file.
"""
is_builtin_aero(::AbstractAeroModel) = false

"""
    aero_hash_id(mode::AbstractAeroModel) -> Tuple

Structural fields of `mode` that change the generated equations and therefore
must enter the model-cache key. Return only fields that alter the equation
structure, never runtime-mutable values (those are read live via registered
getters). Defaults to an empty tuple.
"""
aero_hash_id(::AbstractAeroModel) = ()

"""
    NameRef = Union{Int, Symbol}

A reference to another component, either by symbolic name
(`:ground`) or integer index (`1`).

## Name resolution

Components reference each other by name or index at construction
time. These are stored in `_ref` fields (e.g. `point_refs`,
`wing_ref`). During [`SystemStructure`](@ref) construction,
`assign_indices_and_resolve!` maps every ref to a numeric index
via `build_name_dict` (name → vector position) and stores the
result in the corresponding `_idx` fields (e.g. `point_idxs`,
`wing_idx`).

Each component has a `name` field (`const`, set once at
construction) that identifies it for lookup. The type includes
`Nothing` for forward-compatibility but no public constructor
produces `name=nothing`; a nothing-named component would simply
be unreferenceable by name (only by vector index).
"""
const NameRef = Union{Int, Symbol}

# ==================== POINT ==================== #

"""
    mutable struct Point

A point mass, representing a node in the mass-spring system.

$(TYPEDFIELDS)
"""
mutable struct Point
    "Index in the points vector (assigned by SystemStructure)."
    idx::Int64
    "Name used for lookup by other components' `_ref` fields."
    const name::Union{Int, Symbol, Nothing}
    "Resolved transform index (filled by SystemStructure)."
    transform_idx::Int64
    "Resolved wing index (filled by SystemStructure)."
    wing_idx::Int64
    "Resolved rigid-body index for a body-anchored point (filled by SystemStructure). 0 = not anchored."
    body_idx::Int64
    "Raw transform reference (name or idx). 0 = no transform."
    const transform_ref::Union{Int, Symbol}
    "Raw wing reference (name or idx). 0 = no wing."
    const wing_ref::Union{Int, Symbol}
    "Raw rigid-body reference for anchoring (name or idx). 0 = not anchored to a body."
    const body_ref::Union{Int, Symbol}
    "Position in CAD frame [m]."
    const pos_cad::KVec3
    "Position relative to wing COM in principal frame [m]."
    const pos_b::KVec3
    "Anchor offset in the rigid body's body frame [m] (body-anchored points)."
    const anchor_b::KVec3
    "Position in world frame [m] (updated during simulation)."
    const pos_w::KVec3
    "Velocity in world frame [m/s] (updated during simulation)."
    const vel_w::KVec3
    "External force applied to the point, world frame [N] (settable). Was `disturb`."
    const ext_force_w::KVec3
    "Net force on the point [N] (updated during simulation)."
    const force::KVec3
    "Aerodynamic force in body frame [N] (PARTICLE_DYNAMICS WING points)."
    const aero_force_b::KVec3
    "Total drag force in world frame [N], including the point's own drag and any segment drag contributions assigned to it."
    const drag_force::KVec3
    "Apparent velocity in body frame [m/s] (VSM per-point)."
    const va_b::KVec3
    "Dynamics type (DYNAMIC, STATIC, WING, BODY_STATIC)."
    const type::DynamicsType
    "User-provided mass [kg]."
    extra_mass::SimFloat
    "Total mass [kg]: extra_mass + segment contributions (computed during simulation)."
    total_mass::SimFloat
    "Per-axis damping in body frame [N·s/m]."
    body_frame_damping::KVec3
    "Per-axis damping in world frame [N·s/m]."
    world_frame_damping::KVec3
    "Cross-sectional area for drag [m²]."
    area::SimFloat
    "Drag coefficient [-]."
    drag_coeff::SimFloat
    "If true, constrain point to a sphere."
    fix_sphere::Bool
    "If true, dynamically freeze point position."
    fix_static::Bool
end

# Deprecated `point.disturb` forwards to `ext_force_w` (external world-frame force).
Base.getproperty(point::Point, sym::Symbol) =
    getfield(point, sym === :disturb ? :ext_force_w : sym)

"""
    Point(name, pos_cad, type; wing=1, transform=1, ...)

Constructs a `Point` object, which can be of four different [`DynamicsType`](@ref)s:
- `STATIC`: The point does not move. ``\\ddot{\\mathbf{r}} = \\mathbf{0}``
- `DYNAMIC`: The point moves according to Newton's second law. ``\\ddot{\\mathbf{r}} = \\mathbf{F}/m``
- `WING`: The point has a static position in the rigid body wing frame. ``\\mathbf{r}_w = \\mathbf{r}_{wing} + \\mathbf{R}_{b\\rightarrow w} \\mathbf{r}_b``
- `BODY_STATIC`: The point is static in a [`Body`](@ref)'s body frame
  (pass `body`); it rides the body and feeds its net force and moment into it.

# Arguments
- `name::Union{Int, Symbol}`: Name/identifier for the point (e.g., `:kcu`, `:le_1`, or `1` for legacy).
- `pos_cad::KVec3`: Position of the point in the CAD frame.
- `type::DynamicsType`: Dynamics type of the point (`STATIC`, `DYNAMIC`, etc.).
  Pass `BODY_STATIC` together with `body` to anchor the point to a rigid body.

# Keyword Arguments
- `wing::Union{Int, Symbol}=1`: Reference to the wing (name or index).
- `transform::Union{Int, Symbol}=1`: Reference to the transform (name or index).
- `body::Union{Int, Symbol}`: Reference to a [`Body`](@ref) to anchor the
  point to (requires `type = BODY_STATIC`). The point then rides the body
  kinematically and feeds its net force (and the moment about the body COM)
  into the body. Defaults to no anchor.
- `anchor_b::KVec3`: Anchor offset in the body frame [m] (used with `body`).
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
    body=nothing, anchor_b=nothing,
    extra_mass=0.0, body_frame_damping=nothing, world_frame_damping=nothing,
    area=0.0, drag_coeff=0.0,
    fix_sphere=false, fix_static=false
)
    if type == BODY_STATIC
        isnothing(body) && error("Point $name: BODY_STATIC requires a `body` " *
            "reference to a Body.")
    elseif !isnothing(body)
        error("Point $name: `body` is only valid with type BODY_STATIC.")
    end
    # transform 0 means no transform; a body-anchored point has no wing (wing_ref 0).
    body_ref = isnothing(body) ? 0 : body
    wing_ref = isnothing(wing) ? (type == BODY_STATIC ? 0 : 1) : wing
    transform_ref = isnothing(transform) ? 0 : transform
    anchor = isnothing(anchor_b) ? zeros(KVec3) : KVec3(anchor_b...)
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

    # idx, transform_idx, wing_idx, body_idx are placeholders - resolved by SystemStructure
    Point(0, name, 0, 0, 0, transform_ref, wing_ref, body_ref,
        KVec3(pos_cad...), zeros(KVec3), anchor, zeros(KVec3),
        vel, zeros(KVec3), zeros(KVec3), zeros(KVec3), zeros(KVec3), zeros(KVec3),
        type, extra_mass, 0.0,
        bf_damp, wf_damp, area, drag_coeff,
        fix_sphere, fix_static)
end

# ==================== TWIST_SURFACE ==================== #

"""
    mutable struct TwistSurface

A set of bridle lines that share the same twist angle and trailing edge angle.

$(TYPEDFIELDS)
"""
mutable struct TwistSurface
    "Index in the twist_surfaces vector (assigned by SystemStructure)."
    idx::Int64
    "Name used for lookup by other components' `_ref` fields."
    const name::Union{Int, Symbol, Nothing}
    "Resolved point indices (filled by SystemStructure)."
    point_idxs::Vector{Int64}
    "Raw point references (names or indices)."
    const point_refs::Vector{NameRef}
    "Leading edge position in body frame [m] (from closest VSM panel)."
    le_pos::KVec3
    "Chord vector in body frame [m] (from closest VSM panel)."
    chord::KVec3
    "Spanwise vector in local panel frame (from closest VSM panel)."
    y_airf::KVec3
    "Dynamics type (DYNAMIC or STATIC)."
    const type::DynamicsType
    "Chordwise rotation point fraction (0=LE, 1=TE)."
    moment_frac::SimFloat
    "Damping coefficient for twist dynamics [N·m·s/rad]."
    damping::SimFloat
    "Current twist angle [rad]."
    twist::SimFloat
    "Current twist angular velocity [rad/s]."
    twist_ω::SimFloat
    "Tether force contribution [N]."
    tether_force::SimFloat
    "Tether moment contribution [N·m]."
    tether_moment::SimFloat
    "Aerodynamic moment [N·m]."
    aero_moment::SimFloat
    "Indices of VSM unrefined sections in this twist_surface."
    unrefined_section_idxs::Vector{Int64}
    "Surface area [m²] (flat-plate sections; `NaN` when unused)."
    area::SimFloat
end

"""
    TwistSurface(name, points, type, moment_frac; damping=50.0)

Constructs a `TwistSurface` object representing a collection of points on a
kite body that share a common twist deformation.

TwistSurface geometry (le_pos, chord, y_airf) is computed later by SystemStructure
using the closest VSM panel to the twist_surface's mean point position.

# Arguments
- `name::Union{Int, Symbol}`: Name/identifier for the twist_surface.
- `points::Vector`: References to points (names or indices).
- `type::DynamicsType`: DYNAMIC or STATIC.
- `moment_frac::SimFloat`: Chordwise rotation point (0=LE, 1=TE).

# Keyword Arguments
- `damping::SimFloat=50.0`: Damping coefficient for twist dynamics.
- `x_airf=nothing`: Chord-direction reference (body frame). When given, stored as
  the `chord` field — twist is measured relative to it. Defaults to auto-derived
  from the closest VSM panel during SystemStructure construction.
- `y_airf=nothing`: Spanwise reference (body frame). Auto-derived when omitted.
- `area=NaN`: Surface area [m²] for flat-plate (`AeroPlate`) sections.
- `twist=0.0`: Initial twist angle [rad] (prescribed input for `STATIC` sections).

# Returns
- `TwistSurface`: A new `TwistSurface` object. The `idx` and `point_idxs` are resolved by SystemStructure.
  When `x_airf`/`y_airf` are omitted the geometry fields (le_pos, chord, y_airf) are
  computed during SystemStructure construction from the closest VSM panel.
"""
function TwistSurface(name, points, type, moment_frac;
                      damping=50.0, x_airf=nothing, y_airf=nothing,
                      area=NaN, twist=0.0)
    point_refs = Vector{NameRef}([p isa Integer ? Int(p) : Symbol(p) for p in points])
    chord_vec = isnothing(x_airf) ? zeros(KVec3) : KVec3(x_airf)
    y_vec = isnothing(y_airf) ? zeros(KVec3) : KVec3(y_airf)
    TwistSurface(0, name, Int64[], point_refs,
          zeros(KVec3), chord_vec, y_vec,
          type, moment_frac, damping,
          SimFloat(twist), 0.0, 0.0, 0.0, 0.0,
          Int64[], SimFloat(area))
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
    "Index in the segments vector (assigned by SystemStructure)."
    idx::Int64
    "Name used for lookup by other components' `_ref` fields."
    const name::Union{Int, Symbol, Nothing}
    "Resolved endpoint indices (filled by SystemStructure)."
    point_idxs::Tuple{Int64, Int64}
    "Raw endpoint references (names or indices)."
    const point_refs::Tuple{NameRef, NameRef}
    "Stiffness per unit length [N]. Effective k = unit_stiffness/length [N/m]."
    unit_stiffness::SimFloat
    "Damping per unit length [N·s]. Effective c = unit_damping/length [N·s/m]."
    unit_damping::SimFloat
    "Rest (unstretched) length [m]."
    l0::SimFloat
    "Compressive/tensile stiffness ratio (0-1). 0 = no compression stiffness."
    compression_frac::SimFloat
    "Segment diameter [m]."
    diameter::SimFloat
    "Material density [kg/m³]."
    density::SimFloat
    "Current length [m] (updated during simulation)."
    len::SimFloat
    "Current force [N] (updated during simulation)."
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

# Keyword Arguments
- `density::SimFloat=NaN`: Material density [kg/m³] used for mass.
"""
function Segment(name, point_i, point_j, unit_stiffness, unit_damping, diameter;
    l0=zero(SimFloat), compression_frac=0.1, density=NaN
)
    p1 = point_i isa Integer ? Int(point_i) : Symbol(point_i)
    p2 = point_j isa Integer ? Int(point_j) : Symbol(point_j)
    Segment(0, name, (0, 0), (p1, p2), unit_stiffness, unit_damping, l0, compression_frac,
        diameter, SimFloat(density), zero(SimFloat), zero(SimFloat))
end

"""
    Segment(name, set, point_i, point_j; l0, compression_frac,
            diameter_mm, unit_stiffness, unit_damping)

Constructs a `Segment` using settings for material properties.

# Arguments
- `name::Union{Int, Symbol}`: Name/identifier for the segment.
- `set::Settings`: The settings object containing material properties.
- `point_i`, `point_j`: References to the two endpoint points
  (names or indices).

# Keyword Arguments
- `l0::SimFloat=zero(SimFloat)`: Unstretched length [m].
  Calculated from point positions if zero.
- `compression_frac::SimFloat=0.0`: Compressive/tensile stiffness
  ratio (0-1). 0 = no compression stiffness.
- `diameter_mm::Float64=NaN`: Tether diameter [mm]. If `NaN`,
  uses `set.d_tether`.
- `unit_stiffness::Float64=NaN`: Stiffness per unit length [N].
  Effective k = unit_stiffness/length.
- `unit_damping::Float64=NaN`: Damping per unit length [N·s].
  Effective c = unit_damping/length.
- `density::Float64=NaN`: Material density [kg/m³]. If `NaN`,
  uses `set.rho_tether`.
"""
function Segment(name, set, point_i, point_j;
    l0=zero(SimFloat), compression_frac=0.0,
    diameter_mm=NaN, unit_stiffness=NaN,
    unit_damping=NaN, density=NaN
)
    p1 = point_i isa Integer ? Int(point_i) : Symbol(point_i)
    p2 = point_j isa Integer ? Int(point_j) : Symbol(point_j)

    # Set default diameter from settings if not specified
    if isnan(diameter_mm)
        diameter_mm = set.d_tether
    end
    # Convert diameter from mm to m
    diameter_m = 0.001 * diameter_mm

    # Compute unit_stiffness if not provided
    if isnan(unit_stiffness)
        unit_stiffness = set.e_tether * (diameter_m/2)^2 * π
    end

    # Compute unit_damping if not provided
    if isnan(unit_damping)
        if hasproperty(set, :rel_damping) &&
                set.rel_damping != 0.0
            unit_damping = set.rel_damping * unit_stiffness
        elseif hasproperty(set, :unit_damping) &&
                hasproperty(set, :unit_stiffness) &&
                set.unit_damping != 0.0
            unit_damping = (set.unit_damping /
                set.unit_stiffness) * unit_stiffness
        else
            @warn "Segment $(name): unit_damping is zero " *
                "(no rel_damping or unit_damping in settings)."
            unit_damping = 0.0
        end
    end

    if isnan(density)
        density = set.rho_tether
    end

    Segment(0, name, (0, 0), (p1, p2),
        unit_stiffness, unit_damping, l0,
        compression_frac, diameter_m, SimFloat(density),
        zero(SimFloat), zero(SimFloat))
end

"""Deprecated: SegmentType parameter removed."""
function Segment(_name, _set, _p_i, _p_j, _type::SegmentType;
             _kw...)
    error("Segment `type` (SegmentType) parameter removed. " *
          "Use `unit_stiffness`, `unit_damping`, " *
          "`diameter_mm` kwargs, or a YAML material.")
end

# ==================== PULLEY ==================== #

"""
    mutable struct Pulley

A pulley described by two segments with the common point of the segments being the pulley.

$(TYPEDFIELDS)
"""
mutable struct Pulley
    "Index in the pulleys vector (assigned by SystemStructure)."
    idx::Int64
    "Name used for lookup by other components' `_ref` fields."
    const name::Union{Int, Symbol, Nothing}
    "Resolved segment indices (filled by SystemStructure)."
    segment_idxs::Tuple{Int64, Int64}
    "Raw segment references (names or indices)."
    const segment_refs::Tuple{NameRef, NameRef}
    "Dynamics type (DYNAMIC)."
    const type::DynamicsType
    "Sum of connected segment lengths [m]."
    sum_len::SimFloat
    "Current pulley length [m] (updated during simulation)."
    len::SimFloat
    "Current pulley velocity [m/s] (updated during simulation)."
    vel::SimFloat
end

"""
    Pulley(name, segment_i, segment_j, type)

Constructs a `Pulley` object that enforces length redistribution between two segments.

# Arguments
- `name::Union{Int, Symbol}`: Name/identifier for the pulley.
- `segment_i`, `segment_j`: References to the two segments (names or indices).
- `type::DynamicsType`: Dynamics type (`DYNAMIC`).
"""
function Pulley(name, segment_i, segment_j, type)
    s1 = segment_i isa Integer ? Int(segment_i) : Symbol(segment_i)
    s2 = segment_j isa Integer ? Int(segment_j) : Symbol(segment_j)
    return Pulley(0, name, (0, 0), (s1, s2), type, 0.0, 0.0, 0.0)
end

# ==================== TETHER ==================== #

"""
    mutable struct Tether

A collection of segments forming a flexible line.

Can be constructed two ways:
- **Route 1** (explicit segments): Provide segment references directly.
- **Route 2** (auto-generation): Provide start/end points and `n_segments`;
  intermediate points and segments are created by `expand_auto_tethers!`.

$(TYPEDFIELDS)
"""
mutable struct Tether
    "Index in the tethers vector (assigned by SystemStructure)."
    idx::Int64
    "Name used for lookup by other components' `_ref` fields."
    const name::Union{Int, Symbol, Nothing}
    "Resolved segment indices (filled by SystemStructure)."
    segment_idxs::Vector{Int64}
    "Raw segment references (names or indices)."
    const segment_refs::Vector{NameRef}
    "Resolved start point index (filled by SystemStructure)."
    start_point_idx::Int64
    "Raw start point reference. Nothing for Route 1."
    const start_point_ref::Union{NameRef, Nothing}
    "Resolved end point index (filled by SystemStructure)."
    end_point_idx::Int64
    "Raw end point reference. Nothing for Route 1."
    const end_point_ref::Union{NameRef, Nothing}
    "Number of segments (Route 2 only)."
    const n_segments::Int64
    "Stiffness per unit length [N]. NaN = derive from Settings."
    const unit_stiffness::SimFloat
    "Damping per unit length [N·s]. NaN = derive from Settings."
    const unit_damping::SimFloat
    "Tether diameter [m]. NaN = derive from Settings."
    const diameter::SimFloat
    "Material density [kg/m³]. NaN = derive from Settings."
    const density::SimFloat
    "Current stretched length [m] (updated during simulation)."
    stretched_len::SimFloat
    """Unstretched tether length [m] (sum of segment l0).
    ODE state variable. Segment l0 = len / n_segments."""
    len::SimFloat
    """Initial stretched standoff [m] — the placed point
    geometry (Σ segment norms). Drives placement of root
    tethers. `nothing` = use the geometric (CAD) length,
    i.e. no scaling."""
    init_stretched_len::Union{SimFloat, Nothing}
    """Target initial spring force [N], default 0. `reinit!`
    solves the unstretched `len` from the placed stretched
    length: `len = stretched · (1 − force/unit_stiffness)`.
    Mutually exclusive with `init_stretch_frac`."""
    init_tether_force::Union{SimFloat, Nothing}
    """Initial unstretched/stretched length fraction. `reinit!`
    sets `len = init_stretch_frac · stretched`; 0.9 gives 10%
    pre-stretch, 1.0 no tension, >1.0 slack. Must be positive.
    Mutually exclusive with `init_tether_force`."""
    init_stretch_frac::Union{SimFloat, Nothing}
end

function Base.setproperty!(t::Tether, name::Symbol, x)
    if name === :init_stretch_frac
        isnothing(x) || setfield!(t, :init_tether_force, nothing)
        setfield!(t, :init_stretch_frac, x)
    elseif name === :init_tether_force
        isnothing(x) || setfield!(t, :init_stretch_frac, nothing)
        setfield!(t, :init_tether_force, x)
    else
        setfield!(t, name, x)
    end
end

"""
    Tether(name, segments, stretched_length=nothing;
           start_point=nothing, end_point=nothing,
           tether_force=nothing, stretch_frac=nothing)

Route 1: Construct a `Tether` from explicit segment references.

# Arguments
- `name::Union{Int, Symbol}`: Name/identifier for the tether.
- `segments::Vector`: References to segments (names or indices).
- `stretched_length=nothing`: Stretched standoff [m] (placed point
  geometry). Drives placement of root tethers. `nothing` = use the
  geometric length.

# Keyword Arguments
- `start_point=nothing`: Optional start point ref.
- `end_point=nothing`: Optional end point ref.
- `tether_force=nothing`: Target initial spring force [N], default 0.
- `stretch_frac=nothing`: Initial `len/stretched` fraction. Mutually
  exclusive with `tether_force`.
"""
function Tether(name, segments::AbstractVector, stretched_length=nothing;
                start_point=nothing, end_point=nothing,
                winch_point=nothing, tether_force=nothing,
                stretch_frac=nothing)
    if !isnothing(winch_point)
        error("`winch_point` moved from Tether to " *
              "Winch. Use Tether(name, segments, " *
              "len) and pass winch_point to the " *
              "Winch constructor.")
    end
    init_force, init_frac =
        resolve_tether_init(name, tether_force, stretch_frac)
    segment_refs = Vector{NameRef}(name_ref.(segments))
    init_stretched = opt_simfloat(stretched_length)
    return Tether(0, name, Int64[], segment_refs,
                  0, name_ref(start_point), 0, name_ref(end_point),
                  length(segments),
                  NaN, NaN, NaN, NaN, 0.0,
                  0.0, init_stretched, init_force, init_frac)
end

name_ref(::Nothing) = nothing
name_ref(x::Integer) = Int(x)
name_ref(x) = Symbol(x)

opt_simfloat(::Nothing) = nothing
opt_simfloat(x) = SimFloat(x)

function resolve_tether_init(name, tether_force, stretch_frac)
    if !isnothing(tether_force) && !isnothing(stretch_frac)
        error("Tether $name: set only one of `tether_force` and " *
              "`stretch_frac`.")
    end
    !isnothing(stretch_frac) && return nothing, SimFloat(stretch_frac)
    !isnothing(tether_force) && return SimFloat(tether_force), nothing
    return SimFloat(0.0), nothing
end

"""
    Tether(name, stretched_length=nothing;
           start_point, end_point, n_segments,
           unit_stiffness=NaN, unit_damping=NaN,
           diameter=NaN, tether_force=nothing, stretch_frac=nothing)

Route 2: Construct a `Tether` for auto-generation of intermediate
points and segments by `expand_auto_tethers!`.

# Arguments
- `name::Union{Int, Symbol}`: Name/identifier for the tether.
- `stretched_length=nothing`: Stretched standoff [m] (placed point
  geometry). Drives placement of root tethers. `nothing` = use the
  geometric length.

# Keyword Arguments
- `start_point`: Reference to the start point (required).
- `end_point`: Reference to the end point (required).
- `n_segments::Int`: Number of segments to generate (required).
- `unit_stiffness::Float64=NaN`: Per-unit-length stiffness [N].
  NaN = derive from Settings during auto-expansion.
- `unit_damping::Float64=NaN`: Per-unit-length damping [N·s].
  NaN = derive from Settings during auto-expansion.
- `diameter::Float64=NaN`: Tether diameter [m].
  NaN = derive from Settings during auto-expansion.
- `density::Float64=NaN`: Material density [kg/m³].
  NaN = derive from Settings during auto-expansion.
- `tether_force=nothing`: Target initial spring force [N], default 0.
- `stretch_frac=nothing`: Initial `len/stretched` fraction. Mutually
  exclusive with `tether_force`.
"""
function Tether(name, stretched_length=nothing;
                start_point, end_point, n_segments,
                unit_stiffness=NaN, unit_damping=NaN,
                diameter=NaN, density=NaN, tether_force=nothing,
                stretch_frac=nothing)
    init_force, init_frac =
        resolve_tether_init(name, tether_force, stretch_frac)
    seg_refs = Vector{NameRef}(
        [Symbol("$(name)_seg_$i") for i in 1:n_segments])
    init_stretched = opt_simfloat(stretched_length)
    return Tether(0, name, Int64[], seg_refs,
                  0, name_ref(start_point), 0, name_ref(end_point),
                  Int64(n_segments),
                  Float64(unit_stiffness),
                  Float64(unit_damping),
                  Float64(diameter), Float64(density), 0.0,
                  0.0, init_stretched, init_force, init_frac)
end

# ==================== WINCH ==================== #

"""
    abstract type AbstractWinchModel

Selects the winch motor dynamics. Each concrete model carries its own
parameter fields and adds a `winch_component(model, sys_struct, idx; name,
params)` method building the MTK subsystem (mirrors [`AbstractAeroModel`](@ref)).
The model lives in [`Winch`](@ref)`.model`; common drum parameters
(`gear_ratio`, `drum_radius`, `f_coulomb`, `c_vf`, `inertia_total`) stay on
the `Winch` struct. New models = new struct + a few methods; see
[`TorqueWinch`](@ref).
"""
abstract type AbstractWinchModel end

"""
    mutable struct Winch

A set of tethers (or a single tether) connected to a winch mechanism.

$(TYPEDFIELDS)
"""
mutable struct Winch
    "Index in the winches vector (assigned by SystemStructure)."
    idx::Int64
    "Name used for lookup by other components' `_ref` fields."
    const name::Union{Int, Symbol, Nothing}
    "Resolved tether indices (filled by SystemStructure)."
    tether_idxs::Vector{Int64}
    "Raw tether references (names or indices)."
    const tether_refs::Vector{NameRef}
    "Resolved winch point index (filled by SystemStructure)."
    winch_point_idx::Int64
    "Raw winch point reference (name or index)."
    const winch_point_ref::NameRef
    "Initial reel-out velocity [m/s]. Applied on reinit!."
    init_vel::SimFloat
    "Current reel-out velocity [m/s]. ODE state variable."
    vel::SimFloat
    "Current winch acceleration [m/s²] from motor dynamics."
    acc::SimFloat
    """Abstract setpoint passed to the winch component as the
    `set_value` connector. Interpretation is the component's
    choice (e.g. motor torque, current, set velocity, set length).
    The default component treats it as motor torque [N·m]."""
    set_value::SimFloat
    """Brake input in [0, 1]. The outer integrator freezes
    `winch_vel` and `tether_len` when `> 0.5`; custom components
    may interpret intermediate values as a continuous brake."""
    brake::SimFloat
    """If true, reel-out velocity is prescribed externally rather
    than integrated from motor dynamics: winch acceleration is forced
    to 0 (ignoring `model`). Set the velocity via `winch.vel`."""
    speed_controlled::Bool
    "Force vector at winch point [N]."
    const force::KVec3
    "Gear ratio [-]."
    gear_ratio::SimFloat
    "Drum radius [m]."
    drum_radius::SimFloat
    "Coulomb friction force [N]."
    f_coulomb::SimFloat
    "Viscous friction coefficient [N·s/m]."
    c_vf::SimFloat
    "Total rotational inertia [kg·m²]."
    inertia_total::SimFloat
    "Current friction force [N] (updated during simulation)."
    friction::SimFloat
    """Winch motor dynamics model carrying its own parameters.
    Defaults to [`TorqueWinch`](@ref). See [`AbstractWinchModel`](@ref)."""
    model::AbstractWinchModel
end

"""
    Winch(name, set, tethers; winch_point, ...)

Constructs a `Winch` object that controls tether length through
torque or speed regulation.

# Arguments
- `name::Union{Int, Symbol}`: Name/identifier for the winch.
- `set::Settings`: Settings object for winch parameters.
- `tethers::Vector`: References to tethers connected to this
  winch (names or indices).

# Keyword Arguments
- `winch_point`: Reference to the ground attachment point
  (name or index). Required.
- `init_vel::SimFloat=0.0`: Initial reel-out rate [m/s].
- `brake=0.0`: Brake input in [0, 1]. `> 0.5` engages a hard
  freeze on `winch_vel` and `tether_len` at the outer integrator.
- `speed_controlled::Bool=false`: If true, prescribe reel-out
  velocity via `winch.vel` instead of integrating motor dynamics;
  winch acceleration is forced to 0, ignoring `model`.
- `friction_epsilon::SimFloat=6.0`: Coulomb-friction smoothing width;
  forwarded into the default [`TorqueWinch`](@ref) (ignored when an
  explicit non-default `model` is passed).
- `model::AbstractWinchModel=TorqueWinch()`: Winch motor dynamics
  model. See [`AbstractWinchModel`](@ref) for plugging in your own.
"""
function Winch(name, set::Settings, tethers;
               winch_point,
               init_vel=0.0, brake=0.0, speed_controlled=false,
               friction_epsilon=6.0,
               model::AbstractWinchModel=TorqueWinch(; friction_epsilon))
    tether_refs = Vector{NameRef}(
        [t isa Integer ? Int(t) : Symbol(t) for t in tethers])
    wp = winch_point isa Integer ? Int(winch_point) :
         Symbol(winch_point)
    return Winch(0, name, Int64[], tether_refs, 0, wp,
                 init_vel, 0.0,
                 0.0, 0.0,
                 SimFloat(brake), speed_controlled, zeros(KVec3),
                 set.gear_ratio, set.drum_radius,
                 set.f_coulomb, set.c_vf,
                 set.inertia_total, zero(SimFloat),
                 model)
end

"""
    Winch(name, tethers, gear_ratio, drum_radius, f_coulomb,
          c_vf, inertia_total; winch_point, ...)

Constructs a `Winch` by directly providing physical parameters.

# Arguments
- `name::Union{Int, Symbol}`: Name/identifier for the winch.
- `tethers::Vector`: References to tethers (names or indices).
- `gear_ratio`, `drum_radius`, `f_coulomb`, `c_vf`,
  `inertia_total`: Physical parameters.

# Keyword Arguments
- `winch_point`: Reference to ground attachment point. Required.
- `init_vel::SimFloat=0.0`: Initial reel-out rate [m/s].
"""
function Winch(name, tethers, gear_ratio, drum_radius,
               f_coulomb, c_vf, inertia_total;
               winch_point,
               init_vel=0.0, brake=0.0, speed_controlled=false,
               friction_epsilon=6.0,
               model::AbstractWinchModel=TorqueWinch(; friction_epsilon))
    tether_refs = Vector{NameRef}(
        [t isa Integer ? Int(t) : Symbol(t) for t in tethers])
    wp = winch_point isa Integer ? Int(winch_point) :
         Symbol(winch_point)
    return Winch(0, name, Int64[], tether_refs, 0, wp,
                 init_vel, 0.0,
                 0.0, 0.0,
                 SimFloat(brake), speed_controlled, zeros(KVec3),
                 gear_ratio, drum_radius, f_coulomb,
                 c_vf, inertia_total, zero(SimFloat),
                 model)
end

# ==================== TRANSFORM ==================== #

"""
    mutable struct Transform

Describes the spatial transformation (position and orientation) of system components
relative to a base reference point.

$(TYPEDFIELDS)
"""
mutable struct Transform
    "Index in the transforms vector (assigned by SystemStructure)."
    idx::Int64
    "Name used for lookup by other components' `_ref` fields."
    const name::Union{Int, Symbol, Nothing}
    "Resolved wing index (filled by SystemStructure)."
    wing_idx::Union{Int64, Nothing}
    "Raw wing reference (name or index). Nothing = uses rot_point."
    const wing_ref::Union{NameRef, Nothing}
    "Resolved rotation point index (filled by SystemStructure)."
    rot_point_idx::Union{Int64, Nothing}
    "Raw rotation point reference. Nothing = uses wing."
    const rot_point_ref::Union{NameRef, Nothing}
    "Resolved base point index (filled by SystemStructure)."
    base_point_idx::Union{Int64, Nothing}
    "Raw base point reference."
    const base_point_ref::Union{NameRef, Nothing}
    "Resolved base transform index (filled by SystemStructure)."
    base_transform_idx::Union{Int64, Nothing}
    "Raw base transform reference. Nothing = uses base_pos."
    const base_transform_ref::Union{NameRef, Nothing}
    "Elevation angle [rad]."
    elevation::SimFloat
    "Azimuth angle [rad]."
    azimuth::SimFloat
    "Heading angle [rad]."
    heading::SimFloat
    "Angular velocity in elevation direction [rad/s]."
    elevation_vel::SimFloat
    "Angular velocity in azimuth direction [rad/s]."
    azimuth_vel::SimFloat
    "Angular velocity around radial axis [rad/s]."
    turn_rate::SimFloat
    "Base position [m]. Nothing = derived from base_transform."
    base_pos::Union{KVec3, Nothing}
end

# Helper to convert ref to NameRef or nothing
to_ref(::Nothing) = nothing
to_ref(x::Integer) = Int(x)
to_ref(x) = Symbol(x)

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

    wing_ref = to_ref(wing)
    rot_point_ref = to_ref(rot_point)
    base_point_ref = to_ref(base_point)
    base_transform_ref = to_ref(base_transform)

    Transform(0, name, nothing, wing_ref, nothing, rot_point_ref,
              nothing, base_point_ref, nothing, base_transform_ref,
              elevation, azimuth, heading, elevation_vel, azimuth_vel, turn_rate,
              isnothing(base_pos) ? nothing : KVec3(base_pos...))
end

"""
    get_rot_pos(transform::Transform, bodies, points)

Get the world position of the rotating object (body or point).
"""
function get_rot_pos(transform::Transform, bodies, points)
    wing_idx = transform.wing_idx
    if !isnothing(wing_idx)
        return bodies[something(wing_idx)].pos_w
    end
    rot_point_idx = transform.rot_point_idx
    if !isnothing(rot_point_idx)
        return points[something(rot_point_idx)].pos_w
    end
    error("Transform #$(transform.idx): " *
        "neither wing_idx nor rot_point_idx is set")
end

"""
    get_rot_pos_cad(transform::Transform, bodies, points)

Get the CAD-frame position of the rotating object (body or point).
Used by `get_base_pos` to compute the translation offset for
chained transforms.
"""
function get_rot_pos_cad(transform::Transform, bodies, points)
    wing_idx = transform.wing_idx
    if !isnothing(wing_idx)
        return bodies[something(wing_idx)].pos_cad
    end
    rot_point_idx = transform.rot_point_idx
    if !isnothing(rot_point_idx)
        return points[something(rot_point_idx)].pos_cad
    end
    error("Transform #$(transform.idx): " *
        "neither wing_idx nor rot_point_idx is set")
end

"""
    get_base_pos(transform, transforms, bodies, points)

Get `(base_pos, curr_base_pos)` for a transform.

For chained transforms (`base_transform`): returns the parent's
current world position and CAD position, so
`T = base_pos - curr_base_pos` shifts child points by the same
displacement the parent transform applied.

For direct transforms (`base_pos` + `base_point`): returns the
user-specified position and the base point's current position.
"""
function get_base_pos(transform::Transform,
        transforms, bodies, points)
    base_transform_idx = transform.base_transform_idx
    if !isnothing(base_transform_idx)
        base_tf = transforms[something(
            base_transform_idx)]
        rot_pos = get_rot_pos(base_tf, bodies, points)
        rot_pos_cad = get_rot_pos_cad(
            base_tf, bodies, points)
        return rot_pos, rot_pos_cad
    end
    curr_base_pos = points[something(
        transform.base_point_idx)].pos_w
    base_pos = transform.base_pos
    if !isnothing(base_pos)
        return something(base_pos), curr_base_pos
    end
    error("Transform #$(transform.idx): neither " *
        "base_pos nor base_transform_idx is set")
end
