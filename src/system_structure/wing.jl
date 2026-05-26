# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

"""
Wing types and VSM-related wing code.

This file contains:
- AbstractWing, BaseWing, VSMWing structs
- Wing constructors and helper functions
- VSM panel adjustment functions
"""

# ==================== ABSTRACT WING ==================== #

"""
    abstract type AbstractWing

Abstract base type for all wing implementations.

Concrete subtypes must implement rigid body dynamics and provide a reference frame
for attached points and groups.
"""
abstract type AbstractWing end

# ==================== BASE WING ==================== #

"""
    mutable struct BaseWing <: AbstractWing

A rigid wing body that can have multiple groups of points attached to it.

The wing provides a rigid body reference frame for attached points and groups.
Points with `type == WING` move rigidly with the wing body according to the
wing's orientation matrix `R_b_to_w` and position `pos_w`.

# Special Properties
The wing's orientation can be accessed as a rotation matrix or a quaternion:
```julia
R_matrix = wing.R_b_to_w
wing.R_b_to_w = R_matrix

quat = wing.Q_b_to_w
wing.Q_b_to_w = quat
```

$(TYPEDFIELDS)
"""
mutable struct BaseWing <: AbstractWing
    idx::Int64  # Assigned by SystemStructure based on vector position
    const name::Union{Int, Symbol, Nothing}  # Name/identifier (Int for backwards compat)

    # Structural information - resolved indices
    group_idxs::Vector{Int64}  # Resolved by SystemStructure from group_refs
    transform_idx::Int64       # Resolved by SystemStructure from transform_ref

    # Structural information - raw references
    const group_refs::Vector{NameRef}   # Raw references to groups (names or indices)
    const transform_ref::NameRef        # Raw reference to transform (name or idx)

    # Geometry
    const R_b_to_c::Matrix{SimFloat}       # Body frame → CAD (from ref points)
    const R_p_to_c::Matrix{SimFloat}       # Principal frame → CAD (from inertia diag)
    const R_b_to_p::Matrix{SimFloat}       # Body → principal (= R_p_to_c' * R_b_to_c, constant)
    const pos_cad::KVec3                # Body origin in CAD frame
    const com_offset_b::KVec3           # COM offset from body origin in body frame
    const inertia_principal::KVec3
    const dynamics_type::WingType
    aero_mode::AeroMode

    # Principal frame ODE state (RIGID_DYNAMICS dynamics)
    const com_w::KVec3                  # COM position in world frame
    const com_vel::KVec3                # COM velocity in world frame
    const Q_p_to_w::Vector{SimFloat}       # Principal frame quaternion (length 4)
    const ω_p::KVec3                    # Angular velocity in principal frame

    # Body frame output (algebraic for RIGID_DYNAMICS, from ref points)
    const Q_b_to_w::Vector{SimFloat}
    const ω_b::KVec3
    const pos_w::KVec3                  # Body origin world position
    const vel_w::KVec3                  # Body origin velocity
    const acc_w::KVec3                  # Body origin acceleration

    # Derived variables and parameters, updated during simulation
    const wind_disturb::KVec3
    drag_frac::SimFloat
    const va_b::KVec3 # apparent wind in body frame
    const v_wind::KVec3 # wind velocity in world frame at the wing
    const aero_force_b::KVec3 # aerodynamic force in body frame
    const aero_moment_b::KVec3 # aerodynamic moment in body frame
    const tether_moment::KVec3 # tether moment in world frame
    const tether_force::KVec3 # tether force in world frame
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
    y_damping::SimFloat
    angular_damping::SimFloat
    z_disturb::SimFloat
    mass::SimFloat  # Total mass of wing (sum of WING point masses if set.mass is zero)
end

function Base.getproperty(wing::BaseWing, sym::Symbol)
    if sym == :R_b_to_w
        return quaternion_to_rotation_matrix(getfield(wing, :Q_b_to_w))
    elseif sym == :R_p_to_w
        return quaternion_to_rotation_matrix(getfield(wing, :Q_p_to_w))
    else
        return getfield(wing, sym)
    end
end

function Base.setproperty!(wing::BaseWing, sym::Symbol, value)
    if sym == :R_b_to_w
        if value isa AbstractMatrix
            getfield(wing, :Q_b_to_w) .= rotation_matrix_to_quaternion(value)
        else
            error("Cannot set R_b_to_w with non-matrix value of type $(typeof(value))")
        end
    elseif hasfield(BaseWing, sym)
        setfield!(wing, sym, value)
    else
        error("BaseWing has no field `$(sym)`")
    end
end

# ==================== WEIGHTED REF POINTS ==================== #

"""
    WeightedRefPoints

Weighted combination of reference points for wing frame
definition. Supports single points, equal-weight averaging,
and arbitrary weight combinations.

# Fields
- `refs`: Unresolved names/indices (filled at construction)
- `ids`: Resolved point indices (filled by `resolve!`)
- `weights`: Normalized weights (sum to 1.0)
"""
mutable struct WeightedRefPoints
    const refs::Vector{NameRef}
    ids::Vector{Int64}
    const weights::Vector{Float64}
end

"""Single point from a Symbol ref."""
WeightedRefPoints(ref::Symbol) =
    WeightedRefPoints(NameRef[ref], Int64[], [1.0])

"""
    WeightedRefPoints(id::Integer)

Single point from a resolved index. Stores in `ids`
(not `refs`) so `resolve!` is a no-op.
"""
WeightedRefPoints(id::Integer) =
    WeightedRefPoints(NameRef[], Int64[Int64(id)], [1.0])

"""
    WeightedRefPoints(refs::AbstractVector)

Equal-weight average of multiple ref points, or weighted
if elements are `(name, weight)` tuples.

Supports:
- `[:le, :te]` → equal-weight average
- `[(:le, 0.7), (:te, 0.3)]` → weighted combination
"""
function WeightedRefPoints(refs::AbstractVector)
    isempty(refs) && error(
        "WeightedRefPoints requires at least one " *
        "reference point, got empty vector")
    if refs[1] isa Tuple
        names = NameRef[_to_name_ref(t[1]) for t in refs]
        weights = Float64[Float64(t[2]) for t in refs]
        _validate_weights!(weights)
        return WeightedRefPoints(names, Int64[], weights)
    end
    names = NameRef[_to_name_ref(v) for v in refs]
    n = length(names)
    WeightedRefPoints(names, Int64[], fill(1.0 / n, n))
end

# ==================== VSM WING ==================== #

"""
    mutable struct VSMWing <: AbstractWing

A wing that uses the Vortex Step Method (VSM) for aerodynamic computations.

This struct extends the base wing functionality with VSM-specific aerodynamic
modeling capabilities, including vortex wake computations and aerodynamic loads.

$(TYPEDFIELDS)
"""
mutable struct VSMWing{BA<:VortexStepMethod.BodyAerodynamics,
                       W<:VortexStepMethod.AbstractWing,
                       SL<:VortexStepMethod.Solver} <: AbstractWing
    # Base wing functionality
    base::BaseWing

    # VSM aerodynamics
    vsm_aero::BA
    vsm_wing::W
    vsm_solver::SL

    # Aerodynamic linearization state (RIGID_DYNAMICS)
    # aero_y: operating point inputs
    #   [alpha, beta, ω1, ω2, ω3, twist...]
    # aero_x: baseline wind-axis coefficients
    #   [CL, CD, CS, CM1, CM2, CM3, cm_1..cm_n]
    # aero_jac: dense Jacobian d(aero_x)/d(aero_y)
    aero_y::Vector{SimFloat}
    aero_x::Vector{SimFloat}
    aero_jac::Matrix{SimFloat}

    # PARTICLE_DYNAMICS-specific fields (Nothing for RIGID_DYNAMICS wings)
    point_to_vsm_point::Union{Nothing, Dict{Int64, Tuple{Int64, Symbol}}}
    wing_segments::Union{Nothing, Vector{Tuple{Int64, Int64}}}

    # Orientation reference points (WeightedRefPoints carry
    # both refs and resolved ids, with weights).
    # Z-axis: Normal to wing plane
    # Y-axis: Spanwise, X = Y × Z (chord)
    z_ref_points::Union{Nothing,
        Tuple{WeightedRefPoints, WeightedRefPoints}}
    y_ref_points::Union{Nothing,
        Tuple{WeightedRefPoints, WeightedRefPoints}}

    # Origin point - RESOLVED index
    # Defines wing.pos_w = pos[:, origin_idx]
    origin_idx::Union{Nothing, Int64}

    # KCU origin point - RAW reference (name or idx)
    const origin_ref::Union{Nothing, NameRef}

    # Additional aerodynamic force scale to compensate chord length errors (PARTICLE_DYNAMICS)
    aero_scale_chord::SimFloat

    # Body frame z-axis offset for VSM aerodynamics (RIGID_DYNAMICS only)
    # Shifts VSM panel positions in positive z direction (body frame)
    # to adjust moment arm for improved stability
    aero_z_offset::SimFloat

    function VSMWing(base::BaseWing, vsm_aero,
                     vsm_wing, vsm_solver,
                     aero_y, aero_x, aero_jac,
                     point_to_vsm_point, wing_segments,
                     z_ref_points, y_ref_points,
                     origin_idx, origin_ref,
                     aero_scale_chord, aero_z_offset)
        new{typeof(vsm_aero), typeof(vsm_wing),
            typeof(vsm_solver)}(
            base, vsm_aero, vsm_wing, vsm_solver,
            aero_y, aero_x, aero_jac,
            point_to_vsm_point, wing_segments,
            z_ref_points, y_ref_points,
            origin_idx, origin_ref,
            aero_scale_chord, aero_z_offset)
    end
end

# Delegate property access to base wing for VSMWing
const VSM_WING_OWN_FIELDS = (
    :base, :vsm_aero, :vsm_wing, :vsm_solver,
    :aero_y, :aero_x, :aero_jac,
    :point_to_vsm_point, :wing_segments,
    :z_ref_points, :y_ref_points,
    :origin_idx, :origin_ref,
    :aero_scale_chord, :aero_z_offset)

function Base.getproperty(wing::VSMWing, sym::Symbol)
    if sym in VSM_WING_OWN_FIELDS
        return getfield(wing, sym)
    elseif sym == :vsm_aoa
        # Compute mean angle of attack from VSM solver solution
        solver = getfield(wing, :vsm_solver)
        return mean(solver.sol.alpha_dist)
    else
        return getproperty(getfield(wing, :base), sym)
    end
end

function Base.setproperty!(wing::VSMWing, sym::Symbol, value)
    if sym in VSM_WING_OWN_FIELDS
        setfield!(wing, sym, value)
    else
        setproperty!(getfield(wing, :base), sym, value)
    end
end

# ==================== CONSTRUCTORS ==================== #

# Helper to convert to NameRef
_to_name_ref(x::Integer) = Int(x)
_to_name_ref(x) = Symbol(x)

"""
    BaseWing(name, groups, R_b_to_c, pos_cad, inertia_principal; transform=nothing, y_damping=150.0, dynamics_type=RIGID_DYNAMICS)

Constructs a `BaseWing` object representing a rigid body reference frame.

# Arguments
- `name::Union{Int, Symbol}`: Name/identifier for the wing (e.g., `:main_wing` or `1` for legacy).
- `groups::Vector`: References to groups attached to this wing (names or indices).
- `R_b_to_c::Matrix{SimFloat}`: Rotation matrix from body frame to CAD frame.
- `pos_cad::KVec3`: Position of wing body origin in CAD frame.
- `inertia_principal::KVec3`: Principal moments of inertia [Ixx, Iyy, Izz] in principal frame.

# Keyword Arguments
- `transform=nothing`: Reference to the transform (name or index). Defaults to 1 if not specified.
- `y_damping::SimFloat=150.0`: Damping coefficient for y-axis (pitch) rotation.
- `angular_damping::SimFloat=0.0`: Isotropic angular damping on all 3 rotation axes.
- `dynamics_type::WingType=RIGID_DYNAMICS`: Wing aerodynamic model type.

# Returns
- `BaseWing`: A new base wing object. The `idx`, `group_idxs`, and `transform_idx` are resolved by SystemStructure.
"""
function BaseWing(name, groups::AbstractVector, R_b_to_c::AbstractMatrix,
                  pos_cad, inertia_principal;
                  transform=nothing, y_damping=150.0,
                  angular_damping=0.0,
                  dynamics_type::Union{Nothing,WingType}=nothing,
                  aero_mode::Union{Nothing,AeroMode}=nothing,
                  wing_type::Union{Nothing,WingType}=nothing)
    # Handle deprecated wing_type keyword
    if !isnothing(wing_type)
        if !isnothing(dynamics_type)
            error("Cannot specify both `wing_type` and `dynamics_type`; `wing_type` is deprecated, use `dynamics_type`.")
        end
        @warn "Keyword argument `wing_type` is deprecated; use `dynamics_type` instead."
        dynamics_type = wing_type
    end
    # Apply defaults now that dynamics_type is resolved
    isnothing(dynamics_type) && (dynamics_type = RIGID_DYNAMICS)
    isnothing(aero_mode) && (aero_mode = dynamics_type == RIGID_DYNAMICS ?
        AERO_LINEARIZED : AERO_DIRECT)
    # Convert groups to NameRef vector
    group_refs = Vector{NameRef}([_to_name_ref(g) for g in groups])
    # Handle nothing - default to transform 1
    tf = isnothing(transform) ? 1 : transform
    transform_ref = _to_name_ref(tf)

    # idx, group_idxs, transform_idx are placeholders - resolved by SystemStructure
    return BaseWing(0, name,
        # Structural information - resolved (placeholders)
        Int64[], 0,
        # Structural information - raw references
        group_refs, transform_ref,
        # Geometry
        R_b_to_c, Matrix{SimFloat}(I, 3, 3),  # R_p_to_c placeholder
        Matrix{SimFloat}(I, 3, 3),         # R_b_to_p placeholder
        pos_cad, zeros(KVec3),  # com_offset_b placeholder
        inertia_principal, dynamics_type,
        aero_mode,
        # Principal frame ODE state
        zeros(KVec3), zeros(KVec3),  # com_w, com_vel
        zeros(4), zeros(KVec3),      # Q_p_to_w, ω_p
        # Body frame output
        zeros(4), zeros(KVec3),      # Q_b_to_w, ω_b
        zeros(KVec3), zeros(KVec3), zeros(KVec3),  # pos_w, vel_w, acc_w
        # Derived variables and parameters, updated during simulation
        zeros(KVec3), one(SimFloat),
        zeros(KVec3), zeros(KVec3), zeros(KVec3), zeros(KVec3),
        zeros(KVec3), zeros(KVec3),
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        zeros(KVec3), zeros(KVec3), 0.0, 0.0, false,
        y_damping, angular_damping, 0.0,
        0.0)  # mass initialized to 0, set by SystemStructure
end

"""Warn and normalize if weights don't sum to 1."""
function _validate_weights!(weights::Vector{Float64})
    s = sum(weights)
    s > 0 || error(
        "Ref point weights sum to $s; " *
        "all weights must be positive")
    if !isapprox(s, 1.0; atol=1e-6)
        @warn "Ref point weights sum to $s, " *
            "normalizing to 1.0"
        weights ./= s
    end
end

"""
    create_vsm_wing(set::Settings, vsm_set::VortexStepMethod.VSMSettings; prn=true, kwargs...)

Create a `Wing` geometry object from the settings provided.

This function checks for .obj and .dat files in the model directory.
If found, it uses `VortexStepMethod.ObjWing(obj_path, dat_path)` to load the wing.
Otherwise, it falls back to loading from `aero_geometry.yaml`.

Reads geometry from the `Settings` object and initializes the `Wing` object
from `VortexStepMethod.jl`.
"""
function create_vsm_wing(set::Settings, vsm_set::VortexStepMethod.VSMSettings; prn=true, kwargs...)
    # Check for .obj and .dat files in the model directory
    model_dir = get_data_path()
    obj_path = joinpath(model_dir, set.model)
    dat_path = joinpath(model_dir, set.foil_file)

    if isfile(obj_path) && isfile(dat_path)
        # Use ObjWing constructor (default path)
        prn && @info "Loading wing from .obj/.dat files"

        if set.physical_model == "simple_ram"
            n_unrefined_sections = 2
        else
            n_unrefined_sections = 4
        end

        return VortexStepMethod.ObjWing(obj_path, dat_path;
            mass=set.mass, crease_frac=set.crease_frac, n_unrefined_sections,
            align_to_principal=true, prn, kwargs...
        )
    end

    # Fallback: load from aero_geometry.yaml using provided vsm_set
    prn && @info "Using provided VSMSettings for wing creation"
    # Resolve relative geometry_file paths against data dir
    for ws in vsm_set.wings
        gf = ws.geometry_file
        if !isempty(gf) && !isabspath(gf)
            ws.geometry_file = joinpath(model_dir, basename(gf))
        end
    end
    return VortexStepMethod.Wing(vsm_set; kwargs...)
end

"""
    VSMWing(name, set, groups, vsm_set; transform=nothing, y_damping=150.0, ...)

Constructs a `VSMWing` object with Vortex Step Method aerodynamics.
Creates vsm_wing, vsm_aero, and vsm_solver internally.

# Arguments
- `name::Union{Int, Symbol}`: Name/identifier for the wing.
- `set::Settings`: Settings object for VSM configuration.
- `groups::Vector`: References to groups (names or indices). Used by both
  `RIGID_DYNAMICS` and `PARTICLE_DYNAMICS` wings during
  `match_aero_sections_to_structure!` for LE/TE panel identification.
  For `PARTICLE_DYNAMICS` wings, `group_idxs` is cleared from the wing
  after section matching (groups remain in the `SystemStructure`).
- `vsm_set::VortexStepMethod.VSMSettings`: VSM settings for wing creation.

# Keyword Arguments
- `transform=nothing`: Reference to the transform (name or index). Defaults to 1 if not specified.
- `R_b_to_c::Matrix{SimFloat}`: Rotation matrix body→CAD.
- `pos_cad::KVec3`: Position of wing COM in CAD frame.
- `y_damping::SimFloat=150.0`: Lateral damping coefficient.
- `dynamics_type::WingType=RIGID_DYNAMICS`: Aerodynamic model type.
- `point_to_vsm_point`: 1:1 structural point to VSM point mapping (PARTICLE_DYNAMICS only).
- `wing_segments`: LE/TE pairs (populated for all VSM wing types by
  `match_aero_sections_to_structure!`).
- `z_ref_points`: Chord direction reference points (PARTICLE_DYNAMICS only, names or indices).
- `y_ref_points`: Span direction reference points (PARTICLE_DYNAMICS only, names or indices).
- `origin`: Reference to origin point (PARTICLE_DYNAMICS only, name or index).
- `aero_z_offset::SimFloat=0.0`: Body frame z-offset for VSM panels (RIGID_DYNAMICS only).

# Returns
- `VSMWing`: A new VSM wing object. References are resolved by SystemStructure.
"""
function VSMWing(name, set::Settings,
                 groups::AbstractVector,
                 vsm_set::VortexStepMethod.VSMSettings;
                 R_b_to_c::Union{Nothing,AbstractMatrix}=nothing,
                 pos_cad::Union{Nothing,AbstractVector}=nothing,
                 transform=nothing, y_damping=150.0,
                 angular_damping=0.0,
                 inertia_diag=nothing,
                 dynamics_type::Union{Nothing,WingType}=nothing,
                 aero_mode::Union{Nothing,AeroMode}=nothing,
                 wing_type::Union{Nothing,WingType}=nothing,
                 point_to_vsm_point::Union{Nothing, Dict{Int64, Tuple{Int64, Symbol}}}=nothing,
                 wing_segments::Union{Nothing, Vector{Tuple{Int64, Int64}}}=nothing,
                 z_ref_points=nothing,
                 y_ref_points=nothing,
                 origin=nothing,
                 aero_scale_chord::SimFloat=0.0,
                 aero_z_offset::SimFloat=0.0)

    # Handle deprecated wing_type keyword
    if !isnothing(wing_type)
        if !isnothing(dynamics_type)
            error("Cannot specify both `wing_type` and `dynamics_type`; `wing_type` is deprecated, use `dynamics_type`.")
        end
        @warn "Keyword argument `wing_type` is deprecated; use `dynamics_type` instead."
        dynamics_type = wing_type
    end
    # Apply defaults now that dynamics_type is resolved
    isnothing(dynamics_type) && (dynamics_type = RIGID_DYNAMICS)
    isnothing(aero_mode) && (aero_mode = dynamics_type == RIGID_DYNAMICS ?
        AERO_LINEARIZED : AERO_DIRECT)

    # Validation
    if dynamics_type == PARTICLE_DYNAMICS
        @assert !isnothing(origin)
            "PARTICLE_DYNAMICS wings require origin to define KCU position"
        if !isnothing(pos_cad)
            @warn "Wing '$name': pos_cad is unused for " *
                "PARTICLE_DYNAMICS wings (position comes from " *
                "origin point)"
            pos_cad = nothing
        end
    else
        @assert isnothing(point_to_vsm_point)
            "RIGID_DYNAMICS wings: no point_to_vsm_point"
        # origin, z_ref_points, y_ref_points are now
        # accepted for RIGID_DYNAMICS wings (body frame
        # defined by structural ref points)
    end

    # Convert ref points to WeightedRefPoints
    z_ref = isnothing(z_ref_points) ? nothing :
        (WeightedRefPoints(z_ref_points[1]),
         WeightedRefPoints(z_ref_points[2]))
    y_ref = isnothing(y_ref_points) ? nothing :
        (WeightedRefPoints(y_ref_points[1]),
         WeightedRefPoints(y_ref_points[2]))
    origin_ref = isnothing(origin) ? nothing :
        _to_name_ref(origin)

    # Create VSM wing, aero, and solver
    vsm_wing = create_vsm_wing(set, vsm_set; prn=false,
        sort_sections=false)
    vsm_aero = VortexStepMethod.BodyAerodynamics([vsm_wing])
    vsm_solver = VortexStepMethod.Solver(vsm_aero, vsm_set)

    # Placeholders — overwritten by SystemStructure
    # from point masses (RIGID_DYNAMICS) or ref points (PARTICLE_DYNAMICS)
    isnothing(R_b_to_c) && (R_b_to_c = Matrix{SimFloat}(I, 3, 3))
    isnothing(pos_cad) && (pos_cad = zeros(KVec3))
    inertia_vec = isnothing(inertia_diag) ?
        ones(MVector{3, SimFloat}) : inertia_diag

    base = BaseWing(name, groups, R_b_to_c, pos_cad,
                    inertia_vec; transform, y_damping,
                    angular_damping, dynamics_type, aero_mode)

    # Size aero state vectors based on wing type
    # For RIGID_DYNAMICS: placeholder sizes using n_unrefined
    # as group count proxy; resized in SystemStructure
    # after groups are resolved.
    if dynamics_type == PARTICLE_DYNAMICS
        nx = 0
        ny = 0
    else
        n_groups_est = vsm_wing.n_unrefined_sections
        nx = 6 + n_groups_est
        ny = 5 + n_groups_est
    end

    return VSMWing(base, vsm_aero, vsm_wing, vsm_solver,
                   zeros(SimFloat, ny), zeros(SimFloat, nx),
                   zeros(SimFloat, nx, ny),
                   point_to_vsm_point, wing_segments,
                   z_ref, y_ref,
                   nothing, origin_ref,
                   aero_scale_chord, aero_z_offset)
end

"""
    VSMWing(name, vsm_aero, vsm_wing, vsm_solver, groups, R_b_to_c, pos_cad; transform=nothing)

Legacy constructor accepting pre-created VSM objects directly.
Kept for backward compatibility with predefined structures.

# Arguments
- `name::Union{Int, Symbol}`: Wing name/identifier
- `vsm_aero`: Pre-created BodyAerodynamics
- `vsm_wing`: Pre-created VortexStepMethod.Wing
- `vsm_solver`: Pre-created Solver
- `groups`: References to groups (names or indices)
- `R_b_to_c`: Rotation matrix body→CAD
- `pos_cad`: Position in CAD frame

# Keyword Arguments
- `transform=nothing`: Reference to the transform. Defaults to 1 if not specified.

# Returns
- `VSMWing`: Wing with RIGID_DYNAMICS type
"""
function VSMWing(name, vsm_aero, vsm_wing, vsm_solver,
                 groups::AbstractVector,
                 R_b_to_c::AbstractMatrix,
                 pos_cad::AbstractVector;
                 transform=nothing)
    # Placeholder inertia — overwritten by SystemStructure
    inertia_vec = ones(MVector{3, SimFloat})
    base = BaseWing(name, groups, R_b_to_c, pos_cad,
                    inertia_vec; transform)
    # Placeholder aero arrays — resized by SystemStructure
    n_groups_est = vsm_wing.n_unrefined_sections
    nx = 6 + n_groups_est
    ny = 5 + n_groups_est
    return VSMWing(base, vsm_aero, vsm_wing, vsm_solver,
        zeros(SimFloat, ny), zeros(SimFloat, nx),
        zeros(SimFloat, nx, ny),
        nothing, nothing,
        nothing, nothing,  # z/y_ref_points
        nothing, nothing,  # origin_idx and origin_ref
        0.0, 0.0)
end

"""
    Wing(name, vsm_aero, vsm_wing, vsm_solver, groups, R_b_to_c, pos_cad; transform=1)

Constructs a `VSMWing` object (backward compatibility constructor).

This is a convenience constructor that creates a VSMWing for backward compatibility
with existing code. New code should use `VSMWing(...)` directly.

# Arguments
- `name::Union{Int, Symbol}`: Name/identifier for the wing.
- `vsm_aero`, `vsm_wing`, `vsm_solver`: Vortex Step Method components.
- `groups::Vector`: References to groups attached to this wing (names or indices).
- `R_b_to_c::Matrix{SimFloat}`: Rotation matrix from body frame to CAD frame.
- `pos_cad::KVec3`: Position of wing center of mass in CAD frame.

# Keyword Arguments
- `transform::Union{Int, Symbol}=1`: Reference to the transform.
- `y_damping::SimFloat=150.0`: Damping coefficient for lateral motion.

# Returns
- `VSMWing`: A new VSM wing object.
"""
function SymbolicAWEModels.Wing(name, vsm_aero, vsm_wing, vsm_solver, groups, R_b_to_c,
                                pos_cad; kwargs...)
    return VSMWing(name, vsm_aero, vsm_wing, vsm_solver, groups, R_b_to_c, pos_cad; kwargs...)
end

# ==================== PLATE SURFACE ==================== #

"""
    struct PlateSurface

A flat aerodynamic plate defined by orientation vectors, area,
and a center-of-pressure WING point. Internal to PlateWing.

$(TYPEDFIELDS)
"""
mutable struct PlateSurface
    "Name identifier for this surface."
    name::Union{Symbol, Nothing}
    "Chord direction in body frame (unit vector)."
    x_airf::KVec3
    "Span direction in body frame (unit vector)."
    y_airf::KVec3
    "Plate area [m²]."
    area::SimFloat
    "Raw reference to center-of-pressure WING point."
    point_ref::NameRef
    "Resolved point index (filled by SystemStructure)."
    point_idx::Int64
    "Twist angle [rad] (mutable control input)."
    twist::SimFloat
    "Current AoA [deg] (updated by update_sys_struct!)."
    aoa::SimFloat
end

"""
    PlateSurface(name, x_airf, y_airf, area, point;
                 twist=0.0)

Construct a PlateSurface with the given geometry. The
`point_idx` is resolved later by SystemStructure.
"""
function PlateSurface(name, x_airf, y_airf, area, point;
                      twist=0.0)
    ref = point isa Integer ? Int(point) : Symbol(point)
    PlateSurface(
        isnothing(name) ? nothing : Symbol(name),
        KVec3(x_airf), KVec3(y_airf), area,
        ref, 0,
        twist, 0.0)
end

# ==================== PLATE WING ==================== #

"""
    mutable struct PlateWing <: AbstractWing

A wing with flat-plate CL/CD aerodynamics. Each PlateSurface
computes lift and drag from angle-of-attack lookup tables.
Twist is set directly on each PlateSurface.

Supports both RIGID_DYNAMICS (rigid body) and PARTICLE_DYNAMICS (point mass)
wing dynamics via BaseWing.dynamics_type.

$(TYPEDFIELDS)
"""
mutable struct PlateWing <: AbstractWing
    "Base wing functionality."
    base::BaseWing
    "Plate surfaces (one per aerodynamic plate)."
    surfaces::Vector{PlateSurface}
    "Z-axis reference points for body frame."
    z_ref_points::Union{Nothing,
        Tuple{WeightedRefPoints, WeightedRefPoints}}
    "Y-axis reference points for body frame."
    y_ref_points::Union{Nothing,
        Tuple{WeightedRefPoints, WeightedRefPoints}}
    "Resolved origin point index."
    origin_idx::Union{Nothing, Int64}
    "Raw origin point reference."
    const origin_ref::Union{Nothing, NameRef}
    "CL lookup: callable(alpha_deg) → CL."
    calc_cl::Any
    "CD lookup: callable(alpha_deg) → CD."
    calc_cd::Any
    "Drag correction factor (0.93 for KPS4)."
    drag_corr::SimFloat
    "Pitch moment coefficient."
    cmq::SimFloat
    "Steering moment coefficient."
    smc::SimFloat
    "Mean aerodynamic chord [m]."
    cord_length::SimFloat
end

# Delegate property access to base wing for PlateWing
const PLATE_WING_OWN_FIELDS = (
    :base, :surfaces,
    :z_ref_points, :y_ref_points,
    :origin_idx, :origin_ref,
    :calc_cl, :calc_cd,
    :drag_corr, :cmq, :smc, :cord_length)

function Base.getproperty(wing::PlateWing, sym::Symbol)
    if sym in PLATE_WING_OWN_FIELDS
        return getfield(wing, sym)
    else
        return getproperty(getfield(wing, :base), sym)
    end
end

function Base.setproperty!(wing::PlateWing, sym::Symbol, value)
    if sym in PLATE_WING_OWN_FIELDS
        setfield!(wing, sym, value)
    else
        setproperty!(getfield(wing, :base), sym, value)
    end
end

"""
    PlateWing(name, surfaces, calc_cl, calc_cd;
              dynamics_type=PARTICLE_DYNAMICS, transform=nothing,
              y_damping=150.0, angular_damping=0.0, drag_corr=0.93,
              cmq=1.0, smc=1.0, cord_length=1.0,
              z_ref_points=nothing, y_ref_points=nothing,
              origin=nothing)

Construct a PlateWing with flat-plate aerodynamics.

# Arguments
- `name`: Wing name/identifier.
- `surfaces`: Vector of PlateSurface definitions.
- `calc_cl`: CL lookup callable(alpha_deg) → CL.
- `calc_cd`: CD lookup callable(alpha_deg) → CD.

# Keyword Arguments
- `dynamics_type`: `RIGID_DYNAMICS` or `PARTICLE_DYNAMICS` (default).
- `transform`: Reference to transform (name or index).
- `y_damping`: Damping coefficient for y-axis (pitch) rotation.
- `angular_damping`: Angular damping coefficient.
- `drag_corr`: Drag correction factor.
- `cmq`: Pitch moment coefficient.
- `smc`: Steering moment coefficient.
- `cord_length`: Mean aerodynamic chord [m].
- `z_ref_points`, `y_ref_points`: Body frame references.
- `origin`: Origin point reference.
"""
function PlateWing(name, surfaces::Vector{PlateSurface},
                   calc_cl, calc_cd;
                   dynamics_type::WingType=PARTICLE_DYNAMICS,
                   transform=nothing,
                   y_damping=150.0,
                   angular_damping=0.0,
                   drag_corr=0.93,
                   cmq=1.0, smc=1.0, cord_length=1.0,
                   z_ref_points=nothing,
                   y_ref_points=nothing,
                   origin=nothing)
    # PlateWing has no groups
    base = BaseWing(name, NameRef[], Matrix{SimFloat}(I, 3, 3),
                    zeros(KVec3), ones(MVector{3, SimFloat});
                    transform, y_damping, angular_damping,
                    dynamics_type, aero_mode=AERO_PLATE)

    z_ref = isnothing(z_ref_points) ? nothing :
        (WeightedRefPoints(z_ref_points[1]),
         WeightedRefPoints(z_ref_points[2]))
    y_ref = isnothing(y_ref_points) ? nothing :
        (WeightedRefPoints(y_ref_points[1]),
         WeightedRefPoints(y_ref_points[2]))
    origin_ref = isnothing(origin) ? nothing :
        _to_name_ref(origin)

    PlateWing(base, surfaces,
              z_ref, y_ref,
              nothing, origin_ref,
              calc_cl, calc_cd,
              drag_corr, cmq, smc, cord_length)
end

"""
    plate_alpha(wing::PlateWing, surf::PlateSurface)

Compute current AoA [deg] from body-frame apparent wind and
twist. Requires `va_b` to be up to date.
"""
function plate_alpha(wing::PlateWing, surf::PlateSurface)
    tw = surf.twist
    ct, st = cos(tw), sin(tw)
    x_tw = ct * surf.x_airf + st * (surf.y_airf × surf.x_airf)
    z_tw = x_tw × surf.y_airf
    v_tan = wing.va_b ⋅ x_tw
    v_norm = wing.va_b ⋅ z_tw
    rad2deg(atan(v_norm, v_tan))
end

# ==================== HELPER FUNCTIONS ==================== #

"""
    adjust_vsm_panels_to_origin!(vsm_wing, origin_offset)

Adjust VSM panel positions when body frame origin changes.

When RIGID_DYNAMICS wings are loaded from YAML, the panel positions in aero_geometry.yaml
are specified in an absolute body frame. However, the body frame origin is adjusted
to the mean position of all WING points. This function updates all panel positions
to be relative to the new origin by subtracting the offset.

# Arguments
- `vsm_wing`: VortexStepMethod.Wing with sections to adjust
- `origin_offset`: Vector [x, y, z] to subtract from panel positions
"""
function adjust_vsm_panels_to_origin!(vsm_wing, origin_offset)
    for section in vsm_wing.refined_sections
        section.LE_point .-= origin_offset
        section.TE_point .-= origin_offset
    end
    for section in vsm_wing.non_deformed_sections
        section.LE_point .-= origin_offset
        section.TE_point .-= origin_offset
    end
    for section in vsm_wing.unrefined_sections
        section.LE_point .-= origin_offset
        section.TE_point .-= origin_offset
    end
end

"""
    rotate_vsm_sections!(vsm_wing, R)

Rotate all VSM section LE/TE points by rotation matrix `R`.

Used during initialization to transform sections from CAD
frame to body frame. After the first step, `update_vsm!()`
updates positions from `pos_b` (already in body frame).
"""
function rotate_vsm_sections!(vsm_wing, R)
    for section in vsm_wing.refined_sections
        section.LE_point .= R * section.LE_point
        section.TE_point .= R * section.TE_point
    end
    for section in vsm_wing.non_deformed_sections
        section.LE_point .= R * section.LE_point
        section.TE_point .= R * section.TE_point
    end
    for section in vsm_wing.unrefined_sections
        section.LE_point .= R * section.LE_point
        section.TE_point .= R * section.TE_point
    end
end

"""
    apply_aero_z_offset!(vsm_wing, aero_z_offset)

Apply z-axis offset to VSM panel positions in body frame.

For RIGID_DYNAMICS wings, this shifts the aerodynamic center of pressure
in the positive z-direction (body frame) to adjust the moment arm.
This is applied AFTER the COM adjustment.

# Arguments
- `vsm_wing`: VortexStepMethod.Wing with sections to adjust
- `aero_z_offset`: Distance to shift panels in +z direction [m]
"""
function apply_aero_z_offset!(vsm_wing, aero_z_offset)
    if abs(aero_z_offset) > 1e-10
        offset_vec = [0.0, 0.0, aero_z_offset]
        for section in vsm_wing.refined_sections
            section.LE_point .+= offset_vec
            section.TE_point .+= offset_vec
        end
        for section in vsm_wing.non_deformed_sections
            section.LE_point .+= offset_vec
            section.TE_point .+= offset_vec
        end
        for section in vsm_wing.unrefined_sections
            section.LE_point .+= offset_vec
            section.TE_point .+= offset_vec
        end
    end
end

"""
    calc_pos(wing::Wing, gamma, frac)

Calculate a position on the kite based on spanwise (`gamma`) and chordwise (`frac`) parameters.
"""
function calc_pos(wing::VortexStepMethod.Wing, gamma, frac)
    le_pos = [wing.le_interp[i](gamma) for i in 1:3]
    chord = [wing.te_interp[i](gamma) for i in 1:3] .- le_pos
    pos = le_pos .+ chord .* frac
    return pos
end

"""
    cad_to_body_frame(wing::Wing, pos)

Transform a position from the CAD frame to the wing's body frame.
"""
function cad_to_body_frame(wing::VortexStepMethod.Wing, pos)
    return wing.R_cad_body * (pos + wing.T_cad_body)
end
