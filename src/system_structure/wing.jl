# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

"""
Wing type and VSM-related wing code.

This file contains:
- AbstractWing, Wing structs (abstract aero types + VSMEngine live in types.jl;
  concrete aero modes live in src/aero_modes/)
- Wing constructors and helper functions
- VSM panel adjustment functions
"""

# ==================== ABSTRACT WING ==================== #

"""
    abstract type AbstractWing

Abstract base type for all wing implementations.

Concrete subtypes must implement rigid body dynamics and provide a reference frame
for attached points and twist_surfaces.
"""
abstract type AbstractWing end

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

"""Single point from a String ref (converted to Symbol)."""
WeightedRefPoints(ref::AbstractString) =
    WeightedRefPoints(Symbol(ref))

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
        names = NameRef[to_name_ref(entry[1]) for entry in refs]
        weights = Float64[Float64(entry[2]) for entry in refs]
        validate_weights!(weights)
        return WeightedRefPoints(names, Int64[], weights)
    end
    names = NameRef[to_name_ref(v) for v in refs]
    n = length(names)
    WeightedRefPoints(names, Int64[], fill(1.0 / n, n))
end

# ==================== VSM ENGINE ==================== #

# ==================== WING ==================== #

"""
    mutable struct Wing <: AbstractWing

A wing body that can have multiple twist_surfaces of points attached to it.

The wing provides a body reference frame for attached points and twist_surfaces.
Points with `type == WING` move with the wing body according to the wing's
orientation matrix `R_b_to_w` and position `pos_w`. Its `dynamics_type` selects
rigid-body (`RIGID_DYNAMICS`) or per-particle (`PARTICLE_DYNAMICS`) behaviour, and
its [`aero`](@ref AbstractAeroModel) field selects the aerodynamic model. When the
mode is a VSM mode ([`AbstractVSMAero`](@ref)) its [`VSMEngine`](@ref) fields
(`vsm_wing`, `aero_x`, …) are forwarded through the wing.

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
mutable struct Wing <: AbstractWing
    idx::Int64  # Assigned by SystemStructure based on vector position
    const name::Union{Int, Symbol, Nothing}  # Name/identifier (Int for backwards compat)

    # Structural information - resolved indices
    twist_surface_idxs::Vector{Int64}  # Resolved by SystemStructure from twist_surface_refs
    transform_idx::Int64       # Resolved by SystemStructure from transform_ref

    # Structural information - raw references
    const twist_surface_refs::Vector{NameRef}   # Raw references to twist_surfaces (names or indices)
    const transform_ref::NameRef        # Raw reference to transform (name or idx)

    # Geometry
    const R_b_to_c::Matrix{SimFloat}       # Body frame → CAD (from ref points)
    const R_p_to_c::Matrix{SimFloat}       # Principal frame → CAD (from inertia diag)
    const R_b_to_p::Matrix{SimFloat}       # Body → principal (= R_p_to_c' * R_b_to_c, constant)
    const pos_cad::KVec3                # Body origin in CAD frame
    const com_offset_b::KVec3           # COM offset from body origin in body frame
    const inertia_principal::KVec3
    const dynamics_type::WingType
    aero::AbstractAeroModel

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
    "Whether in-group (twist_surface) points contribute their moment to the wing body."
    group_points_moment::Bool
    y_damping::SimFloat
    angular_damping::SimFloat
    z_disturb::SimFloat
    mass::SimFloat  # Total mass of wing (sum of WING point masses if set.mass is zero)

    # Body-frame reference points (define R_b_to_w / pos_w from structural points)
    z_ref_points::Union{Nothing, Tuple{WeightedRefPoints, WeightedRefPoints}}
    y_ref_points::Union{Nothing, Tuple{WeightedRefPoints, WeightedRefPoints}}
    origin::Union{Nothing, WeightedRefPoints}
end

function Base.getproperty(wing::Wing, sym::Symbol)
    if sym === :R_b_to_w
        return quaternion_to_rotation_matrix(getfield(wing, :Q_b_to_w))
    elseif sym === :R_p_to_w
        return quaternion_to_rotation_matrix(getfield(wing, :Q_p_to_w))
    elseif sym === :vsm_aoa
        return mean(wing_vsm_engine(wing, sym).vsm_solver.sol.alpha_dist)
    elseif sym in VSM_ENGINE_FIELDS
        return getproperty(wing_vsm_engine(wing, sym), sym)
    else
        return getfield(wing, sym)
    end
end

function Base.setproperty!(wing::Wing, sym::Symbol, value)
    if sym === :R_b_to_w
        if value isa AbstractMatrix
            getfield(wing, :Q_b_to_w) .= rotation_matrix_to_quaternion(value)
        else
            error("Cannot set R_b_to_w with non-matrix value of type $(typeof(value))")
        end
    elseif sym in VSM_ENGINE_FIELDS
        setproperty!(wing_vsm_engine(wing, sym), sym, value)
    elseif hasfield(Wing, sym)
        setfield!(wing, sym, value)
    else
        error("Wing has no field `$(sym)`")
    end
end

# Return the wing's VSM engine (so VSM-field access gives a clear error on
# wings whose aero mode carries no engine, e.g. AeroNone/AeroPlate).
function wing_vsm_engine(wing::Wing, sym::Symbol)
    engine = vsm_engine(getfield(wing, :aero))
    engine === nothing && error(
        "Wing $(getfield(wing, :name)) aero mode " *
        "$(typeof(getfield(wing, :aero))) has no VSM engine; " *
        "field `$sym` unavailable.")
    return engine
end

"""
    count_aero_log_points(wings) -> Int

Total extra `SysState` log slots contributed by the wings' aero modes
([`n_aero_log_points`](@ref)): 4 per VSM panel for VSM modes, 4 per section
quad for flat-plate wings, 0 by default.
"""
function count_aero_log_points(wings)
    total = 0
    for wing in wings
        total += n_aero_log_points(getfield(wing, :aero), wing)
    end
    return total
end

# ==================== CONSTRUCTORS ==================== #

# Helper to convert to NameRef
to_name_ref(x::Integer) = Int(x)
to_name_ref(x) = Symbol(x)

"""Warn and normalize if weights don't sum to 1."""
function validate_weights!(weights::Vector{Float64})
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
    Wing(name, twist_surfaces, R_b_to_c, pos_cad, inertia_principal;
         transform=nothing, y_damping=150.0, angular_damping=0.0,
         dynamics_type=RIGID_DYNAMICS, aero=nothing,
         z_ref_points=nothing, y_ref_points=nothing, origin=nothing, vsm=nothing)

Low-level `Wing` constructor. The `idx`, `twist_surface_idxs`, and `transform_idx`
are resolved by `SystemStructure`.

# Arguments
- `name::Union{Int, Symbol}`: Name/identifier (e.g. `:main_wing`).
- `twist_surfaces::Vector`: References to attached twist_surfaces (names or indices).
- `R_b_to_c::Matrix{SimFloat}`: Rotation matrix from body frame to CAD frame.
- `pos_cad::KVec3`: Position of wing body origin in CAD frame.
- `inertia_principal::KVec3`: Principal moments of inertia `[Ixx, Iyy, Izz]`.

# Keyword Arguments
- `transform`: Reference to the transform (name or index). Defaults to 1.
- `y_damping`, `angular_damping`: Damping coefficients.
- `dynamics_type::WingType`: `RIGID_DYNAMICS` (default) or `PARTICLE_DYNAMICS`.
- `aero::AbstractAeroModel`: Aerodynamic model (defaults by `dynamics_type`).
- `group_points_moment::Bool=true`: When `false`, in-group (twist_surface) points
  add no moment to the wing body; their force still contributes. Runtime-switchable.
- `z_ref_points`, `y_ref_points`, `origin`: Body-frame reference points (raw refs).
  A VSM engine, when needed, lives in the `aero` mode (built by the VSM constructors).
"""
function Wing(name, twist_surfaces::AbstractVector, R_b_to_c::AbstractMatrix,
              pos_cad, inertia_principal;
              transform=nothing, y_damping=150.0,
              angular_damping=0.0,
              dynamics_type::Union{Nothing,WingType}=nothing,
              aero::Union{Nothing,AbstractAeroModel}=nothing,
              wing_type::Union{Nothing,WingType}=nothing,
              group_points_moment::Bool=true,
              z_ref_points=nothing, y_ref_points=nothing, origin=nothing)
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
    isnothing(aero) && (aero = dynamics_type == RIGID_DYNAMICS ?
        AeroLinearized() : AeroDirect())
    twist_surface_refs = Vector{NameRef}([to_name_ref(twist_surface) for twist_surface in twist_surfaces])
    transform_value = isnothing(transform) ? 1 : transform
    transform_ref = to_name_ref(transform_value)

    z_ref = isnothing(z_ref_points) ? nothing :
        (WeightedRefPoints(z_ref_points[1]), WeightedRefPoints(z_ref_points[2]))
    y_ref = isnothing(y_ref_points) ? nothing :
        (WeightedRefPoints(y_ref_points[1]), WeightedRefPoints(y_ref_points[2]))
    origin_rp = isnothing(origin) ? nothing : WeightedRefPoints(origin)

    # idx, twist_surface_idxs, transform_idx are placeholders - resolved by SystemStructure
    return Wing(0, name,
        # Structural information - resolved (placeholders)
        Int64[], 0,
        # Structural information - raw references
        twist_surface_refs, transform_ref,
        # Geometry
        R_b_to_c, Matrix{SimFloat}(I, 3, 3),  # R_p_to_c placeholder
        Matrix{SimFloat}(I, 3, 3),         # R_b_to_p placeholder
        pos_cad, zeros(KVec3),  # com_offset_b placeholder
        inertia_principal, dynamics_type,
        aero,
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
        group_points_moment,
        y_damping, angular_damping, 0.0,
        0.0,  # mass initialized to 0, set by SystemStructure
        z_ref, y_ref, origin_rp)
end

"""
    create_vsm_wing(set::Settings, vsm_set::VortexStepMethod.VSMSettings; prn=true, sort_sections=true)

Create a `VortexStepMethod.Wing` geometry object from the settings provided.

This function checks for .obj and .dat files in the model directory.
If found, it uses `VortexStepMethod.ObjWing(obj_path, dat_path)` to load the wing.
Otherwise, it falls back to loading from `aero_geometry.yaml`.
"""
function create_vsm_wing(set::Settings, vsm_set::VortexStepMethod.VSMSettings; prn=true, sort_sections=true)
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
            mass=1.0, crease_frac=set.crease_frac, n_unrefined_sections,
            align_to_principal=false, prn
        )
    end

    # Fallback: load from aero_geometry.yaml using provided vsm_set
    prn && @info "Using provided VSMSettings for wing creation"
    # Resolve relative geometry_file paths against data dir
    for wing_settings in vsm_set.wings
        geometry_file = wing_settings.geometry_file
        if !isempty(geometry_file) && !isabspath(geometry_file)
            wing_settings.geometry_file = joinpath(model_dir, basename(geometry_file))
        end
    end
    return VortexStepMethod.Wing(vsm_set; sort_sections)
end

"""
    build_vsm_engine(set, vsm_set, dynamics_type; point_to_vsm_point=nothing,
                     wing_segments=nothing, aero_scale_chord=0.0, aero_z_offset=0.0)

Build a [`VSMEngine`](@ref): create the VortexStepMethod `vsm_wing`/`vsm_aero`/
`vsm_solver` and size the linearization state vectors. Aero-state sizes are
placeholders for `RIGID_DYNAMICS` (using `n_unrefined_sections` as the
twist_surface-count proxy) and resized by `SystemStructure` once twist_surfaces
are resolved.
"""
function build_vsm_engine(set::Settings, vsm_set::VortexStepMethod.VSMSettings,
                          dynamics_type::WingType;
                          point_to_vsm_point=nothing, wing_segments=nothing,
                          aero_scale_chord=0.0, aero_z_offset=0.0)
    vsm_wing = create_vsm_wing(set, vsm_set; prn=false, sort_sections=false)
    vsm_aero = VortexStepMethod.BodyAerodynamics([vsm_wing])
    vsm_solver = VortexStepMethod.Solver(vsm_aero, vsm_set)

    if dynamics_type == PARTICLE_DYNAMICS
        num_aero_outputs = 0
        num_aero_inputs = 0
    else
        n_twist_surfaces_est = vsm_wing.n_unrefined_sections
        num_aero_outputs = 6 + n_twist_surfaces_est
        num_aero_inputs = 5 + n_twist_surfaces_est
    end

    return VSMEngine(vsm_aero, vsm_wing, vsm_solver,
        zeros(SimFloat, num_aero_inputs),
        zeros(SimFloat, num_aero_outputs),
        zeros(SimFloat, num_aero_outputs, num_aero_inputs),
        point_to_vsm_point, wing_segments,
        SimFloat(aero_scale_chord), SimFloat(aero_z_offset))
end

"""
    VSMWing(name, set, twist_surfaces, vsm_set; transform=nothing, y_damping=150.0, ...)

Construct a [`Wing`](@ref) with Vortex Step Method aerodynamics. Builds the
[`VSMEngine`](@ref) (`vsm_wing`/`vsm_aero`/`vsm_solver`) internally and attaches
it to the wing.

# Arguments
- `name::Union{Int, Symbol}`: Name/identifier for the wing.
- `set::Settings`: Settings object for VSM configuration.
- `twist_surfaces::Vector`: References to twist_surfaces (names or indices).
- `vsm_set`: VSM settings for engine creation. Required for VSM-backed aero
  modes ([`AbstractVSMAero`](@ref)); may be `nothing` for engine-less modes
  like [`AeroNone`](@ref).

# Keyword Arguments
- `transform=nothing`: Reference to the transform. Defaults to 1.
- `R_b_to_c`, `pos_cad`, `inertia_diag`: Geometry placeholders (resolved later).
- `y_damping`, `angular_damping`: Damping coefficients.
- `dynamics_type::WingType=RIGID_DYNAMICS`: Aerodynamic model type.
- `aero::AbstractAeroModel`: Aerodynamic model (defaults by `dynamics_type`).
- `group_points_moment::Bool=true`: When `false`, in-group (twist_surface) points
  add no moment to the wing body; their force still contributes. Runtime-switchable.
- `point_to_vsm_point`, `wing_segments`: VSM structural↔panel maps.
- `z_ref_points`, `y_ref_points`, `origin`: Body-frame references.
- `aero_scale_chord`, `aero_z_offset`: VSM force/panel adjustments.
"""
function VSMWing(name, set::Settings,
                 twist_surfaces::AbstractVector,
                 vsm_set::Union{Nothing, VortexStepMethod.VSMSettings};
                 R_b_to_c::Union{Nothing,AbstractMatrix}=nothing,
                 pos_cad::Union{Nothing,AbstractVector}=nothing,
                 transform=nothing, y_damping=150.0,
                 angular_damping=0.0,
                 inertia_diag=nothing,
                 dynamics_type::Union{Nothing,WingType}=nothing,
                 aero::Union{Nothing,AbstractAeroModel}=nothing,
                 wing_type::Union{Nothing,WingType}=nothing,
                 group_points_moment::Bool=true,
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
    isnothing(dynamics_type) && (dynamics_type = RIGID_DYNAMICS)

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
    end

    # Resolve the aero mode and, when it is VSM-backed, build and attach the
    # engine. Engine-less modes (AeroNone/AeroPlate) need no vsm_set.
    isnothing(aero) && (aero = dynamics_type == RIGID_DYNAMICS ?
        AeroLinearized() : AeroDirect())
    if aero isa AbstractVSMAero
        isnothing(vsm_set) && error(
            "Wing '$name': aero mode $(typeof(aero)) needs VSM geometry " *
            "but no vsm_set was provided.")
        aero.engine = build_vsm_engine(set, vsm_set, dynamics_type;
            point_to_vsm_point, wing_segments, aero_scale_chord, aero_z_offset)
    end

    # Placeholders — overwritten by SystemStructure
    isnothing(R_b_to_c) && (R_b_to_c = Matrix{SimFloat}(I, 3, 3))
    isnothing(pos_cad) && (pos_cad = zeros(KVec3))
    inertia_vec = isnothing(inertia_diag) ?
        ones(MVector{3, SimFloat}) : inertia_diag

    return Wing(name, twist_surfaces, R_b_to_c, pos_cad, inertia_vec;
        transform, y_damping, angular_damping, dynamics_type, aero,
        group_points_moment, z_ref_points, y_ref_points, origin)
end

"""
    VSMWing(name, vsm_aero, vsm_wing, vsm_solver, twist_surfaces, R_b_to_c, pos_cad; transform=nothing)

Construct a `RIGID_DYNAMICS` [`Wing`](@ref) from pre-created VSM objects. Kept for
backward compatibility with predefined structures.
"""
function VSMWing(name, vsm_aero, vsm_wing, vsm_solver,
                 twist_surfaces::AbstractVector,
                 R_b_to_c::AbstractMatrix,
                 pos_cad::AbstractVector;
                 transform=nothing)
    inertia_vec = ones(MVector{3, SimFloat})
    n_twist_surfaces_est = vsm_wing.n_unrefined_sections
    num_aero_outputs = 6 + n_twist_surfaces_est
    num_aero_inputs = 5 + n_twist_surfaces_est
    engine = VSMEngine(vsm_aero, vsm_wing, vsm_solver,
        zeros(SimFloat, num_aero_inputs),
        zeros(SimFloat, num_aero_outputs),
        zeros(SimFloat, num_aero_outputs, num_aero_inputs),
        nothing, nothing, SimFloat(0.0), SimFloat(0.0))
    return Wing(name, twist_surfaces, R_b_to_c, pos_cad, inertia_vec;
        transform, aero=AeroLinearized(engine))
end

"""
    Wing(name, vsm_aero, vsm_wing, vsm_solver, twist_surfaces, R_b_to_c, pos_cad; transform=1)

Backward-compatibility constructor: builds a VSM [`Wing`](@ref) from pre-created
VSM objects (delegates to [`VSMWing`](@ref)).
"""
function Wing(name, vsm_aero, vsm_wing, vsm_solver, twist_surfaces, R_b_to_c,
              pos_cad; kwargs...)
    return VSMWing(name, vsm_aero, vsm_wing, vsm_solver, twist_surfaces, R_b_to_c, pos_cad; kwargs...)
end

# ==================== PLATE WING ==================== #

"""
    PlateWing(name, twist_surfaces, calc_cl, calc_cd;
              dynamics_type=PARTICLE_DYNAMICS, transform=nothing,
              y_damping=150.0, angular_damping=0.0, drag_corr=0.93,
              z_ref_points=nothing, y_ref_points=nothing, origin=nothing)

Construct a flat-plate [`Wing`](@ref) (no VSM engine; `vsm === nothing`). Each
flat-plate section is a 1-point `FIXED` [`TwistSurface`](@ref) carrying the
section's body-frame reference frame, area, and prescribed twist; the shared
polar lookups live on the wing's [`AeroPlate`](@ref) `aero` model. Supports both
`RIGID_DYNAMICS` and `PARTICLE_DYNAMICS`.

# Arguments
- `name`: Wing name/identifier.
- `twist_surfaces`: References (names or indices) to the wing's flat-plate
  sections — each a 1-point `FIXED` [`TwistSurface`](@ref).
- `calc_cl`: CL lookup callable(alpha_deg) → CL.
- `calc_cd`: CD lookup callable(alpha_deg) → CD.

# Keyword Arguments
- `dynamics_type`: `RIGID_DYNAMICS` or `PARTICLE_DYNAMICS` (default).
- `transform`: Reference to transform (name or index).
- `y_damping`, `angular_damping`: Damping coefficients.
- `drag_corr`: Drag correction factor (stored on the `AeroPlate` model).
- `z_ref_points`, `y_ref_points`, `origin`: Body-frame references.
"""
function PlateWing(name, twist_surfaces::AbstractVector,
                   calc_cl, calc_cd;
                   dynamics_type::WingType=PARTICLE_DYNAMICS,
                   transform=nothing,
                   y_damping=150.0,
                   angular_damping=0.0,
                   drag_corr=0.93,
                   z_ref_points=nothing,
                   y_ref_points=nothing,
                   origin=nothing)
    return Wing(name, twist_surfaces, Matrix{SimFloat}(I, 3, 3),
                zeros(KVec3), ones(MVector{3, SimFloat});
                transform, y_damping, angular_damping, dynamics_type,
                aero=AeroPlate(calc_cl, calc_cd; drag_corr),
                z_ref_points, y_ref_points, origin)
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
frame to body frame. After the first step, `refresh_aero!()`
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
