# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# ==================== ENUMS ==================== #

"""
    DynamicsType `DYNAMIC` `STATIC` `BODY_STATIC` `KINEMATIC`

Enumeration for the dynamic model governing a point's motion, a rigid body's
motion, or a station's twist. `DYNAMIC` quantities carry differential state
and are solved by the dynamics; the others are *prescribed* — no state, held
constant within a step, but mutable between steps. A wing's aerodynamic-surface
points are ordinary `DYNAMIC` (particle wing) or `BODY_STATIC` (rigid wing) points;
their wing membership comes from station membership (`point.is_wing_node`).

# Elements
- `DYNAMIC`: Solved by the dynamics. A point moves by Newton's second law; a
  rigid body integrates its 6-DOF state; a station twist solves its
  equilibrium.
- `STATIC`: Prescribed, no state. A point's position is welded in the **world**
  frame; a rigid body is clamped to the world; a station twist is a
  prescribed control input read live via the registered getter.
- `BODY_STATIC`: A point rides a body — static in the rigid
  body's body frame; it feeds its net force (and the moment about the body COM)
  into the body.
- `KINEMATIC`: A twist source whose deflection is prescribed by geometry (e.g. a
  flap hinge between two bodies).
"""
@enum DynamicsType begin
    DYNAMIC
    STATIC
    BODY_STATIC
    KINEMATIC
end

"""
    WingType `RIGID_DYNAMICS` `PARTICLE_DYNAMICS`

Enumeration for the structural representation of a wing.

# Elements
- `RIGID_DYNAMICS`: Wing uses quaternion-based rigid body dynamics, with the
  deformation confined to the twist degrees of freedom of its
  spanwise stations. Aerodynamic forces/moments are applied to the wing
  center of mass; its structural points are `BODY_STATIC` and ride the body.
- `PARTICLE_DYNAMICS`: Wing has no rigid body constraint. Its structural points
  are ordinary `DYNAMIC` particles and the aerodynamic loads are applied to them
  directly, so the wing deforms with the structure.
"""
@enum WingType begin
    RIGID_DYNAMICS
    PARTICLE_DYNAMICS
end

# Backwards-compatible deprecated aliases for the previous WingType names.
Base.@deprecate_binding QUATERNION RIGID_DYNAMICS
Base.@deprecate_binding REFINE PARTICLE_DYNAMICS

"""
    PrincipalFrameMethod

Strategy for computing the principal frame from an inertia tensor.

- `EIGEN_DECOMP`: full 3-axis eigendecomposition + permutation search
  ([`principal_frame`](@ref)). General-purpose; correct for any body.
- `Y_ROTATION`: closed-form rotation about Y only
  ([`calc_inertia_y_rotation`](@ref)). Use for wings symmetric about the
  XZ-plane where the generic permutation search is ambiguous when two
  principal moments are close.
"""
@enum PrincipalFrameMethod begin
    EIGEN_DECOMP
    Y_ROTATION
end

"""
    abstract type WingDynamics

Tag of a wing's dynamics, which a wing's aero component dispatches on
([`wing_dynamics`](@ref)) instead of branching on its `dynamics_type`. Subtypes:
[`RigidDynamics`](@ref), [`ParticleDynamics`](@ref).
"""
abstract type WingDynamics end

"""Dynamics tag for a `RIGID_DYNAMICS` wing (integrated 6-DOF body)."""
struct RigidDynamics <: WingDynamics end

"""Dynamics tag for a `PARTICLE_DYNAMICS` wing (frame fitted from points)."""
struct ParticleDynamics <: WingDynamics end

"""The [`WingDynamics`](@ref) of wing `w` of `sam`, which its aero component dispatches
on."""
wing_dynamics(sam, w) = wing_dynamics_type(sam.definition.wings[w]) == RIGID_DYNAMICS ?
    RigidDynamics() : ParticleDynamics()

"""
    AbstractAeroModel

Supertype for a wing's aerodynamic model. The concrete subtype selects, by
dispatch, the [`aero_component`](@ref) builder that emits the wing's aero
equations. Built-in subtypes: [`AeroNone`](@ref), [`AeroDirect`](@ref),
[`AeroLinearized`](@ref), [`ContinuousAero`](@ref), [`AeroPressure`](@ref),
[`AeroPlate`](@ref). Subtype it and add an
`aero_component` method to plug in custom aerodynamics; see the VSM coupling
documentation for a worked example with a live-updating field.
"""
abstract type AbstractAeroModel end

"""
    UnsteadyAero(; apparent_mass=0.0, wagner=false, wagner_gains=(0.165, 0.335),
                 wagner_rates=(0.0455, 0.3))

Unsteady corrections a wing adds on top of its frozen-circulation VSM forces, held
by the wing's [`VSMEngine`](@ref) and off by default. `apparent_mass` scales the
entrained-air inertia [`apply_apparent_mass!`](@ref) puts on the wing nodes, `1`
being the thin-plate value and `0` disabling it. `wagner` enables the two-state
Wagner lift lag, whose indicial response is `1 - A1·exp(-b1·s) - A2·exp(-b2·s)` in
semi-chords travelled, with gains `(A1, A2)` and rates `(b1, b2)`; the defaults are
R. T. Jones' fit to Wagner's function. The four lag constants are registered
parameters, so retuning them is a sync and not a rebuild; only `wagner` itself is
structural, because it is what adds or removes the two states.

$(TYPEDFIELDS)
"""
mutable struct UnsteadyAero
    "Scale on the thin-plate entrained-air mass; 0 leaves the nodes' inertia alone."
    apparent_mass::SimFloat
    "Whether the wing carries the two-state Wagner lift lag."
    wagner::Bool
    "Wagner indicial gains `(A1, A2)`, the share of lift each lag state withholds."
    wagner_gains::KVec2
    "Wagner indicial rates `(b1, b2)`, per semi-chord travelled."
    wagner_rates::KVec2
end

UnsteadyAero(; apparent_mass=0.0, wagner=false,
             wagner_gains=(0.165, 0.335), wagner_rates=(0.0455, 0.3)) =
    UnsteadyAero(SimFloat(apparent_mass), wagner,
                 KVec2(wagner_gains...), KVec2(wagner_rates...))

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
- `unsteady`: the wing's [`UnsteadyAero`](@ref) corrections.
- `station_le_b`, `station_chord`, `station_y_airf`: the leading edge from the body
  origin [m], chord vector [m] and spanwise axis of the VSM section nearest each of the
  wing's stations, body frame, zero for a station whose row gives its chord.
- `station_sections`, `station_panels`: the unrefined sections and the panels each
  station takes ([`compute_spatial_station_mapping!`](@ref)).
- `station_body_mass`: each station's share of the mass of the wing's body [kg].

The station vectors are indexed by a station's position among its wing's stations.
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
    unsteady::UnsteadyAero
    station_le_b::Vector{KVec3}
    station_chord::Vector{KVec3}
    station_y_airf::Vector{KVec3}
    station_sections::Vector{Vector{Int}}
    station_panels::Vector{Vector{Int}}
    station_body_mass::Vector{SimFloat}
end


"""
    abstract type AbstractVSMAero <: AbstractAeroModel

Aero modes backed by a [`VSMEngine`](@ref) (stored in the `engine` field). VSM
operations dispatch on this supertype; the engine's fields (`vsm_wing`, `aero_x`,
…) are forwarded from the mode. Built-in subtypes: [`AeroDirect`](@ref),
[`AeroLinearized`](@ref), [`ContinuousAero`](@ref) and [`AeroPressure`](@ref).
"""
abstract type AbstractVSMAero <: AbstractAeroModel end

"""
    station_vsm_value(mode, field, k)

What the aero mode `mode` of a wing worked out of the wing's `k`-th station: the `k`-th
of its VSM engine's `field` ([`VSMEngine`](@ref)), `nothing` for a mode without an
engine or a station it worked out nothing of.
"""
station_vsm_value(::AbstractAeroModel, field, k) = nothing
function station_vsm_value(mode::AbstractVSMAero, field, k)
    values = getfield(getfield(mode, :engine), field)
    return k <= length(values) ? values[k] : nothing
end

# VSM engine fields forwarded from a VSM aero mode to its `engine`.
const VSM_ENGINE_FIELDS = (
    :vsm_aero, :vsm_wing, :vsm_solver, :aero_y, :aero_x, :aero_jac,
    :point_to_vsm_point, :wing_segments, :aero_scale_chord, :aero_z_offset,
    :unsteady)

"""
    AeroSpec{M}

The VSM aero mode `M` a structure's wing is flown with, and the options `M` is built with:
what a structure holds of a mode a model builds with its engine ([`model_aero`](@ref)).
"""
struct AeroSpec{M, O <: NamedTuple} <: AbstractAeroModel
    options::O
end
AeroSpec{M}(; options...) where {M} =
    AeroSpec{M, typeof(values(options))}(values(options))

"""
    model_aero(mode, wing, n_stations, set, vsm_set)

The aero mode a model flies the wing whose `wings` row is `wing`, with `n_stations`
stations, with: `mode`, or for an [`AeroSpec`](@ref) its mode with an engine of its own
built from `set` and `vsm_set`, unless there is no `set` to build it from.
"""
model_aero(mode::AbstractAeroModel, wing, n_stations, set, vsm_set) = mode
function model_aero(spec::AeroSpec{M}, wing, n_stations, set, vsm_set) where {M}
    isnothing(set) && return spec
    isnothing(vsm_set) && error("Wing $(wing.name): aero mode $(nameof(M)) needs VSM " *
        "geometry, but the model has no `vsm_set` and the project names no `vsm_settings`.")
    options = spec.options
    unsteady = haskey(options, :unsteady) ? copy_unsteady(options.unsteady) : UnsteadyAero()
    return M(build_vsm_engine(set, vsm_set, wing_dynamics_type(wing); n_stations,
        aero_scale_chord=something(wing.model.aero_scale_chord, 0.0),
        aero_z_offset=something(wing.model.aero_z_offset, 0.0), unsteady);
        Base.structdiff(options, NamedTuple{(:unsteady,)})...)
end

"""A new [`UnsteadyAero`](@ref) with the values of `unsteady`, so a model's engine shares
none with another's."""
copy_unsteady(unsteady::UnsteadyAero) =
    UnsteadyAero((copy(getfield(unsteady, name)) for name in fieldnames(UnsteadyAero))...)

"""The aero mode wing row `wing` is flown with: the one `aero` gives under its name, else
the one its `aero_mode` column names, else `AeroLinearized` for a rigid wing and
`AeroDirect` for a particle wing."""
function wing_aero_spec(wing, aero)
    haskey(aero, wing.name) && return aero[wing.name]
    isnothing(wing.model.aero_mode) || return parse_aero_mode(wing.model.aero_mode)
    return wing_dynamics_type(wing) == RIGID_DYNAMICS ? AeroLinearized() : AeroDirect()
end

"""The [`WingType`](@ref) wing row `wing`, or its model, names, `RIGID_DYNAMICS` where it
names none."""
wing_dynamics_type(wing) = wing_dynamics_type(wing.model)
wing_dynamics_type(model::WingModel) =
    parse_wing_type(something(model.dynamics_type, "RIGID_DYNAMICS"))

function Base.getproperty(mode::AbstractVSMAero, sym::Symbol)
    sym in VSM_ENGINE_FIELDS && return getproperty(getfield(mode, :engine), sym)
    return getfield(mode, sym)
end

function Base.setproperty!(mode::AbstractVSMAero, sym::Symbol, value)
    sym in VSM_ENGINE_FIELDS && return setproperty!(getfield(mode, :engine), sym, value)
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

# ==================== POINT ==================== #

"""A value one per axis: a vector column as it is, one number on every axis, zero unset."""
per_axis(values::Vector{Float64}) = values
per_axis(value::Float64) = SVector(value, value, value)
per_axis(::Nothing) = zero(SVector{3, Float64})

"""The [`DynamicsType`](@ref) the `type` column of `row` names, by KiteGeometry's
`parse_dynamics_type`, refusing one it leaves unset."""
function row_dynamics_type(row)
    isnothing(row.type) && error("$(row.name) needs its `type`")
    name = KiteGeometry.parse_dynamics_type(row.type)
    return only(type for type in instances(DynamicsType) if string(type) == name)
end

"""The dynamics type of the twist of `station`, from its `type` column."""
station_type(station) = row_dynamics_type(station)

"""The dynamics type of `pulley`, from its `type` column."""
pulley_type(pulley) = row_dynamics_type(pulley)

"""The dynamics type of `point`: `STATIC` where it is held, `BODY_STATIC` where it rides a
body or a tube, `DYNAMIC` otherwise."""
function point_type(point)
    point.fixed && return STATIC
    return isnothing(point.body) && isnothing(point.tube) ? DYNAMIC : BODY_STATIC
end

"""A per-axis damping as a column of one value per axis, zero for `nothing`."""
axis_column(::Nothing) = zeros(3)
axis_column(value::Real) = fill(Float64(value), 3)
axis_column(values::AbstractVector) = Vector{Float64}(values)

# ==================== SEGMENT ==================== #

"""
    segment_unit_damping(row, material, set)

The damping per unit length [N·s] of segment `row` of `material` in a run with settings
`set`: its `unit_damping`, else its `damping_per_stiffness` times its stiffness, else
`set.rel_damping` times its stiffness.
"""
function segment_unit_damping(row::KiteGeometry.Segment{SegmentSpring}, material, set)
    (; name, unit_damping, damping_per_stiffness) = row
    material.law == "linear" || haskey(material.parameters, "unit_stiffness") ||
        error("segment $name: material law $(material.law) has no force law here; " *
              "give its unit_stiffness in Julia")
    isnothing(unit_damping) || isnothing(damping_per_stiffness) ||
        error("segment $name: give either `unit_damping` or `damping_per_stiffness`, " *
              "not both")
    isnothing(unit_damping) || return unit_damping
    stiffness = material_unit_stiffness(material, row.diameter)
    linear = stiffness isa Real
    if !isnothing(damping_per_stiffness)
        linear || error("segment $name: `damping_per_stiffness` needs a linear " *
                        "`unit_stiffness`")
        return damping_per_stiffness * stiffness
    end
    linear || error("segment $name: give `unit_damping` for a nonlinear `unit_stiffness`")
    isnothing(set) && return NaN
    set.rel_damping == 0 && @warn "segment $name: unit_damping is zero (no rel_damping " *
                                  "in the settings)."
    return set.rel_damping * stiffness
end

"""
    material_unit_stiffness(material, diameter, read_held=identity)

The stiffness per unit length of `material` at `diameter` [m]: KiteGeometry's
`unit_stiffness` for a `linear` material [N], else the number or force law `F(ε)` it holds,
read through `read_held`.
"""
@inline function material_unit_stiffness(material, diameter,
                                         read_held::F=identity) where {F}
    stiffness = KiteGeometry.unit_stiffness(material, diameter)
    return isnothing(stiffness) ? read_held(material.parameters["unit_stiffness"]) :
           stiffness
end

"""`value` as a `Float64`."""
to_float64(value) = value isa Float64 ? value : Float64(value)::Float64


# ==================== WINCH ==================== #

"""
    abstract type AbstractWinchModel

Selects the winch motor dynamics. Each concrete model carries its own
parameter fields and adds a `winch_component(model, sam, idx; name,
params)` method building the MTK subsystem (mirrors [`AbstractAeroModel`](@ref)).
The model is a winch row's `model.dynamics` ([`WinchModel`](@ref)); the common drum
parameters (`gear_ratio`, `drum_radius`, `coulomb_friction`, `viscous_coefficient`,
`inertia_total`) are the row's and its `WinchModel`'s. New models = new struct + a few
methods; see [`TorqueWinch`](@ref).
"""
abstract type AbstractWinchModel end

