# Copyright (c) 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

"""
    PointModel

What SymbolicAWEModels holds of a point beyond its KiteGeometry `points` row: the columns
only it reads, each `nothing` where the row leaves it unset.

$(TYPEDFIELDS)
"""
@kwdef mutable struct PointModel <: KiteGeometry.AbstractModel
    "Damping per unit mass on the body axes [1/s], one value or one per axis."
    body_frame_damping::Union{Nothing, Float64, Vector{Float64}} = nothing
    "Damping per unit mass on the world axes [1/s], one value or one per axis."
    world_frame_damping::Union{Nothing, Float64, Vector{Float64}} = nothing
    "Drag area [m²]."
    area::Union{Nothing, Float64} = nothing
    "Drag coefficient."
    drag_coeff::Union{Nothing, Float64} = nothing
    "Whether the point is held on a sphere."
    fix_sphere::Union{Nothing, Bool} = nothing
    "Whether the point is frozen in place."
    fix_static::Union{Nothing, Bool} = nothing
    "Where the point sits in its body's frame [m]."
    anchor_b::Union{Nothing, Vector{Float64}} = nothing
    "Entrained-air mass resisting acceleration without adding weight [kg]."
    apparent_mass::Union{Nothing, Float64} = nothing
    "The `wings` row whose frame a bridle or pulley point is damped in."
    wing_idx::Union{Nothing, KiteGeometry.NameRef} = nothing
end

KiteGeometry.model_columns(::Type{PointModel}) = fieldnames(PointModel)
KiteGeometry.model_references(::Type{PointModel}) = (wing_idx=:wings,)

"""
    SegmentSpring

The spring of a segment, the model of a `segments` row: its damping and how it carries
compression as a document authors them, each `nothing` where the row leaves it unset.
A run resolves its damping with [`segment_unit_damping`](@ref).

$(TYPEDFIELDS)
"""
@kwdef mutable struct SegmentSpring <: KiteGeometry.AbstractModel
    "Damping per unit length [N·s]."
    unit_damping::Union{Nothing, Float64} = nothing
    "Damping per unit length over stiffness per unit length [s]."
    damping_per_stiffness::Union{Nothing, Float64} = nothing
    "Compressive over tensile stiffness (0-1)."
    compression_frac::Union{Nothing, Float64} = nothing
    "Fraction of the damping still acting under compression (0-1)."
    compression_damping_frac::Union{Nothing, Float64} = nothing
end

"""
    StationModel

What SymbolicAWEModels holds of a station beyond its KiteGeometry `spanwise_stations` row:
the columns only it reads, each `nothing` where the row leaves it unset.

$(TYPEDFIELDS)
"""
@kwdef mutable struct StationModel <: KiteGeometry.AbstractModel
    "Dynamics type of the station's twist, such as `DYNAMIC` or `STATIC`."
    type::Union{Nothing, String} = nothing
    "Chordwise twist axis as a fraction of the chord (0 leading edge, 1 trailing edge)."
    moment_frac::Union{Nothing, Float64} = nothing
    "Twist damping [N·m·s/rad]."
    damping::Union{Nothing, Float64} = nothing
    "Twist restoring stiffness [N·m/rad]."
    stiffness::Union{Nothing, Float64} = nothing
    "The bodies the station twists with, by name or index."
    bodies::Union{Nothing, Vector{KiteGeometry.NameRef}} = nothing
    "The main and the flap body of a body flap, by name or index."
    flap_bodies::Union{Nothing, Vector{KiteGeometry.NameRef}} = nothing
    "The three points of a point flap, by name or index."
    flap_points::Union{Nothing, Vector{KiteGeometry.NameRef}} = nothing
    "Unit flap-hinge axis."
    flap_axis::Union{Nothing, Vector{Float64}} = nothing
    "Chord vector of a flat-plate section, its wing's frame [m]."
    chord::Union{Nothing, Vector{Float64}} = nothing
    "Spanwise vector of a flat-plate section, its local panel frame."
    y_airf::Union{Nothing, Vector{Float64}} = nothing
    "Surface area of a flat-plate section [m²]."
    area::Union{Nothing, Float64} = nothing
    "Twist angle a run starts from, held for a `STATIC` station [rad]."
    twist::Union{Nothing, Float64} = nothing
    "Reference chord directions `[main, flap]` of a body flap, each body's frame."
    flap_chord_refs::Union{Nothing, Vector{Vector{Float64}}} = nothing
end

KiteGeometry.model_columns(::Type{StationModel}) = fieldnames(StationModel)
KiteGeometry.model_references(::Type{StationModel}) =
    (bodies=:bodies, flap_bodies=:bodies, flap_points=:points)

"""
    WingModel

The columns of a KiteGeometry `wings` row only SymbolicAWEModels reads, each `nothing`
where the row leaves it unset.

$(TYPEDFIELDS)
"""
@kwdef struct WingModel <: KiteGeometry.AbstractModel
    "Name of the aerodynamic model, such as `AeroLinearized`."
    aero_mode::Union{Nothing, String} = nothing
    "Damping of the spin relative to the wing's `small_earth_frame`, about the body axes
    [1/s], one value or one per axis."
    angular_damping::Union{Nothing, Float64, Vector{Float64}} = nothing
    "Scale on the VSM panel chord."
    aero_scale_chord::Union{Nothing, Float64} = nothing
    "Offset of the VSM panels along the body z-axis [m]."
    aero_z_offset::Union{Nothing, Float64} = nothing
    "How the principal frame is found: `EIGEN_DECOMP` or `Y_ROTATION`."
    principal_frame_method::Union{Nothing, String} = nothing
    "`PARTICLE_DYNAMICS` or `RIGID_DYNAMICS`."
    dynamics_type::Union{Nothing, String} = nothing
    "The points, weighted, whose position is the wing's origin."
    origin_idx::Any = nothing
    "The two weighted point groups the wing's z-axis runs between."
    z_ref_points::Any = nothing
    "The two weighted point groups the wing's y-axis runs between."
    y_ref_points::Any = nothing
    "Origin of the wing's body frame in the authoring frame [m]."
    pos_cad::Union{Nothing, Vector{Float64}} = nothing
    "Centre of the wing's own mass in the authoring frame [m]."
    com::Union{Nothing, Vector{Float64}} = nothing
    "Inertia per unit mass [m²], a 3×3 matrix or `[Ixx, Iyy, Izz, Ixy, Ixz, Iyz]`."
    unit_inertia::Any = nothing
    "Drag correction of a flat-plate wing."
    drag_corr::Union{Nothing, Float64} = nothing
    "The stations the wing lists, by name or index."
    stations::Union{Nothing, Vector{KiteGeometry.NameRef}} = nothing
    "Whether the points of the wing's stations add their moment to it."
    group_points_moment::Union{Nothing, Bool} = nothing
end

"""
    PulleyModel

What SymbolicAWEModels holds of a pulley beyond its KiteGeometry `pulleys` row: the columns
only it reads, each `nothing` where the row leaves it unset.

$(TYPEDFIELDS)
"""
@kwdef mutable struct PulleyModel <: KiteGeometry.AbstractModel
    "Dynamics type, such as `DYNAMIC` or `STATIC`."
    type::Union{Nothing, String} = nothing
    "Fraction of line tension the sheave passes on (0-1)."
    efficiency::Union{Nothing, Float64} = nothing
    "Damping on rope travel [N·s/m]."
    damping::Union{Nothing, Float64} = nothing
    "Whether the rope split is frozen."
    brake::Union{Nothing, Bool} = nothing
    "Friction smoothing width [m/s]."
    friction_epsilon::Union{Nothing, Float64} = nothing
end

KiteGeometry.model_columns(::Type{PulleyModel}) =
    (:type, :efficiency, :damping, :brake, :friction_epsilon)

"""
    WinchModel{W}

What SymbolicAWEModels holds of a winch beyond its KiteGeometry `winches` row: the columns
only it reads, each `nothing` where the row leaves it unset, and the motor dynamics `W` (an
[`AbstractWinchModel`](@ref)).
A property it has no field for is the dynamics'; `friction_epsilon` is the dynamics' where
they have one.

$(TYPEDFIELDS)
"""
@kwdef mutable struct WinchModel{W} <: KiteGeometry.AbstractModel
    "Reel-out speed at the start [m/s]."
    init_vel::Union{Nothing, Float64} = nothing
    "Brake setting."
    brake::Union{Nothing, Float64} = nothing
    "Whether the winch follows a set speed rather than a set torque."
    speed_controlled::Union{Nothing, Bool} = nothing
    "Friction smoothing width [m/s]."
    friction_epsilon::Union{Nothing, Float64} = nothing
    "Coulomb friction force [N]."
    coulomb_friction::Union{Nothing, Float64} = nothing
    "Viscous friction coefficient [N·s/m]."
    viscous_coefficient::Union{Nothing, Float64} = nothing
    "Total rotational inertia [kg·m²]."
    inertia_total::Union{Nothing, Float64} = nothing
    "Motor dynamics, a [`TorqueWinch`](@ref) by default."
    dynamics::W = torque_winch(friction_epsilon)
end

WinchModel(init_vel, brake, speed_controlled, friction_epsilon, coulomb_friction,
           viscous_coefficient, inertia_total) =
    WinchModel(; init_vel, brake, speed_controlled, friction_epsilon, coulomb_friction,
               viscous_coefficient, inertia_total)

KiteGeometry.model_columns(::Type{<:WinchModel}) =
    (:init_vel, :brake, :speed_controlled, :friction_epsilon, :coulomb_friction,
     :viscous_coefficient, :inertia_total)

@inline function Base.getproperty(model::WinchModel{W}, property::Symbol) where {W}
    (hasfield(W, property) || !hasfield(WinchModel{W}, property)) &&
        return getproperty(getfield(model, :dynamics), property)
    return getfield(model, property)
end

function Base.setproperty!(model::WinchModel{W}, property::Symbol, value) where {W}
    hasfield(WinchModel{W}, property) &&
        setfield!(model, property, convert(fieldtype(typeof(model), property), value))
    (hasfield(W, property) || !hasfield(WinchModel{W}, property)) &&
        setproperty!(getfield(model, :dynamics), property, value)
    return value
end

"""
    TimoshenkoTubeModel

The columns of a KiteGeometry `tubes` row with model `timoshenko`, each `nothing` where the
row leaves it unset: a 2-node corotational Timoshenko beam element, whose stiffness couples
each node's transverse displacement to its rotation (shear reduction
`Φ = 12·EI/(k·GA·L²)`); a chain of tubes forms a beam. A rigidity is a number or, in a
structure built in Julia, a callable of its strain or curvature returning the effective
rigidity there; one left unset follows from the tube's material.

$(TYPEDFIELDS)
"""
@kwdef mutable struct TimoshenkoTubeModel <: KiteGeometry.AbstractModel
    "Axial rigidity [N]."
    EA::Any = nothing
    "Shear rigidity [N]."
    GA::Any = nothing
    "Torsional rigidity [N·m²]."
    GJ::Any = nothing
    "Bending rigidity about the y-axis [N·m²]."
    EIy::Any = nothing
    "Bending rigidity about the z-axis [N·m²]."
    EIz::Any = nothing
    "Shear correction coefficient, 5/6 when unset."
    shear_coeff::Union{Nothing, Float64} = nothing
    "Rayleigh damping, damping over stiffness [s]; 0 when unset."
    damping::Union{Nothing, Float64} = nothing
    "Tube end in the first body's frame, from its frame origin [m]."
    anchor_a::Union{Nothing, Vector{Float64}} = nothing
    "Tube end in the second body's frame, from its frame origin [m]."
    anchor_b::Union{Nothing, Vector{Float64}} = nothing
end
KiteGeometry.model_columns(::Type{TimoshenkoTubeModel}) =
    (:EA, :GA, :GJ, :EIy, :EIz, :shear_coeff, :damping, :anchor_a, :anchor_b)


"""
    ElasticTubeModel

The columns of a KiteGeometry `tubes` row with model `elastic`, each `nothing` where the
row leaves it unset: a lumped 6-DOF spring between the tube's two ends, their relative pose
in the first body's frame split into axial, shear, torsion and bending springs. A
stiffness is a number (force `k·Δ`) or, in a structure built in Julia, a callable of its
displacement or angle (force `f(Δ)`); every stiffness must be set.

$(TYPEDFIELDS)
"""
@kwdef mutable struct ElasticTubeModel <: KiteGeometry.AbstractModel
    "Axial stiffness [N/m]."
    stiffness_axial::Any = nothing
    "Shear stiffness [N/m]."
    stiffness_shear::Any = nothing
    "Torsional stiffness [N·m/rad]."
    stiffness_torsion::Any = nothing
    "Bending stiffness [N·m/rad]."
    stiffness_bending::Any = nothing
    "Rayleigh damping, damping over stiffness [s]; 0 when unset."
    damping::Union{Nothing, Float64} = nothing
    "Tube end in the first body's frame, from its frame origin [m]."
    anchor_a::Union{Nothing, Vector{Float64}} = nothing
    "Tube end in the second body's frame, from its frame origin [m]."
    anchor_b::Union{Nothing, Vector{Float64}} = nothing
end
KiteGeometry.model_columns(::Type{ElasticTubeModel}) =
    (:stiffness_axial, :stiffness_shear, :stiffness_torsion, :stiffness_bending,
     :damping, :anchor_a, :anchor_b)

"""
    rigidity_fields(model) -> Tuple

The columns of the tube row model `model` holding a rigidity or stiffness.
"""
rigidity_fields(::TimoshenkoTubeModel) = (:EA, :GA, :GJ, :EIy, :EIz)
rigidity_fields(::ElasticTubeModel) =
    (:stiffness_axial, :stiffness_shear, :stiffness_torsion, :stiffness_bending)


KiteGeometry.model_references(::Type{WingModel}) = (
    origin_idx=KiteGeometry.WeightedRefs(:points),
    z_ref_points=(KiteGeometry.WeightedRefs(:points), KiteGeometry.WeightedRefs(:points)),
    y_ref_points=(KiteGeometry.WeightedRefs(:points), KiteGeometry.WeightedRefs(:points)),
    stations=:spanwise_stations)

"""What frames and weighs the body of a wing whose row carries `model`, for KiteGeometry's
`complete_wings!`: its dynamics type, reference points and mesh mass properties."""
KiteGeometry.wing_spec(model::WingModel) = KiteGeometry.WingSpec(;
    particle=wing_dynamics_type(model) == PARTICLE_DYNAMICS, origin=model.origin_idx,
    z_ref_points=model.z_ref_points, y_ref_points=model.y_ref_points, com=model.com,
    unit_inertia=unit_inertia_tensor(model.unit_inertia))

"""Per-unit-mass inertia [m²] as a 3×3 matrix, given as one or as `[Ixx, Iyy, Izz, Ixy,
Ixz, Iyz]`; `nothing` for none."""
unit_inertia_tensor(::Nothing) = nothing
unit_inertia_tensor(inertia::AbstractMatrix) = inertia
unit_inertia_tensor(inertia::AbstractVector) = first(inertia) isa AbstractVector ?
    reduce(vcat, permutedims.(inertia)) : ObjAdapter.unit_inertia_matrix(inertia)

"""
    register_definition_models()

Register SymbolicAWEModels' models of each KiteGeometry block, so `load_authoring` and
`load_structure` read its columns into them: the default model of points, segments,
stations, wings, pulleys, bodies and winches, and the tube models `timoshenko`, the
default, and `elastic`.
"""
function register_definition_models()
    register_model!(:points, "sam_point", PointModel,
                    ("N*s/m", "N*s/m", "m^2", "-", "-", "-", "m", "kg", "-");
                    default=true)
    register_model!(:segments, "spring", SegmentSpring, ("N*s", "s", "-", "-");
                    default=true)
    register_model!(:spanwise_stations, "sam_station", StationModel,
                    ("-", "-", "N*m*s/rad", "N*m/rad", "-", "-", "-", "-", "m", "-", "m^2",
                     "rad", "-");
                    default=true)
    register_model!(:wings, "sam_wing", WingModel,
                    ("-", "1/s", "-", "m", "-", "-", "-", "-", "-", "m", "m", "m^2", "-",
                     "-", "-");
                    default=true)
    register_model!(:pulleys, "sam_pulley", PulleyModel,
                    ("-", "-", "N*s/m", "-", "m/s"); default=true)
    register_model!(:bodies, "sam_body", BodyModel,
                    ("-", "1/s", "1/s", "1/s", "-", "-", "N", "N", "N*m", "-", "kg", "m");
                    default=true)
    register_model!(:winches, "sam_winch", WinchModel,
                    ("m/s", "-", "-", "m/s", "N", "N*s/m", "kg*m^2");
                    default=true)
    register_model!(:tubes, "timoshenko", TimoshenkoTubeModel,
                    ("N", "N", "N*m^2", "N*m^2", "N*m^2", "-", "s", "m", "m");
                    default=true)
    register_model!(:tubes, "elastic", ElasticTubeModel,
                    ("N/m", "N/m", "N*m/rad", "N*m/rad", "s", "m", "m"))
    return nothing
end
