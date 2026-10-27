# Copyright (c) 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

"""
    PointModel

What SymbolicAWEModels holds of a point beyond its KiteGeometry `points` row: the columns
only it reads, each `nothing` where the row leaves it unset, then what it derives when it
builds the structure, which no document holds.

$(TYPEDFIELDS)
"""
@kwdef mutable struct PointModel <: KiteGeometry.AbstractModel
    "Dynamics type, such as `DYNAMIC`, `STATIC` or `BODY_STATIC`."
    type::Union{Nothing, String} = nothing
    "Damping in the body frame [N·s/m], one value or one per axis."
    body_frame_damping::Union{Nothing, Float64, Vector{Float64}} = nothing
    "Damping in the world frame [N·s/m], one value or one per axis."
    world_frame_damping::Union{Nothing, Float64, Vector{Float64}} = nothing
    "Drag area [m²]."
    area::Union{Nothing, Float64} = nothing
    "Drag coefficient."
    drag_coeff::Union{Nothing, Float64} = nothing
    "The tube the point rides on, by name or index."
    tube::Union{Nothing, KiteGeometry.NameRef} = nothing
    "Whether the point is held on a sphere."
    fix_sphere::Union{Nothing, Bool} = nothing
    "Whether the point is frozen in place."
    fix_static::Union{Nothing, Bool} = nothing
    "The wing the point belongs to, by name or index."
    wing_idx::Union{Nothing, KiteGeometry.NameRef} = nothing
    "Where the point sits in its body's frame [m]."
    anchor_b::Union{Nothing, Vector{Float64}} = nothing
    "Index in the structure's points."
    idx::Int64 = 0
    "The transform placing the point, by name or index; 0 for none."
    transform_ref::KiteGeometry.NameRef = 0
    "Index of the transform placing the point; 0 for none."
    transform_idx::Int64 = 0
    "Index of the wing the point belongs to; 0 for none."
    resolved_wing_idx::Int64 = 0
    "Index of the body the point rides; 0 for none."
    resolved_body_idx::Int64 = 0
    "Index of the tube the point rides; 0 for none."
    resolved_tube_idx::Int64 = 0
    "Where the point sits in the frame of the body it rides [m]: `anchor_b`, else where its
    position puts it."
    resolved_anchor_b::KVec3 = zeros(KVec3)
    "Undeformed position relative to the wing's centre of mass, principal frame [m]."
    pos_undeformed_b::KVec3 = zeros(KVec3)
    "External force on the point, world frame [N]."
    ext_force_w::KVec3 = zeros(KVec3)
    "Own mass plus half of each of its segments' masses at rest length [kg]."
    total_mass::SimFloat = 0.0
    "Entrained-air mass resisting acceleration without adding weight [kg]."
    apparent_mass::SimFloat = 0.0
    "Whether the point is a structural node of a wing: a member of one of its stations."
    is_wing_node::Bool = false
    "Position `s ∈ [0,1]` along the beam element of the tube the point rides."
    beam_frac::SimFloat = 0.0
    "Offset off that beam's centerline, rest element frame [m]."
    beam_offset_b::KVec3 = zeros(KVec3)
end

PointModel(type, body_frame_damping, world_frame_damping, area, drag_coeff, tube,
           fix_sphere, fix_static, wing_idx, anchor_b) =
    PointModel(; type, body_frame_damping, world_frame_damping, area, drag_coeff, tube,
               fix_sphere, fix_static, wing_idx, anchor_b,
               resolved_anchor_b=KVec3(something(anchor_b, zeros(3))))

KiteGeometry.model_columns(::Type{PointModel}) =
    (:type, :body_frame_damping, :world_frame_damping, :area, :drag_coeff, :tube,
     :fix_sphere, :fix_static, :wing_idx, :anchor_b)

"""
    SegmentModel

What SymbolicAWEModels holds of a segment beyond its KiteGeometry `segments` row: the
columns of KiteGeometry's segment model `spring`, then what it derives when it builds the
structure, which no document holds.

$(TYPEDFIELDS)
"""
@kwdef mutable struct SegmentModel <: KiteGeometry.AbstractModel
    "Damping per unit length [N·s]."
    unit_damping::Float64
    "Compressive over tensile stiffness; 0 for a segment that goes slack."
    compression_frac::Float64 = 0.1
    "Fraction of `unit_damping` still acting under compression."
    compression_damping_frac::Float64 = 1.0
    "Index in the structure's segments."
    idx::Int64 = 0
    "Stiffness per unit length: a number [N], or a callable force law `F(ε)` [N]."
    unit_stiffness::Any = NaN
    "Material density [kg/m³]."
    density::SimFloat = NaN
end

SegmentModel(unit_damping, compression_frac, compression_damping_frac) =
    SegmentModel(; unit_damping, compression_frac, compression_damping_frac)

KiteGeometry.model_columns(::Type{SegmentModel}) =
    (:unit_damping, :compression_frac, :compression_damping_frac)

"""
    StationModel

What SymbolicAWEModels holds of a station beyond its KiteGeometry `spanwise_stations` row:
the columns only it reads, each `nothing` where the row leaves it unset, then what it
derives when it builds the structure, which no document holds.

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
    "Index in the structure's stations."
    idx::Int64 = 0
    "Indices of the station's points."
    point_idxs::Vector{Int64} = Int64[]
    "Index of the body whose twist the station carries; 0 for none."
    resolved_wing_idx::Int64 = 0
    "Indices of the bodies the station twists with."
    body_idxs::Vector{Int64} = Int64[]
    "Indices of the main and the flap body of a body flap."
    flap_body_idxs::Vector{Int64} = Int64[]
    "Indices of the three points of a point flap."
    flap_point_idxs::Vector{Int64} = Int64[]
    "Leading edge position, body frame [m]."
    le_pos::KVec3 = zeros(KVec3)
    "Chord vector, body frame [m]."
    chord::KVec3 = zeros(KVec3)
    "Spanwise vector, local panel frame."
    y_airf::KVec3 = zeros(KVec3)
    "Twist angle a run starts from, held for a `STATIC` station [rad]."
    twist::SimFloat = 0.0
    "Indices of the unrefined VSM sections of the station."
    unrefined_section_idxs::Vector{Int64} = Int64[]
    "Indices of the VSM panels whose twist moment the station takes."
    panel_idxs::Vector{Int64} = Int64[]
    "Share of its wing's mass on no point, by the area of its panels [kg]."
    body_mass::SimFloat = 0.0
    "Surface area of a flat-plate section [m²]; `NaN` unused."
    area::SimFloat = NaN
    "Reference chord directions `[main, flap]` of a flap, each body's frame."
    flap_chord_refs::Vector{KVec3} = KVec3[]
    "Deflection of the flap in its initial pose [rad]."
    flap_rest_delta::SimFloat = 0.0
end

StationModel(type, moment_frac, damping, stiffness, bodies, flap_bodies, flap_points,
             flap_axis) =
    StationModel(; type, moment_frac, damping, stiffness, bodies, flap_bodies,
                 flap_points, flap_axis)

KiteGeometry.model_columns(::Type{StationModel}) =
    (:type, :moment_frac, :damping, :stiffness, :bodies, :flap_bodies, :flap_points,
     :flap_axis)

"""
    WingModel

The columns of a KiteGeometry `wings` row only SymbolicAWEModels reads, each `nothing`
where the row leaves it unset.

$(TYPEDFIELDS)
"""
struct WingModel <: KiteGeometry.AbstractModel
    "Name of the aerodynamic model, such as `AeroLinearized`."
    aero_mode::Union{Nothing, String}
    "Spin damping about the body axes [1/s], one value or one per axis."
    angular_damping::Union{Nothing, Float64, Vector{Float64}}
    "Scale on the VSM panel chord."
    aero_scale_chord::Union{Nothing, Float64}
    "Offset of the VSM panels along the body z-axis [m]."
    aero_z_offset::Union{Nothing, Float64}
    "How the principal frame is found: `EIGEN_DECOMP` or `Y_ROTATION`."
    principal_frame_method::Union{Nothing, String}
    "`PARTICLE_DYNAMICS` or `RIGID_DYNAMICS`."
    dynamics_type::Union{Nothing, String}
    "The points, weighted, whose position is the wing's origin."
    origin_idx::Any
    "The two weighted point groups the wing's z-axis runs between."
    z_ref_points::Any
    "The two weighted point groups the wing's y-axis runs between."
    y_ref_points::Any
    "Origin of the wing's body frame in the authoring frame [m]."
    pos_cad::Union{Nothing, Vector{Float64}}
    "Centre of the wing's own mass in the authoring frame [m]."
    com::Union{Nothing, Vector{Float64}}
    "Inertia per unit mass [m²], a 3×3 matrix or `[Ixx, Iyy, Izz, Ixy, Ixz, Iyz]`."
    unit_inertia::Any
    "Drag correction of a flat-plate wing."
    drag_corr::Union{Nothing, Float64}
    "The stations the wing lists, by name or index."
    stations::Union{Nothing, Vector{KiteGeometry.NameRef}}
end

"""
    PulleyModel

What SymbolicAWEModels holds of a pulley beyond its KiteGeometry `pulleys` row: the columns
only it reads, each `nothing` where the row leaves it unset, then what it derives when it
builds the structure, which no document holds.

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
    "Index in the structure's pulleys."
    idx::Int64 = 0
    "Sum of the lengths of the pulley's two segments [m]."
    sum_len::SimFloat = 0.0
end

PulleyModel(type, efficiency, damping, brake, friction_epsilon) =
    PulleyModel(; type, efficiency, damping, brake, friction_epsilon)

KiteGeometry.model_columns(::Type{PulleyModel}) =
    (:type, :efficiency, :damping, :brake, :friction_epsilon)

"""
    TetherModel

What SymbolicAWEModels holds of a tether beyond its KiteGeometry `tethers` row: the columns
only it reads, each `nothing` where the row leaves it unset, then what it derives when it
builds the structure, which no document holds.

$(TYPEDFIELDS)
"""
@kwdef mutable struct TetherModel <: KiteGeometry.AbstractModel
    "Length the tether is stretched to when placed [m]."
    init_stretched_length::Union{Nothing, Float64} = nothing
    "Tension the placed tether carries, which sets its rest length [N]."
    init_tether_force::Union{Nothing, Float64} = nothing
    "Strain of the placed tether, which sets its rest length."
    init_stretch_frac::Union{Nothing, Float64} = nothing
    "Index in the structure's tethers."
    idx::Int64 = 0
    "Indices of the tether's segments."
    segment_idxs::Vector{Int64} = Int64[]
    "Index of the tether's start point; 0 for none."
    start_point_idx::Int64 = 0
    "Index of the tether's end point; 0 for none."
    end_point_idx::Int64 = 0
    "Number of segments a tether given by its end points is split into."
    n_segments::Int64 = 0
    "Stiffness per unit length of the segments it generates: a number [N] or a callable
    force law `F(ε)` [N]; `NaN` from the settings."
    unit_stiffness::Any = NaN
    "Damping per unit length of the segments it generates [N·s]; `NaN` from the settings."
    unit_damping::SimFloat = NaN
    "Diameter [m]; `NaN` from the settings."
    diameter::SimFloat = NaN
    "Material density [kg/m³]; `NaN` from the settings."
    density::SimFloat = NaN
    "Young's modulus [Pa], the diameter-independent form of `unit_stiffness`."
    youngs_modulus::SimFloat = NaN
    "Damping per unit stiffness [s], the diameter-independent form of `unit_damping`."
    damping_per_stiffness::SimFloat = NaN
    "Compressive over tensile stiffness of the segments it generates."
    compression_frac::SimFloat = 0.1
    "Fraction of `unit_damping` still acting under compression, per segment."
    compression_damping_frac::SimFloat = 1.0
    "Unstretched length where [`place!`](@ref) puts it: the sum of its segments' `l0` [m]."
    len::SimFloat = 0.0
end

TetherModel(init_stretched_length, init_tether_force, init_stretch_frac) =
    TetherModel(; init_stretched_length, init_tether_force, init_stretch_frac)

KiteGeometry.model_columns(::Type{TetherModel}) =
    (:init_stretched_length, :init_tether_force, :init_stretch_frac)

"""
    WinchModel

The columns of a KiteGeometry `winches` row only SymbolicAWEModels reads, each `nothing`
where the row leaves it unset.

$(TYPEDFIELDS)
"""
struct WinchModel <: KiteGeometry.AbstractModel
    "Reel-out speed at the start [m/s]."
    init_vel::Union{Nothing, Float64}
    "Brake setting."
    brake::Union{Nothing, Float64}
    "Whether the winch follows a set speed rather than a set torque."
    speed_controlled::Union{Nothing, Bool}
    "Friction smoothing width [m/s]."
    friction_epsilon::Union{Nothing, Float64}
end

"""
    TimoshenkoTubeModel

The columns of a KiteGeometry `tubes` row with model `timoshenko`, a
[`TimoshenkoTube`](@ref), each `nothing` where the row leaves it unset.

$(TYPEDFIELDS)
"""
struct TimoshenkoTubeModel <: KiteGeometry.AbstractModel
    "Axial rigidity [N]."
    EA::Union{Nothing, Float64}
    "Shear rigidity [N]."
    GA::Union{Nothing, Float64}
    "Torsional rigidity [N·m²]."
    GJ::Union{Nothing, Float64}
    "Bending rigidity about the y-axis [N·m²]."
    EIy::Union{Nothing, Float64}
    "Bending rigidity about the z-axis [N·m²]."
    EIz::Union{Nothing, Float64}
    "Shear correction coefficient."
    shear_coeff::Union{Nothing, Float64}
    "Rayleigh damping, damping over stiffness [s]."
    damping::Union{Nothing, Float64}
    "Tube end in the first body's frame [m]."
    anchor_a::Union{Nothing, Vector{Float64}}
    "Tube end in the second body's frame [m]."
    anchor_b::Union{Nothing, Vector{Float64}}
end

"""
    ElasticTubeModel

The columns of a KiteGeometry `tubes` row with model `elastic`, an
[`ElasticTube`](@ref), each `nothing` where the row leaves it unset.

$(TYPEDFIELDS)
"""
struct ElasticTubeModel <: KiteGeometry.AbstractModel
    "Axial stiffness [N/m]."
    stiffness_axial::Union{Nothing, Float64}
    "Shear stiffness [N/m]."
    stiffness_shear::Union{Nothing, Float64}
    "Torsional stiffness [N·m/rad]."
    stiffness_torsion::Union{Nothing, Float64}
    "Bending stiffness [N·m/rad]."
    stiffness_bending::Union{Nothing, Float64}
    "Rayleigh damping, damping over stiffness [s]."
    damping::Union{Nothing, Float64}
    "Tube end in the first body's frame [m]."
    anchor_a::Union{Nothing, Vector{Float64}}
    "Tube end in the second body's frame [m]."
    anchor_b::Union{Nothing, Vector{Float64}}
end

"""
    register_definition_models()

Register SymbolicAWEModels' models of each KiteGeometry block, so `load_authoring` and
`load_structure` read its columns into them: the default model of points, segments, stations,
wings, pulleys, bodies, tethers and winches, and the tube models `timoshenko`, the default, and
`elastic`.
"""
function register_definition_models()
    register_model!(:points, "sam_point", PointModel,
                    ("-", "N*s/m", "N*s/m", "m^2", "-", "-", "-", "-", "-", "m");
                    default=true)
    register_model!(:segments, "sam_segment", SegmentModel, ("N*s", "-", "-");
                    default=true)
    register_model!(:spanwise_stations, "sam_station", StationModel,
                    ("-", "-", "N*m*s/rad", "N*m/rad", "-", "-", "-", "-");
                    default=true)
    register_model!(:wings, "sam_wing", WingModel,
                    ("-", "1/s", "-", "m", "-", "-", "-", "-", "-", "m", "m", "m^2", "-",
                     "-");
                    default=true)
    register_model!(:pulleys, "sam_pulley", PulleyModel,
                    ("-", "-", "N*s/m", "-", "m/s"); default=true)
    register_model!(:bodies, "sam_body", BodyModel,
                    ("-", "-", "1/s", "1/s", "1/s", "-", "-", "N", "N", "N*m", "-");
                    default=true)
    register_model!(:tethers, "sam_tether", TetherModel, ("m", "N", "-"); default=true)
    register_model!(:winches, "sam_winch", WinchModel, ("m/s", "-", "-", "m/s");
                    default=true)
    register_model!(:tubes, "timoshenko", TimoshenkoTubeModel,
                    ("N", "N", "N*m^2", "N*m^2", "N*m^2", "-", "s", "m", "m");
                    default=true)
    register_model!(:tubes, "elastic", ElasticTubeModel,
                    ("N/m", "N/m", "N*m/rad", "N*m/rad", "s", "m", "m"))
    return nothing
end
