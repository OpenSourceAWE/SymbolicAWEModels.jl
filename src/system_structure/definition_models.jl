# Copyright (c) 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

"""
    PointModel

The columns of a KiteGeometry `points` row only SymbolicAWEModels reads, each `nothing`
where the row leaves it unset.

$(TYPEDFIELDS)
"""
struct PointModel <: KiteGeometry.AbstractModel
    "Dynamics type, such as `DYNAMIC`, `STATIC` or `BODY_STATIC`."
    type::Union{Nothing, String}
    "Damping in the body frame [N·s/m], one value or one per axis."
    body_frame_damping::Union{Nothing, Float64, Vector{Float64}}
    "Damping in the world frame [N·s/m], one value or one per axis."
    world_frame_damping::Union{Nothing, Float64, Vector{Float64}}
    "Drag area [m²]."
    area::Union{Nothing, Float64}
    "Drag coefficient."
    drag_coeff::Union{Nothing, Float64}
    "The tube the point rides on, by name or index."
    tube::Union{Nothing, KiteGeometry.NameRef}
    "Whether the point is held on a sphere."
    fix_sphere::Union{Nothing, Bool}
    "Whether the point is frozen in place."
    fix_static::Union{Nothing, Bool}
    "The wing the point belongs to, by name or index."
    wing_idx::Union{Nothing, KiteGeometry.NameRef}
    "Where the point sits in its body's frame [m]."
    anchor_b::Union{Nothing, Vector{Float64}}
end

"""
    StationModel

The columns of a KiteGeometry `spanwise_stations` row only SymbolicAWEModels reads, each
`nothing` where the row leaves it unset.

$(TYPEDFIELDS)
"""
struct StationModel <: KiteGeometry.AbstractModel
    "Dynamics type of the station's twist, such as `DYNAMIC` or `STATIC`."
    type::Union{Nothing, String}
    "Chordwise twist axis as a fraction of the chord (0 leading edge, 1 trailing edge)."
    moment_frac::Union{Nothing, Float64}
    "Twist damping [N·m·s/rad]."
    damping::Union{Nothing, Float64}
    "Twist restoring stiffness [N·m/rad]."
    stiffness::Union{Nothing, Float64}
    "The bodies the station twists with, by name or index."
    bodies::Union{Nothing, Vector{KiteGeometry.NameRef}}
    "The main and the flap body of a body flap, by name or index."
    flap_bodies::Union{Nothing, Vector{KiteGeometry.NameRef}}
    "The three points of a point flap, by name or index."
    flap_points::Union{Nothing, Vector{KiteGeometry.NameRef}}
    "Unit flap-hinge axis."
    flap_axis::Union{Nothing, Vector{Float64}}
end

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
end

"""
    PulleyModel

The columns of a KiteGeometry `pulleys` row only SymbolicAWEModels reads, each `nothing`
where the row leaves it unset.

$(TYPEDFIELDS)
"""
struct PulleyModel <: KiteGeometry.AbstractModel
    "Dynamics type, such as `DYNAMIC` or `STATIC`."
    type::Union{Nothing, String}
    "Fraction of line tension the sheave passes on (0-1)."
    efficiency::Union{Nothing, Float64}
    "Damping on rope travel [N·s/m]."
    damping::Union{Nothing, Float64}
    "Whether the rope split is frozen."
    brake::Union{Nothing, Bool}
    "Friction smoothing width [m/s]."
    friction_epsilon::Union{Nothing, Float64}
end

"""
    BodyModel

The columns of a KiteGeometry `bodies` row only SymbolicAWEModels reads, each `nothing`
where the row leaves it unset.

$(TYPEDFIELDS)
"""
struct BodyModel <: KiteGeometry.AbstractModel
    "Dynamics type, `DYNAMIC` or `STATIC`."
    type::Union{Nothing, String}
    "The wing whose body frame the body's damping acts in, by name or index."
    wing::Union{Nothing, KiteGeometry.NameRef}
    "Spin damping about the body axes [1/s], one value or one per axis."
    angular_damping::Union{Nothing, Float64, Vector{Float64}}
    "Damping per unit mass on the world axes [1/s], one value or one per axis."
    world_frame_damping::Union{Nothing, Float64, Vector{Float64}}
    "Damping per unit mass on the body axes [1/s], one value or one per axis."
    body_frame_damping::Union{Nothing, Float64, Vector{Float64}}
    "Whether the body is held on a sphere."
    fix_sphere::Union{Nothing, Bool}
    "Whether the body is frozen in place."
    fix_static::Union{Nothing, Bool}
    "External force, world frame [N]."
    ext_force_w::Union{Nothing, Vector{Float64}}
    "External force, body frame [N]."
    ext_force_b::Union{Nothing, Vector{Float64}}
    "External moment, body frame [N·m]."
    ext_moment_b::Union{Nothing, Vector{Float64}}
    "How the principal frame is found: `EIGEN_DECOMP` or `Y_ROTATION`."
    principal_frame_method::Union{Nothing, String}
end

"""
    TetherModel

The columns of a KiteGeometry `tethers` row only SymbolicAWEModels reads, each `nothing`
where the row leaves it unset: how the tether is placed.

$(TYPEDFIELDS)
"""
struct TetherModel <: KiteGeometry.AbstractModel
    "Length the tether is stretched to when placed [m]."
    init_stretched_length::Union{Nothing, Float64}
    "Tension the placed tether carries, which sets its rest length [N]."
    init_tether_force::Union{Nothing, Float64}
    "Strain of the placed tether, which sets its rest length."
    init_stretch_frac::Union{Nothing, Float64}
end

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
`load_structure` read its columns into them: the default model of points, stations, wings,
pulleys, bodies, tethers and winches, and the tube models `timoshenko`, the default, and
`elastic`.
"""
function register_definition_models()
    register_model!(:points, "sam_point", PointModel,
                    ("-", "N*s/m", "N*s/m", "m^2", "-", "-", "-", "-", "-", "m");
                    default=true)
    register_model!(:spanwise_stations, "sam_station", StationModel,
                    ("-", "-", "N*m*s/rad", "N*m/rad", "-", "-", "-", "-");
                    default=true)
    register_model!(:wings, "sam_wing", WingModel,
                    ("-", "1/s", "-", "m", "-", "-", "-", "-", "-", "m", "m", "m^2", "-");
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
