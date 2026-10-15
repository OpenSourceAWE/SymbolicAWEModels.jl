# Copyright (c) 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

"""
    PointModel

The columns of a KiteGeometry `points` row only SymbolicAWEModels reads, each `nothing`
where the row leaves it unset.

$(TYPEDFIELDS)
"""
struct PointModel <: KiteGeometry.AbstractModel
    "Damping in the body frame [N·s/m], one value or one per axis."
    body_frame_damping::Union{Nothing, Float64, Vector{Float64}}
    "Damping in the world frame [N·s/m], one value or one per axis."
    world_frame_damping::Union{Nothing, Float64, Vector{Float64}}
    "The tube the point rides on, by name or index."
    tube::Union{Nothing, KiteGeometry.NameRef}
    "Whether the point is held on a sphere."
    fix_sphere::Union{Nothing, Bool}
    "Whether the point is frozen in place."
    fix_static::Union{Nothing, Bool}
end

"""
    StationModel

The columns of a KiteGeometry `stations` row only SymbolicAWEModels reads, each `nothing`
where the row leaves it unset.

$(TYPEDFIELDS)
"""
struct StationModel <: KiteGeometry.AbstractModel
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
end

"""
    PulleyModel

The columns of a KiteGeometry `pulleys` row only SymbolicAWEModels reads, each `nothing`
where the row leaves it unset.

$(TYPEDFIELDS)
"""
struct PulleyModel <: KiteGeometry.AbstractModel
    "Damping on rope travel [N·s/m]."
    damping::Union{Nothing, Float64}
    "Whether the rope split is frozen."
    brake::Union{Nothing, Bool}
    "Friction smoothing width [m/s]."
    friction_epsilon::Union{Nothing, Float64}
end

"""
    register_definition_models()

Register [`PointModel`](@ref), [`StationModel`](@ref), [`WingModel`](@ref) and
[`PulleyModel`](@ref) as the default models of their KiteGeometry blocks, so
`load_authoring` and `load_structure` read SymbolicAWEModels' columns into them.
"""
function register_definition_models()
    register_model!(:points, "sam_point", PointModel,
                    ("N*s/m", "N*s/m", "-", "-", "-"); default=true)
    register_model!(:stations, "sam_station", StationModel,
                    ("-", "N*m*s/rad", "N*m/rad", "-", "-", "-", "-"); default=true)
    register_model!(:wings, "sam_wing", WingModel, ("-", "1/s", "-", "m", "-");
                    default=true)
    register_model!(:pulleys, "sam_pulley", PulleyModel, ("N*s/m", "-", "m/s");
                    default=true)
    return nothing
end
