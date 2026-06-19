```@meta
CurrentModule = SymbolicAWEModels
```

## Introduction

The [`SystemStructure`](@ref) provides a flexible framework for defining mechanical
systems using discrete mass-spring-damper models. It serves as input to the
[`SymbolicAWEModel`](@ref), which automatically generates symbolic differential
algebraic equations from the structural definition.

See [Building a system using Julia](tutorial_julia.md) and
[Building a system using YAML](tutorial_yaml.md) for tutorials on creating systems.

## Public enumerations

```@docs
DynamicsType
WingType
```

## Aerodynamic models

```@docs
AbstractAeroModel
AeroNone
AbstractVSMAero
AeroDirect
AeroLinearized
AeroPlate
ContinuousAero
aero_component
is_builtin_aero
aero_hash_id
```

## Winch models

```@docs
AbstractWinchModel
DefaultWinchModel
```

## Core model type

```@docs
SymbolicAWEModel
SymbolicAWEModel(set::Settings, sys_struct::SystemStructure; kwargs...)
```

## System structure and components

```@docs
SystemStructure
SystemStructure(name, set; points, groups, segments, pulleys, tethers, winches, wings, transforms)
Point
Point(name, pos_cad, type; wing, transform, extra_mass, body_frame_damping, world_frame_damping, fix_sphere)
TwistSurface
TwistSurface(name, points, type, moment_frac; damping=50.0)
Segment
Segment(name, set, point_i, point_j; l0, compression_frac, diameter_mm, unit_stiffness, unit_damping)
Segment(name, point_i, point_j, unit_stiffness, unit_damping, diameter; l0, compression_frac)
Pulley
Pulley(name, segment_i, segment_j, type)
Tether
Tether(name, segments::AbstractVector, stretched_length; start_point, end_point, tether_force, stretch_frac)
Tether(name, stretched_length; start_point, end_point, n_segments, unit_stiffness, unit_damping, diameter, tether_force, stretch_frac)
Winch
Winch(name, set::Settings, tethers; winch_point, init_vel, brake)
Winch(name, tethers, gear_ratio, drum_radius, f_coulomb, c_vf, inertia_total; winch_point, init_vel, brake)
AbstractWing
Wing
VSMEngine
VSMWing
PlateWing
create_plate_interpolations
RigidBody
RigidBody(name; mass, inertia_principal, pos, vel, Q_b_to_w, ω_b, com_offset_b, R_b_to_p, angular_damping, ext_force_w, ext_moment_b)
ElasticJoint
ElasticJoint(name, body_a, body_b; anchor_a, anchor_b, stiffness_axial, stiffness_shear, stiffness_torsion, stiffness_bending, damping_trans, damping_rot)
Transform
Transform(name, elevation, azimuth, heading; base_point, base_pos, base_transform, wing, rot_point)
```

## Indexing

```@docs
NamedCollection
NameRef
```

## System state

```@docs
SysState
update_sys_state!
update_from_sysstate!
```
