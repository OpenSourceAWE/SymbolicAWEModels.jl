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
SegmentType
PrincipalFrameMethod
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
AeroPressure
aero_component
is_builtin_aero
aero_hash_id
```

## Winch models

```@docs
AbstractWinchModel
TorqueWinch
CascadedLengthWinch
```

## Core model type

```@docs
SymbolicAWEModel
SymbolicAWEModel(set::Settings, sys_struct::SystemStructure; kwargs...)
```

## Backends

```@docs
ModelBackend
MonolithBackend
KernelBackend
BackendUnsupportedError
default_backend
default_backend!
```

## System structure and components

```@docs
SystemStructure
SystemStructure(name, set; points, twist_surfaces, segments, pulleys, tethers, winches, wings, transforms, bodies, elastic_joints, timoshenko_joints)
Point
Point(name, pos_cad, type; wing, transform, extra_mass, body_frame_damping, world_frame_damping, fix_sphere)
TwistSurface
TwistSurface(name, points, type, moment_frac; damping=50.0)
Segment
Segment(name, set, point_i, point_j; l0, compression_frac, compression_damping_frac, diameter_mm, unit_stiffness, unit_damping, density, youngs_modulus, damping_per_stiffness)
Segment(name, point_i, point_j, unit_stiffness, unit_damping, diameter; l0, compression_frac, compression_damping_frac)
Pulley
Pulley(name, segment_i, segment_j, type; efficiency, damping, brake, friction_epsilon)
Tether
Tether(name, segments::AbstractVector, stretched_length; start_point, end_point, tether_force, stretch_frac)
Tether(name, stretched_length; start_point, end_point, n_segments, unit_stiffness, unit_damping, diameter, tether_force, stretch_frac)
Winch
Winch(name, set::Settings, tethers; winch_point, init_vel, brake)
Winch(name, tethers, gear_ratio, drum_radius, coulomb_friction, viscous_coefficient, inertia_total; winch_point, init_vel, brake)
AbstractWing
RigidWing
ParticleWing
Wing
VSMEngine
VSMWing
PlateWing
create_plate_interpolations
Body
Body(name; mass, inertia_principal, pos, vel, Q_b_to_w, ω_b, com_offset_b, R_b_to_p, angular_damping, ext_force_w, ext_moment_b)
ElasticJoint
ElasticJoint(name, body_a, body_b; anchor_a, anchor_b, stiffness_axial, stiffness_shear, stiffness_torsion, stiffness_bending, damping_trans, damping_rot)
TimoshenkoJoint
TimoshenkoJoint(name, body_a, body_b; anchor_a, anchor_b, EA, GA, GJ, EIy, EIz, shear_coeff, damping_trans, damping_rot, rest_length)
Transform
Transform(name, elevation, azimuth, heading; base_point, base_pos, base_transform, wing, rot_point)
```

## Indexing

```@docs
NamedCollection
NameRef
WeightedRefPoints
```

## Inflated-tube rigidity

```@docs
TubeRigidityLaw
TUBE_SHEAR_COEFF
TUBE_POISSON_RATIO
```

## System state

```@docs
SysState
update_sys_state!
update_from_sysstate!
```
