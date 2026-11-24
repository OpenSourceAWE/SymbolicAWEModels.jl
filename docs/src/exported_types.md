```@meta
CurrentModule = SymbolicAWEModels
```

## Introduction

The `SystemDefinition` provides a flexible framework for defining mechanical
systems using discrete mass-spring-damper models. It serves as input to the
[`SymbolicAWEModel`](@ref), which automatically generates symbolic differential
algebraic equations from the structural definition.

See [Building a system using Julia](tutorial_julia.md) and
[Building a system using YAML](tutorial_yaml.md) for tutorials on creating systems.

## Public enumerations

```@docs
DynamicsType
WingType
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
UnsteadyAero
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
SymbolicAWEModel(set::Settings, sys_struct::SystemDefinition; kwargs...)
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

A system is built from KiteGeometry's keyword rows (`Point`, `Segment`, `SpanwiseStation`,
`Pulley`, `Tether`, `Winch`, `Wing`, `Body`, `Tube`, `Material`), each with one of the
models below in its `model` column. Every model field defaults to `nothing`, the run's
default. [`WingModel`](@ref) and [`BodyModel`](@ref) are documented in
[Private API — definition](private_definition.md).

```@docs
system_definition(::AbstractString, ::Any)
PointModel
StationModel
SegmentSpring
PulleyModel
WinchModel
VSMEngine
create_plate_interpolations
TimoshenkoTubeModel
ElasticTubeModel
AbstractTubeModel
TimoshenkoTube
TimoshenkoTube(; EA, GA, GJ, EIy, EIz, shear_coeff, damping)
ElasticTube
ElasticTube(; stiffness_axial, stiffness_shear, stiffness_torsion, stiffness_bending, damping)
```

## Inflated-tube rigidity

```@docs
TubeRigidityLaw
TUBE_SHEAR_COEFF
TUBE_POISSON_RATIO
```

## System state

```@docs
SystemState
SysState
update_sys_state!
update_from_sysstate!
```
