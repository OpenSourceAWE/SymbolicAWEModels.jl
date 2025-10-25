```@meta
CurrentModule = SymbolicAWEModels
```

## Introduction

The [`SystemStructure`](https://www.google.com/search?q=%40ref) provides a flexible framework for defining the physical
structure of airborne wind energy (AWE) systems using discrete mass-spring-damper models.
This structure can represent many different AWE system configurations, from simple
single-line kites to complex multi-wing systems with intricate bridle networks.

The [`SystemStructure`](https://www.google.com/search?q=%40ref) serves as input to the [`SymbolicAWEModel`](https://www.google.com/search?q=%40ref), which is
based on ModelingToolkit and automatically generates symbolic differential algebraic
equations from the structural definition.

## Workflow

1.  Define system components ([`Point`](https://www.google.com/search?q=%40ref), [`Segment`](https://www.google.com/search?q=%40ref), [`Group`](https://www.google.com/search?q=%40ref), etc.)
2.  Assemble into a [`SystemStructure`](https://www.google.com/search?q=%40ref)
3.  Pass to [`SymbolicAWEModel`](https://www.google.com/search?q=%40ref) for automatic MTK model generation
4.  Simulate the resulting symbolic model

## Public enumerations

```@docs
SegmentType
DynamicsType
```

## Core Model Type

This is the main struct that defines any complete simulation model.

```@docs
SymbolicAWEModel
SymbolicAWEModel(set::Settings, sys_struct::SystemStructure; kwargs...)
SymbolicAWEModel(set::Settings; kwargs...)
SymbolicAWEModel(set::Settings, name::String; kwargs...)
```

## System structure and components

```@docs
SystemStructure
SystemStructure(name, set; points, groups, segments, pulleys, tethers, winches, wings, transforms)
Point
Point(idx, pos_cad, type; wing_idx, vel_w, transform_idx, mass, bridle_damping, fix_sphere)
Group
Group(idx, point_idxs, vsm_wing::RamAirWing, gamma, type, moment_frac)
Group(idx, point_idxs, le_pos, chord, y_airf, type, moment_frac)
Segment
Segment(idx, set, point_idxs, type; l0, compression_frac, axial_stiffness, axial_damping)
Segment(idx, point_idxs, axial_stiffness, axial_damping, diameter; l0, compression_frac)
Pulley
Pulley(idx, segment_idxs, type)
Tether
Tether(idx, segment_idxs, winch_idx)
Winch
Winch(idx, set::Settings, tether_idxs; tether_len, tether_vel, brake)
Winch(idx, tether_idxs, gear_ratio, drum_radius, f_coulomb, c_vf, inertia_total; tether_len, tether_vel, brake)
AbstractWing
BaseWing
VSMWing
Wing
Wing(idx, vsm_aero, vsm_wing, vsm_solver, group_idxs, R_b_c, pos_cad; transform_idx)
Transform
Transform(idx, elevation, azimuth, heading; base_point_idx, base_pos, base_transform_idx, wing_idx, rot_point_idx)
Transform(idx, set, base_point_idx; kwargs...)
```

## System state

```@docs
SysState
update_sys_state!
update_from_sysstate!
```

