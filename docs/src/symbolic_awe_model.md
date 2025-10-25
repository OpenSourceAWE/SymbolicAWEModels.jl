```@meta
CurrentModule = SymbolicAWEModels
```

# SymbolicAWEModel

The [`SymbolicAWEModel`](@ref) is the main type that encapsulates a complete symbolic simulation model of an Airborne Wind Energy (AWE) system. It uses ModelingToolkit.jl to automatically generate symbolic differential algebraic equations from a `SystemStructure` definition.

## Overview

A [`SymbolicAWEModel`](@ref) contains:
- **SystemStructure**: Physical topology of the system (points, segments, tethers, wings, winches)
- **ODEProblem**: Compiled ModelingToolkit differential equations
- **Integrator**: OrdinaryDiffEq integrator for time-stepping
- **Linear models**: Simplified state-space representations for control
- **Settings**: Configuration from YAML files

## Creating a SymbolicAWEModel

There are several ways to create a [`SymbolicAWEModel`](@ref):

### 1. From Settings (using default physical model)

```julia
using SymbolicAWEModels
set = Settings("system.yaml")
sam = SymbolicAWEModel(set)
```

This creates a model using the `physical_model` specified in the settings file (e.g., `"ram"`, `"simple_ram"`, `"4_attach_ram"`).

### 2. From Settings with explicit model name

```julia
sam = SymbolicAWEModel(set, "ram")  # Detailed ram-air kite
sam = SymbolicAWEModel(set, "simple_ram")  # Simplified model
sam = SymbolicAWEModel(set, "4_attach_ram")  # Most detailed model
```

### 3. From custom SystemStructure

```julia
# Create custom system structure (see System Structure page)
sys_struct = SystemStructure("my_model", set; points, segments, tethers, winches, wings)
sam = SymbolicAWEModel(set, sys_struct)
```

## Available Physical Models

The package provides several predefined physical models via factory functions:

- **`"ram"`** (standard): Ram-air kite with simplified bridle system
- **`"simple_ram"`**: Single-segment tethers, no bridle system
- **`"4_attach_ram"`** (detailed): Full ram-air kite with 4-point bridle attachment, deformable wing groups, pulleys, and 3 winches
- **`"tether"`**: Helper model for tether property calculations

See the [System Structure](system_structure.md#Predefined-System-Structures) page for the factory functions that create these models.

## Basic Workflow

```julia
# 1. Create model
set = Settings("system.yaml")
sam = SymbolicAWEModel(set, "ram")

# 2. Initialize
init!(sam)

# 3. Find steady state
find_steady_state!(sam)

# 4. Simulate
(log, _) = sim_oscillate!(sam; dt=0.05, total_time=10.0)

# 5. Visualize (requires GLMakie)
using GLMakie
plot(sam.sys_struct, log)
```

## Type Documentation

```@docs
SymbolicAWEModel
SymbolicAWEModel(set::Settings, sys_struct::SystemStructure; kwargs...)
SymbolicAWEModel(set::Settings; kwargs...)
SymbolicAWEModel(set::Settings, name::String; kwargs...)
```

## Wing Types

The package provides specialized wing types for aerodynamic modeling. Base types (`AbstractWing`, `BaseWing`) are defined in [KiteUtils.jl](https://opensourceawe.github.io/KiteUtils.jl/dev/).

```@docs
VSMWing
Wing
Wing(idx, vsm_aero, vsm_wing, vsm_solver, group_idxs, R_b_c, pos_cad; transform_idx)
```

## See Also

- [System Structure](system_structure.md) - Physical topology definition
- [State Management](state_management.md) - Working with system state
- [Simulation Functions](exported_functions.md) - Running simulations
- [Custom Model Tutorial](tutorial_system_structure.md) - Step-by-step examples
