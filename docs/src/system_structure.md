```@meta
CurrentModule = SymbolicAWEModels
```

# System Structure

The `SystemStructure` type defines the physical topology of Airborne Wind Energy (AWE) systems using discrete mass-spring-damper models. It describes the structure using points, segments, tethers, wings, winches, and other components.

## Important Note

The `SystemStructure` type and its components (`Point`, `Segment`, `Tether`, `Winch`, `Pulley`, `Group`, `Transform`, etc.) are **defined in the KiteUtils.jl package**, not in SymbolicAWEModels.jl.

**For comprehensive documentation of SystemStructure and all its components, please refer to:**

**[KiteUtils.jl System Structure Documentation](https://opensourceawe.github.io/KiteUtils.jl/stable/system_structure/)**

This page documents only the SymbolicAWEModels-specific extensions and usage patterns.

## Overview

The [`SystemStructure`](@ref) serves as input to the [`SymbolicAWEModel`](@ref), which uses ModelingToolkit.jl to automatically generate symbolic differential algebraic equations from the structural definition.

A `SystemStructure` can represent many different AWE system configurations:
- Simple single-line kites
- Complex multi-wing systems
- Intricate bridle networks with pulleys
- Systems with multiple winches

## Workflow

1. Define system components ([`Point`](@ref), [`Segment`](@ref), [`Group`](@ref), etc.)
2. Assemble into a [`SystemStructure`](@ref)
3. Pass to [`SymbolicAWEModel`](@ref) for automatic symbolic model generation
4. Simulate using standard simulation functions

## Predefined System Structures

SymbolicAWEModels provides factory functions to create common system topologies:

```@docs
create_ram_sys_struct
create_simple_ram_sys_struct
create_tether_sys_struct
```

## Example: Creating a Custom Structure

```julia
using SymbolicAWEModels

set = Settings("system.yaml")

# Define components
points = Point[]
segments = Segment[]
tethers = Tether[]
winches = Winch[]

# Add static ground point
push!(points, Point(1, zeros(3), STATIC))

# Add dynamic points and segments
for i in 1:10
    pos = [0.0, 0.0, i * 5.0]
    push!(points, Point(i+1, pos, DYNAMIC))
    push!(segments, Segment(i, set, (i, i+1), BRIDLE))
end

# Create structure
sys_struct = SystemStructure("custom_model", set;
                             points=points,
                             segments=segments)

# Use in simulation
sam = SymbolicAWEModel(set, sys_struct)
init!(sam)
```

For more detailed examples, see the [Custom Model Tutorial](@ref).

## Component Types Defined in KiteUtils

The following types are documented in [KiteUtils.jl](https://opensourceawe.github.io/KiteUtils.jl/stable/system_structure/):

- **`SystemStructure`** - Main container type
- **`Point`** - Mass particles with position/velocity
- **`Segment`** - Spring-damper connections between points
- **`Tether`** - Collection of segments forming a line
- **`Pulley`** - Constraint enforcing equal tension
- **`Group`** - Deformable wing section with twist dynamics
- **`Winch`** - Motor/generator for tether control
- **`Transform`** - Orientation and positioning transformations

### Enumerations

- **`DynamicsType`** - Point dynamics: `STATIC`, `QUASI_STATIC`, `DYNAMIC`, `WING`
- **`SegmentType`** - Segment types: `BRIDLE`, `TETHER`, etc.

## SymbolicAWEModels-Specific Types

The following wing-related types are defined in SymbolicAWEModels:

```@docs
VSMWing
Wing
```

## See Also

- **[KiteUtils System Structure Docs](https://opensourceawe.github.io/KiteUtils.jl/stable/system_structure/)** - Complete component documentation
- [SymbolicAWEModel](@ref) - Using structures in simulations
- [Custom Model Tutorial](@ref) - Step-by-step examples
