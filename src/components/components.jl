# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

"""
# Modern ModelingToolkit Components for AWE Systems

This module provides a component-based architecture for building Airborne Wind Energy
(AWE) system models using ModelingToolkit.jl's modern `@mtkmodel` and `@connector` approach.

## Overview

The component library enables modular, reusable, and composable modeling of AWE systems.
Instead of writing monolithic equation-generating functions, users can:

1. **Instantiate components** (points, segments, winches, wings)
2. **Connect them** using symbolic connections
3. **Build systems** hierarchically
4. **Leverage MTK optimizations** (structural_simplify, code generation)

## Component Hierarchy

### Connectors
- `MechanicalNode`: 3D mechanical connection point (position, velocity, force)
- `ControlSignal`: Scalar control input connector

### Physical Components
- `PointMass`: Dynamic, static, or quasi-static point mass
- `TetherSegment`: Elastic tether segment with aerodynamic drag
- `PulleyComponent`: Length redistribution constraint between segments
- `WinchComponent`: Ground-based motor/drum/gearbox system
- `RigidWing`: 6-DOF wing with quaternion dynamics (STUB - to be implemented)
- `WingGroup`: Wing section with twist dynamics (STUB - to be implemented)

**Note**: Component names have `Component` suffix or different names (e.g., `PointMass` instead of `Point`)
to avoid conflicts with legacy `SystemStructure` types.

## Usage Example

### Simple Tether System
```julia
using ModelingToolkit, SymbolicAWEModels
using ModelingToolkit: t_nounits as t, D_nounits as D

# Create components
@named ground = MassPoint(dynamics_type=4, fixed_pos=[0, 0, 0], mass=1.0)
@named kite = MassPoint(dynamics_type=1, mass=5.0)
@named tether = SpringDamperSegment(
    l0 = 100.0,
    axial_stiffness = 1.2e6,
    axial_damping = 500.0,
    diameter = 0.004
)

# Connect components
eqs = [
    connect(ground.node, tether.node1)
    connect(tether.node2, kite.node)
]

# Build and simplify system
@named sys = ODESystem(eqs, t)
simple_sys = structural_simplify(sys)

# Create problem and simulate
prob = ODEProblem(simple_sys, [], (0.0, 10.0))
sol = solve(prob, Tsit5())
```

### Winch-Controlled Tether
```julia
@named winch = WinchComponent(gear_ratio=110.0, drum_radius=0.1632)
@named ground_point = PointMass(dynamics_type=4, fixed_pos=[0,0,0])
@named tether1 = TetherSegment(l0=50.0, ...)

# Control winch torque
eqs = [
    connect(ground_point.node, tether1.node1)
    winch.ctrl.value ~ 50.0  # 50 Nm motor torque
    winch.tether_force ~ tether1.spring_force  # Link tension
    winch.tether_len ~ tether1.l0  # Control segment length
]
```

## Design Principles

### 1. Modularity
Each component encapsulates its own physics and is self-contained. Users can:
- Mix and match components
- Create custom variants (e.g., different wing models)
- Compose complex systems from simple parts

### 2. Hybrid Connection Strategy
- **Physical connectors** (`@connector`): For mechanical interfaces where Kirchhoff's
  laws apply (force balance at connection points)
- **Direct symbolic linking**: For control signals, aerodynamic coupling, and
  derived quantities

### 3. Backward Compatibility
The legacy `SystemStructure` approach in `src/system_structure.jl` and
`src/generate_system.jl` remains fully functional. Users can choose:
- Legacy: `SymbolicAWEModel(set, "ram")`
- Components: `SymbolicAWEModel(set, "ram"; use_components=true)` (future)

## Component Implementation Status

✅ **Complete:**
- `MechanicalNode` connector
- `ControlSignal` connector
- `MassPoint` with all dynamics types
- `SpringDamperSegment` with drag
- `Pulley` with dynamic/quasi-static modes
- `Winch` with full motor/friction dynamics

🚧 **Stubs (to be implemented):**
- `RigidWing`: Quaternion-based 6-DOF wing
- `DeformableGroup`: Wing twist dynamics

## Architecture Benefits

### Performance
- `structural_simplify` eliminates redundant equations
- Better code generation from ModelingToolkit
- Potential for automatic parallelization

### Maintainability
- Physics equations live in component definitions (declarative)
- Easy to test components independently
- Clear separation of concerns

### Extensibility
- Add new components without modifying existing code
- Users can define custom components in their own code
- Component library can grow independently

## Migration Path

For users of the legacy system:

1. **Current code still works** - No breaking changes
2. **Experiment with components** - Try building simple systems
3. **Gradual adoption** - Mix legacy and component approaches
4. **Full migration** - When ready, switch to component-based builder

## Developer Notes

### Adding a New Component

1. Create file in `src/components/<name>.jl`
2. Define `@mtkmodel <Name> begin ... end`
3. Add comprehensive docstring with equations
4. Include in `src/components/components.jl`
5. Add tests in `test/test_components.jl`

### Component Checklist
- [ ] Clear parameter descriptions
- [ ] Documented state variables
- [ ] Mathematical equations in docstring
- [ ] Usage example
- [ ] Unit tests
- [ ] Integration test in system context

## References

- [ModelingToolkit.jl Documentation](https://docs.sciml.ai/ModelingToolkit/stable/)
- [MTK Components Tutorial](https://docs.sciml.ai/ModelingToolkit/stable/tutorials/acausal_components/)
- [Original generate_system.jl](../generate_system.jl) - Legacy approach
- [SystemStructure](../system_structure.jl) - Legacy data structures
"""

# Include all component definitions
include("mass_point.jl")
include("segment.jl")
include("pulley.jl")
include("winch.jl")
include("rigid_wing.jl")
include("wing_group.jl")
