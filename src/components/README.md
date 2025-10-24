# Modern ModelingToolkit Components for AWE Systems

This directory contains a **component-based architecture** for building Airborne Wind Energy (AWE) system models using ModelingToolkit.jl's modern `@mtkmodel` and `@connector` approach.

## 🎯 Goals

1. **Modularity**: Build complex AWE systems from composable, reusable components
2. **Performance**: Leverage MTK's `structural_simplify` for automatic equation optimization
3. **Extensibility**: Add new physics models without modifying existing code
4. **Backward Compatibility**: Coexists with legacy `SystemStructure` approach

## 📂 Directory Structure

```
src/components/
├── README.md              # This file
├── components.jl          # Main include file with documentation
├── connectors.jl          # @connector definitions
├── mass_point.jl          # @mtkmodel MassPoint
├── segment.jl             # @mtkmodel SpringDamperSegment
├── pulley.jl              # @mtkmodel Pulley
├── winch.jl               # @mtkmodel Winch
├── wing.jl                # @mtkmodel RigidWing (stub)
├── group.jl               # @mtkmodel DeformableGroup (stub)
└── examples.jl            # Example system builders
```

## 🧩 Component Library

### Connectors

#### `MechanicalNode`
3D mechanical connection point with position, velocity, and force.
- **Variables**: `pos(t)[1:3]`, `vel(t)[1:3]`, `force(t)[1:3]`
- **Flow variable**: `force` (Kirchhoff's law enforced at connections)

#### `ControlSignal`
Scalar control input (e.g., winch torque, steering angle).
- **Variables**: `value(t)`

### Physical Components

#### `MassPoint`
Point mass with configurable dynamics:
- **DYNAMIC** (type=1): Newton's 2nd law, `F = ma`
- **STATIC** (type=4): Fixed position
- **QUASI_STATIC** (type=2): Force balance, `F = 0`
- **WING** (type=3): Attached to wing rigid body

**Key parameters**: `mass`, `dynamics_type`, `fixed_pos`, `fix_sphere`, `bridle_damping`

#### `SpringDamperSegment`
Elastic tether segment with:
- Axial spring-damper (Hooke's law + damping)
- Asymmetric compression stiffness
- Aerodynamic drag
- Wind shear effects

**Key parameters**: `l0`, `axial_stiffness`, `axial_damping`, `diameter`, `cd_tether`

#### `Pulley`
Length redistribution constraint between two segments:
- Enforces `l1 + l2 = const`
- Dynamic or quasi-static modes
- Accounts for tether mass inertia

**Key parameters**: `sum_len`, `dynamics_type`, `damping`

#### `Winch`
Ground-based motor/drum/gearbox system:
- Motor torque control
- Coulomb + viscous friction
- Rotational inertia
- Brake functionality

**Key parameters**: `gear_ratio`, `drum_radius`, `f_coulomb`, `inertia_total`

#### `RigidWing` (STUB)
6-DOF wing with quaternion orientation (to be implemented).
- Quaternion kinematics
- Euler rotation equations
- VortexStepMethod aerodynamics integration

#### `DeformableGroup` (STUB)
Wing section twist dynamics (to be implemented).
- Twist angle state
- Bridle/aerodynamic moment coupling

## 🚀 Quick Start

### Example 1: Simple Pendulum

```julia
using SymbolicAWEModels, ModelingToolkit, OrdinaryDiffEq

# Create system
sys = simple_pendulum_system()

# Simplify and solve
simple_sys = structural_simplify(sys)
prob = ODEProblem(simple_sys, [], (0.0, 10.0))
sol = solve(prob, Tsit5())

# Plot results
using Plots
plot(sol, vars=[sys.kite.node.pos...])
```

### Example 2: Custom System

```julia
using ModelingToolkit
using ModelingToolkit: t_nounits as t

# Define components
@named ground = MassPoint(dynamics_type=4, fixed_pos=[0,0,0])
@named kite = MassPoint(dynamics_type=1, mass=5.0)
@named tether = SpringDamperSegment(
    l0 = 100.0,
    axial_stiffness = 1.2e6,
    diameter = 0.004
)

# Connect
eqs = [
    connect(ground.node, tether.node1)
    connect(tether.node2, kite.node)
]

# Build
@named sys = ODESystem(eqs, t;
    systems = [ground, kite, tether]
)

simple_sys = structural_simplify(sys)
# ... solve as usual
```

## 🔗 Connection Strategies

### Physical Connectors (Recommended for mechanical interfaces)

```julia
connect(point1.node, segment.node1)
connect(segment.node2, point2.node)
```

**What this does:**
- Enforces position continuity: `point1.pos ~ segment.pos1`
- Enforces velocity continuity: `point1.vel ~ segment.vel1`
- Enforces force balance: `point1.force + segment.force1 ~ 0` (Flow property)

### Direct Symbolic Linking (For control/coupling)

```julia
eqs = [
    winch.ctrl.value ~ 50.0  # Set motor torque
    winch.tether_force ~ tether.spring_force  # Link tension
    pulley.seg1_force ~ seg1.spring_force  # Link forces
]
```

## 📐 Component Design Patterns

### 1. State Variables
Each component declares its differential variables:
```julia
@mtkmodel MyComponent begin
    @variables begin
        pos(t)[1:3]  # State
        vel(t)[1:3]  # State derivative
        acc(t)[1:3]  # Computed from dynamics
    end
end
```

### 2. Parameters
Physical constants and configuration:
```julia
@parameters begin
    mass = 5.0, [description = "Mass [kg]"]
    stiffness = 1e6, [description = "Stiffness [N/m]"]
end
```

### 3. Connectors
Interfaces to other components:
```julia
@components begin
    node = MechanicalNode()
    ctrl = ControlSignal()
end
```

### 4. Equations
Physics captured declaratively:
```julia
@equations begin
    D(pos) ~ vel
    D(vel) ~ acc
    acc ~ force / mass
end
```

### 5. Conditional Logic
Type-dependent behavior:
```julia
@equations begin
    if dynamics_type == 1  # DYNAMIC
        D(pos) ~ vel
        acc ~ force / mass
    elseif dynamics_type == 4  # STATIC
        pos ~ fixed_pos
        vel ~ zeros(3)
    end
end
```

## 🧪 Testing Strategy

### Unit Tests (Component Level)
Each component should be tested in isolation:

```julia
@testset "MassPoint DYNAMIC" begin
    @named pt = MassPoint(dynamics_type=1, mass=5.0)
    # Test equation count, variable names, etc.
end
```

### Integration Tests (System Level)
Test complete systems:

```julia
@testset "Simple pendulum simulation" begin
    sys = simple_pendulum_system()
    simple_sys = structural_simplify(sys)
    prob = ODEProblem(simple_sys, [], (0.0, 1.0))
    sol = solve(prob, Tsit5())
    @test SciMLBase.successful_retcode(sol)
end
```

## 🔧 Development Workflow

### Adding a New Component

1. **Create file**: `src/components/my_component.jl`
2. **Define component**:
   ```julia
   @mtkmodel MyComponent begin
       @parameters begin
           # ... parameters
       end
       @variables begin
           # ... variables
       end
       @components begin
           # ... connectors
       end
       @equations begin
           # ... physics
       end
   end
   ```
3. **Add documentation**: Complete docstring with equations
4. **Include in `components.jl`**: `include("my_component.jl")`
5. **Write tests**: Add to `test/test_components.jl`
6. **Create example**: Add usage example to `examples.jl`

### Component Checklist
- [ ] Clear docstring with mathematical equations
- [ ] All parameters have descriptions and units
- [ ] Variables have descriptions
- [ ] Equations match physics documentation
- [ ] Usage example included
- [ ] Unit test written
- [ ] Integration test in system context

## 📊 Performance Considerations

### Compilation Time
- First run will be slow (~minutes) due to MTK code generation
- Subsequent runs are fast (compiled code is cached)
- Use `PackageCompiler.jl` for production deployments

### Runtime Performance
- `structural_simplify` eliminates redundant equations
- MTK generates optimized Julia code
- Use stiff solvers (e.g., `Rosenbrock23()`, `FBDF()`) for tether dynamics

### Memory Usage
- Components are type-stable and stack-allocated where possible
- Large systems may benefit from sparse matrix representations

## 🎓 Learning Resources

### ModelingToolkit.jl Documentation
- [MTK Components Tutorial](https://docs.sciml.ai/ModelingToolkit/stable/tutorials/acausal_components/)
- [MTK Language Guide](https://docs.sciml.ai/ModelingToolkit/stable/basics/MTKLanguage/)
- [Connectors and Composition](https://docs.sciml.ai/ModelingToolkit/stable/basics/MTKModel_Connector/)

### AWE System Modeling
- Original implementation: `../generate_system.jl` (1400 lines of equation generation)
- System structure: `../system_structure.jl` (legacy data structures)
- Predefined models: `../predefined_structures.jl` (factory functions)

## 🚧 Implementation Status

### ✅ Completed (Phase 1)
- [x] Connectors (`MechanicalNode`, `ControlSignal`)
- [x] `MassPoint` with all dynamics types
- [x] `SpringDamperSegment` with drag
- [x] `Pulley` with dynamic/quasi-static modes
- [x] `Winch` with motor/friction dynamics
- [x] Example systems (pendulum, double tether, winch control)

### 🚧 Stubs (Phase 2 - To Be Implemented)
- [ ] `RigidWing`: Full 6-DOF quaternion dynamics
- [ ] `DeformableGroup`: Wing twist mechanics
- [ ] VSM integration: Aerodynamic force/moment calculation
- [ ] Attachment point transformations (body → world frame)

### 📋 Planned (Phase 3 - Future)
- [ ] System builder: `build_component_system(set, "ram")`
- [ ] Backward-compatible constructor: `SymbolicAWEModel(set; use_components=true)`
- [ ] State conversion: Components ↔ `SysState`
- [ ] Visualization integration: Makie plotting
- [ ] Advanced components: `ElasticWing`, `ActiveBridle`, etc.

## 🤝 Contributing

To contribute a new component:

1. Fork the repository
2. Create a branch: `git checkout -b feature/my-component`
3. Implement component following the patterns above
4. Add tests and examples
5. Submit pull request with clear description

## 📝 License

Copyright (c) 2025 Bart van de Lint
SPDX-License-Identifier: MPL-2.0

---

**Questions?** Open an issue on the [GitHub repository](https://github.com/OpenSourceAWE/SymbolicAWEModels.jl)
