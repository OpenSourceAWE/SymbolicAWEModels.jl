```@meta
CurrentModule = SymbolicAWEModels
```

# SymbolicAWEModels.jl

Documentation for [SymbolicAWEModels.jl](https://github.com/OpenSourceAWE/SymbolicAWEModels.jl).

## What is SymbolicAWEModels.jl?

SymbolicAWEModels.jl is a **compiler** for mechanical systems, built for
**Airborne Wind Energy** (AWE) modelling. It takes a structural description
of a system — defined in Julia code or a YAML file — and compiles it into
an efficient ODE problem using [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl).

The compilation pipeline works as follows:

```
 Define Components         Assemble             Compile            Simulate
┌──────────────────┐   ┌──────────────┐    ┌─────────────────┐    ┌────────────┐
│ Point, Segment,  │──▶│ System       │───▶│ SymbolicAWE     │───▶│ init!()    │
│ Wing, Winch, ... │   │ Structure    │    │ Model           │    │ next_step! │
│                  │   │              │    │ (symbolic eqs → │    │ sim!()     │
│ Julia or YAML    │   │ (resolves    │    │  ODEProblem)    │    │            │
│                  │   │  references) │    │                 │    │            │
└──────────────────┘   └──────────────┘    └─────────────────┘    └────────────┘
```

The first compilation is slow (minutes) as ModelingToolkit generates and simplifies the
symbolic equations. The result is cached to a binary file, making subsequent runs fast
(seconds).

## Quick start

Install Julia using [juliaup](https://github.com/JuliaLang/juliaup):

```bash
curl -fsSL https://install.julialang.org | sh
juliaup add release
juliaup default release
```

Create a project and add SymbolicAWEModels:

```bash
mkdir my_project && cd my_project
julia --project="."
```

```julia
using Pkg
pkg"add SymbolicAWEModels"
```

### Minimal example (Julia)

```julia
using SymbolicAWEModels
using GLMakie

set = Settings("system.yaml")
set.v_wind = 0.0

# Define components using symbolic names
points = [
    Point(:anchor, [0, 0, 0], STATIC),
    Point(:mass, [0, 0, -50], DYNAMIC; extra_mass=1.0),
]
segments = [Segment(:spring, :anchor, :mass,
    614600.0, 473.0, 0.004)]
transforms = [Transform(:tf, deg2rad(-80), 0.0, 0.0;
    base_pos=[0, 0, 50], base_point=:anchor, rot_point=:mass)]

# Assemble and compile
sys = SystemStructure("pendulum", set; points, segments, transforms)
sam = SymbolicAWEModel(set, sys)
init!(sam)

# Simulate
for _ in 1:100
    next_step!(sam)
end

# Visualize the result
plot(sam.sys_struct)
```

For the full tutorial, see [Building a System using Julia](tutorial_julia.md).
For YAML-based model definition, see [Building a System using YAML](tutorial_yaml.md).

## What can it model?

SymbolicAWEModels provides building blocks for flexible mechanical systems:

- [`Point`](@ref) **masses** — static, dynamic, or quasi-static nodes
- [`Segment`](@ref) **spring-dampers** — with per-unit-length stiffness, damping, and drag
- [`Tether`](@ref)s — collections of segments controlled by a winch
- [`Winch`](@ref)es — torque-controlled motors with Coulomb and viscous friction
- [`Pulley`](@ref)s — equal-tension constraints between segments
- [`Wing`](@ref AbstractWing)s — rigid body quaternion dynamics with aerodynamic forces from the
  [Vortex Step Method](https://github.com/Albatross-Kite-Transport/VortexStepMethod.jl)
- [`Group`](@ref)s — twist degrees of freedom for aeroelastic coupling
- [`Transform`](@ref)s — spherical coordinate positioning of components

These components can be combined to model a wide range of systems, from simple
hanging masses to complex kite power systems with multiple tethers, bridles,
and wings.

## Ecosystem

Key related packages:
- [RamAirKite.jl](https://github.com/OpenSourceAWE/RamAirKite.jl) — ram air kite model
- [V3Kite.jl](https://github.com/OpenSourceAWE/V3Kite.jl) — TU Delft V3 kite model
- [KiteUtils.jl](https://github.com/OpenSourceAWE/KiteUtils.jl) — shared types and utilities
- [VortexStepMethod.jl](https://github.com/Albatross-Kite-Transport/VortexStepMethod.jl) — aerodynamic solver
- [AtmosphericModels.jl](https://github.com/aenarete/AtmosphericModels.jl) — wind profiles
- [KiteModels.jl](https://github.com/ufechner7/KiteModels.jl) — non-symbolic, predefined kite models
- [KiteSimulators.jl](https://github.com/aenarete/KiteSimulators.jl) — meta-package
- [KiteControllers.jl](https://github.com/aenarete/KiteControllers.jl) — control algorithms

Visualisation uses the built-in GLMakie extension
(`ext/SymbolicAWEModelsMakieExt.jl`) — just `using GLMakie` to enable
plotting.

## See also
- [Research Fechner](https://research.tudelft.nl/en/publications/?search=Fechner+wind&pageSize=50&ordering=rating&descending=true) for the scientific background of the winches and tethers

## Questions?
If you have questions or problems, please submit an
[issue](https://github.com/OpenSourceAWE/SymbolicAWEModels.jl/issues/new) or start a
[discussion](https://github.com/OpenSourceAWE/SymbolicAWEModels.jl/discussions/new/choose).
The Julia community is also very helpful:
[Julia Discourse](https://discourse.julialang.org/).

Authors: Bart van de Lint (bart@vandelint.net), Uwe Fechner (uwe.fechner.msc@gmail.com), Jelle Poland
