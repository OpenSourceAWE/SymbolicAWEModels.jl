```@meta
CurrentModule = SymbolicAWEModels
```

# Compilation pipeline

SymbolicAWEModels works like a compiler: it takes a structural description and
transforms it through several stages into an efficient numerical ODE solver.
This page explains each stage.

## Overview

```
    Stage 1            Stage 2              Stage 3               Stage 4
┌──────────────┐   ┌──────────────┐   ┌──────────────────┐   ┌──────────────┐
│ Component    │──▶│ System       │──▶│ Symbolic Eqs     │──▶│ ODEProblem   │
│ Definition   │   │ Structure    │   │ ModelingToolkit  │   │ + Integrator │
│              │   │              │   │                  │   │              │
│ Point()      │   │ resolve refs │   │ point_eqs!()     │   │ init!()      │
│ Segment()    │   │ validate     │   │ segment_eqs!()   │   │ cache to .bin│
│ Wing()       │   │ compute COM  │   │ wing_eqs!()      │   │              │
│ Winch()      │   │              │   │ winch_eqs!()     │   │              │
│ ...          │   │              │   │ vsm_eqs!()       │   │              │
└──────────────┘   └──────────────┘   └──────────────────┘   └──────────────┘
                                            │
                                            ▼
                                  structural_simplify()
```

## Stage 1: component definition

Components are created using constructors ([`Point`](@ref), [`Segment`](@ref), etc.)
with symbolic name references. At this stage, references are unresolved — a segment
knows it connects `:anchor` to `:mass`, but doesn't know their numeric indices yet.

```julia
points = [
    Point(:anchor, [0, 0, 0], STATIC),
    Point(:mass, [0, 0, -50], DYNAMIC; extra_mass=1.0),
]
# :anchor and :mass are just names here, not yet resolved
segments = [Segment(:spring, :anchor, :mass,
    614600.0, 473.0, 0.004)]
```

Alternatively, components can be parsed from a YAML file using
[`load_sys_struct_from_yaml`](@ref), which calls the same constructors internally.

## Stage 2: SystemStructure assembly

The [`SystemStructure`](@ref) constructor takes the component vectors and:

1. **Assigns indices** — each component gets an `idx` based on its position in the
   vector (1, 2, 3, ...)
2. **Resolves references** — symbolic names like `:anchor` are mapped to indices via
   `assign_indices_and_resolve!()`
3. **Computes derived properties**:
   - Segment `l0` from point positions (if zero)
   - Wing center of mass from mass-weighted point centroids
   - Inertia tensor from point masses
   - VSM panel geometry adjustments
4. **Validates** — checks for NaN masses, zero stiffness, invalid pulley constraints,
   etc. via `validate_sys_struct()`

```julia
sys = SystemStructure("my_model", set; points, segments, transforms)
# All references now resolved: segments[1].point_idxs == (1, 2)
```

## Stage 3: symbolic equation generation

`create_sys!()` generates the full set of differential-algebraic equations (DAEs)
using ModelingToolkit.jl. It calls specialized equation builders for each subsystem:

| Function | Source file | Purpose |
|----------|-----------|---------|
| `point_eqs!()` | `src/generate_system/point_eqs.jl` | Newton's law for each point mass |
| `segment_eqs!()` | `src/generate_system/segment_eqs.jl` | Spring-damper forces with drag |
| `wing_eqs!()` | `src/generate_system/wing_eqs.jl` | Quaternion dynamics, angular momentum |
| `winch_eqs!()` | `src/generate_system/winch_eqs.jl` | Motor dynamics, Coulomb/viscous friction |
| `tether_eqs!()` | `src/generate_system/tether_eqs.jl` | Tether length kinematics |
| `pulley_eqs!()` | `src/generate_system/pulley_eqs.jl` | Equal-tension constraints |
| `group_eqs!()` | `src/generate_system/group_eqs.jl` | Twist deformation dynamics |
| `scalar_eqs!()` | `src/generate_system/scalar_eqs.jl` | Winch dynamics, kinematics |
| `vsm_eqs!()` | `src/generate_system/vsm_eqs.jl` | Linearized aerodynamics from VSM |

After generating all equations, `structural_simplify()` from ModelingToolkit reduces
the DAE system by eliminating algebraic constraints and identifying the minimal set
of independent variables.

## Stage 4: compilation and caching

[`init!`](@ref) creates the `ODEProblem` from the simplified symbolic system and
initializes the ODE integrator. This stage is expensive on first run because Julia
JIT-compiles the generated code.

The compiled system is serialized to a binary cache file
(`model_<julia_ver>_<name>_<wing_type>_...bin`) in the data directory. On subsequent
runs, the cached model is deserialized instead of recompiled, reducing startup from
minutes to seconds.

Force a rebuild by deleting the cache file or passing `remake_cache=true`.

## Stage 5: time-stepping

Once compiled, the simulation loop consists of:

1. **`next_step!(sam)`** — advances the ODE integrator by one time step
2. **`update_sys_struct!()`** — copies the integrator state back to the mutable
   component structs (point positions, wing orientation, etc.)
3. **`update_vsm!()`** — periodically calls the Vortex Step Method to update
   aerodynamic forces (controlled by `vsm_interval`)

```julia
using KiteUtils: next_step!

for i in 1:1000
    next_step!(sam; set_values=[torque])
end
```

[`sim!`](@ref) wraps this loop with a matrix of control inputs.

## Runtime parameter changes

The symbolic system uses **registered functions** (`@register_symbolic`) that read
from the mutable component structs at ODE evaluation time. This means many parameters
can be changed at runtime without recompiling:

- Winch parameters: `inertia_total`, `f_coulomb`, `c_vf`, `gear_ratio`
- Segment properties: `l0` (via tether/winch control)
- Wing damping: `body_frame_damping`, `world_frame_damping`
- VSM state: `vsm_jac`, `vsm_x`, `vsm_y` (updated by `update_vsm!()`)

Since registered functions read directly from the structs, changes take effect
instantly on the next ODE evaluation — no `init!` or `remake` call is needed.

Changes that **do** require recompilation (rebuilding the symbolic system):
- Adding or removing components (points, segments, wings)
- Changing the system topology (which points connect to which)
- Changing dynamics types (STATIC ↔ DYNAMIC)

These require creating a new [`SystemStructure`](@ref) and [`SymbolicAWEModel`](@ref).
