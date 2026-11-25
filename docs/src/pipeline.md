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
│ Point rows   │   │ resolve refs │   │ point_eqs!()     │   │ init!()      │
│ Segment rows │   │ validate     │   │ segment_eqs!()   │   │ cache to .bin│
│ Wing rows    │   │ compute COM  │   │ wing_eqs!()      │   │              │
│ Winch rows   │   │              │   │ winch_eqs!()     │   │              │
│ ...          │   │              │   │ aero_eqs!()      │   │              │
└──────────────┘   └──────────────┘   └──────────────────┘   └──────────────┘
                                            │
                                            ▼
                                        mtkcompile
```

## Stage 1: component definition

Components are KiteGeometry rows (`Point`, `Segment`, etc.), built with keywords, with
symbolic name references and a SymbolicAWEModels model ([`PointModel`](@ref),
[`SegmentSpring`](@ref), etc.) in their `model` column. At this stage, references are unresolved — a segment
knows it connects `:anchor` to `:mass`, but doesn't know their numeric indices yet.

```julia
materials = KiteGeometry.MaterialTable(KiteGeometry.Material[], Dict{String, Any}())
points = [
    Point(; name=:anchor, pos_ENU=[0, 0, 0], fixed=true),
    Point(; name=:mass, pos_ENU=[0, 0, -50], extra_mass=1.0),
]
# :anchor and :mass are just names here, not yet resolved
segments = [KiteGeometry.line_segment(:spring, (:anchor, :mass), set, materials;
    diameter=0.004, unit_stiffness=614600.0, model=SegmentSpring(; unit_damping=473.0))]
```

Alternatively, components can be parsed from a YAML file using
`load_authoring`, which builds the same rows.

## Stage 2: SystemDefinition assembly

[`system_definition`](@ref) takes the component vectors and:

1. **Builds the definition** — the rows and the `materials` they name become a
   KiteGeometry `SystemDefinition`; a winch without a `gear_ratio` or `drum_radius`
   takes `set`'s
2. **Completes and places the wings** — KiteGeometry's `complete_wings!` puts every
   wing's body in its design pose and gives it `set.mass` where it carries none, and
   KiteGeometry's `place!` places the whole
3. **Computes derived properties**:
   - Segment `l0` from point positions (if zero)
   - Wing center of mass and inertia tensor — from the `.obj` mesh when one is
     supplied, otherwise from the mass-weighted point masses
   - Body and principal frames ([`PrincipalFrameMethod`](@ref))
   - VSM panel geometry adjustments
4. **Validates** — checks for NaN masses, zero stiffness, invalid pulley constraints,
   etc. via `validate_sys_struct()`

```julia
sys = system_definition("my_model", set; points, segments, materials)
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
| `station_eqs!()` | `src/generate_system/station_eqs.jl` | Twist deformation dynamics |
| `rigid_body_eqs!()` | `src/generate_system/rigid_body_eqs.jl` | 6-DOF body dynamics in the principal frame |
| `body_eqs!()` | `src/generate_system/body_eqs.jl` | Body-frame pose outputs |
| `elastic_tube_eqs!()` | `src/generate_system/elastic_tube_eqs.jl` | [`ElasticTube`](@ref) 6-DOF springs |
| `timoshenko_tube_eqs!()` | `src/generate_system/timoshenko_tube_eqs.jl` | Corotational beam elements |
| `scalar_eqs!()` | `src/generate_system/scalar_eqs.jl` | Winch dynamics, kinematics |
| `aero_eqs!()` | `src/generate_system/aero_eqs.jl` | Wires each wing's aero component (mode-agnostic) |

After generating all equations, `mtkcompile` from ModelingToolkit reduces the
DAE system by eliminating algebraic constraints and identifying the minimal set
of independent variables.

## Stage 4: compilation and caching

[`init!`](@ref) creates the `ODEProblem` from the simplified symbolic system and
initializes the ODE integrator. This stage is expensive on first run because Julia
JIT-compiles the generated code.

The compiled system is serialized to a binary cache file
(`model_<pkg_ver>_<julia_ver>_<name>_<dynamics_type>_<aero_tag>_...bin`) in the data
directory. On subsequent runs, the cached model is deserialized instead of
recompiled, reducing startup from minutes to seconds. The package version is part
of the filename, so upgrading invalidates stale caches automatically.

Force a rebuild by deleting the cache file or passing `init!(sam; remake=true)`.
A model carrying a custom aero or winch component rebuilds by default — its
equations are not captured by the model hash.

## Stage 5: time-stepping

Once compiled, the simulation loop consists of:

1. **`next_step!(sam)`** — advances the ODE integrator by one time step
2. **`update_sys_struct!()`** — copies the integrator state back to the mutable
   component structs (point positions, wing orientation, etc.)
3. **`refresh_aero!()`** — periodically calls the Vortex Step Method to update
   aerodynamic forces (controlled by `vsm_interval`)

```julia
using KiteUtils: next_step!

for i in 1:1000
    next_step!(sam; set_values=[torque])
end
```

[`sim!`](@ref) wraps this loop with a matrix of control inputs.

## Runtime parameter changes

Numeric fields of the component structs are emitted as **flat MTK parameters** and
synced from the live `SystemDefinition` once per step, so many parameters
can be changed at runtime without recompiling:

- Winch parameters: `inertia_total`, `coulomb_friction`, `viscous_coefficient`, `gear_ratio`
- Segment properties: `l0` (via tether/winch control)
- Point damping: `body_frame_damping`, `world_frame_damping`
- Aero state: `aero_jac`, `aero_x`, `aero_y`, `aero_force_b` (updated by
  `refresh_aero!()`)

Mutate the struct field between steps and it is picked up at the next sync — no
`init!` or `remake` call is needed. Flat parameters compile to direct buffer
loads, which is why the generated RHS is allocation-free.

Changes that **do** require recompilation (rebuilding the symbolic system):
- Adding or removing components (points, segments, wings, bodies, tubes)
- Changing the system topology (which points connect to which)
- Changing dynamics types (STATIC ↔ DYNAMIC)
- Changing an aero mode, or any structural field it reports via `aero_hash_id`

These require creating a new `SystemDefinition` and [`SymbolicAWEModel`](@ref).
