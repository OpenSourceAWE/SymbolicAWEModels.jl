```@meta
CurrentModule = SymbolicAWEModels
```

# Building a system using YAML

This tutorial explains how to define mechanical systems using YAML configuration files.
YAML is the recommended approach for complex models with many components, since it
separates geometry data from simulation code.

## Overview

The YAML workflow has three steps:

1. **Write a YAML file** — define materials, points, segments, and other components in
   a structured text file
2. **Load with [`load_sys_struct_from_yaml`](@ref)** — parses the YAML and calls the
   same Julia constructors used in the [Julia tutorial](tutorial_julia.md)
3. **Compile and simulate** — same as the Julia path: [`SymbolicAWEModel`](@ref) →
   [`init!`](@ref) → [`next_step!`](@ref)

The YAML loader does as little as possible: it parses YAML, converts string enum values,
resolves material references, and calls constructors. All defaults and derived
calculations happen in the component constructors.

## YAML file structure

A YAML file can contain any of these top-level blocks:

| Block | Purpose |
|-------|---------|
| `materials` | Material property lookup table |
| `points` | Point masses (nodes in the system) |
| `segments` | Spring-damper connections |
| `pulleys` | Equal-tension constraints |
| `groups` | Deformable wing sections |
| `tethers` | Winch-controlled segment groups |
| `winches` | Torque-controlled motors |
| `wings` | Aerodynamic bodies |
| `transforms` | Spherical coordinate positioning |

Each block uses a **headers + data** format:

```yaml
points:
  headers: [name, pos_cad, type, extra_mass]
  data:
    - [anchor, [0, 0, 0], STATIC, 0.0]
    - [mass, [0, 0, -50], DYNAMIC, 1.0]
```

The `headers` row defines column names. Each `data` row is a list of values matching
those headers. Missing trailing columns default to `nothing`.

Alternatively, you can use a **dict format** where each row is a dictionary:

```yaml
points:
  data:
    - {name: anchor, pos_cad: [0, 0, 0], type: STATIC}
    - {name: mass, pos_cad: [0, 0, -50], type: DYNAMIC, extra_mass: 1.0}
```

## Minimal example

Here is a complete YAML file for a simple two-point tether:

```yaml
# simple_tether.yaml

points:
  headers: [idx, pos_cad, type, wing_idx, transform_idx, extra_mass]
  data:
    - [1, [0, 0, 0], STATIC, nothing, 1, 0.0]
    - [2, [0, 0, -50], DYNAMIC, nothing, 1, 1.0]

segments:
  headers: [idx, point_i, point_j, l0, diameter_mm,
            unit_stiffness, unit_damping, compression_frac]
  data:
    - [1, 1, 2, 50.0, 5.0, 100000, 50.0, 0.001]

transforms:
  headers: [idx, elevation, azimuth, heading,
            base_pos, base_point_idx, rot_point_idx]
  data:
    - [1, -80, 0, 0, [0, 0, 50], 1, 2]
```

Load and simulate:

```julia
using SymbolicAWEModels

set = Settings("system.yaml")
set.v_wind = 0.0

sys = load_sys_struct_from_yaml("simple_tether.yaml";
    system_name="simple_tether", set=set)
sam = SymbolicAWEModel(set, sys)
init!(sam)

for _ in 1:100
    next_step!(sam)
end
```

!!! note "Transform angles"
    Transform elevation and azimuth values in YAML are specified in **degrees**
    (converted automatically), unlike the Julia constructor which takes radians.

## Materials and references

Materials allow you to define physical properties once and reference them across
segments. When a segment's `unit_stiffness` column contains a string instead of a
number, it is treated as a material reference.

```yaml
materials:
  headers: [name, youngs_modulus, density, damping_per_stiffness]
  data:
    - [dyneema, 55e9, 724, 0.00077]

segments:
  headers: [idx, point_i, point_j, l0, diameter_mm,
            unit_stiffness, unit_damping, compression_frac]
  data:
    # 'dyneema' in unit_stiffness triggers material lookup
    - [1, 1, 2, 5.0, 5.0, dyneema, nothing, 0.01]
    # Explicit stiffness (no material lookup)
    - [2, 2, 3, 5.0, 1.0, 100000, 50.0, 0.01]
```

When a material is referenced, derived properties are calculated automatically:

- **`unit_stiffness`** = `youngs_modulus * pi * (diameter_mm/2000)^2`
- **`unit_damping`** = `damping_per_stiffness * unit_stiffness`

## Component reference

### Points

```yaml
points:
  headers: [idx, pos_cad, type, wing_idx, transform_idx,
            extra_mass, body_frame_damping, world_frame_damping,
            area, drag_coeff]
```

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `idx` | Int | required | Point identifier |
| `pos_cad` | [x,y,z] | required | Position in CAD frame [m] |
| `type` | String | required | `STATIC`, `DYNAMIC`, `QUASI_STATIC`, or `WING` |
| `wing_idx` | Int/nothing | 1 | Wing this point belongs to |
| `transform_idx` | Int/nothing | nothing | Transform for initial positioning |
| `extra_mass` | Float | 0.0 | Additional mass [kg] |
| `body_frame_damping` | Float | 0.0 | Damping in body frame [Ns/m] |
| `world_frame_damping` | Float | 0.0 | Damping in world frame [Ns/m] |
| `area` | Float | 0.0 | Cross-sectional area for drag [m^2] |
| `drag_coeff` | Float | 0.0 | Drag coefficient |

### Segments

Two formats are supported:

**With explicit stiffness:**
```yaml
segments:
  headers: [idx, point_i, point_j, l0, diameter_mm,
            unit_stiffness, unit_damping, compression_frac]
  data:
    - [1, 1, 2, 5.0, 5.0, 100000, 50.0, 0.01]
```

**With material reference:**
```yaml
segments:
  headers: [idx, point_i, point_j, l0, diameter_mm,
            unit_stiffness, unit_damping, compression_frac]
  data:
    - [1, 1, 2, 5.0, 5.0, dyneema, nothing, 0.01]
```

| Field | Type | Description |
|-------|------|-------------|
| `idx` | Int | Segment identifier |
| `point_i`, `point_j` | Int | Endpoint point indices |
| `l0` | Float | Unstretched length [m] (0 = calculate from points) |
| `diameter_mm` | Float | Diameter [mm] |
| `unit_stiffness` | Float/String | Per-unit-length stiffness [N], or material name |
| `unit_damping` | Float/nothing | Per-unit-length damping [Ns], or nothing for auto |
| `compression_frac` | Float | Compressive/tensile stiffness ratio (0-1) |

### Pulleys

```yaml
pulleys:
  headers: [idx, segment_i, segment_j, type]
  data:
    - [1, 3, 4, DYNAMIC]
```

### Tethers

**Route 1** (explicit segments):
```yaml
tethers:
  headers: [idx, segment_idxs]
  data:
    - [1, [1, 2, 3]]
```

**Route 2** (auto-generated segments):
```yaml
tethers:
  headers: [name, start_point, end_point, n_segments]
  data:
    - [main, kite, ground, 5]
```

### Winches

```yaml
winches:
  headers: [idx, tether_idxs, winch_point]
  data:
    - [1, [1], ground]
```

### Transforms

```yaml
transforms:
  headers: [idx, elevation, azimuth, heading,
            base_pos, base_point_idx, rot_point_idx]
  data:
    - [1, -80, 0, 0, [0, 0, 50], 1, 2]
```

## Loading workflow

The full loading workflow for a model with aerodynamics:

```julia
using SymbolicAWEModels, VortexStepMethod

set_data_path("data/2plate_kite")
set = Settings("system.yaml")
vsm_set = VortexStepMethod.VSMSettings(
    joinpath(get_data_path(), "vsm_settings.yaml"); data_prefix=false)

struc_yaml = joinpath(get_data_path(),
    "quat_struc_geometry.yaml")
sys = load_sys_struct_from_yaml(struc_yaml;
    system_name="2plate_kite",
    set=set,
    vsm_set=vsm_set)

sam = SymbolicAWEModel(set, sys)
init!(sam)
```

![2-plate kite structure](assets/2plate_kite_structure.png)

After compilation, a cache file (`model_*.bin`) is saved. Subsequent loads skip the
expensive symbolic compilation and deserialize the cached model instead. Force a
rebuild with `remake_cache=true`.

## YAML vs Julia

| Aspect | YAML | Julia constructors |
|--------|------|--------------------|
| **Best for** | Complex models, data from CAD/measurements | Simple models, programmatic generation |
| **Readability** | Easy to scan geometry at a glance | Better for computed geometry (loops, formulas) |
| **Material refs** | Built-in: reference by name | Manual: pass stiffness/damping directly |
| **Version control** | Clean diffs for parameter changes | Code diffs mix logic and parameters |

Both paths produce the same [`SystemStructure`](@ref) type and are equally
capable. They can be freely mixed — for example, load a YAML model and then
modify component fields in Julia before simulation.
