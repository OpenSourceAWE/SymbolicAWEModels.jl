# VSM coupling

This document explains how SymbolicAWEModels couples with the
Vortex Step Method (VSM) for aerodynamic force computation. The
coupling is configured by two orthogonal choices:
[`WingType`](@ref) (structural representation) and
[`AeroMode`](@ref) (force computation strategy).

## Overview

Both wing types use **unrefined sections** as the fundamental
element that maps to VSM geometry:

- **Unrefined section**: A structural element defined by two
  points (leading edge and trailing edge) along the wing span
- **Refined panel**: VSM can subdivide unrefined sections into
  multiple panels for higher aerodynamic fidelity
- **`refined_panel_mapping`**: Maps each refined panel back to
  its parent unrefined section

The VSM solver computes aerodynamic coefficients (cl, cd, cm,
alpha) at both refined and unrefined levels. SymbolicAWEModels
uses the **unrefined-level coefficients** to map forces back to
the structural model.

## Wing types

[`WingType`](@ref) controls the **structural representation** of
the wing — how it deforms and how forces are distributed to
structural degrees of freedom.

### PARTICLE_DYNAMICS

The `PARTICLE_DYNAMICS` wing type creates the most direct coupling between
structure and aerodynamics.

#### Structural model

- Wing structure consists of [`WING`](@ref DynamicsType)-type
  points organized in leading edge (LE) and trailing edge (TE)
  pairs
- Each consecutive pair (point i, point i+1) forms a structural
  segment (strut)
- Points can move independently — the wing can deform
  structurally
- Number of structural segments = (number of WING points) / 2

#### VSM mapping

The structural segments map 1:1 to VSM unrefined sections:

```
Structural points:  [LE₁, TE₁]  [LE₂, TE₂]  [LE₃, TE₃]
                        ↓            ↓            ↓
Unrefined sections:   Sec₁         Sec₂         Sec₃
```

Each VSM section is defined by its LE and TE point positions,
taken directly from the structural point positions. VSM can
subdivide these into refined panels for higher fidelity;
`refined_panel_mapping` maps each refined panel back to its
parent section.

#### Geometry update

Each timestep, structural point positions update VSM section
geometry:

```julia
# For each structural point mapped to a VSM section point:
pos_b = R_b_to_w' * (point.pos_w - wing.origin)
section.LE_point = pos_b  # or section.TE_point
```

This bidirectional coupling allows structural deformation to
affect aerodynamics.

#### Force distribution

Per-panel forces are distributed to structural points:

1. **Panel → section mapping**: Each refined panel maps to its
   parent unrefined section via `refined_panel_mapping`
2. **Moment-preserving LE/TE split** (`compute_aerostruc_loads`):
   Each panel's force and moment are split into LE and TE
   contributions that preserve the total moment about a reference
   point
3. **Accumulation at structural points**: LE/TE forces are
   accumulated at the corresponding structural points via the
   `point_to_vsm_point` mapping

### RIGID_DYNAMICS

The `RIGID_DYNAMICS` wing type uses a rigid body representation with
optional deformable groups.

#### Structural model

- Wing treated as a rigid body with quaternion-based orientation
- No per-point wing structure — aerodynamic forces applied to
  wing center of mass
- Optional **groups** represent deformable sections with twist
  degrees of freedom
- Groups control segment twist angles. With
  `use_prior_polar=true`, group LE/TE positions also define the
  aerodynamic section geometry

#### VSM mapping

VSM still uses unrefined sections, but they don't correspond to
individual structural points:

```
Unrefined sections: [Sec₁,  Sec₂,  Sec₃,  Sec₄,  Sec₅]
                        ╲       ╱       ╲       ╱
Groups:              [─ Group₁ ─]    [─ Group₂ ─]
                     twist DOF θ₁    twist DOF θ₂
```

Multiple unrefined sections can be combined into a single group
for twist control. `compute_spatial_group_mapping!` builds the
mapping automatically: each unrefined section is assigned to the
nearest group centre (Voronoi partition in the body frame), and
the group's single twist DOF then drives every section it owns
as a rigid unit. `n_groups > n_unrefined` is rejected — a twist
DOF without a section to drive would be undefined.

#### Force distribution

Integrated force and moment coefficients are applied to the rigid
body, driving quaternion dynamics. Each group's aerodynamic moment
is the sum of its unrefined section moments, driving the twist
DOF.

## Aero modes

[`AeroMode`](@ref) controls **how aerodynamic forces enter the
ODE system** — orthogonal to the wing type choice.

### AERO_DIRECT

A single VSM solve at the current operating point. The resulting
forces are stored in the wing struct and read by registered
symbolic functions during ODE evaluation:

1. `_vsm_aero_coeffs` (Float64 path) sets the VSM body wind/ω
   from the live wing state and calls `VortexStepMethod.solve!`
2. `_apply_direct_forces!` reconstructs physical forces in the
   wind-axis basis: `F = q∞ · A · (CL · lift + CD · drag + CS · side)`
3. Forces are stored in `wing.aero_force_b` and
   `wing.aero_moment_b` (RIGID_DYNAMICS) or per-point via
   `distribute_panel_forces_to_points!` (PARTICLE_DYNAMICS)
4. Between VSM updates (controlled by `vsm_interval`), forces
   are held constant

For RIGID_DYNAMICS wings, the symbolic equations read the stored
forces via `get_aero_force_override` / `get_aero_moment_override`.
For PARTICLE_DYNAMICS wings, per-point forces are read via
`get_point_aero_force`.

### AERO_LINEARIZED

First-order Taylor expansion of the VSM solver around the last
operating point, so the ODE RHS sees smooth force variations
between VSM updates without a nonlinear solve each call.

State:

- `aero_y = [α, β, ω₁, ω₂, ω₃, θ_g_1, …, θ_g_n]`
- `aero_x = [CL, CD, CS, CMx, CMy, CMz, cm_g_1, …, cm_g_n]`
- `aero_jac = ∂aero_x/∂aero_y` (dense)

Every `vsm_interval` steps, `update_vsm!`:

1. Float64 VSM solve at `y0` to refresh `aero_x` and the
   converged circulation γ₀
2. ForwardDiff Jacobian via a lazily-allocated Dual-shadow
   solver, warm-started from γ₀ (1–2 Picard iters per column)
3. `_safe_vsm_solve!` guards each solve, checking convergence
   and finiteness of both Dual values and partials (plain
   `isfinite(::Dual)` misses partial NaNs)

The ODE then reconstructs forces in the wind-axis basis
(`drag = va/|va|`, `lift = normalize(drag × span)`,
`side = lift × drag`) using
`coef_i = aero_x_0[i] + Σ_j aero_jac[i,j] · Δaero_y[j]`.

### AERO_NONE

Returns zero forces. Useful for debugging rigid body dynamics
without aerodynamic coupling.

### Compatibility

| Wing type | Default aero mode | Supported modes |
|-----------|-------------------|-----------------|
| **RIGID_DYNAMICS** | `AERO_LINEARIZED` | `AERO_LINEARIZED`, `AERO_DIRECT`, `AERO_NONE` |
| **PARTICLE_DYNAMICS** | `AERO_DIRECT` | `AERO_DIRECT`, `AERO_NONE` |

`PARTICLE_DYNAMICS` + `AERO_LINEARIZED` is not yet implemented (raises an
error at runtime).

## Aligning aero sections to structure

When the number of aerodynamic sections differs from the number
of structural LE/TE pairs, `match_aero_sections_to_structure!`
rebuilds the unrefined sections so their geometry matches the
structure. This applies to both wing types and requires
`use_prior_polar=true` on the VortexStepMethod wing.

The steps are:

1. **Find structural LE/TE pairs**: `identify_wing_segments`
   extracts pairs from groups (preferred) or uses a
   consecutive-pair heuristic
2. **Rebuild unrefined sections**: For each structural pair, a
   new `Section` is created with LE/TE positions from the
   structural points (in body frame). Its airfoil data
   (`aero_model`, `aero_data`) is copied from the nearest
   original unrefined section by span index
3. **Re-refine**: `refine!` updates refined panel geometry from
   the rebuilt unrefined sections. Because `use_prior_polar=true`
   and `n_panels` is unchanged, existing refined panel polars are
   preserved — only positions are re-interpolated
4. **Resize linearization state**: For non-PARTICLE_DYNAMICS wings, `aero_y`,
   `aero_x`, and `aero_jac` are resized to match the new group
   count

## Refined panel mapping

Both wing types use `refined_panel_mapping` to handle VSM mesh
refinement:

### Purpose

VSM can subdivide unrefined sections into multiple refined panels
for higher aerodynamic fidelity. The mapping tracks which parent
unrefined section each refined panel belongs to.

### Computation

After VSM refinement, `compute_refined_panel_mapping!` finds the
closest unrefined section for each refined panel by comparing
center positions:

```julia
for each refined_panel in wing.refined_sections
    center = compute_center(refined_panel)
    closest_section = argmin(
        distance(center, unrefined_section_centers))
    refined_panel_mapping[refined_panel_idx] = closest_section
end
```

### Usage

The mapping enables:

1. **Group twist angles**: Applying the correct twist angle from
   groups to refined panels via their parent section
2. **Force distribution (PARTICLE_DYNAMICS)**: Accumulating refined panel
   forces at the structural points of their parent section
3. **Linearization (RIGID_DYNAMICS + AERO_LINEARIZED)**: Propagating
   state perturbations through the correct sections

## Wing type summary

| Aspect | PARTICLE_DYNAMICS | RIGID_DYNAMICS |
|--------|--------|------------|
| **Structural repr.** | Individual WING points | Rigid body + quaternion |
| **Section count** | = structural LE/TE pairs | Independent; optionally rebuilt via `use_prior_polar` |
| **Force distribution** | Per-point moment-preserving LE/TE split | Integrated force/moment on body |
| **Deformation** | Direct: point motion → VSM geometry | Indirect: group twists → sections |
| **Default aero mode** | `AERO_DIRECT` | `AERO_LINEARIZED` |

## Implementation files

- `src/vsm_refine.jl`: Aero-to-structure alignment (all wing
  types), PARTICLE_DYNAMICS force distribution, and geometry updates
- `src/system_structure/types.jl`: Component type definitions
  including `WingType` and `AeroMode` enums
- `src/system_structure/wing.jl`: Wing and VSMWing type
  definitions, group-to-section mapping
- `src/generate_system/vsm_eqs.jl`: Symbolic VSM equation
  generation (all wing type × aero mode combinations)
- `src/generate_system/wing_eqs.jl`: Wing dynamics equation
  generation
- `src/linearize.jl`: VSM update dispatch — linearization
  (RIGID_DYNAMICS) and nonlinear solve (PARTICLE_DYNAMICS)
- VortexStepMethod.jl `src/wing_geometry.jl`:
  `refined_panel_mapping` computation
