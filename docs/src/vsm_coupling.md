```@meta
CurrentModule = SymbolicAWEModels
```

# VSM coupling

This document explains how SymbolicAWEModels couples with the
Vortex Step Method (VSM) for aerodynamic force computation. The
coupling is configured by two orthogonal choices:
[`WingType`](@ref) (structural representation) and
[`AbstractAeroModel`](@ref) (force computation strategy).

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

- Wing structure consists of ordinary [`DYNAMIC`](@ref DynamicsType)
  points organized in leading edge (LE) and trailing edge (TE)
  pairs. Their wing membership comes from twist-surface membership
  (`point.is_wing_node`), not from a dedicated point type
- Each consecutive pair (point i, point i+1) forms a structural
  segment (strut)
- Points can move independently — the wing can deform
  structurally
- Number of structural segments = (number of wing nodes) / 2

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
optional deformable twist surfaces.

#### Structural model

- Wing treated as a rigid body with quaternion-based orientation
- No per-point wing structure — aerodynamic forces applied to
  wing center of mass
- Optional **twist surfaces** represent deformable sections with twist
  degrees of freedom
- Twist surfaces control segment twist angles. With
  `use_prior_polar=true`, their LE/TE positions also define the
  aerodynamic section geometry

#### VSM mapping

VSM still uses unrefined sections, but they don't correspond to
individual structural points:

```
Unrefined sections: [Sec₁,  Sec₂,  Sec₃,  Sec₄,  Sec₅]
                        ╲       ╱       ╲       ╱
Twist surfaces:      [─ Surf₁ ─]    [─ Surf₂ ─]
                     twist DOF θ₁    twist DOF θ₂
```

Multiple unrefined sections can be combined into a single twist
surface for twist control. `compute_spatial_twist_surface_mapping!`
builds the mapping automatically: each unrefined section is assigned
to the nearest twist-surface centre (Voronoi partition in the body
frame), and the surface's single twist DOF then drives every section
it owns
as a rigid unit. More twist surfaces than unrefined sections is
rejected — a twist DOF without a section to drive would be undefined.

#### Force distribution

Integrated force and moment coefficients are applied to the rigid
body, driving quaternion dynamics. Each twist surface's aerodynamic
moment is the sum of its unrefined section moments, driving the twist
DOF.

## Aero modes

[`AbstractAeroModel`](@ref) controls **how aerodynamic forces enter the
ODE system** — orthogonal to the wing type choice.

### `AeroDirect`

A single VSM solve at the current operating point. The resulting
forces are stored in the wing struct and read back as flat
parameters during ODE evaluation:

1. `vsm_aero_coeffs` (Float64 path) sets the VSM body wind/ω
   from the live wing state and calls `VortexStepMethod.solve!`
2. `apply_direct_forces!` reconstructs physical forces in the
   wind-axis basis: `F = q∞ · A · (CL · lift + CD · drag + CS · side)`
3. Forces are stored in `wing.aero_force_b` and
   `wing.aero_moment_b` (RIGID_DYNAMICS) or per-point via
   `distribute_panel_forces_to_points!` (PARTICLE_DYNAMICS)
4. Between VSM updates (controlled by `vsm_interval`), forces
   are held constant

The component reads those buffers through its `params` view —
`params.wings[i].aero_force_b` / `.aero_moment_b` for a rigid wing,
`params.points[j].aero_force_b` for a particle wing — so the values
are synced from the live `SystemStructure` once per step without a
registered-function call in the RHS.

### `AeroLinearized`

First-order Taylor expansion of the VSM solver around the last
operating point, so the ODE RHS sees smooth force variations
between VSM updates without a nonlinear solve each call.

State:

- `aero_y = [α, β, ω₁, ω₂, ω₃, θ_g_1, …, θ_g_n]`
- `aero_x = [CL, CD, CS, CMx, CMy, CMz, cm_g_1, …, cm_g_n]`
- `aero_jac = ∂aero_x/∂aero_y` (dense)

Every `vsm_interval` steps, `refresh_aero!`:

1. Float64 VSM solve at `y0` to refresh `aero_x` and the
   converged circulation γ₀
2. ForwardDiff Jacobian via a lazily-allocated Dual-shadow
   solver, warm-started from γ₀ (1–2 Picard iters per column)
3. `safe_vsm_solve!` guards each solve, checking convergence
   and finiteness of both Dual values and partials (plain
   `isfinite(::Dual)` misses partial NaNs)

The ODE then reconstructs forces in the wind-axis basis
(`drag = va/|va|`, `lift = normalize(drag × span)`,
`side = lift × drag`) using
`coef_i = aero_x_0[i] + Σ_j aero_jac[i,j] · Δaero_y[j]`.

### `ContinuousAero`

Frozen-circulation VSM with the force assembly in the symbolic RHS.
The refresh runs only the circulation solve
(`VortexStepMethod.solve_base!`) and freezes each refined panel's
induced velocity; every RHS step re-derives panel geometry from the
live strut points (frozen mesh-interpolation weights), the effective
angle of attack, the polar coefficients and the lift/drag directions.
Forces therefore respond to wing motion *between* VSM updates —
aerodynamic damping through the changing angle of attack — unlike
`AeroDirect`'s piecewise-constant forces. All per-panel quantities
(`alpha`, `cl`, `q_dyn`, `panel_force`, …) are observable component
variables.

### `AeroPressure`

Continuous like `ContinuousAero` — the per-panel *total* force is
re-derived symbolically each RHS step from the frozen circulation —
but the force is **scattered onto arbitrary structural surface
points** rather than split LE/TE, using the airfoil surface traction
(``-C_p \hat{n} + c_f \hat{s}``) as the placement pattern. VSM owns
the totals; the traction pattern owns only placement and direction.
Each point gets its frozen traction plus an equal share of
(live panel force − frozen pattern net), so the point forces sum to
the live total exactly.

It also supports a live flap deflection ``\delta`` — the signed angle
between two structural bodies about a hinge, modelled by a
`KINEMATIC` [`TwistSurface`](@ref) — which feeds the ``(\alpha,
\delta)`` polars each RHS step.

### `AeroPlate`

Flat-plate lift and drag from a polar table (CL/CD vs α) evaluated by
registered symbolic interpolants, one 1-point `STATIC`
[`TwistSurface`](@ref) per plate. No VSM engine.

### `AeroNone`

Returns zero forces. Useful for debugging rigid body dynamics
without aerodynamic coupling. No VSM engine, so no `vsm_set` or VSM
geometry is needed.

### Compatibility

A mode supports a wing dynamics exactly when it defines the matching
[`aero_component`](@ref) method; anything else errors during the model
build.

| Aero mode | YAML `aero_mode` | `RIGID_DYNAMICS` | `PARTICLE_DYNAMICS` |
|-----------|------------------|:----------------:|:-------------------:|
| [`AeroLinearized`](@ref) | `linearized` | ✔ (default) | — |
| [`AeroDirect`](@ref) | `direct` | ✔ | ✔ (default) |
| [`ContinuousAero`](@ref) | `continuous` | — | ✔ |
| [`AeroPressure`](@ref) | `pressure` | — | ✔ |
| [`AeroPlate`](@ref) | `plate` | — | ✔ |
| [`AeroNone`](@ref) | `none` | ✔ | ✔ |

## Swappable aero components (dispatch)

Each wing carries an `aero::AbstractAeroModel` field. The builder is selected
by dispatch on **both** its type and the wing's dynamics,
[`aero_component`](@ref)`(mode, wing::AbstractWing, sys_struct; name, params)`,
returning a `System` exactly like a winch's [`Winch`](@ref) `model`. The
built-in subtypes [`AeroNone`](@ref), [`AeroDirect`](@ref),
[`AeroLinearized`](@ref) ship their own methods. A mode supports a wing
dynamics by defining the matching method, dispatched on [`RigidWing`](@ref)
and/or [`ParticleWing`](@ref) — rigid, particle, or both. To plug in your own
aerodynamics, subtype `AbstractAeroModel` and add a component builder (per
dynamics you support) and a cache tag:

```julia
struct MyAero <: AbstractAeroModel end

function SymbolicAWEModels.aero_component(::MyAero, wing::RigidWing, sys_struct; name, params)
    # ... build and return a System with the connectors below ...
end
SymbolicAWEModels.aero_mode_tag(::MyAero) = "myaero"

Wing(name, twist_surfaces, R_b_to_c, pos_cad, inertia; aero = MyAero())
```

Everything else is an **optional hook with a working default**, dispatched
on the mode:

- **Lifecycle**: [`setup_aero!`](@ref) (construction),
  [`remake_aero!`](@ref) (settings change), [`validate_aero_structure`](@ref)
  (build-time checks), [`resize_aero_state!`](@ref) (after name resolution),
  [`init_aero_state!`](@ref) (initial operating point).
- **Low-frequency refresh** (every `vsm_interval` steps, orchestrated by
  [`refresh_aero!`](@ref)): [`refresh_rigid_aero!`](@ref) /
  [`refresh_particle_aero!`](@ref).
- **Diagnostics**: [`calc_aoa`](@ref) (default `NaN`),
  [`normalized_inertia`](@ref) — per-unit-mass inertia [m²], scaled by the
  wing's mass at the single consumer (default: normalized point-mass inertia
  from the wing's structural points).
- **Log-point visualization**: [`n_aero_log_points`](@ref) /
  [`write_aero_log_points!`](@ref) / [`read_aero_log_points!`](@ref) /
  [`restore_aero_twist!`](@ref) — extra `SysState` slots for the mode's
  display geometry (defaults: none).
- **Live Makie rendering**: [`plot_wing_aero!`](@ref) /
  [`update_wing_aero_plot!`](@ref) — methods live in the Makie extension,
  so a custom mode draws with full Makie access (default: draws nothing).
- **Traits**: [`couples_to_sections`](@ref) (needs per-section twist
  surfaces; default `false`), [`provides_aero_override`](@ref),
  [`stores_point_force`](@ref), and the cache controls
  [`is_builtin_aero`](@ref) / [`aero_hash_id`](@ref) (see below).

Subtyping [`AbstractVSMAero`](@ref) (a [`VSMEngine`](@ref) in an `engine`
field, exposed via [`vsm_engine`](@ref)) inherits the VSM implementation of
every hook. There are no `isa`/`is_vsm` branches in the pipeline, so a
custom mode is never excluded from a code path it cannot extend.

The returned `System`'s connectors are fixed by the wing's dynamics type
([`RigidWing`](@ref)/[`ParticleWing`](@ref); all quantities in the wing
**body frame**):

- **`RIGID_DYNAMICS`** (`ng = length(wing.twist_surface_idxs)`):
  - inputs: `va[1:3]`, `rho`, `R_b_w[1:3,1:3]`, `omega[1:3]`,
    and — when `ng > 0` — `twist[1:ng]`, `twist_vel[1:ng]`
  - outputs: `force[1:3]`, `moment[1:3]`, `twist_moment[1:ng]`
- **`PARTICLE_DYNAMICS`** (`np` = number of wing nodes):
  - inputs: `point_pos[1:3,1:np]`, `point_vel[1:3,1:np]`,
    `va[1:3,1:np]`, `rho[1:np]`
  - outputs: `point_force[1:3,1:np]`

The wiring layer drives the inputs and reads the outputs; the component is
flattened by `mtkcompile`, so its connectors become inlined unknowns (no
array crosses a registered-function boundary). `validate_aero_component`
checks the contract at build time. A rigid component may additionally expose
an `aero_input` connector vector (as `AeroLinearized` does); it is detected by
name and logged as wing state — no extra method needed.

A custom model returns `false` from [`is_builtin_aero`](@ref) by default, so
its equations bypass the compiled-model cache and `init!` rebuilds (via
`has_custom_component`). Structural fields that change the *generated
equations* must be reported by [`aero_hash_id`](@ref); runtime-mutable fields
must not (they are read live, see below).

### Live-updating fields

Put the mutable value on the mode struct and read it through the build-time
`params` view: the field becomes a flat MTK parameter synced from the live
`SystemStructure` once per step (no `@register_symbolic`, no `psys`):

```julia
mutable struct ConstantLiftAero <: AbstractAeroModel
    CL::Float64                       # live-tunable
end

function SymbolicAWEModels.aero_component(::ConstantLiftAero,
                                          wing::RigidWing, sys_struct; name, params)
    CL = params.wings[wing.idx].aero.CL   # flat param, synced live each step
    # ... use CL in the force equation ...
end
```

Mutate the field between steps — no `remake`, picked up at the next sync:

```julia
sam.sys_struct.wings[1].aero.CL = 0.8
```

A numeric field becomes a scalar/array param; a **callable** field (an
interpolation or polar) becomes a callable param applied as `CL(α)` — see
[`ContinuousPolar`](@ref). A field that changes the equation *structure* is a
compile-time change (`init!(sam; remake=true)`) and belongs in
[`aero_hash_id`](@ref).

!!! note "Zero-allocation RHS"
    The built-in modes generate an allocation-free ODE RHS (asserted by
    `test_bench.jl`). A custom component is not tested in-package; to keep the
    RHS allocation-free, read data through the `params` view (flat params compile
    to direct buffer loads) rather than `@register_symbolic` getters, which box
    array arguments/returns and allocate. Check with the test helper
    `validate_rhs_allocs(sam; max_bytes=0, diagnose=true)` from `test/util.jl`.

## Aligning aero sections to structure

When the number of aerodynamic sections differs from the number
of structural LE/TE pairs, `match_aero_sections_to_structure!`
rebuilds the unrefined sections so their geometry matches the
structure. This applies to both wing types and requires
`use_prior_polar=true` on the VortexStepMethod wing.

The steps are:

1. **Find structural LE/TE pairs**: `identify_wing_segments`
   extracts pairs from twist surfaces (preferred) or uses a
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
   `aero_x`, and `aero_jac` are resized to match the new
   twist-surface count

## Refined panel mapping

Both wing types use `refined_panel_mapping` to handle VSM mesh
refinement:

### Purpose

VSM can subdivide unrefined sections into multiple refined panels
for higher aerodynamic fidelity. The mapping tracks which parent
unrefined section each refined panel belongs to.

### Computation

After VSM refinement, VortexStepMethod's
`compute_refined_panel_mapping!` finds the closest unrefined section
for each refined panel by comparing center positions:

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

1. **Twist-surface twist angles**: Applying the correct twist angle
   from twist surfaces to refined panels via their parent section
2. **Force distribution (PARTICLE_DYNAMICS)**: Accumulating refined panel
   forces at the structural points of their parent section
3. **Linearization (RIGID_DYNAMICS + [`AeroLinearized`](@ref))**:
   Propagating state perturbations through the correct sections

## Wing type summary

| Aspect | PARTICLE_DYNAMICS | RIGID_DYNAMICS |
|--------|--------|------------|
| **Structural repr.** | Individual `DYNAMIC` wing nodes | Rigid body + quaternion |
| **Section count** | = structural LE/TE pairs | Independent; optionally rebuilt via `use_prior_polar` |
| **Force distribution** | Per-point moment-preserving LE/TE split | Integrated force/moment on body |
| **Deformation** | Direct: point motion → VSM geometry | Indirect: twist-surface twists → sections |
| **Default aero mode** | [`AeroDirect`](@ref) | [`AeroLinearized`](@ref) |

## Implementation files

- `src/vsm_refine.jl`: Aero-to-structure alignment (all wing
  types), PARTICLE_DYNAMICS force distribution, and geometry updates
- `src/system_structure/types.jl`: Component type definitions
  including the `WingType` enum and `AbstractAeroModel` types
- `src/system_structure/wing.jl`: Wing and VSMWing definitions,
  twist-surface-to-section mapping
- `src/aero_modes/`: one file per aero mode (`none.jl`, `direct.jl`,
  `linearized.jl`, `continuous.jl`, `pressure.jl`, `plate.jl`), plus
  `common.jl` with the dispatch interface, the connector scaffolding
  and the shared VSM numerics
- `src/generate_system/aero_eqs.jl`: mode-agnostic wiring of the aero
  component into the rest of the system
- `src/generate_system/twist_surface_eqs.jl`: twist DOF equations
- `src/generate_system/wing_eqs.jl`: Wing dynamics equation
  generation
- `src/linearize.jl`: low-frequency refresh dispatch (`refresh_aero!`)
- VortexStepMethod.jl `src/wing_geometry.jl`:
  `refined_panel_mapping` computation
