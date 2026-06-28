# Coordinate Frames

## Overview

SymbolicAWEModels uses three coordinate frames to describe geometry
and dynamics:

- **CAD frame (c)**: where geometry is originally defined
- **Body frame (b)**: attached to the wing, used for aerodynamics
- **World frame (w)**: the simulation frame

The transformation chain is:

```
             R_b_to_c, pos_cad          R_b_to_w, wing.pos_w
CAD frame ──────────────────▶ Body frame ─────────────────▶ World frame
(geometry)                    (wing-attached)               (simulation)
```

Each step involves both a **rotation** and a **translation**:
- **CAD to Body**: rotation `R_b_to_c` and origin shift to
  `pos_cad` (COM for RIGID_DYNAMICS, origin point for PARTICLE_DYNAMICS)
- **Body to World**: rotation `R_b_to_w` (from quaternion state or
  structural points) and translation to `wing.pos_w`

`R_b_to_c` is a constant rotation computed once during
[`SystemStructure`](@ref) construction.
`R_b_to_w` evolves during simulation — from the quaternion state
(RIGID_DYNAMICS) or from deformed point positions (PARTICLE_DYNAMICS).

## CAD Frame

The CAD frame is the coordinate system in which geometry is originally
defined, whether in a YAML file or via Julia constructors.

- Every point stores `pos_cad` — the original design position
- `pos_cad` is never modified by the codebase; it is kept as a
  permanent reference
- There is no imposed convention on orientation or origin — use
  whatever is convenient for your geometry
- Wing `pos_cad` is set to the centre of mass (RIGID_DYNAMICS) or the
  origin point position (PARTICLE_DYNAMICS) during construction
- VSM panel positions start in the CAD frame and are transformed to
  the body frame during construction

## Transform: CAD to World Initial Positioning

A [`Transform`](@ref) repositions CAD-frame geometry into the world
frame for the initial condition. Without a Transform,
`pos_w = pos_cad`.

When a Transform is applied, `reinit!` performs three steps:

1. **Translation**: `pos_w = pos_cad + (base_pos - curr_base_pos)`
2. **Rotation**: spherical repositioning using `elevation` and
   `azimuth` angles around the base point
3. **Heading**: orientation solve for wings (yaw about the radial
   axis)

This lets you place geometry defined in any convenient CAD orientation
into the correct world-frame position (e.g. a kite at 70deg
elevation).

```yaml
transforms:
  - name: tf
    elevation: -80.0      # degrees
    azimuth: 0.0
    heading: 0.0
    base_pos: [0, 0, 50]
    base_point: anchor
    rot_point: tip
```

The `base_point` is the reference point that gets placed at
`base_pos`. The `rot_point` (or `wing`) is what gets rotated to the
specified elevation and azimuth. Transforms can chain: use
`base_transform` instead of `base_pos` to use the already-rotated
`rot_point`/`wing` position of another transform as the base.

See `reinit!` in `transforms.jl`.

## World Frame

The world frame is the simulation-global coordinate system:

- **Origin**: ground station
- **Z-axis**: points up (positive upward)
- **X/Y axes**: define the horizontal plane
- Gravity acts in the `-Z` direction

All simulation quantities (`pos_w`, `vel_w`, forces) and the wind
vector are expressed in the world frame.

## Body Frame — RIGID_DYNAMICS

For `RIGID_DYNAMICS` wings, the body frame is **auto-computed**
as the principal-axis frame of the wing's point masses. This
diagonalizes the XZ-block of the inertia tensor, which simplifies
the rotational equations of motion.

### Algorithm

Given the set of `WING`-type points assigned to this wing:

1. **COM**: mass-weighted centroid in CAD frame
   ``\text{com} = \frac{\sum m_i \, \mathbf{p}_i}{\sum m_i}``
2. **Inertia tensor** ``I_\text{cad}`` about COM from point masses:
   ``I_\text{cad} = \sum m_i \left[
       (\mathbf{r}_i \cdot \mathbf{r}_i)\, \mathbf{I}_3
       - \mathbf{r}_i \mathbf{r}_i^\top \right]``
   where ``\mathbf{r}_i = \mathbf{p}_i - \text{com}``
3. **Y-rotation** ``R_y`` diagonalizes the XZ block:
   ``\theta = \tfrac{1}{2}\arctan\!\left(
       \frac{2\,I_{13}}{I_{11} - I_{33}}\right)``
4. **Body-to-CAD rotation**: ``R_{b \to c} = R_y^\top``
5. **Principal inertia**: ``I_\text{diag} = \text{diag}(
       R_y \, I_\text{cad} \, R_y^\top)``

Point positions in the body frame are:
``\mathbf{p}_b = R_{b \to c}^\top \,
    (\mathbf{p}_\text{cad} - \text{com})``

At runtime, the quaternion state ``Q_{b \to w}`` gives
``R_{b \to w}``, and world positions are recovered as:
``\mathbf{p}_w = \mathbf{wing.pos}_w +
    R_{b \to w} \, \mathbf{p}_b``

See `principal_frame` and the RIGID_DYNAMICS setup block in
`system_structure_core.jl`.

## Body Frame — PARTICLE_DYNAMICS

For `PARTICLE_DYNAMICS` wings, the user defines the body frame by
choosing structural reference points. This gives full control over
the frame orientation, which updates dynamically as the structure
deforms.

### Configuration

```yaml
wings:
  - dynamics_type: PARTICLE_DYNAMICS
    origin_idx: kcu
    z_ref_points: [kcu, le_center]
    y_ref_points: [le_right, le_left]
```

### Algorithm

Given the reference point positions in the world frame:

1. ``\mathbf{z} = \text{normalize}(
       \mathbf{p}_{z2} - \mathbf{p}_{z1})`` — body Z axis
2. ``\mathbf{y}_\text{temp} = \text{normalize}(
       \mathbf{p}_{y2} - \mathbf{p}_{y1})`` — approximate span
3. ``\mathbf{x} = \text{normalize}(
       \mathbf{y}_\text{temp} \times \mathbf{z})``
   — chord direction (orthogonal to Z)
4. ``\mathbf{y} = \mathbf{z} \times \mathbf{x}``
   — span direction (ensures right-handed frame)
5. ``R_{b \to w} = [\mathbf{x} \;\; \mathbf{y} \;\; \mathbf{z}]``
6. Origin = `pos_w[origin_idx]`

Key points:

- `z_ref_points` defines the body Z direction (e.g. kcu to le\_center
  gives a direction roughly along the tether, normal to the wing
  surface)
- `y_ref_points` defines the approximate span direction
- X is derived automatically as the orthogonal chord direction
- The frame is **recomputed each timestep** from current point
  positions, so it tracks structural deformation
- Different reference point choices produce different body frames —
  pick what makes physical sense for your model

See `calc_particle_dynamics_wing_frame` in `transforms.jl`.

## CAD to Body Transformation (VSM Panels)

Both wing types transform VSM panel positions from the CAD frame to
the body frame during [`SystemStructure`](@ref) construction:

1. **Translate**: subtract origin (`adjust_vsm_panels_to_origin!`)
2. **Rotate**: apply ``R_{b \to c}^\top`` to all section LE/TE points
   (`rotate_vsm_sections!`)
3. **Z-offset** (RIGID_DYNAMICS only): apply `aero_z_offset` to shift the
   aerodynamic reference vertically in the body frame
   (`apply_aero_z_offset!`)

After this transformation, all VSM geometry is expressed in the body
frame. During simulation, `R_b_to_w` maps panel positions to the world
frame for aerodynamic calculations.
