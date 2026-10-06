# Coordinate Frames

## Overview

SymbolicAWEModels uses three coordinate frames to describe geometry
and dynamics:

- **Body frame (b, KA)**: attached to the wing, used for aerodynamics
- **Principal frame (p)**: diagonal-inertia frame carrying the
  rigid-body ODE state; a constant rotation off the body frame
- **World frame (w, ENU)**: the simulation frame

The transformation chain is:

```
             R_b_to_w, wing.pos_w
Body frame ─────────────────────▶ World frame
(wing-attached)                   (simulation)
```

`R_b_to_w` evolves during simulation — from the quaternion state
(RIGID_DYNAMICS) or from deformed point positions (PARTICLE_DYNAMICS).

## Initial Pose

Every point and body holds an initial pose in the world frame: `pos_ENU`, and for a
body `Q_KA_to_ENU`, its body→world orientation. It is where a run starts from, and
a run never changes it; what moves is `pos_w` and `Q_b_to_w`.

- Everything derived from geometry is derived from the initial pose, by
  [`init!`](@ref) on every call: each body's mass properties and principal frame,
  each tube's rest geometry and each flap's rest deflection. A run restarted from a
  logged state therefore keeps the rest shape it started with.
- A wing's `pos_ENU` is its centre of mass (RIGID_DYNAMICS without reference
  points) or its weighted `origin` reference position.
- [`reset_to_initial_pose!`](@ref) puts the structure back at it.

## Transform: Authoring Geometry to World

Geometry is authored wherever convenient — the YAML's `pos_cad` column, or the
position passed to a constructor — and a [`Transform`](@ref) moves it into the
world. [`place!`](@ref) applies each Transform in three steps, then makes where the
structure lands its initial pose:

1. **Translation**: `pos_w = pos_ENU + (base_pos - curr_base_pos)`
2. **Rotation**: spherical repositioning using `elevation` and
   `azimuth` angles around the base point
3. **Heading**: orientation solve for wings (yaw about the radial
   axis)

A VSM wing's structure and its aerodynamic geometry must be authored in one frame: construction errors, naming the wing, when a station node or the mesh COM lies outside the sections' bounding box, or when the span from `y_ref_points` is turned away from the sections' span ([`check_aero_frame`](@ref SymbolicAWEModels.check_aero_frame) gives the margins).

Without a Transform, the authored position is the initial pose. This lets you
place geometry defined in any convenient orientation into the correct world-frame
position (e.g. a kite at 70deg elevation).

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

`SystemStructure` places once when it is built; after changing a Transform, call
`place!(sys_struct)` before `init!`. Placing again starts from the initial pose,
and the steps above land on the same pose whatever orientation they start from.

## World Frame

The world frame is the simulation-global coordinate system:

- **Origin**: ground station
- **Z-axis**: points up (positive upward)
- **X/Y axes**: define the horizontal plane
- Gravity acts in the `-Z` direction

All simulation quantities (`pos_w`, `vel_w`, forces) and the wind
vector are expressed in the world frame.

## Body Frame — RIGID_DYNAMICS

For `RIGID_DYNAMICS` wings the body frame is built the same way as
for `PARTICLE_DYNAMICS` — from user-chosen reference points (see
below) — but it is fitted once, in the authored geometry, instead of being
refitted every step, since the wing body is rigid. If a wing declares no
`origin`/`z_ref_points`/`y_ref_points`, the body frame keeps the authored
orientation with its origin at the wing body's own COM.

The wing's mass properties are those of the wing body with every point
it carries (see [Mass of a rigid body](@ref)):

1. **Own part**: `extra_mass` ``m_e`` with inertia ``I_e`` about its own
   COM ``\mathbf{c}_e``, spread like the `.obj` mesh when there is one,
   else like the frame points' `extra_mass`.
2. **Carried points**: each wing node and `BODY_STATIC` rider, as a
   point mass ``m_i`` (its `total_mass`) at its body-frame position
   ``\mathbf{p}_i``.
3. **COM**: ``\text{com\_offset}_b = \frac{m_e \mathbf{c}_e + \sum m_i \mathbf{p}_i}
   {m_e + \sum m_i}``, measured from the body origin (`wing.pos_ENU`,
   the weighted `origin` reference position), with each ``\mathbf{p}_i`` taken
   from the initial pose.
4. **Inertia** about that COM, by the parallel-axis theorem:
   ``I_b = I_e + m_e S(\mathbf{c}_e - \text{com}) + \sum m_i S(\mathbf{p}_i - \text{com})``
   with ``S(\mathbf{r}) = (\mathbf{r} \cdot \mathbf{r})\, \mathbf{I}_3 - \mathbf{r}\mathbf{r}^\top``.

At runtime, the quaternion state gives ``R_{b \to w}``, and world
positions are recovered as
``\mathbf{p}_w = \mathbf{wing.pos}_w + R_{b \to w} \, \mathbf{p}_b``.

See `setup_wing_frame!` and [`update_mass_properties!`](@ref) in
`system_structure_core.jl`.

## Principal Frame — RIGID_DYNAMICS

The rigid-body ODE state (`com_w`, `com_vel`, `Q_p_to_w`,
``\omega_p``) lives in the **principal frame (p)**, where the inertia
tensor is diagonal so the Euler equations have no product-of-inertia
terms. It is a constant rotation off the body frame,
``R_{b \to p} = R_{p \to c}^\top R_{b \to c}``.

The inertia it diagonalises is the body's with the points it carries
([`update_mass_properties!`](@ref)), expressed in the body frame, and
[`PrincipalFrameMethod`](@ref) selects how ``R_{b \to p}`` is found from it:

- `EIGEN_DECOMP` (`principal_frame`) — full 3-axis eigendecomposition
  with a permutation search. General-purpose, correct for any body.
- `Y_ROTATION` (`calc_inertia_y_rotation`) — closed-form rotation
  about Y only, diagonalizing the XZ block:
  ``\theta = \tfrac{1}{2}\arctan\!\left(
      \frac{2\,I_{13}}{I_{11} - I_{33}}\right)``.
  Use it for wings symmetric about the XZ-plane, where the generic
  permutation search is ambiguous when two principal moments are
  close.

The choice is a gauge: it changes the state representation, not the
physics.

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

## Aero Geometry to Body Transformation (VSM Panels)

The VSM aero geometry is read in the frame it was authored in, the same one as the
structure's authored geometry. Both wing types move it into the body frame during
[`SystemStructure`](@ref) construction, from the wing's authored origin and
orientation, which the VSM wing records as `T_cad_body` and `R_cad_body` so that a
rebuilt aero geometry lands the same way:

1. **Translate**: subtract origin (`adjust_vsm_panels_to_origin!`)
2. **Rotate**: apply the inverse of the wing's authored orientation to all section
   LE/TE points (`rotate_vsm_sections!`)
3. **Z-offset** (RIGID_DYNAMICS only): apply `aero_z_offset` to shift the
   aerodynamic reference vertically in the body frame
   (`apply_aero_z_offset!`)

After this transformation, all VSM geometry is expressed in the body
frame. During simulation, `R_b_to_w` maps panel positions to the world
frame for aerodynamic calculations.
