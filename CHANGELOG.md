# CHANGELOG

## Unreleased

### Added
- `Segment(...; compression_damping_frac)` scales the damping term of a segment
  the way `compression_frac` scales the stiffness, live in the right-hand side.
  The default `1.0` is the previous behaviour, where damping is unaffected by
  compression and the force is continuous at `len == l0`, but the damping ratio
  jumps by `1/sqrt(compression_frac)` as the segment goes slack. Setting it equal
  to `compression_frac` gives one damping ratio on both branches, and
  `compression_damping_frac = compression_frac = 0` makes a slack segment carry
  no force at all, at the cost of a force step of
  `unit_damping / l0 · spring_vel` at the crossing. Also readable as a
  `segments` YAML column.
- `Tether(...; compression_frac, compression_damping_frac)` hands both to every
  segment a Route 2 tether generates, so an auto-generated line is no longer
  stuck at the `Segment` default of `compression_frac = 0.1`. Also readable as
  `tethers` YAML columns. A tether without a winch holds its length fixed, so
  this is the way to split one line into several segments without writing every
  segment and intermediate point out by hand.
- A `SysLog` now carries the whole differential state, so a simulation can be
  restarted from any logged step. `update_sys_state!` writes point velocities,
  per-frame body turn rates, twist_surface twist and rate, and pulley length and
  rate alongside the positions and orientations it already wrote, and
  `update_from_sysstate!` reads all of them back — including the body rates,
  which it rebuilds into the principal frame the ODE actually integrates. It no
  longer zeroes point velocities or caps twist surfaces at four. Pass
  `precision=Float64` to `SysState(sam)` or `Logger(sam, steps)` to log a state
  that reproduces `integrator.u`; the default stays `Float32`, so ordinary
  telemetry logs keep their size. Restart with
  `init!(sam; remake=false, reinit_sys=false)`, which seeds the integrator from
  the structure instead of resetting it from the CAD frame.
- `test_init!` round-trips every model the suite builds through a `SysLog` on
  disk and checks `integrator.u` is unchanged, which covers the `DynamicsType`
  combinations that decide whether a component is part of the state at all.
- The YAML loader reads back state that previously had no column: point `vel_w`
  (it was zeroed unconditionally, discarding the column it already accepted),
  pulley `sum_len`/`len`/`vel`, twist_surface `twist`/`twist_vel`, tether `len`,
  winch `vel`/`set_value`, and body `vel`/`Q_b_to_w`/`omega_b`.
- `init!(sam; analytic_jacobian)` gives the solver a Jacobian instead of leaving
  it to differentiate the right-hand side numerically. On the `KernelBackend` it
  is on by default: the right-hand side is a layered composition of small
  components, so each kernel is differentiated once at its own width and the
  results are composed through the constant wiring, rather than differentiating
  the whole model `n_states / chunk` times. On the `MonolithBackend` it is off by
  default and selects MTK's symbolic `jac=true`, which is expensive to build and
  has no derivative for a registered numerical leaf such as the wind profile.
  `nothing` (the default) takes the backend's `default_analytic_jacobian`. The
  choice is part of the serialized model's name, so the two builds are cached
  apart.
- New top-level YAML block `variables`: named values reused across the file.
  A variable name written in any column is replaced by its value (numbers,
  strings, lists), and variables may be defined in terms of each other:
  ```yaml
  variables:
    bridle_comp: 0.01
    dyneema: {youngs_modulus: 55.0e9, damping_per_stiffness: 0.00077,
              density: 724.0}
  ```
  A variable holding a mapping is a *multi-variable*: it fills the columns it
  names at once, so the row carries one entry for the whole group
  (`- [bridle, le_left, kcu, nothing, 1.0, dyneema, bridle_comp]`). The fields
  must match the columns at that position; in a dict row they are merged in
  instead, without overriding the row. Naming a variable after a component is
  an error.
- `youngs_modulus` [Pa] and `damping_per_stiffness` [s] are now ordinary
  `Segment` and `Tether` inputs (YAML columns and constructor keywords). They
  are the diameter-independent forms of `unit_stiffness` and `unit_damping`,
  so one material can be shared by elements of different diameter:
  `unit_stiffness = youngs_modulus * pi * (diameter_mm/2000)^2` and
  `unit_damping = damping_per_stiffness * unit_stiffness`. Giving both forms
  of the same quantity is an error.
- `Pulley` also takes `damping` [N·s/m] and `brake`, both defaulting to off and
  neither a sheave property: `damping` opposes rope travel proportionally to
  speed, to settle a ringing split while debugging, and `brake` freezes the
  split where it is, to isolate whether a problem comes from the rope
  redistributing at all. Both are settable from the constructor and from YAML
  columns of the same names.

### Fixed
- Air density is now clamped at ground level on every path. A rigid wing read
  `calc_rho` at its raw world `z`, so a wing that dipped below `z = 0`
  extrapolated the atmospheric model backwards, while points and segments
  clamped. The clamp now lives in one place, `air_density`, which every
  consumer in both backends reads through.
- `init_stretched_length` placement now carries a beam wing. A `BODY_STATIC`
  point riding a Timoshenko element has `body_idx == 0` and holds its
  association in `joint_idx`, so collecting the bodies to translate by
  `body_idx` alone left every beam body behind while the tether's free end
  moved. The ride constraint then snapped the points back, and the only symptom
  was a non-converged VSM solve much later. Bodies reachable from the moved ones
  through the beam graph now translate too — including bodies that carry no
  point of their own — stopping at `STATIC` bodies, since a beam with a clamped
  end deforms rather than translating.
- The Breukels inflated-tube correlations now error instead of quietly
  returning a negative rigidity when evaluated outside the range they were
  fitted in. The bending slope changes sign below a radius of roughly 38 mm at
  0.25 bar (a 20 mm tube at 0.1 bar gave `EI0 = -94 N·m²` and a negative `EA`),
  and the torsion factor `c2` changes sign for a fat tube above 1 bar. A
  non-positive radius or pressure is rejected as well.
- `set.g_earth` is settable after construction on the `KernelBackend` too. A
  particle point baked the value in while a rigid body and a ride point read the
  registered parameter, so changing gravity on an existing model moved the bodies
  and left the tether points falling at the old rate.
- An aero panel that maps to no twist surface reads a zero deflection on both
  backends. The `KernelBackend` already did; the monolith indexed the
  `twist_surface_delta` array at `0`. Unreachable today, because every panel is
  assigned a nearest flap whenever a wing has one, but live for any mode that
  emits a partially-zero map.

### Changed
- The monolith and the `KernelBackend` evaluate one shared definition of every
  equation, instead of a copy each. Points, segments, pulleys, rigid bodies,
  ref-point wing frames, twist deformation, the ground wind fallback and the aero
  wiring now route through the helpers in `src/components/`, which previously had
  the `KernelBackend` as their only caller. A backend still chooses how equations
  are assembled; it no longer chooses what they are.
- The aero tests are layered, so a new mode inherits its contracts by being
  added to a table rather than by someone remembering to write tests.
  `test_aero_modes.jl` holds what every mode owes and now enrolls
  `AeroPressure`; the new `test_continuous_modes.jl` holds what the continuous
  modes owe — including that per-section inflow is gathered per section and
  never as a wing-wide mean, which is the bug fixed above; and
  `test_continuous_aero.jl` and `test_pressure_aero.jl` keep only what is
  specific to one mode. `AeroPressure` previously sat in none of these layers.
- Nine segment observables (`stiffness`, `damping`, `segment_height`,
  `segment_vel`, `segment_rho`, `wind_vel`, `va`, `area`, `app_perp_vel`) are
  gone, and the point observable `point_drag_force` is renamed
  `point_aero_drag`, where it shadowed the shared function of that name. They
  were intermediates of the monolith's own copy of the spring and drag law,
  which no longer exists; keeping them would have meant re-deriving the force.
  Nothing in the package, the examples or the docs reads them and no exported
  accessor returns a symbolic variable by name, so only code indexing the
  generated system directly is affected. `len`, `spring_force`,
  `spring_force_vec`, `spring_vel`, `unit_vec`, `segment_vec`, `rel_vel`, `l0`
  and `total_drag` are unchanged.

### Breaking
- A segment's damper now resists the rate of change of its extension
  `len - l0`, not the rate its two endpoints separate. The two agree wherever
  `l0` is a fixed parameter, so an ordinary segment is unchanged — but a pulley
  leg and a winched tether member carry `l0` as a state, and the difference is
  that state's own rate. A pulley's rope split therefore appeared in no damper
  at all, while the legs' dampers still drove it: a one-way velocity coupling
  that could feed the split rather than settle it, growing with
  `damping_per_stiffness`, so raising the damping could make such a model
  *less* stable. The split now carries `c_left + c_right`. Reeling likewise no
  longer charges a tether a damper force proportional to reel speed, so a
  reel-out loop is genuinely less damped than the old model made it look —
  expect controllers tuned against it to need retuning.
- A `Pulley` resists rope travel by its `efficiency`, the fraction of line
  tension its sheave passes on. The friction is `(1 − efficiency) ·
  line_tension` carrying the sign of the motion, smoothed over
  `friction_epsilon` [m/s], with `line_tension` the mean of the two leg
  tensions — so the default is the textbook `T_out = 0.95 · T_in`. It replaces
  the hidden `pulley_damp` parameter, which was fixed at 5.0 and could be
  neither read nor set from the `SystemStructure`. `efficiency` defaults to
  0.95, a sealed ball-bearing sheave; published ranges are 0.94–0.97 for those
  and 0.88–0.92 for a bronze bushing. Set 1.0 for an ideal pulley.

  A sheave's losses scale with load, not speed: bearing drag rises with the
  force on the axle and the rope's bending hysteresis with the tension being
  bent. `pulley_damp` was proportional to rope travel and to rope mass, so no
  fixed coefficient reproduces it across models — expect to retune ones that
  leaned on it to settle. A slack pulley is now frictionless, and coasts, but
  it has no tension driving its split either.
- `Winch.f_coulomb` and `Winch.c_vf` are renamed to `Winch.coulomb_friction` and
  `Winch.viscous_coefficient`, matching the `Pulley` fields above and saying
  which of the two is a force [N] and which a coefficient [N·s/m]. The
  `settings.yaml` keys are owned by `KiteUtils` and keep their old names.
- `compression_frac` now defaults to 0.1 in both `Segment` constructors. They
  disagreed — 0.0 from the settings-based one, 0.1 from the direct one — which
  was an oversight rather than a choice. 0.1 lets a rope push back with a tenth
  of its tensile stiffness, and is the value models have actually been built
  against; 0.0 leaves a slack segment with damping but no stiffness at all,
  which a pulley model does not survive.
- The `materials`, `elements` and `segment_properties` YAML blocks were
  removed. A material is now a multi-variable listing `youngs_modulus`,
  `damping_per_stiffness` and `density`, and those are ordinary columns.
  Loading a file with one of the removed blocks errors.
- A segment whose `unit_damping` is missing or `nothing` now takes the
  `Segment` constructor default (derived from the settings) instead of being
  silently forced to zero by the loader.

## v0.13.0 06-08-2026

### Added
- New aero mode `ContinuousAero` (`PARTICLE_DYNAMICS`, YAML
  `aero_mode: continuous`): frozen-circulation VSM with the full force
  assembly in the symbolic RHS. The low-frequency refresh runs only the
  circulation solve (`solve_base!`) and freezes each refined panel's induced
  velocity (`AIC·γ`); every RHS step re-derives panel geometry from the live
  strut points (frozen mesh-interpolation weights), effective angle of
  attack, polar coefficients (registered `Dual`-safe lookups on the panel
  polars), lift/drag directions, and forces. Forces therefore respond to
  wing motion between VSM updates — aerodynamic damping through the changing
  angle of attack — unlike `AeroDirect`'s piecewise-constant forces. All
  per-panel quantities (`alpha`, `cl`, `q_dyn`, `panel_force`, …) are
  observable component variables. The mesh weights enter the model-cache
  hash via `aero_hash_id`.
- New aero mode `AeroPressure` (YAML `aero_mode: pressure`): like
  `ContinuousAero`, the VSM solve stays discrete (every `vsm_interval`) but
  the per-panel force is re-derived symbolically each RHS step from the
  frozen circulation, then *scattered onto arbitrary `BODY_STATIC` surface
  points* via a build-time station→point map. Gives aero damping between
  refreshes.
- Deflecting flap for `AeroPressure`: a live flap deflection `δ` (the signed
  angle between two structural bodies about a hinge, relative to rest) feeds
  the `(α, δ)` polars each RHS step. Modeled by extending `TwistSurface` with
  a `KINEMATIC` twist source and flap fields (`flap_body_idxs`, `flap_axis`,
  `flap_chord_refs`, `flap_rest_delta`) — no new component. Needs
  VortexStepMethod's 3-arg `calculate_cl/cd/cm(panel, α, δ)`.
- Structural `Body` / beam subsystem: a rigid `Body` type with new
  `BODY_STATIC` and `KINEMATIC` `DynamicsType`s, `ElasticJoint`, and
  `TimoshenkoJoint` (a 2-node corotational Timoshenko beam element between
  two `Body`s — the distributed-stiffness counterpart of `ElasticJoint`;
  chains form a beam, closed-form validated for bending with the `PL/kGA`
  shear term, axial and torsion, with nonlinear-stiffness dropoff).
- Beam-anchored points: a `BODY_STATIC` point rides a body rigidly, or a
  joint's corotational-Hermite deformed centerline (`joint_idx`,
  `beam_frac`, `beam_offset_b`, splitting its load to both nodes). A
  body-anchored point's body-frame offset is auto-derived from its
  `pos_cad`, so riding a body needs only `body=`/`joint=` plus a CAD
  position.
- Per-segment material `density` [kg/m³]: each `Segment` and `Tether` carries
  its own `density` (from the YAML `materials` table), replacing the single
  global `set.rho_tether` in mass calculations. Falls back to `set.rho_tether`
  when unset.
- Twist-surface restoring `stiffness` [N·m/rad] for the twist DOF; the
  number of unrefined VSM sections can be read from the YAML.
- `SymbolicAWEModels.ObjAdapter`: mesh mass properties (`center_of_mass`,
  `calculate_inertia_tensor`, `unit_inertia_from_obj`) computed in this
  package instead of VortexStepMethod, plus per-wing `mass`, `com` and
  `unit_inertia` (symmetric per-unit-mass 6-vector `[Ixx,Iyy,Izz,Ixy,Ixz,Iyz]`
  in [m²], scaled by `wing.mass`) as YAML columns and `VSMWing` kwargs. When
  omitted they auto-compute from an `.obj` if one is supplied, else fall back
  to point-mass inertia. `create_vsm_wing` generates the aero YAML from an
  `.obj` at panel resolution and caches it under `<model_dir>/obj_geometry`.
- `PrincipalFrameMethod` (`EIGEN_DECOMP`, `Y_ROTATION`) is exported, selecting
  how a wing's inertia tensor is diagonalized.
- `position_slots(sys_struct)` is exported: the index layout of a `SysState`'s
  `X/Y/Z` arrays (`points`, `panel_corners`, `wings`, `bodies`, `total`), so a
  log consumer no longer has to reconstruct the packing by hand.
- Replay camera pan and zoom controls (#250).
- Makie plotting: beam elements are drawn as instanced see-through tubes
  (`show_beams`), every layer can be toggled from the replay window, and an
  `aero_mapping=true` overlay draws the `AeroPressure` station→point scatter
  map (`m` toggles it live).
- `examples/inflated_beam_fit.jl`: fits a per-joint nonlinear bending law of a
  pressurised tube from the Comer-Levy wrinkled-section theory and validates the
  resulting `Body`/`ElasticJoint` chain as a cantilever.

### Changed
- BREAKING: winch motor dynamics are a type, not a builder function.
  `Winch.model` is now an `AbstractWinchModel` (`TorqueWinch`, the default, or
  `CascadedLengthWinch`) instead of a `Function`, and the builder hook is
  `winch_component(model, sys_struct, winch_idx; name, params)` rather than
  `default_winch_component`. *How to migrate:* wrap a custom builder in a new
  `AbstractWinchModel` subtype with a `winch_component` method; a custom model
  also wants `is_builtin_winch` left at its `false` default so `init!`
  rebuilds the cached model.
- BREAKING: `Wing` and `RigidBody` are one `Body{A, D}` type. A wing is a
  `Body` that carries aero, so `sys.wings` is a filtered view of `sys.bodies`.
  The exported concrete types are `RigidWing{A} = Body{A, RigidDynamics}`
  (6-DOF pose) and `ParticleWing{A} = Body{A, ParticleDynamics}` (pose fitted
  from structural points), replacing direct use of the old wing struct.
- The dependency stack was widened to the newer MTK/Symbolics generation
  (Symbolics 7→8, SymbolicUtils 4→5, SymbolicIndexingInterface 0.3→0.4,
  OrdinaryDiffEqBDF 1→2/3, OrdinaryDiffEqCore 3→4, NonlinearSolve 4→5/6,
  SteadyStateDiffEq 2.5→3/4; `ADTypes` added as a direct dependency).
  `VortexStepMethod` compat is raised to `4` and `KiteUtils` to `0.11.9`.
  `mtkcompile` now keeps body-frame outputs (`body_pos_w`, `body_vel_w`,
  `body_R_b_to_w`) and pulley `l0` in the torn-out state set. Not breaking:
  the package version is part of the `.bin` cache filename, so upgrading
  auto-invalidates stale model caches.
- RHS performance: model parameters are flattened for faster evaluation
  (#234). Rigid-body / beam element handling improved with aero separated
  into its own component (#236); `winch.speed_controlled` restored after the
  winch refactor dropped it.
- BREAKING: the Makie extension is triggered by `MakieControlPlots` instead of
  `Makie`, so `using GLMakie` alone no longer enables `plot`/`replay`/`record`.
  *How to migrate:* add `using MakieControlPlots` next to your Makie backend.

### Fixed
- Rigid-body pose drift: the principal quaternion ODE had no norm constraint, so
  `|Q_p|` drifted (~0.999) and `quaternion_to_rotation_matrix`, which assumes a
  unit quaternion, produced a non-orthonormal `R_b_to_w` (`det = |Q|^6`). Points
  were placed with that skewed `R` while the getters re-orthonormalised, giving a
  ~0.02 pose mismatch under twist. `rigid_body_eqs.jl` now carries a Baumgarte
  term `k(1 − ‖Q‖²)Q` that holds `|Q_p|` on the unit sphere.
- `AeroLinearized`'s refresh writes the operating-point force via
  `apply_direct_forces!` (as `AeroDirect` does), so `wing.aero_force_b` tracks
  the VSM solve without a second state sync.
- `quaternion_to_rotation_matrix` now normalizes internally, dividing each
  product by `sum(abs2, q)` so the result is orthonormal for any nonzero `q`
  (no `sqrt`, no measurable cost). The integrated attitude quaternion only
  meets its unit-norm constraint to solver tolerance — the Baumgarte term in
  `rigid_body_eqs.jl` bounds the drift but, with a time constant of
  `1/(2*quat_norm_gain) = 0.05 s`, cannot remove it within a step. At the
  default `abs_tol = rel_tol = 0.01` this left `R_b_to_w` off orthonormal by
  up to 3e-4 during fast transients, a scale error in every body-to-world
  transform feeding `moment_p`, `wing_pos`, `wing_vel` and `wing_acc`.
  Orthonormality error is now at machine precision (9e-16 worst case over a
  steering ramp, from 3.3e-4). Regenerates the model equations, so cached
  `model_*.bin` files from earlier commits must be rebuilt (`remake=true`).
- Wing body frame no longer inherits the principal frame (#245). A wing that
  declares no `origin`/`z_ref_points`/`y_ref_points` used to get
  `R_b_to_c = R_p_to_c`, so the eigendecomposition's ambiguous axis assignment
  for near-equal principal moments leaked into the *body* frame and flipped it
  ~90° relative to the VSM/CAD frame, causing growing lift/drag oscillation and
  eventual VSM non-convergence during reel-out. The fallback now keeps the CAD
  orientation (`R_b_to_c = I`, origin at the COM), leaving the principal choice
  a pure gauge — asserted by `test_principal_frame_invariance.jl`. For wings
  symmetric about the XZ-plane, `principal_frame_method=Y_ROTATION` selects the
  unique closed-form diagonalization instead of the permutation search.
- Transform placement (`apply_azimuth_elevation!`) is now roll-free and
  zenith/nadir-safe: the current radial is rotated onto the target with a
  minimal rotation (no dependence on the source frame, undefined at the
  zenith), with roll set solely by the heading step and a warning at the
  zenith/nadir where azimuth and heading are undefined. Placement of existing
  models shifts slightly.

### Removed
- BREAKING: the `WING`, `QUASI_STATIC` and `FIXED` `DynamicsType`s.
  `DynamicsType` is now `{DYNAMIC, STATIC, BODY_STATIC, KINEMATIC}`.
  *How to migrate:* a `WING` point becomes `DYNAMIC` on a particle wing, or
  `BODY_STATIC` on a rigid/beam wing (riding a body via `body=`/`wing=`, or a
  beam element via `joint=`); `QUASI_STATIC` becomes `DYNAMIC`; a `FIXED`
  twist surface becomes `STATIC` (`KINEMATIC` is the new geometry-prescribed
  twist source, e.g. a flap hinge). Wing↔point membership is no longer a
  point type: it is carried by `twist_surface.point_idxs`, exposed as
  `point.is_wing_node`, and queried with `wing_frame_member(point, wing_idx)`.
- BREAKING: `auto_create_twist_surfaces!`. A section-coupled (VSM) RIGID
  wing that declared no `twist_surfaces` used to have one `DYNAMIC`
  `TwistSurface` invented per LE/TE section; this silent black box is gone
  and the case now raises a clear error. *How to migrate:* declare explicit
  `twist_surfaces` covering the wing's LE/TE structural sections (see
  `data/2plate_kite/rigid_structural_geometry.yaml`). `AeroNone` wings are
  unaffected (they never coupled to sections).
- The redundant wing-side `point_idxs` list in YAML (non-breaking: the loader
  never read it — wing↔point membership is the point row's `wing_idx`
  column). Drop it from your wing definitions.
- BREAKING: `update_segment_forces!`. Segment forces are written by
  `update_sys_struct!` every step, so the separate refresh call had no
  remaining caller.
- `src/tether_properties.jl` and its `calc_spring_props` / `calc_tether_props`
  / `in_percent_band` helpers (unexported, unused — the one-segment equivalent
  spring fit went away with `copy_to_simple!` in v0.6.0).

## v0.12.0 12-06-2026

### Added
- Winch interface (#210): each `Winch` carries a `model` builder
  (`model(sys_struct, winch_idx; name) -> System`, default
  `default_winch_component`) so custom winch dynamics plug in as subsystems;
  contract checked by `validate_winch_component`. New `speed_controlled` flag
  prescribes reel-out speed directly. See `examples/custom_tape_winch.jl`.
- Swappable per-wing aero modes (#221 and follow-ups): each `Wing` carries an
  `aero::AbstractAeroModel` selecting its aerodynamics by dispatch. Built-ins
  `AeroLinearized` (default for `RIGID_DYNAMICS`), `AeroDirect` (default for
  `PARTICLE_DYNAMICS`), `AeroPlate`, and `AeroNone`, one file each under
  `src/aero_modes/`; chosen via the `aero` kwarg or the YAML `aero_mode`
  column. The mode's `aero_component(mode, sys_struct, wing_idx; name)`
  returns a subsystem wired at a fixed body-frame connector contract per
  `dynamics_type` (RIGID: `va`, `rho`, `R_b_w`, `omega`, `twist`, `twist_vel`
  → `force`, `moment`, `twist_moment`; PARTICLE: per-point `pos`/`vel`/`va`/
  `rho` → `point_force`), validated by `validate_aero_component`. The
  generated RHS stays allocation-free (`test_bench.jl`).
- A custom aero mode needs exactly two methods: `aero_component` and
  `aero_mode_tag` (cache tag). Everything else is an optional hook with a
  working default, dispatched on the mode — lifecycle (`setup_aero!`,
  `remake_aero!`, `validate_aero_structure`, `resize_aero_state!`,
  `init_aero_state!`), low-frequency refresh (`refresh_rigid_aero!`,
  `refresh_particle_aero!`, orchestrated by `refresh_aero!`), diagnostics
  (`calc_aoa`, `normalized_inertia`),
  log-point visualization (`n_aero_log_points`, `write_aero_log_points!`,
  `read_aero_log_points!`, `restore_aero_twist!`), and live Makie rendering
  (`plot_wing_aero!` / `update_wing_aero_plot!`, with methods in the Makie
  extension). There are no `isa`/`is_vsm` branches anywhere in the
  pipeline, so a custom mode is never excluded from a code path it cannot
  extend. VSM state (solver, geometry, linearization buffers) lives in a
  `VSMEngine` carried by `AbstractVSMAero` modes; subtyping it inherits the
  VSM implementation of every hook.
- `normalized_inertia` returns per-unit-mass inertia [m²] for every mode —
  the VSM `ObjWing` mesh tensor is already normalized and is now passed
  through as-is, the default normalizes the WING-point point-mass inertia
  (`normalized_point_inertia`), and the single scaling by `wing.mass`
  happens in `setup_wing_frame!`.
- `has_custom_component(sys_struct)`: `init!` defaults `remake` to rebuild
  automatically when a custom winch/aero component is present (their
  equations are not captured by the model hash). Structural mode fields enter
  the cache key via `aero_hash_id`; the mode tag enters the cache filename.
- Flat-plate wings log a display quad per section (4 corners, square of side
  `sqrt(area)`, structural point at quarter chord) via the log-point hooks,
  so plate geometry shows up in `SysState` logs like VSM panels do.
- New `FIXED` `DynamicsType`: a twist surface whose twist is a prescribed
  control input (no differential state). Flat-plate surfaces use it.

### Changed
- BREAKING: `Group` is renamed to `TwistSurface` throughout (type, YAML
  section, and fields, e.g. `wing.group_idxs` → `wing.twist_surface_idxs`).
  Flat-plate surfaces are now 1-point `FIXED` `TwistSurface`s instead of a
  separate plate type.
- BREAKING: the wing types are merged into one `Wing` struct. `VSMWing` and
  `PlateWing` remain as constructor functions; the polar lookups and drag
  correction of a flat-plate wing live on its `AeroPlate` mode. The
  `wing_type` keyword is deprecated in favour of `dynamics_type`.
- `AeroNone` carries no VSM engine and needs no VSM geometry or `vsm_set`, so
  a pure rigid-body wing builds without any VSM setup (`VSMWing` accepts
  `vsm_set=nothing` for engine-less modes).
- The symbolic aero generation was restructured: `vsm_eqs.jl`, `plate_eqs.jl`
  and `linearize.jl` are replaced by a thin mode-agnostic wiring layer
  (`aero_eqs.jl`), the per-mode files under `src/aero_modes/`, and
  `twist_surface_eqs.jl` (formerly `group_eqs.jl`). `SystemStructure`
  construction is split into `setup_wing_frame!` (mass/inertia/body frame,
  aero-independent) and the mode-dispatched `setup_aero!`.
- The Makie extension is aero-mode agnostic: wing rendering dispatches on
  the aero mode via `plot_wing_aero!(ax, sys, wing, mode)` (and the per-frame
  `update_wing_aero_plot!`), with methods living in the extension so a custom
  mode draws with full Makie access. Flat-plate wings now render their
  section quads in `plot` and `replay` in the VSM panel style (red mesh,
  black borders).
- The transform pipeline (`apply_heading!`, `finalize_transforms!`) no longer
  filters wings by aero mode; flat-plate wings now get heading and frame
  finalization like every other wing.
- Internal renames for readability: leading-underscore function names and
  short abbreviations were removed throughout.

### Removed
- BREAKING: the exported `PlateSurface` type and the `AeroMode` enum with its
  `AERO_NONE`/`AERO_DIRECT`/`AERO_LINEARIZED`/`AERO_PLATE` values, along with
  the `BaseWing` type. Use the `AbstractAeroModel` mode structs and the single
  `Wing` type instead.
- BREAKING: `VSMWing` and `PlateWing` are now constructor functions, not types,
  so `wing isa VSMWing` / `isa PlateWing` errors. Use the exported
  `wing.aero isa AbstractVSMAero` / `wing.aero isa AeroPlate` if you need the
  check, or better, dispatch on the aero mode.
- The dead `SystemStructure` fields `y`, `x`, `jac` (legacy linearization
  buffers; the per-wing state lives in each mode's `VSMEngine`).
- The `exposes_aero_input` trait: the `aero_input` connector is detected by
  name on the built subsystem instead.
- The V3-kite-specific analysis code in the Makie extension:
  `compute_ekf_yaw_and_rate`, `compute_ekf_yaw_and_rate_tension`,
  `calculate_cs`, `calc_ref_area`, `middle_le_to_kcu_dir` and their helpers,
  along with the `plot_cs`, `plot_yaw_rate_paper` and `plot_gk_paper` panels
  (the last hardcoded a V3 segment index). V3Kite carries its own copies.
- The `tape_lengths` kwarg of the multi-panel plot and the hardcoded
  steering reconstruction from `segments[87]`: the `plot_us` and `plot_gk`
  panels now read the logged `syslog.steering` directly (so `steering` must
  be written into the `SysState` before `log!` for these panels to show
  data).
- `set_depower_steering!`, `min_chord_len`, and the
  `SymbolicAWEModel.set_tether_len` field (3-line-kite-specific set-point
  logic with hardcoded tether indices). `calc_side_slip` no longer
  dispatches on the aero mode — it is the same apparent-wind formula for
  every mode and takes just the wing.

### Fixed
- A `DYNAMIC` twist surface without aero sections left
  `twist_surface_aero_moment` unbound and broke `mtkcompile`; the wiring now
  binds the aero component's `twist_moment` for every non-`FIXED`-empty
  surface.
- Makie extension: `wing isa VSMWing` checks in the panel plotting and the
  log-slot lookup threw at runtime (`VSMWing` is a constructor function since
  the wing merge, not a type).

### Compatibility notes
- Plate logs recorded before the quad logging have a different point count
  and will not replay.

## v0.11.1 06-06-2026

### Added
- `init_stretch_frac` (YAML column and `Tether(...; stretch_frac)` kwarg),
  mutually exclusive with `init_tether_force`: `reinit!` derives the
  unstretched `len` as `len = stretch_frac·stretched`. Setting one input
  clears the other. `init_stretch_frac` must be positive: `<1` pre-stretch,
  `1` neutral, `>1` slack.
- `test_twist_alignment.jl`: under group twist the structural strut
  trailing-edge points stay aligned with the deformed VSM panel trailing
  edges for a `RIGID_DYNAMICS` wing.

### Changed
- `VortexStepMethod` compat raised to `3.3.5`.

### Fixed
- Per-group unrefined moment uses the VSM solver field
  `moment_coeff_unrefined_dist`.
- Body-frame camera tracking across animation frames (Makie ext):
  `update_cam!` with explicit up-vector and `PLOT_BODY_PREV_WING_POS`
  to eliminate view drift.

## v0.11.0 02-06-2026

### Breaking
- Tether `init_unstretched_length` (YAML) removed; specifying it errors.
  The unstretched rest length is now *derived*: placement is driven by
  `init_stretched_length` (the standoff / placed point geometry,
  default = geometric) and `init_tether_force` (default 0), and
  `len = stretched·(1 − force/unit_stiffness)`.
- `Tether.init_stretched_len`/`init_unstretched_len` are now
  `Union{SimFloat,Nothing}` (`init_unstretched_len` is derived); `Tether`
  gained `init_tether_force`; the positional length constructor arg
  (now the stretched length) is optional. Serialized models must be
  rebuilt.
- `VSMWing` `origin_idx`/`origin_ref` replaced by
  `origin::WeightedRefPoints` (weighted body-frame origin).
- `update_yaml_from_sys_struct!` and `update_sys_struct_from_yaml!`
  removed (unreliable line-based YAML round-tripping, no longer used).

### Added
- `init_tether_force` (YAML / `Tether(...; tether_force)`, default 0):
  `reinit!` derives every tether's unstretched `len` from the placed
  stretched length, `len = stretched·(1 − force/unit_stiffness)`;
  force 0 gives zero tension.
- `init!`/`reinit!` `apply_tether_lengths` kwarg to skip placement.
- `WeightedRefPoints(::AbstractString)`; `yaml_parse_origin` for
  weighted origin specs.
- Helpers: `apply_tether_init_forces!`, `tether_unit_stiffness`,
  `tether_anchor_free`, `rigid_point_siblings`, `parse_tether_init`,
  `tether_ordered_point_idxs`, `tether_downstream_idxs`,
  `group_tethers_by_overlap`, `apply_cluster_init_stretched_len!`,
  `_wing_log_pos`; `test_tether_init.jl`.

### Changed
- Tether placement honored only on *root* tethers (one endpoint on a
  `STATIC`/winch boundary — the fixed anchor, either end); a tether
  with neither endpoint anchored is an error. Tethers sharing a
  `RIGID_DYNAMICS` wing are treated as one cluster (rigid-body
  connectivity). Multi-root clusters placed by the mean displacement of
  all roots (length + direction), logging `@info` (gated on `prn`).
- Wing position stored in dedicated `SysState` slots; reads via
  `update_from_sysstate!` / `_wing_log_pos` / Makie body-frame arrows
  use `wing.pos_w` directly.
- `build_point_to_vsm_point_mapping` takes a `VSMWing`, using
  body-frame closest-point distances.

### Fixed
- Makie zoom/pan world-camera save/restore (no view drift); body-frame
  zoom distance preserved across mode switches.
- `vsm_refine.jl`: RIGID_DYNAMICS wings always keep their aerodynamic
  panel geometry (mesh- or YAML-defined); section rebuilding from
  structural points is now PARTICLE_DYNAMICS-only. The 2plate aero
  geometry was corrected to match its structural points.
- `get_sys_struct_hash` hashes `wing.origin`.

## v0.10.0 30-05-2026

### Changed
- BREAKING: `WingType` constants `QUATERNION` and `REFINE` are now
  deprecated. Use `RIGID_DYNAMICS` and `PARTICLE_DYNAMICS` instead.
  Deprecated bindings emit a warning and will be removed in a future
  release.
- `DataInterpolations` added as a package dependency (required for
  `PlateWing` polar interpolation).
- `bin/install` now displays an interactive menu to choose Julia
  version (1.11 or 1.12) when no version parameter is provided. The
  currently active Julia version is highlighted as the default. Menu
  is skipped if a version is specified via `--version` or `+X.Y`
  parameters.

### Added
- `PlateWing` and `PlateSurface` types for flat-plate CL/CD lookup
  aerodynamics.
- `AERO_PLATE` aerodynamics mode — evaluates lift and drag from a
  polar table (CL/CD vs α) via registered symbolic interpolants.
- `create_plate_interpolations(alpha_deg, cl_data, cd_data)` — helper
  to build CL and CD interpolation objects (cubic or linear spline)
  for use with `PlateWing`.
- `examples/kps4_comparison.jl` — comparison of a `PlateWing`-based
  rigid-body kite model against the KiteModels kps4 reference.
- `data/kps4/` — YAML settings and system definition for the kps4
  plate model.
- Added missing examples to `examples/menu.jl`:
  `coupled_linearize`, `cosine_steering_trajectory`,
  `kps4_comparison`, `vsm_linearization`, and `sam_tutorial`.

### Fixed
- `init_stretched_len` now works for multi-tether systems. Tethers
  sharing downstream structure are placed to a single effective
  length (the average of several specified values, with a warning),
  and the initial-positioning BFS no longer drags other tethers'
  ground anchors — it stops at `STATIC` points and winch points
  (which may be `DYNAMIC`).
- `bin/create_sys_image`: fixed a bug that prevented deletion of
  stale `.so` files before rebuilding the system image.
- `AUTHORS.md`: corrected contributor entry.
- `examples/kps4_comparison.jl`: fixed soft-scope ambiguity warning
  for `sys_state` inside the simulation loop.
- Multi-log `plot()` legend labels now render correctly as LaTeX.
  Added `lbl()` helper that places the symbol and suffix inside a
  single `$...\text{...}$` math environment, fixing the literal
  `$\gamma$ (SymAWE)` display in the legend.

### Removed
- `examples/makie_polar_plots.jl` — removed (functionality
  superseded).

## v0.9.0 20-05-2026

### Changed
- BREAKING: simplified `AERO_LINEARIZED`. ForwardDiff Jacobian
  over `[α, β, ω, θ_groups]` returning wind-axis coefficients
  `[CL, CD, CS, CM, cm_groups]`. Wing fields and accessors
  renamed `vsm_*` → `aero_*`.
- A RIGID_DYNAMICS wing can now have fewer groups than unrefined
  aero sections (one twist DOF drives several sections via a
  spatial partition). More groups than sections errors.
- Bumped `VortexStepMethod` compat to `3.3.0`.
- License changed from MIT/MPL-2.0 to LGPL-3.0-only. All source
  files updated with REUSE-compliant SPDX headers.
- `bin/install` rewritten: unified menu, optional precompile skip,
  removed `bin/update_manifest` and `bin/create_sys_image2`.
- `bin/create_sys_image` updated with improved comments and options.
- `bin/reuse_lint` made more robust with fallbacks for missing tools.
- Safe `atan`/`smooth_normalize` replacements for `asin`/`normalize`
  in VSM equations and linearisation to avoid NaN at edge cases.

### Added
- `examples/vsm_linearization.jl` — plots the VSM linearisation
  tangents around the operating point.
- `test/util.jl` — shared test utilities for allocation checks across
  all integrators.

## v0.8.3 03-05-2026

### Changed
- VSM solver type is taken from VSM settings instead of being
  hard-coded to `NONLIN`.
- At low apparent wind, aero outputs are zeroed instead of warning
  and skipping. Threshold via new `vsm_min_wind` kwarg (default 0.5)
  on `init!`, `reinit!`, `next_step!`.
- Bumped `VortexStepMethod` compat to `3.2.0`.

## v0.8.2 26-04-2026

### Changed
- Updated the default manifest files.

### Added
- `drag_force` field on `Point` — total drag in world frame (point's
  own aerodynamic drag plus its share of connected segment drag).
  Populated by `update_sys_struct!` each timestep.
- Manifest freshness tests in `test_helpers.jl`: verify that no bare
  `Manifest.toml` exists and that `.default` manifests are at least as
  recent as `Project.toml`.
- CI step to copy version-specific `.default` manifest before build,
  ensuring the correct manifest is used per Julia version.
- Drag-related tests in `test_point.jl`, `test_segment.jl`, and
  `test_wing.jl`.

### Fixed
- Crash with Julia 1.11; `setup_env` updated to fix that.

### Removed
- `plot_recipe.jl` — unused legacy Plots.jl recipe. Visualization is
  handled by `SymbolicAWEModelsMakieExt`.

## v0.8.1 23-04-2026

### Changed
- `SystemStructure.set` field is no longer `const`, allowing change
  after deserialisation.
- Replaced all `@unpack` macro usage with Julia's native destructuring
  syntax `(; a, b) = x`.

### Fixed
- Fixed JETLS warnings across multiple source files.
- `bin/install` now copies `.JETLSConfig.toml.default` to
  `.JETLSConfig.toml` if it does not exist, and warns when an existing
  config differs from the default.
- `bin/install` warning messages now use colored output for visibility.

## v0.8.0 18-04-2026

### Changed
- BREAKING: `SegmentType` positional argument removed from `Segment`
  constructor. Use `unit_stiffness`, `unit_damping`, `diameter_mm`
  kwargs or a YAML material instead. The `SegmentType` enum is kept
  temporarily to produce a helpful deprecation error.
- BREAKING: `winch_point` moved from `Tether` to `Winch`. Pass
  `winch_point` as a keyword to the `Winch` constructor instead.
- BREAKING: Heading calculation changed from wind-perpendicular
  projection to tangential sphere frame. `calc_heading(R_b_to_w,
  wind_norm)` → `calc_heading(R_b_to_w, wing_pos)`.
  `get_heading_components()` removed. `solve_heading_rotation` takes
  `wing_pos` instead of `k, wind_norm`.
- BREAKING: `Tether` struct fields restructured — `winch_point_idx/ref`
  removed, new fields: `start_point_idx/ref`, `end_point_idx/ref`,
  `n_segments`, `unit_stiffness`, `unit_damping`, `diameter`.
- BREAKING: `create_tether()` utility returns a 5-tuple (added
  `ground_point_idx`) and no longer takes a `SegmentType` argument.
- BREAKING: YAML segment format no longer has a `type` column. Existing
  YAML files with a `type` column in segments will raise an error.
- BREAKING: `tether_len` moved from `Winch` to `Tether`. Each tether
  now owns its length as an ODE state variable. Winch-connected tethers
  evolve via `D(tether_len) = winch_vel`; winch-less tethers have
  constant length (`D(tether_len) = 0`).
- BREAKING: `tether_vel` renamed to `winch_vel` and remains on `Winch`.
  `tether_acc` renamed to `winch_acc` in the generated equations.
- BREAKING: `SimpleLinModelWithAttributes` removed. The
  `simple_lin_model` field is no longer part of `SymbolicAWEModel`.
  `simple_linearize!` is no longer exported.
- BREAKING: `sim_oscillate!` and `sim_turn!` removed. Use `sim!` with
  a custom `set_values` matrix instead.
- BREAKING: `update_aero_yaml_from_struc_yaml!` no longer exported.
- BREAKING: `set` field removed from `SymbolicAWEModel`. Settings are
  now read from `sam.sys_struct.set`. The `set_set` setter was removed
  from `ProbWithAttributes` and `LinProbWithAttributes`.
- BREAKING: `get_struct_state` removed from `ProbWithAttributes`.
- Wind equations now use `get_wind_vec` internally instead of
  separate `get_v_wind`, `get_upwind_dir`, and `get_wind_elevation`
  accessors. Not breaking: KiteUtils `Settings` syncs `wind_vec`
  from `v_wind`/`upwind_dir`/`upwind_elevation` automatically when
  `use_wind_vec=false` (the default).
- Tethers no longer require a connected winch. Winch-less tethers use
  constant `l0` from segment properties.
- `compression_frac` description clarified: "Compressive/tensile
  stiffness ratio (0-1). 0 = no compression stiffness."
- `init!`, `next_step!`, `update_sys_state!` are no longer exported
  and must be imported from `KiteUtils`.
- `sim!` now requires `y_op` keyword argument when `lin_model` is
  provided (previously obtained from the removed simple lin model).
- `SerializedModel` type parameters tightened for `defaults` and
  `guesses` fields.
- fixed most `JETLS` warnings for improved robustness and performance.
- Package version is now included in `.bin` cache filenames, so
  upgrading the package automatically invalidates stale cached models.
- the script `bin/run_julia` was updated to work also with Julia 1.12.6

### Added
- Route 2 tether auto-generation: `Tether(name; start_point,
  end_point, n_segments)` automatically creates intermediate points
  and segments, evenly spaced between endpoints. YAML format:
  `headers: [name, start_point, end_point, n_segments, ...]`.
- Route 1 tethers auto-detect `start_point_idx` and `end_point_idx`
  from the first/last segment endpoints.
- Comprehensive docstrings on all `Point`, `Group`, `Segment`,
  `Pulley`, `Tether`, `Winch`, and `Transform` struct fields.
- `WeightedRefPoints` exported for weighted reference point support.
- `init!` keyword `reinit_sys` to optionally skip system structure
  reinitialization.
- New tests: "Route 2 auto-generated tether" and "Tether without
  winch" in `test_tether_winch.jl`.
- New test file `test_tether_init.jl` for tether initialization.
- New test file `test_yaml_weighted_ref.jl` for weighted reference
  point YAML loading.
- Airbag pressurized membrane simulation example (`examples/airbag.jl`).
- the script `bin/install`. Use it after installation from git.
- the script `bin/create_sys_image`. Improves time for first run
  by a factor of 3-5.
- the scripts `bin/install_jetls` and `bin/jetls` to install and run
  `JETLS.jl`, a static code checker for Julia.
- Developer documentation improvements (troubleshooting section for
  segfault issues, updated docs to use GLMakie).

### Fixed
- YAML `calculate_derived_properties!` no longer requires `l0` to
  compute `unit_stiffness` from material properties (needed for
  Route 2 tethers).
- YAML `update_yaml_from_sys_struct!` regex updated for the new
  segment format (no `type` column).
- YAML weighted reference point loading fixed (broken deserialization
  of weighted refs).
- Heading calculation uses tangential sphere frame, fixing drift issues
  with the old wind-perpendicular projection.
- Unknown solver string (e.g. `DFBDF` from default KiteUtils settings)
  no longer throws an error — a warning is emitted and the solver
  falls back to `FBDF`.
- README code examples now include the required
  `SymbolicAWEModels.init_module(; force=false)` call so they work
  correctly on a fresh install.
- README pendulum example also calls `set_data_path("data/base")`
  before loading `Settings`.

### Removed
- `SimpleLinModelWithAttributes` struct and `simple_linearize!`.
- `sim_oscillate!` and `sim_turn!` simulation functions.
- `getstate` and `setstate!` functions from `linearize.jl`.
- `upwind_dir` helper function (replaced by `wind_vec`).
- Branch-specific system images: `bin/create_sys_image` and
  `bin/run_julia` no longer embed the git branch name in the `.so`
  filename. A single `kps-image-<julia_major>.so` is used instead.

### Tests
- README pendulum example and README 2-plate kite example are now
  executed in `test/setup_integration.jl`.

## v0.7.2 18-03-2026

### Added
- `speed_controlled` field on `Winch` — when `true`, tether velocity
  is prescribed externally (`D(tether_vel) = 0`) while length still
  tracks velocity.
- Multi-system `record()` for recording side-by-side SysLog animations
  to video (MP4/GIF/MKV/WebM).
- Makie extension test suite (`test_makie_extension.jl`) covering
  multi-system plot, record, and replay.
- Zenodo metadata (`.zenodo.json`) and `CITATION.cff` for citing the
  package.
- CI: GLMakie tests on Linux via `xvfb-run`, Julia 1.12 test matrix.

### Fixed
- `reposition!()` now uses the analytical `solve_heading_rotation`
  for wind-relative heading, consistent with `reinit!`. Previously
  heading was applied as a relative delta, causing drift.
- `reposition!()` correctly updates PARTICLE_DYNAMICS wings by recalculating
  `R_b_to_w` and `pos_b` from structural points.
- Multi-system `plot()` now passes vector-typed segment colors,
  fixing a crash when `setup_segment_hover_events!` assigned
  `Vector{RGBA}`.
- `init!()` validates that `SystemStructure` uses `VSMWing` type
  before equation generation.
- `sim_reposition!()` passes absolute heading to the transform
  instead of subtracting the current wing heading.
- Typo fixes in README and documentation ("ODE solver" → "ODE
  problem").

### Changed
- `sam_tutorial.jl` example updated: adds WING-type points and uses
  `VSMSettings` with `data_prefix=false`.
- Examples updated to pass `data_prefix=false` to `VSMSettings`.
- 2plate_kite aero geometry TE z-coordinates adjusted.
- `settings.yaml` now includes `sample_freq` field.

## v0.7.1 27-02-2026

### Added
- `update_sys_struct_from_yaml!()` — update a `SystemStructure` in-place
  from a modified YAML file (point `pos_cad` and segment `l0`).
- `segment_cad_length()` and `autocalc_tether_len()` shared helpers,
  replacing duplicated code in the constructor, `reinit!`, and YAML loader.

### Fixed
- `SystemStructure` constructor auto-calculates `winch.tether_len` from
  all connected tethers (was only using the first).

## v0.7.0 DD-02-2026

### Changed
- BREAKING: Julia version requirement raised from 1.10 to 1.11, 1.12.
- `reinit!()` uses a unified code path for all wing types, calling
  `match_aero_sections_to_structure!` and
  `compute_spatial_group_mapping!` during VSM rebuild.
- `test_bench.jl` refactored from ad-hoc benchmarks into a proper
  `@testset` suite with `setup_bench_sam()` helper.
- Added `[workspace]` configuration in `Project.toml` for docs, examples,
  scripts, and test sub-projects.
- Manifest files renamed to `.default` suffix and gitignored.

### Added
- Asymmetric aero/structural section counts: aerodynamic and structural
  meshes can now have different numbers of sections. When counts differ,
  `match_aero_sections_to_structure!()` rebuilds unrefined
  sections from structural LE/TE positions while `use_prior_polar=true`
  preserves existing refined panel polars. Opt-in via
  `use_prior_polar=true` on the VortexStepMethod wing.
- `identify_wing_segments()` — identifies LE/TE pairs from groups
  (preferred) or via a consecutive-pair heuristic.
- `compute_spatial_group_mapping!()` — maps groups to VSM sections by
  spatial proximity, supporting n_groups != n_aero_sections.
- PARTICLE_DYNAMICS wings can now have groups (used for LE/TE pair identification).
- RIGID_DYNAMICS wings can now have `wing_segments` for structural geometry
  locking.
- YAML loader fallback LE/TE detection in
  `update_aero_yaml_from_struc_yaml!()` when no groups are defined
  (consecutive-pair heuristic with x-coordinate check).
- `test_match_aero_sections.jl` — tests geometry matching and polar
  interpolation for both PARTICLE_DYNAMICS and RIGID_DYNAMICS wings, including
  mismatched section counts.
- Helper scripts: `bin/install` (environment setup, Julia version detection)
  and `bin/run_julia` (launcher with system image support).

## v0.6.1 23-02-2026

### Fixed
- Disable VSM auto-sorting of sections (`sort_sections=false`) in all
  VortexStepMethod calls. Auto-sorting silently broke the correspondence
  between VSM sections and structural point indices / group mappings.

## v0.6.0 21-02-2026

### Changed
- Component constructors (`Point`, `Segment`, `Wing`, `Winch`,
  `Transform`) now accept a symbolic `name` (Symbol) as the first
  argument in addition to numeric indices. Numeric `idx` values still
  work. Use e.g. `Point(:kcu, pos, DYNAMIC)`.
- BREAKING: `Segment` constructor takes separate `point_i`, `point_j`
  arguments instead of a `point_idxs` vector.
- BREAKING: Rotation matrix fields renamed from `R_a_b` to `R_a_to_b`
  throughout (e.g. `wing.R_b_w` → `wing.R_b_to_w`).
- BREAKING: `ControlPlotsExt` package extension removed. Visualization is
  now handled entirely by `SymbolicAWEModelsMakieExt`.
- BREAKING: Predefined model factory functions removed
  (`create_ram_sys_struct`, `create_simple_ram_sys_struct`). Build models
  using component constructors or YAML instead.
- BREAKING: Ram air kite and V3 kite models moved to dedicated packages
  ([RamAirKite.jl](https://github.com/OpenSourceAWE/RamAirKite.jl),
  [V3Kite.jl](https://github.com/OpenSourceAWE/V3Kite.jl)).
  Their data directories are removed from this package.
- `src/system_structure.jl` split into modular files under
  `src/system_structure/` (types, core, utilities, transforms, wing,
  named_collection).
- `src/generate_system.jl` split into 13 focused modules under
  `src/generate_system/` (point_eqs, segment_eqs, wing_eqs, group_eqs,
  winch_eqs, pulley_eqs, tether_eqs, scalar_eqs, vsm_eqs, accessors,
  helpers, create_sys).
- Makie extension significantly overhauled with new plotting functions.
- Test suite completely rewritten. The old tests (`test_simulation`,
  `test_linearization`, `test_initialization`, `test_sam`, etc.) tested
  the full assembled kite model as a black box, making failures hard to
  diagnose. The new tests isolate each component with minimal models
  built from constructors, verifying physics against analytical
  solutions:
  - `test_point` — gravity free-fall, damping, quasi-static equilibrium
  - `test_segment` — spring-damper forces, stiffness, drag
  - `test_wing` — RIGID_DYNAMICS and PARTICLE_DYNAMICS wing construction, VSM coupling
  - `test_wing_dynamics` — rigid body torque response, precession,
    angular momentum conservation
  - `test_tether_winch` — reel-out dynamics, Coulomb and viscous
    friction, terminal velocity
  - `test_pulley` — equal-tension constraints, multi-segment pulleys
  - `test_transform` — spherical coordinate positioning
  - `test_quaternion_conversions` — quaternion ↔ rotation matrix
  - `test_quaternion_auto_groups` — auto-generated twist DOFs
  - `test_principal_body_frame` — principal vs body frame separation
  - `test_heading_calculation` — kite heading from tether geometry
  - `test_section_alignment` — VSM section ↔ structural point mapping
  - `test_profile_law` — atmospheric wind profile verification
  - `test_bench` — performance regression tracking
- Complete documentation overhaul with new pages: coordinate_frames,
  vsm_coupling, pipeline, tutorial_julia, tutorial_yaml.
- Data files reorganised: base settings moved to `data/base/`, new
  `data/2plate_kite/` and `data/saddle_form/` model directories added.

### Added
- `NamedCollection` indexing — components support symbolic names
  (e.g. `sys.points[:kcu]`, `sys.segments[:bridle_1]`).
  `SystemStructure` resolves all symbolic references to numeric indices
  automatically via `assign_indices_and_resolve!()`.
- `WingType` enum (`RIGID_DYNAMICS`, `PARTICLE_DYNAMICS`) for explicit wing type
  selection. BREAKING: these names replace the previous `QUATERNION` and
  `REFINE` wing types. Update YAML configs from `type: QUATERNION` /
  `type: REFINE` to `dynamics_type: RIGID_DYNAMICS` / `dynamics_type: PARTICLE_DYNAMICS`,
  and rename the wing `type` field to `dynamics_type`.
  Update any code using the old exported constants. `PARTICLE_DYNAMICS`
  applies per-panel forces directly to structural points for higher
  fidelity aeroelastic coupling.
- `AeroMode` enum (`AERO_NONE`, `AERO_DIRECT`, `AERO_LINEARIZED`) for
  build-time control over aerodynamic computation strategy.
- YAML-based model definition via `load_sys_struct_from_yaml()`,
  `update_yaml_from_sys_struct!()`, and
  `update_aero_yaml_from_struc_yaml!()`.
- PARTICLE_DYNAMICS wing support (`src/vsm_refine.jl`) — structural deformation
  coupled directly to VSM panel geometry with moment-preserving force
  distribution.
- Principal vs body frame separation for RIGID_DYNAMICS wings. Principal
  frame (diagonal inertia) used for Euler equations, body frame (from
  reference points) used for output and VSM coupling.
- Auto-group generation for RIGID_DYNAMICS wings when groups are not
  explicitly provided.
- `record()` for saving simulation replays to MP4.
- `plot_sphere_trajectory`, `plot_body_frame`, `plot_aoa` plotting
  functions.
- `update_segment_forces!`, `set_world_frame_damping`,
  `set_body_frame_damping`, `segment_stretch_stats` utility functions.
- New examples: `hanging_mass`, `catenary_line`, `saddle_form`,
  `coupled_2plate_kite`, `coupled_realtime_visualization`,
  `coupled_linearize`, `coupled_simple_lin_model`,
  `coupled_tether_deflection`, `heading_gate`,
  `cosine_steering_trajectory`, `makie_polar_plots`,
  `static_load_2plate_kite`.
- Benchmark test (`test_bench.jl`) for performance tracking.

### Removed
- `predefined_structures.jl` and factory functions
  (`create_ram_sys_struct`, `create_simple_ram_sys_struct`,
  `create_tether_sys_struct`, `copy_to_simple!`).
- Ram air kite data files, LEI kite directory, `data/kite.obj`.
- Old examples: `ram_air_kite`, `lin_ram_model`, `simple_lin_model`,
  `lin_simple_tuned_model`, `simple_tuned_model`,
  `realtime_visualization`, `reposition`, `tether_props`.
- `SymbolicAWEModelsControlPlotsExt` package extension.
- `src/precompile.jl`.

## v0.5.0 25-08-2024
### Removed
- BREAKING: the Winch struct doesn't have a model field anymore. Instead, all equations are symbolic, and the WinchModels dependency is removed.
### Added
- The function `calc_steady_torque` calculates the torque that will result in zero acceleration.

## v0.4.2 24-08-2024
### Fixed
- Don't write protect manifest

## v0.4.1 13-08-2025
### Fixed
- Update Artifacts.toml.default

## v0.4.0 13-08-2025
### Added
- Structs with attributes for better serialization and code structure (`SimpleLinModelWithAttributes`, `ProbWithAttributes`, `LinProbWithAttributes`, `ControlFuncWithAttributes`).
- `plot_force` option to the plot recipe.
- `model_management.jl` file to better organize the code.
### Changed
- BREAKING: `init_module` function to simplify project setup, replacing `install_examples`, `copy_examples`, `copy_bin` and `copy_model_settings`.
- Major refactoring of the `SymbolicAWEModel` and its initialization process. The `SerializedModel` struct is now much simpler and more robust.
- The `run_julia` script is now much more powerful, with argument parsing for `--copy-manifest` and `--precompile`.
- The precompilation process now uses artifacts instead of downloading files directly.
### Fixed
- URLs in `Artifacts.toml.default`.
- Cross-correlation analysis in tests.
### Removed
- `data/kite.obj` file.
- `copy_examples`, `copy_bin`, `copy_model_settings`, `install_examples` functions.

## v0.3.3 07-08-2025
### Fixed
- Fix non-persistent state bug with `calc_tether_props`

## v0.3.2 07-08-2025
### Fixed
- Fix documentation for sim_oscillate!

## v0.3.1 06-08-2025
### Fixed
- Fix examples and menu

## v0.3.0 06-08-2025
### Changed
- Breaking: sim!, sim_oscillate! and sim_turn! return a tuple (sl, lin_sl) instead of just a sl
### Fixed
- Restrict LinearSolve version to `<3.25.0`
- Fixed `linearize!(sam)` to get updated when the state gets updated
### Added
- Added `lin_simple_tuned_model.jl` example

## v0.2.1 01-08-2025
### Fixed
- Import Pkg

## v0.2.0 01-08-2025
### Added
- Adds simple model and tether model
- Adds `copy_to_simple!` function, which copies the ram model state to the simple model state, uses the tether model to find the equivalent 1-segment spring properties of the tether
- Adds open-loop sim functions `sim!`, `sim_oscillate!`, `sim_turn!`
- Adds plotting function `plot(sys_struct::SystemStructure, sys_log::SysLog)`
- Adds documentation
- Adds new updated tests: test/test_sam.jl
### Fixed
- Fixes documentation
- Fixes the bug where the kite could not have negative position
### Changed
- Improved precompilation
- Breaking: `Segment` constructor has different arguments
### Removed
- Removed `.bin` files from git, will be added as release artifacts

## v0.1.3 18-07-2025
### Changed
- Add interface keyword arguments to `init!`

## v0.1.2 13-07-2025
### Changed
- Update VortexStepMethod.jl

## v0.1.1 13-07-2025
### Added
- Added a simple linearized model
### Changed
- Improved the reinitialization using scalar settings values
- Update KiteUtils and AtmosphericModels

## v0.1.0
- Moved the SymbolicAWEModel from KiteModels.jl to SymbolicAWEModels.jl
