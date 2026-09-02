# CHANGELOG

## Unreleased

### Changed

- BREAKING: a wing's `twist_surface` is now a `station`, in the code and in the
  structural geometry YAML: the `twist_surfaces:` key and the per-wing
  `twist_surfaces: [...]` list are both `stations:`. The entity is the
  structural points lying in one chordwise plane, and twist is one of the things
  it carries, alongside the spanwise stations the aero load and the live polars
  are read on. It had grown three names — `twist_surface`, `station` and `strut`
  — which is how the load path and the deformation path came to compute the same
  blend twice without anyone noticing they were the same quantity.
  `TwistSurface` is `Station`, `station_control_points` says what it returns,
  and the file is `station_eqs.jl`.

### Fixed

- A panel's drag direction no longer re-derives its lift direction inside itself.
  `drag_cross` inlined `lift_cross / |lift_cross|` where the `dir_lift` variable
  already holds exactly that, so normalising `dir_drag` squared and summed three
  copies of it: one panel's drag direction reached 199k expression nodes and a
  wing's system 3.4M. Reading the variable halves the system the monolith backend
  compiles, and leaves the value unchanged.

- Building a wing model no longer constructs the per-section flow curvature rate
  when the wing's solver has the term switched off, which is the default. The
  symbolic inflow reconstruction built it for every refined section and dropped
  it one line later, while the kernel backend already gated the same gather.

- A transform's `elevation_vel` and `azimuth_vel` now place the velocity they
  describe. Both terms were wrong: the elevation term was never rotated by the
  azimuth, so off the `azimuth = 0` meridian it pointed out of the sphere, and
  the azimuth term had the sign the position convention does not have and was
  missing its `cos(elevation)`, which at 70° elevation is a factor of three. The
  field was also a per-component approximation — every component took the wing's
  own spherical velocity scaled by its radius, so the structure was not moving
  rigidly — and it was written before the heading step rotated the positions it
  was measured against, so a placed heading left position and velocity
  disagreeing by metres per second. It is now one rigid rotation about the base,
  `spherical_spin(transform)` crossed into `pos_w - base_pos`, applied after the
  heading; a rigid body takes the matching `ω_b` where it used to be zeroed.
- `reposition!` takes `update_vel` (`false` as before), so a state can be moved
  onto a pose and released already flying on it. Nothing applied the transform
  velocity except `reinit!`, which resets from CAD.

- A panel's aerodynamic load is shared between the two stations it lies between
  rather than rounded onto the nearer one. Rounding is a discontinuity, and two
  panels mirroring each other across the wing land on spanwise places that agree
  to within a few ulps but not exactly, so the nearer station could be a
  different one for each and a whole panel's load moved a station across the
  span. A symmetric wing at a symmetric pose carried yaw and roll that nothing
  cancelled; both are now identically zero. A panel straddling a station also
  took the mean of two weights measured against different station pairs, which
  means nothing averaged.

- `AeroPressure` now places each panel's pitching moment, not only its force.
  The scatter anchored the force to the VSM panel force and left the moment to
  whatever the frozen `Cp` pattern happened to integrate to, so the polar's `Cm`
  never reached the structure; the couple it did place was the flow-curvature
  increment alone. `pressure_couple` now places the deficit between the moment
  the pattern already carries — its own, plus the residual's on its
  share-weighted arm — and the moment `ContinuousAero` would place for the same
  panel, its force at the quarter chord plus `panel_couple`. Both the force and
  the moment now come out exact rather than approximate. The divisor carries the
  airfoil frame's own triple product, the axes not being quite orthogonal on a
  billowed panel.

### Added

- `check_live_polar(mode, wing; panel_idx)` runs one XFoil solve against the
  live polar a panel is flying, on the state the model is in now. It defaults to
  the panel NeuralFoil is least confident about.
  `Live polar off the trained shape range` says the network is extrapolating;
  this says whether the extrapolation is any good, which a confidence cannot.
- `AeroPressure(; live_polars=true)`: instead of reading a tabulated `(α, δ)`
  polar, every panel's polar is regenerated from its deformed shape at each VSM
  solve. The structural points a panel already scatters its load onto are read
  as control points, their offset off the deformed chord line deforms the
  panel's Kulfan fit analytically, NeuralFoil is evaluated on a grid of angles
  about the panel's own α, and those values are written into the panel's own
  polar table. This is for a chord that bends into a shape no single hinge angle
  stands for — the flap angle δ then carries nothing, so it is dropped from the
  generated equations entirely and only α stays live in the RHS. Requires
  VortexStepMethod 4.3.

  The surface traction pattern is regenerated from the same deformed shapes once
  the solve has converged: the contour offset by the deformation's own camber
  increment, `Cp` from one batched NeuralFoil pass at the converged α, and skin
  friction from the flat-plate closure at each panel's own Reynolds. Force and
  placement therefore both follow the deformation — previously a deformation
  changed how hard a panel pulled but not where it pulled, and every traction
  normal was the reference section's. The contour's chord fractions and its
  nodes' assignment to structural points do not move, so the scatter stays
  continuous.

  A live polar spans only the angles it was sampled over and holds its last
  value past them, so bring-up still wants damping enough to keep the solve
  inside that range — on the SK100, `start_world_damping` 300 rather than 20,
  which is what the placed geometry needs before the wing reaches 20 m/s in 10
  ms and slews α 20° in a single step. Past the range the answer is stale rather
  than wrong-signed. Settled at 12 m/s the two agree; at 20 m/s the live polars
  settle to a peak tube bending margin of 1.02 against the tabulated 4.71, with
  3 of 49 joints past the collapse knee rather than 10.
- The continuous VSM modes carry VortexStepMethod's `flow_curvature` moment
  increment, the thin-airfoil `Δcm = -(π/4)·q̂` of a section rotating about its
  own spanwise axis. They re-derive the panel force symbolically instead of
  reading `solve!`'s output, so the term was absent from them while `AeroDirect`
  and `AeroLinearized` had it, and enabling the solver flag silently meant two
  different models. Each panel's rate comes from its sections' own trailing
  minus leading edge apparent wind, so a deforming wing gets the true
  per-section rate — twist and flapping included — rather than a projection of
  one body rate. Off by default, and read from the wing's solver, so no mode can
  have the term while another does not.
- Unsteady aerodynamics on the particle VSM modes, carried by the wing's new
  `UnsteadyAero` (`wing.unsteady`) and off by default. `unsteady_aero(wing)`
  reads a wing's settings back and `apply_apparent_mass!` puts the entrained air
  on the structure. See the "Unsteady aerodynamics" docs page.
  - `apparent_mass` gives the wing's nodes the air they entrain, the thin-plate
    `ρ·π·c²/4` per unit span, spread over them by the weights that carry each
    panel's force. It sits in a new `Point`/`Body` field rather than
    `extra_mass`, because entrained air resists acceleration but has no weight,
    and it lands on whatever integrates the node's translation — the node
    itself, the body that places it, or the two beam bodies its joint spans. It
    is a scale, `1` being the thin-plate value, translational and isotropic.
  - `wagner` adds the two-state Wagner lift lag: circulation reaches its steady
    value over `φ(s) = 1 - A₁·exp(-b₁·s) - A₂·exp(-b₂·s)` semi-chords rather
    than at once. Two states serve the whole wing, driven by its mean apparent
    wind, and the one deficiency shifts every panel's angle of attack before the
    polars are read. `wagner_gains` and `wagner_rates` are registered
    parameters, so retuning a lag syncs instead of rebuilding; only the on/off
    switch is structural.

## v0.15.2 22-08-2026

### Added

- `reinit!(sam, integrator; solver=integrator.alg, kwargs...)`, which defaults
  the solver to the one the integrator was built with. The three-argument method
  takes it as a required argument, and a caller passing a bare `FBDF()` made the
  monolith backend compile its whole right-hand side again at
  `ForwardDiff.Dual`.

### Fixed

- `sim_reposition!` reinitializes with the solver `init!` built, rather than
  rebuilding one and dropping its `linsolve`.
- Every wing node's `aero_force_b`, and the `aero_force_x/y/z` log channels with
  it, stayed zero under the modes that scatter panel loads (`AeroPressure`,
  `ContinuousAero`), so `plot` drew no aero arrows on such a wing. Both backends
  now read the load back out of the model, where it already existed.
- `update_from_sysstate!` set the wing's aero and tether loads to `NaN` to keep
  them unplotted, but `plot` skips a load only with `iszero`, so one `NaN`
  blanked every aero arrow in a replay through the shared adaptive scale. All
  four loads and the per-point forces are restored from the log instead.
- `plot` drops non-finite loads before scaling the rest.

### Changed

- `point.aero_force_b` is the one place a per-point aerodynamic load is read
  from, filled by refresh (`AeroDirect`) or read back out of the model. The
  private `aero_point_forces` and `stores_point_force` are removed.

## v0.15.1 22-08-2026

### Changed

- The `KernelBackend` now defaults to a sparse Jacobian and a
  `KLUFactorization`, through the new `default_sparse` and `default_linsolve`
  dispatches; `init!` takes `sparse=nothing` and a `linsolve` keyword to
  override either. On a 392-point beam the Jacobian is 1301 states square and
  5.15% dense, so the dense factorization was most of a step: 35.5 ms against
  50.2 ms, and 62.1 steps/s against 4.8 with eight models stepping at once,
  because a dense factorization opens a BLAS thread pool per worker over the
  same cores. `sparse` is part of the model bin's name, so a kernel model builds
  once more. The `MonolithBackend` keeps a dense Jacobian and LinearSolve's own
  choice.
- Changing `g_earth`, `wind_vec` or `cd_tether` no longer misses the model cache
  and rebuilds. `get_set_hash` drops the three: each enters the equations as a
  flat parameter that `sync_params!` refreshes from `set`, not as a literal, so
  one build serves a sweep over them.

### Fixed

- `sync_params!` no longer allocates, and is 22.13 ms -> 0.41 ms a call on
  SK100. It ran twice a step and was ~90% of it. Its readers were held in a
  `Vector{Any}`, so each of the 22756 parameters cost a dynamic dispatch, and
  `PathReader` walked its path as values, so a field name was never a
  compile-time constant and every step of the walk boxed both its argument and
  its result. Readers are now grouped by concrete type behind a function barrier
  (`ReaderGroup`) and `PathReader` carries its field names in the type, keeping
  container indices as values so components of one kind still share a compiled
  reader. The remaining allocation is the MTK `setp` setter and the `Any`
  buffers of the monolith's array and callable kinds, a fixed 256 B that does
  not grow with the model.
- `SystemStructure.wings` is a stored field instead of a `getproperty` that ran
  `filter(is_wing, bodies)` and rebuilt a `NamedCollection` on every access.
  Same objects as `bodies`, not copies.
- `test_getter_allocations.jl` covers both backends and asserts `sync_params!`
  against `SYNC_ALLOC_BUDGET`. It previously ran only the default backend, so
  the `KernelBackend` sync paths had no allocation coverage at all.
- `reinit!` cold-starts the initial VSM solve, through a `cold_start` keyword on
  `refresh_aero!` and the refresh chain below it. `VortexStepMethod.solve!`
  warm-starts from the circulation left in `solver.sol`, and past stall the
  iteration has more than one fixed point, so the frozen aero after a state
  restore depended on what had run on that model before: restoring one `.arrow`
  twice with a run in between moved the effective sectional α by 0.29° and the
  frozen panel traction by 84 N. Per-step refreshes still warm-start, so
  stepping costs the same.
- `KernelBackend` reads a fitted (`KINEMATIC`) wing's pose back from the
  weighted blend of its reference points, as `MonolithBackend` does. The
  read-back kept only the first point of each `WeightedRefPoints`, so a wing
  whose `origin_idx`, `z_ref_points` or `y_ref_points` name several points
  reported the wrong `pos_w`/`vel_w`/`R_b_to_w` — and the apparent wind and
  reported scalars derived from them. The dynamics were already correct; only
  the struct read-back was not.
- The `KernelBackend` reads each point's `va_b` out of its `AeroInflowPoint`
  rather than refitting it in `wing_kinematics_from_points!`. The compiled model
  already computed the apparent wind the aero is solved on; the refit overwrote
  it from `set.wind_vec` and a frame fitted from the reference points, and
  disagreed: on the parked V3 it gave `va_b` of 19.26 m/s against the monolith's
  16.44, an AoA of −0.25° against +5.64°, and an L/D of 70-or-`NaN` against 9.5.
  The refit stays as the fallback for a wing with no aero instance.
- A `RIGID_DYNAMICS` wing gets its reported scalars on the `KernelBackend`:
  `heading`, `elevation`, `azimuth`, `course`, `aoa`, `turn_rate` and `turn_acc`
  were only filled for a fitted wing, so on a rigid one they held whatever the
  struct was built with. Read from the body's `frame` output and its angular
  acceleration `alpha_b`, now observed by every body kernel, through the same
  `wing_scalar_kinematics` the monolith's equations use.
- Every anchored (`BODY_STATIC`) point on the `KernelBackend` reports the wind
  at its own height and its net force. `ride_wrench_variables` declared only
  `total_drag`, where `point_eqs!` writes all three for every point, so on a
  `RIGID_DYNAMICS` wing every structural node read back a zero `wind_vec` — the
  full 15.5 m/s of the test case — and a `force` left at whatever the struct was
  built with. Bound once in `ride_wrench_eqs`, so `RideWrench`,
  `HermiteRideWrench` and `TwistNodeWrench` all report them; the last of these
  had `ride_load`'s body inlined and now shares it.
- The `KernelBackend` fills the read-backs a fitted (`KINEMATIC`) wing and its
  points were missing, so both backends scatter the same `SystemStructure` out
  of a solved model: `acc_w` and the
  `elevation_acc`/`azimuth_acc`/`distance_acc` that follow from it — a fitted
  body reports no acceleration of its own, so it is blended from the origin
  reference points, as `wing_eqs!` blends it — `Q_p_to_w` and `ω_p`, aliased to
  the body frame the same way `wing_eqs!` binds them, each point's `total_mass`
  as its segments' rest lengths move, and `va_b` for the points off the wing.
  `test_backend_parity.jl` diffs every real-valued field of every component
  across the two backends, at `init!` and after stepping, so a read-back added
  to one and not the other fails there.
- The live aero modes (`ContinuousAero`, `AeroPressure`) take a panel's chord
  direction between the two bounding sections blended by panel spacing, as
  `VortexStepMethod` does, instead of at their midpoint. On a non-uniformly
  panelled wing that midpoint tilted the whole airfoil frame: the angle of
  attack was off by 1.9e-4 rad and the per-panel force by 0.1%, a systematic
  bias that partly cancels spanwise and so never showed up in the wing total.
  Every per-panel quantity now matches VSM to machine precision. The weights are
  refreshed with the frozen circulation, so they follow the deforming mesh. They
  are a new parameter, so a model cached before this needs `remake=true` once.
- `set.g_earth` reaches a `PARTICLE_DYNAMICS` wing's points. `WingNodePoint`
  called `point_acceleration` without it, so the default keyword baked 9.81 into
  the generated equations while the parameter it registered went unread: gravity
  was frozen at whatever the model was built with. `g_earth` is now a positional
  argument of `point_acceleration` and `point_net_force`, so a caller that omits
  it is a `MethodError` at build. It was only reachable once gravity left the
  set hash; before that a changed `g_earth` rebuilt the model and re-baked the
  literal.
- `test_continuous_modes.jl` checks the per-panel equations by evaluating
  `panel_force_eqs` on the panel's own geometry and inflow, instead of reading
  `aero_N.chord`/`.alpha`/`.panel_force` out of the compiled system. Only the
  `MonolithBackend` builds one, so the testset was skipped on the
  `KernelBackend`, and the names it reached for were codegen artefacts a rename
  would break.

## v0.15.0 21-08-2026

### Added

- `Body.fix_static` freezes a body where it is, like a `Point`'s. A parameter,
  not a `type`, so it can be toggled on a built model. Also a `fix_static`
  column of the `bodies` YAML table.
- `TimoshenkoJoint` and `ElasticJoint` take `damping`, the Rayleigh
  stiffness-proportional coefficient β in seconds (`C = βK`), giving modal ratio
  `ζ = βω/2` per mode. Being built from `K` it inherits the rigid-body null
  space, so it damps every deformation mode — axial, shear, bending, torsion —
  and rigid motion by construction not at all.
- `test_joint_invariance.jl` runs every joint type through the same invariants:
  rigid motion conserves linear and angular momentum, and the logarithmic
  decrement `exp(−2πζ/√(1−ζ²))` pins β to `K`. New joint types inherit the check
  by being added to `joint_cases`.
- `Body` takes `world_frame_damping` and `body_frame_damping`, per-mass [1/s]
  translational damping of the COM velocity resolved on the world and body axes,
  next to the angular `damping` it already had. A beam wing carries most of its
  mass in bodies, and nothing else damps their rigid motion; a per-axis
  body-frame coefficient such as `[0, 0, 20]` resists flapping normal to the
  wing without slowing flight along it. Settable from the `bodies` YAML table.
- `set_world_frame_damping` and `set_body_frame_damping` take a component
  vector, so `set_body_frame_damping(sys.bodies, damping)` reaches the bodies.
  The `SystemStructure` methods still mean the points.
- `Body` takes `wing`, resolved to `wing_idx`/`wing_ref` like a `Point`'s,
  naming the parent wing its `body_frame_damping` resolves against. A body
  without one damps its own velocity on its own axes, as before.
- `set_angular_damping(bodies, damping)` sets the per-axis spin damping the
  `Body` already carried but nothing exposed: `dω/dt -= c .* ω` on the body's
  own axes, so `[0, 20, 0]` resists rotation about `y` alone. It damps
  _absolute_ spin, unlike the joints' Rayleigh `damping`, which has no
  rigid-body effect. Its docstring said N·m·s; the coefficient is applied to
  angular acceleration, so the unit is 1/s.

### Changed

- BREAKING: `plot(::SystemStructure, ...)` now extends `MakieControlPlots.plot`,
  which is the `plot` exported into your session and shadows Makie's, so bare
  `plot(sys_struct)` works instead of erroring. Replace qualified
  `Makie.plot`/`GLMakie.plot` calls with `MakieControlPlots.plot`, or use
  `import GLMakie` instead of `using GLMakie`. `plot!` is unchanged.
- BREAKING: `Body.damping` is now `Body.angular_damping`, as its `Wing`
  constructor kwarg already called it, and the `bodies` YAML key follows. A body
  carries three damping fields and the bare name said nothing about which one it
  was; the other two are `world_frame_damping` and `body_frame_damping`, and the
  joints' Rayleigh β keeps the name `damping`.
- BREAKING: the wing constructors (`Wing`, `VSMWing`, `PlateWing`) and the
  `wings` YAML block drop `y_damping`; `angular_damping` is now the whole
  per-axis vector and defaults to `[0, 150, 0]`, which is what the pair used to
  produce. It was assembled as `[a, a + y_damping, a]`, so one axis was special
  cased and the two knobs summed on it. A scalar is still broadcast to all
  three.
- BREAKING: `TimoshenkoJoint` and `ElasticJoint` no longer accept
  `damping_trans`/`damping_rot`; use `damping` (Rayleigh β, seconds). A dashpot
  on relative node velocity has no rigid-body null space, which is what made the
  bug below possible. YAML `timoshenko_joints`/`elastic_joints` tables must
  rename those two columns to a single `damping`.
- `SystemStructure`'s `diff_vars` property is now `state_vars`, and covers every
  component regardless of `DynamicsType` rather than only the ones the ODE
  integrates. Whether a field is integrated, torn as an iteration variable or
  eliminated is the compiler's choice and has changed across stacks, so a
  `DynamicsType` filter was a standing bet on that choice. It is what
  `validate_sysstate_roundtrip` scrambles, and the bet losing is what let the
  `KINEMATIC` restore below go unnoticed.

### Fixed

- A tether with no winch takes its rest length from `params.tethers[i].len`
  instead of an integrated `tether_len` whose derivative is zero, so writing
  `tether.len` retrims a running simulation and no longer needs a `reinit!` to
  take effect. `MonolithBackend` only: `KernelBackend` already read it as a
  parameter, so the two backends disagreed on whether a mid-run rest-length
  change did anything.
- The Comer-Levy `TubeRigidityLaw` bending branch is C1 across its knee. Knee
  and tail were fitted independently, so `dEI/dκ` jumped from 0 to 29% of
  `EI0/κ_knee` there and the beam Jacobian was discontinuous. The post-knee
  deficit now leaves the knee at the slope the linear branch sets, and its
  exponent is fitted on relative moment error rather than the log deficit (rms
  1.15% → 1.10%). The Breukels fit is unchanged — its knee is a fitting anchor,
  not a tangency.
- A `Body`'s `body_frame_damping` damps its velocity _relative to its parent
  wing_, resolved on the wing's axes, through the same `body_frame_damp_accel` a
  wing node uses — so it damps deformation and leaves rigid flight alone, and it
  means on a body what it already meant on a point. It damped absolute velocity
  on the body's own axes, which for an isotropic coefficient is the same
  operator as `world_frame_damping`. Both backends reach it through the one
  `body_integration`, the kernel declaring `wing_frame`/`wing_velocity` inputs
  for a parented body.
- `store_induced_velocity!` reads `BodyAerodynamics.AIC` as
  `(panel, panel,
  component)`, the layout VortexStepMethod uses so each
  `AIC[:, :, k]` slice is a contiguous BLAS matrix. It still indexed the old
  `(component, panel, panel)` order, so every continuous VSM mode threw a
  `BoundsError` on the first aero refresh.
- The continuous VSM modes sample each section's inflow at the 3/4-chord
  collocation point of thin-airfoil theory instead of at mid-chord. The
  mid-chord average kept a chord twist rate almost entirely out of `alpha`, so
  the aero gave no pitch damping to a chord pivoting near mid-chord and damping
  of the wrong sign to one pivoting between mid- and 3/4-chord. A wing whose
  chord is structurally flexible had no aerodynamic resistance to flapping.
  Inert in steady flight: with no twist rate the LE and TE stations share a
  velocity, so every chordwise fraction gives the same inflow.
- `setproperty!(::Body, ...)` converts to the field type, as Julia's default
  does. Its custom method called `setfield!` raw, so assigning an `SVector` to
  any `KVec3` field of a body threw a `TypeError` where the same assignment on a
  `Point` worked.
- `TimoshenkoJoint` no longer brakes rigid rotation of the bodies it connects.
  Damping resisted the raw relative node velocity `Δv`, whose transverse part is
  by construction the corotational frame's own spin and carries no strain, so a
  chain of elements assembled a spurious rotational damping tensor
  `Σ c_t(|d|²I − ddᵀ)`. On the V3 beam kite that was 15332 N·m·s/rad about yaw,
  giving the wing a 6 ms rotational decay constant and cutting steering
  authority by three orders of magnitude.
- `elastic_joint_wrench` differentiated its deformation measure in the wrong
  frame: the spring uses `Ra'(x_b − x_a)` but the damper used `Ra'(v_b − v_a)`,
  missing the transport term `ω_a × (x_b − x_a)`, so a rigid rotation registered
  as relative anchor velocity whenever the anchors were separated.
- `elastic_joint_wrench`'s wrench pair was not self-equilibrated with separated
  anchors — it injected a net torque `(x_b − x_a) × F`, so angular momentum was
  not conserved. Body A now also carries the transport couple of the transmitted
  force.
- `update_from_sysstate!` now rebuilds the principal-frame state of every body,
  where it skipped `KINEMATIC` ones. Their `com_w`/`com_vel` are not logged, but
  the wing origin point's position alias-eliminates onto `com_w`, so restoring a
  state left that one point at the geometry-file position while every other
  point moved to the logged one. The body frame fitted from those points came
  out tens of degrees off and `init!` then failed DAE initialization
  (`your u0 did not satisfy the initialization requirements`) — on the V3 at 250
  m tether the origin was stranded 148 m from the rest of the kite.

## v0.14.0 13-08-2026

### Added

- `Pulley` gains `efficiency`, the fraction of line tension its sheave passes
  on. Friction is `(1 − efficiency) · line_tension`, carrying the sign of the
  motion, smoothed over `friction_epsilon` [m/s], with `line_tension` the mean
  of the two leg tensions — so the default is the textbook `T_out = 0.95·T_in`.
  It replaces the hidden `pulley_damp`, fixed at 5.0 and neither readable nor
  settable from the `SystemStructure`. A sheave's losses scale with load, not
  speed, and `pulley_damp` scaled with rope travel and rope mass, so no fixed
  coefficient reproduces it — expect to retune models that leaned on it to
  settle. 0.95 is a sealed ball bearing (published 0.94–0.97; 0.88–0.92 for a
  bronze bushing); 1.0 is ideal. A slack pulley is now frictionless and coasts,
  but has no tension driving its split either.
- `Pulley` also takes `damping` [N·s/m] and `brake`, both off by default and
  neither a sheave property: `damping` settles a ringing split while debugging,
  `brake` freezes the split to isolate whether a problem comes from the rope
  redistributing at all. Constructor keywords and YAML columns.
- `Segment(...; compression_damping_frac)` scales the damping term the way
  `compression_frac` scales the stiffness, live in the right-hand side. The
  default `1.0` is the previous behaviour: damping unaffected by compression and
  the force continuous at `len == l0`, but the damping ratio jumping by
  `1/sqrt(compression_frac)` as the segment goes slack. Set it equal to
  `compression_frac` for one damping ratio on both branches;
  `compression_damping_frac = compression_frac = 0` makes a slack segment carry
  no force at all, at the cost of a force step of
  `unit_damping / l0 ·
  spring_vel` at the crossing. Also a `segments` YAML
  column.
- `Tether(...; compression_frac, compression_damping_frac)` hands both to every
  segment a Route 2 tether generates, so an auto-generated line is no longer
  stuck at the `Segment` default of 0.1. Also `tethers` YAML columns. A tether
  without a winch holds its length fixed, so this is how to split one line into
  several segments without writing out every segment and intermediate point.
- `youngs_modulus` [Pa] and `damping_per_stiffness` [s] are ordinary `Segment`
  and `Tether` inputs (YAML columns and constructor keywords) — the
  diameter-independent forms of `unit_stiffness` and `unit_damping`, so one
  material can be shared across diameters:
  `unit_stiffness = youngs_modulus ·
  pi · (diameter_mm/2000)^2`,
  `unit_damping = damping_per_stiffness ·
  unit_stiffness`. Giving both forms
  of the same quantity is an error.
- New top-level YAML block `variables`: named values reused across the file. A
  variable name in any column is replaced by its value (numbers, strings,
  lists), and variables may be defined in terms of each other:
  ```yaml
  variables:
    bridle_comp: 0.01
    dyneema: {
      youngs_modulus: 55.0e9,
      damping_per_stiffness: 0.00077,
      density: 724.0,
    }
  ```
  A variable holding a mapping is a _multi-variable_: it fills the columns it
  names at once, so the row carries one entry for the group
  (`- [bridle, le_left, kcu, nothing, 1.0, dyneema, bridle_comp]`). Its fields
  must match the columns at that position; in a dict row they are merged in
  without overriding the row. Naming a variable after a component is an error.
- A `SysLog` now carries the whole differential state, so a simulation can be
  restarted from any logged step. `update_sys_state!` writes point velocities,
  per-frame body turn rates, twist-surface twist and rate, and pulley length and
  rate alongside what it already wrote; `update_from_sysstate!` reads them back,
  rebuilding the body rates into the principal frame the ODE integrates. It no
  longer zeroes point velocities or caps twist surfaces at four. Pass
  `precision=Float64` to `SysState(sam)` or `Logger(sam, steps)` for a state
  that reproduces `integrator.u`; the default stays `Float32` so telemetry logs
  keep their size. Restart with `init!(sam; remake=false, reinit_sys=false)`,
  which seeds the integrator from the structure instead of the CAD frame.
- The YAML loader reads back state that previously had no column: point `vel_w`
  (zeroed unconditionally, discarding the column it already accepted), pulley
  `sum_len`/`len`/`vel`, twist-surface `twist`/`twist_vel`, tether `len`, winch
  `vel`/`set_value`, and body `vel`/`Q_b_to_w`/`omega_b`.
- `init!(sam; analytic_jacobian)` gives the solver a Jacobian instead of
  differentiating the right-hand side numerically. On by default on the
  `KernelBackend`, where each kernel is differentiated once at its own width and
  composed through the constant wiring, rather than differentiating the whole
  model `n_states / chunk` times. Off by default on the `MonolithBackend`, where
  it selects MTK's symbolic `jac=true` — expensive to build, and with no
  derivative for a registered numerical leaf such as the wind profile. `nothing`
  (the default) takes the backend's `default_analytic_jacobian`. The choice is
  part of the serialized model's name, so the builds cache apart.
- `test_init!` round-trips every model the suite builds through a `SysLog` on
  disk and checks `integrator.u` is unchanged, covering the `DynamicsType`
  combinations that decide whether a component is part of the state at all.

### Fixed

- A segment's damper resists the rate of change of its extension `len - l0`, not
  the rate its endpoints separate. The two agree wherever `l0` is a fixed
  parameter, so an ordinary segment is unchanged — but a pulley leg and a
  winched tether member carry `l0` as a state, and the difference is that
  state's own rate. A pulley's rope split therefore appeared in no damper at all
  while the legs' dampers still drove it: a one-way velocity coupling that could
  feed the split rather than settle it, growing with `damping_per_stiffness`, so
  raising the damping could make a model _less_ stable. The split now carries
  `c_left + c_right`. Reeling likewise no longer charges a tether a damper force
  proportional to reel speed, so a reel-out loop is genuinely less damped than
  the old model made it look — expect controllers tuned against it to need
  retuning.
- A segment whose `unit_damping` is missing or `nothing` takes the `Segment`
  constructor default (derived from the settings) instead of being silently
  forced to zero by the loader.
- Air density is clamped at ground level on every path. A rigid wing read
  `calc_rho` at its raw world `z`, extrapolating the atmospheric model backwards
  below `z = 0` while points and segments clamped. The clamp now lives in
  `air_density`, which every consumer in both backends reads through.
- `init_stretched_length` placement carries a beam wing. A `BODY_STATIC` point
  riding a Timoshenko element has `body_idx == 0` and holds its association in
  `joint_idx`, so collecting bodies by `body_idx` alone left every beam body
  behind while the tether's free end moved; the ride constraint then snapped the
  points back, and the only symptom was a non-converged VSM solve much later.
  Bodies reachable through the beam graph now translate too — including bodies
  carrying no point of their own — stopping at `STATIC` bodies, since a beam
  with a clamped end deforms rather than translates.
- The Breukels inflated-tube correlations error instead of quietly returning a
  negative rigidity outside the range they were fitted in. The bending slope
  changes sign below roughly 38 mm radius at 0.25 bar (a 20 mm tube at 0.1 bar
  gave `EI0 = -94 N·m²` and a negative `EA`), and the torsion factor `c2`
  changes sign for a fat tube above 1 bar. Non-positive radius or pressure is
  rejected too.
- `set.g_earth` is settable after construction on the `KernelBackend` too. A
  particle point baked the value in while a rigid body and a ride point read the
  registered parameter, so changing gravity moved the bodies and left the tether
  points falling at the old rate.
- An aero panel that maps to no twist surface reads a zero deflection on both
  backends; the monolith indexed `twist_surface_delta` at `0`. Unreachable
  today, since every panel is assigned a nearest flap whenever a wing has one,
  but live for any mode that emits a partially-zero map.
- `plot` handles a system with fewer than three winches. The reel-out,
  winch-force and set-torque panels pulled three channels out of every logged
  sample no matter how many the log held, so a one-winch log threw a
  `BoundsError` before the figure was drawn. Each panel now takes its channel
  count from the log.

### Changed

- BREAKING: the `materials`, `elements` and `segment_properties` YAML blocks are
  removed, and loading a file with one errors. A material is now a
  multi-variable listing `youngs_modulus`, `damping_per_stiffness` and
  `density`, all ordinary columns.
- BREAKING: `Winch.f_coulomb` and `Winch.c_vf` are renamed to
  `Winch.coulomb_friction` and `Winch.viscous_coefficient`, matching the
  `Pulley` fields and saying which is a force [N] and which a coefficient
  [N·s/m]. The `settings.yaml` keys are owned by `KiteUtils` and keep their old
  names.
- The monolith and the `KernelBackend` evaluate one shared definition of every
  equation instead of a copy each. Points, segments, pulleys, rigid bodies,
  ref-point wing frames, twist deformation, the ground wind fallback and the
  aero wiring route through the helpers in `src/components/`, which previously
  had the `KernelBackend` as their only caller. A backend still chooses how
  equations are assembled; it no longer chooses what they are.
- `compression_frac` defaults to 0.1 in both `Segment` constructors. They
  disagreed — 0.0 from the settings-based one, 0.1 from the direct one — which
  was an oversight. 0.1 lets a rope push back with a tenth of its tensile
  stiffness and is what models were built against; 0.0 leaves a slack segment
  with damping but no stiffness, which a pulley model does not survive.
- Nine segment observables (`stiffness`, `damping`, `segment_height`,
  `segment_vel`, `segment_rho`, `wind_vel`, `va`, `area`, `app_perp_vel`) are
  gone and the point observable `point_drag_force` is renamed `point_aero_drag`,
  where it shadowed the shared function of that name. They were intermediates of
  the monolith's own copy of the force law, which no longer exists. Nothing in
  the package, examples or docs reads them and no exported accessor returns a
  symbolic variable by name, so only code indexing the generated system directly
  is affected. `len`, `spring_force`, `spring_force_vec`, `spring_vel`,
  `unit_vec`, `segment_vec`, `rel_vel`, `l0` and `total_drag` are unchanged.
- The aero tests are layered, so a new mode inherits its contracts by being
  added to a table rather than by someone remembering to write tests.
  `test_aero_modes.jl` holds what every mode owes and now enrolls
  `AeroPressure`; the new `test_continuous_modes.jl` holds what the continuous
  modes owe — including that per-section inflow is gathered per section and
  never as a wing-wide mean; `test_continuous_aero.jl` and
  `test_pressure_aero.jl` keep only what is specific to one mode.

## v0.13.0 06-08-2026

### Added

- New aero mode `ContinuousAero` (`PARTICLE_DYNAMICS`, YAML
  `aero_mode: continuous`): frozen-circulation VSM with the full force assembly
  in the symbolic RHS. The low-frequency refresh runs only the circulation solve
  (`solve_base!`) and freezes each refined panel's induced velocity (`AIC·γ`);
  every RHS step re-derives panel geometry from the live strut points (frozen
  mesh-interpolation weights), effective angle of attack, polar coefficients
  (registered `Dual`-safe lookups on the panel polars), lift/drag directions,
  and forces. Forces therefore respond to wing motion between VSM updates —
  aerodynamic damping through the changing angle of attack — unlike
  `AeroDirect`'s piecewise-constant forces. All per-panel quantities (`alpha`,
  `cl`, `q_dyn`, `panel_force`, …) are observable component variables. The mesh
  weights enter the model-cache hash via `aero_hash_id`.
- New aero mode `AeroPressure` (YAML `aero_mode: pressure`): like
  `ContinuousAero`, the VSM solve stays discrete (every `vsm_interval`) but the
  per-panel force is re-derived symbolically each RHS step from the frozen
  circulation, then _scattered onto arbitrary `BODY_STATIC` surface points_ via
  a build-time station→point map. Gives aero damping between refreshes.
- Deflecting flap for `AeroPressure`: a live flap deflection `δ` (the signed
  angle between two structural bodies about a hinge, relative to rest) feeds the
  `(α, δ)` polars each RHS step. Modeled by extending `TwistSurface` with a
  `KINEMATIC` twist source and flap fields (`flap_body_idxs`, `flap_axis`,
  `flap_chord_refs`, `flap_rest_delta`) — no new component. Needs
  VortexStepMethod's 3-arg `calculate_cl/cd/cm(panel, α, δ)`.
- Structural `Body` / beam subsystem: a rigid `Body` type with new `BODY_STATIC`
  and `KINEMATIC` `DynamicsType`s, `ElasticJoint`, and `TimoshenkoJoint` (a
  2-node corotational Timoshenko beam element between two `Body`s — the
  distributed-stiffness counterpart of `ElasticJoint`; chains form a beam,
  closed-form validated for bending with the `PL/kGA` shear term, axial and
  torsion, with nonlinear-stiffness dropoff).
- Beam-anchored points: a `BODY_STATIC` point rides a body rigidly, or a joint's
  corotational-Hermite deformed centerline (`joint_idx`, `beam_frac`,
  `beam_offset_b`, splitting its load to both nodes). A body-anchored point's
  body-frame offset is auto-derived from its `pos_cad`, so riding a body needs
  only `body=`/`joint=` plus a CAD position.
- Per-segment material `density` [kg/m³]: each `Segment` and `Tether` carries
  its own `density` (from the YAML `materials` table), replacing the single
  global `set.rho_tether` in mass calculations. Falls back to `set.rho_tether`
  when unset.
- Twist-surface restoring `stiffness` [N·m/rad] for the twist DOF; the number of
  unrefined VSM sections can be read from the YAML.
- `SymbolicAWEModels.ObjAdapter`: mesh mass properties (`center_of_mass`,
  `calculate_inertia_tensor`, `unit_inertia_from_obj`) computed in this package
  instead of VortexStepMethod, plus per-wing `mass`, `com` and `unit_inertia`
  (symmetric per-unit-mass 6-vector `[Ixx,Iyy,Izz,Ixy,Ixz,Iyz]` in [m²], scaled
  by `wing.mass`) as YAML columns and `VSMWing` kwargs. When omitted they
  auto-compute from an `.obj` if one is supplied, else fall back to point-mass
  inertia. `create_vsm_wing` generates the aero YAML from an `.obj` at panel
  resolution and caches it under `<model_dir>/obj_geometry`.
- `PrincipalFrameMethod` (`EIGEN_DECOMP`, `Y_ROTATION`) is exported, selecting
  how a wing's inertia tensor is diagonalized.
- `position_slots(sys_struct)` is exported: the index layout of a `SysState`'s
  `X/Y/Z` arrays (`points`, `panel_corners`, `wings`, `bodies`, `total`), so a
  log consumer no longer has to reconstruct the packing by hand.
- Replay camera pan and zoom controls (#250).
- Makie plotting: beam elements are drawn as instanced see-through tubes
  (`show_beams`), every layer can be toggled from the replay window, and an
  `aero_mapping=true` overlay draws the `AeroPressure` station→point scatter map
  (`m` toggles it live).
- `examples/inflated_beam_fit.jl`: fits a per-joint nonlinear bending law of a
  pressurised tube from the Comer-Levy wrinkled-section theory and validates the
  resulting `Body`/`ElasticJoint` chain as a cantilever.

### Changed

- BREAKING: winch motor dynamics are a type, not a builder function.
  `Winch.model` is now an `AbstractWinchModel` (`TorqueWinch`, the default, or
  `CascadedLengthWinch`) instead of a `Function`, and the builder hook is
  `winch_component(model, sys_struct, winch_idx; name, params)` rather than
  `default_winch_component`. _How to migrate:_ wrap a custom builder in a new
  `AbstractWinchModel` subtype with a `winch_component` method; a custom model
  also wants `is_builtin_winch` left at its `false` default so `init!` rebuilds
  the cached model.
- BREAKING: `Wing` and `RigidBody` are one `Body{A, D}` type. A wing is a `Body`
  that carries aero, so `sys.wings` is a filtered view of `sys.bodies`. The
  exported concrete types are `RigidWing{A} = Body{A, RigidDynamics}` (6-DOF
  pose) and `ParticleWing{A} = Body{A, ParticleDynamics}` (pose fitted from
  structural points), replacing direct use of the old wing struct.
- The dependency stack was widened to the newer MTK/Symbolics generation
  (Symbolics 7→8, SymbolicUtils 4→5, SymbolicIndexingInterface 0.3→0.4,
  OrdinaryDiffEqBDF 1→2/3, OrdinaryDiffEqCore 3→4, NonlinearSolve 4→5/6,
  SteadyStateDiffEq 2.5→3/4; `ADTypes` added as a direct dependency).
  `VortexStepMethod` compat is raised to `4` and `KiteUtils` to `0.11.9`.
  `mtkcompile` now keeps body-frame outputs (`body_pos_w`, `body_vel_w`,
  `body_R_b_to_w`) and pulley `l0` in the torn-out state set. Not breaking: the
  package version is part of the `.bin` cache filename, so upgrading
  auto-invalidates stale model caches.
- RHS performance: model parameters are flattened for faster evaluation (#234).
  Rigid-body / beam element handling improved with aero separated into its own
  component (#236); `winch.speed_controlled` restored after the winch refactor
  dropped it.
- BREAKING: the Makie extension is triggered by `MakieControlPlots` instead of
  `Makie`, so `using GLMakie` alone no longer enables `plot`/`replay`/`record`.
  _How to migrate:_ add `using MakieControlPlots` next to your Makie backend.

### Fixed

- Rigid-body pose drift: the principal quaternion ODE had no norm constraint, so
  `|Q_p|` drifted (~0.999) and `quaternion_to_rotation_matrix`, which assumes a
  unit quaternion, produced a non-orthonormal `R_b_to_w` (`det = |Q|^6`). Points
  were placed with that skewed `R` while the getters re-orthonormalised, giving
  a ~0.02 pose mismatch under twist. `rigid_body_eqs.jl` now carries a Baumgarte
  term `k(1 − ‖Q‖²)Q` that holds `|Q_p|` on the unit sphere.
- `AeroLinearized`'s refresh writes the operating-point force via
  `apply_direct_forces!` (as `AeroDirect` does), so `wing.aero_force_b` tracks
  the VSM solve without a second state sync.
- `quaternion_to_rotation_matrix` now normalizes internally, dividing each
  product by `sum(abs2, q)` so the result is orthonormal for any nonzero `q` (no
  `sqrt`, no measurable cost). The integrated attitude quaternion only meets its
  unit-norm constraint to solver tolerance — the Baumgarte term in
  `rigid_body_eqs.jl` bounds the drift but, with a time constant of
  `1/(2*quat_norm_gain) = 0.05 s`, cannot remove it within a step. At the
  default `abs_tol = rel_tol = 0.01` this left `R_b_to_w` off orthonormal by up
  to 3e-4 during fast transients, a scale error in every body-to-world transform
  feeding `moment_p`, `wing_pos`, `wing_vel` and `wing_acc`. Orthonormality
  error is now at machine precision (9e-16 worst case over a steering ramp, from
  3.3e-4). Regenerates the model equations, so cached `model_*.bin` files from
  earlier commits must be rebuilt (`remake=true`).
- Wing body frame no longer inherits the principal frame (#245). A wing that
  declares no `origin`/`z_ref_points`/`y_ref_points` used to get
  `R_b_to_c = R_p_to_c`, so the eigendecomposition's ambiguous axis assignment
  for near-equal principal moments leaked into the _body_ frame and flipped it
  ~90° relative to the VSM/CAD frame, causing growing lift/drag oscillation and
  eventual VSM non-convergence during reel-out. The fallback now keeps the CAD
  orientation (`R_b_to_c = I`, origin at the COM), leaving the principal choice
  a pure gauge — asserted by `test_principal_frame_invariance.jl`. For wings
  symmetric about the XZ-plane, `principal_frame_method=Y_ROTATION` selects the
  unique closed-form diagonalization instead of the permutation search.
- Transform placement (`apply_azimuth_elevation!`) is now roll-free and
  zenith/nadir-safe: the current radial is rotated onto the target with a
  minimal rotation (no dependence on the source frame, undefined at the zenith),
  with roll set solely by the heading step and a warning at the zenith/nadir
  where azimuth and heading are undefined. Placement of existing models shifts
  slightly.

### Removed

- BREAKING: the `WING`, `QUASI_STATIC` and `FIXED` `DynamicsType`s.
  `DynamicsType` is now `{DYNAMIC, STATIC, BODY_STATIC, KINEMATIC}`. _How to
  migrate:_ a `WING` point becomes `DYNAMIC` on a particle wing, or
  `BODY_STATIC` on a rigid/beam wing (riding a body via `body=`/`wing=`, or a
  beam element via `joint=`); `QUASI_STATIC` becomes `DYNAMIC`; a `FIXED` twist
  surface becomes `STATIC` (`KINEMATIC` is the new geometry-prescribed twist
  source, e.g. a flap hinge). Wing↔point membership is no longer a point type:
  it is carried by `twist_surface.point_idxs`, exposed as `point.is_wing_node`,
  and queried with `wing_frame_member(point, wing_idx)`.
- BREAKING: `auto_create_twist_surfaces!`. A section-coupled (VSM) RIGID wing
  that declared no `twist_surfaces` used to have one `DYNAMIC` `TwistSurface`
  invented per LE/TE section; this silent black box is gone and the case now
  raises a clear error. _How to migrate:_ declare explicit `twist_surfaces`
  covering the wing's LE/TE structural sections (see
  `data/2plate_kite/rigid_structural_geometry.yaml`). `AeroNone` wings are
  unaffected (they never coupled to sections).
- The redundant wing-side `point_idxs` list in YAML (non-breaking: the loader
  never read it — wing↔point membership is the point row's `wing_idx` column).
  Drop it from your wing definitions.
- BREAKING: `update_segment_forces!`. Segment forces are written by
  `update_sys_struct!` every step, so the separate refresh call had no remaining
  caller.
- `src/tether_properties.jl` and its `calc_spring_props` / `calc_tether_props` /
  `in_percent_band` helpers (unexported, unused — the one-segment equivalent
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
  `src/aero_modes/`; chosen via the `aero` kwarg or the YAML `aero_mode` column.
  The mode's `aero_component(mode, sys_struct, wing_idx; name)` returns a
  subsystem wired at a fixed body-frame connector contract per `dynamics_type`
  (RIGID: `va`, `rho`, `R_b_w`, `omega`, `twist`, `twist_vel` → `force`,
  `moment`, `twist_moment`; PARTICLE: per-point `pos`/`vel`/`va`/ `rho` →
  `point_force`), validated by `validate_aero_component`. The generated RHS
  stays allocation-free (`test_bench.jl`).
- A custom aero mode needs exactly two methods: `aero_component` and
  `aero_mode_tag` (cache tag). Everything else is an optional hook with a
  working default, dispatched on the mode — lifecycle (`setup_aero!`,
  `remake_aero!`, `validate_aero_structure`, `resize_aero_state!`,
  `init_aero_state!`), low-frequency refresh (`refresh_rigid_aero!`,
  `refresh_particle_aero!`, orchestrated by `refresh_aero!`), diagnostics
  (`calc_aoa`, `normalized_inertia`), log-point visualization
  (`n_aero_log_points`, `write_aero_log_points!`, `read_aero_log_points!`,
  `restore_aero_twist!`), and live Makie rendering (`plot_wing_aero!` /
  `update_wing_aero_plot!`, with methods in the Makie extension). There are no
  `isa`/`is_vsm` branches anywhere in the pipeline, so a custom mode is never
  excluded from a code path it cannot extend. VSM state (solver, geometry,
  linearization buffers) lives in a `VSMEngine` carried by `AbstractVSMAero`
  modes; subtyping it inherits the VSM implementation of every hook.
- `normalized_inertia` returns per-unit-mass inertia [m²] for every mode — the
  VSM `ObjWing` mesh tensor is already normalized and is now passed through
  as-is, the default normalizes the WING-point point-mass inertia
  (`normalized_point_inertia`), and the single scaling by `wing.mass` happens in
  `setup_wing_frame!`.
- `has_custom_component(sys_struct)`: `init!` defaults `remake` to rebuild
  automatically when a custom winch/aero component is present (their equations
  are not captured by the model hash). Structural mode fields enter the cache
  key via `aero_hash_id`; the mode tag enters the cache filename.
- Flat-plate wings log a display quad per section (4 corners, square of side
  `sqrt(area)`, structural point at quarter chord) via the log-point hooks, so
  plate geometry shows up in `SysState` logs like VSM panels do.
- New `FIXED` `DynamicsType`: a twist surface whose twist is a prescribed
  control input (no differential state). Flat-plate surfaces use it.

### Changed

- BREAKING: `Group` is renamed to `TwistSurface` throughout (type, YAML section,
  and fields, e.g. `wing.group_idxs` → `wing.twist_surface_idxs`). Flat-plate
  surfaces are now 1-point `FIXED` `TwistSurface`s instead of a separate plate
  type.
- BREAKING: the wing types are merged into one `Wing` struct. `VSMWing` and
  `PlateWing` remain as constructor functions; the polar lookups and drag
  correction of a flat-plate wing live on its `AeroPlate` mode. The `wing_type`
  keyword is deprecated in favour of `dynamics_type`.
- `AeroNone` carries no VSM engine and needs no VSM geometry or `vsm_set`, so a
  pure rigid-body wing builds without any VSM setup (`VSMWing` accepts
  `vsm_set=nothing` for engine-less modes).
- The symbolic aero generation was restructured: `vsm_eqs.jl`, `plate_eqs.jl`
  and `linearize.jl` are replaced by a thin mode-agnostic wiring layer
  (`aero_eqs.jl`), the per-mode files under `src/aero_modes/`, and
  `twist_surface_eqs.jl` (formerly `group_eqs.jl`). `SystemStructure`
  construction is split into `setup_wing_frame!` (mass/inertia/body frame,
  aero-independent) and the mode-dispatched `setup_aero!`.
- The Makie extension is aero-mode agnostic: wing rendering dispatches on the
  aero mode via `plot_wing_aero!(ax, sys, wing, mode)` (and the per-frame
  `update_wing_aero_plot!`), with methods living in the extension so a custom
  mode draws with full Makie access. Flat-plate wings now render their section
  quads in `plot` and `replay` in the VSM panel style (red mesh, black borders).
- The transform pipeline (`apply_heading!`, `finalize_transforms!`) no longer
  filters wings by aero mode; flat-plate wings now get heading and frame
  finalization like every other wing.
- Internal renames for readability: leading-underscore function names and short
  abbreviations were removed throughout.

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
- The `exposes_aero_input` trait: the `aero_input` connector is detected by name
  on the built subsystem instead.
- The V3-kite-specific analysis code in the Makie extension:
  `compute_ekf_yaw_and_rate`, `compute_ekf_yaw_and_rate_tension`,
  `calculate_cs`, `calc_ref_area`, `middle_le_to_kcu_dir` and their helpers,
  along with the `plot_cs`, `plot_yaw_rate_paper` and `plot_gk_paper` panels
  (the last hardcoded a V3 segment index). V3Kite carries its own copies.
- The `tape_lengths` kwarg of the multi-panel plot and the hardcoded steering
  reconstruction from `segments[87]`: the `plot_us` and `plot_gk` panels now
  read the logged `syslog.steering` directly (so `steering` must be written into
  the `SysState` before `log!` for these panels to show data).
- `set_depower_steering!`, `min_chord_len`, and the
  `SymbolicAWEModel.set_tether_len` field (3-line-kite-specific set-point logic
  with hardcoded tether indices). `calc_side_slip` no longer dispatches on the
  aero mode — it is the same apparent-wind formula for every mode and takes just
  the wing.

### Fixed

- A `DYNAMIC` twist surface without aero sections left
  `twist_surface_aero_moment` unbound and broke `mtkcompile`; the wiring now
  binds the aero component's `twist_moment` for every non-`FIXED`-empty surface.
- Makie extension: `wing isa VSMWing` checks in the panel plotting and the
  log-slot lookup threw at runtime (`VSMWing` is a constructor function since
  the wing merge, not a type).

### Compatibility notes

- Plate logs recorded before the quad logging have a different point count and
  will not replay.

## v0.11.1 06-06-2026

### Added

- `init_stretch_frac` (YAML column and `Tether(...; stretch_frac)` kwarg),
  mutually exclusive with `init_tether_force`: `reinit!` derives the unstretched
  `len` as `len = stretch_frac·stretched`. Setting one input clears the other.
  `init_stretch_frac` must be positive: `<1` pre-stretch, `1` neutral, `>1`
  slack.
- `test_twist_alignment.jl`: under group twist the structural strut
  trailing-edge points stay aligned with the deformed VSM panel trailing edges
  for a `RIGID_DYNAMICS` wing.

### Changed

- `VortexStepMethod` compat raised to `3.3.5`.

### Fixed

- Per-group unrefined moment uses the VSM solver field
  `moment_coeff_unrefined_dist`.
- Body-frame camera tracking across animation frames (Makie ext): `update_cam!`
  with explicit up-vector and `PLOT_BODY_PREV_WING_POS` to eliminate view drift.

## v0.11.0 02-06-2026

### Breaking

- Tether `init_unstretched_length` (YAML) removed; specifying it errors. The
  unstretched rest length is now _derived_: placement is driven by
  `init_stretched_length` (the standoff / placed point geometry, default =
  geometric) and `init_tether_force` (default 0), and
  `len = stretched·(1 − force/unit_stiffness)`.
- `Tether.init_stretched_len`/`init_unstretched_len` are now
  `Union{SimFloat,Nothing}` (`init_unstretched_len` is derived); `Tether` gained
  `init_tether_force`; the positional length constructor arg (now the stretched
  length) is optional. Serialized models must be rebuilt.
- `VSMWing` `origin_idx`/`origin_ref` replaced by `origin::WeightedRefPoints`
  (weighted body-frame origin).
- `update_yaml_from_sys_struct!` and `update_sys_struct_from_yaml!` removed
  (unreliable line-based YAML round-tripping, no longer used).

### Added

- `init_tether_force` (YAML / `Tether(...; tether_force)`, default 0): `reinit!`
  derives every tether's unstretched `len` from the placed stretched length,
  `len = stretched·(1 − force/unit_stiffness)`; force 0 gives zero tension.
- `init!`/`reinit!` `apply_tether_lengths` kwarg to skip placement.
- `WeightedRefPoints(::AbstractString)`; `yaml_parse_origin` for weighted origin
  specs.
- Helpers: `apply_tether_init_forces!`, `tether_unit_stiffness`,
  `tether_anchor_free`, `rigid_point_siblings`, `parse_tether_init`,
  `tether_ordered_point_idxs`, `tether_downstream_idxs`,
  `group_tethers_by_overlap`, `apply_cluster_init_stretched_len!`,
  `_wing_log_pos`; `test_tether_init.jl`.

### Changed

- Tether placement honored only on _root_ tethers (one endpoint on a
  `STATIC`/winch boundary — the fixed anchor, either end); a tether with neither
  endpoint anchored is an error. Tethers sharing a `RIGID_DYNAMICS` wing are
  treated as one cluster (rigid-body connectivity). Multi-root clusters placed
  by the mean displacement of all roots (length + direction), logging `@info`
  (gated on `prn`).
- Wing position stored in dedicated `SysState` slots; reads via
  `update_from_sysstate!` / `_wing_log_pos` / Makie body-frame arrows use
  `wing.pos_w` directly.
- `build_point_to_vsm_point_mapping` takes a `VSMWing`, using body-frame
  closest-point distances.

### Fixed

- Makie zoom/pan world-camera save/restore (no view drift); body-frame zoom
  distance preserved across mode switches.
- `vsm_refine.jl`: RIGID_DYNAMICS wings always keep their aerodynamic panel
  geometry (mesh- or YAML-defined); section rebuilding from structural points is
  now PARTICLE_DYNAMICS-only. The 2plate aero geometry was corrected to match
  its structural points.
- `get_sys_struct_hash` hashes `wing.origin`.

## v0.10.0 30-05-2026

### Changed

- BREAKING: `WingType` constants `QUATERNION` and `REFINE` are now deprecated.
  Use `RIGID_DYNAMICS` and `PARTICLE_DYNAMICS` instead. Deprecated bindings emit
  a warning and will be removed in a future release.
- `DataInterpolations` added as a package dependency (required for `PlateWing`
  polar interpolation).
- `bin/install` now displays an interactive menu to choose Julia version (1.11
  or 1.12) when no version parameter is provided. The currently active Julia
  version is highlighted as the default. Menu is skipped if a version is
  specified via `--version` or `+X.Y` parameters.

### Added

- `PlateWing` and `PlateSurface` types for flat-plate CL/CD lookup aerodynamics.
- `AERO_PLATE` aerodynamics mode — evaluates lift and drag from a polar table
  (CL/CD vs α) via registered symbolic interpolants.
- `create_plate_interpolations(alpha_deg, cl_data, cd_data)` — helper to build
  CL and CD interpolation objects (cubic or linear spline) for use with
  `PlateWing`.
- `examples/kps4_comparison.jl` — comparison of a `PlateWing`-based rigid-body
  kite model against the KiteModels kps4 reference.
- `data/kps4/` — YAML settings and system definition for the kps4 plate model.
- Added missing examples to `examples/menu.jl`: `coupled_linearize`,
  `cosine_steering_trajectory`, `kps4_comparison`, `vsm_linearization`, and
  `sam_tutorial`.

### Fixed

- `init_stretched_len` now works for multi-tether systems. Tethers sharing
  downstream structure are placed to a single effective length (the average of
  several specified values, with a warning), and the initial-positioning BFS no
  longer drags other tethers' ground anchors — it stops at `STATIC` points and
  winch points (which may be `DYNAMIC`).
- `bin/create_sys_image`: fixed a bug that prevented deletion of stale `.so`
  files before rebuilding the system image.
- `AUTHORS.md`: corrected contributor entry.
- `examples/kps4_comparison.jl`: fixed soft-scope ambiguity warning for
  `sys_state` inside the simulation loop.
- Multi-log `plot()` legend labels now render correctly as LaTeX. Added `lbl()`
  helper that places the symbol and suffix inside a single `$...\text{...}$`
  math environment, fixing the literal `$\gamma$ (SymAWE)` display in the
  legend.

### Removed

- `examples/makie_polar_plots.jl` — removed (functionality superseded).

## v0.9.0 20-05-2026

### Changed

- BREAKING: simplified `AERO_LINEARIZED`. ForwardDiff Jacobian over
  `[α, β, ω, θ_groups]` returning wind-axis coefficients
  `[CL, CD, CS, CM, cm_groups]`. Wing fields and accessors renamed `vsm_*` →
  `aero_*`.
- A RIGID_DYNAMICS wing can now have fewer groups than unrefined aero sections
  (one twist DOF drives several sections via a spatial partition). More groups
  than sections errors.
- Bumped `VortexStepMethod` compat to `3.3.0`.
- License changed from MIT/MPL-2.0 to LGPL-3.0-only. All source files updated
  with REUSE-compliant SPDX headers.
- `bin/install` rewritten: unified menu, optional precompile skip, removed
  `bin/update_manifest` and `bin/create_sys_image2`.
- `bin/create_sys_image` updated with improved comments and options.
- `bin/reuse_lint` made more robust with fallbacks for missing tools.
- Safe `atan`/`smooth_normalize` replacements for `asin`/`normalize` in VSM
  equations and linearisation to avoid NaN at edge cases.

### Added

- `examples/vsm_linearization.jl` — plots the VSM linearisation tangents around
  the operating point.
- `test/util.jl` — shared test utilities for allocation checks across all
  integrators.

## v0.8.3 03-05-2026

### Changed

- VSM solver type is taken from VSM settings instead of being hard-coded to
  `NONLIN`.
- At low apparent wind, aero outputs are zeroed instead of warning and skipping.
  Threshold via new `vsm_min_wind` kwarg (default 0.5) on `init!`, `reinit!`,
  `next_step!`.
- Bumped `VortexStepMethod` compat to `3.2.0`.

## v0.8.2 26-04-2026

### Changed

- Updated the default manifest files.

### Added

- `drag_force` field on `Point` — total drag in world frame (point's own
  aerodynamic drag plus its share of connected segment drag). Populated by
  `update_sys_struct!` each timestep.
- Manifest freshness tests in `test_helpers.jl`: verify that no bare
  `Manifest.toml` exists and that `.default` manifests are at least as recent as
  `Project.toml`.
- CI step to copy version-specific `.default` manifest before build, ensuring
  the correct manifest is used per Julia version.
- Drag-related tests in `test_point.jl`, `test_segment.jl`, and `test_wing.jl`.

### Fixed

- Crash with Julia 1.11; `setup_env` updated to fix that.

### Removed

- `plot_recipe.jl` — unused legacy Plots.jl recipe. Visualization is handled by
  `SymbolicAWEModelsMakieExt`.

## v0.8.1 23-04-2026

### Changed

- `SystemStructure.set` field is no longer `const`, allowing change after
  deserialisation.
- Replaced all `@unpack` macro usage with Julia's native destructuring syntax
  `(; a, b) = x`.

### Fixed

- Fixed JETLS warnings across multiple source files.
- `bin/install` now copies `.JETLSConfig.toml.default` to `.JETLSConfig.toml` if
  it does not exist, and warns when an existing config differs from the default.
- `bin/install` warning messages now use colored output for visibility.

## v0.8.0 18-04-2026

### Changed

- BREAKING: `SegmentType` positional argument removed from `Segment`
  constructor. Use `unit_stiffness`, `unit_damping`, `diameter_mm` kwargs or a
  YAML material instead. The `SegmentType` enum is kept temporarily to produce a
  helpful deprecation error.
- BREAKING: `winch_point` moved from `Tether` to `Winch`. Pass `winch_point` as
  a keyword to the `Winch` constructor instead.
- BREAKING: Heading calculation changed from wind-perpendicular projection to
  tangential sphere frame. `calc_heading(R_b_to_w,
  wind_norm)` →
  `calc_heading(R_b_to_w, wing_pos)`. `get_heading_components()` removed.
  `solve_heading_rotation` takes `wing_pos` instead of `k, wind_norm`.
- BREAKING: `Tether` struct fields restructured — `winch_point_idx/ref` removed,
  new fields: `start_point_idx/ref`, `end_point_idx/ref`, `n_segments`,
  `unit_stiffness`, `unit_damping`, `diameter`.
- BREAKING: `create_tether()` utility returns a 5-tuple (added
  `ground_point_idx`) and no longer takes a `SegmentType` argument.
- BREAKING: YAML segment format no longer has a `type` column. Existing YAML
  files with a `type` column in segments will raise an error.
- BREAKING: `tether_len` moved from `Winch` to `Tether`. Each tether now owns
  its length as an ODE state variable. Winch-connected tethers evolve via
  `D(tether_len) = winch_vel`; winch-less tethers have constant length
  (`D(tether_len) = 0`).
- BREAKING: `tether_vel` renamed to `winch_vel` and remains on `Winch`.
  `tether_acc` renamed to `winch_acc` in the generated equations.
- BREAKING: `SimpleLinModelWithAttributes` removed. The `simple_lin_model` field
  is no longer part of `SymbolicAWEModel`. `simple_linearize!` is no longer
  exported.
- BREAKING: `sim_oscillate!` and `sim_turn!` removed. Use `sim!` with a custom
  `set_values` matrix instead.
- BREAKING: `update_aero_yaml_from_struc_yaml!` no longer exported.
- BREAKING: `set` field removed from `SymbolicAWEModel`. Settings are now read
  from `sam.sys_struct.set`. The `set_set` setter was removed from
  `ProbWithAttributes` and `LinProbWithAttributes`.
- BREAKING: `get_struct_state` removed from `ProbWithAttributes`.
- Wind equations now use `get_wind_vec` internally instead of separate
  `get_v_wind`, `get_upwind_dir`, and `get_wind_elevation` accessors. Not
  breaking: KiteUtils `Settings` syncs `wind_vec` from
  `v_wind`/`upwind_dir`/`upwind_elevation` automatically when
  `use_wind_vec=false` (the default).
- Tethers no longer require a connected winch. Winch-less tethers use constant
  `l0` from segment properties.
- `compression_frac` description clarified: "Compressive/tensile stiffness ratio
  (0-1). 0 = no compression stiffness."
- `init!`, `next_step!`, `update_sys_state!` are no longer exported and must be
  imported from `KiteUtils`.
- `sim!` now requires `y_op` keyword argument when `lin_model` is provided
  (previously obtained from the removed simple lin model).
- `SerializedModel` type parameters tightened for `defaults` and `guesses`
  fields.
- fixed most `JETLS` warnings for improved robustness and performance.
- Package version is now included in `.bin` cache filenames, so upgrading the
  package automatically invalidates stale cached models.
- the script `bin/run_julia` was updated to work also with Julia 1.12.6

### Added

- Route 2 tether auto-generation:
  `Tether(name; start_point,
  end_point, n_segments)` automatically creates
  intermediate points and segments, evenly spaced between endpoints. YAML
  format: `headers: [name, start_point, end_point, n_segments, ...]`.
- Route 1 tethers auto-detect `start_point_idx` and `end_point_idx` from the
  first/last segment endpoints.
- Comprehensive docstrings on all `Point`, `Group`, `Segment`, `Pulley`,
  `Tether`, `Winch`, and `Transform` struct fields.
- `WeightedRefPoints` exported for weighted reference point support.
- `init!` keyword `reinit_sys` to optionally skip system structure
  reinitialization.
- New tests: "Route 2 auto-generated tether" and "Tether without winch" in
  `test_tether_winch.jl`.
- New test file `test_tether_init.jl` for tether initialization.
- New test file `test_yaml_weighted_ref.jl` for weighted reference point YAML
  loading.
- Airbag pressurized membrane simulation example (`examples/airbag.jl`).
- the script `bin/install`. Use it after installation from git.
- the script `bin/create_sys_image`. Improves time for first run by a factor of
  3-5.
- the scripts `bin/install_jetls` and `bin/jetls` to install and run `JETLS.jl`,
  a static code checker for Julia.
- Developer documentation improvements (troubleshooting section for segfault
  issues, updated docs to use GLMakie).

### Fixed

- YAML `calculate_derived_properties!` no longer requires `l0` to compute
  `unit_stiffness` from material properties (needed for Route 2 tethers).
- YAML `update_yaml_from_sys_struct!` regex updated for the new segment format
  (no `type` column).
- YAML weighted reference point loading fixed (broken deserialization of
  weighted refs).
- Heading calculation uses tangential sphere frame, fixing drift issues with the
  old wind-perpendicular projection.
- Unknown solver string (e.g. `DFBDF` from default KiteUtils settings) no longer
  throws an error — a warning is emitted and the solver falls back to `FBDF`.
- README code examples now include the required
  `SymbolicAWEModels.init_module(; force=false)` call so they work correctly on
  a fresh install.
- README pendulum example also calls `set_data_path("data/base")` before loading
  `Settings`.

### Removed

- `SimpleLinModelWithAttributes` struct and `simple_linearize!`.
- `sim_oscillate!` and `sim_turn!` simulation functions.
- `getstate` and `setstate!` functions from `linearize.jl`.
- `upwind_dir` helper function (replaced by `wind_vec`).
- Branch-specific system images: `bin/create_sys_image` and `bin/run_julia` no
  longer embed the git branch name in the `.so` filename. A single
  `kps-image-<julia_major>.so` is used instead.

### Tests

- README pendulum example and README 2-plate kite example are now executed in
  `test/setup_integration.jl`.

## v0.7.2 18-03-2026

### Added

- `speed_controlled` field on `Winch` — when `true`, tether velocity is
  prescribed externally (`D(tether_vel) = 0`) while length still tracks
  velocity.
- Multi-system `record()` for recording side-by-side SysLog animations to video
  (MP4/GIF/MKV/WebM).
- Makie extension test suite (`test_makie_extension.jl`) covering multi-system
  plot, record, and replay.
- Zenodo metadata (`.zenodo.json`) and `CITATION.cff` for citing the package.
- CI: GLMakie tests on Linux via `xvfb-run`, Julia 1.12 test matrix.

### Fixed

- `reposition!()` now uses the analytical `solve_heading_rotation` for
  wind-relative heading, consistent with `reinit!`. Previously heading was
  applied as a relative delta, causing drift.
- `reposition!()` correctly updates PARTICLE_DYNAMICS wings by recalculating
  `R_b_to_w` and `pos_b` from structural points.
- Multi-system `plot()` now passes vector-typed segment colors, fixing a crash
  when `setup_segment_hover_events!` assigned `Vector{RGBA}`.
- `init!()` validates that `SystemStructure` uses `VSMWing` type before equation
  generation.
- `sim_reposition!()` passes absolute heading to the transform instead of
  subtracting the current wing heading.
- Typo fixes in README and documentation ("ODE solver" → "ODE problem").

### Changed

- `sam_tutorial.jl` example updated: adds WING-type points and uses
  `VSMSettings` with `data_prefix=false`.
- Examples updated to pass `data_prefix=false` to `VSMSettings`.
- 2plate_kite aero geometry TE z-coordinates adjusted.
- `settings.yaml` now includes `sample_freq` field.

## v0.7.1 27-02-2026

### Added

- `update_sys_struct_from_yaml!()` — update a `SystemStructure` in-place from a
  modified YAML file (point `pos_cad` and segment `l0`).
- `segment_cad_length()` and `autocalc_tether_len()` shared helpers, replacing
  duplicated code in the constructor, `reinit!`, and YAML loader.

### Fixed

- `SystemStructure` constructor auto-calculates `winch.tether_len` from all
  connected tethers (was only using the first).

## v0.7.0 DD-02-2026

### Changed

- BREAKING: Julia version requirement raised from 1.10 to 1.11, 1.12.
- `reinit!()` uses a unified code path for all wing types, calling
  `match_aero_sections_to_structure!` and `compute_spatial_group_mapping!`
  during VSM rebuild.
- `test_bench.jl` refactored from ad-hoc benchmarks into a proper `@testset`
  suite with `setup_bench_sam()` helper.
- Added `[workspace]` configuration in `Project.toml` for docs, examples,
  scripts, and test sub-projects.
- Manifest files renamed to `.default` suffix and gitignored.

### Added

- Asymmetric aero/structural section counts: aerodynamic and structural meshes
  can now have different numbers of sections. When counts differ,
  `match_aero_sections_to_structure!()` rebuilds unrefined sections from
  structural LE/TE positions while `use_prior_polar=true` preserves existing
  refined panel polars. Opt-in via `use_prior_polar=true` on the
  VortexStepMethod wing.
- `identify_wing_segments()` — identifies LE/TE pairs from groups (preferred) or
  via a consecutive-pair heuristic.
- `compute_spatial_group_mapping!()` — maps groups to VSM sections by spatial
  proximity, supporting n_groups != n_aero_sections.
- PARTICLE_DYNAMICS wings can now have groups (used for LE/TE pair
  identification).
- RIGID_DYNAMICS wings can now have `wing_segments` for structural geometry
  locking.
- YAML loader fallback LE/TE detection in `update_aero_yaml_from_struc_yaml!()`
  when no groups are defined (consecutive-pair heuristic with x-coordinate
  check).
- `test_match_aero_sections.jl` — tests geometry matching and polar
  interpolation for both PARTICLE_DYNAMICS and RIGID_DYNAMICS wings, including
  mismatched section counts.
- Helper scripts: `bin/install` (environment setup, Julia version detection) and
  `bin/run_julia` (launcher with system image support).

## v0.6.1 23-02-2026

### Fixed

- Disable VSM auto-sorting of sections (`sort_sections=false`) in all
  VortexStepMethod calls. Auto-sorting silently broke the correspondence between
  VSM sections and structural point indices / group mappings.

## v0.6.0 21-02-2026

### Changed

- Component constructors (`Point`, `Segment`, `Wing`, `Winch`, `Transform`) now
  accept a symbolic `name` (Symbol) as the first argument in addition to numeric
  indices. Numeric `idx` values still work. Use e.g.
  `Point(:kcu, pos, DYNAMIC)`.
- BREAKING: `Segment` constructor takes separate `point_i`, `point_j` arguments
  instead of a `point_idxs` vector.
- BREAKING: Rotation matrix fields renamed from `R_a_b` to `R_a_to_b` throughout
  (e.g. `wing.R_b_w` → `wing.R_b_to_w`).
- BREAKING: `ControlPlotsExt` package extension removed. Visualization is now
  handled entirely by `SymbolicAWEModelsMakieExt`.
- BREAKING: Predefined model factory functions removed (`create_ram_sys_struct`,
  `create_simple_ram_sys_struct`). Build models using component constructors or
  YAML instead.
- BREAKING: Ram air kite and V3 kite models moved to dedicated packages
  ([RamAirKite.jl](https://github.com/OpenSourceAWE/RamAirKite.jl),
  [V3Kite.jl](https://github.com/OpenSourceAWE/V3Kite.jl)). Their data
  directories are removed from this package.
- `src/system_structure.jl` split into modular files under
  `src/system_structure/` (types, core, utilities, transforms, wing,
  named_collection).
- `src/generate_system.jl` split into 13 focused modules under
  `src/generate_system/` (point_eqs, segment_eqs, wing_eqs, group_eqs,
  winch_eqs, pulley_eqs, tether_eqs, scalar_eqs, vsm_eqs, accessors, helpers,
  create_sys).
- Makie extension significantly overhauled with new plotting functions.
- Test suite completely rewritten. The old tests (`test_simulation`,
  `test_linearization`, `test_initialization`, `test_sam`, etc.) tested the full
  assembled kite model as a black box, making failures hard to diagnose. The new
  tests isolate each component with minimal models built from constructors,
  verifying physics against analytical solutions:
  - `test_point` — gravity free-fall, damping, quasi-static equilibrium
  - `test_segment` — spring-damper forces, stiffness, drag
  - `test_wing` — RIGID_DYNAMICS and PARTICLE_DYNAMICS wing construction, VSM
    coupling
  - `test_wing_dynamics` — rigid body torque response, precession, angular
    momentum conservation
  - `test_tether_winch` — reel-out dynamics, Coulomb and viscous friction,
    terminal velocity
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

- `NamedCollection` indexing — components support symbolic names (e.g.
  `sys.points[:kcu]`, `sys.segments[:bridle_1]`). `SystemStructure` resolves all
  symbolic references to numeric indices automatically via
  `assign_indices_and_resolve!()`.
- `WingType` enum (`RIGID_DYNAMICS`, `PARTICLE_DYNAMICS`) for explicit wing type
  selection. BREAKING: these names replace the previous `QUATERNION` and
  `REFINE` wing types. Update YAML configs from `type: QUATERNION` /
  `type: REFINE` to `dynamics_type: RIGID_DYNAMICS` /
  `dynamics_type: PARTICLE_DYNAMICS`, and rename the wing `type` field to
  `dynamics_type`. Update any code using the old exported constants.
  `PARTICLE_DYNAMICS` applies per-panel forces directly to structural points for
  higher fidelity aeroelastic coupling.
- `AeroMode` enum (`AERO_NONE`, `AERO_DIRECT`, `AERO_LINEARIZED`) for build-time
  control over aerodynamic computation strategy.
- YAML-based model definition via `load_sys_struct_from_yaml()`,
  `update_yaml_from_sys_struct!()`, and `update_aero_yaml_from_struc_yaml!()`.
- PARTICLE_DYNAMICS wing support (`src/vsm_refine.jl`) — structural deformation
  coupled directly to VSM panel geometry with moment-preserving force
  distribution.
- Principal vs body frame separation for RIGID_DYNAMICS wings. Principal frame
  (diagonal inertia) used for Euler equations, body frame (from reference
  points) used for output and VSM coupling.
- Auto-group generation for RIGID_DYNAMICS wings when groups are not explicitly
  provided.
- `record()` for saving simulation replays to MP4.
- `plot_sphere_trajectory`, `plot_body_frame`, `plot_aoa` plotting functions.
- `update_segment_forces!`, `set_world_frame_damping`, `set_body_frame_damping`,
  `segment_stretch_stats` utility functions.
- New examples: `hanging_mass`, `catenary_line`, `saddle_form`,
  `coupled_2plate_kite`, `coupled_realtime_visualization`, `coupled_linearize`,
  `coupled_simple_lin_model`, `coupled_tether_deflection`, `heading_gate`,
  `cosine_steering_trajectory`, `makie_polar_plots`, `static_load_2plate_kite`.
- Benchmark test (`test_bench.jl`) for performance tracking.

### Removed

- `predefined_structures.jl` and factory functions (`create_ram_sys_struct`,
  `create_simple_ram_sys_struct`, `create_tether_sys_struct`,
  `copy_to_simple!`).
- Ram air kite data files, LEI kite directory, `data/kite.obj`.
- Old examples: `ram_air_kite`, `lin_ram_model`, `simple_lin_model`,
  `lin_simple_tuned_model`, `simple_tuned_model`, `realtime_visualization`,
  `reposition`, `tether_props`.
- `SymbolicAWEModelsControlPlotsExt` package extension.
- `src/precompile.jl`.

## v0.5.0 25-08-2024

### Removed

- BREAKING: the Winch struct doesn't have a model field anymore. Instead, all
  equations are symbolic, and the WinchModels dependency is removed.

### Added

- The function `calc_steady_torque` calculates the torque that will result in
  zero acceleration.

## v0.4.2 24-08-2024

### Fixed

- Don't write protect manifest

## v0.4.1 13-08-2025

### Fixed

- Update Artifacts.toml.default

## v0.4.0 13-08-2025

### Added

- Structs with attributes for better serialization and code structure
  (`SimpleLinModelWithAttributes`, `ProbWithAttributes`,
  `LinProbWithAttributes`, `ControlFuncWithAttributes`).
- `plot_force` option to the plot recipe.
- `model_management.jl` file to better organize the code.

### Changed

- BREAKING: `init_module` function to simplify project setup, replacing
  `install_examples`, `copy_examples`, `copy_bin` and `copy_model_settings`.
- Major refactoring of the `SymbolicAWEModel` and its initialization process.
  The `SerializedModel` struct is now much simpler and more robust.
- The `run_julia` script is now much more powerful, with argument parsing for
  `--copy-manifest` and `--precompile`.
- The precompilation process now uses artifacts instead of downloading files
  directly.

### Fixed

- URLs in `Artifacts.toml.default`.
- Cross-correlation analysis in tests.

### Removed

- `data/kite.obj` file.
- `copy_examples`, `copy_bin`, `copy_model_settings`, `install_examples`
  functions.

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

- Breaking: sim!, sim_oscillate! and sim_turn! return a tuple (sl, lin_sl)
  instead of just a sl

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
- Adds `copy_to_simple!` function, which copies the ram model state to the
  simple model state, uses the tether model to find the equivalent 1-segment
  spring properties of the tether
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
