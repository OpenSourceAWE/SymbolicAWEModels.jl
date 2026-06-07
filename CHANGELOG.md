# CHANGELOG

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
