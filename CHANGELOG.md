<!--
SPDX-FileCopyrightText: 2025 Uwe Fechner, Bart van de Lint
SPDX-License-Identifier: MPL-2.0
-->

# v0.7.1 27-02-2026

## Added
- `update_sys_struct_from_yaml!()` — update a `SystemStructure` in-place
  from a modified YAML file (point `pos_cad` and segment `l0`).
- `segment_cad_length()` and `autocalc_tether_len()` shared helpers,
  replacing duplicated code in the constructor, `reinit!`, and YAML loader.

## Fixed
- `SystemStructure` constructor auto-calculates `winch.tether_len` from
  all connected tethers (was only using the first).

# v0.7.0 DD-02-2026

## Changed
- BREAKING: Julia version requirement raised from 1.10 to 1.11, 1.12.
- `reinit!()` uses a unified code path for all wing types, calling
  `match_aero_sections_to_structure!` and
  `compute_spatial_group_mapping!` during VSM rebuild.
- `test_bench.jl` refactored from ad-hoc benchmarks into a proper
  `@testset` suite with `setup_bench_sam()` helper.
- Added `[workspace]` configuration in `Project.toml` for docs, examples,
  scripts, and test sub-projects.
- Manifest files renamed to `.default` suffix and gitignored.

## Added
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
- REFINE wings can now have groups (used for LE/TE pair identification).
- QUATERNION wings can now have `wing_segments` for structural geometry
  locking.
- YAML loader fallback LE/TE detection in
  `update_aero_yaml_from_struc_yaml!()` when no groups are defined
  (consecutive-pair heuristic with x-coordinate check).
- `test_match_aero_sections.jl` — tests geometry matching and polar
  interpolation for both REFINE and QUATERNION wings, including
  mismatched section counts.
- Helper scripts: `bin/install` (environment setup, Julia version detection)
  and `bin/run_julia` (launcher with system image support).

# v0.6.1 23-02-2026

## Fixed
- Disable VSM auto-sorting of sections (`sort_sections=false`) in all
  VortexStepMethod calls. Auto-sorting silently broke the correspondence
  between VSM sections and structural point indices / group mappings.

# v0.6.0 21-02-2026

## Changed
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
  - `test_wing` — QUATERNION and REFINE wing construction, VSM coupling
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

## Added
- `NamedCollection` indexing — components support symbolic names
  (e.g. `sys.points[:kcu]`, `sys.segments[:bridle_1]`).
  `SystemStructure` resolves all symbolic references to numeric indices
  automatically via `assign_indices_and_resolve!()`.
- `WingType` enum (`QUATERNION`, `REFINE`) for explicit wing type
  selection. `REFINE` applies per-panel forces directly to structural
  points for higher fidelity aeroelastic coupling.
- `AeroMode` enum (`AERO_NONE`, `AERO_DIRECT`, `AERO_LINEARIZED`) for
  build-time control over aerodynamic computation strategy.
- YAML-based model definition via `load_sys_struct_from_yaml()`,
  `update_yaml_from_sys_struct!()`, and
  `update_aero_yaml_from_struc_yaml!()`.
- REFINE wing support (`src/vsm_refine.jl`) — structural deformation
  coupled directly to VSM panel geometry with moment-preserving force
  distribution.
- Principal vs body frame separation for QUATERNION wings. Principal
  frame (diagonal inertia) used for Euler equations, body frame (from
  reference points) used for output and VSM coupling.
- Auto-group generation for QUATERNION wings when groups are not
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

## Removed
- `predefined_structures.jl` and factory functions
  (`create_ram_sys_struct`, `create_simple_ram_sys_struct`,
  `create_tether_sys_struct`, `copy_to_simple!`).
- Ram air kite data files, LEI kite directory, `data/kite.obj`.
- Old examples: `ram_air_kite`, `lin_ram_model`, `simple_lin_model`,
  `lin_simple_tuned_model`, `simple_tuned_model`,
  `realtime_visualization`, `reposition`, `tether_props`.
- `SymbolicAWEModelsControlPlotsExt` package extension.
- `src/precompile.jl`.

# v0.5.0 25-08-2024
## Removed
- BREAKING: the Winch struct doesn't have a model field anymore. Instead, all equations are symbolic, and the WinchModels dependency is removed.
## Added
- The function `calc_steady_torque` calculates the torque that will result in zero acceleration.

# v0.4.2 24-08-2024
## Fixed
- Don't write protect manifest

# v0.4.1 13-08-2025
## Fixed
- Update Artifacts.toml.default

# v0.4.0 13-08-2025
## Added
- Structs with attributes for better serialization and code structure (`SimpleLinModelWithAttributes`, `ProbWithAttributes`, `LinProbWithAttributes`, `ControlFuncWithAttributes`).
- `plot_force` option to the plot recipe.
- `model_management.jl` file to better organize the code.
## Changed
- BREAKING: `init_module` function to simplify project setup, replacing `install_examples`, `copy_examples`, `copy_bin` and `copy_model_settings`.
- Major refactoring of the `SymbolicAWEModel` and its initialization process. The `SerializedModel` struct is now much simpler and more robust.
- The `run_julia` script is now much more powerful, with argument parsing for `--copy-manifest` and `--precompile`.
- The precompilation process now uses artifacts instead of downloading files directly.
## Fixed
- URLs in `Artifacts.toml.default`.
- Cross-correlation analysis in tests.
## Removed
- `data/kite.obj` file.
- `copy_examples`, `copy_bin`, `copy_model_settings`, `install_examples` functions.

# v0.3.3 07-08-2025
## Fixed
- Fix non-persistent state bug with `calc_tether_props`

# v0.3.2 07-08-2025
## Fixed
- Fix documentation for sim_oscillate!

# v0.3.1 06-08-2025
## Fixed
- Fix examples and menu

# v0.3.0 06-08-2025
## Changed
- Breaking: sim!, sim_oscillate! and sim_turn! return a tuple (sl, lin_sl) instead of just a sl
## Fixed
- Restrict LinearSolve version to `<3.25.0`
- Fixed `linearize!(sam)` to get updated when the state gets updated
## Added
- Added `lin_simple_tuned_model.jl` example

# v0.2.1 01-08-2025
## Fixed
- Import Pkg

# v0.2.0 01-08-2025
## Added
- Adds simple model and tether model
- Adds `copy_to_simple!` function, which copies the ram model state to the simple model state, uses the tether model to find the equivalent 1-segment spring properties of the tether
- Adds open-loop sim functions `sim!`, `sim_oscillate!`, `sim_turn!`
- Adds plotting function `plot(sys_struct::SystemStructure, sys_log::SysLog)`
- Adds documentation
- Adds new updated tests: test/test_sam.jl
## Fixed
- Fixes documentation
- Fixes the bug where the kite could not have negative position
## Changed
- Improved precompilation
- Breaking: `Segment` constructor has different arguments
## Removed
- Removed `.bin` files from git, will be added as release artifacts

# v0.1.3 18-07-2025
## Changed
- Add interface keyword arguments to `init!`

# v0.1.2 13-07-2025
## Changed
- Update VortexStepMethod.jl

# v0.1.1 13-07-2025
## Added
- Added a simple linearized model
## Changed
- Improved the reinitialization using scalar settings values
- Update KiteUtils and AtmosphericModels

# v0.1.0
- Moved the SymbolicAWEModel from KiteModels.jl to SymbolicAWEModels.jl
