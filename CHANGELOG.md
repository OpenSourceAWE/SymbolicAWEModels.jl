<!--
SPDX-FileCopyrightText: 2025 Uwe Fechner, Bart van de Lint
SPDX-License-Identifier: MPL-2.0
-->

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
