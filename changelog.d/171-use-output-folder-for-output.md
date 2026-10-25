### Added

- `get_output_path` and `set_output_path`, re-exported from KiteUtils, name the folder simulation results go to: loading SymbolicAWEModels sets it to `output` in the working directory, created on first use, and `set_output_path` names another.

### Changed

- Simulation results go to the output folder instead of the working directory or the data folder: the logs `sim!` and `sim_reposition!` save, the logs of every example, a `record` to a relative filename, and the replay viewer's Save button screenshot (#227). `save_log` and `load_log` use the output folder when no `path` is passed, so a log kept in the data folder is read with `path=get_data_path()`.
