### Added

- `system.yaml` is the project file: besides `sim_settings:` it names a model's
  `structural_geometry:`, `aero_geometry:` and `vsm_settings:`, the layout
  `V3Kite.jl` already uses. `project_file(entry)` resolves one against the data
  path and throws an `ArgumentError` naming the entry where the project names
  none, so examples ask the project for a model's files instead of spelling their
  names. A project that names only `sim_settings:` keeps working.

### Changed

- A wing's aero geometry is the project file's `aero_geometry:` where the project
  names one, so the file is named in one place instead of also being the
  `geometry_file:` of every wing in `vsm_settings.yaml`. A `VSMSettings` built
  outside the project still flies its own `geometry_file:`.
- The `kite.struc_geometry_path` and `kite.aero_geometry_path` fields of
  `settings.yaml` are gone from the shipped models and from the documented
  schema. Nothing ever read them; `structural_geometry:` and `aero_geometry:` in
  the project file take their place.
