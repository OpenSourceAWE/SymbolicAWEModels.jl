### Added
- `using SymbolicAWEModels` also loads KiteGeometry and exports its names, except where SymbolicAWEModels has its own of the same name, such as `Point` or `structure_document`.
- KiteGeometry's `load_authoring` and `load_structure` read SymbolicAWEModels' own point, station, wing and pulley columns, such as `body_frame_damping`, `moment_frac` and `aero_z_offset`, instead of refusing them.
- The 2-plate kite's structural YAMLs give the wing's canopy faces.
