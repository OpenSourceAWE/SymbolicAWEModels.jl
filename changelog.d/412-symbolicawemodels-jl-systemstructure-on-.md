### Added
- `using SymbolicAWEModels` also loads KiteGeometry and exports its names, except where SymbolicAWEModels has its own of the same name, such as `Point` or `structure_document`.
- KiteGeometry's `load_authoring` and `load_structure` read SymbolicAWEModels' own point, station, wing, pulley and body columns, such as a point's `type`, `body_frame_damping` and `drag_coeff`, a station's `moment_frac` or a wing's `aero_z_offset`, instead of refusing them. SymbolicAWEModels needs KiteGeometry 0.2.
- The 2-plate kite's structural YAMLs give the wing's canopy faces.
