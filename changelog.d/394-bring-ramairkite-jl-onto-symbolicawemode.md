### Changed
- BREAKING: the positional `VSMWing(name, vsm_aero, vsm_wing, vsm_solver, stations,
  R_b_to_c, pos_cad)` and the `Wing` method taking the same arguments are removed;
  build a VSM wing with `VSMWing(name, set, stations, vsm_set; R_b_to_c, pos_cad)`,
  which also seeds the mesh inertia.
