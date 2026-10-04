### Changed
- BREAKING: the positional `VSMWing(name, vsm_aero, vsm_wing, vsm_solver, stations,
  R_b_to_c, pos_cad)` and the `Wing` method taking the same arguments are removed;
  build a VSM wing with `VSMWing(name, set, stations, vsm_set; R_b_to_c, pos_cad)`.
  A wing with an `.obj` mesh that was built positionally used point-mass inertia;
  the keyword constructor gives it the mesh inertia tensor, so its dynamics change.
