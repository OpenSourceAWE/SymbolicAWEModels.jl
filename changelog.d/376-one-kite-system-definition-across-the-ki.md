### Added
- `place!(sys_struct)` puts a structure in the world by its transforms and makes
  where it lands its initial pose. The `SystemStructure` constructor runs it.
- `init_sys_struct!(sys_struct, set; remake_vsm)` brings a placed structure to the
  start of a run: mass properties and rest geometry derived from the initial pose
  (`init_derived_properties!`), principal-frame state, aero engines and wind.
  `init!` runs it, so a run restarted from a logged state keeps its rest shape.

### Changed
- Requires VortexStepMethod 6.1, whose panel mapping breaks a tie between two
  sections outboard, so a mirror-symmetric wing maps mirror panels to mirror sections.
- BREAKING: `SystemStructure` carries `tubes::Vector{Tube}` in place of
  `elastic_joints` and `timoshenko_joints`. A `Tube` names its two bodies, diameter,
  pressure and law, and its `model` is a `TimoshenkoTube` or an `ElasticTube`, whose
  rigidities come from the law unless given. The authoring YAML's two joint tables
  are one `tubes` table, with `model` and its parameters as extra columns.
- BREAKING: `init!` no longer places the structure; it starts from `sam.sys_struct`
  as it stands. After changing a transform, call `place!(sam.sys_struct)` before
  `init!`. `reinit!(sys_struct, set)` is `place!(sys_struct)`, and `init!` drops
  `reinit_sys`, `reset_vel`, `ignore_l0` and `apply_tether_lengths`. `init!` always
  sets the wind from `set.wind_vec`; `remake_vsm=false` keeps a hand-edited aero.
- BREAKING: `Point` and `Body` hold their initial pose in the world frame,
  `pos_ENU` and (bodies) the quaternion `Q_KA_to_ENU`, in place of the CAD frame
  `pos_cad`, `R_b_to_c` and `R_p_to_c`; the `Wing` constructor takes
  `Q_KA_to_ENU, pos_ENU`, and `VSMWing` the keywords of those names. `reset_to_cad!` is `reset_to_initial_pose!`. The
  authoring YAML keeps its `pos_cad` column.
- BREAKING: a structure document is written and read against awesIO
  `structure_schema.yml` 1.0.0: the system in its initial pose, every position in
  ENU, bodies with `Q_KA_to_ENU` and an `extra_inertia_KA` tensor about their own
  mass centre, a `units` row per table, a `wings` block and one `tubes` table. A
  0.1.0 document is refused. Reading a document with wings is not supported yet
  (#396).

### Fixed
- A structure document is read in the initial pose it gives, without `place!`, so a
  pre-tensioned tether keeps the `l0` the document holds (`SystemStructure(...;
  placed=true)`).
- Placing a structure hung on several tethers with `init_stretched_length` moves it
  straight away from their mean anchor until their mean length is the mean target,
  so tethers that already have their lengths stay and placing again moves nothing.
- `init_pulley_lengths!` sets the two pulley segments' `l0` to the split it writes
  into the pulley, as the model uses them, so the mass of their points matches it.
