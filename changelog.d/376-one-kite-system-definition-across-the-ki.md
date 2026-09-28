### Added
- `place!(sys_struct)` puts a structure in the world from its CAD geometry by its
  transforms. The `SystemStructure` constructor ends with it.
- `init_sys_struct!(sys_struct, set; remake_vsm)` brings a placed structure to the
  start of a run: principal-frame state, aero engines and wind. `init!` runs it.

### Changed
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
- BREAKING: a structure document is written and read against awesIO
  `structure_schema.yml` 1.0.0: the system where it is placed, every position in
  ENU, bodies with `Q_KA_to_ENU` and an `extra_inertia_KA` tensor about their own
  mass centre, a `units` row per table, a `wings` block and one `tubes` table. A
  0.1.0 document is refused. Reading a document with wings is not supported yet
  (#396).
