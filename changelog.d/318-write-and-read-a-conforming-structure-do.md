### Added

- A `SystemStructure` writes and reads as a *structure document* conforming to
  awesIO's `structure_schema.yml`: `structure_document` renders one, and
  `save_structure_document` / `load_structure_document` write and read it as YAML
  or JSON by file extension, both through the one `sys_struct_from_document`
  loader. The document is the resolved structure — points, segments, stations,
  pulleys, tethers, winches, bodies and joints, every reference by name — without
  the transforms that place it in the world.
