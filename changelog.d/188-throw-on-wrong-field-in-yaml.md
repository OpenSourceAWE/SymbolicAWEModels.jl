### Fixed
- A structural YAML naming a block or a row field the loader does not read, or
  giving a row more values than the table has headers, stops the load and says
  what the block accepts. The field was silently dropped and the row silently
  skipped before, so a typo left the component on its default. The `idx` column
  the tutorial showed on points, segments, pulleys, tethers, winches and
  transforms was one of those and now errors: rename it to `name`.
- `anchor_b` on a point row reaches the point. A `BODY_STATIC` point anchored to
  a body took its offset from the geometry even when its row gave one, because
  the loader never passed the column on.
