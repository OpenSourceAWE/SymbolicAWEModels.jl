### Added

- `segment_color` in `plot`, `plot!` and `replay` also takes a vector with one colour per
  segment, or a function `segment -> colour`.
- `segment_role(sys, segment)` names a segment `:winched_tether`, `:unwinched_tether`,
  `:wing` or `:free`, to colour a plot by role.
