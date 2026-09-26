### Changed
- BREAKING: building a VSM wing errors, naming the wing, when its structure and its
  aerodynamic geometry are not in one CAD frame: a station node or the mesh COM
  outside the sections' bounding box grown by 10 % of its largest side, or the span
  from `y_ref_points` more than 5° from the sections' span.
