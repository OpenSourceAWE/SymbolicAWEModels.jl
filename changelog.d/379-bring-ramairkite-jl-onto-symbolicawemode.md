### Fixed
- A VSM wing whose structure and aerodynamic geometry are in different CAD frames
  now errors at load, naming the wing, instead of loading silently (#248). Both
  were always meant to share one CAD frame, as the coordinate-frames page says;
  when they do not, the sections move into the body frame fitted to the structure
  and the aero forces act at the wrong place on the wing. A wing is refused when a
  station node or the mesh COM lies outside the sections' bounding box grown by
  10 % of its largest side, or when the span from `y_ref_points` is more than 10°
  from the sections' span.
