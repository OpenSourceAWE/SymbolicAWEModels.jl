### Added

- `TurbulentWind`, a wind mode in which every point and wing flies in the turbulent
  wind field of AtmosphericModels at its own position, sampled at the start of every
  step and held through it. `set.use_turbulence` scales the turbulence; the mean wind
  keeps `set.upwind_dir` and `set.upwind_elevation`. The settings files carry the
  field extent `environment.grid` it needs.
