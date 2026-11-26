### Added

- `find_equilibrium!(sam)` moves a `KernelBackend` model to the equilibrium of its structure at its current pose, under the loads it carries now with the aerodynamics held at their last solve, and returns whether it converged, the largest load left unbalanced and the moment about the base that holding the pose takes.
