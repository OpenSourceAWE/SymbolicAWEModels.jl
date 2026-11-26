### Added

- `find_equilibrium!(sam)` moves a `KernelBackend` model to its quasi-steady state at rest at its current pose: it solves the aerodynamics on the current shape, settles the structure under them and repeats, each round starting from a quasi-Newton prediction (IQN-ILS), until the forces balance. It returns whether it converged, the rounds and Newton steps it took, the largest load left unbalanced and the moment about the base that holding the pose takes.
