### Added

- `find_quasi_steady!(sam)` finds a kite's quasi-steady state at the pose its transforms hold, on the `KernelBackend`: it relaxes the kite under damping measured against the rigid rotation fitted to its own motion and judges it with all damping off. It returns whether it converged and the residual force and moment it reached, so a pose with no attracting equilibrium is a result rather than a failure. `quasi_steady_residual(sam)` is that judgement at the current state, and `quasi_steady_log(sam)` a one-frame `SysLog` of it.
