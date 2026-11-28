### Changed

- The model cache key carries the tether base point each body damps against, so every cached model with a body recompiles once.

### Fixed

- Angular damping and body-frame damping no longer brake a kite's flight over the sphere about its tether's base point, its turns about the tether, or its reel-out. Both damp motion relative to the kite's small-earth frame: `z` from the base point to the kite, `y` the kite's span made orthogonal to `z`, stretching radially with the tether. `angular_damping` damps a body's spin relative to that frame's rotation, so a yaw coefficient no longer acts on a turn about the tether. `body_frame_damping` of a wing's nodes and of a body with a parent wing damps their velocity relative to the frame's velocity at their position. A kite flying across the sphere therefore trims differently than before. A body that no transform moves damps its absolute spin as before.
- On the kernel backend a fitted (particle) wing had zero spin, so points riding on it did not turn with it. It now spins as its fitted frame turns, and every free point tied to a wing, rigid or fitted, is damped against that wing's frame, as on the monolith.
- A spinning rigid wing (`AeroDirect`, `AeroLinearized`) lost a few percent of its aerodynamic force. Its coefficients were scaled by the panel inflow, which includes the spin, and turned back into force with the wing's apparent wind alone. At a 0.1 rad/s pitch rate the force was 2.6% low.
