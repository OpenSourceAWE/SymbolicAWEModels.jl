### Changed

- BREAKING: Angular damping and body-frame damping no longer brake a kite's flight over the sphere about its tether's base point. `angular_damping` damps a body's spin relative to the rotation `r × v / |r|²` of the wing it flies with over that sphere, and the `body_frame_damping` of a wing's nodes and of a body with a parent wing damps their velocity relative to the wing's velocity turned by that rotation. A kite flying across the sphere therefore trims differently than before. A body that no transform moves damps its absolute spin as before.
- The model cache key carries the tether base point each body damps against, so every cached model with a body recompiles once.

### Fixed

- A spinning rigid wing (`AeroDirect`, `AeroLinearized`) lost a few percent of its aerodynamic force. Its coefficients were scaled by the panel inflow, which includes the spin, and turned back into force with the wing's apparent wind alone. At a 0.1 rad/s pitch rate the force was 2.6% low.
