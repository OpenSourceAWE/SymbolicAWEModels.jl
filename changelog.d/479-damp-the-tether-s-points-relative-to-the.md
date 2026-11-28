### Added

- A tether's inner points can damp their velocity relative to the frame that turns and stretches with the tether's straight line from its anchor to its free end: their `body_frame_damping`, set per point or for all of a tether's inner points by the tether's `body_frame_damping` column or keyword. This damps transverse waves and bouncing in the line and leaves its quasi-steady motion alone. It is off by default, and the model cache key counts the damped tether points.
