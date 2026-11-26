```@meta
CurrentModule = SymbolicAWEModels
```

# Private API — equilibrium solver

The internals of [`find_equilibrium!`](@ref): the coordinates it moves, the forces and
stiffness it reads off the model, and the damped Newton search on them. Not part of
the public API.

```@docs
SymbolicAWEModels.EquilibriumLayout
SymbolicAWEModels.state_slots
SymbolicAWEModels.tether_points
SymbolicAWEModels.inertial_mass
SymbolicAWEModels.twist_inertia
SymbolicAWEModels.quaternion_product
SymbolicAWEModels.quaternion_turn
SymbolicAWEModels.skew
SymbolicAWEModels.coordinate_ranges
SymbolicAWEModels.triples
SymbolicAWEModels.largest_norm
SymbolicAWEModels.EquilibriumProblem
SymbolicAWEModels.pose_rows
SymbolicAWEModels.at_rest
SymbolicAWEModels.generalized_force
SymbolicAWEModels.stiffness
SymbolicAWEModels.move
SymbolicAWEModels.pose_moment
SymbolicAWEModels.unbalanced
SymbolicAWEModels.largest_load
SymbolicAWEModels.newton_step
SymbolicAWEModels.energy_slope
SymbolicAWEModels.line_search
SymbolicAWEModels.limit_step
SymbolicAWEModels.descent_step
SymbolicAWEModels.kite_load
SymbolicAWEModels.settle!
```
