```@meta
CurrentModule = SymbolicAWEModels
```
## Introduction
Most of the functions work on a [`SymbolicAWEModel`](@ref) object.
For this, the variable `s` is used.
Such a variable can be created with the lines:
```julia
set = Settings("system_ram.yaml")
s = SymbolicAWEModel(set)
```
Functions with an "!" as last character of the function name modify one of more of their
parameters, in this context mostly the variable s.

## Input functions
```@docs
set_depower_steering!
set_v_wind_ground!
```

## Output functions
```@docs
unstretched_length
tether_length
calc_height
winch_force
spring_forces
calc_aoa
pos
```

## Simulation interface
```@docs
init_sim!
next_step!
```

## Helper functions
```@docs
copy_examples
copy_bin
```
