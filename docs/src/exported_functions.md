```
@meta
CurrentModule = SymbolicAWEModels
```

# API Reference

This page provides a detailed reference for all public functions exported by the 
`SymbolicAWEModels.jl` package.

## High-Level Simulation Functions

These functions provide convenient wrappers for running common simulation scenarios.

```@docs
sim!
sim_oscillate!
sim_turn!
```

## Low-Level Simulation and Analysis

These functions provide direct control over the simulation and tools for model analysis.

```@docs
init!
next_step!
find_steady_state!
linearize!
copy_to_simple!
simple_linearize!
```

## System Structure Constructors

These functions are used to procedurally generate predefined [`SystemStructure`](@ref) topologies.

```@docs
create_ram_sys_struct
create_tether_sys_struct
create_simple_ram_sys_struct
```

## State Accessor Functions (Getters)

Use these functions to retrieve state information and calculated values from a model instance.

```@docs
winch_force
unstretched_length
tether_length
```

## Utility and Helper Functions

General helper functions for package management and setup.

```@docs
init_module
```
