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
sim_reposition!
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

These functions are used to procedurally generate predefined `SystemStructure` topologies.

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

## Visualization Functions

SymbolicAWEModels provides plotting functionality through a package extension that automatically loads when you import GLMakie. The extension provides the following functions:

### 3D System Visualization

Plot the 3D structure of the kite system with interactive features:
```julia
using GLMakie
plot(sys::SystemStructure; kwargs...)
```

**Keyword Arguments:**
- `size::Tuple=(1200, 800)`: Figure size in pixels
- `margin::Float64=10.0`: Margin around the system in world units
- `segment_color=:black`: Default color for segments
- `highlight_color=:red`: Color for highlighted segments
- `show_points::Bool=true`: Show point markers
- `show_segments::Bool=true`: Show tether segments
- `show_orient::Bool=true`: Show wing orientation axes

**Interactive Features:**
- Hover over segments to highlight them
- Click on a segment to zoom in
- Click in empty space to zoom out
- Rotate, pan, and zoom with mouse

### Time-Series Visualization

Plot simulation results as multi-panel time-series:
```julia
plot(sys::SystemStructure, log::SysLog; kwargs...)
```

**Keyword Arguments:**
- `plot_default::Bool=true`: Enable default plot panels
- `plot_reelout::Bool=plot_default`: Show reel-out velocities
- `plot_aero_force::Bool=plot_default`: Show aerodynamic forces
- `plot_twist::Bool=plot_default`: Show wing twist angles
- `plot_aoa::Bool=plot_default`: Show angle of attack
- `plot_heading::Bool=plot_default`: Show heading angle
- `plot_winch_force::Bool=plot_default`: Show winch forces
- `plot_aero_moment::Bool=false`: Show aerodynamic moments
- `plot_turn_rates::Bool=false`: Show angular velocities
- `plot_elevation::Bool=false`: Show elevation angle
- `plot_azimuth::Bool=false`: Show azimuth angle
- `plot_tether_moment::Bool=false`: Show tether-induced moments
- `plot_set_values::Bool=false`: Show set torque values
- `suffix::String=" - " * sys.name`: Suffix for plot labels
- `size::Tuple=(1200, 800)`: Figure size in pixels

**Example:**
```julia
using SymbolicAWEModels, GLMakie
set = Settings("system.yaml")
sam = SymbolicAWEModel(set, "ram")
init!(sam)
(log, _) = sim_oscillate!(sam)

# Plot only heading and angle of attack
plot(sam.sys_struct, log;
     plot_default=false,
     plot_heading=true,
     plot_aoa=true)
```

!!! note "Automatic Extension Loading"
    You don't need to explicitly load the plotting extension. Simply `using GLMakie` after loading SymbolicAWEModels will automatically make the `plot` functions available.

## Utility and Helper Functions

General helper functions for package management and setup.

```@docs
init_module
```
