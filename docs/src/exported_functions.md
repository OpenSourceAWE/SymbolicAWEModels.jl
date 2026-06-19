```@meta
CurrentModule = SymbolicAWEModels
```

# API reference

This page provides a detailed reference for all public functions exported by
`SymbolicAWEModels.jl`.

## High-level simulation functions

These functions provide convenient wrappers for running common simulation scenarios.

```@docs
sim!
sim_reposition!
```

## Low-level simulation and analysis

These functions provide direct control over the simulation and tools for model analysis.

```@docs
init!
next_step!
find_steady_state!
linearize!
position_slots
```

## YAML loading

```@docs
load_sys_struct_from_yaml
```

## System configuration

```@docs
set_world_frame_damping
set_body_frame_damping
calc_steady_torque
```

## Winch components

The winch motor model is a pluggable [`AbstractWinchModel`](@ref) struct.
Pass one to `Winch(...; model=...)` to override the default
torque-driven motor ([`DefaultWinchModel`](@ref)). The builder is
selected by dispatch on the model type, [`winch_component`](@ref); custom
components must respect the connector contract enforced by
[`validate_winch_component`](@ref).

```@docs
winch_component
is_builtin_winch
validate_winch_component
```

## State accessor functions

Use these functions to retrieve state information and calculated values from a model
instance.

```@docs
winch_force
unstretched_length
tether_length
segment_stretch_stats
```

## Visualization functions

SymbolicAWEModels provides plotting functionality through a package extension that
automatically loads when you import GLMakie.

### 3D system visualization

Plot the 3D structure of the system with interactive features:
```julia
using GLMakie
plot(sys::SystemStructure; kwargs...)
```

**Keyword arguments:**
- `size::Tuple=(1200, 800)`: Figure size in pixels
- `margin::Float64=10.0`: Margin around the system in world units
- `segment_color=:black`: Default color for segments
- `highlight_color=:red`: Color for highlighted segments
- `show_points::Bool=true`: Show point markers
- `show_segments::Bool=true`: Show tether segments
- `show_orient::Bool=true`: Show wing orientation axes

**Interactive features:**
- Hover over segments to highlight them
- Click on a segment to zoom in
- Click in empty space to zoom out
- Rotate, pan, and zoom with mouse

### Time-series visualization

Plot simulation results as multi-panel time-series:
```julia
plot(sys::SystemStructure, log::SysLog; kwargs...)
```

**Keyword arguments:**
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

!!! note "Automatic extension loading"
    Simply `using GLMakie` after loading SymbolicAWEModels to make
    the `plot` functions available.

## Utility and helper functions

```@docs
init_module
```
