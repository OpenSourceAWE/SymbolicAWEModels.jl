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
check_live_polar
```

## YAML loading

```@docs
load_sys_struct_from_yaml
```

## System configuration

```@docs
set_world_frame_damping
set_body_frame_damping
set_angular_damping
calc_steady_torque
```

## Winch components

The winch motor model is a pluggable [`AbstractWinchModel`](@ref) struct.
Pass one to `Winch(...; model=...)` to override the default
torque-driven motor ([`TorqueWinch`](@ref)). The builder is
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

SymbolicAWEModels provides plotting functionality through a package extension
that loads once a Makie backend and `MakieControlPlots` are both available.

### 3D system visualization

Plot the 3D structure of the system with interactive features:
```julia
import GLMakie
using MakieControlPlots
plot(sys::SystemStructure; kwargs...)
```

**Keyword arguments:**
- `vector_scale::Float64=1.0`: Length scale of the force/orientation arrows
- `force_color::Bool=false`: Colour segments by tension instead of `segment_color`
- `segment_color=RGBf(0.25, 0.25, 0.25)`: Default colour for segments
- `relmargin::Float64=0.2`: Margin around the system, as a fraction of its extent
- `body_frame`: Track a wing body frame with the camera (defaults on with `aero_mapping`)
- `zoom`, `pan_horizontal`, `pan_vertical`, `tilt_horizontal`, `tilt_vertical`: Camera placement
- `plot_aero::Bool=true`: Draw the aerodynamic geometry of each wing
- `extra_points`, `extra_groups`, `mesh`: Extra geometry to overlay

Remaining keywords are forwarded to `plot!`, including:
- `show_points`, `show_segments`, `show_orient`, `show_beams`: Layer visibility
- `transparency::Bool=true`: Order-independent transparency; `false` is much faster
- `aero_mapping::Bool=false`: Overlay the [`AeroPressure`](@ref) station→point map
- `linewidth`, `point_size`, `beam_color`, `airfoil_color`, …: Styling

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

**Keyword arguments** — `plot_default::Bool=true` switches the default panel
set on or off as a group; the panels below default to `plot_default`:

- `plot_reelout`: Reel-out velocities of the steering winches
- `plot_aero_force`: z-component of the aerodynamic force
- `plot_aoa`: Angle of attack
- `plot_heading`: Heading and course angles
- `plot_winch_force`: Winch forces

Opt-in panels (all default `false`): `plot_twist`, `plot_turn_rates`,
`plot_turn_radius`, `plot_aero_moment`, `plot_tether_moment`, `plot_tether`,
`plot_tether_actual`, `plot_v_app`, `plot_elevation`, `plot_azimuth`,
`plot_distance`, `plot_yaw_rate`, `plot_cone_angle`, `plot_old_heading`,
`plot_kiteutils_course`, `plot_set_values`.

Appearance: `suffix::String=" - " * sys.name`, `size::Tuple=(1200, 800)`,
`label_fontsize::Int=16`, `ticklabelsize::Int=12`, `legendsize::Int=10`, and
the per-panel limits `aoa_ylims`, `gk_ylims`, `turn_radius_ylims`.

Passing a `Vector{SysLog}` instead of a single log overlays several runs on the
same panels for comparison.

!!! note "Extension loading"
    The `plot` functions become available once both a Makie backend and
    `MakieControlPlots` are loaded — `using GLMakie` on its own is not enough.
    `plot` extends `MakieControlPlots.plot`, the generic of the figure-returning
    plot commands; the scene-mutating `plot!` extends `Makie.plot!`. A Makie
    backend exports a `plot` of its own, so load the backend with `import
    GLMakie` rather than `using GLMakie` to keep the name unambiguous.

## Inflated-tube rigidity laws

Rigidities for the [`TimoshenkoJoint`](@ref)s of a beam wing whose leading edge
and struts are pressurised fabric tubes. [`tube_bending_law`](@ref) and
[`tube_torsion_law`](@ref) come from the empirical Breukels correlations;
[`comer_levy_bending_law`](@ref) is the analytical alternative that stays valid
past collapse but needs the fabric membrane stiffness `E·t`, which
[`breukels_membrane_stiffness`](@ref) can supply.

```@docs
tube_linear_rigidities
tube_bending_law
tube_torsion_law
tube_mass
breukels_tip_force
breukels_collapse_deflection
breukels_membrane_stiffness
comer_levy_bending_stiffness
comer_levy_wrinkling_moment
comer_levy_collapse_moment
comer_levy_bending_law
membrane_linear_rigidities
frame_quaternion
frame_quaternion_xy
```

## Unsteady aerodynamics

```@docs
unsteady_aero
apply_apparent_mass!
```

## Utility and helper functions

```@docs
init_module
```
