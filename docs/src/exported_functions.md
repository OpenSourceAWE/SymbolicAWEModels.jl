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

KiteGeometry's `load_authoring` reads and places an authoring YAML, and the structure
is built from the definition it returns, with the columns SymbolicAWEModels reads held
by the models it registers with KiteGeometry.

## Structure documents

A `SystemDefinition` reads and writes as its awesIO *structure document* through KiteGeometry's `to_yaml`, `from_yaml`, `to_json` and `from_json`: one `headers`/`units`/`data` table per block, conforming to the `structure_schema.yml` KiteGeometry ships. Every position is in the ENU world frame, each body at
the centre of its own mass; the authoring geometry, the transforms and the live state
do not survive a round trip. `save_log(logger, sys, name)` writes a log that carries the document as JSON under the metadata key `topology`.

## System configuration

```@docs
set_world_frame_damping
set_body_frame_damping
set_angular_damping
calc_steady_torque
```

## Placing a SystemDefinition

[`place!`](@ref) puts a structure in the world by KiteGeometry transforms and makes where it lands its initial pose (`pos_ENU`, `Q_KA2ENU`); `system_definition` and `load_authoring` run it. A model, the run of a structure, starts there ([`start_run!`](@ref)), and `init!(sam)` takes the run as it stands, so a run restarted from a logged state resumes. Mass properties and rest geometry are not stored: the model works them out from its definition's initial pose whenever they are read ([`run_value`](@ref)). To move the system for a new run, change its transforms and `place!(sam)` it again. To change only part of it — a tether's unstretched length, say, without moving the kite — change the run before `init!`:

```julia
set_unstretched_length!(sam, 1, 240.0)
init!(sam)
```

```@docs
place!
Pose(::SymbolicAWEModel)
start_run!
reset_to_initial_pose!
remake_wing_aero!
init_wind!
set_unstretched_length!
```

## Wind

Where a model takes its wind from is a [`WindMode`](@ref) it is built with. The default [`ProfileWind`](@ref) scales the ground wind `set.wind_vec` by the height profile `set.profile_law`. [`PerPointWind`](@ref) instead gives every point its own wind vector, settable between steps.

```julia
sys, transforms = load_authoring(path; set)
sam = SymbolicAWEModel(set, sys; transforms, wind_mode=PerPointWind())
init!(sam)
for _ in 1:steps
    for i in eachindex(sam.definition.points)    # a gust field, refreshed every step
        pos = run_value(sam, :points, i, :pos_w)
        set_run_value!(sam, :points, i, :wind_vec, my_wind_field(pos, sam.integrator.t))
    end
    set_run_value!(sam, :wings, 1, :wind_vec, my_wind_field(...))
    next_step!(sam)
end
```

Each point's own drag and apparent wind read its `wind_vec` directly, and a
segment's tether drag the mean of its two endpoints — the same averaging its
drag already applies to their velocities. A wing reads its own `wind_vec` (its body's),
which is the wind a rigid wing's aerodynamics fly in; the per-point winds reach
the VSM panels through the per-point inflow the panels already interpolate.
`init!` seeds every point and wing with `set.wind_vec`, so a model that is never
written to flies in a uniform wind.

```@docs
WindMode
ProfileWind
PerPointWind
per_point_wind
```

## Winch components

The winch motor model is a pluggable [`AbstractWinchModel`](@ref) struct.
Pass one as the `dynamics` of a winch's [`WinchModel`](@ref) to override the default
torque-driven motor ([`TorqueWinch`](@ref)). The builder is
selected by dispatch on the model type, [`winch_component`](@ref); custom
components must respect the connector contract enforced by
[`validate_winch_component`](@ref).

```@docs
winch_component
is_builtin_winch
validate_winch_component
```

## Run values

A model is the run of its structure: every value of a point, segment, station, pulley,
tether, winch, body, wing or tube of the run is read with [`run_value`](@ref) and
written with [`set_run_value!`](@ref), by collection, index or row name, and quantity.
A quantity the integrator holds as a parameter reads its definition's value before
`init!`; write the definition row (`sam.definition.segments[k].rest_length = ...`) to
change it before then.

```julia
pos = run_value(sam, :points, :kite, :pos_w)
set_run_value!(sam, :segments, :steering_left, :rest_length, 10.0)
```

```@docs
run_value
set_run_value!
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

Plot a model, the run of a system, in 3D with interactive features:
```julia
import GLMakie
using MakieControlPlots
plot(sam::SymbolicAWEModel; kwargs...)
```
To plot a structure that has no run yet, plot a model built on it without settings,
`plot(SymbolicAWEModel(nothing, sys))`; building it compiles nothing.

**Keyword arguments:**
- `vector_scale::Float64=1.0`: Length scale of the force/orientation arrows
- `force_color::Bool=false`: Colour segments by tension instead of `segment_color`
- `segment_color=RGBf(0.25, 0.25, 0.25)`: Segment colour — one colour, a vector
  with one per segment, or a function `segment -> colour`
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

To tell the wing from the tethers and the free segments, colour the segments by
[`segment_role`](@ref):
```julia
role_colors = (wing=:black, free=:steelblue, unwinched_tether=:seagreen,
               winched_tether=:darkorange)
plot(sam; segment_color=segment -> role_colors[segment_role(sam.definition, segment)])
```

**Interactive features:**
- Hover over segments to highlight them
- Click on a segment to zoom in
- Click in empty space to zoom out
- Rotate, pan, and zoom with mouse

### Static and vector figures

For a figure with labelled axes, or a vector PDF, draw the system into an axis
of your own with `plot!(ax, sam; kwargs...)`. It takes the drawing keywords of
`plot` (colours, layers, styling) but not its camera ones, draws into any Makie
axis (`Axis3`, `LScene`) and returns a `Dict` of the plot objects by layer
(`:segments`, `:points`, …). With CairoMakie as the backend, `CairoMakie.save` writes
the figure as vector graphics. Here `sam` is any model, such as the one
the [2-Plate Kite](@ref plate-kite-2) example builds:
```julia
using SymbolicAWEModels, CairoMakie, MakieControlPlots
CairoMakie.activate!()

role_colors = (wing=:black, free=:steelblue, unwinched_tether=:seagreen,
               winched_tether=:darkorange)
fig = Figure(size=(600, 500))
ax = Axis3(fig[1, 1]; aspect=:data, xlabel="x [m]", ylabel="y [m]", zlabel="z [m]")
sys = sam.definition
plot!(ax, sam; segment_color=segment -> role_colors[segment_role(sys, segment)],
      show_orient=false, show_wing_frame=false)
Legend(fig[1, 2], [LineElement(color=color) for color in values(role_colors)],
       [string(role) for role in keys(role_colors)])
CairoMakie.save("kite.pdf", fig)
```

### Time-series visualization

Plot simulation results as multi-panel time-series:
```julia
plot(sam::SymbolicAWEModel, log::SysLog; kwargs...)
```

**Keyword arguments** — `plot_default::Bool=true` switches the default panel
set on or off as a group; the panels below default to `plot_default`:

- `plot_reelout`: Reel-out velocities of the steering winches
- `plot_aero_force`: z-component of the aerodynamic force
- `plot_aoa`: Angle of attack
- `plot_heading`: Heading and course angles
- `plot_winch_force`: Winch forces

Opt-in panels (all default `false`): `plot_twist`, `plot_turn_rates`,
`plot_turn_radius`, `plot_aero_moment`, `plot_tether`, `plot_tether_actual`,
`plot_v_app`, `plot_elevation`, `plot_azimuth`, `plot_distance`,
`plot_yaw_rate`, `plot_cone_angle`, `plot_old_heading`,
`plot_kiteutils_course`, `plot_set_values`.

Appearance: `suffixes::Vector{String}` (one label suffix per log, defaulting to
the system names), `size::Tuple=(1200, 800)`, `label_fontsize::Int=16`,
`ticklabelsize::Int=12`, `legendsize::Int=10`, and `turn_radius_ylims`.

Passing a `Vector{SysLog}` instead of a single log overlays several runs on the
same panels for comparison.

!!! note "Extension loading"
    The `plot` functions become available once both a Makie backend and
    `MakieControlPlots` are loaded — `using GLMakie` on its own is not enough.
    `plot` extends `MakieControlPlots.plot`, the generic of the figure-returning
    plot commands; the scene-mutating `plot!` extends `Makie.plot!`. A Makie
    backend exports a `plot` and a `save` of its own, so with both packages
    brought in by `using` those two names are ambiguous: call them qualified,
    as `MakieControlPlots.plot(sam)` and `CairoMakie.save("kite.pdf", fig)`.

## Inflated-tube rigidity laws

Rigidities for the Timoshenko tubes ([`TimoshenkoTubeModel`](@ref)) of a beam wing whose leading edge
and struts are pressurised fabric tubes. [`tube_bending_law`](@ref) and
[`tube_torsion_law`](@ref) come from the empirical Breukels correlations;
[`comer_levy_bending_law`](@ref) is the analytical alternative that stays valid
past collapse but needs the fabric membrane stiffness `E·t`, which
[`breukels_membrane_stiffness`](@ref) can supply.

```@docs
tube_linear_rigidities
tube_law_rigidities
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
segment_role
point_type
body_type
body_dynamics_type
station_type
pulley_type
```
