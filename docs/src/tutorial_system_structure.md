```@meta
CurrentModule = SymbolicAWEModels
```
# Custom Model Tutorial

This tutorial demonstrates how to create custom kite power system models using `SystemStructure`.

!!! info "API Documentation"
    For detailed API documentation, see:
    - [System Structure](@ref) - Component types and reference to KiteUtils documentation
    - [SymbolicAWEModel](@ref) - Main simulation type
    - [State Management](@ref) - Working with system state

A custom `SystemStructure` can be used to create models of kite power systems of almost any configuration:
- Custom amount of tethers
- Custom bridle configurations
- Quasi-static or dynamic point masses
- Different amounts of stiffness, damping and diameter on different tether segments

## Precondition
First, following the [Home](index.md) section up to the installation of the examples. To start Julia,
either use `julia --project=.`, or `./bin/run_julia`. Make sure that
at least `SymbolicAWEModels` version 0.2 is installed by typing `]status` in the Julia REPL. 
If necessary, update `SymbolicAWEModels` by typing `]up`.

## Creating a simple tether

We start by loading the necessary packages and defining settings and parameters.

```julia
using SymbolicAWEModels, VortexStepMethod, Makie

set = Settings("system.yaml")
set.segments = 20
set.l_tether = 50.0
dynamics_type = DYNAMIC
```

Then, we define vectors of the system structure types we are going to use. For this simple example we only need points and segments.

```julia
points = Point[]
segments = Segment[]

points = push!(points, Point(1, zeros(3), STATIC; wing_idx=0))
```

The first point we add is a static point. There are four different `DynamicsType`s to choose from: `STATIC`, `QUASI_STATIC`, `DYNAMIC` and `WING`. `STATIC` just means that the point doesn't move. `DYNAMIC` is a point modeled with acceleration, while `QUASI_STATIC` constrains this acceleration to be zero at all times. A `WING` point is connected to a wing body. (See [KiteUtils DynamicsType docs](https://opensourceawe.github.io/KiteUtils.jl/dev/system_structure/#KiteUtils.DynamicsType) for details.)

Now we can add `DYNAMIC` points and connect them to each other with segments. `BRIDLE` segments don't need to have a tether, because they have a constant unstretched length.
```julia
segment_idxs = Int[]
for i in 1:set.segments
    global points, segments
    point_idx = i+1
    pos = [0.0, 0.0, i * set.l_tether / set.segments]
    if i < set.segments
        push!(points, Point(point_idx, pos, dynamics_type; wing_idx=0))
    else
        push!(points, Point(point_idx, pos, dynamics_type; mass=1.0, wing_idx=0))
    end
    segment_idx = i
    push!(segments, Segment(segment_idx, set, (point_idx-1, point_idx), BRIDLE))
    push!(segment_idxs, segment_idx)
end
```

In order to describe the initial orientation of the structure, we define a `Transform(idx, elevation, azimuth, heading)` with an elevation (-80 degrees), azimuth and heading, and a base position `[0.0, 0.0, 50.0]`.
```julia
transforms = [Transform(1, deg2rad(-80), 0.0, 0.0; 
              base_pos = [0.0, 0.0, 50.0], base_point_idx=points[1].idx,
              rot_point_idx=points[end].idx)]
```

From the points, segments and transform we create a [`SystemStructure(name, set)`](@ref), which can be plotted in 3D to quickly investigate if the model is correct.
```julia
sys_struct = SystemStructure("tether", set; points, segments, transforms)
plot(sys_struct, 0.0)
```
![SystemStructure visualization](assets/tether_sys_struct.png)

If the system looks good, we can easily model it, by first creating a [`SymbolicAWEModel`](@ref), initializing it and stepping through time.
```julia
sam = SymbolicAWEModel(set, sys_struct)

init!(sam; remake=false)
plot(sys_struct, 0.0)  # Create initial plot
for i in 1:80
    next_step!(sam)
    plot(sys_struct, i/set.sample_freq)  # Update plot
end
```
![Tether during simulation](assets/tether_during_sim.png)

# Using a winch and a tether

Let's try to adjust the length of the tether in the last example. To do this we first need to create a set of segments with a common changing `l0`, called a `Tether`.
```julia
set.v_wind = 0.0
tethers = [Tether(1,[segment.idx for segment in segments])]
```
As you can see, we just add all of the segments from the simple tether to our `Tether` struct.
The next step is to create a winch. Each winch can be connected to one or more tethers, so it is possible to connect multiple tethers to the same winch.
The winch is a torque controlled winch.

```julia
using WinchModels
winches = [Winch(1, set, [1])]
```

The plot of the `SystemStructure` should still look the same, so we don't have to plot that. We can just create the system, and simulate it. We call `plot(sys_struct, 0.0)` to create a new plot.

```julia
sys_struct = SystemStructure("winch", set; points, segments, tethers, winches, transforms)
sam = SymbolicAWEModel(set, sys_struct)
init!(sam; remake=false)
ss = SysState(sam)

plot(sys_struct, 0.0)  # Create initial plot
for i in 1:80
    next_step!(sam; set_values=[-20.0])
    update_sys_state!(ss, sam)
    plot(sys_struct, i/set.sample_freq)  # Update plot
end
```

# Using a pulley

First, we need to update some settings. `l_tether` is specified such that the plot window is zoomed in correctly.

```julia
using SymbolicAWEModels, VortexStepMethod, Makie

set = se("system.yaml")
set.v_wind = 10.0
set.l_tether = 5.0
set.abs_tol = 1e-4
set.rel_tol = 1e-4
dynamics_type = DYNAMIC
```

Now we create points and segments in a similar manner as in the last example. The `mass` keyword can be used to specify the mass of the point itself. When `mass=0.0` the mass of the point just consists of the tether/segment mass.

```julia
points = Point[]
segments = Segment[]
pulleys = Pulley[]

push!(points, Point(1, [0.0, 0.0, 2.0], STATIC))
push!(points, Point(2, [2.0, 0.0, 2.0], STATIC))
push!(points, Point(3, [0.1, 0.0, 1.0], DYNAMIC))
push!(points, Point(4, [0.1, 0.0, 0.0], DYNAMIC; mass=0.1))

push!(segments, Segment(1, set, (3,1), BRIDLE))
push!(segments, Segment(2, set, (3,2), BRIDLE))
push!(segments, Segment(3, set, (3,4), BRIDLE))
```

Pulleys can be modeled when three or more `Segment`s are connected to a common `Point`. When creating a pulley, only two segments are specified: these are the segments of the tether moving through the pulley.

```julia
push!(pulleys, Pulley(1, (1,2), DYNAMIC))
```

We can then use a `Transform` to describe the orientation of the initial system. 

```julia
transforms = [Transform(1, -deg2rad(0.0), 0.0, 0.0; base_pos=[1.0, 0.0, 4.0], base_point_idx=1, rot_point_idx=2)]
sys_struct = SymbolicAWEModels.SystemStructure("pulley", set; points, segments, pulleys, transforms)
plot(sys_struct, 0.0)
```

If the plot of the `SystemStructure` looks good, we can continue by creating a [`SymbolicAWEModel`](@ref) and simulating through time.

```julia
sam = SymbolicAWEModel(set, sys_struct)
init!(sam; remake=false)
plot(sys_struct, 0.0)  # Create initial plot
for i in 1:100
    next_step!(sam)
    plot(sys_struct, i/set.sample_freq)  # Update plot
end
```


