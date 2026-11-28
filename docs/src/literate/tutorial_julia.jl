# Copyright (c) 2025 Bart van de Lint, Jelle Poland
# SPDX-License-Identifier: LGPL-3.0-only

#md # ```@meta
#md # CurrentModule = SymbolicAWEModels
#md # ```

# # Building a system using Julia
#
# This tutorial walks through building mechanical systems from scratch
# in Julia. Each component is a KiteGeometry row built with keywords,
# with a SymbolicAWEModels model in its `model` column; the rows are
# assembled into a `SystemDefinition`, which is then compiled into a
# simulation model.
#
# ## Overview
#
# The construction workflow has four steps:
#
# 1. **Define components** — create KiteGeometry `Point` and `Segment`
#    rows, the `Material` table the segments name, a KiteGeometry
#    `Transform`, and optionally `Tether`, `Winch` and `Pulley` rows
# 2. **Assemble** — pass the rows to [`system_definition`](@ref), which
#    builds and places the `SystemDefinition`
# 3. **Run** — wrap it in a [`SymbolicAWEModel`](@ref), the run of the
#    structure; [`init!`](@ref) generates the symbolic equations via
#    ModelingToolkit, compiles an `ODEProblem` and validates the run
# 4. **Simulate** — call [`next_step!`](@ref)
#
# All components use symbolic names (`:anchor`, `:spring`, etc.) and
# reference each other by name. [`system_definition`](@ref) resolves
# these references to numeric indices automatically. A reference can
# also be an integer index.

# ## Step 1: a simple tether
#
# We start with the simplest possible system: a chain of point masses
# connected by spring-damper segments, hanging under gravity.

using SymbolicAWEModels, VortexStepMethod
using SymbolicAWEModels: KiteGeometry
import GLMakie
using MakieControlPlots
GLMakie.activate!(; visible=false)                                #hide
ASSETS = joinpath(@__DIR__, "..", "assets")                        #hide
pkg_root = joinpath(@__DIR__, "..", "..", "..")                    #hide
set_data_path(joinpath(pkg_root, "data", "base"))                 #hide

set = Settings("system.yaml")
set.solver = "FBDF"
set.v_wind = 1.0

n_segments = 20
l_tether = 50.0

# First, define the points. The anchor is fixed in space
# (`fixed=true`). The intermediate points are free (governed by
# Newton's second law). The tip point has some extra mass.

points = Point[]

push!(points, Point(; name=:anchor, pos_ENU=zeros(3), fixed=true))

for i in 1:n_segments
    pos = [0.0, 0.0, i * l_tether / n_segments]
    if i < n_segments
        push!(points, Point(; name=Symbol("p_$i"), pos_ENU=pos))
    else
        push!(points, Point(; name=:tip, pos_ENU=pos,
            extra_mass=1.0))
    end
end

# A point's row gives one of three of the four
# [`DynamicsType`](@ref)s:
# - `STATIC` — `fixed=true`: the point is welded in the world frame
# - `DYNAMIC` — no `fixed`, `body` or `tube`: the point moves
#   according to ``\ddot{\mathbf{r}} = \mathbf{F}/m``
# - `BODY_STATIC` — a `body` (or a `tube`): the point rides that body
#   (or beam element) rigidly and feeds its net force into it
#
# The fourth, `KINEMATIC`, is a `type` of a `SpanwiseStation`'s
# [`StationModel`](@ref) whose deflection is prescribed by geometry,
# not of a point. Further point columns of the run, such as damping
# or drag, go in a [`PointModel`](@ref) in the point's `model` column.
#
# Next, connect the points with segments. Each segment is a
# spring-damper element with explicit stiffness and damping per unit
# length. `KiteGeometry.line_segment` builds a `Segment` row and adds
# the material of its stiffness and density to the `materials` table;
# the damping goes in its [`SegmentSpring`](@ref) model.

unit_stiffness = 614600.0  # [N]
unit_damping = 473.0       # [Ns]
diameter = 0.002           # [m]

materials = KiteGeometry.MaterialTable(KiteGeometry.Material[], Dict{String, Any}())
segments = Segment[]

for i in 1:n_segments
    p_i = points[i].name
    p_j = points[i+1].name
    push!(segments,
        KiteGeometry.line_segment(Symbol("seg_$i"), (p_i, p_j), set, materials;
            diameter, unit_stiffness, model=SegmentSpring(; unit_damping)))
end

# Now assemble everything into a `SystemDefinition`.
#
# A KiteGeometry `Transform` describes the initial orientation of the
# system using spherical coordinates (elevation, azimuth, heading)
# relative to a base position. It names the base point, the point it
# turns and the points it moves by their index; [`place!`](@ref) puts
# the structure where it says. Then visualize it, as a model built on it
# without settings (building one compiles nothing):

sys = system_definition("tether", set; points, segments, materials)
transforms = [Transform(; name=:tf, elevation=deg2rad(-80), azimuth=0.0,
    heading=0.0, base_pos=[0.0, 0.0, 50.0], base_point=1,
    rot_point=[length(points) => 1.0], points=collect(eachindex(points)))]
place!(sys, transforms)

scene = plot(SymbolicAWEModel(nothing, sys))                      #hide
plot(SymbolicAWEModel(nothing, sys))

GLMakie.save(joinpath(ASSETS, "tether_structure.png"), scene)     #hide

#md # ![SystemDefinition visualization](assets/tether_structure.png)
#
# If the system looks correct, compile and simulate:

using KiteUtils: init!, next_step!, update_sys_state!

sam = SymbolicAWEModel(set, sys; transforms)
init!(sam; prn=false, remake=true)                                             #hide
init!(sam)

logger = Logger(sam, 201)
sys_state = SysState(sam)
log!(logger, sys_state)

for i in 1:200
    next_step!(sam)
    update_sys_state!(sys_state, sam)
    log!(logger, sys_state)
end

save_log(logger, "tether_sim")
lg = load_log("tether_sim")
SymbolicAWEModels.record(lg, sam, "tether_sim.gif")
SymbolicAWEModels.record(lg, sam,                                 #hide
    joinpath(ASSETS, "tether_sim.gif"); framerate=20)             #hide

#md # ![Tether simulation](assets/tether_sim.gif)

# ## Step 2: adding a winch
#
# To control tether length, we group the segments into a `Tether`
# row and connect it to a `Winch` row. The winch applies a torque to
# reel the tether in or out; its friction and inertia go in its
# [`WinchModel`](@ref).

## Group all segments into a tether. The points already sit at the
## right spacing, so the tether takes its length from that geometry;
## pass a stretched standoff only to reposition a root tether.
seg_names = [s.name for s in segments]
tethers = [Tether(; name=:main, segments=seg_names)]

## Create a winch connected to the tether
gear_ratio = 1.0
drum_radius = 0.11          # [m]
coulomb_friction = 122.0    # [N]
viscous_coefficient = 30.6  # [Ns/m]
inertia_total = 0.024       # [kgm²]

winches = [Winch(; name=:winch, tethers=[:main], winch_point=:anchor,
    gear_ratio, drum_radius,
    model=WinchModel(; coulomb_friction, viscous_coefficient, inertia_total))]

# Build and simulate with a constant torque of -20 Nm:

sys = system_definition("winch", set; points, segments, materials, tethers, winches)
place!(sys, transforms)

scene = plot(SymbolicAWEModel(nothing, sys))                      #hide
GLMakie.save(joinpath(ASSETS, "winch_structure.png"), scene)      #hide

#md # ![Winch system structure](assets/winch_structure.png)

using KiteUtils: init!, next_step!, update_sys_state!

sam = SymbolicAWEModel(set, sys; transforms)
init!(sam; prn=false, remake=true)                                             #hide
init!(sam)

logger = Logger(sam, 201)
sys_state = SysState(sam)
log!(logger, sys_state)

for i in 1:200
    next_step!(sam; set_values=[-20.0])
    update_sys_state!(sys_state, sam)
    log!(logger, sys_state)
end

save_log(logger, "winch_sim")
lg = load_log("winch_sim")
SymbolicAWEModels.record(lg, sam, "winch_sim.gif")
SymbolicAWEModels.record(lg, sam,                                 #hide
    joinpath(ASSETS, "winch_sim.gif"); framerate=20)              #hide

#md # ![Winch simulation](assets/winch_sim.gif)

# ## Step 3: adding a pulley
#
# A `Pulley` row enforces length redistribution between two
# segments that share a common point. This models a physical pulley
# where the total rope length through the pulley is conserved.

set = Settings("system.yaml")
set.solver = "FBDF"
set.v_wind = 1.0
set.abs_tol = 1e-4
set.rel_tol = 1e-4

# Create points — two static anchor points, and two dynamic points
# forming the pulley system:

points = [
    Point(; name=:left,   pos_ENU=[0.0, 0.0, 2.0], fixed=true),
    Point(; name=:right,  pos_ENU=[2.0, 0.0, 2.0], fixed=true),
    Point(; name=:pulley, pos_ENU=[0.1, 0.0, 1.0]),
    Point(; name=:weight, pos_ENU=[0.1, 0.0, 0.0], extra_mass=0.1),
]

materials = KiteGeometry.MaterialTable(KiteGeometry.Material[], Dict{String, Any}())
spring = SegmentSpring(; unit_damping)
segments = [
    KiteGeometry.line_segment(:to_left,  (:pulley, :left), set, materials;
        diameter, unit_stiffness, model=spring),
    KiteGeometry.line_segment(:to_right, (:pulley, :right), set, materials;
        diameter, unit_stiffness, model=spring),
    KiteGeometry.line_segment(:hanging,  (:pulley, :weight), set, materials;
        diameter, unit_stiffness, model=spring),
]

# The pulley connects the two upper segments. When one gets shorter,
# the other gets longer by the same amount. Its [`PulleyModel`](@ref)
# gives its dynamics type:

pulleys = [Pulley(; name=:pulley, segments=(:to_left, :to_right),
    model=PulleyModel(; type="DYNAMIC"))]

# Build the system:

sys = system_definition("pulley", set; points, segments, materials, pulleys)

scene = plot(SymbolicAWEModel(nothing, sys))                      #hide
GLMakie.save(joinpath(ASSETS, "pulley_structure.png"), scene)     #hide

#md # ![Pulley system structure](assets/pulley_structure.png)

using KiteUtils: init!, next_step!, update_sys_state!

sam = SymbolicAWEModel(set, sys)
init!(sam; prn=false, remake=true)                                             #hide
init!(sam)

logger = Logger(sam, 201)
sys_state = SysState(sam)
log!(logger, sys_state)

for i in 1:200
    next_step!(sam)
    update_sys_state!(sys_state, sam)
    log!(logger, sys_state)
end

save_log(logger, "pulley_sim")
lg = load_log("pulley_sim")
SymbolicAWEModels.record(lg, sam, "pulley_sim.gif")
SymbolicAWEModels.record(lg, sam,                                 #hide
    joinpath(ASSETS, "pulley_sim.gif"); framerate=20)             #hide

#md # ![Pulley simulation](assets/pulley_sim.gif)

# ## Step 4: wings
#
# Wings can be added as KiteGeometry `Wing`, `SpanwiseStation` and
# `Body` rows with a [`WingModel`](@ref), [`StationModel`](@ref) and
# [`BodyModel`](@ref), or via YAML — both paths produce the same
# `SystemDefinition`. See the
# [2-Plate Kite example](examples.md#plate-kite-2) for a complete
# wing system, and [VSM Coupling](vsm_coupling.md) for how
# aerodynamic forces are computed.

# ## Named references
#
# A reference in a row is a symbolic name (`:anchor`, `:spring`) or
# an integer index (`1`, `2`). Components reference each other by
# name:
#
# ```julia
# ## Segment references points by name
# KiteGeometry.line_segment(:spring, (:anchor, :mass), set, materials;
#     diameter=0.002, unit_stiffness=1000.0,
#     model=SegmentSpring(; unit_damping=50.0))
#
# ## Tether references segments by name
# Tether(; name=:main, segments=[:seg1, :seg2])
#
# ## Winch references tethers and winch point by name
# Winch(; name=:winch, tethers=[:main], winch_point=:anchor,
#     gear_ratio=1.0, drum_radius=0.11,
#     model=WinchModel(; coulomb_friction=122.0, viscous_coefficient=30.6,
#         inertia_total=0.024))
#
# ## Pulley references segments by name
# Pulley(; name=:pulley, segments=(:seg_a, :seg_b),
#     model=PulleyModel(; type="DYNAMIC"))
#
# ```
#
# [`system_definition`](@ref) resolves all names to indices. After construction, you can access components by name:
#
# ```julia
# sys.points[:anchor]       # Access point by name
# sys.segments[:spring]     # Access segment by name
# sys.winches[:winch]       # Access winch by name
# ```

# ## Component summary
#
# | KiteGeometry row | Constructor | Model | Purpose |
# |------|------------|-------|---------|
# | `Point` | `Point(; name, pos_ENU, fixed, body, tube, extra_mass, model)` | [`PointModel`](@ref) | Mass node |
# | `Segment` | `KiteGeometry.line_segment(name, (p_i, p_j), set, materials; ...)` | [`SegmentSpring`](@ref) | Spring-damper |
# | `Tether` | `Tether(; name, segments, start_point, end_point, ...)` | — | Winch-controlled segments |
# | `Winch` | `Winch(; name, tethers, winch_point, gear_ratio, drum_radius, model)` | [`WinchModel`](@ref) | Torque-controlled motor |
# | `Pulley` | `Pulley(; name, segments, model)` | [`PulleyModel`](@ref) | Equal-tension constraint |
# | `SpanwiseStation` | `SpanwiseStation(; name, wing, points, model)` | [`StationModel`](@ref) | Wing twist section |
# | `Wing` | `Wing(; name, pos_ENU, Q_KA2ENU, body, model)` | [`WingModel`](@ref) | Lifting surface |
# | `Body` | `Body(; name, pos_ENU, Q_KA2ENU, extra_mass, extra_inertia_KA, fixed, model)` | [`BodyModel`](@ref) | Rigid body |
# | `Tube` | `Tube(; name, bodies, diameter, pressure, material, model)` | [`TimoshenkoTubeModel`](@ref), [`ElasticTubeModel`](@ref) | Inflated tube between bodies |
# | `Material` | `Material(; name, law, parameters)` | — | Stiffness law a segment or tube names |
# | `Transform` | `Transform(; name, elevation, azimuth, heading, ...)` | — | Spherical positioning |
#
# A model's fields all default to `nothing`, the run's default. See
# the [Types](exported_types.md) page for the model types.
