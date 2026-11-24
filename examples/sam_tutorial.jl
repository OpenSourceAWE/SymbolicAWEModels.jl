# SPDX-FileCopyrightText: 2025 Bart van de Lint
#
# SPDX-License-Identifier: LGPL-3.0-only

"""
Progressive tutorial: builds up from a simple tether to a full kite
model, adding a winch, pulley, and wing step by step.
"""

using Pkg
if Base.active_project() != joinpath(@__DIR__, "Project.toml")
    Pkg.activate(@__DIR__)
end

using KiteUtils: init!, next_step!, update_sys_state!
using SymbolicAWEModels, VortexStepMethod
using LinearAlgebra: Diagonal

set_data_path("data/2plate_kite")
set = Settings("system.yaml")
set.segments = 20
set.l_tether = 50.0

# Segment properties (Dyneema, 4mm diameter)
seg_stiffness = 614_600.0  # EA [N]
seg_damping   = 473.0      # [N·s]
seg_diameter  = 0.004      # [m]

# Every segment names its material in this table
materials = KiteGeometry.MaterialTable(KiteGeometry.Material[], Dict{String, Any}())
dyneema(name, ends; rest_length=0.0) = KiteGeometry.line_segment(
    name, ends, set, materials; rest_length, diameter=seg_diameter,
    unit_stiffness=seg_stiffness, model=SegmentSpring(; unit_damping=seg_damping))

# --- STEP 1: Tether ---

points = [Point(; name=Symbol(1), pos_ENU=zeros(3), fixed=true)]
segments = KiteGeometry.Segment[]
l_seg = set.l_tether / set.segments
for i in 1:set.segments
    pos = [0.0, 0.0, i * l_seg]
    kw = i == set.segments ? (; extra_mass=1.0) : (;)
    push!(points, Point(; name=Symbol(i + 1), pos_ENU=pos, kw...))
    push!(segments, dyneema(Symbol(i), (i, i + 1); rest_length=l_seg))
end

transforms = [
    Transform(; name=:main, elevation=deg2rad(-80), azimuth=0, heading=0,
              base_pos=[0, 0, 0], base_point=1,
              rot_point=[length(points) => 1.0], points=collect(eachindex(points))),
]

sys = system_definition("tether", set; points, segments, materials)
place!(sys, transforms)
sam = SymbolicAWEModel(set, sys; transforms)
init!(sam; remake=false)

for _ in 1:80
    next_step!(sam)
end
@info "Tether simulation completed" steps=80

# --- STEP 2: Add a winch ---

set.v_wind = 0.0
n_seg = length(segments)
tethers = [Tether(; name=:main, segments=collect(1:n_seg))]
winches = [Winch(; name=:winch, tethers=[:main], winch_point=1)]

sys = system_definition("winch", set; points, segments, materials, tethers, winches)
place!(sys, transforms)
sam = SymbolicAWEModel(set, sys; transforms)
init!(sam; remake=false)

for _ in 1:80
    next_step!(sam; set_values=[-20.0])
end
@info "Winch simulation completed" steps=80

# --- STEP 3: Add a pulley ---

push!(points, Point(; name=Symbol(22), pos_ENU=[0, 0, set.l_tether + 5], extra_mass=0.1))
push!(points, Point(; name=Symbol(23), pos_ENU=[1, 0, set.l_tether + 5], fixed=true))
push!(segments, dyneema(Symbol(21), (21, 22)))
push!(segments, dyneema(Symbol(22), (21, 23)))
pulleys = [Pulley(; name=Symbol(1), segments=(21, 22), model=PulleyModel(; type="DYNAMIC"))]
transforms = [
    Transform(; name=:main, elevation=deg2rad(-85), azimuth=0, heading=0,
              base_pos=[0, 0, 0], base_point=1, rot_point=[21 => 1.0],
              points=collect(eachindex(points))),
]

sys = system_definition("pulley", set;
    points, segments, materials, tethers, winches, pulleys)
place!(sys, transforms)
sam = SymbolicAWEModel(set, sys; transforms)
init!(sam; remake=false)

for _ in 1:80
    next_step!(sam; set_values=[-10.0])
end
@info "Pulley simulation completed" steps=80

# --- STEP 4: Add a kite ---

vsm_set = VortexStepMethod.VSMSettings(
    project_file("vsm_settings");
    data_prefix=false)
# One station per aero section joins its LE/TE pair to the wing.
wing_z = set.l_tether + 6
# The wing's points ride its body.
stations = SpanwiseStation[]
for (i, y) in enumerate([-1.0, 0.0, 1.0])
    n = length(points)
    push!(points, Point(; name=Symbol(n + 1), pos_ENU=[-0.5, y, wing_z],
                        body=Symbol(1), extra_mass=0.1))
    push!(points, Point(; name=Symbol(n + 2), pos_ENU=[0.5, y, wing_z],
                        body=Symbol(1), extra_mass=0.1))
    push!(stations, SpanwiseStation(; name=Symbol(i), wing=Symbol(1), points=[n + 1, n + 2],
                                    model=StationModel(; type="DYNAMIC", moment_frac=0.25)))
end

wings = [KiteGeometry.Wing(; name=Symbol(1), pos_ENU=zeros(3), body=Symbol(1),
               model=WingModel(; stations=collect(eachindex(stations))))]
# No spin damping: a wing's body is otherwise damped 150 1/s about y.
bodies = [Body(; name=Symbol(1), pos_ENU=zeros(3), extra_mass=0.0,
               extra_inertia_KA=Matrix(Diagonal(ones(3))),
               model=BodyModel(; angular_damping=0.0))]
transforms = [
    Transform(; name=:main, elevation=deg2rad(-85), azimuth=0, heading=0,
              base_pos=[0, 0, 0], base_point=1, rot_point=[21 => 1.0],
              points=collect(eachindex(points)), bodies=[1]),
]

sys = system_definition("wing", set;
    points, spanwise_stations=stations, segments, materials, tethers, winches, pulleys,
    wings, bodies)
place!(sys, transforms)
sam = SymbolicAWEModel(set, sys; transforms, vsm_set)
@info "Wing model created successfully"
