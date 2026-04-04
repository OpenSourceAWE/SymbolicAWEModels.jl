# SPDX-FileCopyrightText: 2025 Bart van de Lint
#
# SPDX-License-Identifier: MPL-2.0

"""
Progressive tutorial: builds up from a simple tether to a full kite
model, adding a winch, pulley, and wing step by step.
"""

using SymbolicAWEModels, VortexStepMethod, LinearAlgebra
using SymbolicAWEModels: Point

set_data_path("data/2plate_kite")
set = Settings("system.yaml")
set.segments = 20
set.l_tether = 50.0

# Segment properties (Dyneema, 4mm diameter)
seg_stiffness = 614_600.0  # EA [N]
seg_damping   = 473.0      # [N·s]
seg_diameter  = 0.004      # [m]

# --- STEP 1: Tether ---

points = [Point(1, zeros(3), STATIC)]
segments = Segment[]
l_seg = set.l_tether / set.segments
for i in 1:set.segments
    pos = [0.0, 0.0, i * l_seg]
    kw = i == set.segments ? (; extra_mass=1.0) : (;)
    push!(points, Point(i + 1, pos, DYNAMIC; kw...))
    push!(segments, Segment(i, i, i + 1,
        seg_stiffness, seg_damping, seg_diameter;
        l0=l_seg))
end

transforms = [
    Transform(1, deg2rad(-80), 0, 0;
              base_pos=[0, 0, 0],
              base_point=1, rot_point=length(points)),
]

sys = SystemStructure("tether", set;
                      points, segments, transforms)
sam = SymbolicAWEModel(set, sys)
init!(sam; remake=false)

for i in 1:80
    next_step!(sam)
end
@info "Tether simulation completed" steps=80

# --- STEP 2: Add a winch ---

set.v_wind = 0.0
n_seg = length(segments)
tethers = [Tether(:main, collect(1:n_seg))]
winches = [Winch(:winch, set, [:main]; winch_point=1)]

sys = SystemStructure("winch", set;
    points, segments, tethers, winches, transforms)
sam = SymbolicAWEModel(set, sys)
init!(sam; remake=false)

for i in 1:80
    next_step!(sam; set_values=[-20.0])
end
@info "Winch simulation completed" steps=80

# --- STEP 3: Add a pulley ---

push!(points, Point(22, [0, 0, set.l_tether + 5], DYNAMIC; extra_mass=0.1))
push!(points, Point(23, [1, 0, set.l_tether + 5], STATIC))
push!(segments, Segment(21, 21, 22,
    seg_stiffness, seg_damping, seg_diameter))
push!(segments, Segment(22, 21, 23,
    seg_stiffness, seg_damping, seg_diameter))
pulleys = [Pulley(1, 21, 22, DYNAMIC)]
transforms = [
    Transform(1, deg2rad(-85), 0, 0;
              base_pos=[0, 0, 0], base_point=1,
              rot_point=21),
]

sys = SystemStructure("pulley", set;
    points, segments, tethers, winches, pulleys, transforms)
sam = SymbolicAWEModel(set, sys)
init!(sam; remake=false)

for i in 1:80
    next_step!(sam; set_values=[-10.0])
end
@info "Pulley simulation completed" steps=80

# --- STEP 4: Add a kite ---

vsm_set = VortexStepMethod.VSMSettings(
    joinpath(get_data_path(), "vsm_settings.yaml");
    data_prefix=false)
vsm_wing = VortexStepMethod.Wing(vsm_set)
vsm_aero = BodyAerodynamics([vsm_wing])
vsm_solver = Solver(vsm_aero;
    solver_type=NONLIN, atol=2e-8, rtol=2e-8)
wings = [SymbolicAWEModels.Wing(1, vsm_aero, vsm_wing,
    vsm_solver, [], I(3), [0.5, 0, set.l_tether + 6])]

# WING-type points: 3 LE/TE pairs matching the 3 aero sections
wing_z = set.l_tether + 6
for (i, y) in enumerate([-1.0, 0.0, 1.0])
    n = length(points)
    push!(points, Point(n + 1, [-0.5, y, wing_z],
        WING; wing=1, transform=1, extra_mass=0.1))
    push!(points, Point(n + 2, [0.5, y, wing_z],
        WING; wing=1, transform=1, extra_mass=0.1))
end

sys = SystemStructure("wing", set;
    points, segments, tethers, winches, pulleys,
    wings, transforms)
sam = SymbolicAWEModel(set, sys)
@info "Wing model created successfully"
