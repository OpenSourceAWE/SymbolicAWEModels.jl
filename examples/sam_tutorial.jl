# SPDX-FileCopyrightText: 2025 Bart van de Lint
#
# SPDX-License-Identifier: MPL-2.0
using SymbolicAWEModels, VortexStepMethod, ControlPlots, LinearAlgebra

set = Settings("system_ram.yaml")
set.segments = 20
set.l_tether = 50.0
dynamics_type = DYNAMIC

points = Point[]
segments = Segment[]

points = push!(points, Point(1, zeros(3), STATIC; wing_idx=0))

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
    push!(segments, Segment(segment_idx, (point_idx-1, point_idx), BRIDLE))
    push!(segment_idxs, segment_idx)
end

transforms = [Transform(1, deg2rad(-80), 0.0, 0.0; 
              base_pos = [0.0, 0.0, 50.0], base_point_idx=points[1].idx,
              rot_point_idx=points[end].idx)]

sys_struct = SystemStructure("tether", set; points, segments, transforms)
plot(sys_struct, 0.0)

sam = SymbolicAWEModel(set, sys_struct)

init_sim!(sam; remake=false)
for i in 1:80
    plot(sam, i/set.sample_freq)
    next_step!(sam)
end

# ADDING A WINCH
set.v_wind = 0.0
tethers = [Tether(1,[segment.idx for segment in segments])]

using WinchModels
wm = TorqueControlledMachine(set)
winches = [Winch(1, wm, [1])]

sys_struct = SystemStructure("winch", set; points, segments, tethers, winches, transforms)
@show set.v_wind
sam = SymbolicAWEModel(set, sys_struct)
init_sim!(sam; remake=false)
ss = SysState(sam)

for i in 1:80
    plot(sam, (i-1)/set.sample_freq)
    next_step!(sam; set_values=[-20.0])
    update_sys_state!(ss, sam)
end
@show ss.l_tether[1]

# ADDING A PULLEY
push!(points, Point(22, [0, 0, set.l_tether+5], DYNAMIC))
push!(points, Point(23, [1, 0, set.l_tether+5], STATIC))
push!(segments, Segment(21, (21,22), BRIDLE))
push!(segments, Segment(22, (21,23), BRIDLE))
pulleys = [Pulley(1, (21,22), DYNAMIC)]
transforms[1].elevation = deg2rad(-85.0)
sys_struct = SystemStructure("pulley", set; points, segments, tethers, winches, pulleys, transforms)
plot(sys_struct, 0.0)

sam = SymbolicAWEModel(set, sys_struct)

init_sim!(sam; remake=false)
for i in 1:80
    plot(sam, i/set.sample_freq)
    next_step!(sam; set_values=[-10.0])
end

# ADDING A KITE
vsm_wing = RamAirWing(set; prn=false)
vsm_aero = BodyAerodynamics([vsm_wing])
vsm_solver = Solver(vsm_aero; solver_type=NONLIN, atol=2e-8, rtol=2e-8)
wings = [SymbolicAWEModels.Wing(1, Group[], I(3), [0.5, 0, set.l_tether+6])]

sys_struct = SystemStructure("wing", set; points, segments, tethers, winches, pulleys, wings, transforms)
plot(sys_struct, 0.0)

sam = SymbolicAWEModel(set, sys_struct, [vsm_aero], [vsm_solver])



