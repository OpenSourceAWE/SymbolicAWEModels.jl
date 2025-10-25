# SPDX-FileCopyrightText: 2025 Bart van de Lint <bart@vandelint.net>
#
# SPDX-License-Identifier: MPL-2.0

"""
    SystemStructure(set::Settings; kwargs...)

Factory function to create a `SystemStructure` for a specific physical model.

This function acts as a dispatcher, calling the appropriate `create_*_sys_struct`
function based on the `physical_model` field in the `Settings` object.

# Arguments
- `set::Settings`: The settings object that defines which model to create.
- `kwargs...`: Keyword arguments passed through to the specific constructor.

# Returns
- `SystemStructure`: The fully constructed system model.
"""
function KiteUtils.SystemStructure(set::Settings; kwargs...)
    func_name = Symbol("create_$(set.physical_model)_sys_struct")
    return getfield(SymbolicAWEModels, func_name)(set; kwargs...)
end

"""
    create_4_attach_ram_sys_struct(set::Settings)

Create a detailed `SystemStructure` for a ram-air kite with a 4-point attachment bridle.

This function procedurally builds a complex kite model. Its key feature is that all four
bridle attachment points on each of the four wing `Group` sections are modeled as
deforming with the group's twist dynamics.

This model includes:
- A flexible wing simulated with 4 deformable `Group` sections.
- A detailed bridle system with multiple segments and `Pulley`s to distribute forces.
- Four main tethers (left/right power, left/right steering) connecting the bridle
  to the ground winches.
- A 3-winch system controlling the tethers.

# Arguments
- `set::Settings`: Configuration parameters defining the kite's geometry and properties.

# Returns
- `SystemStructure`: A new `SystemStructure` object representing the detailed model.
"""
function create_4_attach_ram_sys_struct(set::Settings)
    vsm_wing = RamAirWing(set; prn=false)
    points = Point[]
    groups = Group[]
    segments = Segment[]
    pulleys = Pulley[]
    tethers = Tether[]
    winches = Winch[]
    wings = Wing[]

    attach_points = Point[]
    
    bridle_top_left = [cad_to_body_frame(vsm_wing, set.top_bridle_points[i]) for i in eachindex(set.top_bridle_points)]
    bridle_top_right = [bridle_top_left[i] .* [1, -1, 1] for i in eachindex(set.top_bridle_points)]

    dynamics_type = set.quasi_static ? QUASI_STATIC : DYNAMIC
    z = vsm_wing.R_cad_body[:,3]

    function create_bridle(bridle_top, gammas)
        i_pnt = length(points) # last point idx
        i_seg = length(segments) # last segment idx
        i_pul = length(pulleys) # last pulley idx
        i_grp = length(groups) # last group idx

        # ==================== CREATE DEFORMING WING GROUPS ==================== #
        points = [
            points
            Point(1+i_pnt, calc_pos(vsm_wing, gammas[1], set.bridle_fracs[1]), WING)
            Point(2+i_pnt, calc_pos(vsm_wing, gammas[1], set.bridle_fracs[2]), WING)
            Point(3+i_pnt, calc_pos(vsm_wing, gammas[1], set.bridle_fracs[3]), WING)
            Point(4+i_pnt, calc_pos(vsm_wing, gammas[1], set.bridle_fracs[4]), WING)

            Point(5+i_pnt, calc_pos(vsm_wing, gammas[2], set.bridle_fracs[1]), WING)
            Point(6+i_pnt, calc_pos(vsm_wing, gammas[2], set.bridle_fracs[2]), WING)
            Point(7+i_pnt, calc_pos(vsm_wing, gammas[2], set.bridle_fracs[3]), WING)
            Point(8+i_pnt, calc_pos(vsm_wing, gammas[2], set.bridle_fracs[4]), WING)
        ]
        groups = [
            groups
            Group(1+i_grp, [1+i_pnt, 2+i_pnt, 3+i_pnt, 4+i_pnt], vsm_wing, gammas[1], DYNAMIC, set.bridle_fracs[2])
            Group(2+i_grp, [5+i_pnt, 6+i_pnt, 7+i_pnt, 8+i_pnt], vsm_wing, gammas[2], DYNAMIC, set.bridle_fracs[2])
        ]

        # ==================== CREATE PULLEY BRIDLE SYSTEM ==================== #
        bridle_damping = 1.0
        points = [
            points
            Point(9+i_pnt, bridle_top[1], dynamics_type; bridle_damping)
            Point(10+i_pnt, bridle_top[2], dynamics_type; bridle_damping)
            Point(11+i_pnt, bridle_top[3], dynamics_type; bridle_damping)
            Point(12+i_pnt, bridle_top[4], dynamics_type; bridle_damping)

            Point(13+i_pnt, bridle_top[2] - 1z, dynamics_type; bridle_damping)

            Point(14+i_pnt, bridle_top[1] - 2z, dynamics_type; bridle_damping)
            Point(15+i_pnt, bridle_top[3] - 2z, dynamics_type; bridle_damping)

            Point(16+i_pnt, bridle_top[1] - 4z, dynamics_type; bridle_damping)
            Point(17+i_pnt, bridle_top[3] - 4z, dynamics_type; bridle_damping)
        ]
        segments = [
            segments
            Segment(1+i_seg, set, (1+i_pnt, 9+i_pnt), BRIDLE)
            Segment(2+i_seg, set, (2+i_pnt, 10+i_pnt), BRIDLE)
            Segment(3+i_seg, set, (3+i_pnt, 11+i_pnt), BRIDLE)
            Segment(4+i_seg, set, (4+i_pnt, 12+i_pnt), BRIDLE)

            Segment(5+i_seg, set, (5+i_pnt, 9+i_pnt), BRIDLE)
            Segment(6+i_seg, set, (6+i_pnt, 10+i_pnt), BRIDLE)
            Segment(7+i_seg, set, (7+i_pnt, 11+i_pnt), BRIDLE)
            Segment(8+i_seg, set, (8+i_pnt, 12+i_pnt), BRIDLE)

            Segment(9+i_seg, set, (9+i_pnt, 14+i_pnt), BRIDLE; l0=2)
            Segment(10+i_seg, set, (10+i_pnt, 13+i_pnt), BRIDLE; l0=1)
            Segment(11+i_seg, set, (11+i_pnt, 15+i_pnt), BRIDLE; l0=2)
            Segment(12+i_seg, set, (12+i_pnt, 17+i_pnt), BRIDLE; l0=4)
            
            Segment(13+i_seg, set, (13+i_pnt, 14+i_pnt), BRIDLE; l0=1)
            Segment(14+i_seg, set, (13+i_pnt, 15+i_pnt), BRIDLE; l0=1)
            
            Segment(15+i_seg, set, (14+i_pnt, 16+i_pnt), BRIDLE; l0=2)
            Segment(16+i_seg, set, (15+i_pnt, 16+i_pnt), BRIDLE; l0=2)
            Segment(17+i_seg, set, (15+i_pnt, 17+i_pnt), BRIDLE; l0=2)
        ]
        pulleys = [
            pulleys
            Pulley(1+i_pul, (13+i_seg, 14+i_seg), dynamics_type)
            Pulley(2+i_pul, (16+i_seg, 17+i_seg), dynamics_type)
        ]
        push!(attach_points, points[end-1])
        push!(attach_points, points[end])
        return nothing
    end

    gammas = [-3/4, -1/4, 1/4, 3/4] * vsm_wing.gamma_tip
    create_bridle(bridle_top_left, gammas[[1,2]])
    create_bridle(bridle_top_right, gammas[[3,4]])

    points, segments, tethers, left_power_idx =
        create_tether(1, set, points, segments, tethers, attach_points[1], 
                      POWER_LINE, dynamics_type; z)
    points, segments, tethers, right_power_idx =
        create_tether(2, set, points, segments, tethers, attach_points[3], 
                      POWER_LINE, dynamics_type; z)
    points, segments, tethers, left_steering_idx =
        create_tether(3, set, points, segments, tethers, attach_points[2], 
                      STEERING_LINE, dynamics_type; z)
    points, segments, tethers, right_steering_idx =
        create_tether(4, set, points, segments, tethers, attach_points[4], 
                      STEERING_LINE, dynamics_type; z)

    winches = [
        Winch(1, set, [left_power_idx, right_power_idx])
        Winch(2, set, [left_steering_idx])
        Winch(3, set, [right_steering_idx])
    ]

    vsm_aero = BodyAerodynamics([vsm_wing])
    vsm_solver = Solver(vsm_aero; solver_type=NONLIN, atol=2e-8, rtol=2e-8)
    wings = [VSMWing(1, vsm_aero, vsm_wing, vsm_solver, [1,2,3,4], I(3), zeros(3))]
    transforms = [Transform(1, deg2rad(set.elevation), deg2rad(set.azimuth), deg2rad(set.heading);
                            base_pos= zeros(3), base_point_idx=points[end].idx, wing_idx=1)]
    
    return SystemStructure(set.physical_model, set; points, groups, segments, pulleys, tethers, winches, wings, transforms)
end

"""
    create_ram_sys_struct(set::Settings)

Create a `SystemStructure` for the primary "ram" model with a stability-enhancing bridle.

This is the main, detailed model configuration. It differs from the `4_attach` version
by having each of its four `Group` sections defined by three deforming bridle points and
one statically attached (non-deforming) point. This design improves stability
by ensuring the kite's z-axis remains aligned with the bridle system.

The model features:
- A flexible wing with 4 deformable groups (3 deforming points + 1 static point each).
- A complex bridle system with pulleys.
- Four main tethers and three winches.

# Arguments
- `set::Settings`: Configuration parameters defining the kite's geometry and properties.

# Returns
- `SystemStructure`: A new `SystemStructure` object representing the "ram" model.
"""
function create_ram_sys_struct(set::Settings; d_winch_pos=[zeros(3), zeros(3)])
    vsm_wing = RamAirWing(set; prn=false)
    points = Point[]
    groups = Group[]
    segments = Segment[]
    pulleys = Pulley[]
    tethers = Tether[]
    winches = Winch[]
    wings = Wing[]

    attach_points = Point[]
    
    bridle_top_left = [cad_to_body_frame(vsm_wing, set.top_bridle_points[i]) for i in eachindex(set.top_bridle_points)]
    bridle_top_right = [bridle_top_left[i] .* [1, -1, 1] for i in eachindex(set.top_bridle_points)]

    dynamics_type = set.quasi_static ? QUASI_STATIC : DYNAMIC
    z = vsm_wing.R_cad_body[:,3]

    function create_bridle(bridle_top, gammas)
        i_pnt = length(points) # last point idx
        i_seg = length(segments) # last segment idx
        i_pul = length(pulleys) # last pulley idx
        i_grp = length(groups) # last group idx

        # ==================== CREATE DEFORMING WING GROUPS ==================== #
        points = [
            points
            Point(1+i_pnt, calc_pos(vsm_wing, gammas[1], set.bridle_fracs[1]), WING)
            Point(2+i_pnt, calc_pos(vsm_wing, gammas[1], set.bridle_fracs[3]), WING)
            Point(3+i_pnt, calc_pos(vsm_wing, gammas[1], set.bridle_fracs[4]), WING)

            Point(4+i_pnt, calc_pos(vsm_wing, gammas[2], set.bridle_fracs[1]), WING)
            Point(5+i_pnt, calc_pos(vsm_wing, gammas[2], set.bridle_fracs[3]), WING)
            Point(6+i_pnt, calc_pos(vsm_wing, gammas[2], set.bridle_fracs[4]), WING)
        ]
        groups = [
            groups
            Group(1+i_grp, [1+i_pnt, 2+i_pnt, 3+i_pnt], vsm_wing, gammas[1], DYNAMIC, 0.25)
            Group(2+i_grp, [4+i_pnt, 5+i_pnt, 6+i_pnt], vsm_wing, gammas[2], DYNAMIC, 0.25)
        ]

        # ==================== CREATE PULLEY BRIDLE SYSTEM ==================== #
        bridle_damping = 1.0
        points = [
            points
            Point(7+i_pnt, bridle_top[1], dynamics_type; bridle_damping)
            Point(8+i_pnt, bridle_top[2], WING)
            Point(9+i_pnt, bridle_top[3], dynamics_type; bridle_damping)
            Point(10+i_pnt, bridle_top[4], dynamics_type; bridle_damping)

            Point(11+i_pnt, bridle_top[2] - 1z, dynamics_type; bridle_damping)

            Point(12+i_pnt, bridle_top[1] - 2z, dynamics_type; bridle_damping)
            Point(13+i_pnt, bridle_top[3] - 2z, dynamics_type; bridle_damping)

            Point(14+i_pnt, bridle_top[1] - 4z, dynamics_type; bridle_damping)
            Point(15+i_pnt, bridle_top[3] - 4z, dynamics_type; bridle_damping)
        ]
        segments = [
            segments
            Segment(1+i_seg, set, (1+i_pnt, 7+i_pnt), BRIDLE)
            Segment(2+i_seg, set, (2+i_pnt, 9+i_pnt), BRIDLE)
            Segment(3+i_seg, set, (3+i_pnt, 10+i_pnt), BRIDLE)

            Segment(4+i_seg, set, (4+i_pnt, 7+i_pnt), BRIDLE)
            Segment(5+i_seg, set, (5+i_pnt, 9+i_pnt), BRIDLE)
            Segment(6+i_seg, set, (6+i_pnt, 10+i_pnt), BRIDLE)

            Segment(7+i_seg, set, (7+i_pnt, 12+i_pnt), BRIDLE; l0=2)
            Segment(8+i_seg, set, (8+i_pnt, 11+i_pnt), BRIDLE; l0=1)
            Segment(9+i_seg, set, (9+i_pnt, 13+i_pnt), BRIDLE; l0=2)
            Segment(10+i_seg, set, (10+i_pnt, 15+i_pnt), BRIDLE; l0=4)
            
            Segment(11+i_seg, set, (11+i_pnt, 12+i_pnt), BRIDLE; l0=1)
            Segment(12+i_seg, set, (11+i_pnt, 13+i_pnt), BRIDLE; l0=1)
            
            Segment(13+i_seg, set, (12+i_pnt, 14+i_pnt), BRIDLE; l0=2)
            Segment(14+i_seg, set, (13+i_pnt, 14+i_pnt), BRIDLE; l0=2)
            Segment(15+i_seg, set, (13+i_pnt, 15+i_pnt), BRIDLE; l0=2)
        ]
        pulleys = [
            pulleys
            Pulley(1+i_pul, (11+i_seg, 12+i_seg), dynamics_type)
            Pulley(2+i_pul, (14+i_seg, 15+i_seg), dynamics_type)
        ]
        push!(attach_points, points[end-1])
        push!(attach_points, points[end])
        return nothing
    end

    gammas = [-3/4, -1/4, 1/4, 3/4] * vsm_wing.gamma_tip
    create_bridle(bridle_top_left, gammas[[1,2]])
    create_bridle(bridle_top_right, gammas[[3,4]])

    points, segments, tethers, left_power_idx =
        create_tether(1, set, points, segments, tethers, attach_points[1], 
                      POWER_LINE, dynamics_type; z)
    points, segments, tethers, right_power_idx =
        create_tether(2, set, points, segments, tethers, attach_points[3], 
                      POWER_LINE, dynamics_type; z)
    points, segments, tethers, left_steering_idx =
        create_tether(3, set, points, segments, tethers, attach_points[2], 
                      STEERING_LINE, dynamics_type; z, d_pos=d_winch_pos[1])
    points, segments, tethers, right_steering_idx =
        create_tether(4, set, points, segments, tethers, attach_points[4], 
                      STEERING_LINE, dynamics_type; z, d_pos=d_winch_pos[2])

    winches = [
        Winch(1, set, [left_power_idx, right_power_idx])
        Winch(2, set, [left_steering_idx])
        Winch(3, set, [right_steering_idx])
    ]

    vsm_aero = BodyAerodynamics([vsm_wing])
    vsm_solver = Solver(vsm_aero; solver_type=NONLIN, atol=2e-8, rtol=2e-8)
    wings = [VSMWing(1, vsm_aero, vsm_wing, vsm_solver, [1,2,3,4], I(3), zeros(3))]
    transforms = [Transform(1, deg2rad(set.elevation), deg2rad(set.azimuth), deg2rad(set.heading);
                                      base_pos= zeros(3), base_point_idx=points[end].idx, wing_idx=1)]
    
    return SystemStructure(set.physical_model, set; points, groups, segments, pulleys, tethers, winches, wings, transforms)
end

"""
    create_tether_sys_struct(set::Settings; axial_stiffness, axial_damping)

Create a simplified `SystemStructure` for testing tether dynamics.

This model consists of only four independent tethers, each represented by a dynamic
point mass connected to a fixed ground anchor. It does not include a wing or bridle
system, making it ideal for isolating and analyzing the behavior of the tethers
themselves.

# Arguments
- `set::Settings`: Configuration parameters.

# Keywords
- `axial_stiffness::Vector{Float64}`: Predefined axial stiffness [N] for each tether.
- `axial_damping::Vector{Float64}`: Predefined axial damping [Ns] for each tether.

# Returns
- `SystemStructure`: A new `SystemStructure` representing the 4-tether test system.
"""
function create_tether_sys_struct(set::Settings; 
                                  axial_stiffness=fill(NaN, 4), 
                                  axial_damping=fill(NaN,4))
    points = Point[]
    segments = Segment[]
    tethers = Tether[]
    
    points = [
        Point(1, zeros(3), DYNAMIC; fix_sphere=true)
        Point(2, zeros(3), DYNAMIC; fix_sphere=true)
        Point(3, zeros(3), DYNAMIC; fix_sphere=true)
        Point(4, zeros(3), DYNAMIC; fix_sphere=true)
    ]
    
    points, segments, tethers, left_power_idx =
        create_tether(1, set, points, segments, tethers, points[1], POWER_LINE, DYNAMIC, 
                      axial_stiffness=axial_stiffness[1], axial_damping=axial_damping[1])
    points, segments, tethers, right_power_idx =
        create_tether(2, set, points, segments, tethers, points[2], POWER_LINE, DYNAMIC, 
                      axial_stiffness=axial_stiffness[2], axial_damping=axial_damping[2])
    points, segments, tethers, left_steering_idx =
        create_tether(3, set, points, segments, tethers, points[3], STEERING_LINE, DYNAMIC, 
                      axial_stiffness=axial_stiffness[3], axial_damping=axial_damping[3])
    points, segments, tethers, right_steering_idx =
        create_tether(4, set, points, segments, tethers, points[4], STEERING_LINE, DYNAMIC, 
                      axial_stiffness=axial_stiffness[4], axial_damping=axial_damping[4])
    
    winches = [
        Winch(1, set, [left_power_idx, right_power_idx]; brake=true)
        Winch(2, set, [left_steering_idx]; brake=true)
        Winch(3, set, [right_steering_idx]; brake=true)
    ]
    
    transforms = [Transform(1, deg2rad(set.elevation), deg2rad(set.azimuth), deg2rad(set.heading);
                                      base_pos=zeros(3), base_point_idx=points[end].idx, rot_point_idx=1)]

    return SystemStructure("tether", set; points, segments, tethers, winches, transforms)
end

"""
    create_simple_ram_sys_struct(set::Settings; axial_stiffness, axial_damping)

Create a simplified `SystemStructure` for a ram-air kite with direct tether connections.

This model represents a kite with a flexible wing (2 deformable groups) but simplifies
the bridle by connecting the four main tethers directly to the wing attachment points,
omitting the complex pulley system. Each tether is modeled as a single segment.

# Arguments
- `set::Settings`: Configuration parameters.

# Keywords
- `axial_stiffness::Vector{Float64}`: Predefined axial stiffness [N] for each tether.
- `axial_damping::Vector{Float64}`: Predefined axial damping [Ns] for each tether.

# Returns
- `SystemStructure`: A new `SystemStructure` representing the simplified model.
"""
function create_simple_ram_sys_struct(set::Settings; 
                                      axial_stiffness=fill(NaN, 4), 
                                      axial_damping=fill(NaN,4))
    set.segments = 1
    vsm_wing = RamAirWing(set; prn=false)
    gammas = [-1/2, 1/2] * vsm_wing.gamma_tip
    
    bridle_top_left = [vsm_wing.R_cad_body * (set.top_bridle_points[i] + vsm_wing.T_cad_body) for i in eachindex(set.top_bridle_points)] # cad to kite frame
    bridle_top_right = [bridle_top_left[i] .* [1, -1, 1] for i in eachindex(set.top_bridle_points)]

    points = [
        Point(1, bridle_top_left[2], WING)
        Point(2, bridle_top_right[2], WING)
        Point(3, calc_pos(vsm_wing, gammas[1], set.bridle_fracs[4]), WING)
        Point(4, calc_pos(vsm_wing, gammas[2], set.bridle_fracs[4]), WING)

        Point(5, [0, 0, -set.l_tether], STATIC)
        Point(6, [0, 0, -set.l_tether], STATIC)
        Point(7, [0, 0, -set.l_tether], STATIC)
        Point(8, [0, 0, -set.l_tether], STATIC)
    ]
    groups = [
        Group(1, [3], vsm_wing, gammas[1], DYNAMIC, 0.25)
        Group(2, [4], vsm_wing, gammas[2], DYNAMIC, 0.25)
    ]
    segments = [
        Segment(1, set, (1,5), POWER_LINE; axial_stiffness=axial_stiffness[1],
                axial_damping=axial_damping[1])
        Segment(2, set, (2,6), POWER_LINE; axial_stiffness=axial_stiffness[2],
                axial_damping=axial_damping[2])
        Segment(3, set, (3,7), STEERING_LINE; axial_stiffness=axial_stiffness[3],
                axial_damping=axial_damping[3])
        Segment(4, set, (4,8), STEERING_LINE; axial_stiffness=axial_stiffness[4],
                axial_damping=axial_damping[4])
    ]
    tethers = [
        Tether(1, [1], 5)
        Tether(2, [2], 6)
        Tether(3, [3], 7)
        Tether(4, [4], 8)
    ]
    winches = [
        Winch(1, set, [1,2])
        Winch(2, set, [3])
        Winch(3, set, [4])
    ]
    vsm_aero = BodyAerodynamics([vsm_wing])
    vsm_solver = Solver(vsm_aero; solver_type=NONLIN, atol=2e-8, rtol=2e-8)
    wings = [VSMWing(1, vsm_aero, vsm_wing, vsm_solver, [1,2], I(3), zeros(3))]
    transforms = [
        Transform(1, deg2rad(set.elevation), deg2rad(set.azimuth), deg2rad(set.heading);
                base_pos=zeros(3), base_point_idx=5, wing_idx=1)
    ]

    return SystemStructure(set.physical_model, set; 
        points, groups, segments, tethers, winches, wings, transforms)
end
