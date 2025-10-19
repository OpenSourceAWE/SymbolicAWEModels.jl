# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

# This file contains VSM-specific extensions to the system structure types defined in KiteUtils.jl

"""
    VortexStepMethod.RamAirWing(set::Settings; prn=true, kwargs...)

Create a `RamAirWing` geometry object from the settings provided.

This is a constructor helper that reads the model and foil file paths from the
`Settings` object and initializes the `RamAirWing` object from `VortexStepMethod.jl`.
"""
function VortexStepMethod.RamAirWing(set::Settings; prn=true, kwargs...)
    obj_path = joinpath(dirname(get_data_path()), set.model)
    dat_path = joinpath(dirname(get_data_path()), set.foil_file)
    if set.physical_model == "simple_ram"
        n_groups=2
    else
        n_groups=4
    end
    return RamAirWing(obj_path, dat_path;
        mass=set.mass, crease_frac=set.crease_frac, n_groups,
        align_to_principal=true, prn, kwargs...
    )
end

"""
    Group(idx, point_idxs, vsm_wing::RamAirWing, gamma, type, moment_frac; damping=50.0)

Constructs a `Group` object from a VortexStepMethod RamAirWing geometry.

This constructor extracts the local chord and spanwise vectors from the VSM wing geometry
at the specified spanwise location `gamma`.

# Arguments
- `idx::Int16`: Unique identifier for the group.
- `point_idxs::Vector{Int16}`: Indices of points that move together with this group's twist.
- `vsm_wing::RamAirWing`: Wing geometry object used to extract local chord and spanwise vectors.
- `gamma`: Spanwise parameter (typically -1 to 1) defining the group's location along the wing.
- `type::DynamicsType`: Dynamics type (`DYNAMIC` for time-varying twist, `QUASI_STATIC` for equilibrium).
- `moment_frac::SimFloat`: Chordwise position (0=leading edge, 1=trailing edge) about which the group rotates.

# Keyword Arguments
- `damping::SimFloat=50.0`: Damping coefficient for twist dynamics.

# Returns
- `Group`: A new `Group` object with twist dynamics capability.
"""
function Group(idx, point_idxs, vsm_wing::RamAirWing, gamma, type, moment_frac; damping=50.0)
    le_pos = [vsm_wing.le_interp[i](gamma) for i in 1:3]
    chord = [vsm_wing.te_interp[i](gamma) for i in 1:3] .- le_pos
    y_airf = normalize([vsm_wing.le_interp[i](gamma-0.01) for i in 1:3] - le_pos)
    Group(idx, point_idxs, le_pos, chord, y_airf, type, moment_frac, damping)
end

"""
    mutable struct VSMWing <: AbstractWing

A wing that uses the Vortex Step Method (VSM) for aerodynamic computations.

This struct extends the base wing functionality with VSM-specific aerodynamic
modeling capabilities, including vortex wake computations and aerodynamic loads.

$(TYPEDFIELDS)
"""
mutable struct VSMWing <: AbstractWing
    # Base wing functionality
    base::BaseWing

    # VSM aerodynamics
    vsm_aero::VortexStepMethod.BodyAerodynamics
    vsm_wing::VortexStepMethod.RamAirWing
    vsm_solver::VortexStepMethod.Solver

    # VSM state and linearization
    const vsm_y::Vector{SimFloat}
    const vsm_x::Vector{SimFloat}
    const vsm_jac::Matrix{SimFloat}

    function VSMWing(base::BaseWing, vsm_aero, vsm_wing, vsm_solver, vsm_y, vsm_x, vsm_jac)
        new(base, vsm_aero, vsm_wing, vsm_solver, vsm_y, vsm_x, vsm_jac)
    end
end

# Delegate property access to base wing for VSMWing
function Base.getproperty(wing::VSMWing, sym::Symbol)
    if sym in (:base, :vsm_aero, :vsm_wing, :vsm_solver, :vsm_y, :vsm_x, :vsm_jac)
        return getfield(wing, sym)
    else
        return getproperty(getfield(wing, :base), sym)
    end
end

function Base.setproperty!(wing::VSMWing, sym::Symbol, value)
    if sym in (:base, :vsm_aero, :vsm_wing, :vsm_solver, :vsm_y, :vsm_x, :vsm_jac)
        setfield!(wing, sym, value)
    else
        setproperty!(getfield(wing, :base), sym, value)
    end
end

"""
    VSMWing(idx::Int16, vsm_aero, vsm_wing, vsm_solver, group_idxs::Vector{Int16},
            R_b_c::Matrix{SimFloat}, pos_cad::KVec3; transform_idx=1, y_damping=150.0)

Constructs a `VSMWing` object with Vortex Step Method aerodynamics.

# Arguments
- `idx::Int16`: Unique identifier for the wing.
- `vsm_aero`: VortexStepMethod.BodyAerodynamics object.
- `vsm_wing`: VortexStepMethod.RamAirWing object.
- `vsm_solver`: VortexStepMethod.Solver object.
- `group_idxs::Vector{Int16}`: Indices of groups attached to this wing.
- `R_b_c::Matrix{SimFloat}`: Rotation matrix from body frame to CAD frame.
- `pos_cad::KVec3`: Position of wing center of mass in CAD frame.

# Keyword Arguments
- `transform_idx::Int16=1`: Transform used for initial positioning and orientation.
- `y_damping::SimFloat=150.0`: Damping coefficient for lateral motion.

# Returns
- `VSMWing`: A new VSM wing object.
"""
function VSMWing(idx::Int, vsm_aero, vsm_wing, vsm_solver,
                 group_idxs::AbstractVector, R_b_c::AbstractMatrix, pos_cad::AbstractVector;
                 transform_idx=1, y_damping=150.0)
    base = BaseWing(Int16(idx), convert(Vector{Int16}, group_idxs), R_b_c, convert(KVec3, pos_cad);
                    transform_idx, y_damping)
    ny = length(group_idxs)+3+3
    nx = length(group_idxs)+3+3
    return VSMWing(base, vsm_aero, vsm_wing, vsm_solver,
                   zeros(SimFloat, ny), zeros(SimFloat, nx), zeros(SimFloat, nx, ny))
end

"""
    Wing(idx, vsm_aero, vsm_wing, vsm_solver, group_idxs, R_b_c, pos_cad; transform_idx)

Constructs a `VSMWing` object (backward compatibility constructor).

This is a convenience constructor that creates a VSMWing for backward compatibility
with existing code. New code should use `VSMWing(...)` directly.

# Arguments
- `idx::Int16`: Unique identifier for the wing.
- `vsm_aero`, `vsm_wing`, `vsm_solver`: Vortex Step Method components.
- `group_idxs::Vector{Int16}`: Indices of groups attached to this wing.
- `R_b_c::Matrix{SimFloat}`: Rotation matrix from body frame to CAD frame.
- `pos_cad::KVec3`: Position of wing center of mass in CAD frame.

# Keyword Arguments
- `transform_idx::Int16=1`: Transform used for initial positioning and orientation.
- `y_damping::SimFloat=150.0`: Damping coefficient for lateral motion.

# Returns
- `VSMWing`: A new VSM wing object.
"""
function Wing(idx, vsm_aero, vsm_wing, vsm_solver, group_idxs, R_b_c, pos_cad; kwargs...)
    return VSMWing(idx, vsm_aero, vsm_wing, vsm_solver, group_idxs, R_b_c, pos_cad; kwargs...)
end

# Helper functions that are VSM-specific

"""
    calc_pos(wing::RamAirWing, gamma, frac)

Calculate a position on the kite based on spanwise (`gamma`) and chordwise (`frac`) parameters.
"""
function calc_pos(wing::RamAirWing, gamma, frac)
    le_pos = [wing.le_interp[i](gamma) for i in 1:3]
    chord = [wing.te_interp[i](gamma) for i in 1:3] .- le_pos
    pos = le_pos .+ chord .* frac
    return pos
end

"""
    cad_to_body_frame(wing::RamAirWing, pos)

Transform a position from the CAD frame to the wing's body frame.
"""
function cad_to_body_frame(wing::RamAirWing, pos)
    return wing.R_cad_body * (pos + wing.T_cad_body)
end

"""
    find_axis_point(P, l, v=[0,0,1])

Calculate the coordinates of a point `Q` that lies on a line defined by vector `v`
and is at a distance `l` from a given point `P`.
"""
function find_axis_point(P, l, v=[0,0,1])
    # Compute discriminant
    D = (v ⋅ P)^2 - norm(v)^2 * (norm(P)^2 - l^2)
    D < 0 && error("No real solution: l is too small or parameters invalid")
    # Solve quadratic for t, choose solution for negative direction
    t = (v ⋅ P - √D) / norm(v)^2
    # Compute point Q = t * v
    return [t * v[1], t * v[2], t * v[3]]
end

"""
    create_tether(tether_idx, set, points, segments, tethers, attach_point, type, dynamics_type; z, axial_stiffness, axial_damping)

Procedurally create a multi-segment tether.

This function builds a tether from a specified number of segments, connecting a given
`attach_point` on the kite to a new anchor point on the ground.
"""
function create_tether(tether_idx, set, points, segments, tethers, attach_point,
                       type, dynamics_type; z=[0,0,1], axial_stiffness=NaN,
                       axial_damping=NaN, d_pos=zeros(3))
    winch_pos = find_axis_point(attach_point.pos_cad, set.l_tether, z) .+ d_pos
    dir = winch_pos - attach_point.pos_cad
    segment_idxs = Int16[]
    winch_idx = 0
    for i in 1:set.segments
        frac = i / set.segments
        pos = attach_point.pos_cad + frac * dir
        point_idx = length(points)+1 # last point idx
        segment_idx = length(segments)+1 # last segment idx
        if i == 1
            last_idx = attach_point.idx
        else
            last_idx = point_idx-1
        end
        if i == set.segments
            points = [points; Point(point_idx, pos, STATIC)]
            winch_idx = points[end].idx
        else
            points = [points; Point(point_idx, pos, dynamics_type)]
        end
        segments = [segments; Segment(segment_idx, set, (last_idx, point_idx), type;
                                      axial_stiffness, axial_damping)]
        push!(segment_idxs, segment_idx)
    end
    tethers = [tethers; Tether(tether_idx, segment_idxs, winch_idx)]
    return points, segments, tethers, tethers[end].idx
end

# Rotation helper functions
function rotate_around_x(vec, angle::T) where T
    result = zeros(T, 3)
    result[1] = vec[1]
    result[2] = cos(angle) * vec[2] - sin(angle) * vec[3]
    result[3] = sin(angle) * vec[2] + cos(angle) * vec[3]
    result
end

function rotate_around_y(vec, angle::T) where T
    result = zeros(T, 3)
    result[1] = cos(angle) * vec[1] + sin(angle) * vec[3]
    result[2] = vec[2]
    result[3] = -sin(angle) * vec[1] + cos(angle) * vec[3]
    result
end

function rotate_around_z(vec, angle::T) where T
    result = zeros(T, 3)
    result[1] = cos(angle) * vec[1] - sin(angle) * vec[2]
    result[2] = sin(angle) * vec[1] + cos(angle) * vec[2]
    result[3] = vec[3]
    result
end

function calc_R_t_w(pos)
    r = norm(pos)
    elevation = asin(pos[3] / r)
    azimuth = atan(pos[2], pos[1])
    R_t_w = zeros(3, 3)
    R_t_w[:, 1] = [cos(elevation) * cos(azimuth), cos(elevation) * sin(azimuth), sin(elevation)]
    R_t_w[:, 2] = [-sin(azimuth), cos(azimuth), 0]
    R_t_w[:, 3] = [-sin(elevation) * cos(azimuth), -sin(elevation) * sin(azimuth), cos(elevation)]
    return R_t_w
end

"""
    get_rot_pos(transform::Transform, wings, points)

Get the position of the rotating object (wing or point) for a given transform.
"""
function get_rot_pos(transform::Transform, wings, points)
    if !isnothing(transform.wing_idx)
        return wings[transform.wing_idx].pos_w
    elseif !isnothing(transform.rot_point_idx)
        return points[transform.rot_point_idx].pos_w
    end
end

"""
    get_base_pos(transform::Transform, wings, points)

Get the base position for a given transform, resolving chained transforms if necessary.
"""
function get_base_pos(transform::Transform, wings, points)
    curr_base_pos = points[transform.base_point_idx].pos_cad
    if !isnothing(transform.base_pos)
        return transform.base_pos, curr_base_pos
    elseif !isnothing(transform.base_transform_idx)
        base_transform = transforms[transform.base_transform_idx]
        return get_rot_pos(base_transform, wings, points), curr_base_pos
    end
end

"""
    apply_heading(vec, R_t_w, curr_R_t_w, heading)

Apply a heading rotation to a vector.
"""
function apply_heading(vec, R_t_w, curr_R_t_w, heading)
    vec_along_z = rotate_around_z(curr_R_t_w' * vec, heading)
    return R_t_w * vec_along_z
end

"""
    reinit!(transforms::Vector{Transform}, sys_struct::SystemStructure)

Apply the initial spatial transformations to all components in a `SystemStructure`.

This function iterates through all transforms and applies the specified translation
and rotation to position and orient the kite system components correctly in the
world frame at the beginning of a simulation.
"""
function reinit!(transforms::Vector{Transform}, sys_struct::SystemStructure)
    @unpack points, wings = sys_struct
    for transform in transforms
        # ==================== TRANSLATE ==================== #
        base_pos, curr_base_pos = get_base_pos(transform, wings, points)
        T = base_pos - curr_base_pos
        for point in points
            if point.transform_idx == transform.idx
                point.pos_w .= point.pos_cad .+ T
                point.vel_w .= 0.0
            end
        end
        for wing in wings
            if wing.transform_idx == transform.idx
                wing.pos_w .= wing.pos_cad .+ T
                wing.vel_w .= 0.0
            end
        end

        # ==================== ROTATE ==================== #
        curr_rot_pos = get_rot_pos(transform, wings, points)
        curr_R_t_w = calc_R_t_w(curr_rot_pos - base_pos)
        transform_pos = rotate_around_z(rotate_around_y([1,0,0], -transform.elevation), -transform.azimuth)
        R_t_w = calc_R_t_w(transform_pos)

        for point in points
            if point.transform_idx == transform.idx
                vec = point.pos_w - base_pos
                point.pos_w .= base_pos + apply_heading(vec, R_t_w, curr_R_t_w, transform.heading)
            end
            if point.type == WING
                wing = wings[point.wing_idx]
                point.pos_b .= wing.R_b_c' * (point.pos_cad - wing.pos_cad) # TODO: test this
            end
        end
        for wing in wings
            if wing.transform_idx == transform.idx
                vec = wing.pos_w - base_pos
                wing.pos_w .= base_pos + apply_heading(vec, R_t_w, curr_R_t_w, transform.heading)
                R_b_w = zeros(3,3)
                for i in 1:3
                    R_b_w[:, i] .= apply_heading(wing.R_b_c[:, i], R_t_w, curr_R_t_w, transform.heading)
                end
                wing.R_b_w = R_b_w
            end
        end
    end
end

"""
    reposition!(transforms::Vector{Transform}, sys_struct::SystemStructure)

Update the system's spatial orientation based on its current position, preserving velocities.

This function adjusts the orientation of all components in the `SystemStructure` without
altering their dynamic state. Unlike `reinit!`, it uses the current world positions (`pos_w`)
as the starting point for rotations, rather than resetting from the CAD coordinates.

This function is useful for making real-time adjustments to the system's pose during a simulation.
Crucially, it **preserves the existing velocities (`vel_w`) of all points and wings**.

NOTE: the transform.heading is applied relative to the current heading of the system.

# Arguments
- `sys_struct::SystemStructure`: The system model to update.
"""
function reposition!(transforms::Vector{Transform}, sys_struct::SystemStructure)
    @unpack points, wings = sys_struct
    for transform in transforms
        # Get the current positions of the base and the rotating object
        base_pos = points[transform.base_point_idx].pos_w
        rot_pos = get_rot_pos(transform, wings, points)

        # Calculate the current orientation in spherical coordinates
        curr_rel_pos = rot_pos - base_pos
        curr_R_t_w = calc_R_t_w(curr_rel_pos)
        transform_pos = rotate_around_z(rotate_around_y([1,0,0], -transform.elevation), -transform.azimuth)
        R_t_w = calc_R_t_w(transform_pos)

        # Apply the rotation to all relevant points
        for point in points
            if point.transform_idx == transform.idx
                vec = point.pos_w - base_pos
                point.pos_w .= base_pos + apply_heading(vec, R_t_w, curr_R_t_w, transform.heading)
            end
        end

        # Apply the rotation to all relevant wings
        for wing in wings
            if wing.transform_idx == transform.idx
                # Rotate the wing's position
                vec = wing.pos_w - base_pos
                wing.pos_w .= base_pos + apply_heading(vec, R_t_w, curr_R_t_w, transform.heading)

                # Rotate the wing's orientation matrix
                R_b_w = zeros(3,3)
                current_R_b_w = wing.R_b_w
                for i in 1:3
                    R_b_w[:, i] .= apply_heading(current_R_b_w[:, i], R_t_w, curr_R_t_w, transform.heading)
                end
                wing.R_b_w = R_b_w
            end
        end
    end
end

"""
    reinit!(sys_struct::SystemStructure, set::Settings)

Re-initialize a `SystemStructure` from a `Settings` object.

This function resets various component states (e.g., winch lengths, group twists,
pulley positions) to their initial values as defined in the `Settings` object. It
is typically called before starting a new simulation run.
"""
function reinit!(sys_struct::SystemStructure, set::Settings)
    @unpack points, groups, segments, pulleys, tethers, winches, wings, transforms = sys_struct

    for segment in segments
        @assert (0 < segment.diameter < 1)
    end

    for winch in winches
        winch.tether_len = set.l_tethers[winch.idx]
        winch.tether_vel    = set.v_reel_outs[winch.idx]
    end

    (length(groups) > 0) && (first_moment_frac = groups[1].moment_frac)
    for group in groups
        group.twist = 0.0
        group.twist_ω = 0.0
        @assert group.moment_frac ≈ first_moment_frac "All group.moment_frac must be the same."
    end

    for transform in transforms
        transform.elevation = deg2rad(set.elevations[transform.idx])
        transform.azimuth   = deg2rad(set.azimuths[transform.idx])
        transform.heading   = deg2rad(set.headings[transform.idx])
    end

    for segment in segments
        (segment.l0 ≈ 0) && (segment.l0 = norm(points[segment.point_idxs[1]].pos_cad - points[segment.point_idxs[2]].pos_cad))
        @assert (segment.l0 > 0)
    end

    for pulley in pulleys
        segment1, segment2 = segments[pulley.segment_idxs[1]], segments[pulley.segment_idxs[2]]
        pulley.sum_len = segment1.l0 + segment2.l0
        pulley.len = segment1.l0
        pulley.vel = 0.0
        @assert !(pulley.sum_len ≈ 0)
    end

    reinit!(transforms, sys_struct)
    for wing in wings
        if isa(wing, VSMWing)
            wing.vsm_y .= 0.0
            wing.vsm_y[1:3] .= wing.R_b_w' * [set.v_wind, 0., 0.]
        end
    end

    return nothing
end

"""
    copy!(sys1::SystemStructure, sys2::SystemStructure)

Copy the dynamic state from one `SystemStructure` (`sys1`) to another (`sys2`).

This function is designed to transfer the state (positions, velocities, etc.) between
two system models, which can have different levels of fidelity. For example, it can
copy the state from a detailed multi-segment tether model (`sys1`) to a simplified
single-segment model (`sys2`).

The function handles several cases:
- If `sys1` and `sys2` have the same structure, it performs a direct copy of all point states.
- If `sys2` is a simplified (1-segment per tether) version of `sys1`, it copies the
  positions and velocities of the tether endpoints.
- It also copies the state of wings, groups, winches, and pulleys where applicable.
"""
function copy!(sys1::SystemStructure, sys2::SystemStructure)
    simple = false

    # copy point pos and vel
    if length(sys1.points) > 0
        if length(sys1.points) == length(sys2.points)
            for (point1, point2) in zip(sys1.points, sys2.points)
                point2.pos_w .= point1.pos_w
                point2.vel_w .= point1.vel_w
                point2.disturb .= point1.disturb
            end
        # if different number of points, copy only the tether points
        elseif length(sys1.tethers) > 1 && length(sys1.tethers) == length(sys2.tethers)
            for (tether1, tether2) in zip(sys1.tethers, sys2.tethers)
                if length(tether1.segment_idxs) == length(tether2.segment_idxs)
                    # copy the points of the segments of the tethers
                    for (segment_idx1, segment_idx2) in zip(tether1.segment_idxs, tether2.segment_idxs)
                        point_idxs1 = sys1.segments[segment_idx1].point_idxs
                        point_idxs2 = sys2.segments[segment_idx2].point_idxs
                        for (point_idx1, point_idx2) in zip(point_idxs1, point_idxs2)
                            sys2.points[point_idx2].pos_w .= sys1.points[point_idx1].pos_w
                            sys2.points[point_idx2].vel_w .= sys1.points[point_idx1].vel_w
                            sys2.points[point_idx2].disturb .= sys1.points[point_idx1].disturb
                        end
                    end
                elseif length(tether2.segment_idxs) == 1
                    # copy the first and last point of the tether
                    point_idxs1 = [sys1.segments[tether1.segment_idxs[1]].point_idxs[1],
                                   sys1.segments[tether1.segment_idxs[end]].point_idxs[2]]
                    point_idxs2 = sys2.segments[tether2.segment_idxs[1]].point_idxs
                    sys2.points[point_idxs2[1]].pos_w .= sys1.points[point_idxs1[1]].pos_w
                    sys2.points[point_idxs2[2]].pos_w .= sys1.points[point_idxs1[2]].pos_w
                    sys2.points[point_idxs2[1]].vel_w .= sys1.points[point_idxs1[1]].vel_w
                    sys2.points[point_idxs2[2]].vel_w .= sys1.points[point_idxs1[2]].vel_w
                    sys2.points[point_idxs2[1]].disturb .= sys1.points[point_idxs1[1]].disturb
                    sys2.points[point_idxs2[2]].disturb .= sys1.points[point_idxs1[2]].disturb
                    simple = true
                end
            end
        end
    end

    # copy twist and twist_ω of groups
    if length(sys1.groups) > 1 && length(sys1.groups) == length(sys2.groups)
        for (group1, group2) in zip(sys1.groups, sys2.groups)
            group2.twist = group1.twist
            group2.twist_ω = group1.twist_ω
        end
    end

    # copy winch tether lengths and velocities
    if length(sys1.winches) > 1 && length(sys1.winches) == length(sys2.winches)
        for (winch2, winch1) in zip(sys2.winches, sys1.winches)
            if !simple
                winch2.tether_len = winch1.tether_len
                winch2.tether_vel = winch1.tether_vel
            else
                winch2.tether_len = 0.0
                for tether_idx in winch1.tether_idxs
                    tether2 = sys2.tethers[tether_idx]
                    segment2 = sys2.segments[tether2.segment_idxs[1]]
                    point_idxs2 = segment2.point_idxs
                    slen = norm(sys2.points[point_idxs2[1]].pos_w .-
                                        sys2.points[point_idxs2[2]].pos_w)
                    stiffness = segment2.axial_stiffness / slen
                    nt = length(winch1.tether_idxs)
                    winch2.tether_len += (slen - norm(winch1.force)/stiffness/nt) / nt
                end
                winch2.tether_vel = winch1.tether_vel
            end
        end
    end

    # copy pulley lengths and velocities
    if length(sys1.pulleys) > 1 && length(sys1.pulleys) == length(sys2.pulleys)
        for (pulley1, pulley2) in zip(sys1.pulleys, sys2.pulleys)
            pulley2.len = pulley1.len
            pulley2.vel = pulley1.vel
        end
    end

    # copy wing positions and velocities
    if length(sys1.wings) > 1 && length(sys1.wings) == length(sys2.wings)
        for (wing1, wing2) in zip(sys1.wings, sys2.wings)
            wing2.pos_w .= wing1.pos_w
            wing2.vel_w .= wing1.vel_w
            wing2.ω_b .= wing1.ω_b
            wing2.Q_b_w .= wing1.Q_b_w
        end
    end
end

"""
    update_from_sysstate!(sys::SystemStructure, ss::SysState)

Update the dynamic state of a `SystemStructure` from a `SysState` snapshot.

This function copies the state variables that are present in `SysState` (such as point
positions, wing orientations, winch lengths, and twist angles) into an existing `SystemStructure`.
Fields that cannot be populated from `SysState` (such as aerodynamic forces, moments, and
segment forces) are set to `NaN` to prevent them from being plotted.

This is useful for visualizing a `SysLog` by extracting individual `SysState` snapshots
and applying them to a `SystemStructure` for plotting with the Makie extension.

# Arguments
- `sys::SystemStructure`: The system structure to update (must already exist with correct topology).
- `ss::SysState`: The state snapshot to copy from.

# Example
```julia
# Load a system log
lg = load_log(...)

# Create a SystemStructure with the same topology
sys = SystemStructure(se(), "ram")

# Update from a specific time step
update_from_sysstate!(sys, lg.syslog[100])

# Plot the system at that time step
plot(sys)
```

# Notes
- The `SystemStructure` must have been created with the same model configuration as the
  simulation that generated the `SysLog`.
- Aerodynamic and force fields are set to `NaN` and will not be plotted.
- The number of points in `sys` must match the parametric type `P` of `SysState{P}`.
"""
function update_from_sysstate!(sys::SystemStructure, ss::SysState{P}) where P
    @unpack points, groups, winches, wings = sys

    # Verify compatibility
    if length(points) != P
        error("SystemStructure has $(length(points)) points but SysState has $P points")
    end

    # Update point positions (X, Y, Z from SysState)
    for point in points
        point.pos_w[1] = ss.X[point.idx]
        point.pos_w[2] = ss.Y[point.idx]
        point.pos_w[3] = ss.Z[point.idx]
        # Set velocity to zero (not available in basic SysState)
        point.vel_w .= 0.0
        # Set forces to NaN (not available in SysState)
        point.force .= NaN
    end

    # Update wing state if wings exist
    if length(wings) > 0 && length(wings) == 1  # Currently only support single-wing systems
        wing = wings[1]

        # Copy orientation quaternion
        wing.Q_b_w .= ss.orient

        # Copy spherical coordinates
        wing.elevation = Float64(ss.elevation)
        wing.azimuth = Float64(ss.azimuth)
        wing.heading = Float64(ss.heading)

        # Compute wing position from average of points (if wings exist)
        # For a typical system, the wing COM is near the bridle attachment
        wing.pos_w .= [mean(ss.X), mean(ss.Y), mean(ss.Z)]

        # Copy velocity if available in vel_kite
        wing.vel_w .= ss.vel_kite

        # Set angular velocity to NaN (turn_rates in SysState, but need conversion)
        wing.ω_b .= ss.turn_rates

        # Set aerodynamic quantities to NaN (to prevent plotting)
        wing.aero_force_b .= NaN
        wing.aero_moment_b .= NaN
        wing.tether_force .= NaN
        wing.tether_moment .= NaN
        wing.va_b .= NaN
        wing.v_wind .= ss.v_wind_kite
        wing.aoa = Float64(ss.AoA)
        wing.course = Float64(ss.course)
        wing.acc_w .= 0.0
        wing.turn_rate .= ss.turn_rates
        wing.turn_acc .= 0.0
    end

    # Update group twist angles
    n_groups = min(length(groups), 4)  # SysState stores up to 4 twist angles
    for i in 1:n_groups
        if i <= length(groups)
            groups[i].twist = Float64(ss.twist_angles[i])
            groups[i].twist_ω = 0.0  # Not available in SysState
            # Set forces/moments to NaN
            groups[i].tether_force = NaN
            groups[i].tether_moment = NaN
            groups[i].aero_moment = NaN
        end
    end

    # Update winch state
    n_winches = min(length(winches), 4)  # SysState stores up to 4 winches
    for i in 1:n_winches
        if i <= length(winches)
            winches[i].tether_len = Float64(ss.l_tether[i])
            winches[i].tether_vel = Float64(ss.v_reelout[i])
            winches[i].force .= NaN  # Force not directly available
            winches[i].friction = NaN
            winches[i].tether_acc = 0.0
            winches[i].set_value = Float64(ss.set_torque[i])
        end
    end

    # Update segments - set forces to NaN (not available in SysState)
    for segment in sys.segments
        p1 = points[segment.point_idxs[1]]
        p2 = points[segment.point_idxs[2]]
        segment.len = norm(p1.pos_w - p2.pos_w)
        segment.force = NaN  # Not available in SysState
    end

    # Update global wind vector
    sys.wind_vec_gnd .= ss.v_wind_gnd

    return nothing
end
