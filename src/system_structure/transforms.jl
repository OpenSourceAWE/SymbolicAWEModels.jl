# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

"""
Transform functions for heading calculation and spatial positioning.

This file contains:
- Heading calculation functions (calc_heading, apply_heading, etc.)
- reinit! and reposition! functions for applying transforms

Note: The Transform struct and its constructors are defined in types.jl
"""

# ==================== HEADING CALCULATION ==================== #

"""
    apply_heading(vec, R_t_to_w, curr_R_t_to_w, heading)

Apply a heading rotation to a vector.
"""
function apply_heading(vec, R_t_to_w, curr_R_t_to_w, heading)
    vec_along_z = rotate_around_z(curr_R_t_to_w' * vec, heading)
    return R_t_to_w * vec_along_z
end

"""
    wrap_to_pi(angle)

Wrap angle to [-π, π] range.
"""
function wrap_to_pi(angle)
    return mod(angle + π, 2π) - π
end

"""
    calc_heading(R_b_to_w, wind_norm)

Calculate heading angle from body-to-world rotation matrix and wind direction.
Heading is the angle of the body x-axis projected onto a wind-perpendicular plane.
"""
function calc_heading(R_b_to_w, wind_norm)
    e_x = R_b_to_w[:, 1]
    # Project -e_x onto plane perpendicular to wind
    minus_e_x = -e_x
    proj_on_wind = dot(minus_e_x, wind_norm) * wind_norm
    e_x_perp = minus_e_x - proj_on_wind
    # Heading components in wind-perpendicular plane
    wind_cross_z = [wind_norm[2], -wind_norm[1], 0]
    heading_x = dot(e_x_perp, wind_cross_z)
    heading_z = e_x_perp[3]
    heading = atan(heading_x, heading_z)
    return wrap_to_pi(heading)
end

"""
    get_heading_components(e_x, k, θ, wind_norm)

Get heading_y and heading_z components for body x-axis rotated by θ around k.
"""
function get_heading_components(e_x, k, θ, wind_norm)
    e_x_rot = rotate_v_around_k(e_x, k, θ)
    minus_ex = -e_x_rot
    proj_on_wind = dot(minus_ex, wind_norm) * wind_norm
    e_x_perp = minus_ex - proj_on_wind
    wind_cross_z = [wind_norm[2], -wind_norm[1], 0.0]
    hy = dot(e_x_perp, wind_cross_z)
    hz = e_x_perp[3]
    return hy, hz
end

"""
    solve_heading_rotation(R_b_to_w, target_heading, k, wind_norm)

Analytical solution for heading rotation angle.

The heading components vary with rotation angle θ as:
  hy(θ) = A1*sin(θ) + B1*cos(θ) + C1  (same form for hz)

The equation hy*cos(h) - hz*sin(h) = 0 gives: A*sin(θ) + B*cos(θ) + C = 0

Solution: θ = atan2(A, B) - acos(-C / √(A² + B²))
"""
function solve_heading_rotation(R_b_to_w, target_heading, k, wind_norm)
    k = normalize(k)
    e_x = R_b_to_w[:, 1]

    # Extract coefficients by sampling at θ = 0, π/2, π
    hy_0, hz_0 = get_heading_components(e_x, k, 0.0, wind_norm)
    hy_90, hz_90 = get_heading_components(e_x, k, π/2, wind_norm)
    hy_180, hz_180 = get_heading_components(e_x, k, π, wind_norm)

    C1 = (hy_0 + hy_180) / 2
    B1 = hy_0 - C1
    A1 = hy_90 - C1

    C2 = (hz_0 + hz_180) / 2
    B2 = hz_0 - C2
    A2 = hz_90 - C2

    ch = cos(target_heading)
    sh = sin(target_heading)

    A = A1 * ch - A2 * sh
    B = B1 * ch - B2 * sh
    C = C1 * ch - C2 * sh

    r = sqrt(A^2 + B^2)

    if r < 1e-10
        return 0.0
    end

    base_angle = atan(A, B)
    arg = clamp(-C / r, -1.0, 1.0)
    delta = acos(arg)

    return base_angle - delta
end

# ==================== REFINE WING FRAME CALCULATION ==================== #

"""
    get_ref_position_from_points(points::AbstractVector{Point}, ref::Int64)
    get_ref_position_from_points(points::AbstractVector{Point}, refs::Vector{Int64})

Helper to get position (single point or average of multiple).
Used for REFINE wing reference point calculations.
"""
get_ref_position_from_points(points::AbstractVector{Point}, ref::Int64) = points[ref].pos_w
function get_ref_position_from_points(points::AbstractVector{Point}, refs::Vector{Int64})
    n = length(refs)
    return sum(points[idx].pos_w for idx in refs) / n
end

"""
    calc_refine_wing_frame(points::Vector{Point}, z_ref_points, y_ref_points, origin_idx)

Calculate R_b_to_w rotation matrix and origin position for a REFINE wing from structural point positions.

This function implements the same R_b_to_w calculation logic as used in `generate_system.jl`
for REFINE wings, ensuring consistency between initialization (`reinit!`) and simulation
(symbolic equations).

# Algorithm
1. Extract reference point positions (with averaging if vectors provided)
2. Calculate Z-axis (normal to wing): normalized vector from z_p1 to z_p2
3. Calculate X-axis (chord direction): Y_temp × Z, where Y_temp is from y_p1 to y_p2
4. Calculate Y-axis (spanwise): Z × X (ensures orthogonality and right-handed system)
5. Extract origin position from origin_idx point

# Arguments
- `points::AbstractVector{Point}`: All structural points (must have pos_w initialized)
- `z_ref_points::Tuple{Union{Int64, Vector{Int64}}, Union{Int64, Vector{Int64}}}`:
  Reference points defining Z-axis (normal direction)
- `y_ref_points::Tuple{Union{Int64, Vector{Int64}}, Union{Int64, Vector{Int64}}}`:
  Reference points defining Y-axis (spanwise direction)
- `origin_idx::Int64`: Point index defining wing origin (KCU position)

# Returns
- `R_b_to_w::Matrix{SimFloat}`: 3x3 rotation matrix from body frame to world frame
- `origin::KVec3`: Origin position in world frame
"""
function calc_refine_wing_frame(
    points::AbstractVector{Point},
    z_ref_points::Tuple{Union{Int64, Vector{Int64}}, Union{Int64, Vector{Int64}}},
    y_ref_points::Tuple{Union{Int64, Vector{Int64}}, Union{Int64, Vector{Int64}}},
    origin_idx::Int64
)
    # Extract reference point positions (with averaging if vectors provided)
    z_p1, z_p2 = z_ref_points
    y_p1, y_p2 = y_ref_points

    pos_z1 = get_ref_position_from_points(points, z_p1)
    pos_z2 = get_ref_position_from_points(points, z_p2)
    pos_y1 = get_ref_position_from_points(points, y_p1)
    pos_y2 = get_ref_position_from_points(points, y_p2)

    # Build rotation matrix from structural geometry
    # Z direction (normal to wing, normalized)
    z_axis = normalize(pos_z2 - pos_z1)

    # Y temp direction (not necessarily orthogonal yet)
    y_temp = normalize(pos_y2 - pos_y1)

    # X = Y_temp × Z (chord direction, orthogonal to Z)
    x_axis = normalize(y_temp × z_axis)

    # Y = Z × X (ensure orthogonality and right-handed system)
    y_axis = z_axis × x_axis

    # Construct rotation matrix [x y z]
    R_b_to_w = hcat(x_axis, y_axis, z_axis)

    # Extract origin position
    origin = points[origin_idx].pos_w

    return R_b_to_w, origin
end

# ==================== HELPERS ==================== #

"""
    init_untransformed_components!(points, wings, update_vel;
                                   filter=true)

Initialize `pos_w = pos_cad` for untransformed components.
When `filter=true` (default), only components with
`transform_idx == 0` are initialized. When `filter=false`,
all components are initialized (used when no transforms exist).
"""
function init_untransformed_components!(
    points, wings, update_vel::Bool; filter::Bool=true
)
    for point in points
        filter && point.transform_idx != 0 && continue
        point.pos_w .= point.pos_cad
        update_vel && (point.vel_w .= 0.0)
    end
    for wing in wings
        filter && wing.transform_idx != 0 && continue
        # Body frame output
        wing.pos_w .= wing.pos_cad
        wing.Q_b_to_w .= rotation_matrix_to_quaternion(
            wing.R_b_to_c)
        if update_vel
            wing.vel_w .= 0.0
            wing.ω_b .= 0.0
        end
    end
end

"""
    init_principal_frame!(wings, points)

Compute principal frame ODE state from body frame.
Must be called after body frame (pos_w, R_b_to_w,
vel_w, ω_b) is fully initialized.

Sets: com_w, Q_p_to_w, com_vel, ω_p (derived from body
frame), and pos_b for QUATERNION wing points (body
frame, relative to COM).
"""
function init_principal_frame!(wings, points)
    for wing in wings
        R_b_to_w = wing.R_b_to_w
        # COM position in world frame
        wing.com_w .= wing.pos_w .+
            R_b_to_w * wing.com_offset_b
        # Principal frame quaternion:
        # R_p_to_w = R_b_to_w * R_b_to_c' * R_p_to_c
        R_p_to_w = R_b_to_w * wing.R_b_to_c' * wing.R_p_to_c
        wing.Q_p_to_w .= rotation_matrix_to_quaternion(
            R_p_to_w)
        # Derive principal velocities from body frame
        ω_w = R_b_to_w * wing.ω_b
        r_com_w = R_b_to_w * wing.com_offset_b
        wing.com_vel .= wing.vel_w .+
            cross(ω_w, r_com_w)
        wing.ω_p .= R_p_to_w' * ω_w
        # pos_b: offset from COM in body frame
        wing.wing_type != QUATERNION && continue
        com_cad = wing.pos_cad .+
            wing.R_b_to_c * wing.com_offset_b
        for point in points
            if point.type == WING &&
               point.wing_idx == wing.idx
                point.pos_b .= wing.R_b_to_c' *
                    (point.pos_cad - com_cad)
            end
        end
    end
end

# ==================== REINIT! ==================== #

"""
    reinit!(transforms::AbstractVector{Transform}, sys_struct::SystemStructure)

Apply the initial spatial transformations to all components in a `SystemStructure`.

This function iterates through all transforms and applies the specified translation
and rotation to position and orient the kite system components correctly in the
world frame at the beginning of a simulation.

If transforms is empty, simply initializes pos_w = pos_cad for all components.
"""
function reinit!(transforms::AbstractVector{Transform}, sys_struct::SystemStructure;
                 update_vel::Bool=true)
    @unpack points, wings = sys_struct

    # No transforms: init ALL as pos_w = pos_cad
    if isempty(transforms)
        init_untransformed_components!(
            points, wings, update_vel; filter=false)
        init_principal_frame!(wings, points)
        return
    end

    # Initialize untransformed components (transform_idx == 0)
    init_untransformed_components!(points, wings, update_vel)

    # Apply transforms
    for transform in transforms
        # Warn if turn_rate is not zero (not yet implemented)
        if transform.turn_rate != 0.0
            @warn "Transform #$(transform.idx): turn_rate = $(rad2deg(transform.turn_rate))°/s is not zero, " *
                  "but turn_rate dynamics are not yet implemented. This field will be ignored."
        end

        # ==================== TRANSLATE ==================== #
        base_pos, curr_base_pos = get_base_pos(transform, transforms, wings, points)
        T = base_pos - curr_base_pos
        for point in points
            if point.transform_idx == transform.idx
                point.pos_w .= point.pos_cad .+ T
                update_vel && (point.vel_w .= 0.0)
            end
        end
        for wing in wings
            if wing.transform_idx == transform.idx
                wing.pos_w .= wing.pos_cad .+ T
                update_vel && (wing.vel_w .= 0.0)
            end
        end

        # ==================== ROTATE (azimuth/elevation only) ==================== #
        curr_rot_pos = get_rot_pos(transform, wings, points)
        rel_pos = curr_rot_pos - base_pos

        # Error if wing/rot position coincides with base (cannot define rotation)
        if norm(rel_pos) < 1e-6
            error("Transform #$(transform.idx): Wing/rot position and base position " *
                  "overlap at $(base_pos). Cannot define elevation/azimuth rotation. " *
                  "Use transform_idx: 0 to skip transforms, or adjust positions.")
        else
            curr_R_t_to_w = calc_R_t_to_w(rel_pos)
        end

        transform_pos = rotate_around_z(rotate_around_y([1,0,0], -transform.elevation),
                                        -transform.azimuth)
        R_t_to_w = calc_R_t_to_w(transform_pos)

        # Compute velocity components from spherical coordinate motion
        elev = transform.elevation
        azim = transform.azimuth
        rot_pos = get_rot_pos(transform, wings, points)
        r_rot = norm(rot_pos - base_pos)
        vel_spherical = rotate_around_y([0, 0, r_rot * transform.elevation_vel], -elev) +
                        rotate_around_z([0, r_rot * transform.azimuth_vel, 0], -azim)

        # First apply only azimuth/elevation rotation (heading=0)
        for point in points
            if point.transform_idx == transform.idx
                vec = point.pos_w - base_pos
                point.pos_w .= base_pos + apply_heading(vec, R_t_to_w, curr_R_t_to_w, 0.0)
                if update_vel
                    point.vel_w .= norm(point.pos_w - base_pos) / norm(rot_pos - base_pos) *
                                   vel_spherical
                end
            end
        end
        for wing in wings
            if wing.transform_idx == transform.idx
                vec = wing.pos_w - base_pos
                wing.pos_w .= base_pos + apply_heading(vec, R_t_to_w, curr_R_t_to_w, 0.0)
                if update_vel
                    wing.vel_w .= norm(wing.pos_w - base_pos) / norm(rot_pos - base_pos) *
                                  vel_spherical
                    wing.ω_b .= 0.0
                end
            end
        end

        # ==================== APPLY HEADING VIA NONLINEAR SOLVE ==================== #
        # For each wing in this transform, calculate R_b_to_w and solve for heading
        for wing in wings
            if wing.transform_idx == transform.idx
                # Calculate R_b_to_w depending on wing type
                if !isnothing(wing.z_ref_points)
                    # Ref points available (REFINE or
                    # QUATERNION with ref points)
                    R_b_to_w, _ = calc_refine_wing_frame(
                        points, wing.z_ref_points,
                        wing.y_ref_points,
                        wing.origin_idx)
                else
                    # No ref points: rotate R_b_to_c
                    R_b_to_w = zeros(3, 3)
                    for i in 1:3
                        R_b_to_w[:, i] .= apply_heading(
                            wing.R_b_to_c[:, i],
                            R_t_to_w, curr_R_t_to_w, 0.0)
                    end
                end

                # Solve for the rotation angle that
                # achieves target heading
                k = normalize(wing.pos_w)
                wind_norm = normalize(
                    sys_struct.wind_vec_gnd)
                delta_heading = solve_heading_rotation(
                    R_b_to_w, transform.heading,
                    k, wind_norm)

                # Apply the solved rotation to all
                # points in this transform
                for point in points
                    if point.transform_idx ==
                       transform.idx
                        point.pos_w .= rotate_v_around_k(
                            point.pos_w,
                            k, delta_heading)
                    end
                end

                # Apply rotation to wing position and orientation
                wing.pos_w .= rotate_v_around_k(wing.pos_w, k, delta_heading)
                for i in 1:3
                    R_b_to_w[:, i] .= rotate_v_around_k(R_b_to_w[:, i], k, delta_heading)
                end
                wing.R_b_to_w = R_b_to_w
            end
        end
    end

    # Calculate pos_b for REFINE wing points after all
    # transforms are complete
    for wing in wings
        if wing.wing_type == REFINE
            R_b_to_w, origin = calc_refine_wing_frame(
                points, wing.z_ref_points,
                wing.y_ref_points, wing.origin_idx)

            wing.R_b_to_w = R_b_to_w
            wing.pos_w .= origin

            for point in points
                if point.type == WING &&
                   point.wing_idx == wing.idx
                    point.pos_b .= R_b_to_w' *
                        (point.pos_w - origin)
                end
            end
        end
    end

    # Compute principal frame state from body frame
    # for all wings (QUATERNION pos_b also set here)
    init_principal_frame!(wings, points)
end

"""
    reposition!(transforms::AbstractVector{Transform},
                sys_struct::SystemStructure)

Update the system's spatial orientation based on its current
position, preserving velocities.

Unlike `reinit!`, uses current world positions (`pos_w`) as
the starting point rather than resetting from CAD coordinates.
Heading is wind-relative (absolute), consistent with `reinit!`:
the same analytical `solve_heading_rotation` is used.

Preserves existing velocities (`vel_w`, `ω_b`) of all points
and wings.

# Arguments
- `sys_struct::SystemStructure`: The system model to update.
"""
function reposition!(
    transforms::AbstractVector{Transform},
    sys_struct::SystemStructure
)
    @unpack points, wings = sys_struct
    for transform in transforms
        base_pos = points[transform.base_point_idx].pos_w
        rot_pos = get_rot_pos(transform, wings, points)
        curr_rel_pos = rot_pos - base_pos

        if norm(curr_rel_pos) < 1e-6
            error(
                "Transform #$(transform.idx): " *
                "Wing/rot position and base position " *
                "overlap at $(base_pos). " *
                "Cannot define rotation.")
        else
            curr_R_t_to_w = calc_R_t_to_w(curr_rel_pos)
        end

        transform_pos = rotate_around_z(
            rotate_around_y(
                [1,0,0], -transform.elevation),
            -transform.azimuth)
        R_t_to_w = calc_R_t_to_w(transform_pos)

        # ========= PHASE 1: AZIMUTH/ELEVATION ========= #
        for point in points
            if point.transform_idx == transform.idx
                vec = point.pos_w - base_pos
                point.pos_w .= base_pos +
                    apply_heading(vec, R_t_to_w,
                                  curr_R_t_to_w, 0.0)
            end
        end
        for wing in wings
            if wing.transform_idx == transform.idx
                vec = wing.pos_w - base_pos
                wing.pos_w .= base_pos +
                    apply_heading(vec, R_t_to_w,
                                  curr_R_t_to_w, 0.0)
            end
        end

        # ======== PHASE 2: HEADING VIA SOLVE ========== #
        for wing in wings
            if wing.transform_idx == transform.idx
                # R_b_to_w after azim/elev rotation
                if !isnothing(wing.z_ref_points)
                    R_b_to_w, _ = calc_refine_wing_frame(
                        points, wing.z_ref_points,
                        wing.y_ref_points,
                        wing.origin_idx)
                else
                    R_b_to_w = zeros(3, 3)
                    for i in 1:3
                        R_b_to_w[:, i] .= apply_heading(
                            wing.R_b_to_w[:, i],
                            R_t_to_w, curr_R_t_to_w,
                            0.0)
                    end
                end

                # Analytical heading solve
                k = normalize(wing.pos_w)
                wind_norm = normalize(
                    sys_struct.wind_vec_gnd)
                delta_heading = solve_heading_rotation(
                    R_b_to_w, transform.heading,
                    k, wind_norm)

                # Apply to all points in this transform
                for point in points
                    if point.transform_idx ==
                       transform.idx
                        point.pos_w .= rotate_v_around_k(
                            point.pos_w,
                            k, delta_heading)
                    end
                end

                # Apply to wing position and orientation
                wing.pos_w .= rotate_v_around_k(
                    wing.pos_w, k, delta_heading)
                for i in 1:3
                    R_b_to_w[:, i] .= rotate_v_around_k(
                        R_b_to_w[:, i],
                        k, delta_heading)
                end
                wing.R_b_to_w = R_b_to_w
            end
        end
    end

    # REFINE wings: recalculate R_b_to_w and pos_b
    # from structural points
    for wing in wings
        if wing.wing_type == REFINE
            R_b_to_w, origin = calc_refine_wing_frame(
                points, wing.z_ref_points,
                wing.y_ref_points, wing.origin_idx)

            wing.R_b_to_w = R_b_to_w
            wing.pos_w .= origin

            for point in points
                if point.type == WING &&
                   point.wing_idx == wing.idx
                    point.pos_b .= R_b_to_w' *
                        (point.pos_w - origin)
                end
            end
        end
    end

    # Compute principal frame state from body frame
    init_principal_frame!(wings, points)
end
