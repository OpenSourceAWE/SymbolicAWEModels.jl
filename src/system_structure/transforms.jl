# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

"""
Transform functions for heading calculation and spatial positioning.

This file contains:
- Heading calculation functions (calc_heading, apply_heading, etc.)
- reinit! and reposition! functions for applying transforms

Note: The Transform struct and its constructors are defined in types.jl
"""

function finalize_transforms! end

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
    calc_heading(R_b_to_w, wing_pos)

Calculate heading angle using the tangential sphere frame.

Projects the body x-axis onto the tangent plane of the tether
sphere at `wing_pos`. Heading is measured from the elevation
direction (x_t, away from zenith) toward the azimuthal direction
(y_t). Heading = 0 when the kite nose points toward the ground
station.
"""
function calc_heading(R_b_to_w, wing_pos)
    R_t_to_w = calc_R_t_to_w(wing_pos)
    e_x = R_b_to_w[:, 1]
    e_x_t = R_t_to_w' * e_x
    return atan(e_x_t[2], e_x_t[1])
end

"""
    solve_heading_rotation(R_b_to_w, target_heading, wing_pos)

Calculate the rotation angle around the radial axis needed to
achieve `target_heading`.

With the tangential sphere heading, rotating around the radial
axis by θ simply shifts the heading by θ, so the solution is
`target_heading - current_heading`.
"""
function solve_heading_rotation(
    R_b_to_w, target_heading, wing_pos,
)
    current = calc_heading(R_b_to_w, wing_pos)
    return wrap_to_pi(target_heading - current)
end

# ==================== PARTICLE_DYNAMICS WING FRAME CALCULATION ==================== #

"""
    get_ref_position_from_points(points, ref_pt; field=:pos_w)

Weighted position from structural points. `field` selects
which point coordinate to read (`:pos_w` or `:pos_cad`).
"""
function get_ref_position_from_points(
    points::AbstractVector{Point},
    ref_pt::WeightedRefPoints;
    field::Symbol=:pos_w,
)
    pos = zero(KVec3)
    for (idx, w) in zip(ref_pt.ids, ref_pt.weights)
        pos += w * getproperty(points[idx], field)
    end
    return pos
end

"""
    calc_particle_dynamics_wing_frame(points, z_ref_points, y_ref_points, origin)

Calculate R_b_to_w rotation matrix and origin position
from structural point positions.

# Algorithm
1. Weighted ref point positions
2. Z-axis (normal): z_p1 → z_p2
3. X-axis (chord): Y_temp × Z
4. Y-axis (span): Z × X (orthogonal, right-handed)
5. Origin from weighted `origin` ref points
"""
function calc_particle_dynamics_wing_frame(
    points::AbstractVector{Point},
    z_ref_points::Tuple{WeightedRefPoints,
                        WeightedRefPoints},
    y_ref_points::Tuple{WeightedRefPoints,
                        WeightedRefPoints},
    origin::WeightedRefPoints,
)
    z_p1, z_p2 = z_ref_points
    y_p1, y_p2 = y_ref_points

    pos_z1 = get_ref_position_from_points(points, z_p1)
    pos_z2 = get_ref_position_from_points(points, z_p2)
    pos_y1 = get_ref_position_from_points(points, y_p1)
    pos_y2 = get_ref_position_from_points(points, y_p2)

    # Z direction (normal to wing, normalized)
    z_axis = normalize(pos_z2 - pos_z1)

    # Y temp direction (not necessarily orthogonal)
    y_temp = normalize(pos_y2 - pos_y1)

    # X = Y_temp × Z (chord, orthogonal to Z)
    x_axis = normalize(y_temp × z_axis)

    # Y = Z × X (orthogonal, right-handed)
    y_axis = z_axis × x_axis

    R_b_to_w = hcat(x_axis, y_axis, z_axis)
    origin_pos = get_ref_position_from_points(points, origin)

    return R_b_to_w, origin_pos
end

# ==================== HELPERS ==================== #

"""
    copy_cad_to_world!(points, bodies; update_vel=true)

Copy CAD geometry to world frame for ALL points and bodies.
Sets `pos_w = pos_cad` (and `Q_b_to_w` to initial CAD orientation
for bodies). Must be called before `reinit!(transforms, ...)`.
"""
function copy_cad_to_world!(points, bodies; update_vel::Bool=true)
    for point in points
        point.pos_w .= point.pos_cad
        update_vel && (point.vel_w .= 0.0)
    end
    for body in bodies
        body.pos_w .= body.pos_cad
        body.Q_b_to_w .= rotation_matrix_to_quaternion(body.R_b_to_c)
        if update_vel
            body.vel_w .= 0.0
            body.ω_b .= 0.0
        end
    end
end

"""
    min_rotation(curr_dir, target_dir) -> (axis, angle)

Axis and angle of the minimal rotation taking unit vector `curr_dir` onto
`target_dir`. The rotation lies in the plane the two directions span, so it adds
no roll about either — the tangential (heading) orientation is left untouched.
Falls back to an arbitrary perpendicular axis when the directions are anti-parallel.
"""
function min_rotation(curr_dir, target_dir)
    cosang = clamp(dot(curr_dir, target_dir), -1.0, 1.0)
    if cosang > 1 - 1e-12
        return [0.0, 0.0, 1.0], 0.0
    elseif cosang < -1 + 1e-12
        ref = abs(curr_dir[1]) < 0.9 ? [1.0, 0.0, 0.0] : [0.0, 1.0, 0.0]
        return normalize(cross(curr_dir, ref)), Float64(π)
    end
    return normalize(cross(curr_dir, target_dir)), acos(cosang)
end

"""
    spherical_spin(transform) -> Vector

Angular velocity [rad/s], in world axes, of the rigid rotation about the
transform's base that carries `elevation_vel` and `azimuth_vel`. Both axes are
the ones the transform's own position convention turns about: `-y` of the
tangential frame for elevation, world `-z` for azimuth. A component's velocity is
then `spin` crossed into `pos_w - base_pos`, which moves the whole system rigidly
rather than approximating each component's own spherical velocity.
"""
function spherical_spin(transform)
    azim = transform.azimuth
    return transform.elevation_vel .* [-sin(azim), -cos(azim), 0.0] .-
           transform.azimuth_vel .* [0.0, 0.0, 1.0]
end

"""
    apply_azimuth_elevation!(transform, points, bodies, base_pos)

Apply the azimuth/elevation rotation of a single transform to all components in
it (points and bodies). Rotates the current radial onto the target radial by the
minimal (roll-free) rotation, so placement never depends on the source frame —
which is undefined when the components start at the zenith. Roll about the radial
is set afterwards by the heading step, which is well-defined at the target
elevation/azimuth. Returns `(curr_R_t_to_w, R_t_to_w)` for use in that step.
"""
function apply_azimuth_elevation!(transform, points, bodies, base_pos)
    curr_rot_pos = get_rot_pos(transform, bodies, points)
    rel_pos = curr_rot_pos - base_pos

    if norm(rel_pos) < 1e-6
        error("Transform #$(transform.idx): Wing/rot position and base " *
              "position overlap at $(base_pos). Cannot define " *
              "elevation/azimuth rotation. Use transform_idx: 0 to skip " *
              "transforms, or adjust positions.")
    end
    curr_R_t_to_w = calc_R_t_to_w(rel_pos)

    transform_pos = rotate_around_z(
        rotate_around_y([1, 0, 0], -transform.elevation), -transform.azimuth)
    R_t_to_w = calc_R_t_to_w(transform_pos)

    if abs(abs(transform.elevation) - π / 2) < 1e-6
        @warn "Transform #$(transform.idx): elevation = " *
              "$(round(rad2deg(transform.elevation); digits=2))° is at the " *
              "zenith/nadir, where azimuth and heading are undefined."
    end

    axis, angle = min_rotation(normalize(rel_pos), normalize(transform_pos))

    for point in points
        point.transform_idx == transform.idx || continue
        point.pos_w .= base_pos .+ rotate_v_around_k(point.pos_w .- base_pos, axis, angle)
    end
    for body in bodies
        body.transform_idx == transform.idx || continue
        body.pos_w .= base_pos .+ rotate_v_around_k(body.pos_w .- base_pos, axis, angle)
        R_b = quaternion_to_rotation_matrix(body.Q_b_to_w)
        for i in 1:3
            R_b[:, i] .= rotate_v_around_k(R_b[:, i], axis, angle)
        end
        body.Q_b_to_w .= rotation_matrix_to_quaternion(R_b)
    end

    return curr_R_t_to_w, R_t_to_w
end

"""
    apply_spherical_velocity!(transform, points, bodies, base_pos)

Put every component of `transform` on the rigid rotation [`spherical_spin`](@ref)
describes, as `spin x (pos_w - base_pos)`. Runs after the heading step, the
positions it reads being the placed ones; a particle wing's `ω_b` is not part of
its state, so only a rigid body takes the rate.
"""
function apply_spherical_velocity!(transform, points, bodies, base_pos)
    spin = spherical_spin(transform)
    for point in points
        point.transform_idx == transform.idx || continue
        point.vel_w .= cross(spin, point.pos_w .- base_pos)
    end
    for body in bodies
        body.transform_idx == transform.idx || continue
        body.vel_w .= cross(spin, body.pos_w .- base_pos)
        body.ω_b .= body.dynamics_type == RIGID_DYNAMICS ?
            quaternion_to_rotation_matrix(body.Q_b_to_w)' * spin : 0.0
    end
end

"""
    apply_heading!(transform, points, bodies,
                    curr_R_t_to_w, R_t_to_w, base_pos)

Apply heading rotation to all components in a single transform.
Rotates around the radial axis through `base_pos` (not the origin).
Uses the reference body's `R_b_to_w` for the no-ref-points orientation source.
After `copy_cad_to_world!`, this equals `R_b_to_c` (for
`reinit!`), or the current world orientation (for `reposition!`).
Bodies in the transform rotate with the same heading delta; a transform
without a body target applies no heading (matching point behavior).
"""
function apply_heading!(transform, points, bodies,
                         curr_R_t_to_w, R_t_to_w, base_pos)
    reference_body = heading_reference_body(transform, bodies)
    isnothing(reference_body) && return

    if !isnothing(reference_body.z_ref_points)
        R_b_to_w, _ = calc_particle_dynamics_wing_frame(
            points, reference_body.z_ref_points,
            reference_body.y_ref_points, reference_body.origin)
    else
        R_b_to_w = quaternion_to_rotation_matrix(reference_body.Q_b_to_w)
    end

    rel_pos = reference_body.pos_w - base_pos
    delta_heading = solve_heading_rotation(
        R_b_to_w, transform.heading, rel_pos)
    k = normalize(rel_pos)

    for point in points
        point.transform_idx == transform.idx || continue
        point.pos_w .= base_pos .+ rotate_v_around_k(
            point.pos_w .- base_pos, k, delta_heading)
    end
    for body in bodies
        body.transform_idx == transform.idx || continue
        body.pos_w .= base_pos .+ rotate_v_around_k(
            body.pos_w .- base_pos, k, delta_heading)
        R_b = quaternion_to_rotation_matrix(body.Q_b_to_w)
        for i in 1:3
            R_b[:, i] .= rotate_v_around_k(R_b[:, i], k, delta_heading)
        end
        body.Q_b_to_w .= rotation_matrix_to_quaternion(R_b)
    end
end

"""
    heading_reference_body(transform, bodies) -> Body or nothing

The single body whose frame defines the transform's heading: the wing/rot body
(`transform.wing_idx`) when set, otherwise the first body in the transform. The
heading rotation is applied once about the radial through it — not once per body,
which would compound into a spurious rotation for multi-body (beam) transforms.
"""
function heading_reference_body(transform, bodies)
    if !isnothing(transform.wing_idx)
        return bodies[something(transform.wing_idx)]
    end
    for body in bodies
        body.transform_idx == transform.idx && return body
    end
    return nothing
end

"""
    finalize_transforms!(points, bodies)

Finalize transforms: update PARTICLE_DYNAMICS body frames from structural
point positions, then compute principal frame ODE state for every body
(RIGID_DYNAMICS bodies re-derived from the transformed `pos_w`/`Q_b_to_w`).
"""
function finalize_transforms!(points, bodies)
    for body in bodies
        body.dynamics_type == PARTICLE_DYNAMICS || continue
        (isnothing(body.z_ref_points) || isnothing(body.y_ref_points) ||
         isnothing(body.origin)) && continue
        R_b_to_w, origin = calc_particle_dynamics_wing_frame(
            points, body.z_ref_points, body.y_ref_points, body.origin)
        body.R_b_to_w = R_b_to_w
        body.pos_w .= origin
        for point in points
            if point.is_wing_node && point.wing_idx == body.idx
                point.pos_undeformed_b .= R_b_to_w' * (point.pos_w - origin)
                point.pos_b .= point.pos_undeformed_b
            end
        end
    end
    init_principal_frame!(bodies, points)
end

"""
    refresh_deformed_positions!(points, bodies) -> nothing

Update every PARTICLE_DYNAMICS wing node's `pos_b` from its live world
position, measured in its wing's current frame about the wing origin — the same
reference `finalize_transforms!` builds `pos_undeformed_b` against, so the two
differ only by the deformation the points have picked up.

Called once per state sync, so both backends report it from one definition. A
RIGID_DYNAMICS wing's nodes keep the offset they were constructed with, that
being what places them.
"""
function refresh_deformed_positions!(points, bodies)
    for body in bodies
        body.dynamics_type == PARTICLE_DYNAMICS || continue
        isnothing(body.origin) && continue
        frame = body.R_b_to_w'
        for point in points
            if point.is_wing_node && point.wing_idx == body.idx
                point.pos_b .= frame * (point.pos_w - body.pos_w)
            end
        end
    end
    return nothing
end

"""
    init_principal_frame!(bodies, points)

Compute principal frame ODE state from body frame.
Must be called after body frame (`pos_w`, `R_b_to_w`,
`vel_w`, `ω_b`) is fully initialized.

Sets: `com_w`, `Q_p_to_w`, `com_vel`, `ω_p` (derived from body
frame), and `pos_b` for RIGID_DYNAMICS wing points (body
frame, relative to COM).
"""
function init_principal_frame!(bodies, points)
    for body in bodies
        if body.dynamics_type == RIGID_DYNAMICS
            init_principal_state!(body)
            # pos_b: wing nodes' offset from COM in body frame
            com_cad = body.pos_cad .+ body.R_b_to_c * body.com_offset_b
            for point in points
                if point.is_wing_node && point.wing_idx == body.idx
                    point.pos_undeformed_b .=
                        body.R_b_to_c' * (point.pos_cad - com_cad)
                    point.pos_b .= point.pos_undeformed_b
                end
            end
        else
            # PARTICLE: R_b_to_p is identity, so derive R_p_to_w from R_b_to_c/R_p_to_c.
            R_b_to_w = body.R_b_to_w::Matrix{SimFloat}
            body.com_w .= body.pos_w .+ R_b_to_w * body.com_offset_b
            R_p_to_w = R_b_to_w * body.R_b_to_c' * body.R_p_to_c
            body.Q_p_to_w .= rotation_matrix_to_quaternion(R_p_to_w)
            ω_w = R_b_to_w * body.ω_b
            body.com_vel .= body.vel_w .+ cross(ω_w, R_b_to_w * body.com_offset_b)
            body.ω_p .= R_p_to_w' * ω_w
        end
    end
end

# ==================== REINIT! ==================== #

"""
    reinit!(transforms::AbstractVector{Transform}, sys_struct::SystemStructure;
            update_vel=true)

Apply transforms to all components in a `SystemStructure`.

Expects `pos_w` to already be set (via `copy_cad_to_world!` and optionally
`apply_tether_init_stretched_lens!` from `reinit!(sys_struct, set; ...)`).
Applies: translate (from pos_w) → azimuth/elevation → heading.
"""
function reinit!(transforms::AbstractVector{Transform}, sys_struct::SystemStructure;
                 update_vel::Bool=true)
    (; points, bodies) = sys_struct

    if isempty(transforms)
        finalize_transforms!(points, bodies)
        return
    end

    for transform in transforms
        if transform.turn_rate != 0.0
            @warn "Transform #$(transform.idx): turn_rate = " *
                  "$(rad2deg(transform.turn_rate))°/s is not zero, " *
                  "but turn_rate dynamics are not yet implemented. " *
                  "This field will be ignored."
        end

        # ==================== TRANSLATE ==================== #
        base_pos, curr_base_pos = get_base_pos(transform, transforms, bodies, points)
        T = base_pos - curr_base_pos
        for point in points
            point.transform_idx == transform.idx || continue
            point.pos_w .= point.pos_w .+ T
            update_vel && (point.vel_w .= 0.0)
        end
        # Wings are bodies, handled by the bodies loop.
        for body in bodies
            body.transform_idx == transform.idx || continue
            body.pos_w .= body.pos_w .+ T
            update_vel && (body.vel_w .= 0.0)
        end

        # ==================== ROTATE + HEADING ==================== #
        curr_R_t_to_w, R_t_to_w = apply_azimuth_elevation!(
            transform, points, bodies, base_pos)
        apply_heading!(transform, points, bodies,
            curr_R_t_to_w, R_t_to_w, base_pos)
        update_vel && apply_spherical_velocity!(
            transform, points, bodies, base_pos)
    end

    finalize_transforms!(points, bodies)
end

"""
    reposition!(transforms::AbstractVector{Transform},
                sys_struct::SystemStructure; update_vel=false)

Update the system's spatial orientation based on its current
position, preserving velocities. `update_vel` instead overwrites
them with the rigid rotation each transform's `elevation_vel` and
`azimuth_vel` describe, which is how a state is placed already
moving on its arc rather than at rest.

Unlike `reinit!`, uses current world positions (`pos_w`) as
the starting point (no reset from CAD coordinates, no tether
length scaling). Heading uses the tangential sphere frame,
consistent with `reinit!`.
"""
function reposition!(
    transforms::AbstractVector{Transform},
    sys_struct::SystemStructure;
    update_vel::Bool=false
)
    (; points, bodies) = sys_struct
    for transform in transforms
        base_pos = if !isnothing(
                transform.base_transform_idx)
            base_tf = transforms[something(
                transform.base_transform_idx)]
            get_rot_pos(base_tf, bodies, points)
        else
            points[something(
                transform.base_point_idx)].pos_w
        end
        curr_R_t_to_w, R_t_to_w =
            apply_azimuth_elevation!(
                transform, points, bodies, base_pos)
        apply_heading!(transform, points, bodies,
            curr_R_t_to_w, R_t_to_w, base_pos)
        update_vel && apply_spherical_velocity!(
            transform, points, bodies, base_pos)
    end
    finalize_transforms!(points, bodies)
end
