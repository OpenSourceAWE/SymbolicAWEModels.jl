# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

module SymbolicAWEModelsMakieExt

using Makie
using UnPack
using LinearAlgebra
using StaticArrays
using Statistics
using Printf
using Dates
using LaTeXStrings
using KiteUtils
using KiteUtils: SysLog
using SymbolicAWEModels
using VortexStepMethod
using GeometryBasics

# Global storage for plot observables (for time-based plotting)
const PLOT_GEOMETRY_OBS = Ref{Union{Nothing, Observable}}(nothing)  # Single trigger observable
const PLOT_SEGMENT_COLORS_OBS = Ref{Union{Nothing, Observable}}(nothing)  # Separate for force coloring
const PLOT_SCENE = Ref{Union{Nothing, Scene}}(nothing)
const PLOT_RELEVANT_PLOTS = Ref{Union{Nothing, Vector}}(nothing)
const PLOT_SYSTEM_STRUCTURE = Ref{Union{Nothing, SystemStructure}}(nothing)
const PLOT_VECTOR_SCALE = Ref{Float64}(1.0)
const PLOT_FORCE_COLOR = Ref{Bool}(false)
const PLOT_SEGMENT_COLOR = Ref{Symbol}(:black)
const PLOT_ZOOMED_IN = Ref{Bool}(false)
const PLOT_ZOOM_RELMARGIN = Ref{Float64}(0.2)
const PLOT_ZOOM_SEGMENT_IDX = Ref{Int}(-1)  # Which segment we're zoomed into (-1 = none)
const PLOT_ZOOM_BODY_IDX = Ref{Int}(-1)     # Which rigid body we're zoomed into (-1 = none)
const PLOT_BODY_FRAME = Ref{Bool}(false)  # Whether body frame tracking is active
const PLOT_CAMERA_DISTANCE = Ref{Union{Nothing, Float64}}(nothing)  # Stored camera distance
const PLOT_PREV_BODY_FRAME = Ref{Bool}(false)  # Previous body frame state
const PLOT_PREV_ZOOMED_IN = Ref{Bool}(false)  # Previous zoomed state
const PLOT_PREV_SEGMENT_IDX = Ref{Int}(-1)  # Previous segment index
const PLOT_WORLD_EYEPOS = Ref{Union{Nothing, Vec3f}}(nothing)
const PLOT_WORLD_LOOKAT = Ref{Union{Nothing, Vec3f}}(nothing)
const PLOT_BODY_DISTANCE = Ref{Union{Nothing, Float64}}(nothing)
const PLOT_BODY_PREV_WING_POS = Ref{Union{Nothing, Vec3f}}(nothing)

"""
    wing_log_pos(sl, sys, wing, k)

World position of `wing` at log frame `k`, read from its dedicated
slot appended after the structural points and panel corners in the
`SysLog`. Older logs lack that slot; in that case fall back to the
mean of the wing's `WING`-type structural points (or all points).
"""
function wing_log_pos(sl, sys, wing, k)
    i = SymbolicAWEModels.position_slots(sys).wings[wing.idx]
    if i <= length(sl.X[k])
        return SVector{3, Float64}(
            sl.X[k][i], sl.Y[k][i], sl.Z[k][i])
    end
    idxs = [p.idx for p in sys.points
            if p.type == WING && p.wing_idx == wing.idx]
    isempty(idxs) && (idxs = eachindex(sl.X[k]))
    return SVector{3, Float64}(
        mean(@view sl.X[k][idxs]),
        mean(@view sl.Y[k][idxs]),
        mean(@view sl.Z[k][idxs]))
end

# Multi-system plotting support
const PLOT_MULTI_SYSTEMS = Ref{Union{Nothing, Vector{<:SystemStructure}}}(nothing)
const PLOT_MULTI_GEOMETRY_OBS = Ref{Union{Nothing, Vector{Observable}}}(nothing)

"""
    calculate_segment_force_colors(segments, segment_color)

Calculate segment colors based on their force values.
Maps forces from green (low) to red (high) using linear interpolation.

# Arguments
- `segments`: Collection of segments with force field
- `segment_color`: Default color to use when all forces are equal

# Returns
- Vector of RGBAf colors, one per segment
"""
function calculate_segment_force_colors(segments, segment_color)
    forces = [seg.force for seg in segments]
    max_force = maximum(forces)
    min_force = minimum(forces)
    force_range = max_force - min_force

    return [begin
        if force_range > 0
            normalized_force = (seg.force - min_force) / force_range
            # Interpolate from green (0,1,0) to red (1,0,0)
            RGBAf(normalized_force, 1.0 - normalized_force, 0.0, 1.0)
        else
            to_color(segment_color)
        end
    end for seg in segments]
end

function SymbolicAWEModels.plot_wing_aero!(ax, sys, wing,
        mode::SymbolicAWEModels.AbstractAeroModel;
        use_observables=false, geometry_obs=nothing)
    return nothing
end

function SymbolicAWEModels.plot_wing_aero!(ax, sys, wing,
        mode::SymbolicAWEModels.AbstractVSMAero;
        use_observables=false, geometry_obs=nothing)
    return plot!(ax, mode.vsm_aero; R_b_w=wing.R_b_to_w,
                 T_b_w=aero_plot_translation(wing), use_observables)
end

function SymbolicAWEModels.plot_wing_aero!(ax, sys, wing, mode::AeroPlate;
        use_observables=false, geometry_obs=nothing)
    quad_vertices() = [Point3f(corner)
        for twist_surface_idx in wing.twist_surface_idxs
        for corner in SymbolicAWEModels.plate_corners(
            sys.twist_surfaces[twist_surface_idx],
            sys.points[
                sys.twist_surfaces[twist_surface_idx].point_idxs[1]].pos_w,
            wing.R_b_to_w)]
    initial = quad_vertices()
    isempty(initial) && return nothing

    nan_point = Point3f(NaN, NaN, NaN)
    quad_borders(verts) = [point
        for quad_num in 1:(length(verts) ÷ 4)
        for point in (verts[4quad_num-3], verts[4quad_num-2],
                      verts[4quad_num-1], verts[4quad_num],
                      verts[4quad_num-3], nan_point)]
    faces = GLTriangleFace[]
    for quad_num in 1:(length(initial) ÷ 4)
        base = 4 * (quad_num - 1)
        push!(faces, GLTriangleFace(base + 1, base + 2, base + 3))
        push!(faces, GLTriangleFace(base + 1, base + 3, base + 4))
    end

    animate = use_observables && !isnothing(geometry_obs)
    vertices = animate ?
        @lift(begin; $geometry_obs; quad_vertices(); end) : initial
    p = mesh!(ax, vertices, faces; color=(:red, 0.2), transparency=true)
    borders = animate ?
        @lift(quad_borders($vertices)) : quad_borders(initial)
    lines!(ax, borders; color=:black, transparency=true)
    return p
end

SymbolicAWEModels.update_wing_aero_plot!(wing,
    ::SymbolicAWEModels.AbstractAeroModel) = nothing

SymbolicAWEModels.update_wing_aero_plot!(wing,
    mode::SymbolicAWEModels.AbstractVSMAero) =
    plot!(mode.vsm_aero; R_b_w=wing.R_b_to_w,
          T_b_w=aero_plot_translation(wing))

"""
    aero_plot_translation(wing)

World-frame translation for the wing's VSM panels with `aero_z_offset`
removed, so panels render at the structural pose instead of the shifted
aero pose.
"""
aero_plot_translation(wing) =
    wing.pos_w .- wing.R_b_to_w * [0.0, 0.0, wing.aero_z_offset]

"""
    body_frame_arrows(frames, scale)

Build flat `(origins, directions)` arrays of RGB body-frame axes (x→red,
y→green, z→blue) for a list of `(position, rotation_matrix)` frames, each axis
scaled by `scale`. Shared by the wing and rigid-body renderers.
"""
function body_frame_arrows(frames, scale)
    origins = Point3f[]
    directions = Vec3f[]
    for (position, R) in frames
        origin_pos = Point3f(position)
        for i in 1:3
            push!(origins, origin_pos)
            push!(directions, Vec3f(R[:, i]) * scale)
        end
    end
    return origins, directions
end

"""World-frame `(position, R_b_to_w)` per wing, orientation from the quaternion."""
wing_frames(sys) = [(wing.pos_w,
    SymbolicAWEModels.quaternion_to_rotation_matrix(wing.Q_b_to_w))
    for wing in sys.wings]

"""World-frame `(position, R_b_to_w)` per rigid body, orientation from the quaternion."""
rigid_body_frames(sys) = [(body.pos_w,
    SymbolicAWEModels.quaternion_to_rotation_matrix(body.Q_b_to_w))
    for body in sys.bodies]

"""
    body_joint_spokes(sys) -> Vector{Point3f}

Flat list of segment endpoints (consecutive pairs, for `linesegments!`): for each
elastic joint, a rigid spoke from each connected body's origin to its anchor on
that joint. This shows a body's rigid extent and the actual joint connectivity,
rather than assuming a rod shape or a fixed body order. (Point→body spokes would
slot in here the same way once points can attach to rigid bodies.)
"""
function body_joint_spokes(sys)
    points = Point3f[]
    for joint in sys.elastic_joints
        body_a = sys.bodies[joint.body_a_idx]
        body_b = sys.bodies[joint.body_b_idx]
        R_a = SymbolicAWEModels.quaternion_to_rotation_matrix(body_a.Q_b_to_w)
        R_b = SymbolicAWEModels.quaternion_to_rotation_matrix(body_b.Q_b_to_w)
        push!(points, Point3f(body_a.pos_w),
              Point3f(body_a.pos_w .+ R_a * joint.anchor_a_b))
        push!(points, Point3f(body_b.pos_w),
              Point3f(body_b.pos_w .+ R_b * joint.anchor_b_b))
    end
    return points
end

"""
    spoke_body_indices(sys) -> Vector{Int}

Body index owning each spoke, in the same order as `body_joint_spokes` emits them
(per joint: body_a's spoke then body_b's). Used to highlight a hovered body's
spokes.
"""
function spoke_body_indices(sys)
    idxs = Int[]
    for joint in sys.elastic_joints
        push!(idxs, joint.body_a_idx, joint.body_b_idx)
    end
    return idxs
end

"""
    plot_frame_axes!(ax, plots, key, sys, n, frames_of, scale, geometry_obs)

Plot RGB body-frame axes for `n` entities, reading frames via `frames_of(sys)`.
Static when `geometry_obs === nothing`, otherwise observable-driven from
`PLOT_SYSTEM_STRUCTURE[]`. No-op when `n == 0`.
"""
function plot_frame_axes!(ax, plots, key, sys, n, frames_of, scale, geometry_obs;
                          facets=8)
    n == 0 && return nothing
    axis_colors = repeat([:red, :green, :blue], n)
    if isnothing(geometry_obs)
        origins, directions = body_frame_arrows(frames_of(sys), scale)
        plots[key] = arrows3d!(ax, origins, directions; color=axis_colors, quality=facets)
    else
        origins_dirs = @lift begin
            $geometry_obs  # Trigger dependency
            body_frame_arrows(frames_of(PLOT_SYSTEM_STRUCTURE[]), PLOT_VECTOR_SCALE[])
        end
        plots[key] = arrows3d!(ax, @lift($origins_dirs[1]),
            @lift($origins_dirs[2]); color=axis_colors, quality=facets)
    end
    return nothing
end

function Makie.plot!(ax, sys::SystemStructure;
                     point_color = :darkred, segment_color = :black,
                     wing_colors = Makie.wong_colors(), vector_scale = 1.0,
                     show_points = true, show_segments = true, show_orient = true,
                     force_color = false,
                     plot_vsm = true, plot_aero = true,
                     extra_points = nothing,
                     extra_groups = nothing,
                     mesh = nothing,
                     # `facets` sets the orientation-arrow quality; `linewidth`
                     # sizes the segment/rigid-body lines and `point_size` the
                     # point-node markers.
                     facets::Int = 8,
                     linewidth = 2,
                     point_size = 9,
                     # Optional observable for real-time updates
                     geometry_obs = nothing)

    plots = Dict{Symbol, Any}()

    # === Plot Segments ===
    if show_segments && !isempty(sys.segments)
        if isnothing(geometry_obs)
            # Static plotting: build segment points once
            lineseg_points = Point3f[]
            for seg in sys.segments
                p1 = sys.points[seg.point_idxs[1]].pos_w
                p2 = sys.points[seg.point_idxs[2]].pos_w
                push!(lineseg_points, Point3f(p1))
                push!(lineseg_points, Point3f(p2))
            end
        else
            # Dynamic plotting: compute from PLOT_SYSTEM_STRUCTURE when triggered
            lineseg_points = @lift begin
                $geometry_obs  # Trigger dependency
                sys_ref = PLOT_SYSTEM_STRUCTURE[]
                points = Point3f[]
                for seg in sys_ref.segments
                    p1 = sys_ref.points[seg.point_idxs[1]].pos_w
                    p2 = sys_ref.points[seg.point_idxs[2]].pos_w
                    push!(points, Point3f(p1))
                    push!(points, Point3f(p2))
                end
                points
            end
        end

        num_segments = length(sys.segments)

        # Calculate segment colors based on force if requested
        if force_color
            seg_colors = Observable(calculate_segment_force_colors(sys.segments, segment_color))
        else
            seg_colors = Observable(fill(to_color(segment_color), num_segments))
        end

        plots[:segments] = linesegments!(ax, lineseg_points; color=seg_colors,
                                         linewidth, transparency=true,
                                         label="Segments")
        plots[:segment_colors_obs] = seg_colors
    end

    # === Plot Points ===
    if show_points
        if isnothing(geometry_obs)
            # Static plotting
            point_positions = [Point3f(p.pos_w) for p in sys.points]
        else
            # Dynamic plotting: compute from PLOT_SYSTEM_STRUCTURE when triggered
            point_positions = @lift begin
                $geometry_obs  # Trigger dependency
                [Point3f(p.pos_w) for p in PLOT_SYSTEM_STRUCTURE[].points]
            end
        end
        plots[:points] = scatter!(ax, point_positions, color=point_color,
                                  markersize=point_size, label="Points",
                                  transparency=true)
    end

    # === Plot Rigid Bodies (standalone, e.g. a beam) ===
    if !isempty(sys.bodies)
        # Grey spokes from each body origin to its joint anchors: show the rigid
        # extent and the true joint connectivity (works for chains, stars, Ys).
        if !isempty(sys.elastic_joints)
            if isnothing(geometry_obs)
                spoke_points = body_joint_spokes(sys)
            else
                spoke_points = @lift begin
                    $geometry_obs  # Trigger dependency
                    body_joint_spokes(PLOT_SYSTEM_STRUCTURE[])
                end
            end
            spoke_colors = Observable(fill(to_color(RGBf(0.32, 0.32, 0.32)),
                                           length(spoke_body_indices(sys))))
            plots[:body_chain] = linesegments!(ax, spoke_points; color=spoke_colors,
                                               linewidth, label="Rigid bodies")
            plots[:body_chain_colors_obs] = spoke_colors
        end

        # RGB body-frame axes per body, from the orientation quaternion.
        if show_orient
            plot_frame_axes!(ax, plots, :bodies, sys, length(sys.bodies),
                             rigid_body_frames, vector_scale, geometry_obs; facets)
        end
    end

    # === Plot Wings (RGB body-frame axes from the orientation quaternion) ===
    if show_orient
        plot_frame_axes!(ax, plots, :wings, sys, length(sys.wings),
                         wing_frames, vector_scale, geometry_obs; facets)
    end

    # === Plot VSM Aerodynamics ===
    # VSM panels - use observables if we're in dynamic mode
    if plot_vsm && !isempty(sys.wings)
        plots[:vsm] = []
        use_obs = !isnothing(geometry_obs)
        for wing in sys.wings
            p = SymbolicAWEModels.plot_wing_aero!(ax, sys, wing, wing.aero;
                use_observables=use_obs, geometry_obs)
            isnothing(p) || push!(plots[:vsm], p)
        end
    end

    # === Plot Aero Forces ===
    if plot_aero
        if isnothing(geometry_obs)
            # Static plotting: build aero force arrows
            aero_origins = Point3f[]
            aero_forces_raw = Vec3f[]

            for wing in sys.wings
                if wing.dynamics_type == RIGID_DYNAMICS
                    # For RIGID_DYNAMICS wings, use wing.aero_force_b
                    if !iszero(wing.aero_force_b)
                        aero_force_w = wing.R_b_to_w * wing.aero_force_b
                        push!(aero_origins, Point3f(wing.pos_w))
                        push!(aero_forces_raw, Vec3f(aero_force_w))
                    end
                elseif wing.dynamics_type == SymbolicAWEModels.PARTICLE_DYNAMICS
                    # For PARTICLE_DYNAMICS wings, plot both point forces and total wing force
                    # Plot individual point forces
                    for point in sys.points
                        if point.type == WING && point.wing_idx == wing.idx
                            if !iszero(point.aero_force_b)
                                aero_force_w = wing.R_b_to_w * point.aero_force_b
                                push!(aero_origins, Point3f(point.pos_w))
                                push!(aero_forces_raw, Vec3f(aero_force_w))
                            end
                        end
                    end
                    # Also plot total wing aero force at wing center when vector_scale > 0
                    if vector_scale > 0 && !iszero(wing.aero_force_b)
                        aero_force_w = wing.R_b_to_w * wing.aero_force_b
                        push!(aero_origins, Point3f(wing.pos_w))
                        push!(aero_forces_raw, Vec3f(aero_force_w))
                    end
                end
            end

            if !isempty(aero_origins)
                # Calculate adaptive force scale
                max_force = maximum(norm.(aero_forces_raw))
                # Scale forces to be similar size as vector_scale
                force_scale = vector_scale / max_force
                aero_directions = [f * force_scale for f in aero_forces_raw]

                plots[:aero_forces] = arrows3d!(ax, aero_origins, aero_directions,
                                               color=:magenta,
                                               label="Aero Forces")
            end
        else
            # Dynamic plotting: compute from PLOT_SYSTEM_STRUCTURE when triggered
            aero_origins_dirs = @lift begin
                $geometry_obs  # Trigger dependency
                sys_ref = PLOT_SYSTEM_STRUCTURE[]
                scale = PLOT_VECTOR_SCALE[]
                origins = Point3f[]
                forces_raw = Vec3f[]

                for wing in sys_ref.wings
                    if wing.dynamics_type == RIGID_DYNAMICS
                        if !iszero(wing.aero_force_b)
                            aero_force_w = wing.R_b_to_w * wing.aero_force_b
                            push!(origins, Point3f(wing.pos_w))
                            push!(forces_raw, Vec3f(aero_force_w))
                        end
                    elseif wing.dynamics_type == SymbolicAWEModels.PARTICLE_DYNAMICS
                        for point in sys_ref.points
                            if point.type == WING && point.wing_idx == wing.idx
                                if !iszero(point.aero_force_b)
                                    aero_force_w = wing.R_b_to_w * point.aero_force_b
                                    push!(origins, Point3f(point.pos_w))
                                    push!(forces_raw, Vec3f(aero_force_w))
                                end
                            end
                        end
                        if scale > 0 && !iszero(wing.aero_force_b)
                            aero_force_w = wing.R_b_to_w * wing.aero_force_b
                            push!(origins, Point3f(wing.pos_w))
                            push!(forces_raw, Vec3f(aero_force_w))
                        end
                    end
                end

                # Calculate adaptive force scale
                directions = Vec3f[]
                if !isempty(forces_raw)
                    max_force = maximum(norm.(forces_raw))
                    force_scale = scale / max_force
                    directions = [f * force_scale for f in forces_raw]
                end
                (origins, directions)
            end

            plots[:aero_forces] = arrows3d!(ax, @lift($aero_origins_dirs[1]), @lift($aero_origins_dirs[2]),
                                           color=:magenta)
        end
    end

    # === Calculate arrow scale from data bounding box ===
    data_char_length = 20.0  # fallback for empty systems
    if !isempty(sys.points)
        all_x = [p.pos_w[1] for p in sys.points]
        all_y = [p.pos_w[2] for p in sys.points]
        all_z = [p.pos_w[3] for p in sys.points]

        xr, yr, zr = extrema(all_x), extrema(all_y), extrema(all_z)
        data_char_length = max(
            xr[2] - xr[1], yr[2] - yr[1], zr[2] - zr[1],
            1.0,  # minimum to avoid zero
        )
    elseif !isempty(sys.bodies)
        all_x = [b.pos_w[1] for b in sys.bodies]
        all_y = [b.pos_w[2] for b in sys.bodies]
        all_z = [b.pos_w[3] for b in sys.bodies]
        xr, yr, zr = extrema(all_x), extrema(all_y), extrema(all_z)
        data_char_length = max(xr[2]-xr[1], yr[2]-yr[1], zr[2]-zr[1], 1.0)
    end

    axis_length = data_char_length * 0.2

    # === Plot Global Axes ===
    begin
        origins = [Point3f(0, 0, 0), Point3f(0, 0, 0), Point3f(0, 0, 0)]
        directions = [
            Vec3f(axis_length, 0, 0),
            Vec3f(0, axis_length, 0),
            Vec3f(0, 0, axis_length),
        ]
        axis_colors = [:red, :green, :blue]
        plots[:global_axes] = arrows3d!(ax, origins, directions;
                                        shaftradius=0.02,
                                        tipradius=0.08,
                                        tiplength=0.2,
                                        color=axis_colors,
                                        label="Global Axes")
    end

    # === Plot Mesh (e.g., OBJ file) ===
    if !isnothing(mesh)
        # Transform mesh vertices to world frame using wing orientation
        if !isempty(sys.wings)
            wing = sys.wings[1]
            R_b_w = wing.R_b_to_w
            T_b_w = wing.pos_w
            old_vertices = GeometryBasics.coordinates(mesh)
            new_vertices = [Point3f(R_b_w * Vector(v) + T_b_w) for v in old_vertices]
            transformed_mesh = GeometryBasics.Mesh(new_vertices, GeometryBasics.faces(mesh))
        else
            transformed_mesh = mesh
        end
        plots[:mesh] = mesh!(ax, transformed_mesh;
                            color=(:lightblue, 0.3),
                            transparency=true)
    end

    # === Plot Extra Points (e.g., from external CSV) ===
    if !isnothing(extra_points)
        extra_positions = [Point3f(p...) for p in extra_points]

        # Draw lines connecting points within each group
        if !isnothing(extra_groups)
            for (_gname, indices) in extra_groups
                for i in 1:(length(indices)-1)
                    p1 = extra_positions[indices[i]]
                    p2 = extra_positions[indices[i+1]]
                    lines!(ax, [p1[1], p2[1]], [p1[2], p2[2]], [p1[3], p2[3]];
                           color=(:blue, 0.6), linewidth=2, transparency=true)
                end
            end
        end

        plots[:extra_points] = scatter!(ax, extra_positions, color=:blue,
                                        markersize=10, label="Extra Points",
                                        transparency=true)
    end

    return plots
end

"""
    SymbolicAWEModels.update_plot_observables!(sys::SystemStructure)

Trigger plot updates by updating the geometry observable.

The SystemStructure should already be updated via `update_from_sysstate!`.
This function simply triggers the observable which causes Makie to recompute
all geometry from `PLOT_SYSTEM_STRUCTURE[]` via `@lift` expressions.

# Example
```julia
# Create initial plot
scene = plot(sys_struct)

# In simulation loop:
for step in 1:steps
    next_step!(sam; ...)
    update_plot_observables!(sam.sys_struct)
    sleep(0.001)  # Allow Makie to process updates
end
```
"""
function SymbolicAWEModels.update_plot_observables!(sys::SystemStructure)
    # Trigger geometry observable - this causes all @lift expressions to recompute
    if !isnothing(PLOT_GEOMETRY_OBS[])
        PLOT_GEOMETRY_OBS[][] = time()  # Use timestamp as trigger value
    end

    # Update segment colors if force coloring is enabled
    if !isnothing(PLOT_SEGMENT_COLORS_OBS[]) && PLOT_FORCE_COLOR[]
        PLOT_SEGMENT_COLORS_OBS[][] = calculate_segment_force_colors(sys.segments, PLOT_SEGMENT_COLOR[])
    end

    # Update per-mode aero plots (VSM panel meshes; quad plots update
    # through the geometry observable on their own)
    if !isnothing(PLOT_GEOMETRY_OBS[])
        for wing in sys.wings
            SymbolicAWEModels.update_wing_aero_plot!(wing, wing.aero)
        end
    end

    if !isnothing(PLOT_SCENE[]) && !isnothing(PLOT_RELEVANT_PLOTS[]) && !isnothing(PLOT_SYSTEM_STRUCTURE[])
        scene = PLOT_SCENE[]
        relevant_plots = PLOT_RELEVANT_PLOTS[]
        stored_sys = PLOT_SYSTEM_STRUCTURE[]

        mode_changed = (PLOT_PREV_BODY_FRAME[] != PLOT_BODY_FRAME[] ||
                        PLOT_PREV_ZOOMED_IN[] != PLOT_ZOOMED_IN[] ||
                        PLOT_PREV_SEGMENT_IDX[] != PLOT_ZOOM_SEGMENT_IDX[])

        apply_zoom_mode!(scene, relevant_plots, stored_sys; mode_changed)
    end

    return nothing
end

"""
    point_line_segment_distance(p, a, b)

Calculate the minimum 2D distance from a point `p` to a line segment defined
by points `a` and `b`.
"""
function point_line_segment_distance(p, a, b)
    # Vector from a to b
    ab = b - a
    # Vector from a to p
    ap = p - a
    # Length squared of the segment.
    len_sq = dot(ab, ab)
    # If the segment is just a point, return distance from p to a.
    if len_sq ≈ 0.0
        return norm(ap)
    end
    # Project ap onto ab to find the closest point on the infinite line.
    # t is the normalized position of the projection.
    t = dot(ap, ab) / len_sq
    # Clamp t to the [0, 1] range to stay on the segment.
    t_clamped = clamp(t, 0.0, 1.0)
    # Calculate the closest point on the segment.
    closest_point = a + t_clamped * ab
    # Return the distance from p to that closest point.
    return norm(p - closest_point)
end

const R_ENU2NED = @SMatrix [0.0 1.0 0.0;
                            1.0 0.0 0.0;
                            0.0 0.0 -1.0]

function gradient_uniform(y::AbstractVector{<:Real}, ts::Real)
    n = length(y)
    grad = Vector{Float64}(undef, n)
    if n == 0
        return grad
    elseif n == 1
        grad[1] = 0.0
        return grad
    end
    grad[1] = (y[2] - y[1]) / ts
    for i in 2:(n - 1)
        grad[i] = (y[i + 1] - y[i - 1]) / (2 * ts)
    end
    grad[n] = (y[n] - y[n - 1]) / ts
    return grad
end

function moving_average_same(x::AbstractVector{<:Real}, window::Int)
    n = length(x)
    if window <= 1 || n == 0
        return Float64.(x)
    end
    left = window ÷ 2
    right = window - 1 - left
    padded = Vector{Float64}(undef, n + left + right)
    padded[1:left] .= 0.0
    padded[(left + 1):(left + n)] .= x
    padded[(left + n + 1):end] .= 0.0
    out = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        s = 0.0
        for k in 0:(window - 1)
            s += padded[i + k]
        end
        out[i] = s / window
    end
    return out
end

"""
    compute_turn_radius(sl_in, sys::SystemStructure; smooth_window=10, eps=1e-12)

Compute signed turn radius from velocity and smoothed acceleration using the
instantaneous center of rotation (ICR). The sign follows the body x-axis
projected into the world xy-plane.

# Returns
- `(radius, time)`: Tuple of signed turn radius [m] and time vector [s]
- `nothing` if required data is missing
"""
function compute_turn_radius(sl_in, _sys::SystemStructure; smooth_window=10, eps=1e-12)
    sl = hasproperty(sl_in, :syslog) ? sl_in.syslog : sl_in
    n = length(sl.time)
    if n < 2 || isempty(sl.vel_kite) || isempty(sl.orient)
        return nothing
    end
    if length(sl.vel_kite) < n || length(sl.orient) < n
        return nothing
    end

    ts = mean(diff(sl.time))
    ts = isfinite(ts) && ts > eps ? ts : eps

    v_x = Vector{Float64}(undef, n)
    v_y = Vector{Float64}(undef, n)
    v_z = Vector{Float64}(undef, n)
    @inbounds for k in eachindex(v_x)
        v = sl.vel_kite[k]
        v_x[k] = v[1]
        v_y[k] = v[2]
        v_z[k] = v[3]
    end

    a_x = gradient_uniform(v_x, ts)
    a_y = gradient_uniform(v_y, ts)
    a_z = gradient_uniform(v_z, ts)

    if smooth_window > 1
        a_x = moving_average_same(a_x, smooth_window)
        a_y = moving_average_same(a_y, smooth_window)
        a_z = moving_average_same(a_z, smooth_window)
    end

    radius = Vector{Float64}(undef, n)
    @inbounds for k in eachindex(radius)
        v = SVector{3, Float64}(v_x[k], v_y[k], v_z[k])
        a = SVector{3, Float64}(a_x[k], a_y[k], a_z[k])
        v_norm = norm(v)
        if !isfinite(v_norm) || v_norm <= eps
            radius[k] = NaN
            continue
        end
        v_hat = v / v_norm
        a_t = dot(a, v_hat) * v_hat
        omega = cross(a - a_t, v) / (v_norm^2)
        omega_norm = norm(omega)
        if !isfinite(omega_norm) || omega_norm <= eps
            radius[k] = NaN
            continue
        end
        icr = cross(v, omega) / (omega_norm^2)
        R_b_w = SymbolicAWEModels.quaternion_to_rotation_matrix(sl.orient[k])
        e_x = SVector{3, Float64}(R_b_w[:, 1])
        det = e_x[1] * icr[2] - e_x[2] * icr[1]
        if !isfinite(det) || abs(det) <= eps
            radius[k] = NaN
        else
            radius[k] = -(det < 0 ? -1.0 : 1.0) * norm(icr)
        end
    end

    return radius, sl.time
end


"""
    Makie.plot(sys::SystemStructure, lg::SysLog; kwargs...)

Create a multi-panel plot of key simulation results from a `SysLog`.

# Arguments
- `sys::SystemStructure`: The system structure, used to get component counts (e.g., number of groups).
- `lg::SysLog`: The simulation log data to be plotted.

# Keyword Arguments
- `plot_default::Bool=true`: Defaults to true, enabling all plot panels. If false, all panels are disabled.
- `plot_turn_rates::Bool=false`: Show the panel with the wing's angular velocities (ω_x, ω_y, ω_z).
- `plot_turn_radius::Bool=false`: Show signed turn radius derived from velocity and acceleration.
- `plot_reelout::Bool=plot_default`: Show the panel with the reel-out velocities of the steering winches.
- `plot_aero_force::Bool=plot_default`: Show the panel with the z-component of aerodynamic force.
- `plot_aero_moment::Bool=false`: Show the panel with the y-component of aerodynamic moment.
- `plot_tether_moment::Bool=false`: Show the panel with the y-component of tether-induced moment.
- `plot_tether::Bool=false`: Show the panel with winched tether length.
- `plot_tether_actual::Bool=false`: Show the panel with actual tether length from nodal positions.
- `plot_twist::Bool=false`: Show the panel with the twist angles for each wing group.
- `plot_v_app`::Bool=false`: Show the panel with the apparent wind speed at the wing.
- `plot_aoa::Bool=plot_default`: Show the panel with the angle of attack.
- `aoa_ylims::Union{Nothing,Tuple}=nothing`: Y-axis limits for the AoA panel (`nothing` leaves autoscaling).
- `plot_heading::Bool=plot_default`: Show the panel with the kite's heading and course angles.
- `plot_kiteutils_course::Bool=false`: Also plot course calculated using KiteUtils.calc_course.
- `gk_ylims::Union{Nothing,Tuple}=(0.0, 10.0)`: Y-axis limits for the gk panel (`nothing` leaves autoscaling).
- `plot_elevation::Bool=false`: Show the panel with the kite's elevation angle.
- `plot_azimuth::Bool=false`: Show the panel with the kite's azimuth angle.
- `plot_distance::Bool=false`: Show the panel with the kite distance from origin (norm of position).
- `plot_yaw_rate::Bool=false`: Show yaw rate `dψ/dt` derived from the wind-referenced heading.
- `turn_radius_ylims::Union{Nothing,Tuple}=(0.0, 100.0)`: Y-axis limits for the turn-radius panel (`nothing` leaves autoscaling).
- `plot_cone_angle::Bool=false`: Show the panel with the cone angle (angle between wind vector and normalized kite position).
- `plot_old_heading::Bool=false`: Show the old heading calculated from orientation quaternion (angle between -R_b_w[:,1] and -R_v_w[:,1]).
- `plot_winch_force::Bool=plot_default`: Show the panel with the winch forces.
- `plot_set_values::Bool=false`: Show the panel with the set torque values.
- `suffix::String=" - " * sys.name`: Suffix to append to plot labels.
- `size::Tuple=(1200, 800)`: Figure size in pixels.

# Example
```julia
# Plot only the angle of attack and heading
plot(model.sys_struct, log, plot_reelout=false, plot_aero_force=false, plot_twist=false, plot_winch_force=false)
```
"""
function Makie.plot(sys::SystemStructure, lg::SysLog{N}; kwargs...) where N
    # Wrapper that uses the vector version
    return Makie.plot([sys], [lg]; kwargs...)
end

"""
    Makie.plot(sys::SystemStructure, logs::Vector{SysLog}; kwargs...)

Create a multi-panel plot comparing multiple simulation logs on the same figure.

This method allows plotting multiple syslogs (e.g., from PARTICLE_DYNAMICS and RIGID_DYNAMICS models)
on the same panels for direct comparison. Each log's traces are labeled with its
corresponding system name.

# Arguments
- `sys::SystemStructure`: The system structure (can be from any of the models).
- `logs::Vector{SysLog}`: Vector of simulation logs to compare.

# Keyword Arguments
Same as the single-syslog version. See `Makie.plot(sys::SystemStructure, lg::SysLog)` for details.

# Example
```julia
# Compare PARTICLE_DYNAMICS vs RIGID_DYNAMICS models
plot(sys_struct, [syslog_refine, syslog_quat];
     plot_turn_rates=true, plot_azimuth=false,
     plot_heading=false, plot_v_app=false, plot_aoa=false,
     plot_default=false, plot_aero_force=false)
```
"""
function Makie.plot(sys::SystemStructure, logs::Vector{<:SysLog}; kwargs...)
    # Wrapper that creates a vector of SystemStructure with the same sys for each log
    syss = [sys for _ in logs]
    return Makie.plot(syss, logs; kwargs...)
end

function Makie.plot(syss::Vector{<:SystemStructure}, logs::Vector{<:SysLog};
                   plot_default=true,
                   plot_reelout=plot_default,
                   plot_aero_force=plot_default,
                   plot_twist=false,
                   plot_us=false,
                   plot_gk=false,
                   gk_ylims=(0.0, 15.0),
                   plot_yaw_rate=false,
                   turn_radius_ylims=(0.0, 70.0),
                   plot_v_app=false,
                   plot_kite_vel=false,
                   plot_aoa=plot_default,
                   plot_sideslip=false,
                   plot_heading=plot_default,
                   plot_course=true,
                   plot_kiteutils_course=false,
                   plot_aero_moment=false,
                   plot_turn_rates=false,
                   plot_turn_radius=false,
                   plot_elevation=false,
                   plot_azimuth=false,
                   plot_wind=false,
                   plot_tether_moment=false,
                   plot_tether_actual=false,
                   plot_winch_force=plot_default,
                   plot_set_values=false,
                   plot_distance=false,
                   plot_cone_angle=false,
                   plot_old_heading=false,
                   plot_tether=false,
                   setpoints=nothing,  # Dict with keys matching plot names
                   ylims=nothing,  # Dict with keys matching plot names, values are (min, max)
                   suffixes=nothing,
                   size=(1200, 800),
                   show_legend=true,
                   ticklabelsize=12,
                   compact_labels=false,
                   legend_position=:rt,
                   legendsize=10)

    # Build list of panels to plot by combining data from all logs
    panels = []

    # Generate suffixes: use custom if provided, otherwise use system names
    # Use empty suffix when only one log (legend not shown anyway)
    actual_suffixes = if length(logs) == 1
        [""]
    elseif isnothing(suffixes)
        [" (" * sys.name * ")" for sys in syss]
    else
        [" (" * s * ")" for s in suffixes]
    end

    # Helper to get setpoint for a given key
    function get_setpoint(key::Symbol)
        isnothing(setpoints) && return nothing
        return get(setpoints, key, nothing)
    end

    # Helper to get ylim for a given key
    function get_ylim(key::Symbol)
        isnothing(ylims) && return nothing
        return get(ylims, key, nothing)
    end

    # Helper: combine a LaTeX symbol with a plain-text suffix inside one math
    # environment so MathTeXEngine renders it consistently.
    # e.g. lbl(L"\gamma", " (SymAWE)") → LaTeXString("$\gamma\text{ (SymAWE)}$")
    function lbl(var::LaTeXString, suffix::String)
        isempty(suffix) && return var
        inner = replace(String(var), r"^\$|\$$" => "")
        return LaTeXString("\$" * inner * "\\text{" * suffix * "}\$")
    end

    # Helper to add setpoint data to panel arrays and compute error
    function add_setpoint!(all_data, all_labels, all_times, all_linestyles,
                          setpoint_val, time_vec, label_base, suffix;
                          actual_data=nothing, error_info=nothing, skip_ref_suffix=false)
        isnothing(setpoint_val) && return
        push!(all_data, setpoint_val)
        label = skip_ref_suffix ? label_base * suffix : label_base * " ref" * suffix
        push!(all_labels, label)
        push!(all_times, time_vec)
        push!(all_linestyles, :dash)
        # Compute error at last timestep if actual_data provided
        if !isnothing(actual_data) && !isnothing(error_info) && length(actual_data) > 0
            actual_last = actual_data[end]
            ref_last = setpoint_val[end]
            abs_err = actual_last - ref_last
            rel_err = abs(ref_last) > 1e-6 ? 100.0 * abs_err / abs(ref_last) : NaN
            push!(error_info, (label=label_base, abs=abs_err, rel=rel_err))
        end
    end

    if plot_yaw_rate
        all_data = []
        all_labels = []
        all_times = []
        for (i, lg) in enumerate(logs)
            sl = lg.syslog
            suffix = actual_suffixes[i]

            # Calculate heading rate from diff for quaternion wings
            heading_unwrapped = copy(sl.heading)
            for j in 2:lastindex(heading_unwrapped)
                while heading_unwrapped[j] - heading_unwrapped[j-1] > π
                    heading_unwrapped[j] -= 2π
                end
                while heading_unwrapped[j] - heading_unwrapped[j-1] < -π
                    heading_unwrapped[j] += 2π
                end
            end
            heading_rate = diff(rad2deg.(heading_unwrapped)) ./ diff(sl.time)
            # Filter out high jumps (> 200 °/s)
            for j in eachindex(heading_rate)
                if abs(heading_rate[j]) > 200
                    heading_rate[j] = j > 1 ? heading_rate[j-1] : 0.0
                end
            end
            push!(all_data, heading_rate)
            push!(all_labels, lbl(L"\dot{\psi}", suffix))
            push!(all_times, sl.time[1:end-1])

            # Also plot z-component (yaw rate ω_z) if available
            turn_rates_rad = hcat(sl.turn_rates...)
            if !all(iszero, turn_rates_rad)
                omega_z_deg = rad2deg.(turn_rates_rad[3, :])
                push!(all_data, omega_z_deg)
                push!(all_labels, lbl(L"\omega_z", suffix))
                push!(all_times, sl.time)
            end
        end
        push!(panels, (
            data = all_data,
            labels = all_labels,
            times = all_times,
            ylabel = L"\omega \; [°/s]"
        ))
    end

    if plot_turn_rates
        all_data = []
        all_labels = []
        all_times = []
        for (i, lg) in enumerate(logs)
            sl = lg.syslog
            suffix = " - " * syss[i].name

            # View-frame angular velocity components (ω_v) from log turn_rates
            turn_rates_rad = hcat(sl.turn_rates...)
            if !all(iszero, turn_rates_rad)
                push!(all_data, rad2deg.(turn_rates_rad[1, :]))
                push!(all_labels, "ω_v,x")
                push!(all_times, sl.time)

                push!(all_data, rad2deg.(turn_rates_rad[2, :]))
                push!(all_labels, "ω_v,y")
                push!(all_times, sl.time)

                push!(all_data, rad2deg.(turn_rates_rad[3, :]))
                push!(all_labels, "ω_v,z")
                push!(all_times, sl.time)
            elseif !isempty(sl.orient)
                # Fallback for PARTICLE_DYNAMICS logs: reconstruct ω_v from quaternions
                sys_i = syss[i]
                wing = sys_i.wings[1]
                n = length(sl.time)
                ωx = Vector{Float64}(undef, n - 1)
                ωy = Vector{Float64}(undef, n - 1)
                ωz = Vector{Float64}(undef, n - 1)
                @inbounds for k in 1:(n - 1)
                    dt = sl.time[k + 1] - sl.time[k] + eps()
                    R1 = SymbolicAWEModels.quaternion_to_rotation_matrix(sl.orient[k])
                    R2 = SymbolicAWEModels.quaternion_to_rotation_matrix(sl.orient[k + 1])
                    R_rel = R2 * R1'
                    trR = clamp((R_rel[1, 1] + R_rel[2, 2] + R_rel[3, 3] - 1) / 2, -1.0, 1.0)
                    angle = acos(trR)
                    if angle < 1e-9
                        axis = SVector{3, Float64}(0.0, 0.0, 0.0)
                    else
                        denom = 2 * sin(angle) + eps()
                        axis = SVector{3, Float64}(
                            (R_rel[3, 2] - R_rel[2, 3]) / denom,
                            (R_rel[1, 3] - R_rel[3, 1]) / denom,
                            (R_rel[2, 1] - R_rel[1, 2]) / denom,
                        )
                    end
                    ω_w = (angle / dt) .* axis
                    pos_w = wing_log_pos(sl, sys_i, wing, k)
                    e_x = SVector{3, Float64}(R1[:, 1])
                    R_v_w = SymbolicAWEModels.calc_R_v_to_w(pos_w, e_x)
                    ω_v = R_v_w' * ω_w
                    ωx[k], ωy[k], ωz[k] = ω_v
                end
                push!(all_data, rad2deg.(ωx))
                push!(all_labels, "ω_v,x")
                push!(all_times, sl.time[1:end-1])

                push!(all_data, rad2deg.(ωy))
                push!(all_labels, "ω_v,y")
                push!(all_times, sl.time[1:end-1])

                push!(all_data, rad2deg.(ωz))
                push!(all_labels, "ω_v,z")
                push!(all_times, sl.time[1:end-1])
            end
        end
        if !isempty(all_data)
            push!(panels, (
                data = all_data,
                labels = all_labels,
                times = all_times,
                ylabel = "ω_v [°/s]\nturn-rate"
            ))
        end
    end

    if plot_turn_radius
        all_data = []
        all_labels = []
        all_times = []
        for (i, lg) in enumerate(logs)
            suffix = actual_suffixes[i]
            result = compute_turn_radius(lg, syss[i])
            if result === nothing
                @warn "Missing velocity/orientation data in syslog; skipping plot_turn_radius for $suffix"
                continue
            end
            radius, times = result
            push!(all_data, radius)
            push!(all_labels, "r_turn" * suffix)
            push!(all_times, times)
        end
        if !isempty(all_data)
            push!(panels, (
                data = all_data,
                labels = all_labels,
                times = all_times,
                ylabel = "turn radius [m]",
                ylims = turn_radius_ylims
            ))
        end
    end

    if plot_reelout
        all_data = []
        all_labels = []
        all_times = []
        for (i, lg) in enumerate(logs)
            sl = lg.syslog
            suffix = actual_suffixes[i]
            v_ro = [[sl.v_reelout[i][j] for i in eachindex(sl.v_reelout)] for j in 1:3]
            for j in 1:3
                # Only plot if non-zero or if it's index 1
                if j == 1 || !all(iszero, v_ro[j])
                    push!(all_data, v_ro[j])
                    push!(all_labels, lbl(L"v_{ro,%$j}", suffix))
                    push!(all_times, sl.time)
                end
            end
        end
        if !isempty(all_data)
            push!(panels, (
                data = all_data,
                labels = all_labels,
                times = all_times,
                ylabel = L"v_{ro} \; [m/s]"
            ))
        end
    end

    if plot_tether
        all_data = []
        all_labels = []
        all_times = []
        all_linestyles = Symbol[]
        all_errors = []
        sp = get_setpoint(:l_tether)
        for (i, lg) in enumerate(logs)
            sl = lg.syslog
            suffix = actual_suffixes[i]
            l_tether = [length(sl.l_tether[j]) > 0 ? sl.l_tether[j][1] : 0.0
                        for j in eachindex(sl.l_tether)]
            push!(all_data, l_tether)
            push!(all_labels, lbl(L"l_t", suffix))
            push!(all_times, sl.time)
            push!(all_linestyles, :solid)
            add_setpoint!(all_data, all_labels, all_times, all_linestyles,
                         sp, sl.time, "l_tether", suffix;
                         actual_data=l_tether, error_info=all_errors)
        end
        push!(panels, (
            data = all_data, labels = all_labels, times = all_times,
            linestyles = all_linestyles, error_info = all_errors,
            ylabel = L"l_t \; [m]"
        ))
    end

    if plot_tether_actual
        all_data = []
        all_labels = []
        all_times = []
        for (i, lg) in enumerate(logs)
            sl = lg.syslog
            suffix = actual_suffixes[i]
            sys = syss[i]
            if isempty(sys.tethers)
                @warn "No tethers found in system; skipping plot_tether_actual for $suffix"
                continue
            end
            if isempty(sl.X) || isempty(sl.Y) || isempty(sl.Z)
                @warn "Missing position data in syslog; skipping plot_tether_actual for $suffix"
                continue
            end
            tether = sys.tethers[1]
            seg_idxs = tether.segment_idxs
            if isempty(seg_idxs)
                @warn "Tether has no segments; skipping plot_tether_actual for $suffix"
                continue
            end
            segs = sys.segments
            n_samples = length(sl.time)
            actual_len = Vector{Float64}(undef, n_samples)
            for j in 1:n_samples
                X = sl.X[j]
                Y = sl.Y[j]
                Z = sl.Z[j]
                total = 0.0
                for seg_idx in seg_idxs
                    seg = segs[seg_idx]
                    p1, p2 = seg.point_idxs
                    dx = X[p2] - X[p1]
                    dy = Y[p2] - Y[p1]
                    dz = Z[p2] - Z[p1]
                    total += sqrt(dx * dx + dy * dy + dz * dz)
                end
                actual_len[j] = total
            end
            push!(all_data, actual_len)
            push!(all_labels, "l_tether_actual" * suffix)
            push!(all_times, sl.time)
        end
        if !isempty(all_data)
            push!(panels, (
                data = all_data,
                labels = all_labels,
                times = all_times,
                ylabel = "actual tether length [m]"
            ))
        end
    end

    if plot_aero_force
        all_data = []
        all_labels = []
        all_times = []
        for (i, lg) in enumerate(logs)
            sl = lg.syslog
            suffix = actual_suffixes[i]
            aero_force_z = [sl.aero_force_b[i][3] for i in eachindex(sl.aero_force_b)]
            push!(all_data, aero_force_z)
            push!(all_labels, lbl(L"F_{aero,z}", suffix))
            push!(all_times, sl.time)
        end
        push!(panels, (
            data = all_data,
            labels = all_labels,
            times = all_times,
            ylabel = L"F_{a,z} \; [N]"
        ))
    end

    if plot_aero_moment
        all_data = []
        all_labels = []
        all_times = []
        for (i, lg) in enumerate(logs)
            sl = lg.syslog
            suffix = actual_suffixes[i]
            aero_moment_z = [sl.aero_moment_b[i][3] for i in eachindex(sl.aero_moment_b)]
            push!(all_data, aero_moment_z)
            push!(all_labels, lbl(L"M_{aero,z}", suffix))
            push!(all_times, sl.time)
        end
        push!(panels, (
            data = all_data,
            labels = all_labels,
            times = all_times,
            ylabel = L"M_{a,z} \; [Nm]"
        ))
    end

    if plot_tether_moment
        all_data = []
        all_labels = []
        all_times = []
        for (i, lg) in enumerate(logs)
            sl = lg.syslog
            suffix = actual_suffixes[i]
            tether_moment_y = [sl.tether_induced_moment[i][2] for i in eachindex(sl.tether_induced_moment)]
            push!(all_data, tether_moment_y)
            push!(all_labels, lbl(L"M_{tether,y}", suffix))
            push!(all_times, sl.time)
        end
        push!(panels, (
            data = all_data,
            labels = all_labels,
            times = all_times,
            ylabel = L"M_{t,y} \; [Nm]"
        ))
    end

    if plot_twist
        all_data = []
        all_labels = []
        all_times = []
        for (i, lg) in enumerate(logs)
            sl = lg.syslog
            suffix = actual_suffixes[i]
            n_groups = length(sl.twist_angles[1])
            if n_groups > 0
                twist_deg = rad2deg.(hcat(sl.twist_angles...))
                for j in 1:n_groups
                    push!(all_data, twist_deg[j, :])
                    push!(all_labels, lbl(L"\beta_%$j", suffix))
                    push!(all_times, sl.time)
                end
            end
        end
        if !isempty(all_data)
            push!(panels, (
                data = all_data,
                labels = all_labels,
                times = all_times,
                ylabel = L"\beta \; [°]"
            ))
        end
    end

    if plot_us
        all_data = []
        all_labels = []
        all_times = []
        for (i, lg) in enumerate(logs)
            sl = lg.syslog
            suffix = actual_suffixes[i]
            push!(all_data, collect(sl.steering))
            push!(all_labels, lbl(L"u_s", suffix))
            push!(all_times, sl.time)
        end
        push!(panels, (
            data = all_data,
            labels = all_labels,
            times = all_times,
            ylabel = L"u_s \; [\%]"
        ))
    end

    if plot_gk
        all_data = []
        all_labels = []
        all_times = []
        for (i, lg) in enumerate(logs)
            sl = lg.syslog
            suffix = actual_suffixes[i]

            us_pct = collect(sl.steering)

            # Calculate heading rate from diff for quaternion wings
            heading_unwrapped = copy(sl.heading)
            for j in 2:lastindex(heading_unwrapped)
                while heading_unwrapped[j] - heading_unwrapped[j-1] > π
                    heading_unwrapped[j] -= 2π
                end
                while heading_unwrapped[j] - heading_unwrapped[j-1] < -π
                    heading_unwrapped[j] += 2π
                end
            end
            heading_rate = diff(rad2deg.(heading_unwrapped)) ./ diff(sl.time)
            v_app = sl.v_app[2:end]
            us_seg = us_pct[2:end]

            # calculate gk, guarding against zero steering (use percentage threshold)
            gk = similar(heading_rate)
            @inbounds for k in eachindex(gk)
                gk[k] = abs(us_seg[k]) > 1.0 ? heading_rate[k] / (v_app[k] * us_seg[k] / 100.0) : NaN
            end

            @info "turn-rate $(heading_rate[end])"
            @info "v_app $(v_app[end])"
            @info "us_seg $(us_seg[end])%"
            @info "gk $(gk[end])"
            @info "alpha-vsm $(rad2deg(sl.AoA[end]))"

            # Average over the last 50 seconds (or full log if shorter)
            t_end = sl.time[end]
            t_start = max(sl.time[1], t_end - 50)
            window_mean(vals, times) = begin
                idxs = findall(t -> t >= t_start, times)
                isempty(idxs) && return NaN
                window = vals[idxs]
                finite_vals = filter(isfinite, window)
                isempty(finite_vals) ? NaN : mean(finite_vals)
            end

            heading_rate_avg = window_mean(heading_rate, sl.time[2:end])
            v_app_avg = window_mean(v_app, sl.time[2:end])
            us_seg_avg = window_mean(us_seg, sl.time[2:end])
            gk_avg = window_mean(gk, sl.time[2:end])
            cs_us_avg = isempty(cs_over_us_vec) ? NaN : window_mean(cs_over_us_vec, sl.time[2:end])
            aoa_avg = window_mean(rad2deg.(sl.AoA), sl.time)

            @info "averages over last 50 seconds for log $(i) $(syss[i].name): \n"
            @info "turn-rate $(heading_rate_avg)"
            @info "v_app $(v_app_avg)"
            @info "us_seg $(us_seg_avg)"
            @info "gk ~10 $(gk_avg)"
            @info "cs/us ~0.1 $(abs(cs_us_avg))"
            @info "alpha-vsm $(aoa_avg)"

            # @info "--- resolving alpha mystery ---"
            # # ss.AoA = atan(wing.va_b[3], wing.va_b[1]) # version-1 
            # #---> ss.AoA = wing.vsm_solver.sol.alpha_dist[length(wing.vsm_solver.sol.alpha_dist) ÷ 2 + (length(wing.vsm_solver.sol.alpha_dist) % 2)] # version-2, likely with induction
            # # ss.AoA =wing.vsm_aero.alpha_uncorrected[length(wing.vsm_solver.sol.alpha_dist) ÷ 2 + (length(wing.vsm_solver.sol.alpha_dist) % 2)] # version-3, hopefully without induction
            # @info "alpha VSM (with induction?) $(rad2deg(sl.AoA[end])) deg"

            # # computing alpha geometrically
            # # Report final geometric AoA using hardcoded mid-panel corners (world frame)
            # last_state = sl[end]
            # X = last_state.X; Y = last_state.Y; Z = last_state.Z
            # # Mid-panel corners: 10,11,12,13 (11/13 front; 10/12 back)
            # back = 0.5 .* ([X[10], Y[10], Z[10]] .+ [X[12], Y[12], Z[12]])
            # front = 0.5 .* ([X[11], Y[11], Z[11]] .+ [X[13], Y[13], Z[13]])

            # delta_z = front[3] - back[3]
            # delta_x = front[1] - back[1]
            # aoa_wrt_horizontal = -rad2deg(atan(delta_z, delta_x))
            # @info "alpha wrt horizontal $(round(aoa_wrt_horizontal, digits=2)) deg"

            # mid_panel_vector = front .- back
            # mid_panel_vector_unit = mid_panel_vector / (norm(mid_panel_vector) + 1e-12)
            # # @info "mid-panel vector" mid_panel_vector_unit=round.(mid_panel_vector_unit, digits=5)
            
            # # wind vector in world frame
            # v_wind = sl.v_wind_kite[end]
            # # @info "v_wind_kite" v_wind=round.(v_wind, digits=5)
            # v_wind_unit = v_wind / (norm(v_wind) + 1e-12)

            # # compute angle v_a and vector_mid_panel
            # vel_KCU = sl.vel_kite[end]
            # # @info "vel_KCU" vel_kite=round.(vel_KCU, digits=5)
            # va_kcu = vel_KCU - v_wind
            # va_kcu_unit = va_kcu / (norm(va_kcu) + 1e-12)


            # # Flip chord direction so it points into the incoming flow (front -> back)
            # cos_theta = dot(-mid_panel_vector_unit, va_kcu_unit)
            # alpha_KCU = rad2deg(acos(clamp(cos_theta, -1.0, 1.0)))
            # @info "KCU" va_kcu=round.(va_kcu, digits=5) va_kcu_norm=norm(va_kcu) alpha_KCU=round(alpha_KCU, digits=2)
            # # @info "alpha wing (v_app_KCU) $(round(alpha_KCU, digits=2))"
            
            # # compute wing v_a
            # min1 = sl[end - 1]
            # last_state = sl[end]

            # X_last = last_state.X; Y_last = last_state.Y; Z_last = last_state.Z
            # X_min1 = min1.X; Y_min1 = min1.Y; Z_min1 = min1.Z

            # dt_last_to_min1 = last_state.time - min1.time + 1e-12
            # va_wing = SVector{3,Float64}(
            #     (X_last[1] - X_min1[1]) / (dt_last_to_min1) - v_wind[1],
            #     (Y_last[1] - Y_min1[1]) / (dt_last_to_min1) - v_wind[2],
            #     (Z_last[1] - Z_min1[1]) / (dt_last_to_min1) - v_wind[3],
            # )
            # # @info "v_app wing" va_wing=round.(va_wing, digits=5)
            # va_wing_unit = va_wing / (norm(va_wing) + 1e-12)
            

            # # Use the same convention: chord points front -> back, apparent wind approaches from front
            # cos_theta_wing = dot(-mid_panel_vector_unit, va_wing_unit)
            # alpha_wing = rad2deg(acos(clamp(cos_theta_wing, -1.0, 1.0)))
            # @info "WING" va_wing=round.(va_wing, digits=5) va_wing_norm=norm(va_wing) alpha_wing=round(alpha_wing, digits=2)

            # # computing lift and drag using the total aero force "aero_force_b"
            # # SysLog stores orientation as a quaternion; rebuild R_b_w on the fly
            # R_b_w = SymbolicAWEModels.quaternion_to_rotation_matrix(sl.orient[end])
            # F_aero_b = sl.aero_force_b[end]
            # F_aero_world = R_b_w * F_aero_b
            # # Decompose aero force into drag (opposing apparent wind) and lift (perpendicular)
            # drag_dir = -va_wing_unit               # drag acts against the flow
            # drag = -dot(F_aero_world, va_wing_unit)  # positive magnitude
            # drag_vec = drag * drag_dir
            # lift_vec = F_aero_world - dot(F_aero_world, va_wing_unit) * va_wing_unit
            # lift = norm(lift_vec)
            # lift_dir = lift > 1e-12 ? lift_vec / lift : zeros(3)
            # @info "Aero VSM forces" lift=round(lift, digits=2) drag=round(drag, digits=2) L_over_D=round(lift / (drag + 1e-12), digits=2)

            # # Aero forces of tethers
            # tether_force_w = sl.tether_induced_force[end]
            # drag_tether = -dot(tether_force_w, va_wing_unit)
            # tether_lift_vec = tether_force_w - dot(tether_force_w, va_wing_unit) * va_wing_unit
            # tether_lift = norm(tether_lift_vec)
            # @info "Aero tether forces" lift=round(tether_lift, digits=2) drag=round(drag_tether, digits=2) L_over_D=round(tether_lift / (drag_tether + 1e-12), digits=2)

            # # Total aero forces (wing + tether)
            # total_drag = drag + drag_tether
            # total_lift_vec = lift_vec + tether_lift_vec
            # total_lift = norm(total_lift_vec)
            # total_angle = rad2deg(acos(clamp(dot(total_lift_vec / (total_lift + 1e-12), drag_dir), -1.0, 1.0)))
            # @info "Aero total forces" lift=round(total_lift, digits=2) drag=round(total_drag, digits=2) angle_lift_to_drag=round(total_angle, digits=2) L_over_D=round(total_lift / (total_drag + 1e-12), digits=2)

            push!(all_data, gk)
            push!(all_labels, lbl(L"g_k", suffix))
            push!(all_times, sl.time[2:end])
        end
        push!(panels, (
            data = all_data,
            labels = all_labels,
            times = all_times,
            ylabel = L"g_k \; [-]"
        ))
    end

    if plot_v_app
        all_data = []
        all_labels = []
        all_times = []
        all_linestyles = Symbol[]
        all_errors = []
        sp = get_setpoint(:v_app)
        for (i, lg) in enumerate(logs)
            sl = lg.syslog
            suffix = actual_suffixes[i]
            push!(all_data, sl.v_app)
            push!(all_labels, lbl(L"v_a", suffix))
            push!(all_times, sl.time)
            push!(all_linestyles, :solid)
            add_setpoint!(all_data, all_labels, all_times, all_linestyles,
                         sp, sl.time, "v_app", suffix;
                         actual_data=sl.v_app, error_info=all_errors)
        end
        push!(panels, (
            data = all_data, labels = all_labels, times = all_times,
            linestyles = all_linestyles, error_info = all_errors,
            ylabel = L"v_a \; [m/s]", ylim = get_ylim(:v_app)
        ))
    end

    if plot_kite_vel
        all_data = []
        all_labels = []
        all_times = []
        all_linestyles = Symbol[]
        all_errors = []
        sp = get_setpoint(:v_kite)
        for (i, lg) in enumerate(logs)
            sl = lg.syslog
            suffix = actual_suffixes[i]
            v_kite_norm = [norm(v) for v in sl.vel_kite]
            push!(all_data, v_kite_norm)
            push!(all_labels, lbl(L"v_k", suffix))
            push!(all_times, sl.time)
            push!(all_linestyles, :solid)
            add_setpoint!(all_data, all_labels, all_times, all_linestyles,
                         sp, sl.time, "|v_kite|", suffix;
                         actual_data=v_kite_norm, error_info=all_errors)
        end
        push!(panels, (
            data = all_data, labels = all_labels, times = all_times,
            linestyles = all_linestyles, error_info = all_errors,
            ylabel = L"v_k \; [m/s]"
        ))
    end

    if plot_aoa
        all_data = []
        all_labels = []
        all_times = []
        all_linestyles = Symbol[]
        all_errors = []
        sp = get_setpoint(:aoa)
        for (i, lg) in enumerate(logs)
            sl = lg.syslog
            suffix = actual_suffixes[i]
            aoa_deg = rad2deg.(sl.AoA)
            push!(all_data, aoa_deg)
            push!(all_labels, lbl(L"\alpha", suffix))
            push!(all_times, sl.time)
            push!(all_linestyles, :solid)
            if !isnothing(sp)
                sp_deg = rad2deg.(sp)
                push!(all_data, sp_deg)
                push!(all_labels, lbl(L"\alpha_{ref}", suffix))
                push!(all_times, sl.time)
                push!(all_linestyles, :dash)
                # Compute error
                abs_err = aoa_deg[end] - sp_deg[end]
                rel_err = abs(sp_deg[end]) > 1e-6 ? 100.0 * abs_err / abs(sp_deg[end]) : NaN
                push!(all_errors, (label="α", abs=abs_err, rel=rel_err))
            end
        end
        push!(panels, (
            data = all_data, labels = all_labels, times = all_times,
            linestyles = all_linestyles, error_info = all_errors,
            ylabel = L"\alpha \; [°]", ylim = get_ylim(:aoa)
        ))
    end

    if plot_sideslip
        all_data = []
        all_labels = []
        all_times = []
        for (i, lg) in enumerate(logs)
            sl = lg.syslog
            suffix = actual_suffixes[i]
            sideslip_deg = rad2deg.(sl.side_slip)
            push!(all_data, sideslip_deg)
            push!(all_labels, lbl(L"\beta", suffix))
            push!(all_times, sl.time)
        end
        push!(panels, (
            data = all_data,
            labels = all_labels,
            times = all_times,
            ylabel = L"\beta \; [°]"
        ))
    end

    if plot_heading
        all_data = []
        all_labels = []
        all_times = []
        all_linestyles = Symbol[]
        all_errors = []
        sp = get_setpoint(:heading)
        for (i, lg) in enumerate(logs)
            sl = lg.syslog
            suffix = actual_suffixes[i]
            heading_deg = rad2deg.(sl.heading)
            push!(all_data, heading_deg)
            push!(all_labels, lbl(L"\psi", suffix))
            push!(all_times, sl.time)
            push!(all_linestyles, :solid)
            if plot_course
                course_deg = rad2deg.(sl.course)
                push!(all_data, course_deg)
                push!(all_labels, lbl(L"\chi", suffix))
                push!(all_times, sl.time)
                push!(all_linestyles, :solid)
            end
            if plot_kiteutils_course
                course_kiteutils_deg = [rad2deg(KiteUtils.calc_course(sl.vel_kite[i], sl.elevation[i], sl.azimuth[i]))
                                        for i in eachindex(sl.orient)]
                push!(all_data, course_kiteutils_deg)
                push!(all_labels, lbl(L"\chi_{KU}", suffix))
                push!(all_times, sl.time)
                push!(all_linestyles, :solid)
            end
            if !isnothing(sp)
                sp_deg = rad2deg.(sp)
                push!(all_data, sp_deg)
                push!(all_labels, lbl(L"\psi_{ref}", suffix))
                push!(all_times, sl.time)
                push!(all_linestyles, :dash)
                # Compute error
                abs_err = heading_deg[end] - sp_deg[end]
                rel_err = abs(sp_deg[end]) > 1e-6 ? 100.0 * abs_err / abs(sp_deg[end]) : NaN
                push!(all_errors, (label="ψ", abs=abs_err, rel=rel_err))
            end
        end
        ylabel_heading = plot_course ? L"\psi, \chi \; [°]" : L"\psi \; [°]"
        push!(panels, (
            data = all_data, labels = all_labels, times = all_times,
            linestyles = all_linestyles, error_info = all_errors,
            ylabel = ylabel_heading, ylim = get_ylim(:heading)
        ))
    end

    if plot_old_heading
        all_data = []
        all_labels = []
        all_times = []
        for (i, lg) in enumerate(logs)
            sl = lg.syslog
            suffix = actual_suffixes[i]
            # Calculate old heading from orientation quaternion
            old_heading_rad = similar(sl.heading)
            for i in eachindex(sl.orient)
                R_b_w = SymbolicAWEModels.quaternion_to_rotation_matrix(sl.orient[i])
                R_v_w = [1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 -1.0]
                v1 = -R_b_w[:, 1]
                v2 = -R_v_w[:, 1]
                old_heading_rad[i] = atan(v1[2] - v2[2], v1[1] - v2[1])
            end
            old_heading_deg = rad2deg.(old_heading_rad)
            push!(all_data, old_heading_deg)
            push!(all_labels, lbl(L"\psi_{old}", suffix))
            push!(all_times, sl.time)
        end
        push!(panels, (
            data = all_data,
            labels = all_labels,
            times = all_times,
            ylabel = L"\psi_{old} \; [°]"
        ))
    end

    if plot_distance
        all_data = []
        all_labels = []
        all_times = []
        for (i, lg) in enumerate(logs)
            sl = lg.syslog
            suffix = actual_suffixes[i]
            sys_i = syss[i]
            wing = sys_i.wings[1]
            distance = [norm(wing_log_pos(sl, sys_i, wing, j)) for j in eachindex(sl.X)]
            push!(all_data, distance)
            push!(all_labels, lbl(L"r", suffix))
            push!(all_times, sl.time)
        end
        push!(panels, (
            data = all_data,
            labels = all_labels,
            times = all_times,
            ylabel = L"r \; [m]"
        ))
    end

    if plot_cone_angle
        all_data = []
        all_labels = []
        all_times = []
        for (i, lg) in enumerate(logs)
            sl = lg.syslog
            suffix = actual_suffixes[i]
            sys_i = syss[i]
            wing = sys_i.wings[1]
            # Assuming wind_vec_gnd is available in syslog
            cone_angle_rad = similar(sl.heading)
            for j in eachindex(sl.X)
                pos = wing_log_pos(sl, sys_i, wing, j)
                pos_norm = normalize(pos)
                wind_norm = normalize([sl.v_wind_gnd[1], sl.v_wind_gnd[2], sl.v_wind_gnd[3]])
                cone_angle_rad[j] = acos(clamp(dot(pos_norm, wind_norm), -1.0, 1.0))
            end
            cone_angle_deg = rad2deg.(cone_angle_rad)
            push!(all_data, cone_angle_deg)
            push!(all_labels, lbl(L"\theta_c", suffix))
            push!(all_times, sl.time)
        end
        push!(panels, (
            data = all_data,
            labels = all_labels,
            times = all_times,
            ylabel = L"\theta_c \; [°]"
        ))
    end

    if plot_elevation
        all_data = []
        all_labels = []
        all_times = []
        all_linestyles = Symbol[]
        all_errors = []
        sp = get_setpoint(:elevation)
        for (i, lg) in enumerate(logs)
            sl = lg.syslog
            suffix = actual_suffixes[i]
            elevation_deg = rad2deg.(sl.elevation)
            push!(all_data, elevation_deg)
            push!(all_labels, lbl(L"\gamma", suffix))
            push!(all_times, sl.time)
            push!(all_linestyles, :solid)
            if !isnothing(sp)
                sp_deg = rad2deg.(sp)
                push!(all_data, sp_deg)
                push!(all_labels, lbl(L"\gamma_{ref}", suffix))
                push!(all_times, sl.time)
                push!(all_linestyles, :dash)
                # Compute error
                abs_err = elevation_deg[end] - sp_deg[end]
                rel_err = abs(sp_deg[end]) > 1e-6 ? 100.0 * abs_err / abs(sp_deg[end]) : NaN
                push!(all_errors, (label="elev", abs=abs_err, rel=rel_err))
            end
        end
        push!(panels, (
            data = all_data,
            labels = all_labels,
            times = all_times,
            linestyles = all_linestyles,
            error_info = all_errors,
            ylabel = L"\gamma \; [°]"
        ))
    end

    if plot_azimuth
        all_data = []
        all_labels = []
        all_times = []
        all_linestyles = Symbol[]
        all_errors = []
        sp = get_setpoint(:azimuth)
        for (i, lg) in enumerate(logs)
            sl = lg.syslog
            suffix = actual_suffixes[i]
            azimuth_deg = rad2deg.(sl.azimuth)
            push!(all_data, azimuth_deg)
            push!(all_labels, lbl(L"\phi", suffix))
            push!(all_times, sl.time)
            push!(all_linestyles, :solid)
            if !isnothing(sp)
                sp_deg = rad2deg.(sp)
                push!(all_data, sp_deg)
                push!(all_labels, lbl(L"\phi_{ref}", suffix))
                push!(all_times, sl.time)
                push!(all_linestyles, :dash)
                # Compute error
                abs_err = azimuth_deg[end] - sp_deg[end]
                rel_err = abs(sp_deg[end]) > 1e-6 ? 100.0 * abs_err / abs(sp_deg[end]) : NaN
                push!(all_errors, (label="azim", abs=abs_err, rel=rel_err))
            end
        end
        push!(panels, (
            data = all_data,
            labels = all_labels,
            times = all_times,
            linestyles = all_linestyles,
            error_info = all_errors,
            ylabel = L"\phi \; [°]"
        ))
    end

    if plot_wind
        all_data = []
        all_labels = []
        all_times = []
        for (i, lg) in enumerate(logs)
            sl = lg.syslog
            suffix = actual_suffixes[i]
            # v_wind_gnd is a vector, take norm
            wind_speed = [norm(sl.v_wind_gnd[j]) for j in eachindex(sl.v_wind_gnd)]
            push!(all_data, wind_speed)
            push!(all_labels, lbl(L"v_w", suffix))
            push!(all_times, sl.time)
        end
        push!(panels, (
            data = all_data,
            labels = all_labels,
            times = all_times,
            ylabel = L"v_w \; [m/s]"
        ))
    end

    if plot_winch_force
        all_data = []
        all_labels = []
        all_times = []
        all_linestyles = Symbol[]
        all_errors = []
        sp = get_setpoint(:winch_force)
        for (i, lg) in enumerate(logs)
            sl = lg.syslog
            suffix = actual_suffixes[i]
            winch_force = [[sl.winch_force[i][j] for i in eachindex(sl.winch_force)] for j in 1:3]
            for j in 1:3
                if j == 1 || !all(iszero, winch_force[j])
                    push!(all_data, winch_force[j])
                    if compact_labels
                        push!(all_labels, "Winch force" * suffix)
                    else
                        push!(all_labels, lbl(L"F_{t,%$j}", suffix))
                    end
                    push!(all_times, sl.time)
                    push!(all_linestyles, :solid)
                end
            end
            # Add setpoint only for first winch (use winch_force[1] for error)
            setpoint_label = compact_labels ? "Ref. winch force" : "F_winch ref"
            add_setpoint!(all_data, all_labels, all_times, all_linestyles,
                         sp, sl.time, setpoint_label, suffix;
                         actual_data=winch_force[1], error_info=all_errors, skip_ref_suffix=compact_labels)
        end
        if !isempty(all_data)
            ylabel_winch = compact_labels ? L"[N]" : L"F_t \; [N]"
            push!(panels, (
                data = all_data, labels = all_labels, times = all_times,
                linestyles = all_linestyles, error_info = all_errors,
                ylabel = ylabel_winch, ylim = get_ylim(:winch_force)
            ))
        end
    end

    if plot_set_values
        all_data = []
        all_labels = []
        all_times = []
        for (i, lg) in enumerate(logs)
            sl = lg.syslog
            suffix = actual_suffixes[i]
            set_values = [[sl.set_torque[i][j] for i in eachindex(sl.set_torque)] for j in 1:3]
            for j in 1:3
                # Only plot if non-zero or if it's index 1
                if j == 1 || !all(iszero, set_values[j])
                    push!(all_data, set_values[j])
                    push!(all_labels, lbl(L"\tau_%$j", suffix))
                    push!(all_times, sl.time)
                end
            end
        end
        if !isempty(all_data)
            push!(panels, (
                data = all_data,
                labels = all_labels,
                times = all_times,
                ylabel = L"\tau \; [Nm]"
            ))
        end
    end

    # Check if there's anything to plot
    if isempty(panels)
        error("No plot sections enabled. Enable at least one plot panel.")
    end

    # Create figure with subplots
    n_panels = length(panels)
    fig = Figure(size=size)

    axes = []
    label_fontsize = 16
    for (i, panel) in enumerate(panels)
        # Share x-axis with first subplot
        if i == 1
            ax = Axis(fig[i, 1], ylabel=panel.ylabel, ylabelsize=label_fontsize,
                      xticklabelsize=ticklabelsize, yticklabelsize=ticklabelsize)
        else
            ax = Axis(fig[i, 1], ylabel=panel.ylabel, ylabelsize=label_fontsize,
                      xticklabelsvisible=false,
                      xticklabelsize=ticklabelsize, yticklabelsize=ticklabelsize)
            linkxaxes!(axes[1], ax)
        end

        # Plot each data series in this panel
        linestyles = haskey(panel, :linestyles) ? panel.linestyles : fill(:solid, length(panel.data))
        for (j, (data_series, label, time_vec, ls)) in enumerate(zip(panel.data, panel.labels, panel.times, linestyles))
            if length(data_series) != length(time_vec)
                @warn "Skipping trace '$label': data length $(length(data_series)) != time length $(length(time_vec))"
                continue
            end
            lines!(ax, time_vec, data_series, label=label, linestyle=ls)
        end

        # Apply ylims if specified
        if haskey(panel, :ylim) && !isnothing(panel.ylim)
            Makie.ylims!(ax, panel.ylim...)
        end

        # Add legend if multiple traces and show_legend is true
        if show_legend && length(panel.data) > 1
            axislegend(ax, position=legend_position, labelsize=legendsize, patchsize=(10, 5))
        end

        # Display error info as text annotation with white background if present
        if haskey(panel, :error_info) && !isempty(panel.error_info)
            err_texts = String[]
            for err in panel.error_info
                abs_str = @sprintf("%.2f", err.abs)
                rel_str = isnan(err.rel) ? "N/A" : @sprintf("%.1f%%", err.rel)
                push!(err_texts, "Δ=$(abs_str) ($(rel_str))")
            end
            err_text = join(err_texts, "\n")
            text!(ax, 0.98, 0.02; text=err_text, space=:relative,
                  align=(:right, :bottom), fontsize=12, color=:black,
                  glowwidth=4.0, glowcolor=(:white, 1.0))
        end

        push!(axes, ax)
    end

    # Add x-label to bottom subplot
    axes[end].xlabel = L"t \; [s]"
    axes[end].xlabelsize = label_fontsize
    axes[end].xticklabelsvisible = true

    Makie.resize_to_layout!(fig)
    return fig
end

function zoom_out!(scene, cam, plots, distance=nothing; relmargin=0.2)
    # --- ROBUST ZOOM OUT ---
    # 1. Get the current camera viewing direction vector
    inv_view_matrix = inv(cam.view[])
    cam_dir_vec = normalize(Vec3f(inv_view_matrix[1, 3],
                                  inv_view_matrix[2, 3],
                                  inv_view_matrix[3, 3]))
    # 2. Get the scene's bounding box and its center (the new target)
    #    Expand to include origin so global axes arrows are always visible
    bbox = data_limits(plots)
    lo = min.(bbox.origin, Vec3f(0, 0, 0))
    hi = max.(bbox.origin .+ bbox.widths, Vec3f(0, 0, 0))
    bbox = Rect3f(lo, hi .- lo)
    center = lo .+ (hi .- lo) ./ 2

    # 3. Calculate distance only if not provided
    if isnothing(distance)
        # Calculate the distance needed to see the whole box
        radius = norm(bbox.widths) / 2.0
        fov_rad = 2 * atan(1 / cam.projection[][2, 2])
        distance = radius / tan(fov_rad / 2.0) * (1 + relmargin)
    end

    # 4. Calculate the new camera position
    new_eyepos = center + cam_dir_vec * distance
    # 5. Update the camera to the new "fit-all" view
    update_cam!(scene, new_eyepos, center)

    return distance
end

function zoom_in!(scene, cam, sys, segment_idx, distance=nothing)
    # --- ZOOM IN ON SEGMENT ---
    # Get current segment endpoints
    seg = sys.segments[segment_idx]
    p1_w = sys.points[seg.point_idxs[1]].pos_w
    p2_w = sys.points[seg.point_idxs[2]].pos_w

    # Calculate segment center
    center = (p1_w + p2_w) / 2.0f0

    # Calculate distance only if not provided
    if isnothing(distance)
        segment_len = norm(p2_w - p1_w)
        distance = segment_len * 1.5 + 2.0
    end

    # Get camera direction
    inv_view_matrix = inv(cam.view[])
    cam_dir_vec = normalize(Vec3f(inv_view_matrix[1, 3],
                                  inv_view_matrix[2, 3],
                                  inv_view_matrix[3, 3]))

    # Calculate new camera position
    new_eyepos = center + distance * cam_dir_vec

    # Update camera
    update_cam!(scene, new_eyepos, center)

    return distance
end

"""
    zoom_body_frame!(scene, cam, sys, distance=nothing)

Set camera to body-frame view tracking the wing orientation.
Camera positioned behind the kite looking forward along body x-axis.

# Arguments
- `scene`: The Makie scene
- `cam`: The camera object
- `sys`: SystemStructure with wing to track
- `distance`: Optional fixed distance to use (if nothing, calculates from geometry)

# Returns
- Camera distance used (for storage/reuse)
"""
function zoom_body_frame!(scene, cam, sys, distance=nothing)
    if isempty(sys.wings)
        @warn "No wings in system, cannot use body frame view"
        return nothing
    end

    wing = sys.wings[1]
    kite_pos = wing.pos_w
    R_b_w = wing.R_b_to_w

    # Calculate distance only if not provided
    if isnothing(distance)
        # Calculate characteristic system length
        if !isempty(sys.points)
            all_x = [p.pos_w[1] for p in sys.points]
            all_y = [p.pos_w[2] for p in sys.points]
            all_z = [p.pos_w[3] for p in sys.points]

            xlims = extrema(all_x)
            ylims = extrema(all_y)
            zlims = extrema(all_z)

            char_length = max(xlims[2] - xlims[1], ylims[2] - ylims[1], zlims[2] - zlims[1])
        else
            char_length = 10.0
        end
        distance = char_length * 0.1
    end

    # Camera position: kite_pos - R_b_w * [distance, 0, 0]
    # Places camera in front of kite (negative x in body frame), looking along +x axis
    cam_offset_body = [-distance, 0.0, 0.0]
    cam_offset_world = R_b_w * cam_offset_body
    cam_pos = kite_pos + cam_offset_world

    update_cam!(scene, Vec3f(cam_pos), Vec3f(kite_pos), Vec3f(R_b_w[:, 3]))

    return distance
end

function get_cam_eyepos(cam)
    inv_view = inv(cam.view[])
    return Vec3f(inv_view[1, 4], inv_view[2, 4], inv_view[3, 4])
end

"""
    zoom_in_body!(scene, cam, sys, body_idx, distance=nothing)

Center the camera on rigid body `body_idx`. Without `distance`, derives one from
the body's spoke reach (joint anchors) so the body fills the view.
"""
function zoom_in_body!(scene, cam, sys, body_idx, distance=nothing)
    body = sys.bodies[body_idx]
    center = Vec3f(body.pos_w)
    if isnothing(distance)
        R = SymbolicAWEModels.quaternion_to_rotation_matrix(body.Q_b_to_w)
        reach = 0.0
        for joint in sys.elastic_joints
            joint.body_a_idx == body_idx &&
                (reach = max(reach, norm(R * joint.anchor_a_b)))
            joint.body_b_idx == body_idx &&
                (reach = max(reach, norm(R * joint.anchor_b_b)))
        end
        distance = reach * 4.0 + 2.0
    end
    inv_view = inv(cam.view[])
    cam_dir = normalize(Vec3f(inv_view[1, 3], inv_view[2, 3], inv_view[3, 3]))
    update_cam!(scene, center + distance * cam_dir, center)
    return distance
end

function apply_zoom_mode!(scene, relevant_plots, stored_sys; mode_changed::Bool)
    cam = scene.camera

    if mode_changed
        if PLOT_PREV_BODY_FRAME[] && !isempty(stored_sys.wings)
            wing = stored_sys.wings[1]
            PLOT_BODY_DISTANCE[] = norm(get_cam_eyepos(cam) - Vec3f(wing.pos_w))
        elseif !PLOT_PREV_ZOOMED_IN[] && !PLOT_PREV_BODY_FRAME[]
            cc = Makie.cameracontrols(scene)
            PLOT_WORLD_EYEPOS[] = Vec3f(cc.eyeposition[])
            PLOT_WORLD_LOOKAT[] = Vec3f(cc.lookat[])
        end
    end

    if PLOT_BODY_FRAME[]
        wing_pos = isempty(stored_sys.wings) ? nothing :
                   Vec3f(stored_sys.wings[1].pos_w)
        if mode_changed || isnothing(PLOT_BODY_PREV_WING_POS[])
            dist = zoom_body_frame!(scene, cam, stored_sys, PLOT_BODY_DISTANCE[])
        else
            current_dist = norm(get_cam_eyepos(cam) - PLOT_BODY_PREV_WING_POS[])
            dist = zoom_body_frame!(scene, cam, stored_sys, current_dist)
        end
        PLOT_BODY_DISTANCE[] = dist
        PLOT_CAMERA_DISTANCE[] = dist
        PLOT_BODY_PREV_WING_POS[] = wing_pos
    elseif PLOT_ZOOMED_IN[] && PLOT_ZOOM_SEGMENT_IDX[] > 0
        dist = zoom_in!(scene, cam, stored_sys, PLOT_ZOOM_SEGMENT_IDX[], PLOT_CAMERA_DISTANCE[])
        PLOT_CAMERA_DISTANCE[] = dist
    elseif PLOT_ZOOMED_IN[] && PLOT_ZOOM_BODY_IDX[] > 0
        dist = zoom_in_body!(scene, cam, stored_sys, PLOT_ZOOM_BODY_IDX[], PLOT_CAMERA_DISTANCE[])
        PLOT_CAMERA_DISTANCE[] = dist
    elseif !PLOT_ZOOMED_IN[]
        if mode_changed
            if !isnothing(PLOT_WORLD_EYEPOS[]) && !isnothing(PLOT_WORLD_LOOKAT[])
                update_cam!(scene, PLOT_WORLD_EYEPOS[], PLOT_WORLD_LOOKAT[])
            else
                dist = zoom_out!(scene, cam, relevant_plots, nothing;
                                 relmargin=PLOT_ZOOM_RELMARGIN[])
                PLOT_CAMERA_DISTANCE[] = dist
            end
        end
    end

    if mode_changed
        PLOT_PREV_BODY_FRAME[] = PLOT_BODY_FRAME[]
        PLOT_PREV_ZOOMED_IN[] = PLOT_ZOOMED_IN[]
        PLOT_PREV_SEGMENT_IDX[] = PLOT_ZOOM_SEGMENT_IDX[]
    end

    return nothing
end

"""
    hover_ref_label(name, idx)

Hover label for a point/segment: its `name` ref (e.g. `bridle_1`) when set,
falling back to the integer `idx` for auto-generated components with no name.
"""
hover_ref_label(name, idx) = isnothing(name) ? string(idx) : string(name)

"""
    setup_segment_hover_events!(scene, systems::Vector{<:SystemStructure},
                                 segment_color_obs::Vector, all_plots;
                                 segment_colors, highlight_color=:yellow,
                                 force_color=false, relmargin=0.2)

Add hover labels and click-to-zoom for segments across multiple systems.

`segment_color_obs[i]` is system `i`'s per-segment color observable (the one
driving the `linesegments!` colors); the highlight is applied by writing to it.
Segment endpoints are read from the live `SystemStructure`, not the plot.

For multi-system, labels show "sys_idx:seg_ref" and "sys_idx:pt_ref" format.
For single system, labels show just "seg_ref" and "pt_ref".
"""
function setup_segment_hover_events!(scene, systems::Vector{<:SystemStructure},
                                      segment_color_obs::Vector, all_plots;
                                      segment_colors::Vector,
                                      highlight_color=:yellow,
                                      force_color=false,
                                      relmargin=0.2)
    # Track hover state per system
    last_hovered = Ref((-1, -1))  # (sys_idx, seg_idx)
    zoomed_in = Ref(false)

    # Create label observables
    segment_label = Observable("")
    segment_label_pos = Observable(Point2f(0, 0))
    segment_label_visible = Observable(false)
    text!(scene, segment_label, position=segment_label_pos, space=:pixel,
          fontsize=14, color=:white, strokecolor=:black, strokewidth=1,
          align=(:center, :center), visible=segment_label_visible, transparency=true)

    point1_label = Observable("")
    point1_label_pos = Observable(Point2f(0, 0))
    point1_label_visible = Observable(false)
    text!(scene, point1_label, position=point1_label_pos, space=:pixel,
          fontsize=14, color=:white, strokecolor=:black, strokewidth=1,
          align=(:center, :center), visible=point1_label_visible, transparency=true)

    point2_label = Observable("")
    point2_label_pos = Observable(Point2f(0, 0))
    point2_label_visible = Observable(false)
    text!(scene, point2_label, position=point2_label_pos, space=:pixel,
          fontsize=14, color=:white, strokecolor=:black, strokewidth=1,
          align=(:center, :center), visible=point2_label_visible, transparency=true)

    # --- Hover handler: find closest segment across ALL systems ---
    on(events(scene).mouseposition, priority=2) do mp
        min_dist = Inf
        closest = (-1, -1)  # (sys_idx, seg_idx)
        mouse_pos_2d = Point2f(mp)
        margin_px = 30.0

        for (sys_i, sys) in enumerate(systems)
            for seg_i in eachindex(sys.segments)
                seg = sys.segments[seg_i]
                p1_2d = Makie.project(scene, Point3f(sys.points[seg.point_idxs[1]].pos_w))
                p2_2d = Makie.project(scene, Point3f(sys.points[seg.point_idxs[2]].pos_w))
                dist = point_line_segment_distance(mouse_pos_2d, p1_2d, p2_2d)
                if dist < min_dist
                    min_dist = dist
                    closest = (sys_i, seg_i)
                end
            end
        end

        hover = (min_dist < margin_px) ? closest : (-1, -1)
        if hover != last_hovered[]
            # Reset all segment colors (per-segment obs drives the mesh vertices)
            for (sys_i, sys) in enumerate(systems)
                base_color = segment_colors[sys_i]
                if force_color
                    new_colors = calculate_segment_force_colors(sys.segments, base_color)
                else
                    new_colors = fill(to_color(base_color), length(sys.segments))
                end
                if hover[1] == sys_i && hover[2] != -1
                    new_colors[hover[2]] = to_color(highlight_color)
                end
                segment_color_obs[sys_i][] = new_colors
            end

            if hover != (-1, -1)
                sys_i, seg_i = hover
                sys = systems[sys_i]
                seg = sys.segments[seg_i]
                p1_3d = sys.points[seg.point_idxs[1]].pos_w
                p2_3d = sys.points[seg.point_idxs[2]].pos_w

                # Multi-system label format: "1:seg_ref" or just "seg_ref" if single system
                prefix = length(systems) > 1 ? "$(sys_i):" : ""
                p1 = sys.points[seg.point_idxs[1]]
                p2 = sys.points[seg.point_idxs[2]]

                mid_point_2d = Makie.project(scene, (p1_3d + p2_3d) / 2)
                segment_label[] = "$(prefix)$(hover_ref_label(seg.name, seg_i))"
                segment_label_pos[] = mid_point_2d + Point2f(20, 0)
                segment_label_visible[] = true

                p1_2d = Makie.project(scene, p1_3d)
                p2_2d = Makie.project(scene, p2_3d)
                point1_label[] = "$(prefix)$(hover_ref_label(p1.name, seg.point_idxs[1]))"
                point1_label_pos[] = p1_2d + Point2f(20, 0)
                point1_label_visible[] = true
                point2_label[] = "$(prefix)$(hover_ref_label(p2.name, seg.point_idxs[2]))"
                point2_label_pos[] = p2_2d + Point2f(20, 0)
                point2_label_visible[] = true
            else
                segment_label_visible[] = false
                point1_label_visible[] = false
                point2_label_visible[] = false
            end
            last_hovered[] = hover
        end
    end

    # --- Camera movement: update label positions ---
    on(scene.camera.view, priority=1) do _
        sys_i, seg_i = last_hovered[]
        if sys_i != -1 && seg_i != -1
            sys = systems[sys_i]
            seg = sys.segments[seg_i]
            p1_3d = sys.points[seg.point_idxs[1]].pos_w
            p2_3d = sys.points[seg.point_idxs[2]].pos_w

            mid_point_2d = Makie.project(scene, (p1_3d + p2_3d) / 2)
            segment_label_pos[] = mid_point_2d + Point2f(20, 0)

            p1_2d = Makie.project(scene, p1_3d)
            p2_2d = Makie.project(scene, p2_3d)
            point1_label_pos[] = p1_2d + Point2f(20, 0)
            point2_label_pos[] = p2_2d + Point2f(20, 0)
        end
    end

    # --- Click-to-zoom handler ---
    on(events(scene).mousebutton, priority=1) do event
        if event.button == Mouse.left && event.action == Mouse.press
            cam = scene.camera
            sys_i, seg_i = last_hovered[]

            if sys_i != -1 && seg_i != -1
                sys = systems[sys_i]
                seg = sys.segments[seg_i]
                p1_w = sys.points[seg.point_idxs[1]].pos_w
                p2_w = sys.points[seg.point_idxs[2]].pos_w

                center = (p1_w + p2_w) / 2.0f0
                segment_len = norm(p2_w - p1_w)
                dist_heuristic = segment_len * 1.5 + 2.0

                if !PLOT_ZOOMED_IN[] && !PLOT_BODY_FRAME[]
                    cc = Makie.cameracontrols(scene)
                    PLOT_WORLD_EYEPOS[] = Vec3f(cc.eyeposition[])
                    PLOT_WORLD_LOOKAT[] = Vec3f(cc.lookat[])
                end

                inv_view_matrix = inv(cam.view[])
                cam_dir = normalize(Vec3f(inv_view_matrix[1,3], inv_view_matrix[2,3], inv_view_matrix[3,3]))
                new_eyepos = center + dist_heuristic * cam_dir

                update_cam!(scene, new_eyepos, center)
                PLOT_CAMERA_DISTANCE[] = dist_heuristic
                zoomed_in[] = true
                PLOT_ZOOMED_IN[] = true
                PLOT_ZOOM_SEGMENT_IDX[] = seg_i
                PLOT_ZOOM_BODY_IDX[] = -1
                PLOT_PREV_ZOOMED_IN[] = true
                PLOT_PREV_SEGMENT_IDX[] = seg_i

                sleep(0.01)
                p1_3d = sys.points[seg.point_idxs[1]].pos_w
                p2_3d = sys.points[seg.point_idxs[2]].pos_w
                mid_point_2d = Makie.project(scene, (p1_3d + p2_3d) / 2)
                segment_label_pos[] = mid_point_2d + Point2f(20, 0)
                p1_2d = Makie.project(scene, p1_3d)
                p2_2d = Makie.project(scene, p2_3d)
                point1_label_pos[] = p1_2d + Point2f(20, 0)
                point2_label_pos[] = p2_2d + Point2f(20, 0)

                return Consume(true)
            elseif zoomed_in[] && sys_i == -1
                zoomed_in[] = false
                PLOT_ZOOMED_IN[] = false
                PLOT_ZOOM_SEGMENT_IDX[] = -1
                if !isnothing(PLOT_SYSTEM_STRUCTURE[])
                    apply_zoom_mode!(scene, all_plots, PLOT_SYSTEM_STRUCTURE[]; mode_changed=true)
                else
                    zoom_out!(scene, cam, all_plots, nothing; relmargin)
                end

                return Consume(true)
            end
        end
        return Consume(false)
    end

    return nothing
end

"""
    setup_body_zoom_events!(scene, sys; relmargin=0.2,
                            spoke_colors_obs=nothing, highlight_color=:red)

Hover labels and click-to-zoom for standalone rigid bodies, picking on the grey
joint spokes (origin→anchor segments). Hovering a body shows its name and, when
`spoke_colors_obs` (the per-spoke color observable of `plots[:body_chain]`) is
given, recolors that body's spokes to `highlight_color`. Clicking zooms in and
tracks it (`PLOT_ZOOM_BODY_IDX`, followed each frame by `apply_zoom_mode!`);
clicking empty space restores the prior view. No-op when the system has no rigid
bodies. Mirrors `setup_segment_hover_events!`; the two don't clash because
segment-less models (e.g. a beam) never install the segment one.
"""
function setup_body_zoom_events!(scene, sys; relmargin=0.2,
                                 spoke_colors_obs=nothing, highlight_color=:red)
    isempty(sys.bodies) && return nothing
    last_body = Ref(-1)
    base_spoke_color = to_color(RGBf(0.32, 0.32, 0.32))
    spoke_owner = spoke_body_indices(sys)

    body_label = Observable("")
    body_label_pos = Observable(Point2f(0, 0))
    body_label_visible = Observable(false)
    text!(scene, body_label, position=body_label_pos, space=:pixel,
          fontsize=14, color=:white, strokecolor=:black, strokewidth=1,
          align=(:center, :center), visible=body_label_visible, transparency=true)

    # Min 2D distance from the cursor to any of body `b`'s joint spokes.
    function spoke_distance(b, mouse_2d)
        body = sys.bodies[b]
        R = SymbolicAWEModels.quaternion_to_rotation_matrix(body.Q_b_to_w)
        origin_2d = Makie.project(scene, Point3f(body.pos_w))
        dist = Inf
        for joint in sys.elastic_joints
            anchor = if joint.body_a_idx == b
                joint.anchor_a_b
            elseif joint.body_b_idx == b
                joint.anchor_b_b
            else
                continue
            end
            anchor_2d = Makie.project(scene, Point3f(body.pos_w .+ R * anchor))
            dist = min(dist, point_line_segment_distance(mouse_2d, origin_2d, anchor_2d))
        end
        return dist
    end

    on(events(scene).mouseposition, priority=2) do mp
        mouse_2d = Point2f(mp)
        margin_px = 30.0
        best, best_dist = -1, Inf
        for b in eachindex(sys.bodies)
            d = spoke_distance(b, mouse_2d)
            d < best_dist && ((best, best_dist) = (b, d))
        end
        hovered = best_dist < margin_px ? best : -1
        if hovered != last_body[]
            if !isnothing(spoke_colors_obs)
                spoke_colors_obs[] = [owner == hovered ? to_color(highlight_color) :
                                      base_spoke_color for owner in spoke_owner]
            end
            if hovered != -1
                body = sys.bodies[hovered]
                body_label[] = hover_ref_label(body.name, hovered)
                body_label_pos[] =
                    Makie.project(scene, Point3f(body.pos_w)) + Point2f(20, 0)
                body_label_visible[] = true
            else
                body_label_visible[] = false
            end
            last_body[] = hovered
        end
    end

    on(scene.camera.view, priority=1) do _
        if last_body[] != -1
            body = sys.bodies[last_body[]]
            body_label_pos[] =
                Makie.project(scene, Point3f(body.pos_w)) + Point2f(20, 0)
        end
    end

    on(events(scene).mousebutton, priority=1) do event
        if event.button == Mouse.left && event.action == Mouse.press
            cam = scene.camera
            if last_body[] != -1
                if !PLOT_ZOOMED_IN[] && !PLOT_BODY_FRAME[]
                    cc = Makie.cameracontrols(scene)
                    PLOT_WORLD_EYEPOS[] = Vec3f(cc.eyeposition[])
                    PLOT_WORLD_LOOKAT[] = Vec3f(cc.lookat[])
                end
                dist = zoom_in_body!(scene, cam, sys, last_body[])
                PLOT_CAMERA_DISTANCE[] = dist
                PLOT_ZOOMED_IN[] = true
                PLOT_PREV_ZOOMED_IN[] = true
                PLOT_ZOOM_BODY_IDX[] = last_body[]
                PLOT_ZOOM_SEGMENT_IDX[] = -1
                return Consume(true)
            elseif PLOT_ZOOMED_IN[]
                PLOT_ZOOMED_IN[] = false
                PLOT_ZOOM_BODY_IDX[] = -1
                PLOT_ZOOM_SEGMENT_IDX[] = -1
                if !isnothing(PLOT_WORLD_EYEPOS[]) && !isnothing(PLOT_WORLD_LOOKAT[])
                    update_cam!(scene, PLOT_WORLD_EYEPOS[], PLOT_WORLD_LOOKAT[])
                end
                return Consume(true)
            end
        end
        return Consume(false)
    end
    return nothing
end

function plot_with_panes(sys::SystemStructure;
                    size = (1200, 800),
                    relmargin = 0.2,
                    segment_color = :black,
                    highlight_color = :red,
                    force_color = false,
                    body_frame = false,
                    extra_points = nothing,
                    extra_groups = nothing,
                    mesh = nothing,
                    kwargs...)
    # Use LScene for advanced camera controls
    scene = Scene(; camera=cam3d!, show_axis = false, size, zoommode = :free, samples = 16)
    plots = plot!(scene, sys; segment_color, force_color,
                  extra_points, extra_groups, mesh, kwargs...)
    
    relevant_plots = AbstractPlot[]
    if haskey(plots, :segments)
        push!(relevant_plots, plots[:segments])
    end
    if haskey(plots, :points)
        push!(relevant_plots, plots[:points])
    end
    if haskey(plots, :wings)
        # plots[:wings] can be either an array (static) or a single plot (observables)
        if plots[:wings] isa AbstractArray
            append!(relevant_plots, plots[:wings])
        else
            push!(relevant_plots, plots[:wings])
        end
    end
    # Standalone rigid bodies (so the camera frames a body-only scene, e.g. a beam)
    haskey(plots, :body_chain) && push!(relevant_plots, plots[:body_chain])
    haskey(plots, :bodies) && push!(relevant_plots, plots[:bodies])

    # --- Event Handling for Segments (using extracted reusable function) ---
    if haskey(plots, :segments)
        setup_segment_hover_events!(scene, [sys], [plots[:segment_colors_obs]], relevant_plots;
                                    segment_colors=[segment_color],
                                    highlight_color, force_color, relmargin)
    end
    # Hover labels + highlight + click-to-zoom for standalone rigid bodies.
    setup_body_zoom_events!(scene, sys; relmargin, highlight_color,
                            spoke_colors_obs=get(plots, :body_chain_colors_obs, nothing))

    # Set initial camera position
    cam = scene.camera
    if body_frame
        initial_distance = zoom_body_frame!(scene, cam, sys)
    else
        update_cam!(scene, Vec3f(-100, -100, 100), Vec3f(0, 0, 0))
        initial_distance = zoom_out!(scene, cam, relevant_plots, nothing; relmargin)
    end

    return scene, plots, relevant_plots, initial_distance
end

# Public API function - creates scene with observables for dynamic updates
function Makie.plot(sys::SystemStructure;
                    vector_scale=1.0,
                    force_color=false,
                    segment_color=:black,
                    plot_aero=true,
                    relmargin=0.2,
                    body_frame=false,
                    extra_points=nothing,
                    extra_groups=nothing,
                    mesh=nothing,
                    kwargs...)
    # Store SystemStructure globally FIRST so @lift expressions can access it
    PLOT_SYSTEM_STRUCTURE[] = sys
    PLOT_VECTOR_SCALE[] = vector_scale
    PLOT_FORCE_COLOR[] = force_color
    PLOT_SEGMENT_COLOR[] = segment_color
    PLOT_BODY_FRAME[] = body_frame

    # Create single geometry trigger observable
    geometry_obs = Observable(0.0)

    # Create scene with observables using internal function
    scene, plots, relevant_plots, initial_distance = plot_with_panes(sys;
                geometry_obs,
                vector_scale,
                force_color,
                segment_color,
                plot_aero,
                relmargin,
                body_frame,
                extra_points,
                extra_groups,
                mesh,
                kwargs...)

    # Get the segment colors observable from the plots if available
    segment_colors_obs = nothing
    if haskey(plots, :segment_colors_obs)
        segment_colors_obs = plots[:segment_colors_obs]
    end

    # Store observables and settings globally
    PLOT_GEOMETRY_OBS[] = geometry_obs
    PLOT_SEGMENT_COLORS_OBS[] = segment_colors_obs
    PLOT_SCENE[] = scene
    PLOT_RELEVANT_PLOTS[] = relevant_plots
    # SystemStructure and settings already stored above before plot creation
    PLOT_ZOOMED_IN[] = false  # Initialize zoom state (not zoomed in)
    PLOT_ZOOM_RELMARGIN[] = relmargin  # Store relmargin for auto-updates
    PLOT_ZOOM_SEGMENT_IDX[] = -1  # No segment zoomed initially
    PLOT_CAMERA_DISTANCE[] = initial_distance  # Preserve initial zoom during replay

    return scene
end

"""
    build_geometry_observables(sys::SystemStructure, trigger::Observable)

Build geometry observables for a SystemStructure that update when trigger is notified.
Used for multi-system plotting where each system needs its own observable.

# Arguments
- `sys::SystemStructure`: The system structure to build observables for
- `trigger::Observable`: Observable that triggers geometry recalculation

# Returns
- Dict with keys :segment_points, :point_positions, :wing_data
"""
function build_geometry_observables(sys::SystemStructure, trigger::Observable)
    obs = Dict{Symbol, Any}()

    # Segments: Vector of Point3f pairs
    obs[:segment_points] = @lift begin
        $trigger
        points = Point3f[]
        for seg in sys.segments
            push!(points, Point3f(sys.points[seg.point_idxs[1]].pos_w))
            push!(points, Point3f(sys.points[seg.point_idxs[2]].pos_w))
        end
        points
    end

    # Points: Vector of Point3f
    obs[:point_positions] = @lift begin
        $trigger
        [Point3f(p.pos_w) for p in sys.points]
    end

    # Wing axes: origins and directions
    obs[:wing_data] = @lift begin
        $trigger
        origins = Point3f[]
        directions = Vec3f[]
        for wing in sys.wings
            origin = Point3f(wing.pos_w)
            R = wing.R_b_to_w
            for i in 1:3
                push!(origins, origin)
                push!(directions, Vec3f(R[:, i]))
            end
        end
        (origins, directions)
    end

    return obs
end

"""
    Makie.plot(syss::Vector{<:SystemStructure}; colors=Makie.wong_colors(), ...)

Plot multiple SystemStructures in a single scene with different colors.

# Arguments
- `syss::Vector{<:SystemStructure}`: Vector of system structures to plot

# Keyword Arguments
- `colors`: Color palette for systems (default: wong_colors())
- `vector_scale::Real=1.0`: Scale factor for wing orientation arrows
- `relmargin::Real=0.2`: Relative margin for camera zoom
- `use_observables::Bool=false`: If true, create observables for dynamic updates (used by replay)
- All other kwargs passed to plot!

# Returns
- Scene with all systems plotted

# Example
```julia
scene = plot([sys1, sys2])  # Plot two systems with different colors
```
"""
function Makie.plot(syss::Vector{<:SystemStructure};
                    colors=Makie.wong_colors(),
                    vector_scale=1.0,
                    relmargin=0.2,
                    use_observables=false,
                    highlight_color=:yellow,
                    force_color=false,
                    kwargs...)
    scene = Scene(; camera=cam3d!, show_axis=false, size=(1200, 800))

    PLOT_MULTI_SYSTEMS[] = syss
    all_plots = AbstractPlot[]
    segment_color_obs = Observable[]  # Per-segment color obs per system, for hover
    segment_colors_list = Any[]       # Track base colors for hover reset

    if use_observables
        # Create trigger observables for each system
        triggers = [Observable(0.0) for _ in syss]
        PLOT_MULTI_GEOMETRY_OBS[] = triggers

        for (i, sys) in enumerate(syss)
            color = colors[mod1(i, length(colors))]
            push!(segment_colors_list, color)
            obs = build_geometry_observables(sys, triggers[i])

            # Observable per-segment colors so the hover handler can highlight
            # (setup_segment_hover_events! writes Vector{RGBA} into this obs).
            seg_colors = Observable(fill(to_color(color), length(sys.segments)))
            seg_plot = linesegments!(scene, obs[:segment_points];
                                     color=seg_colors, linewidth=2,
                                     transparency=true)
            push!(all_plots, seg_plot)
            push!(segment_color_obs, seg_colors)

            pt_plot = scatter!(scene, obs[:point_positions];
                              color=color, transparency=true)
            push!(all_plots, pt_plot)

            # Wing arrows
            if !isempty(sys.wings)
                arrows3d!(scene,
                          @lift($(obs[:wing_data])[1]),
                          @lift($(obs[:wing_data])[2] .* $vector_scale);
                          color=repeat([:red, :green, :blue], length(sys.wings)))
            end
        end
    else
        # Static mode: reuse existing plot!
        for (i, sys) in enumerate(syss)
            color = colors[mod1(i, length(colors))]
            push!(segment_colors_list, color)
            plots = plot!(scene, sys;
                          segment_color=color,
                          point_color=color,
                          vector_scale,
                          geometry_obs=nothing,
                          kwargs...)
            if haskey(plots, :segments)
                push!(all_plots, plots[:segments])
            end
            # One obs per system (aligned with `systems`); empty for segment-less.
            push!(segment_color_obs, get(plots, :segment_colors_obs,
                                         Observable(RGBAf[])))
            haskey(plots, :points) && push!(all_plots, plots[:points])
        end
    end

    # Add hover/click events for segments
    if !isempty(segment_color_obs)
        setup_segment_hover_events!(scene, syss, segment_color_obs, all_plots;
                                    segment_colors=segment_colors_list,
                                    highlight_color, force_color, relmargin)
    end

    zoom_out!(scene, scene.camera, all_plots; relmargin)
    return scene
end

"""
    Makie.plot!(sys::SystemStructure; vector_scale=1.0)

Update the currently displayed SystemStructure plot with new data from `sys`.

This function follows standard Makie conventions: `plot!` with `!` mutates the existing
scene by updating its observables. Must be called after an initial `plot(sys)` has created
the scene and observables.

# Arguments
- `sys::SystemStructure`: The system structure with updated state to display

# Keyword Arguments
- `vector_scale::Real=1.0`: Scale factor for wing orientation arrows

# Returns
- `nothing` (mutates existing scene via observables)

# Example
```julia
# Create initial plot
scene = plot(sys_struct)

# In simulation loop, update the plot
for i in 1:100
    next_step!(sam)
    plot!(sys_struct)  # Updates observables
    sleep(0.01)
end
```

# See Also
- `plot(::SystemStructure)`: Create a new scene (non-mutating)
- `update_plot_observables!`: Lower-level observable update function
"""
function Makie.plot!(sys::SystemStructure; vector_scale=1.0)
    # Check if geometry observable exists
    if isnothing(PLOT_GEOMETRY_OBS[]) || isnothing(PLOT_SCENE[])
        error("No existing plot to update. Call plot(sys) first to create a scene.")
    end

    # Check if the scene is still active
    scene = PLOT_SCENE[]
    scene_has_display = false
    try
        if hasfield(typeof(scene), :events) &&
           hasfield(typeof(scene.events), :window_open)
            scene_has_display = scene.events.window_open[]
        else
            scene_has_display = true
        end
    catch
        scene_has_display = false
    end

    if !scene_has_display
        error("Plot window has been closed. Call plot(sys) to create a new scene.")
    end

    # Update vector scale if changed
    if vector_scale != PLOT_VECTOR_SCALE[]
        PLOT_VECTOR_SCALE[] = vector_scale
    end

    # Trigger observable update
    update_plot_observables!(sys)

    return nothing
end

"""
    record(lg::SysLog, sys::SystemStructure, filename::String; framerate=30, kwargs...)

Record a SysLog animation to a video or GIF file. The output format is
determined by the file extension (`.mp4`, `.gif`, `.mkv`, `.webm`).

# Arguments
- `lg::SysLog`: The simulation log to record
- `sys::SystemStructure`: The system structure matching the log's topology
- `filename::String`: Output filename (e.g., `"sim.mp4"`, `"sim.gif"`)

# Keyword Arguments
- `framerate::Int=30`: Framerate (frames per second)
- `vector_scale::Real=1.0`: Scale factor for wing orientation arrows
- All other keyword arguments are passed through to the SystemStructure plot function

# Returns
- The Scene object used for recording

# Example
```julia
record(log, sys_struct, "output.mp4")
record(log, sys_struct, "output.gif"; framerate=20)
record(log, sys_struct, "sim.mp4"; framerate=60, vector_scale=0.3)
```
"""
function SymbolicAWEModels.record(lg::SysLog, sys::SystemStructure, filename::String;
                                   framerate::Int=30,
                                   vector_scale::Real=0.2,
                                   kwargs...)
    n_frames = length(lg.syslog)
    n_frames == 0 && error("Empty SysLog provided for recording")

    println("Recording video to: $filename")
    println("Framerate: $framerate fps")
    println("Total frames: $n_frames")

    # Initialize with first state
    update_from_sysstate!(sys, lg.syslog[1])

    # Create initial plot with observables (following standard Makie pattern)
    scene = plot(sys; vector_scale, kwargs...)

    # Record video by stepping through all frames
    Makie.record(scene, filename, 1:n_frames; framerate=framerate) do frame_num
        # Update system state
        update_from_sysstate!(sys, lg.syslog[frame_num])

        # Update plot observables (mutating, following standard Makie pattern)
        plot!(sys; vector_scale)

        # Critical: Yield to GLMakie's render thread to process observable updates
        # Without this, the record function captures the same (first) frame repeatedly
        sleep(0.001)

        # Print progress every 10% or at the end
        if frame_num % max(1, div(n_frames, 10)) == 0 || frame_num == n_frames
            progress_pct = round(100 * frame_num / n_frames, digits=1)
            println("  Progress: $progress_pct% ($frame_num/$n_frames frames)")
        end
    end

    println("Video saved successfully!")
    return scene
end

"""
    record(logs::Vector{<:SysLog}, syss::Vector{<:SystemStructure},
           filename::String; framerate=30, colors=Makie.wong_colors(),
           vector_scale=0.2, kwargs...)

Record a multi-system SysLog animation to a video or GIF file.
The output format is determined by the file extension
(`.mp4`, `.gif`, `.mkv`, `.webm`).

# Arguments
- `logs::Vector{<:SysLog}`: Simulation logs (one per system)
- `syss::Vector{<:SystemStructure}`: System structures
    (must match logs length)
- `filename::String`: Output filename
    (e.g., `"sim.mp4"`, `"sim.gif"`)

# Keyword Arguments
- `framerate::Int=30`: Framerate (frames per second)
- `colors`: Color palette for distinguishing systems
    (default: wong_colors())
- `vector_scale::Real=0.2`: Scale factor for wing arrows
- All other keyword arguments are passed through to `plot`

# Returns
- The Scene object used for recording

# Example
```julia
record([log1, log2], [sys1, sys2], "output.mp4")
record([log1, log2], [sys1, sys2], "output.mp4";
       framerate=60, colors=[:red, :blue])
```
"""
function SymbolicAWEModels.record(
        logs::Vector{<:SysLog},
        syss::Vector{<:SystemStructure},
        filename::String;
        framerate::Int=30,
        colors=Makie.wong_colors(),
        vector_scale::Real=0.2,
        kwargs...)
    length(logs) == length(syss) ||
        error("logs and systems must have same length")
    n_frames = minimum(length(lg.syslog) for lg in logs)
    n_frames == 0 &&
        error("Empty SysLog provided for recording")

    println("Recording multi-system video to: $filename")
    println("Framerate: $framerate fps")
    println("Total frames: $n_frames")
    println("Systems: $(length(syss))")

    # Initialize all systems with first state
    for (sys, lg) in zip(syss, logs)
        update_from_sysstate!(sys, lg.syslog[1])
    end

    # Create plot with observables for dynamic updates
    scene = plot(syss; colors, vector_scale,
                 use_observables=true, kwargs...)

    # Record video by stepping through all frames
    Makie.record(scene, filename, 1:n_frames;
                 framerate=framerate) do frame_num
        # Update each system's state
        for (sys, lg) in zip(syss, logs)
            frame = min(frame_num, length(lg.syslog))
            update_from_sysstate!(sys, lg.syslog[frame])
        end
        # Trigger all geometry observables
        if !isnothing(PLOT_MULTI_GEOMETRY_OBS[])
            for obs in PLOT_MULTI_GEOMETRY_OBS[]
                obs[] = time()
            end
        end

        # Yield to render thread for observable updates
        sleep(0.001)

        # Print progress every 10% or at the end
        if frame_num % max(1, div(n_frames, 10)) == 0 ||
                frame_num == n_frames
            pct = round(100 * frame_num / n_frames, digits=1)
            println("  Progress: $pct% " *
                    "($frame_num/$n_frames frames)")
        end
    end

    println("Video saved successfully!")
    return scene
end

"""
    setup_replay_controls!(scene, n_frames, update_frame!, get_time, get_dt;
                           replay_speed=1.0, autoplay=false, loop=false)

Create replay UI controls (slider, buttons, animation loop) for a scene.
This is a reusable function that can be used by both single and multi-system replay.

# Arguments
- `scene`: Parent scene to overlay UI on
- `n_frames`: Total number of frames
- `update_frame!(idx)`: Callback to update geometry for frame idx
- `get_time(idx)`: Callback to get simulation time for frame idx
- `get_dt(idx)`: Callback to get time delta between frames (for playback timing)

# Keyword Arguments
- `replay_speed::Real=1.0`: Replay speed factor
- `autoplay::Bool=false`: Start playing automatically
- `loop::Bool=false`: Loop playback continuously

# Returns
- `(ui_scene, frame_idx, is_playing)`: UI scene and observables for external control
"""
function setup_replay_controls!(scene, n_frames, update_frame!, get_time, get_dt;
                                 replay_speed=1.0, autoplay=false, loop=false)
    # Create pixel-space subscene for UI controls overlay
    ui_scene = Scene(scene, viewport=scene.viewport, clear=false, camera=campixel!)

    # Create observable for current frame index
    frame_idx = Observable(1)

    # UI Layout constants (pixel coordinates from bottom-left)
    ui_height = 120
    ui_margin = 20
    button_height = 30
    button_width = 80
    slider_height = 20

    # Get scene size
    scene_width = Observable(scene.viewport[].widths[1])
    scene_height = Observable(scene.viewport[].widths[2])

    # Update scene size on resize
    on(scene.viewport) do area
        scene_width[] = area.widths[1]
        scene_height[] = area.widths[2]
    end

    # --- Slider Implementation ---
    slider_y = ui_margin + button_height + 10
    slider_x_start = ui_margin
    slider_width_obs = @lift($(scene_width) - 2 * ui_margin)

    # Slider track
    slider_track_rect = @lift(Rect2f(slider_x_start, slider_y, $(slider_width_obs), slider_height))
    poly!(ui_scene, slider_track_rect, color=RGBAf(0.3, 0.3, 0.3, 0.8))

    # Slider thumb position
    slider_thumb_x = @lift(slider_x_start + ($(frame_idx) - 1) / max(1, n_frames - 1) * $(slider_width_obs))
    slider_thumb_pos = @lift(Point2f($(slider_thumb_x), slider_y + slider_height / 2))
    scatter!(ui_scene, slider_thumb_pos, color=:white, markersize=15)

    # --- Button Implementation ---
    button_y = ui_margin

    # Play/Pause button
    is_playing = Observable(autoplay)
    play_button_x = ui_margin
    play_button_rect = Rect2f(play_button_x, button_y, button_width, button_height)
    play_button_color = @lift($(is_playing) ? RGBAf(0.8, 0.3, 0.3, 0.8) : RGBAf(0.3, 0.8, 0.3, 0.8))
    poly!(ui_scene, play_button_rect, color=play_button_color)
    play_button_label = @lift($(is_playing) ? "Pause" : "Play")
    text!(ui_scene, play_button_label, position=Point2f(play_button_x + button_width/2, button_y + button_height/2),
          align=(:center, :center), fontsize=14, color=:white)

    # Step backward button
    step_back_x = play_button_x + button_width + 10
    step_back_rect = Rect2f(step_back_x, button_y, button_width, button_height)
    poly!(ui_scene, step_back_rect, color=RGBAf(0.4, 0.4, 0.4, 0.8))
    text!(ui_scene, "<", position=Point2f(step_back_x + button_width/2, button_y + button_height/2),
          align=(:center, :center), fontsize=14, color=:white)

    # Step forward button
    step_forward_x = step_back_x + button_width + 10
    step_forward_rect = Rect2f(step_forward_x, button_y, button_width, button_height)
    poly!(ui_scene, step_forward_rect, color=RGBAf(0.4, 0.4, 0.4, 0.8))
    text!(ui_scene, ">", position=Point2f(step_forward_x + button_width/2, button_y + button_height/2),
          align=(:center, :center), fontsize=14, color=:white)

    # Body frame toggle button
    body_frame_button_x = step_forward_x + button_width + 10
    body_frame_button_rect = Rect2f(body_frame_button_x, button_y, button_width, button_height)
    body_frame_obs = Observable(PLOT_BODY_FRAME[])
    body_frame_button_color = @lift($(body_frame_obs) ? RGBAf(0.3, 0.6, 0.8, 0.8) : RGBAf(0.5, 0.5, 0.5, 0.8))
    poly!(ui_scene, body_frame_button_rect, color=body_frame_button_color)
    body_frame_button_label = @lift($(body_frame_obs) ? "Body" : "World")
    text!(ui_scene, body_frame_button_label, position=Point2f(body_frame_button_x + button_width/2, button_y + button_height/2),
          align=(:center, :center), fontsize=14, color=:white)

    # Save button
    save_button_x = body_frame_button_x + button_width + 10
    save_button_rect = Rect2f(save_button_x, button_y, button_width, button_height)
    poly!(ui_scene, save_button_rect, color=RGBAf(0.2, 0.7, 0.3, 0.8))
    text!(ui_scene, "Save", position=Point2f(save_button_x + button_width/2, button_y + button_height/2),
          align=(:center, :center), fontsize=14, color=:white)

    # Info label
    info_label_x = save_button_x + button_width + 20
    info_text = @lift("Frame: $($(frame_idx))/$n_frames | Time: $(@sprintf("%.2f", get_time($(frame_idx)))) s | Speed: $(replay_speed)x")
    text!(ui_scene, info_text, position=Point2f(info_label_x, button_y + button_height/2),
          align=(:left, :center), fontsize=14, color=:white, strokecolor=:black, strokewidth=1)

    # Combined mouse event handling for slider and buttons
    slider_dragging = Ref(false)

    on(events(ui_scene).mousebutton, priority = 2) do event
        if event.button == Mouse.left
            mp = events(ui_scene).mouseposition[]

            if event.action == Mouse.press
                # Check if click is on slider
                track_rect = slider_track_rect[]
                if mp[2] >= track_rect.origin[2] && mp[2] <= track_rect.origin[2] + track_rect.widths[2] &&
                   mp[1] >= track_rect.origin[1] && mp[1] <= track_rect.origin[1] + track_rect.widths[1]
                    slider_dragging[] = true
                    # Update frame immediately
                    rel_pos = clamp((mp[1] - slider_x_start) / slider_width_obs[], 0.0, 1.0)
                    new_idx = round(Int, 1 + rel_pos * (n_frames - 1))
                    frame_idx[] = clamp(new_idx, 1, n_frames)
                    update_frame!(frame_idx[])
                    return Consume(true)
                end

                # Check play button
                if mp[1] >= play_button_rect.origin[1] && mp[1] <= play_button_rect.origin[1] + play_button_rect.widths[1] &&
                   mp[2] >= play_button_rect.origin[2] && mp[2] <= play_button_rect.origin[2] + play_button_rect.widths[2]
                    is_playing[] = !is_playing[]
                    return Consume(true)
                end

                # Check step back button
                if mp[1] >= step_back_rect.origin[1] && mp[1] <= step_back_rect.origin[1] + step_back_rect.widths[1] &&
                   mp[2] >= step_back_rect.origin[2] && mp[2] <= step_back_rect.origin[2] + step_back_rect.widths[2]
                    new_idx = max(1, frame_idx[] - 1)
                    frame_idx[] = new_idx
                    update_frame!(new_idx)
                    return Consume(true)
                end

                # Check step forward button
                if mp[1] >= step_forward_rect.origin[1] && mp[1] <= step_forward_rect.origin[1] + step_forward_rect.widths[1] &&
                   mp[2] >= step_forward_rect.origin[2] && mp[2] <= step_forward_rect.origin[2] + step_forward_rect.widths[2]
                    new_idx = min(n_frames, frame_idx[] + 1)
                    frame_idx[] = new_idx
                    update_frame!(new_idx)
                    return Consume(true)
                end

                # Check body frame toggle button
                if mp[1] >= body_frame_button_rect.origin[1] && mp[1] <= body_frame_button_rect.origin[1] + body_frame_button_rect.widths[1] &&
                   mp[2] >= body_frame_button_rect.origin[2] && mp[2] <= body_frame_button_rect.origin[2] + body_frame_button_rect.widths[2]
                    PLOT_BODY_FRAME[] = !PLOT_BODY_FRAME[]
                    body_frame_obs[] = PLOT_BODY_FRAME[]
                    if PLOT_BODY_FRAME[]
                        PLOT_ZOOMED_IN[] = false
                        PLOT_ZOOM_SEGMENT_IDX[] = -1
                    end
                    if !isnothing(PLOT_SCENE[]) && !isnothing(PLOT_RELEVANT_PLOTS[]) && !isnothing(PLOT_SYSTEM_STRUCTURE[])
                        apply_zoom_mode!(PLOT_SCENE[], PLOT_RELEVANT_PLOTS[],
                                         PLOT_SYSTEM_STRUCTURE[]; mode_changed=true)
                    end
                    return Consume(true)
                end

                # Check save button
                if mp[1] >= save_button_rect.origin[1] && mp[1] <= save_button_rect.origin[1] + save_button_rect.widths[1] &&
                   mp[2] >= save_button_rect.origin[2] && mp[2] <= save_button_rect.origin[2] + save_button_rect.widths[2]
                    # Generate filename with timestamp and simulation time
                    sim_time = get_time(frame_idx[])
                    timestamp = Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS")
                    filename = "replay_$(timestamp)_t$(round(sim_time, digits=2))s_frame$(frame_idx[]).png"
                    # Hide UI scene, save, then restore
                    ui_scene.visible[] = false
                    save(filename, scene; px_per_unit=2)
                    ui_scene.visible[] = true
                    @info "Saved high-res screenshot" filename
                    return Consume(true)
                end
            elseif event.action == Mouse.release
                if slider_dragging[]
                    slider_dragging[] = false
                    return Consume(true)
                end
            end
        end
        return Consume(false)
    end

    on(events(ui_scene).mouseposition, priority=2) do mp
        if slider_dragging[]
            rel_pos = clamp((mp[1] - slider_x_start) / slider_width_obs[], 0.0, 1.0)
            new_idx = round(Int, 1 + rel_pos * (n_frames - 1))
            frame_idx[] = clamp(new_idx, 1, n_frames)
            update_frame!(frame_idx[])
        end
    end

    # Animation loop with replay speed
    @async begin
        while true
            if is_playing[]
                if frame_idx[] < n_frames
                    new_idx = frame_idx[] + 1
                    frame_idx[] = new_idx
                    update_frame!(new_idx)
                    # Calculate sleep time based on actual time difference and replay speed
                    dt = get_dt(frame_idx[])
                    sleep(max(0.01, dt / replay_speed))
                elseif loop
                    frame_idx[] = 1
                    update_frame!(1)
                else
                    is_playing[] = false  # Stop at end
                end
            end
            sleep(0.02)  # Check state frequently
        end
    end

    return ui_scene, frame_idx, is_playing
end

"""
    replay(lg::SysLog, sys::SystemStructure; replay_speed=1.0, autoplay=false, loop=false, kwargs...)

Replay a SysLog with interactive 3D visualization and playback controls.

# Arguments
- `lg::SysLog`: The simulation log to replay
- `sys::SystemStructure`: The system structure matching the log's topology

# Keyword Arguments
# - `replay_speed::Real=1.0`: Replay speed factor (1.0 = real-time, 2.0 = 2x speed, etc.)
# - `autoplay::Bool=false`: Start playing automatically when opened
# - `loop::Bool=false`: Loop playback continuously
# - `vector_scale::Real=1.0`: Scale factor for wing orientation arrows
# - `show_panes::Bool=true`: Show gray background reference panes (set `false` for white-only background)
# - All other keyword arguments are passed through to the SystemStructure plot function

# Returns
- A Scene with interactive controls overlaid on the 3D visualization

# UI Controls
- **Slider**: Drag to scrub through time
- **Play/Pause button**: Toggle playback (green when stopped, red when playing)
- **< button**: Step backward one frame
- **> button**: Step forward one frame
- **Info display**: Shows current frame, time, and playback speed

# Example
```julia
# Create interactive replay viewer
scene = replay(log, sys_struct)

# Auto-play at 2x speed with looping
scene = replay(log, sys_struct, replay_speed=2.0, autoplay=true, loop=true)

# Replay with custom visualization settings
scene = replay(log, sys_struct, replay_speed=0.5, vector_scale=0.3)
```

# See Also
- `record`: For saving replay to MP4 video file
"""
function SymbolicAWEModels.replay(lg::SysLog, sys::SystemStructure;
                      replay_speed=1.0,
                      autoplay=false,
                      loop=false,
                      vector_scale=1.0,
                      show_panes::Bool=true,
                      kwargs...)

    n_frames = length(lg.syslog)
    n_frames == 0 && error("Empty SysLog provided for replay")

    # Initialize with first state
    update_from_sysstate!(sys, lg.syslog[1])

    # Create initial plot which sets up observables and scene
    scene = plot(sys; vector_scale, kwargs...)

    # Define callbacks for UI controls
    function update_frame!(idx)
        update_from_sysstate!(sys, lg.syslog[idx])
        plot!(sys; vector_scale)
    end

    get_time(idx) = lg.syslog[idx].time
    get_dt(idx) = idx > 1 ? lg.syslog[idx].time - lg.syslog[idx - 1].time : 0.05

    # Setup replay controls using shared function
    setup_replay_controls!(scene, n_frames, update_frame!, get_time, get_dt;
                           replay_speed, autoplay, loop)

    return scene
end

"""
    replay(logs::Vector{<:SysLog}, syss::Vector{<:SystemStructure}; colors=Makie.wong_colors(), ...)

Replay multiple SysLogs with their corresponding SystemStructures simultaneously.

# Arguments
- `logs::Vector{<:SysLog}`: Vector of simulation logs to replay
- `syss::Vector{<:SystemStructure}`: Vector of system structures (must match logs length)

# Keyword Arguments
- `colors`: Color palette for distinguishing systems (default: wong_colors())
- `replay_speed::Real=1.0`: Replay speed factor
- `autoplay::Bool=false`: Start playing automatically
- `loop::Bool=false`: Loop playback continuously
- `vector_scale::Real=1.0`: Scale factor for wing orientation arrows

# Returns
- A Scene with interactive controls overlaid on the 3D visualization

# Example
```julia
# Replay two simulations side by side
scene = replay([log1, log2], [sys1, sys2])

# With custom colors
scene = replay([log1, log2], [sys1, sys2], colors=[:red, :blue])
```
"""
function SymbolicAWEModels.replay(logs::Vector{<:SysLog}, syss::Vector{<:SystemStructure};
                      colors=Makie.wong_colors(),
                      replay_speed=1.0,
                      autoplay=false,
                      loop=false,
                      vector_scale=1.0,
                      kwargs...)

    length(logs) == length(syss) || error("logs and systems must have same length")
    n_frames = minimum(length(lg.syslog) for lg in logs)
    n_frames == 0 && error("Empty SysLog provided for replay")

    # Initialize all systems with first state
    for (sys, lg) in zip(syss, logs)
        update_from_sysstate!(sys, lg.syslog[1])
    end

    # Create plot with observables enabled for dynamic updates
    scene = plot(syss; colors, vector_scale, use_observables=true, kwargs...)

    # Define callbacks for UI controls
    function update_frame!(idx)
        # Update each system's state
        for (sys, lg) in zip(syss, logs)
            frame = min(idx, length(lg.syslog))
            update_from_sysstate!(sys, lg.syslog[frame])
        end
        # Trigger all geometry observables
        if !isnothing(PLOT_MULTI_GEOMETRY_OBS[])
            for obs in PLOT_MULTI_GEOMETRY_OBS[]
                obs[] = time()
            end
        end
    end

    # Use first log for timing
    get_time(idx) = logs[1].syslog[min(idx, length(logs[1].syslog))].time
    function get_dt(idx)
        if idx > 1
            lg = logs[1]
            i1 = min(idx, length(lg.syslog))
            i0 = min(idx - 1, length(lg.syslog))
            return lg.syslog[i1].time - lg.syslog[i0].time
        end
        return 0.05
    end

    # Setup replay controls using shared function
    setup_replay_controls!(scene, n_frames, update_frame!, get_time, get_dt;
                           replay_speed, autoplay, loop)

    return scene
end

"""
    plot_sphere_trajectory(logs; radius=1.0, colors=nothing, labels=nothing, ...)

Plot elevation-azimuth trajectories as lines on a sphere surface.

# Arguments
- `logs`: Vector of SysLog objects (or single SysLog)
- `radius`: Sphere radius (default 1.0)
- `colors`: Optional vector of colors for each trajectory
- `labels`: Optional vector of labels for legend
- `sphere_alpha`: Transparency of reference sphere (default 0.2)
- `linewidth`: Line width for trajectories (default 2.0)
- `size`: Figure size tuple (default (800, 800))
"""
function SymbolicAWEModels.plot_sphere_trajectory(lg::SysLog; kwargs...)
    SymbolicAWEModels.plot_sphere_trajectory([lg]; kwargs...)
end

function SymbolicAWEModels.plot_sphere_trajectory(logs::Vector{<:SysLog};
    radius=1.0,
    colors=nothing,
    labels=nothing,
    sphere_alpha=0.2,
    linewidth=2.0,
    size=(800, 800)
)
    fig = Figure(; size)
    ax = LScene(fig[1, 1], show_axis=false)
    # Disable perspective (orthographic projection)
    ax.scene.camera.projection[] = Makie.orthographicprojection(
        -2f0, 2f0, -2f0, 2f0, -10f0, 10f0)

    # Draw semi-transparent sphere
    sphere_mesh = Sphere(Point3f(0, 0, 0), Float32(radius))
    mesh!(ax, sphere_mesh, color=(:gray, sphere_alpha), transparency=true)

    # Draw origin and reference axes
    axis_len = radius * 1.3
    scatter!(ax, [Point3f(0, 0, 0)], color=:black, markersize=10)
    # X axis: azimuth=0 (east)
    lines!(ax, [Point3f(0, 0, 0), Point3f(axis_len, 0, 0)],
           color=:red, linewidth=3)
    text!(ax, Point3f(axis_len * 1.1, 0, 0), text="X (az=0°)",
          fontsize=14, color=:red)
    # Y axis: azimuth=-90° (north)
    lines!(ax, [Point3f(0, 0, 0), Point3f(0, axis_len, 0)],
           color=:green, linewidth=3)
    text!(ax, Point3f(0, axis_len * 1.1, 0), text="Y (az=-90°)",
          fontsize=14, color=:green)
    # Z axis: elevation=90° (up)
    lines!(ax, [Point3f(0, 0, 0), Point3f(0, 0, axis_len)],
           color=:blue, linewidth=3)
    text!(ax, Point3f(0, 0, axis_len * 1.1), text="Z (el=90°)",
          fontsize=14, color=:blue)

    # Default colors if not provided
    default_colors = [:blue, :red, :green, :orange, :purple, :cyan]
    actual_colors = isnothing(colors) ? default_colors : colors

    # Plot each trajectory
    for (i, lg) in enumerate(logs)
        sl = lg.syslog
        elevation = sl.elevation
        azimuth = sl.azimuth

        # Convert to Cartesian (ENU frame, azimuth_east convention)
        x = radius .* cos.(elevation) .* cos.(azimuth)
        y = -radius .* cos.(elevation) .* sin.(azimuth)
        z = radius .* sin.(elevation)

        color = actual_colors[mod1(i, length(actual_colors))]
        label = isnothing(labels) ? "trajectory $i" : labels[i]
        lines!(ax, x, y, z; color, linewidth, label)
    end

    # Add legend if multiple logs
    if length(logs) > 1 || !isnothing(labels)
        Legend(fig[1, 2], ax)
    end

    return fig
end

"""
    plot_body_frame(sys_struct; extra_points=nothing, dir=:side)

Plot wing points in 2D body frame coordinates.
Updates pos_b for PARTICLE_DYNAMICS wing points and shows all WING-type points.

# Arguments
- `sys_struct::SystemStructure`: System structure to plot
- `extra_points`: Optional vector of (x,y,z) tuples to overlay
- `dir::Symbol`: Viewing direction (:side, :front, or :top, default :side)

# Returns
- Figure with 2D scatter plot of wing points in body frame
"""
function SymbolicAWEModels.plot_body_frame(sys_struct::SystemStructure;
                         extra_points=nothing,
                         dir::Symbol=:front,
                         point_size=10,
                         extra_point_size=6,
                         figsize=(800, 600))
    (; points, wings, segments) = sys_struct

    # Update pos_b for PARTICLE_DYNAMICS wing points based on current wing orientation
    for wing in wings
        if wing.dynamics_type == SymbolicAWEModels.PARTICLE_DYNAMICS
            R_w_b = wing.R_b_to_w'  # transpose to get world-to-body
            for point in points
                if point.wing_idx == wing.idx
                    point.pos_b .= R_w_b * (point.pos_w - wing.pos_w)
                end
            end
        end
    end

    # Collect WING points
    wing_points = [p for p in points if p.type == SymbolicAWEModels.WING]

    # Extract coordinates and depth based on viewing direction
    # :front = looking in +x direction (shows y-z plane)
    # :side = looking in +y direction (shows x-z plane)
    # :top = looking in -z direction (shows x-y plane)
    if dir == :top
        coords = [(p.pos_b[1], p.pos_b[2]) for p in wing_points]
        depths = [-p.pos_b[3] for p in wing_points]  # -z: lower z = further
        xlabel, ylabel = "x [m]", "y [m]"
    elseif dir == :side
        coords = [(p.pos_b[1], p.pos_b[3]) for p in wing_points]
        depths = [p.pos_b[2] for p in wing_points]  # +y: higher y = further
        xlabel, ylabel = "x [m]", "z [m]"
    else  # :front (default)
        coords = [(p.pos_b[2], p.pos_b[3]) for p in wing_points]
        depths = [p.pos_b[1] for p in wing_points]  # +x: higher x = further
        xlabel, ylabel = "y [m]", "z [m]"
    end

    x_vals = [c[1] for c in coords]
    y_vals = [c[2] for c in coords]

    # Compute colors based on depth (closer = red, further = green)
    if !isempty(depths)
        d_min, d_max = extrema(depths)
        d_range = d_max - d_min
        if d_range > 0
            norm_d = [(d - d_min) / d_range for d in depths]
            # Closer (low depth) = red, further (high depth) = green
            point_colors = [RGBf(1.0 - nd, nd, 0.0) for nd in norm_d]
        else
            point_colors = fill(RGBf(0.5, 0.5, 0.0), length(depths))
        end
    else
        point_colors = [:red]
    end

    # Create figure
    fig = Figure(size=figsize)
    ax = Axis(fig[1, 1];
              xlabel, ylabel,
              title="Wing Points (Body Frame)",
              aspect=DataAspect())

    # Helper to get 2D coords from body frame position
    function get_2d(pos_b)
        if dir == :top
            return (pos_b[1], pos_b[2])
        elseif dir == :side
            return (pos_b[1], pos_b[3])
        else  # :front
            return (pos_b[2], pos_b[3])
        end
    end

    # Plot segments as thin gray lines
    wing_point_idxs = Set(p.idx for p in wing_points)
    for seg in segments
        # Only plot segments between wing points
        from_idx, to_idx = seg.point_idxs
        if from_idx in wing_point_idxs && to_idx in wing_point_idxs
            p1 = points[from_idx]
            p2 = points[to_idx]
            c1 = get_2d(p1.pos_b)
            c2 = get_2d(p2.pos_b)
            lines!(ax, [c1[1], c2[1]], [c1[2], c2[2]];
                   color=:gray, linewidth=0.5)
        end
    end

    # Plot wing points with depth-based coloring
    scatter!(ax, x_vals, y_vals;
             markersize=point_size,
             color=point_colors)

    # Add point indices as labels, placed away from other points
    for (i, p) in enumerate(wing_points)
        px, py = coords[i]

        # Calculate direction away from nearby points
        away_x, away_y = 0.0, 0.0
        for (j, _) in enumerate(wing_points)
            if i != j
                ox, oy = coords[j]
                dx, dy = px - ox, py - oy
                dist = sqrt(dx^2 + dy^2) + 0.01
                # Weight by inverse distance (closer points push harder)
                away_x += dx / dist^2
                away_y += dy / dist^2
            end
        end

        # Normalize and determine offset direction
        away_len = sqrt(away_x^2 + away_y^2)
        if away_len > 0
            away_x /= away_len
            away_y /= away_len
        else
            away_x, away_y = 1.0, 1.0
        end

        # Set offset and alignment based on direction
        offset_dist = 12
        offset = (offset_dist * sign(away_x), offset_dist * sign(away_y))
        align_x = away_x >= 0 ? :left : :right
        align_y = away_y >= 0 ? :bottom : :top

        text!(ax, px, py;
              text=string(p.idx),
              fontsize=12,
              align=(align_x, align_y),
              offset=offset)
    end

    # Plot extra points if provided
    if !isnothing(extra_points)
        # Transform extra points to body frame
        wing = wings[1]
        R_w_b = wing.R_b_to_w'

        extra_body = [R_w_b * (collect(p) - wing.pos_w) for p in extra_points]

        if dir == :top
            extra_coords = [(p[1], p[2]) for p in extra_body]
            extra_depths = [-p[3] for p in extra_body]
        elseif dir == :side
            extra_coords = [(p[1], p[3]) for p in extra_body]
            extra_depths = [p[2] for p in extra_body]
        else  # :front
            extra_coords = [(p[2], p[3]) for p in extra_body]
            extra_depths = [p[1] for p in extra_body]
        end

        ex_x = [c[1] for c in extra_coords]
        ex_y = [c[2] for c in extra_coords]

        # Depth-based colors for extra points (same red-to-green as sim)
        if !isempty(extra_depths)
            ed_min, ed_max = extrema(extra_depths)
            ed_range = ed_max - ed_min
            if ed_range > 0
                norm_ed = [(d - ed_min) / ed_range for d in extra_depths]
                extra_colors = [RGBf(1.0 - nd, nd, 0.0) for nd in norm_ed]
            else
                extra_colors = fill(RGBf(0.5, 0.5, 0.0), length(extra_depths))
            end
        else
            extra_colors = [:red]
        end

        scatter!(ax, ex_x, ex_y;
                 markersize=extra_point_size + 6,
                 color=extra_colors,
                 marker=:cross)
    end

    # Auto-zoom to fit all WING points with margin for labels
    if !isempty(x_vals)
        x_min, x_max = extrema(x_vals)
        y_min, y_max = extrema(y_vals)
        # Larger margin to accommodate point index labels
        margin_x = 0.15 * (x_max - x_min) + 0.3
        margin_y = 0.15 * (y_max - y_min) + 0.3
        limits!(ax, x_min - margin_x, x_max + margin_x,
                    y_min - margin_y, y_max + margin_y)
    end

    # Custom legend
    legend_elements = [
        MarkerElement(color=:black, marker=:circle, markersize=10),
    ]
    legend_labels = ["sim"]

    if !isnothing(extra_points)
        push!(legend_elements,
              MarkerElement(color=:black, marker=:cross, markersize=12))
        push!(legend_labels, "photogrammetry")
    end

    # Depth color scale as two entries
    push!(legend_elements,
          MarkerElement(color=RGBf(1, 0, 0), marker=:rect, markersize=10))
    push!(legend_labels, "close")
    push!(legend_elements,
          MarkerElement(color=RGBf(0, 1, 0), marker=:rect, markersize=10))
    push!(legend_labels, "far")

    Legend(fig[1, 2], legend_elements, legend_labels)

    return fig
end

"""
    plot_aoa(sys_struct; wing_idx=1, figsize=(800, 400))

Plot the angle of attack distribution along the wing span.

# Arguments
- `sys_struct::SystemStructure`: System structure containing the wing with VSM solver
- `wing_idx::Int`: Index of the wing to plot (default: 1)
- `figsize`: Figure size tuple (default: (800, 400))

# Returns
- Figure with angle of attack vs panel index
"""
function SymbolicAWEModels.plot_aoa(sys_struct::SystemStructure;
                                    wing_idx::Int=1,
                                    figsize=(800, 400))
    wing = sys_struct.wings[wing_idx]

    # Get alpha distribution (in radians)
    alpha_dist = wing.vsm_solver.sol.alpha_dist

    # Use panel indices for x-axis
    panel_indices = 1:length(alpha_dist)

    # Convert to degrees for plotting
    alpha_deg = rad2deg.(alpha_dist)

    fig = Figure(size=figsize)
    ax = Axis(fig[1, 1];
              xlabel="Panel index",
              ylabel="Angle of attack [deg]",
              title="Angle of Attack Distribution")

    lines!(ax, panel_indices, alpha_deg; linewidth=2)
    scatter!(ax, panel_indices, alpha_deg; markersize=6)

    return fig
end

end
