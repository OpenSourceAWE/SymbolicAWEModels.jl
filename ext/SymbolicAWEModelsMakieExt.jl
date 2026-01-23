# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

module SymbolicAWEModelsMakieExt

using Makie
using UnPack
using LinearAlgebra
using StaticArrays
using Statistics
using Printf
using KiteUtils
using KiteUtils: SysLog
using SymbolicAWEModels
using VortexStepMethod

# Global storage for plot observables (for time-based plotting)
const PLOT_GEOMETRY_OBS = Ref{Union{Nothing, Observable}}(nothing)  # Single trigger observable
const PLOT_SEGMENT_COLORS_OBS = Ref{Union{Nothing, Observable}}(nothing)  # Separate for force coloring
const PLOT_SCENE = Ref{Union{Nothing, Scene}}(nothing)
const PLOT_BACKGROUND_PANES = Ref{Union{Nothing, Vector}}(nothing)
const PLOT_MARGIN = Ref{Float64}(10.0)
const PLOT_RELEVANT_PLOTS = Ref{Union{Nothing, Vector}}(nothing)
const PLOT_SYSTEM_STRUCTURE = Ref{Union{Nothing, SystemStructure}}(nothing)
const PLOT_VECTOR_SCALE = Ref{Float64}(1.0)
const PLOT_FORCE_COLOR = Ref{Bool}(false)
const PLOT_SEGMENT_COLOR = Ref{Symbol}(:black)
const PLOT_ZOOMED_IN = Ref{Bool}(false)
const PLOT_ZOOM_RELMARGIN = Ref{Float64}(0.2)
const PLOT_ZOOM_SEGMENT_IDX = Ref{Int}(-1)  # Which segment we're zoomed into (-1 = none)
const PLOT_BODY_FRAME = Ref{Bool}(false)  # Whether body frame tracking is active
const PLOT_CAMERA_DISTANCE = Ref{Union{Nothing, Float64}}(nothing)  # Stored camera distance
const PLOT_PREV_BODY_FRAME = Ref{Bool}(false)  # Previous body frame state
const PLOT_PREV_ZOOMED_IN = Ref{Bool}(false)  # Previous zoomed state
const PLOT_PREV_SEGMENT_IDX = Ref{Int}(-1)  # Previous segment index

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

function Makie.plot!(ax, sys::SystemStructure;
                     point_color = :darkred, segment_color = :black,
                     wing_colors = Makie.wong_colors(), vector_scale = 1.0,
                     show_points = true, show_segments = true, show_orient = true,
                     show_panes = true, margin = 10.0, force_color = false,
                     plot_vsm = true, plot_aero = true,
                     # Optional observable for real-time updates
                     geometry_obs = nothing,
                     pane_observables = nothing)

    plots = Dict{Symbol, Any}()

    # === Plot Segments ===
    if show_segments
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

        plots[:segments] = linesegments!(ax, lineseg_points, color=seg_colors,
                                         linewidth=2, label="Segments", transparency=true)
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
        plots[:points] = scatter!(ax, point_positions, color=point_color, label="Points",
                                  transparency=true)
    end

    # === Plot Wings ===
    if show_orient
        if isnothing(geometry_obs)
            # Static plotting: create separate arrows for each wing
            plots[:wings] = []
            for (i, wing) in enumerate(sys.wings)
                wing_pos = Point3f(wing.pos_w)
                R = wing.R_b_w
                scale = vector_scale
                origins = [wing_pos, wing_pos, wing_pos]
                directions = [Vec3f(R[:, 1]) * scale, Vec3f(R[:, 2]) * scale, Vec3f(R[:, 3]) * scale]

                axis_colors = [:red, :green, :blue]
                p = arrows3d!(ax, origins, directions, color=axis_colors, label="Wing $i Axes")
                push!(plots[:wings], p)
            end
        else
            # Dynamic plotting: compute from PLOT_SYSTEM_STRUCTURE when triggered
            wing_origins_dirs = @lift begin
                $geometry_obs  # Trigger dependency
                sys_ref = PLOT_SYSTEM_STRUCTURE[]
                scale = PLOT_VECTOR_SCALE[]
                origins = Point3f[]
                directions = Vec3f[]
                for wing in sys_ref.wings
                    wing_pos = Point3f(wing.pos_w)
                    R = wing.R_b_w
                    # Add three arrow vectors for each axis (x, y, z in body frame)
                    for i in 1:3
                        push!(origins, wing_pos)
                        push!(directions, Vec3f(R[:, i]) * scale)
                    end
                end
                (origins, directions)
            end
            axis_colors = repeat([:red, :green, :blue], length(sys.wings))
            plots[:wings] = arrows3d!(ax, @lift($wing_origins_dirs[1]), @lift($wing_origins_dirs[2]),
                                     color=axis_colors)
        end
    end

    # === Plot VSM Aerodynamics ===
    # VSM panels - use observables if we're in dynamic mode
    if plot_vsm && !isempty(sys.wings)
        plots[:vsm] = []
        use_obs = !isnothing(geometry_obs)  # If geometry observable exists, use it for VSM too
        for (i, wing) in enumerate(sys.wings)
            p = plot!(ax, wing.vsm_aero; R_b_w=wing.R_b_w, T_b_w=wing.pos_w, use_observables=use_obs)
            push!(plots[:vsm], p)
        end
    end

    # === Plot Aero Forces ===
    if plot_aero
        if isnothing(geometry_obs)
            # Static plotting: build aero force arrows
            aero_origins = Point3f[]
            aero_forces_raw = Vec3f[]

            for wing in sys.wings
                if wing.wing_type == QUATERNION
                    # For QUATERNION wings, use wing.aero_force_b
                    if !iszero(wing.aero_force_b)
                        aero_force_w = wing.R_b_w * wing.aero_force_b
                        push!(aero_origins, Point3f(wing.pos_w))
                        push!(aero_forces_raw, Vec3f(aero_force_w))
                    end
                elseif wing.wing_type == SymbolicAWEModels.REFINE
                    # For REFINE wings, plot both point forces and total wing force
                    # Plot individual point forces
                    for point in sys.points
                        if point.type == WING && point.wing_idx == wing.idx
                            if !iszero(point.aero_force_b)
                                aero_force_w = wing.R_b_w * point.aero_force_b
                                push!(aero_origins, Point3f(point.pos_w))
                                push!(aero_forces_raw, Vec3f(aero_force_w))
                            end
                        end
                    end
                    # Also plot total wing aero force at wing center when vector_scale > 0
                    if vector_scale > 0 && !iszero(wing.aero_force_b)
                        aero_force_w = wing.R_b_w * wing.aero_force_b
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
                    if wing.wing_type == QUATERNION
                        if !iszero(wing.aero_force_b)
                            aero_force_w = wing.R_b_w * wing.aero_force_b
                            push!(origins, Point3f(wing.pos_w))
                            push!(forces_raw, Vec3f(aero_force_w))
                        end
                    elseif wing.wing_type == SymbolicAWEModels.REFINE
                        for point in sys_ref.points
                            if point.type == WING && point.wing_idx == wing.idx
                                if !iszero(point.aero_force_b)
                                    aero_force_w = wing.R_b_w * point.aero_force_b
                                    push!(origins, Point3f(point.pos_w))
                                    push!(forces_raw, Vec3f(aero_force_w))
                                end
                            end
                        end
                        if scale > 0 && !iszero(wing.aero_force_b)
                            aero_force_w = wing.R_b_w * wing.aero_force_b
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

    # === Calculate system scale for axes and panes ===
    xlims, ylims, zlims = (-10, 10), (-10, 10), (-10, 10) # Default limits
    if !isempty(sys.points)
        all_x = [p.pos_w[1] for p in sys.points]
        all_y = [p.pos_w[2] for p in sys.points]
        all_z = [p.pos_w[3] for p in sys.points]

        xlims_data = extrema(all_x)
        ylims_data = extrema(all_y)
        zlims_data = extrema(all_z)

        xlims = (xlims_data[1] - margin, xlims_data[2] + margin)
        ylims = (ylims_data[1] - margin, ylims_data[2] + margin)
        zlims = (zlims_data[1] - margin, zlims_data[2] + margin)
    end

    # Calculate characteristic length for scaling arrows
    char_length = max(xlims[2] - xlims[1], ylims[2] - ylims[1], zlims[2] - zlims[1])
    axis_length = char_length * 0.15

    # === Plot Global Axes ===
    begin
        origins = [Point3f(0, 0, 0), Point3f(0, 0, 0), Point3f(0, 0, 0)]
        directions = [Vec3f(axis_length, 0, 0), Vec3f(0, axis_length, 0), Vec3f(0, 0, axis_length)]
        axis_colors = [:red, :green, :blue]
        # Use fixed radii so arrows don't get fatter with longer tethers
        plots[:global_axes] = arrows3d!(ax, origins, directions;
                                        shaftradius=0.02,
                                        tipradius=0.08,
                                        tiplength=0.2,
                                        color=axis_colors,
                                        label="Global Axes")
    end

    # === Plot Background Panes ===
    if show_panes
        pane_color = RGBAf(0.95, 0.95, 0.95, 0.3)
        pane_extent = 10000.0f0  # Large value for "infinite" extent

        if isnothing(pane_observables)
            # Static plotting: create new observables
            xz_pane_obs = Observable(Rect3(Vec3f(xlims[1] - pane_extent, ylims[2], zlims[1]),
                                            Vec3f(xlims[2] - xlims[1] + pane_extent, 0.01, pane_extent)))
            yz_pane_obs = Observable(Rect3(Vec3f(xlims[2], ylims[1] - pane_extent, zlims[1]),
                                            Vec3f(0.01, ylims[2] - ylims[1] + pane_extent, pane_extent)))
            xy_pane_obs = Observable(Rect3(Vec3f(xlims[1] - pane_extent, ylims[1] - pane_extent, zlims[1]),
                                            Vec3f(xlims[2] - xlims[1] + pane_extent, ylims[2] - ylims[1] + pane_extent, 0.01)))
            plots[:pane_observables] = [xz_pane_obs, yz_pane_obs, xy_pane_obs]
        else
            # Dynamic plotting: use provided observables and update them
            xz_pane_obs, yz_pane_obs, xy_pane_obs = pane_observables
            xz_pane_obs[] = Rect3(Vec3f(xlims[1] - pane_extent, ylims[2], zlims[1]),
                                  Vec3f(xlims[2] - xlims[1] + pane_extent, 0.01, pane_extent))
            yz_pane_obs[] = Rect3(Vec3f(xlims[2], ylims[1] - pane_extent, zlims[1]),
                                  Vec3f(0.01, ylims[2] - ylims[1] + pane_extent, pane_extent))
            xy_pane_obs[] = Rect3(Vec3f(xlims[1] - pane_extent, ylims[1] - pane_extent, zlims[1]),
                                  Vec3f(xlims[2] - xlims[1] + pane_extent, ylims[2] - ylims[1] + pane_extent, 0.01))
            plots[:pane_observables] = pane_observables
        end

        # Create mesh plots for the panes
        xz_pane = mesh!(ax, plots[:pane_observables][1], color=pane_color)
        yz_pane = mesh!(ax, plots[:pane_observables][2], color=pane_color)
        xy_pane = mesh!(ax, plots[:pane_observables][3], color=pane_color)
        plots[:panes] = [xz_pane, yz_pane, xy_pane]
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

    # Update VSM panel meshes
    if !isnothing(PLOT_GEOMETRY_OBS[])
        for wing in sys.wings
            plot!(wing.vsm_aero; R_b_w=wing.R_b_w, T_b_w=wing.pos_w)
        end
    end

    # Auto-update camera to keep geometry centered (runs every frame during replay)
    # Re-runs the appropriate zoom function to track moving geometry
    if !isnothing(PLOT_SCENE[]) && !isnothing(PLOT_RELEVANT_PLOTS[]) && !isnothing(PLOT_SYSTEM_STRUCTURE[])
        scene = PLOT_SCENE[]
        relevant_plots = PLOT_RELEVANT_PLOTS[]
        stored_sys = PLOT_SYSTEM_STRUCTURE[]

        # Detect mode change
        mode_changed = (PLOT_PREV_BODY_FRAME[] != PLOT_BODY_FRAME[] ||
                        PLOT_PREV_ZOOMED_IN[] != PLOT_ZOOMED_IN[] ||
                        PLOT_PREV_SEGMENT_IDX[] != PLOT_ZOOM_SEGMENT_IDX[])

        if mode_changed
            PLOT_CAMERA_DISTANCE[] = nothing  # Force recalculation
            # Update prev state
            PLOT_PREV_BODY_FRAME[] = PLOT_BODY_FRAME[]
            PLOT_PREV_ZOOMED_IN[] = PLOT_ZOOMED_IN[]
            PLOT_PREV_SEGMENT_IDX[] = PLOT_ZOOM_SEGMENT_IDX[]
        end

        # Call zoom functions with stored distance (preserves manual zoom)
        if PLOT_BODY_FRAME[]
            # Body frame mode: continuously track wing orientation
            dist = zoom_body_frame!(scene, scene.camera, stored_sys, PLOT_CAMERA_DISTANCE[])
            PLOT_CAMERA_DISTANCE[] = dist
        elseif PLOT_ZOOMED_IN[] && PLOT_ZOOM_SEGMENT_IDX[] > 0
            # When zoomed in, keep segment centered as it moves
            dist = zoom_in!(scene, scene.camera, stored_sys, PLOT_ZOOM_SEGMENT_IDX[], PLOT_CAMERA_DISTANCE[])
            PLOT_CAMERA_DISTANCE[] = dist
        elseif !PLOT_ZOOMED_IN[]
            # When not zoomed in, keep full view centered on geometry
            dist = zoom_out!(scene, scene.camera, relevant_plots, PLOT_CAMERA_DISTANCE[]; relmargin=PLOT_ZOOM_RELMARGIN[])
            PLOT_CAMERA_DISTANCE[] = dist
        end
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

normalize_safe(v::SVector{3, Float64}; eps=1e-12) = begin
    n = norm(v)
    if !isfinite(n) || n < eps
        return (SVector{3, Float64}(NaN, NaN, NaN), n)
    end
    return (v / n, n)
end

function bridle_frame_va(va_enu::SVector{3, Float64},
                         tension_enu::SVector{3, Float64};
                         pos_enu::Union{Nothing, SVector{3, Float64}}=nothing,
                         eps=1e-12,
                         ey_prev::Union{Nothing, SVector{3, Float64}}=nothing)
    # Use provided tension direction; higher-level code enforces sign continuity
    ez, nT = normalize_safe(tension_enu; eps=eps)
    if !isfinite(nT) || nT < eps
        return (SMatrix{3, 3, Float64}(ntuple(_ -> NaN, 9)), ey_prev)
    end

    ey_raw = cross(ez, -va_enu)
    ey, nEy = normalize_safe(ey_raw; eps=eps)
    # Keep ey continuous to avoid frame flips when va grazes ez
    if ey_prev !== nothing && isfinite(nEy) && nEy >= eps && dot(ey, ey_prev) < 0
        ey = -ey
    end
    if !isfinite(nEy) || nEy < eps
        if ey_prev !== nothing
            ey = ey_prev
        else
            tmp = abs(ez[1]) < 0.9 ? SVector(1.0, 0.0, 0.0) : SVector(0.0, 1.0, 0.0)
            ey, _ = normalize_safe(cross(ez, tmp); eps=eps)
        end
    end

    ex = cross(ey, ez)
    DCM = @SMatrix [ex[1] ey[1] ez[1];
                    ex[2] ey[2] ez[2];
                    ex[3] ey[3] ez[3]]
    return DCM, ey
end

function yaw_from_dcm_ned(DCM::SMatrix{3, 3, Float64}; eps=1e-9)
    x = clamp(-DCM[3, 1], -1.0, 1.0)
    pitch = asin(x)
    if abs(abs(x) - 1.0) < eps
        yaw = atan(-DCM[1, 2], DCM[2, 2])
        roll = 0.0
        return roll, pitch, yaw
    end
    roll = atan(DCM[3, 2], DCM[3, 3])
    yaw = atan(DCM[2, 1], DCM[1, 1])
    return roll, pitch, yaw
end

function unwrap_phase!(vals::AbstractVector{<:Real}; period=2π, thresh=π)
    if isempty(vals)
        return vals
    end
    offset = 0.0
    prev = vals[1]
    for i in 2:length(vals)
        δ = vals[i] - prev
        if δ > thresh
            offset -= period
        elseif δ < -thresh
            offset += period
        end
        prev = vals[i]
        vals[i] += offset
    end
    return vals
end

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
function compute_turn_radius(sl_in, sys::SystemStructure; smooth_window=10, eps=1e-12)
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
    @inbounds for k in 1:n
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
    @inbounds for k in 1:n
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
    midle_to_kcu_dir(sl, k; eps=1e-12)

Compute the unit vector from the mid leading-edge (avg of points 12 & 14)
to the KCU/bridle hub (point 1) for sample `k` of a syslog entry.
Returns `nothing` if the required points are unavailable or degenerate.
"""
function midle_to_kcu_dir(sl, k; eps=1e-12)
    Xk = sl.X[k]; Yk = sl.Y[k]; Zk = sl.Z[k]
    if length(Xk) < 14 || length(Yk) < 14 || length(Zk) < 14
        return nothing
    end
    p1 = SVector{3, Float64}(Xk[1], Yk[1], Zk[1])
    ple12 = SVector{3, Float64}(Xk[12], Yk[12], Zk[12])
    ple14 = SVector{3, Float64}(Xk[14], Yk[14], Zk[14])
    p_le_mid = (ple12 + ple14) / 2
    dir = p1 - p_le_mid
    n = norm(dir)
    return n > eps ? dir / n : nothing
end

function calc_ref_area(sys::SystemStructure)
    isempty(sys.wings) && return NaN
    wing = sys.wings[1]
    hasproperty(wing, :vsm_aero) || return NaN
    panels = wing.vsm_aero.panels
    isempty(panels) && return NaN
    return sum(p.chord * p.width for p in panels)
end

function calculate_cs(sl_in, sys; rho=1.225, eps=1e-12)
    sl = hasproperty(sl_in, :syslog) ? sl_in.syslog : sl_in
    s_ref = calc_ref_area(sys)
    if !isfinite(s_ref) || s_ref <= eps
        return Float64[], Float64[]
    end

    n = length(sl.time)
    if n == 0 ||
       length(sl.vel_kite) < n || length(sl.v_wind_kite) < n ||
       length(sl.orient) < n || length(sl.aero_force_b) < n ||
       length(sl.X) < n || length(sl.Y) < n || length(sl.Z) < n
        return Float64[], Float64[]
    end

    cs = Vector{Float64}(undef, n)

    @inbounds for k in 1:n
        v_kite = sl.vel_kite[k]
        v_wind = sl.v_wind_kite[k]
        v_a = v_kite - v_wind
        v_a_norm = norm(v_a)
        if v_a_norm <= eps
            cs[k] = NaN
            continue
        end
        drag_dir = -v_a / v_a_norm

        up_dir = midle_to_kcu_dir(sl, k; eps=eps)
        if up_dir === nothing
            cs[k] = NaN
            continue
        end
        up_dir = -up_dir

        side_raw = cross(drag_dir, up_dir)
        side_norm = norm(side_raw)
        if side_norm <= eps
            cs[k] = NaN
            continue
        end
        side_dir = side_raw / side_norm

        R_b_w = SymbolicAWEModels.quaternion_to_rotation_matrix(sl.orient[k])
        F_aero_b = sl.aero_force_b[k]
        F_aero_w = R_b_w * F_aero_b
        F_side = dot(F_aero_w, side_dir)
        cs[k] = F_side / (0.5 * rho * v_a_norm^2 * s_ref)
    end

    return cs, sl.time
end

"""
    compute_ekf_yaw_and_rate(sl_in, sys::SystemStructure; eps=1e-12)

Compute EKF-style yaw and yaw rate using a velocity-based tangent frame.

This matches the Python reference implementation which uses the kite's velocity
direction in the tangent plane (perpendicular to position) to define heading.
This approach is more robust than using apparent wind × tension because it
avoids singularities when these vectors align.

# Returns
- `(yaw, yaw_rate)`: Tuple of unwrapped yaw angles [rad] and yaw rates [deg/s]
- `nothing` if required data is missing
"""
function compute_ekf_yaw_and_rate(sl_in, sys::SystemStructure; eps=1e-12)
    # Accept either a SysLog wrapper or the raw syslog StructVector
    sl = hasproperty(sl_in, :syslog) ? sl_in.syslog : sl_in
    n = length(sl.time)
    if n < 2 || isempty(sl.vel_kite)
        return nothing
    end
    if length(sys.wings) == 0 || length(sl.X) < n || length(sl.Y) < n || length(sl.Z) < n
        return nothing
    end
    
    kite_idx = sys.wings[1].origin_idx
    yaw = Vector{Float64}(undef, n)
    nan_count = 0
    
    # Use velocity-based tangent frame (same as HeadingGate/sphere method)
    # This is more robust than tension × apparent wind
    @inbounds for k in 1:n
        pos = SVector{3, Float64}(sl.X[k][kite_idx], sl.Y[k][kite_idx], sl.Z[k][kite_idx])
        vel = SVector{3, Float64}(sl.vel_kite[k])
        
        npos = norm(pos)
        nvel = norm(vel)
        
        if npos > eps && nvel > eps
            # Define tangent plane perpendicular to radial direction
            radial = pos / npos
            
            # Project velocity into tangent plane
            tang_vel = vel - (vel ⋅ radial) * radial
            ntang = norm(tang_vel)
            
            if ntang > eps
                tang_vel_unit = tang_vel / ntang
                
                # Build local "up" frame at kite position
                # up_z = radial (points away from origin)
                # up_y = perpendicular in horizontal plane
                # up_x = up_y × up_z
                up_z = radial
                up_y_raw = SVector(-pos[2], pos[1], 0.0)
                nup_y = norm(up_y_raw)
                
                if nup_y > eps
                    up_y = up_y_raw / nup_y
                    up_x = up_z × up_y
                    nup_x = norm(up_x)
                    
                    if nup_x > eps
                        up_x = up_x / nup_x
                        # Re-orthonormalize
                        up_y = up_z × up_x
                        
                        # Rotation matrix from world to "up" frame
                        R_up = @SMatrix [up_x[1] up_y[1] up_z[1];
                                         up_x[2] up_y[2] up_z[2];
                                         up_x[3] up_y[3] up_z[3]]
                        
                        # Express tangent velocity in "up" frame
                        heading_vec = R_up' * tang_vel_unit
                        
                        # Yaw from x-y components in "up" frame
                        yaw[k] = atan(heading_vec[2], heading_vec[1])
                        continue
                    end
                end
            end
        end
        
        # Fallback: copy previous or NaN
        yaw[k] = k > 1 ? yaw[k-1] : NaN
        nan_count += 1
    end
    
    if nan_count > 0
        @info "compute_ekf_yaw_and_rate: $nan_count samples with degenerate geometry"
    end
    
    # Unwrap phase
    yaw_unwrapped = copy(yaw)
    unwrap_phase!(yaw_unwrapped)
    
    # Compute yaw rate
    ts = mean(diff(sl.time))
    ts = isfinite(ts) && ts > eps ? ts : eps
    
    yaw_rate = gradient_uniform(yaw_unwrapped, ts)
    # Match Python's smoothing window
    yaw_rate = moving_average_same(yaw_rate, 10)
    
    @info "EKF yaw rate stats" mean=mean(yaw_rate) std=std(yaw_rate) min=minimum(yaw_rate) max=maximum(yaw_rate)
    
    return yaw_unwrapped, rad2deg.(yaw_rate)
end


"""
    compute_ekf_yaw_and_rate_tension(sl_in, sys::SystemStructure; eps=1e-12)

Alternative implementation using tension × apparent wind frame.

This version uses the bridle frame (tension direction as z-axis, apparent wind
to define y-axis). It's more physically motivated but can be less robust when
apparent wind aligns with tension.

# Returns  
- `(yaw, yaw_rate)`: Tuple of unwrapped yaw angles [rad] and yaw rates [deg/s]
- `nothing` if required data is missing
"""
function compute_ekf_yaw_and_rate_tension(sl_in, sys::SystemStructure; eps=1e-12)
    sl = hasproperty(sl_in, :syslog) ? sl_in.syslog : sl_in
    n = length(sl.time)
    if n < 2 || isempty(sl.vel_kite) || isempty(sl.v_wind_kite) || isempty(sl.tether_induced_force)
        return nothing
    end
    if length(sys.wings) == 0 ||
       length(sl.vel_kite) < n || length(sl.v_wind_kite) < n ||
       length(sl.tether_induced_force) < n || length(sl.X) < n ||
       length(sl.Y) < n || length(sl.Z) < n
        return nothing
    end
    
    kite_idx = sys.wings[1].origin_idx
    yaw = Vector{Float64}(undef, n)
    ey_prev = nothing
    ex_prev = nothing
    tension_prev = nothing
    missing_tension = 0
    missing_pos = 0
    nan_yaw = 0
    flip_count = 0
    t_min = Inf
    t_max = -Inf
    t_sum = 0.0
    t_cnt = 0
    
    @inbounds for k in 1:n
        v_kite = SVector{3, Float64}(sl.vel_kite[k])
        v_wind = SVector{3, Float64}(sl.v_wind_kite[k])
        pos = SVector{3, Float64}(sl.X[k][kite_idx], sl.Y[k][kite_idx], sl.Z[k][kite_idx])
        tension_raw = SVector{3, Float64}(sl.tether_induced_force[k])

        # Prefer geometry-based bridle direction
        tension_dir = midle_to_kcu_dir(sl, k; eps=eps)
        
        # Enforce continuity on tension direction
        if tension_dir !== nothing && all(isfinite, tension_dir) && norm(tension_dir) >= eps
            if tension_prev !== nothing && dot(tension_dir, tension_prev) < 0
                tension_dir = -tension_dir
                flip_count += 1
            end
        else
            if tension_prev !== nothing
                tension_dir = tension_prev
            else
                tension_dir = tension_raw
                nT_raw = norm(tension_raw)
                if !isfinite(nT_raw) || nT_raw < eps
                    npos = norm(pos)
                    if isfinite(npos) && npos >= eps
                        tension_dir = -pos / npos
                        missing_tension += 1
                    else
                        missing_pos += 1
                        yaw[k] = k > 1 ? yaw[k - 1] : NaN
                        continue
                    end
                else
                    t_min = min(t_min, nT_raw)
                    t_max = max(t_max, nT_raw)
                    t_sum += nT_raw
                    t_cnt += 1
                    if tension_prev !== nothing && dot(tension_dir, tension_prev) < 0
                        tension_dir = -tension_dir
                        flip_count += 1
                    end
                end
            end
        end
        tension_prev = tension_dir

        va = v_wind - v_kite
        DCM_b2w, ey_prev = bridle_frame_va(va, tension_dir; pos_enu=pos, ey_prev=ey_prev, eps=eps)
        
        # Check ex continuity to detect frame flips
        ex_curr = SVector{3, Float64}(DCM_b2w[1, 1], DCM_b2w[2, 1], DCM_b2w[3, 1])
        if ex_prev !== nothing && dot(ex_curr, ex_prev) < 0
            # Flip ex and ez (keep ey)
            DCM_b2w = @SMatrix [-DCM_b2w[1,1]  DCM_b2w[1,2] -DCM_b2w[1,3];
                                -DCM_b2w[2,1]  DCM_b2w[2,2] -DCM_b2w[2,3];
                                -DCM_b2w[3,1]  DCM_b2w[3,2] -DCM_b2w[3,3]]
            ex_curr = SVector{3, Float64}(DCM_b2w[1, 1], DCM_b2w[2, 1], DCM_b2w[3, 1])
        end
        ex_prev = ex_curr
        
        # Convert to NED and extract yaw
        DCM_b2ned = R_ENU2NED * DCM_b2w
        yaw_raw = atan(DCM_b2ned[2, 1], DCM_b2ned[1, 1])
        
        # Unwrap-style continuity
        if k == 1
            yaw[k] = yaw_raw
        else
            diff = yaw_raw - yaw[k-1]
            while diff > π
                diff -= 2π
            end
            while diff < -π
                diff += 2π
            end
            yaw[k] = yaw[k-1] + diff
        end
        
        nan_yaw += isfinite(yaw[k]) ? 0 : 1
    end
    
    if missing_tension > 0 || missing_pos > 0 || nan_yaw > 0 || flip_count > 0
        t_mean = t_cnt > 0 ? t_sum / t_cnt : NaN
        @info "compute_ekf_yaw_and_rate_tension stats" missing_tension missing_pos nan_yaw flip_count t_min t_max t_mean
    end
    
    ts = mean(diff(sl.time))
    ts = isfinite(ts) && ts > eps ? ts : eps
    
    yaw_rate = gradient_uniform(yaw, ts)
    yaw_rate = moving_average_same(yaw_rate, 10)
    
    if missing_tension > 0 || missing_pos > 0
        @info "EKF yaw rate (tension) stats" mean=mean(yaw_rate) std=std(yaw_rate) min=minimum(yaw_rate) max=maximum(yaw_rate)
    end
    
    return yaw, rad2deg.(yaw_rate)
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
- `plot_yaw_rate_paper::Bool=false`: Plot EKF bridle yaw rate ψ̇ derived from the ENU tension/apparent wind frame (paper yaw-rate).
- `yaw_rate_paper_ylims::Union{Nothing,Tuple}=nothing`: Y-axis limits for the paper yaw-rate panel (`nothing` leaves autoscaling).
- `plot_elevation::Bool=false`: Show the panel with the kite's elevation angle.
- `plot_azimuth::Bool=false`: Show the panel with the kite's azimuth angle.
- `plot_distance::Bool=false`: Show the panel with the kite distance from origin (norm of position).
- `plot_yaw_rate::Bool=false`: Show yaw rate `dψ/dt` derived from the wind-referenced heading.
- `plot_gk_paper::Bool=false`: Plot gk using paper-style ψ̇ and reconstructed steering command.
- `plot_cs::Bool=false`: Plot the side-force coefficient C_S based on apparent wind and logged aero force.
- `cs_ylims::Union{Nothing,Tuple}=(-0.5, 0.5)`: Y-axis limits for the C_S panel (`nothing` leaves autoscaling).
- `turn_radius_ylims::Union{Nothing,Tuple}=(0.0, 100.0)`: Y-axis limits for the turn-radius panel (`nothing` leaves autoscaling).
- `yaw_rate_paper_compare::Bool=false`: Log std/offset comparisons between yaw definitions in the paper panel.
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

This method allows plotting multiple syslogs (e.g., from REFINE and QUATERNION models)
on the same panels for direct comparison. Each log's traces are labeled with its
corresponding system name.

# Arguments
- `sys::SystemStructure`: The system structure (can be from any of the models).
- `logs::Vector{SysLog}`: Vector of simulation logs to compare.

# Keyword Arguments
Same as the single-syslog version. See `Makie.plot(sys::SystemStructure, lg::SysLog)` for details.

# Example
```julia
# Compare REFINE vs QUATERNION models
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

function Makie.plot(syss::Vector{SystemStructure}, logs::Vector{<:SysLog};
                   plot_default=false,
                   plot_reelout=plot_default,
                   plot_aero_force=plot_default,
                   plot_twist=false,
                   plot_us=false,
                   plot_gk=false,
                   gk_ylims=(0.0, 15.0),
                   plot_yaw_rate=false,
                   plot_yaw_rate_paper=false,
                   yaw_rate_paper_ylims=(-90.0, 90.0),
                   yaw_rate_paper_compare=false,
                   plot_gk_paper=false,
                   plot_cs=false,
                   cs_ylims=(-0.0125, 0.0125),
                   turn_radius_ylims=(0.0, 70.0),
                   plot_v_app=false,
                   plot_kite_vel=false,
                   plot_aoa=plot_default,
                   aoa_ylims=(0.0, 15.0),
                   plot_heading=plot_default,
                   plot_kiteutils_course=false,
                   plot_aero_moment=false,
                   plot_turn_rates=false,
                   plot_turn_radius=false,
                   plot_elevation=false,
                   plot_azimuth=false,
                   plot_tether_moment=false,
                   plot_winch_force=plot_default,
                   plot_set_values=false,
                   plot_distance=false,
                   plot_cone_angle=false,
                   plot_old_heading=false,
                   plot_tether=false,
                   plot_tether_actual=false,
                   heading_setpoint=nothing,
                   tether_len_setpoint=nothing,
                   tape_lengths=nothing,
                   suffixes=nothing,
                   size=(1200, 800))

    # Build list of panels to plot by combining data from all logs
    panels = []

    # Generate suffixes: use custom if provided, otherwise use system names
    actual_suffixes = if isnothing(suffixes)
        [" - " * sys.name for sys in syss]
    else
        [" - " * s for s in suffixes]
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
            for j in 2:length(heading_unwrapped)
                while heading_unwrapped[j] - heading_unwrapped[j-1] > π
                    heading_unwrapped[j] -= 2π
                end
                while heading_unwrapped[j] - heading_unwrapped[j-1] < -π
                    heading_unwrapped[j] += 2π
                end
            end
            heading_rate = diff(rad2deg.(heading_unwrapped)) ./ diff(sl.time)
            push!(all_data, heading_rate)
            push!(all_labels, "ψ̇" * suffix)
            push!(all_times, sl.time[1:end-1])
        end
        push!(panels, (
            data = all_data,
            labels = all_labels,
            times = all_times,
            ylabel = "dψ/dt [°/s]\nheading⟂wind"
        ))
    end

    if plot_yaw_rate_paper
        all_data = []
        all_labels = []
        all_times = []
        compare_results = []
        for (i, lg) in enumerate(logs)
            sl = lg.syslog
            suffix = " - " * syss[i].name
            n = length(sl.time)
            if isempty(sl.orient) || isempty(sl.X)
                @warn "Missing orient or position data in syslog; skipping yaw_rate_paper for $suffix"
                continue
            end

            # Compute heading rate from diff for quaternion wings
            heading_unwrapped = copy(sl.heading)
            for j in 2:length(heading_unwrapped)
                while heading_unwrapped[j] - heading_unwrapped[j-1] > π
                    heading_unwrapped[j] -= 2π
                end
                while heading_unwrapped[j] - heading_unwrapped[j-1] < -π
                    heading_unwrapped[j] += 2π
                end
            end
            yaw_rate_perp_wind = Base.diff(rad2deg.(heading_unwrapped)) ./ Base.diff(sl.time)

            yaw_ekf = nothing
            yaw_rate_EKF = nothing
            if !isempty(sl.vel_kite) && !isempty(sl.v_wind_kite) && !isempty(sl.tether_induced_force)
                ekf = compute_ekf_yaw_and_rate(lg, syss[i])
                if ekf !== nothing
                    yaw_ekf, yaw_rate_EKF = ekf
                end
            end

            # Tangential sphere heading: use velocity projected into tangent plane
            yaw_sphere = Vector{Float64}(undef, n)
            kite_idx = syss[i].wings[1].origin_idx
            @inbounds for k in 1:n
                pos = SVector{3, Float64}(sl.X[k][kite_idx], sl.Y[k][kite_idx], sl.Z[k][kite_idx])
                vel = SVector{3, Float64}(sl.vel_kite[k])
                
                if norm(pos) > 1e-9 && norm(vel) > 1e-9
                    radial = pos / norm(pos)
                    tang_vel = vel - (vel ⋅ radial) * radial
                    if norm(tang_vel) > 1e-9
                        tang_vel /= norm(tang_vel)
                        up_z = radial
                        up_y_raw = SVector(-pos[2], pos[1], 0.0)
                        if norm(up_y_raw) > 1e-9
                            up_y = up_y_raw / norm(up_y_raw)
                            up_x = up_z × up_y
                            if norm(up_x) > 1e-9
                                up_x /= norm(up_x)
                                up_y = up_z × up_x  # re-orthonormalize
                                R_up = @SMatrix [up_x[1] up_y[1] up_z[1];
                                                 up_x[2] up_y[2] up_z[2];
                                                 up_x[3] up_y[3] up_z[3]]
                                heading_vec = R_up' * tang_vel
                                yaw_sphere[k] = atan(heading_vec[2], heading_vec[1])
                                continue
                            end
                        end
                    end
                end
                yaw_sphere[k] = k > 1 ? yaw_sphere[k-1] : NaN
            end

            yaw_sphere_unwrapped = copy(yaw_sphere)
            unwrap_phase!(yaw_sphere_unwrapped)
            yaw_rate_sphere = Base.diff(rad2deg.(yaw_sphere_unwrapped)) ./ Base.diff(sl.time)

            # Averages over the last 33 seconds for quick inspection
            t_rate_diff = sl.time[1:end-1]
            window_start = sl.time[end] - 10
            mask_diff = t_rate_diff .>= window_start
            mask_diff = any(mask_diff) ? mask_diff : trues(length(t_rate_diff))
            mask_full = sl.time .>= window_start
            mask_full = any(mask_full) ? mask_full : trues(length(sl.time))
            avg_heading = mean(yaw_rate_perp_wind[mask_diff])
            avg_ekf = yaw_rate_EKF === nothing ? NaN : mean(yaw_rate_EKF[mask_full])
            avg_sphere = mean(yaw_rate_sphere[mask_diff])
            @info "yaw_rate_paper averages (last 33 s)" suffix avg_heading avg_ekf avg_sphere

            push_trace!(data, label, time) = begin
                push!(all_data, data); push!(all_labels, label); push!(all_times, time)
            end
            push_trace!(yaw_rate_perp_wind, "dψ/dt (heading⟂wind)" * suffix, sl.time[1:end-1])
            if yaw_rate_EKF !== nothing
                push_trace!(yaw_rate_EKF, "dψ/dt (EKF yaw)" * suffix, sl.time)
                # push_trace!(yaw_ekf, "ψ (EKF yaw)" * suffix, sl.time)
            end
            # push_trace!(yaw_rate_sphere, "dψ/dt (HeadingGate)" * suffix, sl.time[1:end-1])

            if yaw_rate_paper_compare
                yaw_cmp = yaw_ekf === nothing ? copy(heading_unwrapped) : copy(yaw_ekf)
                yaw_sph = copy(yaw_sphere)
                for k in 2:n
                    if !isnan(yaw_cmp[k]) && !isnan(yaw_cmp[k - 1])
                        dψ = yaw_cmp[k] - yaw_cmp[k - 1]
                        if dψ > π
                            yaw_cmp[k] -= 2π
                        elseif dψ < -π
                            yaw_cmp[k] += 2π
                        end
                    end
                    dψ = yaw_sph[k] - yaw_sph[k - 1]
                    if dψ > π
                        yaw_sph[k] -= 2π
                    elseif dψ < -π
                        yaw_sph[k] += 2π
                    end
                end
                mask = .!isnan.(yaw_cmp) .& .!isnan.(yaw_sph)
                yaw_cmp_f = yaw_cmp[mask]
                yaw_sph_f = yaw_sph[mask]
                t_f = sl.time[mask]
                if isempty(t_f)
                    push!(compare_results, (suffix=suffix,
                                            std_diff=NaN,
                                            std_diff_deg=NaN,
                                            std_rate=NaN,
                                            std_rate_deg=NaN))
                else
                    diff_yaw = yaw_cmp_f .- yaw_sph_f
                    diff_yaw .-= mean(diff_yaw)
                    std_diff = std(diff_yaw)
                    rate_cmp = Base.diff(yaw_cmp_f) ./ Base.diff(t_f)
                    rate_sph = Base.diff(yaw_sph_f) ./ Base.diff(t_f)
                    std_rate = std(rate_cmp .- rate_sph)
                    push!(compare_results, (suffix=suffix,
                                            std_diff=std_diff,
                                            std_diff_deg=rad2deg(std_diff),
                                            std_rate=std_rate,
                                            std_rate_deg=rad2deg(std_rate)))
                end
            end
        end
        if !isempty(all_data)
            push!(panels, (
                data = all_data,
                labels = all_labels,
                times = all_times,
                ylabel = "dψ/dt paper\n",
                ylims = yaw_rate_paper_ylims
            ))
            if yaw_rate_paper_compare && !isempty(compare_results)
                for res in compare_results
                    @info "yaw_rate_paper comparison (EKF vs sphere)" res...
                end
            end
        end
    end

    if plot_gk_paper
        all_data = []
        all_labels = []
        all_times = []
        for (i, lg) in enumerate(logs)
            sl = lg.syslog
            suffix = " - " * syss[i].name
            n = length(sl.time)
            if isempty(sl.vel_kite) || isempty(sl.v_wind_kite)
                @warn "Missing vel_kite or v_wind_kite in syslog; skipping gk_paper for $suffix"
                continue
            end
            # Paper yaw-rate from apparent wind
            yaw = Vector{Float64}(undef, n)
            @inbounds for k in 1:n
                v = sl.vel_kite[k]
                w = sl.v_wind_kite[k]
                va_enu = w .- v
                va_ned = SVector{3, Float64}(va_enu[2], va_enu[1], -va_enu[3])
                yaw[k] = atan(va_ned[2], va_ned[1])
            end
            for k in 2:n
                dψ = yaw[k] - yaw[k - 1]
                if dψ > π
                    yaw[k] -= 2π
                elseif dψ < -π
                    yaw[k] += 2π
                end
            end
            yaw_rate = diff(rad2deg.(yaw)) ./ diff(sl.time)

            # Reconstruct steering command from segment 87 (same calibration as controller)
            seg_left = syss[i].segments[87]
            p_i, p_j = seg_left.point_idxs
            xs = sl.X; ys = sl.Y; zs = sl.Z
            steering_len = similar(sl.time)
            @inbounds for k in eachindex(sl.time)
                p1 = SVector{3, Float64}(xs[k][p_i], ys[k][p_i], zs[k][p_i])
                p2 = SVector{3, Float64}(xs[k][p_j], ys[k][p_j], zs[k][p_j])
                steering_len[k] = norm(p2 - p1)
            end
            steering_l0 = steering_len[1]
            us_cmd = similar(steering_len)
            @inbounds for k in eachindex(us_cmd)
                δ = steering_len[k] - steering_l0
                us_cmd[k] = abs(δ) > 1e-6 ? δ / 1.4 : 0.0
            end
            us_seg = us_cmd[2:end]  # align with diff-based yaw_rate

            v_app = sl.v_app[2:end]  # apparent wind magnitude

            gk = similar(yaw_rate)
            @inbounds for k in eachindex(gk)
                gk[k] = abs(us_seg[k]) > 1e-8 ? yaw_rate[k] / (v_app[k] * us_seg[k]) : NaN
            end

            push!(all_data, gk)
            push!(all_labels, "gk_paper" * suffix)
            push!(all_times, sl.time[2:end])
        end
        if !isempty(all_data)
            push!(panels, (
                data = all_data,
                labels = all_labels,
                times = all_times,
                ylabel = "gk (paper ψ̇)\n[s/m]"
            ))
        end
    end

    if plot_cs
        all_data = []
        all_labels = []
        all_times = []
        for (i, lg) in enumerate(logs)
            suffix = actual_suffixes[i]
            cs_vals, cs_time = calculate_cs(lg, syss[i])
            if isempty(cs_vals)
                @warn "Missing data for C_S; skipping plot_cs for $suffix"
                continue
            end
            push!(all_data, cs_vals)
            push!(all_labels, "C_S" * suffix)
            push!(all_times, cs_time)
        end
        if !isempty(all_data)
            push!(panels, (
                data = all_data,
                labels = all_labels,
                times = all_times,
                ylabel = "C_S [-]\nside force",
                ylims = cs_ylims
            ))
        end
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
                # Fallback for REFINE logs: reconstruct ω_v from quaternions
                kite_idx = syss[i].wings[1].origin_idx
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
                    pos_w = SVector{3, Float64}(sl.X[k][kite_idx], sl.Y[k][kite_idx], sl.Z[k][kite_idx])
                    e_x = SVector{3, Float64}(R1[:, 1])
                    R_v_w = SymbolicAWEModels.calc_R_v_w(pos_w, e_x)
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
                    push!(all_labels, "v_ro,$j" * suffix)
                    push!(all_times, sl.time)
                end
            end
        end
        if !isempty(all_data)
            push!(panels, (
                data = all_data,
                labels = all_labels,
                times = all_times,
                ylabel = "reel-out speed [m/s]"
            ))
        end
    end

    if plot_tether
        all_data = []
        all_labels = []
        all_times = []
        for (i, lg) in enumerate(logs)
            sl = lg.syslog
            suffix = actual_suffixes[i]
            # Extract tether length (first element if vector, otherwise scalar)
            l_tether = [length(sl.l_tether[j]) > 0 ? sl.l_tether[j][1] : 0.0
                        for j in eachindex(sl.l_tether)]
            push!(all_data, l_tether)
            push!(all_labels, "l_tether" * suffix)
            push!(all_times, sl.time)

            # Add setpoint if provided
            if !isnothing(tether_len_setpoint)
                is_multi_setpoint = (tether_len_setpoint isa Vector &&
                                     length(tether_len_setpoint) > 0 &&
                                     tether_len_setpoint[1] isa AbstractVector)

                if is_multi_setpoint
                    if i <= length(tether_len_setpoint) && !isnothing(tether_len_setpoint[i])
                        push!(all_data, tether_len_setpoint[i])
                        push!(all_labels, "l_sp" * suffix)
                        push!(all_times, sl.time)
                    end
                else
                    push!(all_data, tether_len_setpoint)
                    push!(all_labels, "l_sp" * suffix)
                    push!(all_times, sl.time)
                end
            end
        end
        push!(panels, (
            data = all_data,
            labels = all_labels,
            times = all_times,
            ylabel = "winched tether length [m]"
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

    if !isnothing(tape_lengths)
        all_data = []
        all_labels = []
        all_times = []
        for (i, tl) in enumerate(tape_lengths)
            suffix = actual_suffixes[i]
            if hasproperty(tl, :right_steering) && !isempty(tl.right_steering)
                push!(all_data, tl.right_steering)
                push!(all_labels, "L_right" * suffix)
                push!(all_times, tl.time)
            end
            if hasproperty(tl, :depower) && !isempty(tl.depower)
                push!(all_data, tl.depower)
                push!(all_labels, "L_depower" * suffix)
                push!(all_times, tl.time)
            end
        end
        if !isempty(all_data)
            push!(panels, (
                data = all_data,
                labels = all_labels,
                times = all_times,
                ylabel = "tape length [m]"
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
            push!(all_labels, "F_aero,z" * suffix)
            push!(all_times, sl.time)
        end
        push!(panels, (
            data = all_data,
            labels = all_labels,
            times = all_times,
            ylabel = "F_aero,z [N]"
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
            push!(all_labels, "M_aero,z" * suffix)
            push!(all_times, sl.time)
        end
        push!(panels, (
            data = all_data,
            labels = all_labels,
            times = all_times,
            ylabel = "M_aero,z [Nm]"
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
            push!(all_labels, "M_tether,y" * suffix)
            push!(all_times, sl.time)
        end
        push!(panels, (
            data = all_data,
            labels = all_labels,
            times = all_times,
            ylabel = "M_tether,y [Nm]"
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
                    push!(all_labels, "β_$j" * suffix)
                    push!(all_times, sl.time)
                end
            end
        end
        if !isempty(all_data)
            push!(panels, (
                data = all_data,
                labels = all_labels,
                times = all_times,
                ylabel = "twist [°]"
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
            # Reconstruct steering tape length of segment 87 from logged node positions
            seg_left = syss[i].segments[87]
            p_i, p_j = seg_left.point_idxs
            xs = sl.X
            ys = sl.Y
            zs = sl.Z
            n = length(sl.time)
            steering_len = zeros(Float64, n)
            @inbounds for k in 1:n
                p1 = SVector{3,Float64}(xs[k][p_i], ys[k][p_i], zs[k][p_i])
                p2 = SVector{3,Float64}(xs[k][p_j], ys[k][p_j], zs[k][p_j])
                steering_len[k] = norm(p2 - p1)
            end
            steering_l0 = 1.6  # HARDCODED V3 baseline
            # Convert segment length change to steering command (same calibration as controller)
            us_cmd = similar(steering_len)
            @inbounds for k in eachindex(us_cmd)
                δ = steering_len[k] - steering_l0
                us_cmd[k] = abs(δ) > 1e-6 ? δ / 1.4 : 0.0
            end
            us_seg = us_cmd[2:end]

            push!(all_data, us_seg)
            push!(all_labels, "us_cmd" * suffix)
            push!(all_times, sl.time[2:end])  # align with us_seg length
        end
        push!(panels, (
            data = all_data,
            labels = all_labels,
            times = all_times,
            ylabel = "u_s (-)", # \nsteering"
        ))
    end

    if plot_gk
        all_data = []
        all_labels = []
        all_times = []
        for (i, lg) in enumerate(logs)
            sl = lg.syslog
            suffix = actual_suffixes[i]

            # Reconstruct steering tape length of segment 87 from logged node positions
            seg_left = syss[i].segments[87]
            p_i, p_j = seg_left.point_idxs
            xs = sl.X
            ys = sl.Y
            zs = sl.Z
            n = length(sl.time)
            steering_len = zeros(Float64, n)
            @inbounds for k in 1:n
                p1 = SVector{3,Float64}(xs[k][p_i], ys[k][p_i], zs[k][p_i])
                p2 = SVector{3,Float64}(xs[k][p_j], ys[k][p_j], zs[k][p_j])
                steering_len[k] = norm(p2 - p1)
            end
            steering_l0 = 1.6  #TODO: hardcoded V3 baseline
            # Convert segment length change to steering command (same calibration as controller)
            us = similar(steering_len)
            @inbounds for k in eachindex(us)
                δ = steering_len[k] - steering_l0
                us[k] = abs(δ) > 1e-6 ? δ / 1.4 : 0.0
            end

            # Calculate heading rate from diff for quaternion wings
            heading_unwrapped = copy(sl.heading)
            for j in 2:length(heading_unwrapped)
                while heading_unwrapped[j] - heading_unwrapped[j-1] > π
                    heading_unwrapped[j] -= 2π
                end
                while heading_unwrapped[j] - heading_unwrapped[j-1] < -π
                    heading_unwrapped[j] += 2π
                end
            end
            heading_rate = diff(rad2deg.(heading_unwrapped)) ./ diff(sl.time)
            v_app = sl.v_app[2:end]
            us_seg = us[2:end]

            # calculate gk, guarding against zero steering
            gk = similar(heading_rate)
            @inbounds for k in eachindex(gk)
                gk[k] = abs(us_seg[k]) > 1e-8 ? heading_rate[k] / (v_app[k] * us_seg[k]) : NaN
            end
            
            # force ylimits from 0 to 10 for gk axis
            cs_vals = Float64[]
            cs_over_us_vec = Float64[]
            try
                cs_vals, _ = calculate_cs(lg, syss[i])
                n_us = length(us_seg)
                n_cs = length(cs_vals)
                if n_us > 0 && n_cs >= n_us + 1
                    cs_over_us_vec = similar(us_seg)
                    @inbounds for k in eachindex(us_seg)
                        cs_over_us_vec[k] = abs(us_seg[k]) > 1e-8 ? cs_vals[k + 1] / us_seg[k] : NaN
                    end
                end
            catch err
                @warn "Failed to compute cs/us" exception=(err, catch_backtrace())
            end

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
            # # ss.AoA =wing.vsm_aero.alpha_uncorrected[length(wing.vsm_solver.sol.alpha_dist) ÷ 2 + (length(wing.vsm_solver.sol.alpha_dist) % 2)] # version-3, hopefullu without induction
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
            push!(all_labels, "gk" * suffix)
            push!(all_times, sl.time[2:end])
        end
            push!(panels, (
                data = all_data,
                labels = all_labels,
                times = all_times,
                ylabel = "gk [°/m]\n(dψ/dt)/va·us",
                ylims = gk_ylims
            ))
        end

    if plot_v_app
        all_data = []
        all_labels = []
        all_times = []
        for (i, lg) in enumerate(logs)
            sl = lg.syslog
            suffix = actual_suffixes[i]
            v_app = sl.v_app
            push!(all_data, v_app)
            push!(all_labels, "v_app" * suffix)
            push!(all_times, sl.time)
        end
        push!(panels, (
            data = all_data,
            labels = all_labels,
            times = all_times,
            ylabel = "v_a (ms^-1)",
        ))
    end

    if plot_kite_vel
        all_data = []
        all_labels = []
        all_times = []
        for (i, lg) in enumerate(logs)
            sl = lg.syslog
            suffix = actual_suffixes[i]
            v_kite_norm = [norm(v) for v in sl.vel_kite]
            push!(all_data, v_kite_norm)
            push!(all_labels, "|v_kite|" * suffix)
            push!(all_times, sl.time)
        end
        push!(panels, (
            data = all_data,
            labels = all_labels,
            times = all_times,
            ylabel = "kite velocity [m/s]"
        ))
    end

    if plot_aoa
        all_data = []
        all_labels = []
        all_times = []
        for (i, lg) in enumerate(logs)
            sl = lg.syslog
            suffix = actual_suffixes[i]
            aoa_deg = rad2deg.(sl.AoA)
            push!(all_data, aoa_deg)
            push!(all_labels, "α" * suffix)
            push!(all_times, sl.time)
        end
        push!(panels, (
            data = all_data,
            labels = all_labels,
            times = all_times,
            ylabel = "α_w (°)",#\nangle of attack",
            ylims = aoa_ylims
        ))
    end

    if plot_heading
        all_data = []
        all_labels = []
        all_times = []
        for (i, lg) in enumerate(logs)
            sl = lg.syslog
            suffix = actual_suffixes[i]
            heading_deg = rad2deg.(sl.heading)
            course_deg = rad2deg.(sl.course)
            push!(all_data, heading_deg)
            push!(all_labels, "ψ" * suffix)
            push!(all_times, sl.time)
            push!(all_data, course_deg)
            push!(all_labels, "χ" * suffix)
            push!(all_times, sl.time)
            if plot_kiteutils_course
                course_kiteutils_deg = [rad2deg(KiteUtils.calc_course(sl.orient[i])) for i in eachindex(sl.orient)]
                push!(all_data, course_kiteutils_deg)
                push!(all_labels, "χ_KU" * suffix)
                push!(all_times, sl.time)
            end
            # Add setpoint if provided
            if !isnothing(heading_setpoint)
                # Check if heading_setpoint is a vector of vectors (multiple logs)
                # or a vector of numbers (single log)
                is_multi_setpoint = (heading_setpoint isa Vector &&
                                     length(heading_setpoint) > 0 &&
                                     heading_setpoint[1] isa AbstractVector)

                if is_multi_setpoint
                    # Multiple setpoints (one per log)
                    if i <= length(heading_setpoint) && !isnothing(heading_setpoint[i])
                        setpoint_deg = rad2deg.(heading_setpoint[i])
                        push!(all_data, setpoint_deg)
                        push!(all_labels, "ψ_sp" * suffix)
                        push!(all_times, sl.time)
                    end
                else
                    # Single setpoint for all logs (or single log)
                    setpoint_deg = rad2deg.(heading_setpoint)
                    push!(all_data, setpoint_deg)
                    push!(all_labels, "ψ_sp" * suffix)
                    push!(all_times, sl.time)
                end
            end
        end
        push!(panels, (
            data = all_data,
            labels = all_labels,
            times = all_times,
            ylabel = "heading/course [°]"
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
            push!(all_labels, "ψ_old" * suffix)
            push!(all_times, sl.time)
        end
        push!(panels, (
            data = all_data,
            labels = all_labels,
            times = all_times,
            ylabel = "old heading [°]"
        ))
    end

    if plot_distance
        all_data = []
        all_labels = []
        all_times = []
        for (i, lg) in enumerate(logs)
            sl = lg.syslog
            suffix = actual_suffixes[i]
            # Get wing origin index
            kite_idx = syss[i].wings[1].origin_idx
            distance = [norm([sl.X[j][kite_idx], sl.Y[j][kite_idx], sl.Z[j][kite_idx]]) for j in eachindex(sl.X)]
            push!(all_data, distance)
            push!(all_labels, "r" * suffix)
            push!(all_times, sl.time)
        end
        push!(panels, (
            data = all_data,
            labels = all_labels,
            times = all_times,
            ylabel = "distance [m]"
        ))
    end

    if plot_cone_angle
        all_data = []
        all_labels = []
        all_times = []
        for (i, lg) in enumerate(logs)
            sl = lg.syslog
            suffix = actual_suffixes[i]
            # Get wing origin index
            kite_idx = syss[i].wings[1].origin_idx
            # Assuming wind_vec_gnd is available in syslog
            cone_angle_rad = similar(sl.heading)
            for j in eachindex(sl.X)
                pos = [sl.X[j][kite_idx], sl.Y[j][kite_idx], sl.Z[j][kite_idx]]
                pos_norm = normalize(pos)
                wind_norm = normalize([sl.v_wind_gnd[1], sl.v_wind_gnd[2], sl.v_wind_gnd[3]])
                cone_angle_rad[j] = acos(clamp(dot(pos_norm, wind_norm), -1.0, 1.0))
            end
            cone_angle_deg = rad2deg.(cone_angle_rad)
            push!(all_data, cone_angle_deg)
            push!(all_labels, "θ_cone" * suffix)
            push!(all_times, sl.time)
        end
        push!(panels, (
            data = all_data,
            labels = all_labels,
            times = all_times,
            ylabel = "cone angle [°]"
        ))
    end

    if plot_elevation
        all_data = []
        all_labels = []
        all_times = []
        for (i, lg) in enumerate(logs)
            sl = lg.syslog
            suffix = actual_suffixes[i]
            elevation_deg = rad2deg.(sl.elevation)
            push!(all_data, elevation_deg)
            push!(all_labels, "elevation" * suffix)
            push!(all_times, sl.time)
        end
        push!(panels, (
            data = all_data,
            labels = all_labels,
            times = all_times,
            ylabel = "β (°)",#\nelevation"
        ))
    end

    if plot_azimuth
        all_data = []
        all_labels = []
        all_times = []
        for (i, lg) in enumerate(logs)
            sl = lg.syslog
            suffix = actual_suffixes[i]
            azimuth_deg = rad2deg.(sl.azimuth)
            push!(all_data, azimuth_deg)
            push!(all_labels, "azimuth" * suffix)
            push!(all_times, sl.time)
        end
        push!(panels, (
            data = all_data,
            labels = all_labels,
            times = all_times,
            ylabel = "ɸ (°)",#\nazimuth"
        ))
    end

    if plot_winch_force
        all_data = []
        all_labels = []
        all_times = []
        for (i, lg) in enumerate(logs)
            sl = lg.syslog
            suffix = actual_suffixes[i]
            winch_force = [[sl.winch_force[i][j] for i in eachindex(sl.winch_force)] for j in 1:3]
            for j in 1:3
                # Only plot if non-zero or if it's index 1
                if j == 1 || !all(iszero, winch_force[j])
                    push!(all_data, winch_force[j])
                    push!(all_labels, "F_winch,$j" * suffix)
                    push!(all_times, sl.time)
                end
            end
        end
        if !isempty(all_data)
            push!(panels, (
                data = all_data,
                labels = all_labels,
                times = all_times,
                ylabel = "F_t ground (N)",#Winch force [N]"
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
                    push!(all_labels, "T_winch,$j" * suffix)
                    push!(all_times, sl.time)
                end
            end
        end
        if !isempty(all_data)
            push!(panels, (
                data = all_data,
                labels = all_labels,
                times = all_times,
                ylabel = "Set torque [Nm]"
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
    for (i, panel) in enumerate(panels)
        # Share x-axis with first subplot
        if i == 1
            ax = Axis(fig[i, 1], ylabel=panel.ylabel)
        else
            ax = Axis(fig[i, 1], ylabel=panel.ylabel, xticklabelsvisible=false)
            linkxaxes!(axes[1], ax)
        end

        # Plot each data series in this panel
        for (j, (data_series, label, time_vec)) in enumerate(zip(panel.data, panel.labels, panel.times))
            if length(data_series) != length(time_vec)
                @warn "Skipping trace '$label': data length $(length(data_series)) != time length $(length(time_vec))"
                continue
            end
            lines!(ax, time_vec, data_series, label=label)
        end

        # Add legend if multiple traces
        if length(panel.data) > 1
            axislegend(ax, position=:lb)
        end
        if hasproperty(panel, :ylims) && !isnothing(panel.ylims)
            ylims!(ax, panel.ylims...)
        end

        push!(axes, ax)
    end

    # Add x-label to bottom subplot
    axes[end].xlabel = "time [s]"
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
    bbox = data_limits(plots)
    center = bbox.origin .+ bbox.widths ./ 2

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
    R_b_w = wing.R_b_w

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

    # Set camera position and lookat
    update_cam!(scene, Vec3f(cam_pos), Vec3f(kite_pos))

    # Set up vector to align with body z-axis AFTER updating camera
    # This ensures the kite's z-axis always points straight up on screen
    cam.upvector[] = Vec3f(R_b_w[:, 3])

    return distance
end

function _plot_with_panes(sys::SystemStructure;
                    size = (1200, 800),
                    margin = 10.0,
                    relmargin = 0.2,
                    segment_color = :black,
                    highlight_color = :red,
                    force_color = false,
                    body_frame = false,
                    perspective = false,
                    kwargs...)
    # Use LScene for advanced camera controls
    scene = Scene(; camera=cam3d!, show_axis = false, size, zoommode = :free, samples = 16)
    plots = plot!(scene, sys; segment_color, margin, force_color, kwargs...)
    
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

    # --- Event Handling for Segments ---
    if haskey(plots, :segments)
        lineseg_plot = plots[:segments]
        seg_colors_obs = lineseg_plot.color
        last_hovered_idx = Ref(-1)
        zoomed_in = Ref(false)

        # --- Hover Labels ---
        # Segment index label at middle of segment
        segment_label = Observable("")
        segment_label_pos = Observable(Point2f(0, 0))
        segment_label_visible = Observable(false)
        text!(scene, segment_label, position = segment_label_pos, space = :pixel,
              fontsize = 14, color = :white, strokecolor = :black, strokewidth = 1,
              align = (:center, :center), visible = segment_label_visible, transparency = true)

        # Point index labels at segment endpoints
        point1_label = Observable("")
        point1_label_pos = Observable(Point2f(0, 0))
        point1_label_visible = Observable(false)
        text!(scene, point1_label, position = point1_label_pos, space = :pixel,
              fontsize = 14, color = :white, strokecolor = :black, strokewidth = 1,
              align = (:center, :center), visible = point1_label_visible, transparency = true)

        point2_label = Observable("")
        point2_label_pos = Observable(Point2f(0, 0))
        point2_label_visible = Observable(false)
        text!(scene, point2_label, position = point2_label_pos, space = :pixel,
              fontsize = 14, color = :white, strokecolor = :black, strokewidth = 1,
              align = (:center, :center), visible = point2_label_visible, transparency = true)

        # --- Event Handler for Robust Hover Highlighting ---
        on(events(scene).mouseposition, priority = 2) do mp
            # This approach is more robust than `pick` for thin lines.
            # It finds the closest segment in 2D screen space.
            min_dist = Inf
            closest_seg_idx = -1
            mouse_pos_2d = Point2f(mp)
            margin_px = 30.0 # pixel margin for hover detection
            if haskey(plots, :segments)
                seg_points_3d = plots[:segments][1][]
                for i in 1:length(sys.segments)
                    p1_3d = seg_points_3d[2 * i - 1]
                    p2_3d = seg_points_3d[2 * i]
                    p1_2d = Makie.project(scene, p1_3d)
                    p2_2d = Makie.project(scene, p2_3d)
                    dist = point_line_segment_distance(mouse_pos_2d, p1_2d, p2_2d)
                    if dist < min_dist
                        min_dist = dist
                        closest_seg_idx = i
                    end
                end
            end

            hover_idx = (min_dist < margin_px) ? closest_seg_idx : -1
            if hover_idx != last_hovered_idx[]
                num_segments = length(sys.segments)

                # Calculate base colors (either force-based or uniform)
                if force_color
                    new_colors = calculate_segment_force_colors(sys.segments, segment_color)
                else
                    new_colors = fill(to_color(segment_color), num_segments)
                end

                if hover_idx != -1
                    seg = sys.segments[hover_idx]
                    p1_3d = sys.points[seg.point_idxs[1]].pos_w
                    p2_3d = sys.points[seg.point_idxs[2]].pos_w

                    # Show segment index 20px to the right of middle of segment
                    mid_point_3d = (p1_3d + p2_3d) / 2
                    mid_point_2d = Makie.project(scene, mid_point_3d)
                    segment_label[] = string(hover_idx)
                    segment_label_pos[] = mid_point_2d + Point2f(20, 0)
                    segment_label_visible[] = true

                    # Show point indices 20px to the right of segment endpoints
                    p1_2d = Makie.project(scene, p1_3d)
                    p2_2d = Makie.project(scene, p2_3d)
                    point1 = sys.points[seg.point_idxs[1]].idx
                    point2 = sys.points[seg.point_idxs[2]].idx
                    point1_label[] = string(point1)
                    point1_label_pos[] = p1_2d + Point2f(20, 0)
                    point1_label_visible[] = true
                    point2_label[] = string(point2)
                    point2_label_pos[] = p2_2d + Point2f(20, 0)
                    point2_label_visible[] = true

                    new_colors[hover_idx] = to_color(highlight_color)
                else
                    segment_label_visible[] = false
                    point1_label_visible[] = false
                    point2_label_visible[] = false
                end
                seg_colors_obs[] = new_colors
                last_hovered_idx[] = hover_idx
            end
        end

        # --- Event Handler for Camera Movement ---
        on(scene.camera.view, priority = 1) do _
            # Update label positions when camera moves
            hover_idx = last_hovered_idx[]
            if hover_idx != -1
                seg = sys.segments[hover_idx]
                p1_3d = sys.points[seg.point_idxs[1]].pos_w
                p2_3d = sys.points[seg.point_idxs[2]].pos_w

                # Update segment label position
                mid_point_3d = (p1_3d + p2_3d) / 2
                mid_point_2d = Makie.project(scene, mid_point_3d)
                segment_label_pos[] = mid_point_2d + Point2f(20, 0)

                # Update point label positions
                p1_2d = Makie.project(scene, p1_3d)
                p2_2d = Makie.project(scene, p2_3d)
                point1_label_pos[] = p1_2d + Point2f(20, 0)
                point2_label_pos[] = p2_2d + Point2f(20, 0)
            end
        end

        # --- Event Handler for Click-to-Zoom ---
        zoomed_in = Ref(false)
        on(events(scene).mousebutton, priority = 1) do event
            if event.button == Mouse.left && event.action == Mouse.press
                cam = scene.camera
                if !zoomed_in[] || last_hovered_idx[] != -1
                    # --- ZOOM IN --- (This part remains the same)
                    hover_idx = last_hovered_idx[]
                    if hover_idx != -1
                        seg = sys.segments[hover_idx]
                        p1_w = sys.points[seg.point_idxs[1]].pos_w
                        p2_w = sys.points[seg.point_idxs[2]].pos_w
                        
                        center = (p1_w + p2_w) / 2.0f0
                        segment_len = norm(p2_w - p1_w)
                        dist_heuristic = segment_len * 1.5 + 2.0
                        
                        inv_view_matrix = inv(cam.view[])
                        cam_dir_vec = normalize(Vec3f(inv_view_matrix[1, 3], inv_view_matrix[2, 3], inv_view_matrix[3, 3]))
                        new_eyepos = center + dist_heuristic * cam_dir_vec
                        
                        update_cam!(scene, new_eyepos, center)
                        zoomed_in[] = true
                        PLOT_ZOOMED_IN[] = true  # Track global zoom state
                        PLOT_ZOOM_SEGMENT_IDX[] = hover_idx  # Track which segment we're zoomed into

                        # Update label positions after zoom
                        # Small delay to ensure camera update is complete
                        sleep(0.01)
                        p1_3d = sys.points[seg.point_idxs[1]].pos_w
                        p2_3d = sys.points[seg.point_idxs[2]].pos_w
                        
                        # Update segment label position
                        mid_point_3d = (p1_3d + p2_3d) / 2
                        mid_point_2d = Makie.project(scene, mid_point_3d)
                        segment_label_pos[] = mid_point_2d + Point2f(20, 0)
                        
                        # Update point label positions
                        p1_2d = Makie.project(scene, p1_3d)
                        p2_2d = Makie.project(scene, p2_3d)
                        point1_label_pos[] = p1_2d + Point2f(20, 0)
                        point2_label_pos[] = p2_2d + Point2f(20, 0)
                        
                        return Consume(true) # Consume the event
                    end
                else
                    zoom_out!(scene, cam, relevant_plots, nothing; relmargin)
                    zoomed_in[] = false
                    PLOT_ZOOMED_IN[] = false  # Track global zoom state
                    PLOT_ZOOM_SEGMENT_IDX[] = -1  # Clear zoomed segment

                    # Update label positions after zoom out
                    # Small delay to ensure camera update is complete
                    sleep(0.01)
                    hover_idx = last_hovered_idx[]
                    if hover_idx != -1
                        seg = sys.segments[hover_idx]
                        p1_3d = sys.points[seg.point_idxs[1]].pos_w
                        p2_3d = sys.points[seg.point_idxs[2]].pos_w
                        
                        # Update segment label position
                        mid_point_3d = (p1_3d + p2_3d) / 2
                        mid_point_2d = Makie.project(scene, mid_point_3d)
                        segment_label_pos[] = mid_point_2d + Point2f(20, 0)
                        
                        # Update point label positions
                        p1_2d = Makie.project(scene, p1_3d)
                        p2_2d = Makie.project(scene, p2_3d)
                        point1_label_pos[] = p1_2d + Point2f(20, 0)
                        point2_label_pos[] = p2_2d + Point2f(20, 0)
                    end
                    
                    return Consume(true) # Consume the event
                end
            end
            return Consume(false)
        end
    end

    # Extract pane observables from plots (created by plot!())
    pane_observables = haskey(plots, :pane_observables) ? plots[:pane_observables] : nothing

    # Set camera projection type
    cam = scene.camera
    if !perspective
        # Use orthographic projection
        # Get current projection view bounds
        widths = scene.viewport[].widths
        w_half = Float32(widths[1] / 2)
        h_half = Float32(widths[2] / 2)
        cam.projection[] = Makie.orthographicprojection(
            -w_half, w_half,
            -h_half, h_half,
            -10_000f0, 10_000f0
        )
    end

    # Set initial camera position
    if body_frame
        zoom_body_frame!(scene, cam, sys)
    else
        update_cam!(scene, Vec3f(-100, -100, 100), Vec3f(0, 0, 0))
        zoom_out!(scene, cam, relevant_plots, nothing; relmargin)
    end

    # Return scene along with pane_observables, margin, plots dict, and relevant_plots
    # These will be used by time-based plotting
    return scene, pane_observables, margin, plots, relevant_plots
end

# Public API function - creates scene with observables for dynamic updates
function Makie.plot(sys::SystemStructure;
                    vector_scale=1.0,
                    force_color=false,
                    segment_color=:black,
                    plot_aero=true,
                    relmargin=0.2,
                    body_frame=false,
                    perspective=false,
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
    scene, pane_observables, margin, plots, relevant_plots = _plot_with_panes(sys;
                geometry_obs,
                vector_scale,
                force_color,
                segment_color,
                plot_aero,
                relmargin,
                body_frame,
                perspective,
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
    PLOT_BACKGROUND_PANES[] = pane_observables
    PLOT_MARGIN[] = margin
    PLOT_RELEVANT_PLOTS[] = relevant_plots
    # SystemStructure and settings already stored above before plot creation
    PLOT_ZOOMED_IN[] = false  # Initialize zoom state (not zoomed in)
    PLOT_ZOOM_RELMARGIN[] = relmargin  # Store relmargin for auto-updates
    PLOT_ZOOM_SEGMENT_IDX[] = -1  # No segment zoomed initially

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

    # Update background panes based on current point positions
    if !isnothing(PLOT_BACKGROUND_PANES[])
        panes = PLOT_BACKGROUND_PANES[]
        margin = PLOT_MARGIN[]

        # Recalculate limits based on current point positions
        xlims, ylims, zlims = (-10, 10), (-10, 10), (-10, 10)
        if !isempty(sys.points)
            all_x = [p.pos_w[1] for p in sys.points]
            all_y = [p.pos_w[2] for p in sys.points]
            all_z = [p.pos_w[3] for p in sys.points]

            xlims_data = extrema(all_x)
            ylims_data = extrema(all_y)
            zlims_data = extrema(all_z)

            xlims = (xlims_data[1] - margin, xlims_data[2] + margin)
            ylims = (ylims_data[1] - margin, ylims_data[2] + margin)
            zlims = (zlims_data[1] - margin, zlims_data[2] + margin)
        end

        # Update pane observables with infinite extent
        pane_extent = 10000.0f0

        # XZ plane at y_max - extends far in -X and +Z directions
        panes[1][] = Rect3(Vec3f(xlims[1] - pane_extent, ylims[2], zlims[1]),
                           Vec3f(xlims[2] - xlims[1] + pane_extent, 0.01, pane_extent))
        # YZ plane at x_max - extends far in -Y and +Z directions
        panes[2][] = Rect3(Vec3f(xlims[2], ylims[1] - pane_extent, zlims[1]),
                           Vec3f(0.01, ylims[2] - ylims[1] + pane_extent, pane_extent))
        # XY plane at z_min - extends far in -X and -Y directions
        panes[3][] = Rect3(Vec3f(xlims[1] - pane_extent, ylims[1] - pane_extent, zlims[1]),
                           Vec3f(xlims[2] - xlims[1] + pane_extent, ylims[2] - ylims[1] + pane_extent, 0.01))
    end

    return nothing
end

"""
    record(lg::SysLog, sys::SystemStructure, filename::String; framerate=30, kwargs...)

Record a SysLog animation to an MP4 video file.

# Arguments
- `lg::SysLog`: The simulation log to record
- `sys::SystemStructure`: The system structure matching the log's topology
- `filename::String`: Output video filename (e.g., "simulation.mp4")

# Keyword Arguments
- `framerate::Int=30`: Video framerate (frames per second)
- `vector_scale::Real=1.0`: Scale factor for wing orientation arrows
- All other keyword arguments are passed through to the SystemStructure plot function

# Returns
- The Scene object used for recording

# Example
```julia
# Record simulation to MP4
record(log, sys_struct, "output.mp4")

# Record with custom framerate and settings
record(log, sys_struct, "simulation.mp4", framerate=60, vector_scale=0.3)
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

    # Create initial plot which sets up observables and scene (following standard Makie pattern)
    scene = plot(sys; vector_scale, show_panes=show_panes, kwargs...)

    # Create pixel-space subscene for UI controls overlay
    ui_scene = Scene(scene, viewport=scene.viewport, clear=false, camera=campixel!)

    # Create observable for current frame index
    frame_idx = Observable(1)

    # Function to update to a specific frame (using plot! to mutate observables)
    function update_frame!(idx)
        ss = lg.syslog[idx]
        update_from_sysstate!(sys, ss)
        plot!(sys; vector_scale)
    end

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
    slider_thumb_x = @lift(slider_x_start + ($(frame_idx) - 1) / (n_frames - 1) * $(slider_width_obs))
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

    # Info label
    info_label_x = body_frame_button_x + button_width + 20
    info_text = @lift("Frame: $($(frame_idx))/$n_frames | Time: $(@sprintf("%.2f", lg.syslog[$(frame_idx)].time)) s | Speed: $(replay_speed)x")
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
                    # Disable zoom-in mode when switching to body frame
                    if PLOT_BODY_FRAME[]
                        PLOT_ZOOMED_IN[] = false
                        PLOT_ZOOM_SEGMENT_IDX[] = -1
                    end
                    # Force camera update immediately (not just during playback)
                    PLOT_CAMERA_DISTANCE[] = nothing  # Force recalculation on mode change
                    if !isnothing(PLOT_SCENE[]) && !isnothing(PLOT_RELEVANT_PLOTS[]) && !isnothing(PLOT_SYSTEM_STRUCTURE[])
                        scene = PLOT_SCENE[]
                        relevant_plots = PLOT_RELEVANT_PLOTS[]
                        stored_sys = PLOT_SYSTEM_STRUCTURE[]

                        if PLOT_BODY_FRAME[]
                            dist = zoom_body_frame!(scene, scene.camera, stored_sys, PLOT_CAMERA_DISTANCE[])
                            PLOT_CAMERA_DISTANCE[] = dist
                        else
                            dist = zoom_out!(scene, scene.camera, relevant_plots, PLOT_CAMERA_DISTANCE[]; relmargin=PLOT_ZOOM_RELMARGIN[])
                            PLOT_CAMERA_DISTANCE[] = dist
                        end
                    end
                    # Update tracking variables to prevent mode change detection
                    PLOT_PREV_BODY_FRAME[] = PLOT_BODY_FRAME[]
                    PLOT_PREV_ZOOMED_IN[] = PLOT_ZOOMED_IN[]
                    PLOT_PREV_SEGMENT_IDX[] = PLOT_ZOOM_SEGMENT_IDX[]
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
                    if frame_idx[] > 1
                        dt = lg.syslog[frame_idx[]].time - lg.syslog[frame_idx[] - 1].time
                        sleep(max(0.01, dt / replay_speed))
                    else
                        sleep(0.05)
                    end
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

    return scene
end

end
