# Copyright (c) 2025 Bart van de Lint, Jelle Poland
# SPDX-License-Identifier: MPL-2.0

"""
Shared utility functions for coupled examples.
"""

using CSV, DataFrames
using LinearAlgebra
using UnPack

# V3 Kite steering/depower calibration (from KCU documentation)
const V3_STEERING_L0 = 1.4    # Neutral steering tape length (m)
const V3_STEERING_GAIN = 1.4  # Maximum differential (m) at |u_s| = 1
const V3_DEPOWER_L0 = 0.0     # Neutral depower tape length (m)
const V3_DEPOWER_GAIN = 5.0   # Depower range (m) for 0-100%

"""
    steering_percentage_to_lengths(percentage; l0=V3_STEERING_L0, gain=V3_STEERING_GAIN)

Convert steering percentage to left/right tape lengths (m).
Percentage convention: negative = left turn, positive = right turn.
Uses half-gain on each side for symmetric actuation.
"""
function steering_percentage_to_lengths(percentage;
                                        l0=V3_STEERING_L0, gain=V3_STEERING_GAIN)
    u_s = percentage / 100.0
    L_left = l0 - (gain / 2.0) * u_s
    L_right = l0 + (gain / 2.0) * u_s
    return L_left, L_right
end

"""
    csv_steering_percentage_to_lengths(percentage; l0=V3_STEERING_L0, gain=V3_STEERING_GAIN)

Convert CSV steering percentage to left/right tape lengths (m).
Uses opposite sign convention and full gain (matches CSV data format).
"""
function csv_steering_percentage_to_lengths(percentage;
                                            l0=V3_STEERING_L0, gain=V3_STEERING_GAIN)
    u_s = percentage / 100.0
    L_left = l0 + gain * u_s
    L_right = l0 - gain * u_s
    return L_left, L_right
end

"""
    depower_percentage_to_length(percentage; l0=V3_DEPOWER_L0, gain=V3_DEPOWER_GAIN)

Convert depower percentage to tape length (m).
"""
function depower_percentage_to_length(percentage;
                                      l0=V3_DEPOWER_L0, gain=V3_DEPOWER_GAIN)
    u_p = percentage / 100.0
    return l0 + gain * u_p
end

"""
    build_geom_suffix(depower_l0, tip_reduction, te_frac)

Build geometry filename suffix from configuration parameters.
"""
function build_geom_suffix(depower_l0, tip_reduction, te_frac)
    return "depower$(depower_l0)_tip$(tip_reduction)_te$(te_frac)"
end

"""
    load_extra_points(csv_path, sys_struct; body_offset)

Load extra points from CSV and transform from camera frame to simulation frame.
Returns (transformed_points, groups) where groups is Vector of (group_name, indices).

CSV has columns: group, idx_in_group, x, y, z.
Alignment: CSV strut3/strut4 LE centers align with sim points 10, 12.
"""
function load_extra_points(csv_path::String, sys_struct; body_offset=[0.3, 0.0, 0.2])
    df = CSV.read(csv_path, DataFrame)

    # CSV strut centers: strut3[1]/strut4[1] are at TE, [end] are at LE
    strut3 = [[r.x, r.y, r.z] for r in eachrow(df) if r.group == "strut3"]
    strut4 = [[r.x, r.y, r.z] for r in eachrow(df) if r.group == "strut4"]
    csv_le_center = (strut3[end] + strut4[end]) / 2
    csv_te_center = (strut3[1] + strut4[1]) / 2

    # Sim reference: points 10, 12 (center LE)
    sim_p10 = collect(sys_struct.points[10].pos_w)
    sim_p12 = collect(sys_struct.points[12].pos_w)
    sim_le_center = (sim_p10 + sim_p12) / 2

    # Direction vectors
    csv_span = normalize(strut4[end] - strut3[end])

    # CSV basis: y=spanwise, z from wing center geometry, x from cross
    csv_y = csv_span
    csv_wing_center = (csv_le_center + csv_te_center) / 2
    csv_z = normalize(csv_wing_center - csv_y * 0.84/2)
    csv_x = cross(csv_y, csv_z)

    # Sim basis: directly from wing rotation matrix
    R_b_to_w = sys_struct.wings[1].R_b_to_w
    sim_x = R_b_to_w[:, 1]
    sim_y = R_b_to_w[:, 2]
    sim_z = R_b_to_w[:, 3]

    # Rotation: R * csv_basis = sim_basis
    csv_basis = hcat(csv_x, csv_y, csv_z)
    sim_basis = hcat(sim_x, sim_y, sim_z)
    R = sim_basis * csv_basis'

    # Translation: align LE centers
    T = sim_le_center - R * csv_le_center + R_b_to_w * body_offset

    # Transform all points (including camera origin marker at zeros)
    all_pts = [[row.x, row.y, row.z] for row in eachrow(df)]
    push!(all_pts, zeros(3))
    transformed = [Tuple(R * p + T) for p in all_pts]

    # Build group indices (1-based)
    groups = Vector{Tuple{String, Vector{Int}}}()
    current_group = ""
    current_indices = Int[]
    for (i, row) in enumerate(eachrow(df))
        if row.group != current_group
            if !isempty(current_indices)
                push!(groups, (current_group, copy(current_indices)))
            end
            current_group = row.group
            current_indices = [i]
        else
            push!(current_indices, i)
        end
    end
    if !isempty(current_indices)
        push!(groups, (current_group, current_indices))
    end

    return transformed, groups
end

const PLOT_COLORS = [:blue, :green, :orange, :purple, :cyan, :magenta]

"""
    plot_body_frame_local(sys_structs; extra_points, extra_groups, dir, labels)

Plot wing points in 2D body frame. Accepts single sys_struct or vector of them.
Extra points connected per-strut and LE.

# Arguments
- `sys_structs`: System structure or vector of system structures to plot
- `extra_points`: Optional vector of (x,y,z) tuples from load_extra_points
- `extra_groups`: Optional groups from load_extra_points
- `dir::Symbol`: Viewing direction (:side, :front, or :top)
- `labels`: Optional vector of labels for each sys_struct
- `point_idxs`: Optional vector of point indices to plot (default: WING points only)
"""
function plot_body_frame_local(sys_structs;
                               extra_points=nothing,
                               extra_groups=nothing,
                               dir::Symbol=:front,
                               point_size=10,
                               extra_point_size=8,
                               figsize=(800, 600),
                               labels=nothing,
                               point_idxs=nothing,
                               legend=true,
                               title=true,
                               show_point_idxs=true)
    # Normalize to vector
    structs = sys_structs isa Vector ? sys_structs : [sys_structs]
    n_structs = length(structs)

    # Default labels
    if isnothing(labels)
        labels = n_structs == 1 ? ["sim"] : ["sim_$i" for i in 1:n_structs]
    end

    # Set up axis labels
    if dir == :top
        xlabel, ylabel = "x [m]", "y [m]"
    elseif dir == :side
        xlabel, ylabel = "x [m]", "z [m]"
    else  # :front
        xlabel, ylabel = "y [m]", "z [m]"
    end

    fig = Figure(size=figsize)
    ax_title = title ? "Wing Points (Body Frame)" : ""
    ax = Axis(fig[1, 1]; xlabel, ylabel,
              title=ax_title, aspect=DataAspect())

    function get_2d(pos_b)
        if dir == :top
            return (pos_b[1], pos_b[2])
        elseif dir == :side
            return (pos_b[1], pos_b[3])
        else
            return (pos_b[2], pos_b[3])
        end
    end

    # Collect all coordinates for auto-zoom
    all_x_vals = Float64[]
    all_y_vals = Float64[]

    # Plot each sys_struct
    for (s_idx, sys_struct) in enumerate(structs)
        @unpack points, wings, segments = sys_struct
        color = PLOT_COLORS[mod1(s_idx, length(PLOT_COLORS))]

        # Update pos_b for REFINE wing points
        for wing in wings
            if wing.wing_type == SymbolicAWEModels.REFINE
                R_w_b = wing.R_b_to_w'
                for point in points
                    if point.wing_idx == wing.idx
                        point.pos_b .= R_w_b * (point.pos_w - wing.pos_w)
                    end
                end
            end
        end

        # Select points to plot: use point_idxs if provided, otherwise WING points
        if isnothing(point_idxs)
            plot_points = [p for p in points if p.type == SymbolicAWEModels.WING]
        else
            plot_points = [points[i] for i in point_idxs if i <= length(points)]
        end

        # Extract 2D coords based on viewing direction
        if dir == :top
            coords = [(p.pos_b[1], p.pos_b[2]) for p in plot_points]
        elseif dir == :side
            coords = [(p.pos_b[1], p.pos_b[3]) for p in plot_points]
        else  # :front
            coords = [(p.pos_b[2], p.pos_b[3]) for p in plot_points]
        end

        x_vals = [c[1] for c in coords]
        y_vals = [c[2] for c in coords]
        append!(all_x_vals, x_vals)
        append!(all_y_vals, y_vals)

        # Plot segments (skip diagonals 29-46)
        plot_point_idxs = Set(p.idx for p in plot_points)
        for seg in segments
            if 29 <= seg.idx <= 46
                continue
            end
            from_idx, to_idx = seg.point_idxs
            if from_idx in plot_point_idxs && to_idx in plot_point_idxs
                p1 = points[from_idx]
                p2 = points[to_idx]
                c1 = get_2d(p1.pos_b)
                c2 = get_2d(p2.pos_b)
                lines!(ax, [c1[1], c2[1]], [c1[2], c2[2]];
                       color=(color, 0.5), linewidth=3)
            end
        end

        # Plot wing points
        scatter!(ax, x_vals, y_vals;
                 markersize=point_size, color=color, marker=:circle)

        # Add point labels only for first struct
        if show_point_idxs && s_idx == 1
            for (i, p) in enumerate(plot_points)
                px, py = coords[i]
                away_x, away_y = 0.0, 0.0
                for (j, _) in enumerate(plot_points)
                    if i != j
                        ox, oy = coords[j]
                        dx, dy = px - ox, py - oy
                        dist = sqrt(dx^2 + dy^2) + 0.01
                        away_x += dx / dist^2
                        away_y += dy / dist^2
                    end
                end
                away_len = sqrt(away_x^2 + away_y^2)
                if away_len > 0
                    away_x /= away_len
                    away_y /= away_len
                else
                    away_x, away_y = 1.0, 1.0
                end
                offset = (12 * sign(away_x), 12 * sign(away_y))
                align_x = away_x >= 0 ? :left : :right
                align_y = away_y >= 0 ? :bottom : :top
                text!(ax, px, py; text=string(p.idx), fontsize=12,
                      align=(align_x, align_y), offset=offset)
            end
        end
    end

    # Plot extra points with connections (use first struct's wing for transform)
    if !isnothing(extra_points) && !isnothing(extra_groups)
        wing = structs[1].wings[1]
        R_w_b = wing.R_b_to_w'
        extra_body = [R_w_b * (collect(p) - wing.pos_w) for p in extra_points]

        if dir == :top
            extra_coords = [(p[1], p[2]) for p in extra_body]
        elseif dir == :side
            extra_coords = [(p[1], p[3]) for p in extra_body]
        else
            extra_coords = [(p[2], p[3]) for p in extra_body]
        end

        # Draw lines connecting points within each group
        for (gname, indices) in extra_groups
            for i in 1:(length(indices)-1)
                c1 = extra_coords[indices[i]]
                c2 = extra_coords[indices[i+1]]
                lines!(ax, [c1[1], c2[1]], [c1[2], c2[2]];
                       color=(:red, 0.6), linewidth=2)
            end
        end

        # Plot all extra points as circles
        ex_x = [c[1] for c in extra_coords]
        ex_y = [c[2] for c in extra_coords]
        scatter!(ax, ex_x, ex_y;
                 markersize=extra_point_size, color=:red, marker=:circle)
    end

    # Auto-zoom with margin
    if !isempty(all_x_vals)
        x_min, x_max = extrema(all_x_vals)
        y_min, y_max = extrema(all_y_vals)
        margin_x = 0.15 * (x_max - x_min) + 0.3
        margin_y = 0.15 * (y_max - y_min) + 0.3
        limits!(ax, x_min - margin_x, x_max + margin_x,
                    y_min - margin_y, y_max + margin_y)
    end

    # Legend
    if legend
        legend_elements = [
            MarkerElement(color=PLOT_COLORS[mod1(i, length(PLOT_COLORS))],
                          marker=:circle, markersize=10)
            for i in 1:n_structs
        ]
        legend_labels = copy(labels)
        if !isnothing(extra_points)
            push!(legend_elements,
                  MarkerElement(color=:red, marker=:circle, markersize=10))
            push!(legend_labels, "photogrammetry")
        end
        Legend(fig[1, 2], legend_elements, legend_labels)
    end

    return fig
end

"""
    wrap_to_pi(angle)

Wrap angle to [-π, π] range.
"""
function wrap_to_pi(angle)
    return mod(angle + π, 2π) - π
end

"""
    euler_to_quaternion(roll_deg, pitch_deg, yaw_deg)

Convert Euler angles (in degrees) to quaternion.
Converts from NED to ENU frame:
  X_ENU = Y_NED (East)
  Y_ENU = X_NED (North)
  Z_ENU = -Z_NED (Up = -Down)
"""
function euler_to_quaternion(roll_deg, pitch_deg, yaw_deg)
    roll_rad = deg2rad(roll_deg)
    pitch_rad = deg2rad(pitch_deg)
    yaw_rad = deg2rad(yaw_deg)
    rot_ned = RotZYX(yaw_rad, pitch_rad, roll_rad)
    R_ned_to_enu = [0.0 1.0 0.0;
                    1.0 0.0 0.0;
                    0.0 0.0 -1.0]
    rot_enu = R_ned_to_enu * Matrix(rot_ned)
    q = SymbolicAWEModels.rotation_matrix_to_quaternion(rot_enu)
    return q
end

"""
    calc_heading(sys_struct::SystemStructure, R_b_w)

Calculate heading angle from rotation matrix, wrapped to [-π, π].
"""
function calc_heading(sys_struct::SymbolicAWEModels.SystemStructure, R_b_w)
    e_x = R_b_w[:, 1]
    wind_norm = [1,0,0]
    minus_e_x = -e_x
    proj_on_wind = dot(minus_e_x, wind_norm) * wind_norm
    e_x_perp = minus_e_x - proj_on_wind
    wind_cross_z = [wind_norm[2], -wind_norm[1], 0]
    heading_x = dot(e_x_perp, wind_cross_z)
    heading_z = e_x_perp[3]
    heading = atan(heading_x, heading_z)
    return wrap_to_pi(heading)
end

"""
    calc_csv_heading(roll_deg, pitch_deg, yaw_deg, sys_struct)

Calculate heading from CSV Euler angles, wrapped to [-π, π].
"""
function calc_csv_heading(roll_deg, pitch_deg, yaw_deg, sys_struct)
    quat = euler_to_quaternion(roll_deg, pitch_deg, yaw_deg)
    R = SymbolicAWEModels.quaternion_to_rotation_matrix(quat)
    heading = calc_heading(sys_struct, R)
    return wrap_to_pi(heading + π)
end
