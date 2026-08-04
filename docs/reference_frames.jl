# SPDX-FileCopyrightText: 2025 Bart van de Lint
# SPDX-License-Identifier: MIT
#
# Reference frame visualization for kite heading calculation.
# Shows Up frame, Body frame, wind perpendicular plane, and heading projections
# at a single time point.

using GLMakie
using MakieControlPlots
using CairoMakie
using Makie: NoShading
using LinearAlgebra
using LaTeXStrings
using GeometryBasics: Point2f, Point3f, Vec3f, Cylinder

# ============================================================================
# Configuration variables (modify these as needed)
# ============================================================================
tether_length = 50.0
cone_angle_deg = 45.0
t_fraction = 0.8/8  # Fraction of period (1/8 = π/4)

# Output options
save_pdf = true
show_live = true

# ============================================================================
# Frame computation
# ============================================================================

"""Compute all reference frame vectors at a single time point."""
function compute_frame_at_t(t, cone_angle_deg; tether_len=tether_length)
    cone_angle = deg2rad(cone_angle_deg)
    distance = tether_len * cos(cone_angle)
    radius = tether_len * sin(cone_angle)

    # Kite position
    kite_pos = [distance, -radius * sin(t), radius * cos(t)]

    # Body frame: e_z radial, e_y velocity, e_x = e_y × e_z
    e_z = kite_pos / norm(kite_pos)
    velocity = [0.0, cos(t), sin(t)]
    e_y = velocity / norm(velocity)
    e_x = cross(e_y, e_z)

    # Tangent frame: u_z radial, u_y horizontal tangent, u_x = u_y × u_z (-u_x points "up")
    u_z = kite_pos / norm(kite_pos)
    y_vec = [-kite_pos[2], kite_pos[1], 0.0]
    u_y = y_vec / norm(y_vec)
    u_x = cross(u_y, u_z)

    # Wind frame: project -e_x onto wind-perpendicular plane (y-z)
    wind_norm = [1.0, 0.0, 0.0]
    minus_e_x = -e_x
    proj_on_wind = dot(minus_e_x, wind_norm) * wind_norm
    e_x_perp = minus_e_x - proj_on_wind  # -e_x projected onto y-z plane
    e_x_perp_norm = e_x_perp / norm(e_x_perp)

    # Also project u_x onto wind perp plane for comparison
    proj_up_on_wind = dot(u_x, wind_norm) * wind_norm
    up_x_perp = u_x - proj_up_on_wind
    up_x_perp_norm = up_x_perp / norm(up_x_perp)

    return (
        kite_pos = kite_pos,
        e_x = e_x,
        e_y = e_y,
        e_z = e_z,
        u_x = u_x,
        u_y = u_y,
        u_z = u_z,
        minus_e_x = minus_e_x,
        e_x_perp = e_x_perp,
        e_x_perp_norm = e_x_perp_norm,
        up_x_perp = up_x_perp,
        up_x_perp_norm = up_x_perp_norm,
        wind_norm = wind_norm,
    )
end

"""Create arc points between two vectors."""
function arc_between_vectors(origin, v1, v2, radius, n_points=15)
    v1_norm = v1 / norm(v1)
    v2_norm = v2 / norm(v2)

    # Angle between vectors
    angle = acos(clamp(dot(v1_norm, v2_norm), -1, 1))

    # Create arc using spherical interpolation
    points = Vector{Point3f}()
    for i in 0:n_points
        t = i / n_points
        # Slerp between v1 and v2
        theta = angle * t
        if angle > 1e-6
            v = sin((1-t)*angle)/sin(angle) * v1_norm + sin(t*angle)/sin(angle) * v2_norm
        else
            v = v1_norm
        end
        push!(points, Point3f(origin .+ radius * v))
    end
    return points
end

# ============================================================================
# Plotting
# ============================================================================

"""Create sphere guide lines: arcs on XY, XZ, YZ planes, and through kite position."""
function create_sphere_arcs(radius, kite_pos, n=20)
    arcs = Vector{Vector{Point3f}}()

    # Arc on XY plane (z=0): quarter circle from +x to -y
    xy_arc = [Point3f(radius * cos(φ), radius * sin(φ), 0) for φ in range(0, -π/2, length=n)]
    push!(arcs, xy_arc)

    # Arc on XZ plane (y=0): quarter circle from +z to +x
    xz_arc = [Point3f(radius * sin(θ), 0, radius * cos(θ)) for θ in range(0, π/2, length=n)]
    push!(arcs, xz_arc)

    # Arc on YZ plane (x=0): quarter circle from +z to -y
    yz_arc = [Point3f(0, -radius * sin(θ), radius * cos(θ)) for θ in range(0, π/2, length=n)]
    push!(arcs, yz_arc)

    # Meridian arc through kite (vertical great circle in plane containing z-axis and kite)
    φ_kite = atan(kite_pos[2], kite_pos[1])
    meridian_arc = [Point3f(radius * sin(θ) * cos(φ_kite), radius * sin(θ) * sin(φ_kite), radius * cos(θ))
                    for θ in range(0, π/2, length=n)]
    push!(arcs, meridian_arc)

    # Latitude arc through kite (horizontal circle at kite's z-height)
    # Only the portion in visible quadrant (x≥0, y≤0)
    θ_kite = acos(kite_pos[3] / radius)  # polar angle of kite
    lat_radius = radius * sin(θ_kite)    # radius of latitude circle
    z_kite = kite_pos[3]
    lat_arc = [Point3f(lat_radius * cos(φ), lat_radius * sin(φ), z_kite)
               for φ in range(-π/2, 0, length=n)]
    push!(arcs, lat_arc)

    return arcs
end

"""
Plot a flat 2D arrow in 3D space: line shaft + flat triangle head lying on a plane.
The arrow head is perpendicular to both the arrow direction and the plane normal.
"""
function arrow_flat!(ax, origin, direction, tip_width, tip_length, plane_normal;
                     color=:black, linewidth=2)
    origin = collect(origin)
    direction = collect(direction)
    plane_normal = collect(plane_normal)
    dir = normalize(direction)
    arrow_length = norm(direction)
    shaft_end = origin + (arrow_length - tip_length) * dir
    tip_apex = origin + arrow_length * dir

    # Perpendicular direction lying in the plane (perpendicular to both arrow and plane normal)
    perp = cross(dir, plane_normal)
    if norm(perp) < 1e-6
        perp = abs(dir[3]) < 0.9 ? cross(dir, [0.0, 0.0, 1.0]) : cross(dir, [1.0, 0.0, 0.0])
    end
    perp = normalize(perp) * tip_width

    # Draw shaft as line (extend slightly into arrowhead for clean connection)
    shaft_end_visual = origin + (arrow_length - tip_length * 0.7) * dir
    lines!(ax, [Point3f(origin...), Point3f(shaft_end_visual...)], color=color, linewidth=linewidth)

    # Draw flat triangle head
    tri_verts = [
        Point3f((shaft_end - perp)...),
        Point3f((shaft_end + perp)...),
        Point3f(tip_apex...)
    ]
    mesh!(ax, tri_verts, [1 2 3], color=color, shading=NoShading)
end

"""Create combined figure with 3D view and front view (YZ)."""
function plot_reference_frames(; use_glmakie=true)
    if use_glmakie
        GLMakie.activate!()
    else
        CairoMakie.activate!()
    end

    fig = Figure(size=(595*0.9, 420*0.9))
    rowsize!(fig.layout, 1, Auto(0.15))  # Legend row smaller
    fontsize = 15

    # Compute frame data once
    t = t_fraction * 2π
    frame = compute_frame_at_t(t, cone_angle_deg)
    kite_pos = frame.kite_pos
    x, y, z = kite_pos

    # Arrow scaling
    s = tether_length / 10
    scale = 3.0 * s
    arrowscale = 1.5s
    tip_w = 0.2 * arrowscale  # tip width for flat arrows
    tip_l = 0.4 * arrowscale   # tip length for flat arrows
    arrow_lw = 2  # line width for arrow shafts
    # 2D arrow parameters (need larger values for visibility in 2D coordinates)
    s_2d = arrowscale * 3
    shaft_w_2d = 0.12 * s_2d
    tip_w_2d = 0.36 * s_2d
    tip_l_2d = 0.6 * s_2d

    # Projections of -e_x and -u_x onto y-z plane
    minus_e_x = -frame.e_x
    minus_u_x = -frame.u_x
    e_x_yz = [minus_e_x[2], minus_e_x[3]]
    e_x_yz_norm = e_x_yz / norm(e_x_yz)
    u_x_yz = [minus_u_x[2], minus_u_x[3]]
    u_x_yz_norm = u_x_yz / norm(u_x_yz)

    # ===== HORIZONTAL LEGEND (top, spanning both columns) =====
    legend_entries = []
    legend_labels = []
    push!(legend_entries, LineElement(color=:black, linewidth=2))
    push!(legend_labels, L"\textrm{Tether}")
    push!(legend_entries, PolyElement(color=(:lightblue, 0.4)))
    push!(legend_labels, L"\textrm{W plane}")
    push!(legend_entries, PolyElement(color=(:lightgreen, 0.4)))
    push!(legend_labels, L"\textrm{T plane}")
    push!(legend_entries, MarkerElement(marker=:rtriangle, color=:red, markersize=12))
    push!(legend_labels, L"-\hat{e}_x")
    push!(legend_entries, MarkerElement(marker=:rtriangle, color=:limegreen, markersize=12))
    push!(legend_labels, L"-\hat{u}_x")
    push!(legend_entries, MarkerElement(marker=:rtriangle, color=:darkgreen, markersize=12))
    push!(legend_labels, L"\hat{u}_y")
    push!(legend_entries, MarkerElement(marker=:rtriangle, color=:orange, markersize=12))
    push!(legend_labels, L"-\hat{e}_x|_{yz}")
    push!(legend_entries, MarkerElement(marker=:rtriangle, color=:gray, markersize=12))
    push!(legend_labels, L"\hat{y}, \hat{z}")
    push!(legend_entries, MarkerElement(marker=:circle, color=:blue, markersize=10))
    push!(legend_labels, L"O_\textrm{body}")
    push!(legend_entries, LineElement(color=:purple, linewidth=2))
    push!(legend_labels, L"\textrm{T}")
    push!(legend_entries, LineElement(color=:darkorange, linewidth=2))
    push!(legend_labels, L"\textrm{W}")

    Legend(fig[1, 1:2], legend_entries, legend_labels, labelsize=fontsize, orientation=:horizontal, nbanks=2)

    # ===== 3D VIEW (left) =====
    ax3d = Axis3(fig[2, 1],
                 xlabel=L"X \textrm{ (m)}", ylabel=L"Y \textrm{ (m)}", zlabel=L"Z \textrm{ (m)}",
                 aspect=:data, azimuth=-0.3π, elevation=π/6,
                 xlabelsize=fontsize, ylabelsize=fontsize, zlabelsize=fontsize,
                 xticklabelsize=fontsize*0.6, yticklabelsize=fontsize*0.6, zticklabelsize=fontsize*0.6,
                 xlabeloffset=30, ylabeloffset=30, zlabeloffset=30)

    # Tether
    lines!(ax3d, [0, x], [0, y], [0, z], color=:black, linewidth=2)

    # Sphere guide arcs
    sphere_arcs = create_sphere_arcs(tether_length, kite_pos)
    for arc in sphere_arcs
        lines!(ax3d, arc, color=:black, linewidth=1)
    end

    # Wind perp plane
    plane_size = 1.0 * scale
    plane_vertices = [
        Point3f(0, y - plane_size, z - plane_size),
        Point3f(0, y + plane_size, z - plane_size),
        Point3f(0, y + plane_size, z + plane_size),
        Point3f(0, y - plane_size, z + plane_size),
    ]
    plane_faces = [1 2 3; 1 3 4]
    mesh!(ax3d, plane_vertices, plane_faces, color=(:lightblue, 0.4), shading=NoShading)

    # Tangential plane at kite (perpendicular to tether/radial direction)
    # Offset slightly towards origin so it doesn't block arrow views
    tan_plane_center = kite_pos - 0.15 * scale * frame.u_z
    tan_plane_vertices = [
        Point3f((tan_plane_center - plane_size * frame.u_x - plane_size * frame.u_y)...),
        Point3f((tan_plane_center + plane_size * frame.u_x - plane_size * frame.u_y)...),
        Point3f((tan_plane_center + plane_size * frame.u_x + plane_size * frame.u_y)...),
        Point3f((tan_plane_center - plane_size * frame.u_x + plane_size * frame.u_y)...),
    ]
    mesh!(ax3d, tan_plane_vertices, plane_faces, color=(:lightgreen, 0.4), shading=NoShading)

    proj_origin = [scale * 0.1, y, z]

    # Arrows on T plane (tangent plane at kite, normal = u_z)
    t_plane_normal = frame.u_z
    arrow_flat!(ax3d, kite_pos, -frame.e_x * scale, tip_w, tip_l, t_plane_normal;
                color=:red, linewidth=arrow_lw)
    arrow_flat!(ax3d, kite_pos, -frame.u_x * scale, tip_w, tip_l, t_plane_normal;
                color=:limegreen, linewidth=arrow_lw)
    arrow_flat!(ax3d, kite_pos, frame.u_y * scale, tip_w, tip_l, t_plane_normal;
                color=:darkgreen, linewidth=arrow_lw)

    # Arrows on W plane (wind perp plane YZ, normal = x_hat)
    w_plane_normal = [1.0, 0.0, 0.0]
    minus_e_x = -frame.e_x
    e_x_proj = [0.0, minus_e_x[2], minus_e_x[3]]
    e_x_proj_norm = e_x_proj / norm(e_x_proj)
    arrow_flat!(ax3d, proj_origin, e_x_proj_norm * scale, tip_w, tip_l, w_plane_normal;
                color=:orange, linewidth=arrow_lw)
    # World frame unit vectors (y and z) at projection origin
    w_y = [0.0, 1.0, 0.0]
    w_z = [0.0, 0.0, 1.0]
    arrow_flat!(ax3d, proj_origin, w_y * scale, tip_w, tip_l, w_plane_normal;
                color=:gray, linewidth=arrow_lw)
    arrow_flat!(ax3d, proj_origin, w_z * scale, tip_w, tip_l, w_plane_normal;
                color=:gray, linewidth=arrow_lw)

    # Origin markers (simple scatter instead of mesh spheres)
    scatter!(ax3d, [Point3f(0, 0, 0)], markersize=8, color=:black)       # world origin
    scatter!(ax3d, [Point3f(x, y, z)], markersize=8, color=:blue)        # kite/body origin
    scatter!(ax3d, [Point3f(proj_origin...)], markersize=8, color=:gray) # projected origin

    # Arcs
    arc_radius = scale * 0.4
    arc_points = arc_between_vectors(kite_pos, -frame.e_x, -frame.u_x, arc_radius)
    lines!(ax3d, arc_points, color=:purple, linewidth=2)
    arc_points_proj = arc_between_vectors(proj_origin, e_x_proj_norm, w_z, arc_radius)
    lines!(ax3d, arc_points_proj, color=:darkorange, linewidth=2)

    margin = scale * 1.2
    xlims!(ax3d, 0, tether_length)
    ylims!(ax3d, -tether_length, 0)
    zlims!(ax3d, 0, tether_length)

    # ===== FRONT VIEW (YZ, right) =====
    ax_front = Axis(fig[2, 2],
                    xlabel=L"Y \textrm{ (m)}",
                    aspect=DataAspect(),
                    xlabelsize=fontsize,
                    xticklabelsize=fontsize*0.6, yticklabelsize=fontsize*0.6)

    lines!(ax_front, [0, y], [0, z], color=:black, linewidth=2)
    p = poly!(ax_front, Point2f[(y - plane_size, z - plane_size), (y + plane_size, z - plane_size),
                            (y + plane_size, z + plane_size), (y - plane_size, z + plane_size)],
          color=(:lightblue, 0.3))
    translate!(p, 0, 0, -1)  # Push poly behind labels

    arrows2d!(ax_front, [y], [z], [e_x_yz[1] * scale], [e_x_yz[2] * scale],
              color=:red, shaftwidth=shaft_w_2d, tipwidth=tip_w_2d, tiplength=tip_l_2d)
    arrows2d!(ax_front, [y], [z], [u_x_yz[1] * scale], [u_x_yz[2] * scale],
              color=:limegreen, shaftwidth=shaft_w_2d, tipwidth=tip_w_2d, tiplength=tip_l_2d)
    arrows2d!(ax_front, [y], [z], [e_x_yz[1] * scale], [e_x_yz[2] * scale],
              color=:orange, shaftwidth=shaft_w_2d, tipwidth=tip_w_2d, tiplength=tip_l_2d)
    arrows2d!(ax_front, [y], [z], [0.0], [scale],
              color=:gray, shaftwidth=shaft_w_2d, tipwidth=tip_w_2d, tiplength=tip_l_2d)

    scatter!(ax_front, [0], [0], markersize=scale*0.9, color=:black)  # world origin
    scatter!(ax_front, [y], [z], markersize=scale*0.9, color=:blue)   # kite body origin

    angle_body_yz = atan(e_x_yz_norm[1], e_x_yz_norm[2])
    angle_up_yz = atan(u_x_yz_norm[1], u_x_yz_norm[2])
    arc_angles = range(angle_body_yz, angle_up_yz, length=30)
    lines!(ax_front, y .+ arc_radius .* sin.(arc_angles), z .+ arc_radius .* cos.(arc_angles), color=:purple, linewidth=2)
    arc_angles_proj = range(0.0, angle_body_yz, length=30)
    lines!(ax_front, y .+ (arc_radius*0.8) .* sin.(arc_angles_proj), z .+ (arc_radius*0.8) .* cos.(arc_angles_proj), color=:darkorange, linewidth=2)

    xlims!(ax_front, min(y, 0) - margin, 0)
    ylims!(ax_front, 0, z + margin)

    rowgap!(fig.layout, 1, -0)  # Reduce gap between legend and figures
    colgap!(fig.layout, 1, 30)  # Add padding between left and right figures
    return fig
end

# ============================================================================
# Main execution
# ============================================================================

# Combined figure with all views
fig_combined = plot_reference_frames(use_glmakie=false)

if show_live
    GLMakie.activate!()
    fig_live = plot_reference_frames(use_glmakie=true)
    scrn = display(fig_live)
    wait(scrn)
end

if save_pdf
    CairoMakie.activate!()
    save("reference_frames.pdf", fig_combined)
    println("Saved: reference_frames.pdf")
end

println("\nDone!")
GLMakie.activate!()
