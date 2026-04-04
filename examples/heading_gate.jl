# SPDX-FileCopyrightText: 2025 Bart van de Lint
# SPDX-License-Identifier: MIT
#
# Kite heading calculation comparing Tangential Sphere vs Wind Frame methods.
# Assumes: downwind flight, constant tether, zero gravity, straight tether, no sideslip.
# Expects constant turn rate (uniform circular motion).

using Pkg
Pkg.activate(@__DIR__)

using GLMakie
using CairoMakie
using Makie: NoShading
using LinearAlgebra
using LaTeXStrings
using Statistics
using Colors: HSV

# ============================================================================
# Configuration variables (modify these as needed)
# ============================================================================
tether_length = 50.0  # Length of tether
cone_angles = [30.0, 45.0, 60.0, 75.0]  # Cone angles in degrees (comma-separated in original)

# Visualization options
show_body_x = true
show_body_y = false
show_body_z = false
show_up_x = true
show_up_y = false
show_up_z = false

# Output options
save_pdf = true
show_live = false

"""Tangential sphere heading from body velocity direction and up frame basis vectors."""
function tangential_sphere_heading(R_body_y, R_up_x, R_up_y, R_up_z)
    R_up = hcat(R_up_x, R_up_y, R_up_z)
    heading_vec = R_up' * R_body_y
    heading = atan(heading_vec[1], heading_vec[2])
    return heading
end

"""Wind frame heading by projecting body x-axis onto wind-perpendicular plane."""
function wind_frame_heading(R_body_x, wind_vel)
    wind_norm = wind_vel / norm(wind_vel)

    minus_e_x = -R_body_x
    proj_on_wind = dot(minus_e_x, wind_norm) * wind_norm
    e_x_perp = minus_e_x - proj_on_wind

    wind_cross_z = [wind_norm[2], -wind_norm[1], 0.0]
    heading_x = dot(e_x_perp, wind_cross_z)

    heading_z = e_x_perp[3]

    heading = atan(heading_x, heading_z)
    return heading
end

"""Run simulation for given cone angle, returns NamedTuple with time series and error metrics."""
function run_simulation(cone_angle_deg; tether_len=tether_length)
    cone_angle = deg2rad(cone_angle_deg)

    t = range(1e-8, 2π, length=500)  # one period, avoid t=0

    distance = tether_len * cos(cone_angle)
    radius = tether_len * sin(cone_angle)

    x = fill(distance, length(t))
    y = -radius .* sin.(t)
    z = radius .* cos.(t)

    kite_pos = hcat(x, y, z)

    n = length(t)
    R_body_x = zeros(n, 3)
    R_body_y = zeros(n, 3)
    R_body_z = zeros(n, 3)

    R_up_x = zeros(n, 3)
    R_up_y = zeros(n, 3)
    R_up_z = zeros(n, 3)

    heading = zeros(n)

    for i in eachindex(t)
        R_body_z[i, :] = kite_pos[i, :] / norm(kite_pos[i, :])
        velocity = [0.0, -cos(t[i]), -sin(t[i])]
        R_body_y[i, :] = velocity / norm(velocity)
        R_body_x[i, :] = cross(R_body_y[i, :], R_body_z[i, :])

        R_up_z[i, :] = kite_pos[i, :] / norm(kite_pos[i, :])
        y_vec = [-kite_pos[i, 2], kite_pos[i, 1], 0.0]
        R_up_y[i, :] = y_vec / norm(y_vec)
        R_up_x[i, :] = cross(R_up_z[i, :], R_up_y[i, :])

        heading[i] = tangential_sphere_heading(
            R_body_y[i, :], R_up_x[i, :], R_up_y[i, :], R_up_z[i, :])
    end

    wind_vel = [1.0, 0.0, 0.0]
    heading_wind = zeros(n)
    for i in eachindex(t)
        heading_wind[i] = wind_frame_heading(R_body_x[i, :], wind_vel)
    end

    heading_unwrapped = unwrap(heading)
    dt = t[2] - t[1]
    heading_rate = gradient(heading_unwrapped, dt)

    heading_wind_unwrapped = unwrap(heading_wind)
    heading_wind_rate = gradient(heading_wind_unwrapped, dt)

    turn_rate_err_sphere = abs(maximum(heading_rate) - minimum(heading_rate))
    turn_rate_err_wind = abs(maximum(heading_wind_rate) - minimum(heading_wind_rate))

    mean_turn_rate_sphere = mean(abs.(heading_rate))
    mean_turn_rate_wind = mean(abs.(heading_wind_rate))

    turn_rate_rel_err_sphere = turn_rate_err_sphere / mean_turn_rate_sphere
    turn_rate_rel_err_wind = turn_rate_err_wind / mean_turn_rate_wind

    return (
        t = collect(t),
        y = y,
        z = z,
        x = x,
        heading = heading,
        heading_wind = heading_wind,
        heading_rate = heading_rate,
        heading_wind_rate = heading_wind_rate,
        turn_rate_err_sphere = turn_rate_err_sphere,
        turn_rate_err_wind = turn_rate_err_wind,
        turn_rate_rel_err_sphere = turn_rate_rel_err_sphere,
        turn_rate_rel_err_wind = turn_rate_rel_err_wind,
        R_body_x = R_body_x,
        R_body_y = R_body_y,
        R_body_z = R_body_z,
        R_up_x = R_up_x,
        R_up_y = R_up_y,
        R_up_z = R_up_z,
    )
end

"""Unwrap phase angles to avoid discontinuities."""
function unwrap(phase)
    result = copy(phase)
    for i in 2:lastindex(result)
        diff = result[i] - result[i-1]
        if diff > π
            result[i:end] .-= 2π
        elseif diff < -π
            result[i:end] .+= 2π
        end
    end
    return result
end

"""Compute numerical gradient."""
function gradient(y, dt)
    n = length(y)
    grad = zeros(n)
    grad[1] = (y[2] - y[1]) / dt
    grad[end] = (y[end] - y[end-1]) / dt
    for i in 2:lastindex(y)-1
        grad[i] = (y[i+1] - y[i-1]) / (2dt)
    end
    return grad
end

"""Get color for cone angle: green (min) to red (max) using HSV for vibrant gradient."""
cone_color(angle, min_angle, max_angle) = HSV(120 * (1 - (angle - min_angle) / (max_angle - min_angle)), 1.0, 0.9)

"""Plot combined trajectory, heading, and turn rate for multiple cone angles."""
function plot_combined(cone_angles_list, results; use_glmakie=true)
    if use_glmakie
        GLMakie.activate!()
    else
        CairoMakie.activate!()
    end
    fig = Figure(size=(595*0.9, 420*0.9))

    fontsize = 12
    ax1 = Axis(fig[1, 1], ylabel=L"\textrm{pos (m)}", xticklabelsvisible=false,
               ylabelsize=fontsize, xlabelsize=fontsize, xticklabelsize=fontsize, yticklabelsize=fontsize)
    ax2 = Axis(fig[2, 1], ylabel=L"\psi \textrm{ (deg)}", xticklabelsvisible=false,
               ylabelsize=fontsize, xlabelsize=fontsize, xticklabelsize=fontsize, yticklabelsize=fontsize)
    ax3 = Axis(fig[3, 1], xlabel=L"t \textrm{ (s)}", ylabel=L"\dot{\psi} \textrm{ (deg/s)}",
               ylabelsize=fontsize, xlabelsize=fontsize, xticklabelsize=fontsize, yticklabelsize=fontsize)

    rowgap!(fig.layout, 1, 5)  # gap between row 1 and row 2
    rowgap!(fig.layout, 2, 5)  # gap between row 2 and row 3

    linkxaxes!(ax1, ax2, ax3)

    min_a, max_a = extrema(cone_angles_list)

    # Shared angle legend at top with error
    angle_elems = [LineElement(color=cone_color(a, min_a, max_a)) for a in cone_angles_list]
    angle_labels = ["$(Int(a))° ($(round(Int, res.turn_rate_rel_err_sphere * 100))%)"
                    for (a, res) in zip(cone_angles_list, results)]
    Legend(fig[0, 1], angle_elems, angle_labels, L"\textrm{Cone angle (error)}",
           orientation=:horizontal, tellwidth=false, titleposition=:left, halign=1.0,
           labelsize=fontsize, titlesize=fontsize)

    for (angle, res) in zip(cone_angles_list, results)
        c = cone_color(angle, min_a, max_a)
        lines!(ax1, res.t, res.y, color=c)
        lines!(ax1, res.t, res.z, color=c, linestyle=:dash)
    end
    Legend(fig[1, 1],
           [LineElement(color=:black), LineElement(color=:black, linestyle=:dash)],
           [L"y", L"z"],
           tellwidth=false, tellheight=false, halign=:right, valign=:bottom,
           labelsize=fontsize, rowgap=0, padding=(6,6,2,2))

    ref_res = run_simulation(80)

    for (angle, res) in zip(cone_angles_list, results)
        c = cone_color(angle, min_a, max_a)
        lines!(ax2, res.t, rad2deg.(res.heading), color=c)
    end
    lines!(ax2, ref_res.t, rad2deg.(ref_res.heading_wind), linestyle=:dot, color=:gray)
    Legend(fig[2, 1],
           [LineElement(color=:black), LineElement(color=:gray, linestyle=:dot)],
           ["T", "W"],
           tellwidth=false, tellheight=false, halign=:right, valign=:bottom,
           labelsize=fontsize, rowgap=0, padding=(6,6,2,2))

    for (angle, res) in zip(cone_angles_list, results)
        c = cone_color(angle, min_a, max_a)
        lines!(ax3, res.t, rad2deg.(res.heading_rate), color=c)
    end
    lines!(ax3, ref_res.t, rad2deg.(ref_res.heading_wind_rate), linestyle=:dot, color=:gray)
    Legend(fig[3, 1],
           [LineElement(color=:black), LineElement(color=:gray, linestyle=:dot)],
           ["T", "W"],
           tellwidth=false, tellheight=false, halign=:right, valign=:bottom,
           labelsize=fontsize, rowgap=0, padding=(6,6,2,2))

    return fig
end

"""Plot relative turn rate error vs cone angle."""
function plot_error_vs_cone_angle(; use_glmakie=true)
    angles = 0:80
    rel_err_sphere = Float64[]
    rel_err_wind = Float64[]

    for angle in angles
        res = run_simulation(angle)
        push!(rel_err_sphere, res.turn_rate_rel_err_sphere)
        push!(rel_err_wind, res.turn_rate_rel_err_wind)
    end

    if use_glmakie
        GLMakie.activate!()
    else
        CairoMakie.activate!()
    end
    fig = Figure(size=(595*0.9, 420*0.9))
    fontsize = 12

    ax = Axis(fig[1, 1],
              xlabel=L"\textrm{Cone Angle (degrees)}",
              ylabel=L"\textrm{Relative Error}",
              xlabelsize=fontsize, ylabelsize=fontsize,
              xticklabelsize=fontsize, yticklabelsize=fontsize)

    lines!(ax, collect(angles), rel_err_sphere, label=L"\textrm{Tangential Sphere}")
    lines!(ax, collect(angles), rel_err_wind, label=L"\textrm{Wind perp.}")
    axislegend(ax, position=:lt, labelsize=fontsize)

    return fig
end

"""Create 3D trajectory plot for a single cone angle."""
function plot_3d_trajectory(_cone_angle_deg, res; use_glmakie=true)
    if use_glmakie
        GLMakie.activate!()
    else
        CairoMakie.activate!()
    end
    fig = Figure(size=(595*0.9, 420*0.9))
    fontsize = 12

    ax = Axis3(fig[1, 1],
               xlabel=L"X \textrm{ (m)}", ylabel=L"Y \textrm{ (m)}", zlabel=L"Z \textrm{ (m)}",
               aspect=:data, azimuth=-0.3π,
               xlabelsize=fontsize, ylabelsize=fontsize, zlabelsize=fontsize,
               xticklabelsize=fontsize, yticklabelsize=fontsize, zticklabelsize=fontsize)

    x, y, z = res.x, res.y, res.z
    t = res.t

    lines!(ax, x, y, z, color=:dodgerblue, linewidth=2, label=L"\textrm{Trajectory}")
    lines!(ax, [0, x[1]], [0, y[1]], [0, z[1]], color=:black, linewidth=2, label=L"\textrm{Tether}")

    step = 25
    # Arrow dimensions - all scale with tether_length (reference: tether_length=10)
    s = tether_length / 10
    scale = 2.0 * s
    shaft_r = 0.006 * s
    tip_r = 0.018 * s
    tip_l = 0.03 * s

    # Distinct colors for body (red/orange/blue) and up (greens) frames
    color_body_x = :red
    color_body_y = :orange
    color_body_z = :blue
    color_up_x = :limegreen
    color_up_y = :green
    color_up_z = :darkgreen

    first_body_x = true
    first_body_y = true
    first_body_z = true
    first_up_x = true
    first_up_y = true
    first_up_z = true

    quality = 8
    for i in firstindex(t):step:lastindex(t)
        if show_body_x
            label = first_body_x ? L"R_{\textrm{body},x}" : nothing
            first_body_x = false
            arrows3d!(ax, [Point3f(x[i], y[i], z[i])],
                    [Vec3f(res.R_body_x[i, :]...) * scale],
                    color=color_body_x, shaftradius=shaft_r, tipradius=tip_r, tiplength=tip_l, quality=quality, shading=NoShading, label=label)
        end

        if show_body_y
            label = first_body_y ? L"R_{\textrm{body},y}" : nothing
            first_body_y = false
            arrows3d!(ax, [Point3f(x[i], y[i], z[i])],
                    [Vec3f(res.R_body_y[i, :]...) * scale],
                    color=color_body_y, shaftradius=shaft_r, tipradius=tip_r, tiplength=tip_l, quality=quality, shading=NoShading, label=label)
        end

        if show_body_z
            label = first_body_z ? L"R_{\textrm{body},z}" : nothing
            first_body_z = false
            arrows3d!(ax, [Point3f(x[i], y[i], z[i])],
                    [Vec3f(res.R_body_z[i, :]...) * scale],
                    color=color_body_z, shaftradius=shaft_r, tipradius=tip_r, tiplength=tip_l, quality=quality, shading=NoShading, label=label)
        end

        if show_up_x
            label = first_up_x ? L"R_{\textrm{up},x}" : nothing
            first_up_x = false
            arrows3d!(ax, [Point3f(x[i], y[i], z[i])],
                    [Vec3f(res.R_up_x[i, :]...) * scale],
                    color=color_up_x, shaftradius=shaft_r, tipradius=tip_r, tiplength=tip_l, quality=quality, shading=NoShading, label=label)
        end

        if show_up_y
            label = first_up_y ? L"R_{\textrm{up},y}" : nothing
            first_up_y = false
            arrows3d!(ax, [Point3f(x[i], y[i], z[i])],
                    [Vec3f(res.R_up_y[i, :]...) * scale],
                    color=color_up_y, shaftradius=shaft_r, tipradius=tip_r, tiplength=tip_l, quality=4, shading=NoShading, label=label)
        end

        if show_up_z
            label = first_up_z ? L"R_{\textrm{up},z}" : nothing
            first_up_z = false
            arrows3d!(ax, [Point3f(x[i], y[i], z[i])],
                    [Vec3f(res.R_up_z[i, :]...) * scale],
                    color=color_up_z, shaftradius=shaft_r, tipradius=tip_r, tiplength=tip_l, quality=4, shading=NoShading, label=label)
        end
    end

    # Expand limits to include arrows (add scale as margin)
    max_yz = max(maximum(abs.(y)), maximum(abs.(z))) + scale
    xlims!(ax, -scale, x[1] * 1.5 + scale)
    ylims!(ax, -max_yz, max_yz)
    zlims!(ax, -max_yz, max_yz)

    # Create custom legend with arrow-like markers
    legend_entries = []
    legend_labels = []

    # Trajectory and tether
    push!(legend_entries, LineElement(color=:blue, linewidth=2))
    push!(legend_labels, L"\textrm{Trajectory}")
    push!(legend_entries, LineElement(color=:black, linewidth=2))
    push!(legend_labels, L"\textrm{Tether}")

    # Body frame arrows
    if show_body_x
        push!(legend_entries, MarkerElement(marker=:rtriangle, color=color_body_x, markersize=12))
        push!(legend_labels, L"R_{\textrm{body},x}")
    end
    if show_body_y
        push!(legend_entries, MarkerElement(marker=:rtriangle, color=color_body_y, markersize=12))
        push!(legend_labels, L"R_{\textrm{body},y}")
    end
    if show_body_z
        push!(legend_entries, MarkerElement(marker=:rtriangle, color=color_body_z, markersize=12))
        push!(legend_labels, L"R_{\textrm{body},z}")
    end

    # Up frame arrows
    if show_up_x
        push!(legend_entries, MarkerElement(marker=:rtriangle, color=color_up_x, markersize=12))
        push!(legend_labels, L"R_{\textrm{up},x}")
    end
    if show_up_y
        push!(legend_entries, MarkerElement(marker=:rtriangle, color=color_up_y, markersize=12))
        push!(legend_labels, L"R_{\textrm{up},y}")
    end
    if show_up_z
        push!(legend_entries, MarkerElement(marker=:rtriangle, color=color_up_z, markersize=12))
        push!(legend_labels, L"R_{\textrm{up},z}")
    end

    Legend(fig[1, 2], legend_entries, legend_labels, labelsize=fontsize)

    return fig
end

"""Save figure to PDF using CairoMakie."""
function save_figure(basename, plot_func, args...; kwargs...)
    fig = plot_func(args...; use_glmakie=false, kwargs...)
    save("$basename.pdf", fig)
    println("Saved: $basename.pdf")
    return fig
end

# ============================================================================
# Main execution
# ============================================================================

results = []
figs_traj = Dict{Int,Figure}()
fig_combined = nothing
fig_error = nothing

for angle in cone_angles
    res = run_simulation(angle)
    push!(results, res)

    println("\n=== Cone angle: $angle degrees ===")
    println("Turn rate error (Tangential Sphere): $(rad2deg(res.turn_rate_err_sphere)) deg/s, " *
            "relative: $(res.turn_rate_rel_err_sphere)")
    println("Turn rate error (Wind Frame): $(rad2deg(res.turn_rate_err_wind)) deg/s, " *
            "relative: $(res.turn_rate_rel_err_wind)")

    figs_traj[Int(angle)] = plot_3d_trajectory(angle, res; use_glmakie=false)
    if show_live
        scrn = display(figs_traj[Int(angle)])
        wait(scrn)
    end

    if save_pdf
        # Save 3D trajectory as PDF using CairoMakie
        save("trajectory_$(Int(angle)).pdf", figs_traj[Int(angle)])
        println("Saved: trajectory_$(Int(angle)).pdf")
    end
end

open("errors.txt", "w") do f
    write(f, "Relative Turn Rate Errors\n")
    write(f, "=========================\n\n")
    for (angle, res) in zip(cone_angles, results)
        write(f, "Cone angle: $(angle)°\n")
        write(f, "  Tangential Sphere: $(res.turn_rate_rel_err_sphere)\n")
        write(f, "  Wind perp.: $(res.turn_rate_rel_err_wind)\n\n")
    end
end
println("\nSaved: errors.txt")

if length(cone_angles) > 0
    fig_combined = plot_combined(cone_angles, results; use_glmakie=true)
    if show_live
        scrn = display(fig_combined)
        wait(scrn)
    end

    if save_pdf
        save_figure("combined", plot_combined, cone_angles, results)
    end
end

println("\nGenerating error vs cone angle plot (0-80 degrees)...")
fig_error = plot_error_vs_cone_angle(; use_glmakie=true)
if show_live
    scrn = display(fig_error)
    wait(scrn)
end

if save_pdf
    save_figure("error_vs_cone_angle", plot_error_vs_cone_angle)
end

println("\nDone!")
GLMakie.activate!()

