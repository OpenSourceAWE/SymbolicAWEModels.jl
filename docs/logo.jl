# Copyright (c) 2025 Jelle Poland, Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

"""
Create a logo visualization of the V3 kite in Julia's official colors.

2D projection of the kite viewed from the front, split into 3 sections:
- Left: Julia Red (#CB3C33)
- Center: Julia Green (#389826)
- Right: Julia Purple (#9558B2)

Usage:
    julia --project=examples examples/logo.jl
"""

using CairoMakie
using FileIO
using GeometryBasics
using LinearAlgebra
using Statistics
using Colors

# Activate CairoMakie for proper display and PDF/SVG export
CairoMakie.activate!()

# Julia official colors
const JULIA_RED = colorant"#CB3C33"
const JULIA_GREEN = colorant"#389826"
const JULIA_PURPLE = colorant"#9558B2"

# ============= ADJUSTABLE PARAMETERS =============
# Text settings
TEXT_FONTSIZE = 70
# Modern fat round fonts: "Rubik Bold", "Lexend Bold", "Outfit Bold", "Fredoka Bold"
# Install: yay -S ttf-rubik ttf-lexend
TEXT_FONT = "Liberation sans bold"
TEXT_SPACING = 2.7          # Horizontal spacing between words
TEXT_Y_OFFSET = 1.0       # Vertical offset from bottom of kite (negative = below)

# Color split positions (0-1 range along wingspan)
COLOR_SPLIT_LEFT = 0.25     # Below this = red
COLOR_SPLIT_RIGHT = 0.75    # Above this = purple, between = green

# Figure size
FIG_WIDTH = 800
FIG_HEIGHT = 400
# =================================================

# Load the OBJ file
obj_path = joinpath(@__DIR__, "..", "data", "v3", "V3_25.obj")
mesh_3d = load(obj_path)

# Rotation matrices
function rot_x(angle)
    c, s = cos(angle), sin(angle)
    [1 0 0; 0 c -s; 0 s c]
end

function rot_z(angle)
    c, s = cos(angle), sin(angle)
    [c -s 0; s c 0; 0 0 1]
end

# Combined rotation: 90deg around x, then -90deg around z
R = rot_z(-π/2) * rot_x(π/2)

# Transform vertices to world frame
old_vertices = coordinates(mesh_3d)
transformed = [R * Vector(v) + [0, 0, 7.3] for v in old_vertices]

# Project to 2D (front view: use y, z as x, y)
vertices_2d = [Point2f(v[2], v[3]) for v in transformed]

# Find x-range for coloring (wingspan direction)
x_coords = [v[1] for v in vertices_2d]
x_min, x_max = extrema(x_coords)
x_range = x_max - x_min

# Color based on x position
function get_julia_color(x)
    normalized = (x - x_min) / x_range
    if normalized < COLOR_SPLIT_LEFT
        return JULIA_RED
    elseif normalized < COLOR_SPLIT_RIGHT
        return JULIA_GREEN
    else
        return JULIA_PURPLE
    end
end

# Assign color to each vertex based on its x position
vertex_colors = [get_julia_color(v[1]) for v in vertices_2d]

# Create 2D mesh
mesh_faces = faces(mesh_3d)
mesh_2d = GeometryBasics.Mesh(vertices_2d, mesh_faces)

# Create figure with transparent background
fig = Figure(size=(FIG_WIDTH, FIG_HEIGHT), backgroundcolor=:transparent)
ax = Axis(fig[1, 1];
    aspect=DataAspect(),
    backgroundcolor=:transparent)

# Hide all decorations for clean logo look
hidedecorations!(ax)
hidespines!(ax)

# Plot the 2D mesh with vertex colors
mesh!(ax, mesh_2d; color=vertex_colors)

# Add text "Open Source AWE" - each word in a Julia color
text_y = minimum(v[2] for v in vertices_2d) + TEXT_Y_OFFSET
text_x_center = (x_min + x_max) / 2

text!(ax, text_x_center - TEXT_SPACING, text_y; text="Open", color=JULIA_RED,
      fontsize=TEXT_FONTSIZE, font=TEXT_FONT, align=(:center, :top))
text!(ax, text_x_center, text_y; text="Source", color=JULIA_GREEN,
      fontsize=TEXT_FONTSIZE, font=TEXT_FONT, align=(:center, :top))
text!(ax, text_x_center + TEXT_SPACING, text_y; text="AWE", color=JULIA_PURPLE,
      fontsize=TEXT_FONTSIZE, font=TEXT_FONT, align=(:center, :top))

# Display figure
display(fig)

# Save as PNG (small file size) and SVG
save("kite_logo.png", fig; px_per_unit=3, backend=CairoMakie)
save("kite_logo.svg", fig; pt_per_unit=1, backend=CairoMakie)

# Create square version
fig_square = Figure(size=(FIG_WIDTH, FIG_WIDTH), backgroundcolor=:transparent)
ax_square = Axis(fig_square[1, 1];
    aspect=DataAspect(),
    backgroundcolor=:transparent)
hidedecorations!(ax_square)
hidespines!(ax_square)
mesh!(ax_square, mesh_2d; color=vertex_colors)
text!(ax_square, text_x_center - TEXT_SPACING, text_y; text="Open", color=JULIA_RED,
      fontsize=TEXT_FONTSIZE, font=TEXT_FONT, align=(:center, :top))
text!(ax_square, text_x_center, text_y; text="Source", color=JULIA_GREEN,
      fontsize=TEXT_FONTSIZE, font=TEXT_FONT, align=(:center, :top))
text!(ax_square, text_x_center + TEXT_SPACING, text_y; text="AWE", color=JULIA_PURPLE,
      fontsize=TEXT_FONTSIZE, font=TEXT_FONT, align=(:center, :top))
save("kite_logo_square.png", fig_square; px_per_unit=3, backend=CairoMakie)

@info "Logo saved to kite_logo.png, kite_logo.svg, and kite_logo_square.png"

nothing
