# Copyright (c) 2025 Bart van de Lint, Uwe Fechner
# SPDX-License-Identifier: LGPL-3.0-only

# This example demonstrates how to export a system structure to a YAML file
# in the same component-based format that load_sys_struct_from_yaml can read.

using Pkg
if Base.active_project() != joinpath(@__DIR__, "Project.toml")
    Pkg.activate(@__DIR__)
end

using GLMakie
using SymbolicAWEModels
import SymbolicAWEModels: Point  # resolve ambiguity with GLMakie

set_data_path(joinpath(dirname(@__DIR__), "data"))
set = Settings("base/system.yaml")
set.v_wind = 0.0

horizontal_span = 8.0
n_segments = 10
total_length = 10.0

# Segment properties (soft rope, 4mm diameter)
seg_stiffness = 500.0  # EA [N]
seg_damping   = 50.0   # [N·s]
seg_diameter  = 0.004  # [m]

points = [Point(:left, [0, 0, 5], STATIC)]
for i in 1:n_segments - 1
    x = i * horizontal_span / n_segments
    push!(points, Point(Symbol("p$i"), [x, 0, 5], DYNAMIC;
        world_frame_damping=1.0))
end
push!(points, Point(:right, [horizontal_span, 0, 5], STATIC))

l0_seg = total_length / n_segments
segments = Segment[]
for i in 1:n_segments
    push!(segments, Segment(
        Symbol("seg$i"), points[i].name,
        points[i + 1].name,
        seg_stiffness, seg_damping, seg_diameter;
        l0=l0_seg, compression_frac=0.01))
end

sys = SystemStructure("catenary", set; points, segments)

# Export to YAML file
save_sys_struct_to_yaml(sys, joinpath(dirname(@__DIR__), "data", "catenary_export.yaml"))

# Verify by loading it back
sys2 = load_sys_struct_from_yaml(
    joinpath(dirname(@__DIR__), "data", "catenary_export.yaml");
    system_name="catenary_reloaded", set=set)

println("=== Original system: $(sys.name) ===")
println("  Points: $(length(sys.points)), Segments: $(length(sys.segments))")
println("=== Reloaded system: $(sys2.name) ===")
println("  Points: $(length(sys2.points)), Segments: $(length(sys2.segments))")
println("Export/import round-trip successful!")

# Visualize
scene = plot(sys2)