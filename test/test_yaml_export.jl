# SPDX-FileCopyrightText: 2025 Uwe Fechner
# SPDX-License-Identifier: LGPL-3.0-only

# test_yaml_export.jl - Tests for save_sys_struct_to_yaml()
#
# Verifies that:
# 1. A simple catenary system exports and round-trips correctly
# 2. The output uses compact YAML format (inline header arrays)
# 3. Complex kite systems (with wings, tethers, winches, etc.)
#    export and round-trip correctly
# 4. All structural sections round-trip with correct counts

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

using Test
using SymbolicAWEModels
using SymbolicAWEModels: DYNAMIC, STATIC, PARTICLE_DYNAMICS, save_sys_struct_to_yaml, load_sys_struct_from_yaml
using KiteUtils
using LinearAlgebra

# ==================== Helper ==================== #

"""
    count_lines_starting_with(io::IOBuffer, prefix::String) -> Int

Count lines in a buffer that start with a given prefix (after stripping).
"""
function count_lines_starting_with(io::IOBuffer, prefix::String)
    seekstart(io)
    count = 0
    for line in eachline(io)
        stripped = strip(line)
        if startswith(stripped, prefix)
            count += 1
        end
    end
    return count
end

# ==================== Catenary (simple) ==================== #

@testset "YAML Export — Catenary round-trip" begin
    tmpdir = mktempdir()
    yaml_path = joinpath(tmpdir, "catenary.yaml")
    pkg_root = dirname(@__DIR__)

    # Build a simple catenary system
    horizontal_span = 8.0
    n_segments = 10
    total_length = 10.0
    seg_stiffness = 500.0
    seg_damping = 50.0
    seg_diameter = 0.004

    points = [SymbolicAWEModels.Point(:left, [0, 0, 5], STATIC)]
    for i in 1:(n_segments - 1)
        x = i * horizontal_span / n_segments
        push!(points, SymbolicAWEModels.Point(
            Symbol("p$i"), [x, 0, 5], DYNAMIC;
            world_frame_damping=1.0))
    end
    push!(points, SymbolicAWEModels.Point(
        :right, [horizontal_span, 0, 5], STATIC))

    l0_seg = total_length / n_segments
    segments = SymbolicAWEModels.Segment[]
    for i in 1:n_segments
        push!(segments, SymbolicAWEModels.Segment(
            Symbol("seg$i"), points[i].name, points[i + 1].name,
            seg_stiffness, seg_damping, seg_diameter;
            l0=l0_seg, compression_frac=0.01))
    end

    set = Settings(joinpath(pkg_root, "data", "base", "system.yaml"))
    set.v_wind = 0.0
    sys = SymbolicAWEModels.SystemStructure(
        "catenary", set; points, segments)

    # Export
    @test save_sys_struct_to_yaml(sys, yaml_path) === nothing
    @test isfile(yaml_path)

    # Reload
    sys2 = load_sys_struct_from_yaml(
        yaml_path; system_name="reloaded", set=set)

    # Compare counts
    @test length(sys2.points) == length(sys.points)
    @test length(sys2.segments) == length(sys.segments)
    @test length(sys2.wings) == 0
    @test length(sys2.winches) == 0
    @test length(sys2.pulleys) == 0
    @test length(sys2.tethers) == 0
    @test length(sys2.transforms) == 0
    @test length(sys2.twist_surfaces) == 0

    # Compare point names and types
    for i in eachindex(sys.points)
        p1 = sys.points[i]
        p2 = sys2.points[i]
        @test p1.name == p2.name
        @test p1.type == p2.type
        @test p1.pos_cad ≈ p2.pos_cad
    end

    # Compare segment names and endpoint references
    for i in eachindex(sys.segments)
        s1 = sys.segments[i]
        s2 = sys2.segments[i]
        @test s1.name == s2.name
        @test s1.unit_stiffness == s2.unit_stiffness
        @test s1.unit_damping == s2.unit_damping
        @test s1.diameter == s2.diameter
        @test s1.l0 == s2.l0
        @test s1.compression_frac == s2.compression_frac
    end

    # Verify compact YAML format
    yaml_content = read(yaml_path, String)
    @test occursin("headers: [name, pos_cad, type", yaml_content)
    @test occursin("headers: [name, point_i, point_j, l0", yaml_content)
    # Should NOT have expanded header format
    @test !occursin("headers:\n    - \"name\"", yaml_content)
end

# ==================== 2-plate kite (complex) ==================== #

@testset "YAML Export — 2-plate kite round-trip" begin
    pkg_root = dirname(@__DIR__)
    src_data_path = joinpath(pkg_root, "data", "2plate_kite")

    tmpdir = mktempdir()
    data_path = joinpath(tmpdir, "2plate_kite")
    cp(src_data_path, data_path; force=true)

    set_data_path(data_path)
    set = Settings("system.yaml")

    vsm_settings_path = joinpath(data_path, "vsm_settings.yaml")
    vsm_set = SymbolicAWEModels.VortexStepMethod.VSMSettings(
        vsm_settings_path; data_prefix=false)

    # --- Test RIGID_DYNAMICS (quaternion) ---
    @testset "RIGID_DYNAMICS" begin
        rigid_yaml = joinpath(data_path,
            "rigid_structural_geometry.yaml")
        sys = load_sys_struct_from_yaml(
            rigid_yaml; set=set, vsm_set=vsm_set)

        export_path = joinpath(tmpdir, "rigid_export.yaml")
        save_sys_struct_to_yaml(sys, export_path)

        sys2 = load_sys_struct_from_yaml(
            export_path; system_name="reloaded",
            set=set, vsm_set=vsm_set)

        @test length(sys2.points) == length(sys.points)
        @test length(sys2.segments) == length(sys.segments)
        @test length(sys2.wings) == length(sys.wings)
        @test length(sys2.winches) == length(sys.winches)
        @test length(sys2.pulleys) == length(sys.pulleys)
        @test length(sys2.tethers) == length(sys.tethers)
        @test length(sys2.transforms) == length(sys.transforms)
        @test length(sys2.twist_surfaces) ==
              length(sys.twist_surfaces)

        # Compare wing properties
        for i in eachindex(sys.wings)
            w1 = sys.wings[i]
            w2 = sys2.wings[i]
            @test w1.name == w2.name
            @test w1.dynamics_type == w2.dynamics_type
            @test length(w1.twist_surface_idxs) ==
                  length(w2.twist_surface_idxs)
        end

        # Compare transform properties
        for i in eachindex(sys.transforms)
            t1 = sys.transforms[i]
            t2 = sys2.transforms[i]
            @test t1.name == t2.name
            @test t1.elevation ≈ t2.elevation
            @test t1.azimuth ≈ t2.azimuth
            @test t1.heading ≈ t2.heading
        end

        # Compare point names
        for i in eachindex(sys.points)
            @test sys.points[i].name == sys2.points[i].name
        end

        # Verify compact format
        yaml_content = read(export_path, String)
        @test occursin("headers: [name, start_point", yaml_content)
        @test occursin("aero_mode:", yaml_content)
        @test occursin("dynamics_type:", yaml_content)
        @test occursin("origin_idx:", yaml_content)
        @test occursin("z_ref_points:", yaml_content)
        @test occursin("y_ref_points:", yaml_content)
        @test occursin("base_pos:", yaml_content)
    end

    # --- Test PARTICLE_DYNAMICS (refine/particle) ---
    @testset "PARTICLE_DYNAMICS" begin
        particle_yaml = joinpath(data_path,
            "particle_structural_geometry.yaml")
        sys = load_sys_struct_from_yaml(
            particle_yaml; set=set, vsm_set=vsm_set)

        export_path = joinpath(tmpdir, "particle_export.yaml")
        save_sys_struct_to_yaml(sys, export_path)

        sys2 = load_sys_struct_from_yaml(
            export_path; system_name="reloaded_particle",
            set=set, vsm_set=vsm_set)

        @test length(sys2.points) == length(sys.points)
        @test length(sys2.segments) == length(sys.segments)
        @test length(sys2.wings) == length(sys.wings)
        @test length(sys2.winches) == length(sys.winches)
        @test length(sys2.pulleys) == length(sys.pulleys)
        @test length(sys2.tethers) == length(sys.tethers)

        for (w1, w2) in zip(sys.wings, sys2.wings)
            @test w1.dynamics_type == w2.dynamics_type
            @test w1.dynamics_type == PARTICLE_DYNAMICS
        end

        # PARTICLE_DYNAMICS has no twist_surfaces in this config
        @test length(sys2.twist_surfaces) ==
              length(sys.twist_surfaces)

        # Verify dict-format sections
        yaml_content = read(export_path, String)
        @test occursin("name: main_wing", yaml_content)
        @test occursin("PARTICLE_DYNAMICS", yaml_content)
    end
end

# ==================== Error handling ==================== #

@testset "YAML Export — Error handling" begin
    tmpdir = mktempdir()

    # Export to non-existent directory (should fail)
    pkg_root = dirname(@__DIR__)
    set = Settings(joinpath(pkg_root, "data", "base", "system.yaml"))
    p = [SymbolicAWEModels.Point(:test, [0,0,0], STATIC)]
    sys = SymbolicAWEModels.SystemStructure("test", set; points=p)

    bad_path = joinpath(tmpdir, "nonexistent", "out.yaml")
    @test_throws SystemError save_sys_struct_to_yaml(sys, bad_path)
end

# ==================== Empty systems ==================== #

@testset "YAML Export — Minimal system" begin
    tmpdir = mktempdir()
    pkg_root = dirname(@__DIR__)

    # System with just a single point
    set = Settings(joinpath(pkg_root, "data", "base", "system.yaml"))
    points = [SymbolicAWEModels.Point(:origin, [0,0,0], STATIC)]
    sys = SymbolicAWEModels.SystemStructure(
        "minimal", set; points=points)

    export_path = joinpath(tmpdir, "minimal.yaml")
    save_sys_struct_to_yaml(sys, export_path)

    sys2 = load_sys_struct_from_yaml(
        export_path; system_name="minimal2", set=set)

    @test length(sys2.points) == 1
    @test sys2.points[1].name == :origin
    @test length(sys2.segments) == 0
end