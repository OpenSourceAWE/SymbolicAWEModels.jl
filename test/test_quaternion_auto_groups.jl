# Copyright (c) 2025 Bart van de Lint, Jelle Poland
# SPDX-License-Identifier: MPL-2.0

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

# Test auto-creation of groups for QUATERNION wings
using SymbolicAWEModels
using SymbolicAWEModels: VortexStepMethod, WING,
    QUATERNION, REFINE
using Test
using LinearAlgebra

@testset "QUATERNION wing auto-group creation" begin
    # Copy 2plate_kite data to temp directory
    pkg_root = dirname(@__DIR__)
    src_data_path = joinpath(
        pkg_root, "data", "2plate_kite"
    )
    tmpdir = mktempdir()
    data_path = joinpath(tmpdir, "2plate_kite")
    cp(src_data_path, data_path; force=true)
    set_data_path(data_path)

    struc_yaml = joinpath(
        data_path, "quat_struc_geometry.yaml")
    refine_yaml = joinpath(
        data_path, "refine_struc_geometry.yaml")

    set = Settings("system.yaml")
    vsm_set_path = joinpath(
        data_path, "vsm_settings.yaml")
    vsm_set = VortexStepMethod.VSMSettings(
        vsm_set_path; data_prefix=false)

    # ── REFINE: should have 0 groups ──────────────
    sys_refine = load_sys_struct_from_yaml(
        refine_yaml;
        system_name="2plate_refine", set, vsm_set)

    @test length(sys_refine.wings) == 1
    @test sys_refine.wings[1].wing_type == REFINE
    @test length(sys_refine.groups) == 0
    @test length(sys_refine.wings[1].group_idxs) == 0

    # ── QUATERNION with YAML-defined groups ───────
    # quat_struc_geometry.yaml has 3 explicit groups
    # and 7 WING points (6 LE/TE + kcu).
    sys_quat = load_sys_struct_from_yaml(
        struc_yaml;
        system_name="2plate_quat", set, vsm_set,
        wing_type=QUATERNION)

    wing = sys_quat.wings[1]
    @test wing.wing_type == QUATERNION
    @test length(sys_quat.groups) == 3
    @test length(wing.group_idxs) == 3
    @test !isnothing(wing.wing_segments)
    @test length(wing.wing_segments) == 3

    # Geometry was computed from closest VSM panel
    for group in sys_quat.groups
        @test !iszero(group.chord)
        @test !iszero(group.y_airf)
    end

    # ── QUATERNION auto-group creation ────────────
    # Load without explicit groups: auto-group should
    # kick in for the 6 LE/TE WING points (kcu is
    # excluded because it isn't in any group → the
    # group-based path is unavailable and the fallback
    # consecutive-pair heuristic uses only even-count
    # subsets).
    # Easiest approach: load REFINE YAML as QUATERNION
    # (REFINE YAML has 6 WING points, no groups, no
    # kcu WING point).
    sys_auto = load_sys_struct_from_yaml(
        refine_yaml;
        system_name="2plate_auto", set, vsm_set,
        wing_type=QUATERNION)

    wing_auto = sys_auto.wings[1]
    @test wing_auto.wing_type == QUATERNION
    @test length(sys_auto.groups) == 3
    @test length(wing_auto.group_idxs) == 3
    @test !isnothing(wing_auto.wing_segments)
    @test length(wing_auto.wing_segments) == 3

    for group in sys_auto.groups
        @test !iszero(group.chord)
        @test !iszero(group.y_airf)
    end
end
nothing
