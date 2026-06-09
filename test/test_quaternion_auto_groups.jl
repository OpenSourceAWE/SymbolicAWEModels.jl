# Copyright (c) 2025 Bart van de Lint, Jelle Poland
# SPDX-License-Identifier: LGPL-3.0-only

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

# Test auto-creation of twist_surfaces for RIGID_DYNAMICS wings
using SymbolicAWEModels
using SymbolicAWEModels: VortexStepMethod, WING,
    RIGID_DYNAMICS, PARTICLE_DYNAMICS
using Test
using LinearAlgebra

@testset "RIGID_DYNAMICS wing auto-twist_surface creation" begin
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
        data_path, "rigid_structural_geometry.yaml")
    refine_yaml = joinpath(
        data_path, "particle_structural_geometry.yaml")

    set = Settings("system.yaml")
    vsm_set_path = joinpath(
        data_path, "vsm_settings.yaml")
    vsm_set = VortexStepMethod.VSMSettings(
        vsm_set_path; data_prefix=false)

    # ── PARTICLE_DYNAMICS: should have 0 twist_surfaces ──────────────
    sys_refine = load_sys_struct_from_yaml(
        refine_yaml;
        system_name="2plate_refine", set, vsm_set)

    @test length(sys_refine.wings) == 1
    @test sys_refine.wings[1].dynamics_type == PARTICLE_DYNAMICS
    @test length(sys_refine.twist_surfaces) == 0
    @test length(sys_refine.wings[1].twist_surface_idxs) == 0

    # ── RIGID_DYNAMICS with YAML-defined twist_surfaces ───────
    # rigid_structural_geometry.yaml has 3 explicit twist_surfaces
    # and 7 WING points (6 LE/TE + kcu).
    sys_quat = load_sys_struct_from_yaml(
        struc_yaml;
        system_name="2plate_quat", set, vsm_set,
        dynamics_type=RIGID_DYNAMICS)

    wing = sys_quat.wings[1]
    @test wing.dynamics_type == RIGID_DYNAMICS
    @test length(sys_quat.twist_surfaces) == 3
    @test length(wing.twist_surface_idxs) == 3
    @test !isnothing(wing.wing_segments)
    @test length(wing.wing_segments) == 3

    # Geometry was computed from closest VSM panel
    for twist_surface in sys_quat.twist_surfaces
        @test !iszero(twist_surface.chord)
        @test !iszero(twist_surface.y_airf)
    end

    # ── RIGID_DYNAMICS auto-twist_surface creation ────────────
    # Load without explicit twist_surfaces: auto-twist_surface should
    # kick in for the 6 LE/TE WING points (kcu is
    # excluded because it isn't in any twist_surface → the
    # twist_surface-based path is unavailable and the fallback
    # consecutive-pair heuristic uses only even-count
    # subsets).
    # Easiest approach: load PARTICLE_DYNAMICS YAML as RIGID_DYNAMICS
    # (PARTICLE_DYNAMICS YAML has 6 WING points, no twist_surfaces, no
    # kcu WING point).
    sys_auto = load_sys_struct_from_yaml(
        refine_yaml;
        system_name="2plate_auto", set, vsm_set,
        dynamics_type=RIGID_DYNAMICS)

    wing_auto = sys_auto.wings[1]
    @test wing_auto.dynamics_type == RIGID_DYNAMICS
    @test length(sys_auto.twist_surfaces) == 3
    @test length(wing_auto.twist_surface_idxs) == 3
    @test !isnothing(wing_auto.wing_segments)
    @test length(wing_auto.wing_segments) == 3

    for twist_surface in sys_auto.twist_surfaces
        @test !iszero(twist_surface.chord)
        @test !iszero(twist_surface.y_airf)
    end
end
nothing
