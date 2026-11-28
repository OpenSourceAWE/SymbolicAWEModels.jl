# Copyright (c) 2025 Bart van de Lint, Jelle Poland
# SPDX-License-Identifier: LGPL-3.0-only

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

# Test station handling for RIGID_DYNAMICS wings (explicit only)
using SymbolicAWEModels
using SymbolicAWEModels: VortexStepMethod,
    RIGID_DYNAMICS, PARTICLE_DYNAMICS, wing_body, wing_station_idxs
using Test
using LinearAlgebra

@testset "RIGID_DYNAMICS wing station handling" begin
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

    # ── PARTICLE_DYNAMICS: 3 explicit stations (LE/TE membership) ──
    sys_refine, _ = load_authoring(
        refine_yaml;
        name="2plate_refine", set)

    @test length(sys_refine.wings) == 1
    refine_sam = SymbolicAWEModel(set, sys_refine; vsm_set)
    @test body_dynamics_type(refine_sam, wing_body(refine_sam, 1)) == PARTICLE_DYNAMICS
    @test length(sys_refine.spanwise_stations) == 3
    @test length(wing_station_idxs(refine_sam, 1)) == 3

    # ── RIGID_DYNAMICS with YAML-defined stations ───────
    # rigid_structural_geometry.yaml has 3 explicit stations
    # and 7 WING points (6 LE/TE + kcu).
    sys_quat, _ = load_authoring(
        struc_yaml;
        name="2plate_quat", set)

    sam = SymbolicAWEModel(set, sys_quat; vsm_set)
    @test body_dynamics_type(sam, wing_body(sam, 1)) == RIGID_DYNAMICS
    @test length(sys_quat.spanwise_stations) == 3
    @test length(wing_station_idxs(sam, 1)) == 3
    @test !isnothing(sam.aero[1].wing_segments)
    @test length(sam.aero[1].wing_segments) == 3

    # Geometry was computed from closest VSM panel
    for station_idx in eachindex(sam.definition.spanwise_stations)
        @test !iszero(run_value(sam, :spanwise_stations, station_idx, :chord))
        @test !iszero(run_value(sam, :spanwise_stations, station_idx, :y_airf))
    end

end
nothing
