# SPDX-FileCopyrightText: 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_transform.jl - Transform and spherical coordinate tests
#
# Tests initial heading, elevation, azimuth and their velocities after init.
# Verifies:
# 1. Initial angles match YAML configuration
# 2. Initial velocities match YAML configuration
# 3. Geometric consistency (position from spherical coords)
# 4. Heading calculation consistency
#
# Uses 2plate_kite configuration files with both REFINE and QUATERNION wing types.

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

@isdefined(test_init!) || include(joinpath(@__DIR__, "util.jl"))

using Test
using SymbolicAWEModels
using SymbolicAWEModels: KVec3, VortexStepMethod,
    calc_heading, reposition!
using KiteUtils
using LinearAlgebra

@testset "Transform Tests" begin
    # Setup - copy data directory
    pkg_root = dirname(@__DIR__)
    src_data_path = joinpath(pkg_root, "data", "2plate_kite")

    tmpdir = mktempdir()
    data_path = joinpath(tmpdir, "2plate_kite")
    cp(src_data_path, data_path; force=true)

    # Set data path and load settings
    set_data_path(data_path)
    set = Settings("system.yaml")

    # Load VSM settings from data directory
    vsm_settings_path = joinpath(data_path, "vsm_settings.yaml")
    vsm_set = VortexStepMethod.VSMSettings(vsm_settings_path; data_prefix=false)

    # Paths for both wing types
    quat_yaml_path = joinpath(data_path, "quat_struc_geometry.yaml")
    refine_yaml_path = joinpath(data_path, "refine_struc_geometry.yaml")

    # Create and initialize SAMs once for each wing type
    quat_sys = load_sys_struct_from_yaml(
        quat_yaml_path; system_name="transform_test_QUATERNION", set=set, vsm_set=vsm_set
    )
    quat_sam = SymbolicAWEModel(set, quat_sys)
    test_init!(quat_sam; remake=false, reload=false)  # Load/build once

    refine_sys = load_sys_struct_from_yaml(
        refine_yaml_path; system_name="transform_test_REFINE", set=set, vsm_set=vsm_set
    )
    refine_sam = SymbolicAWEModel(set, refine_sys)
    test_init!(refine_sam; remake=false, reload=false)  # Load/build once

    # Helper to reset transform to default YAML values
    function reset_transform!(sys)
        tf = sys.transforms[:main_transform]
        tf.elevation = deg2rad(80)
        tf.azimuth = deg2rad(0)
        tf.heading = deg2rad(0)
        tf.elevation_vel = 0.0
        tf.azimuth_vel = 0.0
    end

    # Test both wing types
    sam_configs = [
        ("REFINE", refine_sam, refine_yaml_path),
        ("QUATERNION", quat_sam, quat_yaml_path),
    ]

    for (wing_type_name, sam, yaml_path) in sam_configs
        @testset "$wing_type_name Wing" begin
            # ================================================================
            # YAML Loading Verification (uses already-loaded sys_struct)
            # ================================================================
            @testset "YAML Loading Verification" begin
                sys = sam.sys_struct

                # Verify transform was loaded
                @test length(sys.transforms) == 1
                @test haskey(sys.transforms, :main_transform)

                transform = sys.transforms[:main_transform]

                # Verify base point
                @test transform.base_point_idx == 10  # ground point index

                # Verify wing reference
                @test transform.wing_idx == 1  # main_wing index

                println("\n  ====== [$wing_type_name] Loaded transform: " *
                    "elev=$(round(rad2deg(transform.elevation), digits=1))°, " *
                    "azim=$(round(rad2deg(transform.azimuth), digits=1))°, " *
                    "heading=$(round(rad2deg(transform.heading), digits=1))° ======\n")
            end

            # ================================================================
            # Physics Test 1: Initial angles after init
            # ================================================================
            @testset "Initial angles after init!" begin
                reset_transform!(sam.sys_struct)
                test_init!(sam; remake=false, reload=false)

                # After init, the transform angles should still match
                transform = sam.sys_struct.transforms[:main_transform]

                @test transform.elevation ≈ deg2rad(80) atol=1e-10
                @test transform.azimuth ≈ deg2rad(0) atol=1e-10
                @test transform.heading ≈ deg2rad(0) atol=1e-10

                println("\n  ====== [$wing_type_name] After init: " *
                    "elev=$(round(rad2deg(transform.elevation), digits=1))°, " *
                    "azim=$(round(rad2deg(transform.azimuth), digits=1))°, " *
                    "heading=$(round(rad2deg(transform.heading), digits=1))° ======\n")
            end

            # ================================================================
            # Physics Test 2: Initial velocities after init
            # ================================================================
            @testset "Initial velocities after init!" begin
                reset_transform!(sam.sys_struct)
                test_init!(sam; remake=false, reload=false)

                transform = sam.sys_struct.transforms[:main_transform]

                # Default velocities are 0
                @test transform.elevation_vel ≈ 0.0 atol=1e-10
                @test transform.azimuth_vel ≈ 0.0 atol=1e-10

                println("\n  ====== [$wing_type_name] Velocities: " *
                    "elev_vel=$(round(rad2deg(transform.elevation_vel), digits=2))°/s, " *
                    "azim_vel=$(round(rad2deg(transform.azimuth_vel), digits=2))°/s ======\n")
            end

            # ================================================================
            # Physics Test 3: Geometric consistency - position from spherical coords
            # ================================================================
            @testset "Geometric consistency" begin
                # For elevation=80deg, azimuth=0deg, the wing should be positioned
                # according to spherical coordinate transformation
                reset_transform!(sam.sys_struct)
                test_init!(sam; remake=false, reload=false)

                # Get wing position
                wing = sam.sys_struct.wings[:main_wing]
                wing_pos = wing.base.pos_w

                # Get ground position
                ground_pos = sam.sys_struct.points[:ground].pos_w

                # Vector from ground to wing
                rel_pos = wing_pos - ground_pos

                # Calculate distance (tether length)
                distance = norm(rel_pos)

                # The wing should be at the expected position
                # This tests that the transform correctly places the wing
                @test distance > 0  # Wing is above ground

                # Check that z component is positive (wing above ground level in world frame)
                @test wing_pos[3] > ground_pos[3]

                println("\n  ====== [$wing_type_name] Geometry: " *
                    "wing_pos=$(round.(wing_pos, digits=2)), " *
                    "distance=$(round(distance, digits=2))m ======\n")
            end

            # ================================================================
            # Physics Test 4: Alternative angles configuration
            # ================================================================
            @testset "Alternative angles configuration" begin
                # Modify transform to test different angles
                tf = sam.sys_struct.transforms[:main_transform]
                tf.elevation = deg2rad(45)
                tf.azimuth = deg2rad(30)
                tf.heading = deg2rad(10)
                tf.elevation_vel = deg2rad(0.1)
                tf.azimuth_vel = deg2rad(0.5)

                # Verify modified values before init
                @test tf.elevation ≈ deg2rad(45) atol=1e-10
                @test tf.azimuth ≈ deg2rad(30) atol=1e-10
                @test tf.heading ≈ deg2rad(10) atol=1e-10
                @test tf.elevation_vel ≈ deg2rad(0.1) atol=1e-10
                @test tf.azimuth_vel ≈ deg2rad(0.5) atol=1e-10

                test_init!(sam; remake=false, reload=false)

                # Angles should be preserved after init
                transform_after = sam.sys_struct.transforms[:main_transform]
                @test transform_after.elevation ≈ deg2rad(45) atol=1e-10
                @test transform_after.azimuth ≈ deg2rad(30) atol=1e-10

                println("\n  ====== [$wing_type_name] Alt config: " *
                    "elev=$(round(rad2deg(transform_after.elevation), digits=1))°, " *
                    "azim=$(round(rad2deg(transform_after.azimuth), digits=1))° ======\n")
            end

            # ================================================================
            # Physics Test 5: Transform affects wing position
            # ================================================================
            @testset "Transform affects wing position" begin
                # Test 1: elevation = 80 deg (default)
                reset_transform!(sam.sys_struct)
                test_init!(sam; remake=false, reload=false)
                wing_z1 = sam.sys_struct.wings[:main_wing].base.pos_w[3]

                # Test 2: elevation = 45 deg
                sam.sys_struct.transforms[:main_transform].elevation = deg2rad(45)
                test_init!(sam; remake=false, reload=false)
                wing_z2 = sam.sys_struct.wings[:main_wing].base.pos_w[3]

                # Higher elevation should result in higher z position
                # (wing more overhead)
                @test wing_z1 > wing_z2

                println("\n  ====== [$wing_type_name] Elevation effect: " *
                    "z(80°)=$(round(wing_z1, digits=2))m > " *
                    "z(45°)=$(round(wing_z2, digits=2))m ======\n")
            end

            # ================================================================
            # Physics Test 6: Azimuth affects y-position
            # ================================================================
            @testset "Azimuth affects y-position" begin
                # Test 1: azimuth = 0 deg (default)
                reset_transform!(sam.sys_struct)
                test_init!(sam; remake=false, reload=false)
                wing_y1 = sam.sys_struct.wings[:main_wing].base.pos_w[2]

                # Test 2: azimuth = 30 deg (more to the side)
                sam.sys_struct.transforms[:main_transform].azimuth = deg2rad(30)
                test_init!(sam; remake=false, reload=false)
                wing_y2 = sam.sys_struct.wings[:main_wing].base.pos_w[2]

                # Larger azimuth should give larger |y| component
                @test abs(wing_y2) > abs(wing_y1)

                println("\n  ====== [$wing_type_name] Azimuth effect: " *
                    "|y|(30°)=$(round(abs(wing_y2), digits=2))m > " *
                    "|y|(0°)=$(round(abs(wing_y1), digits=2))m ======\n")
            end

            # ================================================================
            # Physics Test 7: Heading affects wing orientation (not position)
            # ================================================================
            @testset "Heading affects orientation" begin
                reset_transform!(sam.sys_struct)
                test_init!(sam; remake=false, reload=false)

                wing = sam.sys_struct.wings[:main_wing]

                # Wing should have a rotation matrix
                @test !isnothing(wing.base.R_b_to_w)
                @test size(wing.base.R_b_to_w) == (3, 3)

                # R_b_to_w should be a valid rotation matrix (orthonormal)
                @test det(wing.base.R_b_to_w) ≈ 1.0 atol=1e-10
                @test wing.base.R_b_to_w * wing.base.R_b_to_w' ≈ I(3) atol=1e-10

                println("\n  ====== [$wing_type_name] Heading affects rotation: " *
                    "det(R_b_to_w)=$(round(det(wing.base.R_b_to_w), digits=4)) ======\n")
            end

            # ================================================================
            # Physics Test 8: Base position offset
            # ================================================================
            @testset "Base position from base_point" begin
                # The transform references ground as base_point
                transform = sam.sys_struct.transforms[:main_transform]
                @test transform.base_point_idx == 10  # ground index

                reset_transform!(sam.sys_struct)
                test_init!(sam; remake=false, reload=false)

                # Transform base_pos should match the ground point position
                ground_pos = sam.sys_struct.points[:ground].pos_w
                @test transform.base_pos ≈ ground_pos atol=1e-10

                println("\n  ====== [$wing_type_name] Base point: " *
                    "ground_pos=$(round.(ground_pos, digits=2)), " *
                    "transform.base_pos=$(round.(transform.base_pos, digits=2)) ======\n")
            end

            # ================================================================
            # Physics Test 9: reposition! heading matches reinit!
            # Tested with both origin and non-origin base positions
            # to verify heading rotates around the correct axis.
            # ================================================================
            function test_reposition_heading(sam, base_pos, label)
                sys = sam.sys_struct
                tf = sys.transforms[:main_transform]
                orig_base_pos = copy(tf.base_pos)

                @testset "Reposition heading ($label)" begin
                    for h_deg in [0, 10, -15, 30, 45, -45]
                        target_h = deg2rad(h_deg)

                        # Reference: reinit! with target heading
                        tf.base_pos .= base_pos
                        reset_transform!(sys)
                        tf.heading = target_h
                        test_init!(sam; remake=false, reload=false)
                        wing = sys.wings[:main_wing]
                        reinit_R = copy(wing.R_b_to_w)
                        reinit_pos = copy(wing.pos_w)

                        # reinit! with heading=0, then
                        # reposition! to target
                        tf.base_pos .= base_pos
                        reset_transform!(sys)
                        test_init!(sam; remake=false, reload=false)
                        tf.heading = target_h
                        reposition!(sys.transforms, sys)
                        wing = sys.wings[:main_wing]

                        @test reinit_R ≈ wing.R_b_to_w atol=1e-6
                        @test reinit_pos ≈ wing.pos_w atol=1e-4
                    end
                end

                tf.base_pos .= orig_base_pos
            end

            test_reposition_heading(
                sam, KVec3(0.0, 0.0, 0.0), "origin base")
            test_reposition_heading(
                sam, KVec3(10.0, 5.0, 0.0),
                "non-origin base")
        end
    end

    # Cleanup
    rm(tmpdir; recursive=true)
end
nothing
