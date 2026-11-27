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
# Uses 2plate_kite configuration files with both PARTICLE_DYNAMICS and RIGID_DYNAMICS wing types.

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

@isdefined(test_init!) || include(joinpath(@__DIR__, "util.jl"))

using Test
using SymbolicAWEModels
using SymbolicAWEModels: KVec3, VortexStepMethod, quaternion_to_rotation_matrix,
    write_wing_scalars!
using KiteUtils
using SymbolicAWEModels.KiteGeometry: base_position
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
    rigid_dynamics_yaml_path = joinpath(data_path, "rigid_structural_geometry.yaml")
    particle_dynamics_yaml_path = joinpath(data_path, "particle_structural_geometry.yaml")

    # Create and initialize SAMs once for each wing type
    rigid_dynamics_sys, rigid_dynamics_transforms = load_sys_struct_from_yaml(
        rigid_dynamics_yaml_path; system_name="transform_test_RIGID_DYNAMICS", set=set, vsm_set=vsm_set,
        aero_mode=AeroNone()
    )
    rigid_dynamics_sam = SymbolicAWEModel(set, rigid_dynamics_sys;
                                          transforms=rigid_dynamics_transforms)
    test_init!(rigid_dynamics_sam)

    particle_dynamics_sys, particle_dynamics_transforms = load_sys_struct_from_yaml(
        particle_dynamics_yaml_path; system_name="transform_test_PARTICLE_DYNAMICS", set=set, vsm_set=vsm_set,
        aero_mode=AeroNone()
    )
    particle_dynamics_sam = SymbolicAWEModel(set, particle_dynamics_sys;
                                             transforms=particle_dynamics_transforms)
    test_init!(particle_dynamics_sam)

    # Helper to reset transform to default YAML values
    function reset_transform!(sys)
        tf = sys.transforms[:main_transform]
        tf.elevation = deg2rad(80)
        tf.azimuth = deg2rad(0)
        tf.heading = deg2rad(0)
        tf.elevation_vel = 0.0
        tf.azimuth_vel = 0.0
        tf.turn_rate = 0.0
    end

    # Each member's velocity against each member's position from base_pos: unchanged by
    # a rigid rotation of positions and velocities about base_pos.
    function motion_against_pose(members, base_pos)
        return [dot(member.vel_w, other.pos_w .- base_pos)
                for member in members, other in members]
    end

    # Give every member of `transform` in `state` a distinct velocity and return them.
    function set_member_velocities!(state, transform)
        members = [state.points[transform.points]; state.bodies[transform.bodies]]
        for (i, member) in enumerate(members)
            member.vel_w .= [0.3i, -0.2i, 0.1i]
        end
        return members
    end

    # Test both wing types
    sam_configs = [
        ("PARTICLE_DYNAMICS", particle_dynamics_sam, particle_dynamics_yaml_path),
        ("RIGID_DYNAMICS", rigid_dynamics_sam, rigid_dynamics_yaml_path),
    ]

    for (dynamics_type_name, sam, yaml_path) in sam_configs
        @testset "$dynamics_type_name Wing" begin
            # ================================================================
            # YAML Loading Verification (uses already-loaded sys_struct)
            # ================================================================
            @testset "YAML Loading Verification" begin
                # Verify transform was loaded
                @test length(sam.transforms) == 1
                @test haskey(sam.transforms, :main_transform)

                transform = sam.transforms[:main_transform]

                # Verify base point
                @test transform.base_point == 10  # ground point index

                # Verify wing reference
                @test transform.body == 1  # main_wing index

                println("\n  ====== [$dynamics_type_name] Loaded transform: " *
                    "elev=$(round(rad2deg(transform.elevation), digits=1))°, " *
                    "azim=$(round(rad2deg(transform.azimuth), digits=1))°, " *
                    "heading=$(round(rad2deg(transform.heading), digits=1))° ======\n")
            end

            # ================================================================
            # Physics Test 1: Initial angles after init
            # ================================================================
            @testset "Initial angles after init!" begin
                reset_transform!(sam.state)
                place!(sam; prn=false)
                test_init!(sam; prn=false)

                # After init, the transform angles should still match
                transform = sam.transforms[:main_transform]

                @test transform.elevation ≈ deg2rad(80) atol=1e-10
                @test transform.azimuth ≈ deg2rad(0) atol=1e-10
                @test transform.heading ≈ deg2rad(0) atol=1e-10

                println("\n  ====== [$dynamics_type_name] After init: " *
                    "elev=$(round(rad2deg(transform.elevation), digits=1))°, " *
                    "azim=$(round(rad2deg(transform.azimuth), digits=1))°, " *
                    "heading=$(round(rad2deg(transform.heading), digits=1))° ======\n")
            end

            # ================================================================
            # Physics Test 2: Initial velocities after init
            # ================================================================
            @testset "Initial velocities after init!" begin
                reset_transform!(sam.state)
                place!(sam; prn=false)
                test_init!(sam; prn=false)

                transform = sam.transforms[:main_transform]

                # Default velocities are 0
                @test transform.elevation_vel ≈ 0.0 atol=1e-10
                @test transform.azimuth_vel ≈ 0.0 atol=1e-10

                println("\n  ====== [$dynamics_type_name] Velocities: " *
                    "elev_vel=$(round(rad2deg(transform.elevation_vel), digits=2))°/s, " *
                    "azim_vel=$(round(rad2deg(transform.azimuth_vel), digits=2))°/s ======\n")
            end

            # ================================================================
            # Physics Test 2b: A start velocity is the wing's own rate
            # ================================================================
            @testset "Start velocity round-trips into the wing rates" begin
                sys = sam.state
                tf = sys.transforms[:main_transform]
                reset_transform!(sys)
                tf.elevation = deg2rad(70)
                tf.azimuth = deg2rad(25)
                tf.elevation_vel = deg2rad(3.0)
                tf.azimuth_vel = deg2rad(-2.0)

                place!(sam)

                wing = sys.bodies[:main_wing]
                write_wing_scalars!(wing, sys.points;
                    base_point=tf.base_point)

                @test wing.elevation_vel ≈ tf.elevation_vel atol=1e-8
                @test wing.azimuth_vel ≈ tf.azimuth_vel atol=1e-8

                println("\n  ====== [$dynamics_type_name] Start velocity: " *
                    "elev_vel=$(round(rad2deg(wing.elevation_vel), digits=3))°/s, " *
                    "azim_vel=$(round(rad2deg(wing.azimuth_vel), digits=3))°/s ======\n")

                reset_transform!(sys)
            end

            # ================================================================
            # Physics Test 2c: turn_rate spins the structure as it is placed
            # ================================================================
            @testset "turn_rate advances the placed heading" begin
                sys = sam.state
                tf = sys.transforms[:main_transform]
                turn_rate = deg2rad(6.0)
                dt = 1e-5

                reset_transform!(sys)
                tf.elevation = deg2rad(70)
                tf.azimuth = deg2rad(25)
                tf.heading = deg2rad(15)
                tf.turn_rate = turn_rate
                place!(sam)
                point_pos = [copy(point.pos_w) for point in sys.points]
                point_vel = [copy(point.vel_w) for point in sys.points]
                body_pos = [copy(body.pos_w) for body in sys.bodies]
                body_vel = [copy(body.vel_w) for body in sys.bodies]

                tf.heading += turn_rate * dt
                tf.turn_rate = 0.0
                place!(sam)

                for (idx, point) in enumerate(sys.points)
                    @test (point.pos_w .- point_pos[idx]) ./ dt ≈ point_vel[idx] atol=1e-6
                end
                for (idx, body) in enumerate(sys.bodies)
                    @test (body.pos_w .- body_pos[idx]) ./ dt ≈ body_vel[idx] atol=1e-6
                end

                println("\n  ====== [$dynamics_type_name] Turn rate: " *
                    "$(round(rad2deg(turn_rate), digits=2))°/s moves the fastest " *
                    "point at $(round(maximum(norm, point_vel), digits=3)) m/s ======\n")

                reset_transform!(sys)
            end

            # ================================================================
            # Physics Test 3: Geometric consistency - position from spherical coords
            # ================================================================
            @testset "Geometric consistency" begin
                # For elevation=80deg, azimuth=0deg, the wing should be positioned
                # according to spherical coordinate transformation
                reset_transform!(sam.state)
                place!(sam; prn=false)
                test_init!(sam; prn=false)

                # Get wing position
                wing = sam.state.bodies[:main_wing]
                wing_pos = wing.pos_w

                # Get ground position
                ground_pos = sam.state.points[:ground].pos_w

                # Vector from ground to wing
                rel_pos = wing_pos - ground_pos

                # Calculate distance (tether length)
                distance = norm(rel_pos)

                # The wing should be at the expected position
                # This tests that the transform correctly places the wing
                @test distance > 0  # Wing is above ground

                # Check that z component is positive (wing above ground level in world frame)
                @test wing_pos[3] > ground_pos[3]

                println("\n  ====== [$dynamics_type_name] Geometry: " *
                    "wing_pos=$(round.(wing_pos, digits=2)), " *
                    "distance=$(round(distance, digits=2))m ======\n")
            end

            # ================================================================
            # Physics Test 4: Alternative angles configuration
            # ================================================================
            @testset "Alternative angles configuration" begin
                # Modify transform to test different angles
                tf = sam.transforms[:main_transform]
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

                place!(sam; prn=false)
                test_init!(sam; prn=false)

                # Angles should be preserved after init
                transform_after = sam.transforms[:main_transform]
                @test transform_after.elevation ≈ deg2rad(45) atol=1e-10
                @test transform_after.azimuth ≈ deg2rad(30) atol=1e-10

                println("\n  ====== [$dynamics_type_name] Alt config: " *
                    "elev=$(round(rad2deg(transform_after.elevation), digits=1))°, " *
                    "azim=$(round(rad2deg(transform_after.azimuth), digits=1))° ======\n")
            end

            # ================================================================
            # Physics Test 5: Transform affects wing position
            # ================================================================
            @testset "Transform affects wing position" begin
                # Test 1: elevation = 80 deg (default)
                reset_transform!(sam.state)
                place!(sam; prn=false)
                test_init!(sam; prn=false)
                wing_z1 = sam.state.bodies[:main_wing].pos_w[3]

                # Test 2: elevation = 45 deg; init! alone keeps the placement
                sam.transforms[:main_transform].elevation = deg2rad(45)
                test_init!(sam; prn=false)
                @test sam.state.bodies[:main_wing].pos_w[3] ≈ wing_z1 atol=1e-6
                place!(sam; prn=false)
                test_init!(sam; prn=false)
                wing_z2 = sam.state.bodies[:main_wing].pos_w[3]

                # Higher elevation should result in higher z position
                # (wing more overhead)
                @test wing_z1 > wing_z2

                println("\n  ====== [$dynamics_type_name] Elevation effect: " *
                    "z(80°)=$(round(wing_z1, digits=2))m > " *
                    "z(45°)=$(round(wing_z2, digits=2))m ======\n")
            end

            # ================================================================
            # Physics Test 6: Azimuth affects y-position
            # ================================================================
            @testset "Azimuth affects y-position" begin
                # Test 1: azimuth = 0 deg (default)
                reset_transform!(sam.state)
                place!(sam; prn=false)
                test_init!(sam; prn=false)
                wing_y1 = sam.state.bodies[:main_wing].pos_w[2]

                # Test 2: azimuth = 30 deg (more to the side)
                sam.transforms[:main_transform].azimuth = deg2rad(30)
                place!(sam; prn=false)
                test_init!(sam; prn=false)
                wing_y2 = sam.state.bodies[:main_wing].pos_w[2]

                # Larger azimuth should give larger |y| component
                @test abs(wing_y2) > abs(wing_y1)

                println("\n  ====== [$dynamics_type_name] Azimuth effect: " *
                    "|y|(30°)=$(round(abs(wing_y2), digits=2))m > " *
                    "|y|(0°)=$(round(abs(wing_y1), digits=2))m ======\n")
            end

            # ================================================================
            # Physics Test 7: Heading affects wing orientation (not position)
            # ================================================================
            @testset "Heading affects orientation" begin
                reset_transform!(sam.state)
                place!(sam; prn=false)
                test_init!(sam; prn=false)

                wing = sam.state.bodies[:main_wing]

                # Wing should have a rotation matrix
                @test !isnothing(wing.R_b_to_w)
                @test size(wing.R_b_to_w) == (3, 3)

                # R_b_to_w should be a valid rotation matrix (orthonormal)
                @test det(wing.R_b_to_w) ≈ 1.0 atol=1e-10
                @test wing.R_b_to_w * wing.R_b_to_w' ≈ I(3) atol=1e-10

                println("\n  ====== [$dynamics_type_name] Heading affects rotation: " *
                    "det(R_b_to_w)=$(round(det(wing.R_b_to_w), digits=4)) ======\n")
            end

            # ================================================================
            # Physics Test 8: Base position offset
            # ================================================================
            @testset "Base position from base_point" begin
                # The transform references ground as base_point
                transform = sam.transforms[:main_transform]
                @test transform.base_point == 10  # ground index

                reset_transform!(sam.state)
                place!(sam; prn=false)
                test_init!(sam; prn=false)

                # Transform base_pos should match the ground point position
                ground_pos = sam.state.points[:ground].pos_w
                @test transform.base_pos ≈ ground_pos atol=1e-10

                println("\n  ====== [$dynamics_type_name] Base point: " *
                    "ground_pos=$(round.(ground_pos, digits=2)), " *
                    "transform.base_pos=$(round.(transform.base_pos, digits=2)) ======\n")
            end

            # ================================================================
            # Physics Test 9: reposition! heading matches reinit!
            # Tested with both origin and non-origin base positions
            # to verify heading rotates around the correct axis.
            # ================================================================
            function test_reposition_heading(sam, base_pos, label)
                sys = sam.state
                tf = sys.transforms[:main_transform]
                orig_base_pos = copy(tf.base_pos)

                @testset "Reposition heading ($label)" begin
                    for h_deg in [0, 10, -15, 30, 45, -45]
                        target_h = deg2rad(h_deg)

                        # Reference: reinit! with target heading
                        tf.base_pos .= base_pos
                        reset_transform!(sys)
                        tf.heading = target_h
                        place!(sam; prn=false)
                        test_init!(sam; prn=false)
                        wing = sys.bodies[:main_wing]
                        reinit_R = copy(wing.R_b_to_w)
                        reinit_pos = copy(wing.pos_w)

                        # reinit! with heading=0, then
                        # reposition! to target
                        tf.base_pos .= base_pos
                        reset_transform!(sys)
                        place!(sam; prn=false)
                        test_init!(sam; prn=false)
                        tf.heading = target_h
                        reposition!(sys)
                        wing = sys.bodies[:main_wing]

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

            @testset "Reposition turns velocities with the pose" begin
                sys = sam.state
                tf = sys.transforms[:main_transform]
                reset_transform!(sys)
                place!(sam; prn=false)
                test_init!(sam; prn=false)
                members = set_member_velocities!(sys, tf)
                base_pos = copy(base_position(Pose(sys), sys.transforms, tf))
                before = motion_against_pose(members, base_pos)

                tf.heading = deg2rad(25)
                tf.elevation = deg2rad(60)
                reposition!(sys)

                @test motion_against_pose(members, base_pos) ≈ before atol=1e-8
            end
        end
    end

    # ================================================================
    # Chained Transform Tests (kps4 plate kite with a tilted child transform)
    # ================================================================
    @testset "Chained Transforms" begin
        kite_angle = deg2rad(3.83)
        set_c, sys_c, transforms_c = build_plate_kite(tmpdir; name="chained_test",
            azimuth=deg2rad(10.0), tilt=kite_angle)
        particles = KiteUtils.get_particles(set_c.height_k, set_c.h_bridle,
            set_c.width, set_c.m_k)
        pos_top, pos_right, pos_left = particles[4], particles[5], particles[6]
        sam_c = SymbolicAWEModel(set_c, sys_c; transforms=transforms_c)
        init!(sam_c)

        @testset "get_base_pos returns the parent's current and initial position" begin
            sys = sam_c.state
            tf_child = sys.transforms[:kite_tilt]
            design = Pose([KVec3(point.pos_ENU) for point in sys.points],
                          [KVec3(body.pos_ENU) for body in sys.bodies],
                          [quaternion_to_rotation_matrix(body.Q_KA_to_ENU)
                           for body in sys.bodies])
            base_pos = base_position(Pose(sys), sys.transforms, tf_child)
            curr_base_pos = base_position(design, sys.transforms, tf_child)
            @test base_pos == sys.wings[:plate_wing].pos_w
            @test curr_base_pos == sys.wings[:plate_wing].pos_ENU
        end

        @testset "Child points translated from the authoring geometry" begin
            sys = sam_c.state
            top = sys.points[:top]
            @test !(top.pos_w ≈ pos_top)
            # top should be far from origin (at tether length)
            @test norm(top.pos_w) > 50.0
            @test top.pos_ENU == top.pos_w
        end

        @testset "Child points near parent wing" begin
            sys = sam_c.state
            wing_pos = sys.wings[1].pos_w
            top_pos = sys.points[:top].pos_w
            right_pos = sys.points[:right].pos_w
            left_pos = sys.points[:left].pos_w

            # All child points should be near the wing
            # (within bridle length, ~30m)
            @test norm(top_pos - wing_pos) < 50.0
            @test norm(right_pos - wing_pos) < 50.0
            @test norm(left_pos - wing_pos) < 50.0
        end

        @testset "Distances preserved (rigid body)" begin
            sys = sam_c.state
            # Authoring distances between child points
            cad_dist_top_right = norm(pos_top - pos_right)
            cad_dist_top_left = norm(pos_top - pos_left)
            cad_dist_right_left = norm(pos_right - pos_left)

            # World distances should match CAD distances
            # (transforms are rigid body)
            world_dist_top_right = norm(
                sys.points[:top].pos_w -
                sys.points[:right].pos_w)
            world_dist_top_left = norm(
                sys.points[:top].pos_w -
                sys.points[:left].pos_w)
            world_dist_right_left = norm(
                sys.points[:right].pos_w -
                sys.points[:left].pos_w)

            @test cad_dist_top_right ≈
                world_dist_top_right atol=1e-6
            @test cad_dist_top_left ≈
                world_dist_top_left atol=1e-6
            @test cad_dist_right_left ≈
                world_dist_right_left atol=1e-6
        end

        @testset "reposition! no crash" begin
            sys = sam_c.state
            # Should not crash for chained transforms
            @test_nowarn reposition!(sys)
        end

        @testset "Chained reposition! turns each transform's velocities with it" begin
            sys = sam_c.state
            init!(sam_c; prn=false)
            members = [set_member_velocities!(sys, transform)
                       for transform in sys.transforms]
            pos_before = [[copy(member.pos_w) for member in group] for group in members]
            vel_before = [[copy(member.vel_w) for member in group] for group in members]

            sys.transforms[:main_tf].heading = deg2rad(20)
            sys.transforms[:kite_tilt].elevation += deg2rad(5)
            reposition!(sys)

            for (transform, group, pos, vel) in
                zip(sys.transforms, members, pos_before, vel_before)
                base_pos = base_position(Pose(sys), sys.transforms, transform)
                before = [dot(velocity, position .- base_pos)
                          for velocity in vel, position in pos]
                @test motion_against_pose(group, base_pos) ≈ before atol=1e-8
            end
        end

        @testset "Different child elevation changes pos" begin
            sys = sam_c.state
            tf_child = sys.transforms[:kite_tilt]
            elev = tf_child.elevation - kite_angle

            # Elevation 1
            tf_child.elevation = elev + kite_angle
            place!(sam_c; prn=false)
            init!(sam_c; prn=false)
            top_pos1 = copy(sys.points[:top].pos_w)

            # Elevation 2 (larger tilt)
            tf_child.elevation = elev + 2 * kite_angle
            place!(sam_c; prn=false)
            init!(sam_c; prn=false)
            top_pos2 = copy(sys.points[:top].pos_w)

            # Different elevation should give different pos
            @test !(top_pos1 ≈ top_pos2)
        end
    end

    # Cleanup
    rm(tmpdir; recursive=true)
end
nothing
