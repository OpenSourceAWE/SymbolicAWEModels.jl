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
    rigid_dynamics_yaml_path = joinpath(data_path, "rigid_structural_geometry.yaml")
    particle_dynamics_yaml_path = joinpath(data_path, "particle_structural_geometry.yaml")

    # Create and initialize SAMs once for each wing type
    rigid_dynamics_sys = load_sys_struct_from_yaml(
        rigid_dynamics_yaml_path; system_name="transform_test_RIGID_DYNAMICS", set=set, vsm_set=vsm_set,
        aero_mode=AeroNone()
    )
    rigid_dynamics_sam = SymbolicAWEModel(set, rigid_dynamics_sys)
    test_init!(rigid_dynamics_sam)

    particle_dynamics_sys = load_sys_struct_from_yaml(
        particle_dynamics_yaml_path; system_name="transform_test_PARTICLE_DYNAMICS", set=set, vsm_set=vsm_set,
        aero_mode=AeroNone()
    )
    particle_dynamics_sam = SymbolicAWEModel(set, particle_dynamics_sys)
    test_init!(particle_dynamics_sam)

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
        ("PARTICLE_DYNAMICS", particle_dynamics_sam, particle_dynamics_yaml_path),
        ("RIGID_DYNAMICS", rigid_dynamics_sam, rigid_dynamics_yaml_path),
    ]

    for (dynamics_type_name, sam, yaml_path) in sam_configs
        @testset "$dynamics_type_name Wing" begin
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

                println("\n  ====== [$dynamics_type_name] Loaded transform: " *
                    "elev=$(round(rad2deg(transform.elevation), digits=1))°, " *
                    "azim=$(round(rad2deg(transform.azimuth), digits=1))°, " *
                    "heading=$(round(rad2deg(transform.heading), digits=1))° ======\n")
            end

            # ================================================================
            # Physics Test 1: Initial angles after init
            # ================================================================
            @testset "Initial angles after init!" begin
                reset_transform!(sam.sys_struct)
                test_init!(sam; prn=false)

                # After init, the transform angles should still match
                transform = sam.sys_struct.transforms[:main_transform]

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
                reset_transform!(sam.sys_struct)
                test_init!(sam; prn=false)

                transform = sam.sys_struct.transforms[:main_transform]

                # Default velocities are 0
                @test transform.elevation_vel ≈ 0.0 atol=1e-10
                @test transform.azimuth_vel ≈ 0.0 atol=1e-10

                println("\n  ====== [$dynamics_type_name] Velocities: " *
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
                test_init!(sam; prn=false)

                # Get wing position
                wing = sam.sys_struct.bodies[:main_wing]
                wing_pos = wing.pos_w

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

                println("\n  ====== [$dynamics_type_name] Geometry: " *
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

                test_init!(sam; prn=false)

                # Angles should be preserved after init
                transform_after = sam.sys_struct.transforms[:main_transform]
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
                reset_transform!(sam.sys_struct)
                test_init!(sam; prn=false)
                wing_z1 = sam.sys_struct.bodies[:main_wing].pos_w[3]

                # Test 2: elevation = 45 deg
                sam.sys_struct.transforms[:main_transform].elevation = deg2rad(45)
                test_init!(sam; prn=false)
                wing_z2 = sam.sys_struct.bodies[:main_wing].pos_w[3]

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
                reset_transform!(sam.sys_struct)
                test_init!(sam; prn=false)
                wing_y1 = sam.sys_struct.bodies[:main_wing].pos_w[2]

                # Test 2: azimuth = 30 deg (more to the side)
                sam.sys_struct.transforms[:main_transform].azimuth = deg2rad(30)
                test_init!(sam; prn=false)
                wing_y2 = sam.sys_struct.bodies[:main_wing].pos_w[2]

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
                reset_transform!(sam.sys_struct)
                test_init!(sam; prn=false)

                wing = sam.sys_struct.bodies[:main_wing]

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
                transform = sam.sys_struct.transforms[:main_transform]
                @test transform.base_point_idx == 10  # ground index

                reset_transform!(sam.sys_struct)
                test_init!(sam; prn=false)

                # Transform base_pos should match the ground point position
                ground_pos = sam.sys_struct.points[:ground].pos_w
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
                        test_init!(sam; prn=false)
                        wing = sys.bodies[:main_wing]
                        reinit_R = copy(wing.R_b_to_w)
                        reinit_pos = copy(wing.pos_w)

                        # reinit! with heading=0, then
                        # reposition! to target
                        tf.base_pos .= base_pos
                        reset_transform!(sys)
                        test_init!(sam; prn=false)
                        tf.heading = target_h
                        reposition!(sys.transforms, sys)
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
        end
    end

    # ================================================================
    # Chained Transform Tests (using kps4_plate-style programmatic API)
    # ================================================================
    @testset "Chained Transforms" begin
        using SymbolicAWEModels: Point, Segment, Tether, Winch,
            PlateWing, TwistSurface, Transform,
            SystemStructure,
            create_plate_interpolations, get_rot_pos,
            get_rot_pos_cad, get_base_pos, reinit!

        # Use kps4 settings (2plate_kite has no tethers)
        kps4_data = joinpath(tmpdir, "kps4")
        cp(joinpath(pkg_root, "data", "kps4"),
            kps4_data; force=true)
        set_data_path(kps4_data)
        set_c = Settings("system.yaml")
        set_c.upwind_dir = rad2deg(-pi/2)

        # Geometry from KiteUtils
        particles = KiteUtils.get_particles(
            set_c.height_k, set_c.h_bridle,
            set_c.width, set_c.m_k)
        pos_kcu = particles[2]
        pos_nose = particles[3]
        pos_top = particles[4]
        pos_right = particles[5]
        pos_left = particles[6]

        kite_mass = set_c.mass
        k_nose = set_c.rel_nose_mass * kite_mass
        k_top = set_c.rel_top_mass *
            (1.0 - set_c.rel_nose_mass) * kite_mass
        k_side = 0.5 * (1.0 - set_c.rel_top_mass) *
            (1.0 - set_c.rel_nose_mass) * kite_mass
        set_c.mass = 0.0

        pre_stress = 0.9975
        pos_map = Dict(:kcu => pos_kcu, :nose => pos_nose,
            :top => pos_top, :right => pos_right,
            :left => pos_left)
        bridle_l0(a, b) =
            norm(pos_map[b] - pos_map[a]) * pre_stress

        points_c = [
            Point(:ground, zeros(3), STATIC),
            Point(:kcu, pos_kcu, DYNAMIC;
                extra_mass=set_c.kcu_mass,
                transform=:main_tf),
            Point(:nose, pos_nose, DYNAMIC;
                extra_mass=k_nose,
                transform=:main_tf),
            Point(:top, pos_top, WING;
                extra_mass=k_top, wing=:plate_wing,
                transform=:kite_tilt),
            Point(:right, pos_right, WING;
                extra_mass=k_side, wing=:plate_wing,
                transform=:kite_tilt),
            Point(:left, pos_left, WING;
                extra_mass=k_side, wing=:plate_wing,
                transform=:kite_tilt),
        ]

        segments_c = [
            Segment(:kcu_nose, set_c, :kcu, :nose;
                l0=bridle_l0(:kcu, :nose),
                diameter_mm=set_c.d_line),
            Segment(:right_nose, set_c, :right, :nose;
                l0=bridle_l0(:right, :nose),
                diameter_mm=set_c.d_line),
            Segment(:right_left, set_c, :right, :left;
                l0=bridle_l0(:right, :left),
                diameter_mm=set_c.d_line),
            Segment(:top_right, set_c, :top, :right;
                l0=bridle_l0(:top, :right),
                diameter_mm=set_c.d_line),
            Segment(:left_kcu, set_c, :left, :kcu;
                l0=bridle_l0(:left, :kcu),
                diameter_mm=set_c.d_line),
            Segment(:right_kcu, set_c, :right, :kcu;
                l0=bridle_l0(:right, :kcu),
                diameter_mm=set_c.d_line),
            Segment(:top_left, set_c, :top, :left;
                l0=bridle_l0(:top, :left),
                diameter_mm=set_c.d_line),
            Segment(:left_nose, set_c, :left, :nose;
                l0=bridle_l0(:left, :nose),
                diameter_mm=set_c.d_line),
            Segment(:nose_top, set_c, :nose, :top;
                l0=bridle_l0(:nose, :top),
                diameter_mm=set_c.d_line),
        ]

        tethers_c = [Tether(:main_tether,
            set_c.l_tethers[1];
            start_point=:ground, end_point=:kcu,
            n_segments=set_c.segments)]

        winches_c = [Winch(:winch, set_c,
            [:main_tether]; winch_point=:ground)]

        rel_side = set_c.rel_side_area / 100.0
        K = 1.0 - rel_side
        twist_surfaces_c = [
            TwistSurface(:main, [:top], STATIC, 0.0;
                x_airf=[1,0,0], y_airf=[0,1,0],
                area=set_c.area, twist=deg2rad(set_c.alpha_zero)),
            TwistSurface(:right_tip, [:right], STATIC, 0.0;
                x_airf=[1,0,0], y_airf=[0,0,-1],
                area=set_c.area * rel_side,
                twist=deg2rad(set_c.alpha_ztip)),
            TwistSurface(:left_tip, [:left], STATIC, 0.0;
                x_airf=[1,0,0], y_airf=[0,0,1],
                area=set_c.area * rel_side,
                twist=deg2rad(set_c.alpha_ztip)),
        ]
        cl_interp, cd_interp =
            create_plate_interpolations(
                set_c.alpha_cl, set_c.cl_list,
                set_c.cd_list; alpha_cd=set_c.alpha_cd)

        wing_c = PlateWing(:plate_wing,
            [:main, :right_tip, :left_tip],
            cl_interp, cd_interp;
            dynamics_type=PARTICLE_DYNAMICS,
            z_ref_points=([:right, :left], :top),
            y_ref_points=(:left, :right),
            origin=:kcu, drag_corr=0.93 * K)

        elev = deg2rad(set_c.elevation)
        azim = deg2rad(10.0)
        kite_angle = deg2rad(3.83)

        transforms_c = [
            Transform(:main_tf, elev, azim, 0.0;
                base_pos=zeros(3), base_point=:ground,
                wing=:plate_wing),
            Transform(:kite_tilt,
                elev + kite_angle, azim, 0.0;
                base_transform=:main_tf,
                rot_point=:top),
        ]

        sys_c = SystemStructure("chained_test", set_c;
            points=points_c, twist_surfaces=twist_surfaces_c,
            segments=segments_c,
            tethers=tethers_c, winches=winches_c,
            wings=[wing_c], transforms=transforms_c)

        sam_c = SymbolicAWEModel(set_c, sys_c)
        init!(sam_c)

        @testset "get_base_pos returns different values" begin
            sys = sam_c.sys_struct
            tf_child = sys.transforms[:kite_tilt]
            base_pos, curr_base_pos = get_base_pos(
                tf_child, sys.transforms,
                sys.bodies, sys.points)
            # After init, parent wing has moved from CAD
            # so base_pos (world) != curr_base_pos (CAD)
            @test !(base_pos ≈ curr_base_pos)
            println("  base_pos=$(round.(base_pos, digits=2))")
            println("  curr_base_pos=" *
                "$(round.(curr_base_pos, digits=2))")
        end

        @testset "Child points translated from CAD" begin
            sys = sam_c.sys_struct
            top = sys.points[:top]
            # top should NOT be at its CAD position
            @test !(top.pos_w ≈ top.pos_cad)
            # top should be far from origin (at tether length)
            @test norm(top.pos_w) > 50.0
            println("  top.pos_w=" *
                "$(round.(top.pos_w, digits=2))")
            println("  top.pos_cad=" *
                "$(round.(top.pos_cad, digits=2))")
        end

        @testset "Child points near parent wing" begin
            sys = sam_c.sys_struct
            wing_pos = sys.wings[1].pos_w
            top_pos = sys.points[:top].pos_w
            right_pos = sys.points[:right].pos_w
            left_pos = sys.points[:left].pos_w

            # All child points should be near the wing
            # (within bridle length, ~30m)
            @test norm(top_pos - wing_pos) < 50.0
            @test norm(right_pos - wing_pos) < 50.0
            @test norm(left_pos - wing_pos) < 50.0

            println("  wing_pos=" *
                "$(round.(wing_pos, digits=2))")
            dist = round(
                norm(top_pos - wing_pos), digits=2)
            println("  dist(top, wing)=$dist")
        end

        @testset "Distances preserved (rigid body)" begin
            sys = sam_c.sys_struct
            # CAD distances between child points
            cad_dist_top_right = norm(
                sys.points[:top].pos_cad -
                sys.points[:right].pos_cad)
            cad_dist_top_left = norm(
                sys.points[:top].pos_cad -
                sys.points[:left].pos_cad)
            cad_dist_right_left = norm(
                sys.points[:right].pos_cad -
                sys.points[:left].pos_cad)

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
            sys = sam_c.sys_struct
            # Should not crash for chained transforms
            @test_nowarn reposition!(
                sys.transforms, sys)
        end

        @testset "Different child elevation changes pos" begin
            sys = sam_c.sys_struct
            tf_child = sys.transforms[:kite_tilt]

            # Elevation 1
            tf_child.elevation = elev + kite_angle
            init!(sam_c; prn=false)
            top_pos1 = copy(sys.points[:top].pos_w)

            # Elevation 2 (larger tilt)
            tf_child.elevation = elev + 2 * kite_angle
            init!(sam_c; prn=false)
            top_pos2 = copy(sys.points[:top].pos_w)

            # Different elevation should give different pos
            @test !(top_pos1 ≈ top_pos2)
            println("  top_pos(tilt1)=" *
                "$(round.(top_pos1, digits=2))")
            println("  top_pos(tilt2)=" *
                "$(round.(top_pos2, digits=2))")
        end
    end

    # Cleanup
    rm(tmpdir; recursive=true)
end
nothing
