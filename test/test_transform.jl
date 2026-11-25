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
using SymbolicAWEModels: KVec3, VortexStepMethod, origin_ENU
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
    rigid_dynamics_sys, rigid_dynamics_transforms = load_authoring(
        rigid_dynamics_yaml_path; name="transform_test_RIGID_DYNAMICS", set=set
    )
    rigid_dynamics_sam = SymbolicAWEModel(set, rigid_dynamics_sys;
        transforms=rigid_dynamics_transforms, vsm_set,
        aero=Dict(wing.name => AeroNone() for wing in rigid_dynamics_sys.wings))
    test_init!(rigid_dynamics_sam)

    particle_dynamics_sys, particle_dynamics_transforms = load_authoring(
        particle_dynamics_yaml_path;
        name="transform_test_PARTICLE_DYNAMICS", set=set
    )
    particle_dynamics_sam = SymbolicAWEModel(set, particle_dynamics_sys;
        transforms=particle_dynamics_transforms, vsm_set,
        aero=Dict(wing.name => AeroNone() for wing in particle_dynamics_sys.wings))
    test_init!(particle_dynamics_sam)

    # Helper to reset transform to default YAML values
    function reset_transform!(sam)
        tf = sam.transforms[:main_transform]
        tf.elevation = deg2rad(80)
        tf.azimuth = deg2rad(0)
        tf.heading = deg2rad(0)
        tf.elevation_vel = 0.0
        tf.azimuth_vel = 0.0
        tf.turn_rate = 0.0
    end

    # Test both wing types
    sam_configs = [
        ("PARTICLE_DYNAMICS", particle_dynamics_sam, particle_dynamics_yaml_path),
        ("RIGID_DYNAMICS", rigid_dynamics_sam, rigid_dynamics_yaml_path),
    ]

    for (dynamics_type_name, sam, yaml_path) in sam_configs
        @testset "$dynamics_type_name Wing" begin
            # ================================================================
            # YAML Loading Verification (uses the already-built model)
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
                reset_transform!(sam)
                place!(sam)
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
                reset_transform!(sam)
                place!(sam)
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
            @testset "azimuth rates have the sign of a finite difference of azimuth" begin
                azimuth(t) = 0.3 + 0.4t + 0.2t^2
                position(t) = 100.0 .* [cos(azimuth(t)), sin(azimuth(t)), 0.5]
                t, dt = 1.0, 1e-5
                velocity = (position(t + dt) - position(t - dt)) / 2dt
                acceleration = (position(t + dt) - 2position(t) + position(t - dt)) / dt^2
                rates = SymbolicAWEModels.wing_scalar_kinematics(; rel_pos=position(t),
                    e_x=[1.0, 0.0, 0.0],
                    R_t_to_w=SymbolicAWEModels.sym_calc_R_t_to_w(position(t)),
                    R_v_to_w=[1.0 0 0; 0 1.0 0; 0 0 1.0],
                    R_b_to_w=[1.0 0 0; 0 1.0 0; 0 0 1.0],
                    vel=velocity, acc=acceleration, omega_b=zeros(3), alpha_b=zeros(3),
                    va_b=[1.0, 0.0, 0.0], twist_offset=0.0)
                @test rates.azimuth ≈ azimuth(t) atol=1e-8
                @test rates.azimuth_vel ≈ 0.4 + 0.4t atol=1e-6
                @test sign(rates.azimuth_acc) == sign(0.4)
            end

            @testset "Start velocity round-trips into the wing rates" begin
                tf = sam.transforms[:main_transform]
                reset_transform!(sam)
                tf.elevation = deg2rad(70)
                tf.azimuth = deg2rad(25)
                tf.elevation_vel = deg2rad(3.0)
                tf.azimuth_vel = deg2rad(-2.0)

                place!(sam)

                wing(quantity) = run_value(sam, :wings, :main_wing, quantity)
                rel_pos = wing(:pos_w) - tf.base_pos
                @test atan(rel_pos[2], rel_pos[1]) ≈ tf.azimuth atol=1e-8
                @test wing(:elevation_vel) ≈ tf.elevation_vel atol=1e-8
                @test wing(:azimuth_vel) ≈ tf.azimuth_vel atol=1e-8

                println("\n  ====== [$dynamics_type_name] Start velocity: " *
                    "elev_vel=$(round(rad2deg(wing(:elevation_vel)), digits=3))°/s, " *
                    "azim_vel=$(round(rad2deg(wing(:azimuth_vel)), digits=3))°/s ======\n")

                reset_transform!(sam)
            end

            # ================================================================
            # Physics Test 2c: turn_rate spins the structure as it is placed
            # ================================================================
            @testset "turn_rate advances the placed heading" begin
                tf = sam.transforms[:main_transform]
                turn_rate = deg2rad(6.0)
                dt = 1e-5
                points = eachindex(sam.definition.points)
                bodies = eachindex(sam.definition.bodies)

                reset_transform!(sam)
                tf.elevation = deg2rad(70)
                tf.azimuth = deg2rad(25)
                tf.heading = deg2rad(15)
                tf.turn_rate = turn_rate
                place!(sam)
                point_pos = [run_value(sam, :points, i, :pos_w) for i in points]
                point_vel = [run_value(sam, :points, i, :vel_w) for i in points]
                body_pos = [run_value(sam, :bodies, i, :pos_w) for i in bodies]
                body_vel = [run_value(sam, :bodies, i, :vel_w) for i in bodies]

                tf.heading += turn_rate * dt
                tf.turn_rate = 0.0
                place!(sam)

                for i in points
                    @test (run_value(sam, :points, i, :pos_w) .- point_pos[i]) ./ dt ≈
                          point_vel[i] atol=1e-6
                end
                for i in bodies
                    @test (run_value(sam, :bodies, i, :pos_w) .- body_pos[i]) ./ dt ≈
                          body_vel[i] atol=1e-6
                end

                println("\n  ====== [$dynamics_type_name] Turn rate: " *
                    "$(round(rad2deg(turn_rate), digits=2))°/s moves the fastest " *
                    "point at $(round(maximum(norm, point_vel), digits=3)) m/s ======\n")

                reset_transform!(sam)
            end

            # ================================================================
            # Physics Test 3: Geometric consistency - position from spherical coords
            # ================================================================
            @testset "Geometric consistency" begin
                # For elevation=80deg, azimuth=0deg, the wing should be positioned
                # according to spherical coordinate transformation
                reset_transform!(sam)
                place!(sam)
                test_init!(sam; prn=false)

                # Get wing position
                wing_pos = run_value(sam, :bodies, :main_wing, :pos_w)

                # Get ground position
                ground_pos = run_value(sam, :points, :ground, :pos_w)

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

                place!(sam)
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
                reset_transform!(sam)
                place!(sam)
                test_init!(sam; prn=false)
                wing_z1 = run_value(sam, :bodies, :main_wing, :pos_w)[3]

                # Test 2: elevation = 45 deg; init! alone keeps the placement
                sam.transforms[:main_transform].elevation = deg2rad(45)
                test_init!(sam; prn=false)
                @test run_value(sam, :bodies, :main_wing, :pos_w)[3] ≈ wing_z1 atol=1e-6
                place!(sam)
                test_init!(sam; prn=false)
                wing_z2 = run_value(sam, :bodies, :main_wing, :pos_w)[3]

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
                reset_transform!(sam)
                place!(sam)
                test_init!(sam; prn=false)
                wing_y1 = run_value(sam, :bodies, :main_wing, :pos_w)[2]

                # Test 2: azimuth = 30 deg (more to the side)
                sam.transforms[:main_transform].azimuth = deg2rad(30)
                place!(sam)
                test_init!(sam; prn=false)
                wing_y2 = run_value(sam, :bodies, :main_wing, :pos_w)[2]

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
                reset_transform!(sam)
                place!(sam)
                test_init!(sam; prn=false)

                R_b_to_w = run_value(sam, :bodies, :main_wing, :R_b_to_w)

                # Wing should have a rotation matrix
                @test !isnothing(R_b_to_w)
                @test size(R_b_to_w) == (3, 3)

                # R_b_to_w should be a valid rotation matrix (orthonormal)
                @test det(R_b_to_w) ≈ 1.0 atol=1e-10
                @test R_b_to_w * R_b_to_w' ≈ I(3) atol=1e-10

                println("\n  ====== [$dynamics_type_name] Heading affects rotation: " *
                    "det(R_b_to_w)=$(round(det(R_b_to_w), digits=4)) ======\n")
            end

            # ================================================================
            # Physics Test 8: Base position offset
            # ================================================================
            @testset "Base position from base_point" begin
                # The transform references ground as base_point
                transform = sam.transforms[:main_transform]
                @test transform.base_point == 10  # ground index

                reset_transform!(sam)
                place!(sam)
                test_init!(sam; prn=false)

                # Transform base_pos should match the ground point position
                ground_pos = run_value(sam, :points, :ground, :pos_w)
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
                tf = sam.transforms[:main_transform]
                orig_base_pos = copy(tf.base_pos)
                wing(quantity) = run_value(sam, :bodies, :main_wing, quantity)

                @testset "Reposition heading ($label)" begin
                    for h_deg in [0, 10, -15, 30, 45, -45]
                        target_h = deg2rad(h_deg)

                        # Reference: reinit! with target heading
                        tf.base_pos .= base_pos
                        reset_transform!(sam)
                        tf.heading = target_h
                        place!(sam)
                        test_init!(sam; prn=false)
                        reinit_R = wing(:R_b_to_w)
                        reinit_pos = wing(:pos_w)

                        # reinit! with heading=0, then
                        # reposition! to target
                        tf.base_pos .= base_pos
                        reset_transform!(sam)
                        place!(sam)
                        test_init!(sam; prn=false)
                        tf.heading = target_h
                        reposition!(sam)

                        @test reinit_R ≈ wing(:R_b_to_w) atol=1e-6
                        @test reinit_pos ≈ wing(:pos_w) atol=1e-4
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
        using SymbolicAWEModels: SystemDefinition

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

        points_c = KiteGeometry.Point[
            Point(; name=:ground, pos_ENU=zeros(3), fixed=true),
            Point(; name=:kcu, pos_ENU=pos_kcu,
                extra_mass=set_c.kcu_mass),
            Point(; name=:nose, pos_ENU=pos_nose,
                extra_mass=k_nose),
            Point(; name=:top, pos_ENU=pos_top,
                extra_mass=k_top),
            Point(; name=:right, pos_ENU=pos_right,
                extra_mass=k_side),
            Point(; name=:left, pos_ENU=pos_left,
                extra_mass=k_side),
        ]

        materials_c = KiteGeometry.MaterialTable(KiteGeometry.Material[],
                                                 Dict{String, Any}())
        segments_c = KiteGeometry.Segment[
            KiteGeometry.line_segment(:kcu_nose, (:kcu, :nose), set_c, materials_c;
                rest_length=bridle_l0(:kcu, :nose),
                diameter=0.001 * set_c.d_line, model=SegmentSpring()),
            KiteGeometry.line_segment(:right_nose, (:right, :nose), set_c, materials_c;
                rest_length=bridle_l0(:right, :nose),
                diameter=0.001 * set_c.d_line, model=SegmentSpring()),
            KiteGeometry.line_segment(:right_left, (:right, :left), set_c, materials_c;
                rest_length=bridle_l0(:right, :left),
                diameter=0.001 * set_c.d_line, model=SegmentSpring()),
            KiteGeometry.line_segment(:top_right, (:top, :right), set_c, materials_c;
                rest_length=bridle_l0(:top, :right),
                diameter=0.001 * set_c.d_line, model=SegmentSpring()),
            KiteGeometry.line_segment(:left_kcu, (:left, :kcu), set_c, materials_c;
                rest_length=bridle_l0(:left, :kcu),
                diameter=0.001 * set_c.d_line, model=SegmentSpring()),
            KiteGeometry.line_segment(:right_kcu, (:right, :kcu), set_c, materials_c;
                rest_length=bridle_l0(:right, :kcu),
                diameter=0.001 * set_c.d_line, model=SegmentSpring()),
            KiteGeometry.line_segment(:top_left, (:top, :left), set_c, materials_c;
                rest_length=bridle_l0(:top, :left),
                diameter=0.001 * set_c.d_line, model=SegmentSpring()),
            KiteGeometry.line_segment(:left_nose, (:left, :nose), set_c, materials_c;
                rest_length=bridle_l0(:left, :nose),
                diameter=0.001 * set_c.d_line, model=SegmentSpring()),
            KiteGeometry.line_segment(:nose_top, (:nose, :top), set_c, materials_c;
                rest_length=bridle_l0(:nose, :top),
                diameter=0.001 * set_c.d_line, model=SegmentSpring()),
        ]

        tethers_c = [KiteGeometry.expand_tether!(points_c, segments_c,
            fill(nothing, length(points_c)), :main_tether, :ground, :kcu,
            set_c.segments, set_c, materials_c;
            init=(; init_stretched_length=set_c.l_tethers[1]),
            model=KiteGeometry.NoModel(), segment_model=SegmentSpring())]

        winches_c = [Winch(; name=:winch,
            tethers=[:main_tether], winch_point=:ground)]

        rel_side = set_c.rel_side_area / 100.0
        K = 1.0 - rel_side
        stations_c = [
            SpanwiseStation(; name=:main, wing=:plate_wing, points=[:top],
                model=StationModel(; type="STATIC", moment_frac=0.0,
                    chord=[1.0, 0.0, 0.0], y_airf=[0.0, 1.0, 0.0],
                    area=set_c.area, twist=deg2rad(set_c.alpha_zero))),
            SpanwiseStation(; name=:right_tip, wing=:plate_wing, points=[:right],
                model=StationModel(; type="STATIC", moment_frac=0.0,
                    chord=[1.0, 0.0, 0.0], y_airf=[0.0, 0.0, -1.0],
                    area=set_c.area * rel_side,
                    twist=deg2rad(set_c.alpha_ztip))),
            SpanwiseStation(; name=:left_tip, wing=:plate_wing, points=[:left],
                model=StationModel(; type="STATIC", moment_frac=0.0,
                    chord=[1.0, 0.0, 0.0], y_airf=[0.0, 0.0, 1.0],
                    area=set_c.area * rel_side,
                    twist=deg2rad(set_c.alpha_ztip))),
        ]
        wing_c = KiteGeometry.Wing(; name=:plate_wing, pos_ENU=zeros(3), body=:plate_wing,
            model=WingModel(; aero_mode="AeroPlate", dynamics_type="PARTICLE_DYNAMICS",
                stations=[:main, :right_tip, :left_tip],
                z_ref_points=[[:right, :left], :top], y_ref_points=[:left, :right],
                origin_idx=:kcu, drag_corr=0.93 * K))
        body_c = Body(; name=:plate_wing, pos_ENU=zeros(3), extra_mass=0.0,
                      extra_inertia_KA=Matrix(Diagonal(ones(3))),
                      model=BodyModel(; angular_damping=0.0))

        elev = deg2rad(set_c.elevation)
        azim = deg2rad(10.0)
        kite_angle = deg2rad(3.83)

        sys_c = system_definition("chained_test", set_c;
            points=points_c, spanwise_stations=stations_c,
            segments=segments_c, materials=materials_c,
            tethers=tethers_c, winches=winches_c,
            wings=[wing_c], bodies=[body_c])

        point_idx(name) = KiteGeometry.row_index(sys_c.points, name)
        tether_idxs = [i for (i, point) in enumerate(sys_c.points)
                       if startswith(String(point.name), "main_tether_point")]
        transforms_c = [
            Transform(; name=:main_tf, elevation=elev, azimuth=azim, heading=0.0,
                base_pos=zeros(3), base_point=point_idx(:ground), body=1,
                points=[point_idx(:kcu), point_idx(:nose), tether_idxs...],
                bodies=[1]),
            Transform(; name=:kite_tilt, elevation=elev + kite_angle, azimuth=azim,
                heading=0.0, base_transform=1, rot_point=[point_idx(:top) => 1.0],
                points=point_idx.([:top, :right, :left])),
        ]
        place!(sys_c, transforms_c)

        sam_c = SymbolicAWEModel(set_c, sys_c; transforms=transforms_c, vsm_set)
        init!(sam_c)

        @testset "get_base_pos returns the parent's current and initial position" begin
            sys = sam_c.definition
            tf_child = sam_c.transforms[:kite_tilt]
            design = Pose([KVec3(point.pos_ENU) for point in sys.points],
                          [KVec3(origin_ENU(sam_c, b)) for b in eachindex(sys.bodies)],
                          [SymbolicAWEModels.quaternion_to_rotation_matrix(body.Q_KA2ENU)
                           for body in sys.bodies])
            base_pos = base_position(Pose(sam_c), sam_c.transforms, tf_child)
            curr_base_pos = base_position(design, sam_c.transforms, tf_child)
            wing_idx = KiteGeometry.row_index(sys.bodies, :plate_wing)
            @test base_pos == run_value(sam_c, :bodies, wing_idx, :pos_w)
            @test curr_base_pos == origin_ENU(sam_c, wing_idx)
        end

        @testset "Child points translated from the authoring geometry" begin
            top_pos = run_value(sam_c, :points, :top, :pos_w)
            @test !(top_pos ≈ pos_top)
            # top should be far from origin (at tether length)
            @test norm(top_pos) > 50.0
            @test sam_c.definition.points[:top].pos_ENU == top_pos
        end

        @testset "Child points near parent wing" begin
            wing_pos = run_value(sam_c, :bodies, 1, :pos_w)
            top_pos = run_value(sam_c, :points, :top, :pos_w)
            right_pos = run_value(sam_c, :points, :right, :pos_w)
            left_pos = run_value(sam_c, :points, :left, :pos_w)

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
            pos_w(name) = run_value(sam_c, :points, name, :pos_w)
            # Authoring distances between child points
            cad_dist_top_right = norm(pos_top - pos_right)
            cad_dist_top_left = norm(pos_top - pos_left)
            cad_dist_right_left = norm(pos_right - pos_left)

            # World distances should match CAD distances
            # (transforms are rigid body)
            world_dist_top_right = norm(pos_w(:top) - pos_w(:right))
            world_dist_top_left = norm(pos_w(:top) - pos_w(:left))
            world_dist_right_left = norm(pos_w(:right) - pos_w(:left))

            @test cad_dist_top_right ≈
                world_dist_top_right atol=1e-6
            @test cad_dist_top_left ≈
                world_dist_top_left atol=1e-6
            @test cad_dist_right_left ≈
                world_dist_right_left atol=1e-6
        end

        @testset "reposition! no crash" begin
            # Should not crash for chained transforms
            @test_nowarn reposition!(sam_c)
        end

        @testset "Different child elevation changes pos" begin
            tf_child = sam_c.transforms[:kite_tilt]

            # Elevation 1
            tf_child.elevation = elev + kite_angle
            place!(sam_c)
            init!(sam_c; prn=false)
            top_pos1 = run_value(sam_c, :points, :top, :pos_w)

            # Elevation 2 (larger tilt)
            tf_child.elevation = elev + 2 * kite_angle
            place!(sam_c)
            init!(sam_c; prn=false)
            top_pos2 = run_value(sam_c, :points, :top, :pos_w)

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

@testset "a wing's twist offset is the twist at the middle of its span" begin
    twist(i) = 0.1i
    @test SymbolicAWEModels.mid_span_twist(twist, Int[]) == 0
    @test SymbolicAWEModels.mid_span_twist(twist, [4]) ≈ 0.4
    @test SymbolicAWEModels.mid_span_twist(twist, 3:5) ≈ 0.4
    @test SymbolicAWEModels.mid_span_twist(twist, 3:6) ≈ 0.45
end
