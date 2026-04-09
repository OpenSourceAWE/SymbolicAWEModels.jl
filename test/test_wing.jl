# SPDX-FileCopyrightText: 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

# test_wing.jl - Wing sanity tests
#
# Simple physics sanity checks for both REFINE and QUATERNION wings.
# Uses 2plate_kite configuration. All tests use winch brake engaged
# and loose tolerances to catch "crazy stuff" without being brittle.

using Pkg
Pkg.activate(@__DIR__)

using Test
using SymbolicAWEModels
using SymbolicAWEModels: KVec3, VortexStepMethod
using KiteUtils
using LinearAlgebra

"""
    reset_state!(sam, set)

Reset SAM to default state between tests. Restores gravity, wind,
transform, constraints, damping, and brake. Calls `init!` to
reinitialize the integrator.
"""
function reset_state!(sam, set)
    set.g_earth = 9.81
    set.v_wind = 15.0

    tf = sam.sys_struct.transforms[:main_transform]
    tf.elevation = deg2rad(60)
    tf.azimuth = 0.0
    tf.heading = 0.0
    tf.elevation_vel = 0.0
    tf.azimuth_vel = 0.0

    for point in sam.sys_struct.points
        if point.type == SymbolicAWEModels.DYNAMIC
            point.fix_sphere = false
            point.fix_static = false
        end
    end

    for wing in sam.sys_struct.wings
        wing.fix_sphere = false
    end

    set_world_frame_damping(sam.sys_struct, 0.0)

    sam.sys_struct.winches[:main_winch].brake = true

    init!(sam; remake=false, reload=false, prn=false)
end

@testset "Wing Tests" begin
    # Copy 2plate_kite data to temp directory
    pkg_root = dirname(@__DIR__)
    src_data_path = joinpath(
        pkg_root, "data", "2plate_kite"
    )

    tmpdir = mktempdir()
    data_path = joinpath(tmpdir, "2plate_kite")
    cp(src_data_path, data_path; force=true)

    set_data_path(data_path)
    set = Settings("system.yaml")

    vsm_settings_path = joinpath(
        data_path, "vsm_settings.yaml"
    )
    vsm_set = VortexStepMethod.VSMSettings(
        vsm_settings_path; data_prefix=false
    )

    quat_yaml_path = joinpath(
        data_path, "quat_struc_geometry.yaml"
    )
    refine_yaml_path = joinpath(
        data_path, "refine_struc_geometry.yaml"
    )

    # Build SAMs once (expensive symbolic compilation)
    quat_sys = load_sys_struct_from_yaml(
        quat_yaml_path;
        system_name="wing_test_QUATERNION", set, vsm_set
    )
    quat_sam = SymbolicAWEModel(set, quat_sys)
    init!(quat_sam; remake=true, prn=true)

    refine_sys = load_sys_struct_from_yaml(
        refine_yaml_path;
        system_name="wing_test_REFINE", set, vsm_set
    )
    refine_sam = SymbolicAWEModel(set, refine_sys)
    init!(refine_sam; remake=true, prn=true)

    sam_configs = [
        ("REFINE", refine_sam, SymbolicAWEModels.REFINE),
        ("QUATERNION", quat_sam, SymbolicAWEModels.QUATERNION),
    ]

    for (wtn, sam, expected_wing_type) in sam_configs
        @testset "$wtn Wing" begin
            # ========================================================
            # Test 0: YAML Loading Verification
            # ========================================================
            @testset "YAML loading" begin
                sys = sam.sys_struct

                @test length(sys.wings) == 1
                @test haskey(sys.wings, :main_wing)
                wing = sys.wings[:main_wing]
                @test wing.wing_type == expected_wing_type

                @test haskey(sys.points, :le_left)
                @test haskey(sys.points, :te_left)
                @test haskey(sys.points, :le_center)
                @test haskey(sys.points, :te_center)
                @test haskey(sys.points, :le_right)
                @test haskey(sys.points, :te_right)
                @test haskey(sys.points, :kcu)
                @test haskey(sys.points, :steering_left)
                @test haskey(sys.points, :steering_right)

                @test sys.points[:le_left].type ==
                    SymbolicAWEModels.WING

                @test length(sys.transforms) == 1
                @test haskey(sys.transforms, :main_transform)
            end

            # ========================================================
            # Test 1: No forces → barely moves
            # ========================================================
            @testset "No wind no gravity" begin
                reset_state!(sam, set)
                set.g_earth = 0.0
                set.v_wind = 0.0
                init!(
                    sam; remake=false, reload=false, prn=false
                )
                @show sam.sys_struct.wind_vec_gnd
                @show sam.sys_struct.wings[1].heading

                sys = sam.sys_struct
                wing = sys.wings[:main_wing]
                init_wing_pos = copy(wing.pos_w)
                init_kcu_pos = copy(sys.points[:kcu].pos_w)

                for _ in 1:50
                    next_step!(sam; dt=0.05,
                        vsm_interval=0)
                end

                wing_drift = norm(wing.pos_w - init_wing_pos)
                kcu_drift = norm(
                    sys.points[:kcu].pos_w - init_kcu_pos
                )
                wing_speed = norm(wing.vel_w)

                @test wing_drift < 1.0
                @test kcu_drift < 1.0
                @test wing_speed < 2.0

                println("  [$wtn] No forces: " *
                    "wing_drift=$(round(wing_drift; digits=3))" *
                    "m, kcu_drift=" *
                    "$(round(kcu_drift; digits=3))m, " *
                    "speed=$(round(wing_speed; digits=3))m/s")
            end

            # ========================================================
            # Test 2: Gravity → falls
            # ========================================================
            @testset "Gravity no wind" begin
                reset_state!(sam, set)
                set.v_wind = 0.0
                init!(
                    sam; remake=false, reload=false, prn=false
                )

                wing = sam.sys_struct.wings[:main_wing]
                initial_z = wing.pos_w[3]

                for _ in 1:100
                    next_step!(sam; dt=0.05,
                        vsm_interval=0)
                end

                final_z = wing.pos_w[3]

                @test final_z < initial_z

                println("  [$wtn] Gravity: z: " *
                    "$(round(initial_z; digits=2)) → " *
                    "$(round(final_z; digits=2))m")
            end

            # ========================================================
            # Test 3: Wind+gravity → stable (doesn't blow up)
            # ========================================================
            @testset "Wind+gravity stable" begin
                reset_state!(sam, set)

                wing = sam.sys_struct.wings[:main_wing]
                initial_norm = norm(wing.pos_w)

                for _ in 1:100
                    next_step!(sam; dt=0.001,
                        vsm_interval=1)
                end

                final_norm = norm(wing.pos_w)

                @test final_norm < 2 * initial_norm

                println("  [$wtn] Stable: norm: " *
                    "$(round(initial_norm; digits=1)) → " *
                    "$(round(final_norm; digits=1))m")
            end

            # ========================================================
            # Test 4: Wind → kite moves
            # ========================================================
            @testset "Kite moves" begin
                reset_state!(sam, set)

                wing = sam.sys_struct.wings[:main_wing]
                initial_pos = copy(wing.pos_w)

                for _ in 1:100
                    next_step!(sam; dt=0.001,
                        vsm_interval=1)
                end

                displacement = norm(wing.pos_w - initial_pos)
                speed = norm(wing.vel_w)

                @test displacement > 0.05
                @test speed > 0.05

                println("  [$wtn] Moves: " *
                    "disp=$(round(displacement; digits=2))m," *
                    " speed=$(round(speed; digits=2))m/s")
            end

            # ========================================================
            # Test 5: Wind → kite rotates
            # ========================================================
            @testset "Kite rotates" begin
                reset_state!(sam, set)

                wing = sam.sys_struct.wings[:main_wing]
                initial_Q = copy(wing.Q_b_to_w)

                for _ in 1:100
                    next_step!(sam; dt=0.001,
                        vsm_interval=1)
                end

                q_diff = norm(wing.Q_b_to_w - initial_Q)
                q_norm = norm(wing.Q_b_to_w)

                @test q_diff > 0.001
                @test q_norm ≈ 1.0 atol = 0.1

                println("  [$wtn] Rotates: " *
                    "q_diff=$(round(q_diff; digits=4)), " *
                    "q_norm=$(round(q_norm; digits=4))")
            end

            # ========================================================
            # Test 6: fix_sphere → radial motion only
            # ========================================================
            @testset "fix_sphere" begin
                reset_state!(sam, set)
                set.v_wind = 0.0
                for point in sam.sys_struct.points
                    if point.type == SymbolicAWEModels.DYNAMIC
                        point.fix_sphere = true
                    end
                end
                for wing in sam.sys_struct.wings
                    wing.fix_sphere = true
                end
                init!(
                    sam; remake=false, reload=false, prn=false
                )

                # QUATERNION: fix_sphere constrains com_w
                # (not kcu, which shifts with rotation).
                # REFINE: no separate COM, check kcu directly.
                wing = sam.sys_struct.wings[:main_wing]
                if expected_wing_type ==
                        SymbolicAWEModels.QUATERNION
                    p0 = copy(wing.com_w)
                else
                    p0 = copy(
                        sam.sys_struct.points[:kcu].pos_w
                    )
                end
                r̂0 = normalize(p0)

                for _ in 1:100
                    next_step!(sam; dt=0.05, vsm_interval=0)
                end

                if expected_wing_type == SymbolicAWEModels.QUATERNION
                    p1 = copy(wing.com_w)
                else
                    p1 = copy(sam.sys_struct.points[:kcu].pos_w)
                end

                Δp = p1 - p0
                tangential = norm(Δp - dot(Δp, r̂0) * r̂0)
                @test tangential < 9e-3

                println("  [$wtn] fix_sphere: " *
                    "tangential=" *
                    "$(round(tangential; digits=6))")
            end

            # ========================================================
            # Test 7: fix_static → DYNAMIC points frozen
            # ========================================================
            @testset "fix_static" begin
                reset_state!(sam, set)
                set.v_wind = 0.0
                for point in sam.sys_struct.points
                    if point.type == SymbolicAWEModels.DYNAMIC
                        point.fix_static = true
                    end
                end
                init!(
                    sam; remake=false, reload=false, prn=false
                )

                # QUATERNION: kcu is a WING point (derived
                # from com_w + R * pos_b), not a DYNAMIC
                # point, so fix_static has no effect on it.
                check_names = if expected_wing_type ==
                        SymbolicAWEModels.QUATERNION
                    [:steering_left, :steering_right]
                else
                    [:kcu, :steering_left,
                        :steering_right]
                end
                initial_positions = Dict(
                    n => copy(sam.sys_struct.points[n].pos_w)
                    for n in check_names
                )

                for _ in 1:50
                    next_step!(sam; dt=0.05,
                        vsm_interval=0)
                end

                for name in check_names
                    drift = norm(
                        sam.sys_struct.points[name].pos_w -
                        initial_positions[name]
                    )
                    drift_tol = if expected_wing_type ==
                            SymbolicAWEModels.QUATERNION
                        # QUATERNION keeps wing points on the rigid body;
                        # cross-platform solver differences can cause
                        # millimeter-level residual drift in this test.
                        5e-2
                    else
                        1e-5
                    end
                    @test drift < drift_tol
                end

                println("  [$wtn] fix_static: frozen")
            end

            # ========================================================
            # Test 8: High damping → slow motion
            # ========================================================
            @testset "High damping" begin
                # Undamped run
                reset_state!(sam, set)
                set.v_wind = 0.0
                init!(
                    sam; remake=false, reload=false, prn=false
                )

                wing = sam.sys_struct.wings[:main_wing]
                init_pos = copy(wing.pos_w)

                for _ in 1:100
                    next_step!(sam; dt=0.05,
                        vsm_interval=0)
                end

                undamped_speed = norm(wing.vel_w)
                undamped_drift = norm(wing.pos_w - init_pos)

                # Damped run
                reset_state!(sam, set)
                set.v_wind = 0.0
                set_world_frame_damping(
                    sam.sys_struct, 1000.0
                )
                init!(
                    sam; remake=false, reload=false, prn=false
                )

                init_pos .= wing.pos_w

                for _ in 1:100
                    next_step!(sam; dt=0.05,
                        vsm_interval=0)
                end

                damped_speed = norm(wing.vel_w)
                damped_drift = norm(wing.pos_w - init_pos)

                @test damped_speed < 0.5 * undamped_speed
                @test damped_drift < 0.5 * undamped_drift

                println("  [$wtn] Damping: " *
                    "speed $(round(undamped_speed; digits=2))" *
                    " → $(round(damped_speed; digits=2))m/s," *
                    " drift $(round(undamped_drift; digits=2))" *
                    " → $(round(damped_drift; digits=2))m")
            end
        end
    end

    # Cleanup
    rm(tmpdir; recursive=true)
end
nothing
