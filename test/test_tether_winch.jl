# SPDX-FileCopyrightText: 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_tether_winch.jl - Winch dynamics tests
#
# Tests winch motor dynamics in isolation using zero-stiffness
# and low-stiffness tethers. Each test verifies one physical
# behavior with tight tolerances and swept parameters.

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

@isdefined(test_init!) || include(joinpath(@__DIR__, "util.jl"))

using Test
using SymbolicAWEModels
using SymbolicAWEModels: KVec3
using KiteUtils
using LinearAlgebra
using ModelingToolkit
using ModelingToolkit: t_nounits as tt, D_nounits as DD
using SymbolicAWEModels: AbstractWinchModel

# Custom winch model for Test 9: `set_value` is a target tether length; a PD law
# on `len` tracks it. Drum params come from the flat `params` view (mirrors
# TorqueWinch); `friction_epsilon` is the model's own field. Defined at
# top level because structs cannot live inside `@testset` (a local scope).
if !@isdefined(LengthTrackWinch)
    mutable struct LengthTrackWinch <: AbstractWinchModel
        K_p::Float64
        K_d::Float64
        friction_epsilon::Float64
    end
    LengthTrackWinch(K_p, K_d; friction_epsilon=6.0) =
        LengthTrackWinch(K_p, K_d, friction_epsilon)

    function SymbolicAWEModels.winch_component(model::LengthTrackWinch,
            sam, winch_idx; name, params)
        @variables begin
            vel(tt); len(tt); force(tt); set_value(tt)
            brake(tt); acc(tt); friction(tt)
            ω_motor(tt); err(tt); tau_cmd(tt); tau_net(tt)
        end
        drum_radius         = params.winches[winch_idx].drum_radius
        gear_ratio          = params.winches[winch_idx].gear_ratio
        coulomb_friction    = params.winches[winch_idx].coulomb_friction
        viscous_coefficient = params.winches[winch_idx].viscous_coefficient
        inertia_total       = params.winches[winch_idx].inertia_total
        friction_eps        = params.winches[winch_idx].model.friction_epsilon
        smooth_sign(x, e) = x / sqrt(x * x + e * e)
        ratio = drum_radius / gear_ratio
        eqs = [
            ω_motor  ~ vel / ratio
            friction ~ smooth_sign(ω_motor, friction_eps) *
                       coulomb_friction * ratio +
                       viscous_coefficient * ω_motor * ratio^2
            err      ~ set_value - len
            tau_cmd  ~ model.K_p * err - model.K_d * vel
            tau_net  ~ tau_cmd + ratio * force - friction
            acc      ~ ifelse(brake > 0.5, 0.0,
                              ratio * tau_net / inertia_total)
        ]
        return System(eqs, tt,
            [vel, len, force, set_value, brake, acc, friction,
             ω_motor, err, tau_cmd, tau_net],
            [drum_radius, gear_ratio, coulomb_friction, viscous_coefficient,
             inertia_total, friction_eps]; name)
    end
end

# ================================================================
# Minimal YAML: weight at z=-50 connected to ground via tether
# ================================================================
const WINCH_TEST_YAML = """
points:
  headers: [name, pos_cad, type, wing_idx, transform_idx,
            extra_mass, body_frame_damping, world_frame_damping,
            area, drag_coeff]
  data:
    - [weight, [0.0, 0.0, -50.0], DYNAMIC, nothing, nothing,
       10.0, 0.0, 0.0, 0.0, 0.0]
    - [ground, [0.0, 0.0, 0.0], STATIC, nothing, nothing,
       0.0, 0.0, 0.0, 0.0, 0.0]

segments:
  headers: [name, point_i, point_j, l0, diameter_mm,
            unit_stiffness, unit_damping, compression_frac]
  data:
    - [tether_seg, weight, ground, 50.0, 1.0,
       0.0, 0.0, 0.0]

tethers:
  headers: [name, segment_idxs]
  data:
    - [main_tether, [tether_seg]]

winches:
  headers: [name, tether_idxs, winch_point]
  data:
    - [main_winch, [main_tether], ground]
"""

const WINCH_TEST_SETTINGS = """
system:
    log_file: "data/winch_test"
    g_earth: 9.81

solver:
    solver: "FBDF"
    abs_tol: 0.0001
    rel_tol: 0.0001
    relaxation: 0.6

kite:
    model: ""
    foil_file: "ram_air_kite/ram_air_kite_foil.dat"
    physical_model: "2plate"
    mass: 0.0

tether:
    cd_tether: 0.0
    unit_damping: 350.0
    unit_stiffness: 120000.0
    rho_tether: 724.0
    e_tether: 1000.0
    rel_damping: 0.001
    d_tether: 1.0

winch:
    winch_model: "TorqueControlledMachine"
    max_force: 4000
    v_ro_max: 8.0
    drum_radius: 0.1
    gear_ratio: 1.0
    inertia_total: 0.1
    f_coulomb: 0.0
    c_vf: 0.0

environment:
    rho_0: 1.225
    v_wind: 0.0
    upwind_dir: -90.0
    upwind_elevation: 0.0
    wind_vec: [0.0, 0.0, 0.0]
    h_ref: 6.0
    profile_law: 0
"""

@testset "Winch Tests" begin
    # --- Setup: write YAML files, build model ONCE ---
    tmpdir = mktempdir()
    yaml_path = joinpath(tmpdir, "winch_test_geometry.yaml")
    write(yaml_path, WINCH_TEST_YAML)

    settings_path = joinpath(tmpdir, "settings.yaml")
    write(settings_path, WINCH_TEST_SETTINGS)

    system_yaml = "system:\n  sim_settings: settings.yaml\n"
    write(joinpath(tmpdir, "system.yaml"), system_yaml)

    set_data_path(tmpdir)
    set = Settings("system.yaml")

    sys, transforms = load_authoring(
        yaml_path; name="winch_test", set=set)
    sam = SymbolicAWEModel(set, sys; transforms)
    test_init!(sam)

    winch(quantity) = run_value(sam, :winches, :main_winch, quantity)
    set_winch!(quantity, value) =
        set_run_value!(sam, :winches, :main_winch, quantity, value)
    seg(quantity) = run_value(sam, :segments, :tether_seg, quantity)
    set_seg!(quantity, value) =
        set_run_value!(sam, :segments, :tether_seg, quantity, value)
    r = winch(:drum_radius)   # 0.1 m
    n = winch(:gear_ratio)     # 1.0

    # ============================================================
    # Test 1: Brake holds tether velocity at zero
    # ============================================================
    @testset "Brake holds" begin
        place!(sam)
        test_init!(sam; prn=false)
        set_winch!(:coulomb_friction, 0.0)
        set_winch!(:viscous_coefficient, 0.0)
        set_winch!(:inertia_total, 0.1)
        set_winch!(:brake, true)
        set_seg!(:unit_stiffness, 0.0)
        set_seg!(:unit_damping, 0.0)

        tau_motor = 5.0
        for _ in 1:100
            next_step!(sam; set_values=[tau_motor],
                       dt=0.001, vsm_interval=0)
        end
        @test abs(winch(:vel)) < 1e-10
    end

    # ============================================================
    # Test 2: Acceleration vs inertia
    # Zero stiffness, zero friction. Only motor torque and inertia.
    # Expected: a = (r/n) * tau_motor / I
    # ============================================================
    @testset "Acceleration vs inertia" begin
        tau_motor = 1.0

        for I_test in [0.1, 0.5, 1.0]
            place!(sam)
            test_init!(sam; prn=false)
            set_winch!(:coulomb_friction, 0.0)
            set_winch!(:viscous_coefficient, 0.0)
            set_winch!(:inertia_total, I_test)
            set_winch!(:brake, false)
            set_seg!(:unit_stiffness, 0.0)
            set_seg!(:unit_damping, 0.0)

            dt = 0.001
            n_steps = 20
            for _ in 1:n_steps
                next_step!(sam; set_values=[tau_motor],
                           dt=dt, vsm_interval=0)
            end

            # v(t) = a*t, so a = v / t
            a_measured = winch(:vel) / (n_steps * dt)
            a_expected = (r / n) * tau_motor / I_test

            @test a_measured ≈ a_expected rtol=0.05
        end
    end

    # ============================================================
    # Test 3: Acceleration vs Coulomb friction
    # Zero stiffness, zero viscous friction, small epsilon.
    # After spin-up (so smooth_sign ≈ 1):
    #   a = (r/n)/I * (tau_motor - coulomb_friction * r/n)
    # ============================================================
    @testset "Acceleration vs Coulomb friction" begin
        sam.definition.winches[:main_winch].model.friction_epsilon = 0.01
        tau_motor = 1.0
        I = 0.1

        for f_c in [0.5, 2.0, 5.0]
            place!(sam)
            test_init!(sam; prn=false)
            set_winch!(:viscous_coefficient, 0.0)
            set_winch!(:inertia_total, 0.1)
            set_winch!(:coulomb_friction, f_c)
            set_winch!(:brake, false)
            set_seg!(:unit_stiffness, 0.0)
            set_seg!(:unit_damping, 0.0)

            dt = 0.001
            # Spin up so smooth_sign(ω, 0.01) ≈ 1
            for _ in 1:50
                next_step!(sam; set_values=[tau_motor],
                           dt=dt, vsm_interval=0)
            end
            v0 = winch(:vel)

            # Measure acceleration over next steps
            n_meas = 20
            for _ in 1:n_meas
                next_step!(sam; set_values=[tau_motor],
                           dt=dt, vsm_interval=0)
            end
            v1 = winch(:vel)

            a_measured = (v1 - v0) / (n_meas * dt)
            a_expected = (r / n) / I *
                (tau_motor - f_c * r / n)

            @test a_measured ≈ a_expected rtol=0.05
        end
    end

    # ============================================================
    # Test 4: Terminal velocity vs viscous friction
    # Zero stiffness, zero Coulomb. Motor torque balanced by
    # viscous friction at terminal velocity.
    # Expected: v_term = tau_motor * n / (viscous_coefficient * r)
    # ============================================================
    @testset "Terminal velocity vs viscous friction" begin
        sam.definition.winches[:main_winch].model.friction_epsilon = 6.0
        tau_motor = 1.0
        I = 0.01

        for viscous_test in [50.0, 100.0, 200.0]
            place!(sam)
            test_init!(sam; prn=false)
            set_winch!(:coulomb_friction, 0.0)
            set_winch!(:inertia_total, 0.01  # Small I for fast settling)
            set_winch!(:viscous_coefficient, viscous_test)
            set_winch!(:brake, false)
            set_seg!(:unit_stiffness, 0.0)
            set_seg!(:unit_damping, 0.0)

            # Time constant: I / (viscous_coefficient * (r/n)^2)
            tau = I / (viscous_test * (r / n)^2)
            t_settle = 10 * tau
            dt = 0.001
            n_steps = Int(ceil(t_settle / dt))

            for _ in 1:n_steps
                next_step!(sam; set_values=[tau_motor],
                           dt=dt, vsm_interval=0)
            end

            v_expected = tau_motor * n / (viscous_test * r)
            @test winch(:vel) ≈ v_expected rtol=0.05
        end
    end

    # ============================================================
    # Test 5: Steady-state tether stretch
    # Brake on (l0 constant), weight sags under gravity.
    # At equilibrium: (unit_stiffness / len) * (len - l0) = m*g
    # Exact: extension = m*g*l0 / (unit_stiffness - m*g)
    # ============================================================
    @testset "Steady-state tether stretch" begin
        sam.definition.winches[:main_winch].model.friction_epsilon = 6.0

        unit_k = 50000.0  # Low but not too low [N]
        unit_d = unit_k * 0.01  # Damping for settling
        place!(sam)
        test_init!(sam; prn=false)
        set_winch!(:coulomb_friction, 1.0)
        set_winch!(:viscous_coefficient, 0.5)
        set_winch!(:inertia_total, 0.1)
        set_winch!(:brake, true)
        set_seg!(:unit_stiffness, unit_k)
        set_seg!(:unit_damping, unit_d)

        # Compute total mass at weight point:
        # m = extra_mass + rho * pi * (d/2)^2 * l0 / 2
        l0 = 50.0
        d = seg(:diameter)  # meters
        rho = seg(:density)
        extra_m = sam.definition.points[:weight].extra_mass
        tether_m = rho * π * (d / 2)^2 * l0 / 2
        m = extra_m + tether_m
        g = set.g_earth

        # Run to steady state
        dt = 0.001
        for _ in 1:5000
            next_step!(sam; dt=dt, vsm_interval=0)
        end

        # Exact equilibrium: k * ext = m*g where k = EA/len
        extension = m * g * l0 / (unit_k - m * g)
        expected_z = -(l0 + extension)

        weight_z = run_value(sam, :points, :weight, :pos_w)[3]
        @test weight_z ≈ expected_z rtol=0.01
    end

    # ============================================================
    # Test 6: calc_steady_torque at terminal velocity
    # Stiff tether with heavy weight, viscous friction,
    # small motor torque. At terminal velocity the system
    # is in steady state so calc_steady_torque must equal
    # the motor torque we are applying.
    # ============================================================
    @testset "calc_steady_torque at terminal velocity" begin
        sam.definition.winches[:main_winch].model.friction_epsilon = 6.0

        # Stiff tether to transmit gravity force
        place!(sam)
        test_init!(sam; prn=false)
        set_winch!(:coulomb_friction, 0.0)
        set_winch!(:inertia_total, 0.01)
        set_winch!(:viscous_coefficient, 100.0)
        set_winch!(:brake, false)
        set_seg!(:unit_stiffness, 120000.0)
        set_seg!(:unit_damping, 500.0)

        tau_motor = 0.5

        # Time constant: I / (viscous_coefficient * (r/n)^2)
        tau_tc = winch(:inertia_total) /
            (winch(:viscous_coefficient) * (r / n)^2)
        t_settle = 100 * tau_tc
        dt = 0.1
        n_steps = Int(ceil(t_settle / dt))

        for _ in 1:n_steps
            next_step!(sam; set_values=[tau_motor],
                       dt=dt, vsm_interval=0)
        end

        steady = calc_steady_torque(sam)
        @test run_value(sam, :winches, 1, :acc) ≈ 0.0 atol=0.01
        @test steady[1] ≈ tau_motor rtol=0.01

        println(
            "  calc_steady_torque: " *
            "applied=$(tau_motor), " *
            "computed=$(round(steady[1]; digits=4))"
        )
    end

    # ============================================================
    # Test 7: Route 2 tether auto-generation
    # Same weight-on-string setup but using Route 2 (only
    # start_point, end_point, n_segments — no explicit segments).
    # Verify: correct number of points/segments created, mass
    # falls under gravity to expected steady-state position.
    # ============================================================
    @testset "Route 2 auto-generated tether" begin
        using SymbolicAWEModels: SystemDefinition, SymbolicAWEModel, KiteGeometry, init!

        materials = KiteGeometry.MaterialTable(KiteGeometry.Material[],
                                               Dict{String, Any}())
        points = KiteGeometry.Point[
            Point(; name=:mass, pos_ENU=[0.0, 0.0, -100.0],
                  extra_mass=5.0),
            Point(; name=:anchor, pos_ENU=[0.0, 0.0, 0.0], fixed=true),
        ]
        segments = KiteGeometry.Segment[]
        transforms = [
            Transform(; name=:tf, elevation=deg2rad(-80.0), azimuth=0.0, heading=0.0,
                base_pos=[0, 0, 0], base_point=2, rot_point=[1 => 1.0])
        ]

        # Route 2: auto-generate 4 segments between mass
        # and anchor
        tethers = [KiteGeometry.expand_tether!(points, segments,
            fill(nothing, length(points)), :line, :mass, :anchor, 4,
            set, materials; model=KiteGeometry.NoModel(),
            segment_model=SegmentSpring())]
        winches = [Winch(; name=:winch, tethers=[:line],
            winch_point=:anchor)]

        sys2 = system_definition("route2_test", set;
            points, segments, materials, tethers, winches)

        # 2 original + 3 intermediate = 5 points
        @test length(sys2.points) == 5
        # 4 auto-generated segments
        @test length(sys2.segments) == 4
        # Segment names follow convention
        @test KiteGeometry.row_index(sys2.segments, :line_seg_1) == 1
        @test KiteGeometry.row_index(sys2.segments, :line_seg_4) == 4

        # Tether start/end resolved correctly
        teth = sys2.tethers[:line]
        @test teth.start_point ==
            KiteGeometry.row_index(sys2.points, :mass)
        @test teth.end_point ==
            KiteGeometry.row_index(sys2.points, :anchor)
        @test length(teth.segments) == 4

        # Winch point is anchor
        @test sys2.winches[:winch].winch_point ==
            KiteGeometry.row_index(sys2.points, :anchor)

        # Intermediate points evenly spaced
        for i in 1:3
            pt = sys2.points[Symbol("line_point_$i")]
            expected_z = -100.0 + i * 25.0
            @test pt.pos_ENU[3] ≈ expected_z atol=0.1
        end

        # Segments have correct l0 (100m / 4 = 25m each)
        for i in 1:4
            s = sys2.segments[Symbol("line_seg_$i")]
            @test s.rest_length ≈ 25.0 atol=0.01
        end

        # Build and simulate: mass should settle under gravity
        sam2 = SymbolicAWEModel(set, sys2; transforms)
        test_init!(sam2)

        # Brake on, stiff tether
        place!(sam2)
        test_init!(sam2; prn=false)
        set_run_value!(sam2, :winches, :winch, :brake, true)
        for k in eachindex(sys2.segments)
            set_run_value!(sam2, :segments, k, :unit_stiffness, 50000.0)
            set_run_value!(sam2, :segments, k, :unit_damping, 500.0)
        end

        next_step!(sam2; dt=3.0, vsm_interval=0)

        # Exact equilibrium: each segment carries weight of
        # all mass below its upper end.
        # extension_j = T_j * l0 / (unit_k - T_j)
        l0 = 25.0
        g = set.g_earth
        unit_k = 50000.0
        d = run_value(sam2, :segments, 1, :diameter)
        rho = run_value(sam2, :segments, 1, :density)
        half_seg_m = rho * π * (d / 2)^2 * l0 / 2
        m_mass = 5.0 + half_seg_m         # bottom point
        m_mid = 2 * half_seg_m            # intermediate pts

        # Cumulative weight below each segment (bottom-up)
        cum_mass = [m_mass,
                    m_mass + m_mid,
                    m_mass + 2 * m_mid,
                    m_mass + 3 * m_mid]
        total_len = sum(
            l0 * unit_k / (unit_k - m * g)
            for m in cum_mass)
        expected_z = -total_len

        mass_z = run_value(sam2, :points, :mass, :pos_w)[3]
        @test mass_z ≈ expected_z rtol=0.01
    end

    # ============================================================
    # Test 8: Tether without winch (constant l0)
    # A tether with segments but no winch. Segments should use
    # constant rest length. Verify the system builds and
    # simulates without errors.
    # ============================================================
    @testset "Tether without winch" begin
        using SymbolicAWEModels: SystemDefinition, SymbolicAWEModel, KiteGeometry

        materials = KiteGeometry.MaterialTable(KiteGeometry.Material[],
                                               Dict{String, Any}())
        points = [
            Point(; name=:top, pos_ENU=[0.0, 0.0, -50.0],
                  extra_mass=5.0),
            Point(; name=:bot, pos_ENU=[0.0, 0.0, 0.0], fixed=true),
        ]
        segments = [
            KiteGeometry.line_segment(:seg, (:top, :bot), set, materials;
                rest_length=50.0, diameter=0.001, unit_stiffness=50000.0,
                model=SegmentSpring(; unit_damping=500.0))
        ]
        tethers = [Tether(; name=:free_tether, segments=[:seg],
                          init_stretched_length=50.0)]
        transforms = [
            Transform(; name=:tf, elevation=deg2rad(-80.0), azimuth=0.0, heading=0.0,
                base_pos=[0, 0, 0], base_point=2, rot_point=[1 => 1.0])
        ]

        # No winches — should build fine
        sys3 = system_definition("no_winch", set;
            points, segments, materials, tethers)

        @test length(sys3.winches) == 0
        @test length(sys3.tethers) == 1
        teth = sys3.tethers[:free_tether]
        @test teth.start_point ==
            KiteGeometry.row_index(sys3.points, :top)
        @test teth.end_point ==
            KiteGeometry.row_index(sys3.points, :bot)

        sam3 = SymbolicAWEModel(set, sys3; transforms)
        test_init!(sam3)

        # Settle under gravity. k = unit_stiffness/l0 = 1000 N/m against 5 kg
        # damped by c = unit_damping/l0 = 10 Ns/m, so the transient decays with
        # a 1 s envelope over a 0.45 s period. 5 s leaves under a percent of it,
        # and FBDF picks its own substeps, so dt only has to resolve the period.
        for _ in 1:250
            next_step!(sam3; dt=0.02, vsm_interval=0)
        end

        # Should reach equilibrium (not diverge)
        top_z = run_value(sam3, :points, :top, :pos_w)[3]
        @test isfinite(top_z)
        @test top_z < -49.0  # stretched slightly by gravity

        # Its length is a parameter, so retrimming it mid-run takes effect on
        # the next step, with no reinit! to push the change onto the problem.
        set_run_value!(sam3, :tethers, :free_tether, :len, 48.0)
        for _ in 1:250
            next_step!(sam3; dt=0.02, vsm_interval=0)
        end
        @test run_value(sam3, :segments, :seg, :rest_length) ≈ 48.0 rtol=1e-6
        @test run_value(sam3, :points, :top, :pos_w)[3] > top_z
    end

    # ============================================================
    # Test 9: Custom length-tracking component
    # A custom winch component that interprets `set_value` as
    # target tether length and uses PD on `len` to track it.
    # Exercises the new `len` connector.
    # ============================================================
    @testset "Custom length-tracking component" begin
        using SymbolicAWEModels: SystemDefinition, SymbolicAWEModel, KiteGeometry,
            validate_winch_component

        materials_l = KiteGeometry.MaterialTable(KiteGeometry.Material[],
                                                 Dict{String, Any}())
        points_l = [
            Point(; name=:ground_l, pos_ENU=[0.0, 0.0, 0.0], fixed=true),
            Point(; name=:mass_l, pos_ENU=[0.0, 0.0, -50.0],
                  extra_mass=10.0),
        ]
        segments_l = [
            KiteGeometry.line_segment(:line_l, (:ground_l, :mass_l), set,
                materials_l; rest_length=50.0, diameter=0.005,
                unit_stiffness=50000.0, model=SegmentSpring(; unit_damping=500.0)),
        ]
        tethers_l = [Tether(; name=:main_l, segments=[:line_l],
                            init_stretched_length=50.0)]
        winches_l = [Winch(; name=:winch_l, tethers=[:main_l],
                           winch_point=:ground_l,
                           model=WinchModel(;
                               dynamics=LengthTrackWinch(100.0, 20.0)))]

        sys_l = system_definition("len_track_test", set;
            points=points_l, segments=segments_l, materials=materials_l,
            tethers=tethers_l, winches=winches_l)
        sam_l = SymbolicAWEModel(set, sys_l)
        place!(sam_l)
        test_init!(sam_l; prn=false)


        # Step 1: hold target = 50 for 1s (steady)
        for _ in 1:200
            next_step!(sam_l; set_values=[50.0],
                       dt=0.005, vsm_interval=0)
        end
        @test run_value(sam_l, :tethers, :main_l, :len) ≈ 50.0 atol=0.2

        # Step 2: target = 45 for 2s, expect settle
        for _ in 1:400
            next_step!(sam_l; set_values=[45.0],
                       dt=0.005, vsm_interval=0)
        end
        @test run_value(sam_l, :tethers, :main_l, :len) ≈ 45.0 atol=0.2

        # Step 3: target = 55 for 6s, expect settle. The 10 m of reel-out rings
        # the loop longer than the 5 m step above does; the P law's load droop
        # leaves the same ~0.1 m offset the other two steps settle at.
        for _ in 1:1200
            next_step!(sam_l; set_values=[55.0],
                       dt=0.005, vsm_interval=0)
        end
        @test run_value(sam_l, :tethers, :main_l, :len) ≈ 55.0 atol=0.2

        # Validation: a component missing `len` must fail
        bad_builder = function (sam, winch_idx; name)
            @variables begin
                vel(tt); force(tt); set_value(tt); brake(tt)
                acc(tt); friction(tt)
            end
            eqs = [
                friction ~ 0.0
                acc      ~ 0.0
            ]
            return System(eqs, tt; name)
        end
        bad_subsys = bad_builder(sys_l, 1; name=:bad)
        @test_throws ErrorException validate_winch_component(
            bad_subsys, sys_l.winches[1])
    end

    # ============================================================
    # Test 10: set_unstretched_length! hands a winched tether a new
    # rest length without re-placing the structure
    # ============================================================
    @testset "set_unstretched_length! retrims a winched tether" begin
        place!(sam)
        test_init!(sam; prn=false)
        set_winch!(:brake, true)
        set_seg!(:unit_stiffness, 50000.0)
        set_seg!(:unit_damping, 500.0)

        placed = [run_value(sam, :points, i, :pos_w) for i in eachindex(sys.points)]
        @test unstretched_length(sam) ≈ [50.0]

        set_unstretched_length!(sam, :main_tether, 48.0)
        # A winched tether's length is integrator state, so it reaches `u0`
        # through a fresh `init!`, which keeps the placement.
        init!(sam; remake=false, prn=false)

        @test unstretched_length(sam) ≈ [48.0]
        @test seg(:rest_length) ≈ 48.0
        @test all(run_value(sam, :points, i, :pos_w) ≈ was
                  for (i, was) in enumerate(placed))
    end

    rm(tmpdir; recursive=true)
end
nothing
