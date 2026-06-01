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

# ================================================================
# Minimal YAML: weight at z=-50 connected to ground via tether
# ================================================================
const WINCH_TEST_YAML = """
materials:
  headers: [name, youngs_modulus, density, damping_per_stiffness]
  data:
    - [test_mat, 1000.0, 724, 0.001]

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
    struc_geometry_path: "particle_structural_geometry.yaml"
    aero_geometry_path: "aero_geometry.yaml"
    mass: 0.0
    quasi_static: false

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

    sys = load_sys_struct_from_yaml(
        yaml_path; system_name="winch_test", set=set)
    sam = SymbolicAWEModel(set, sys)
    test_init!(sam)

    winch = sam.sys_struct.winches[:main_winch]
    tether = sam.sys_struct.tethers[:main_tether]
    seg = sam.sys_struct.segments[:tether_seg]
    r = winch.drum_radius   # 0.1 m
    n = winch.gear_ratio     # 1.0

    # ============================================================
    # Test 1: Brake holds tether velocity at zero
    # ============================================================
    @testset "Brake holds" begin
        winch.brake = true
        winch.f_coulomb = 0.0
        winch.c_vf = 0.0
        winch.inertia_total = 0.1
        seg.unit_stiffness = 0.0
        seg.unit_damping = 0.0
        test_init!(sam; prn=false)

        tau_motor = 5.0
        for _ in 1:100
            next_step!(sam; set_values=[tau_motor],
                       dt=0.001, vsm_interval=0)
        end
        @test abs(winch.vel) < 1e-10
    end

    # ============================================================
    # Test 2: Acceleration vs inertia
    # Zero stiffness, zero friction. Only motor torque and inertia.
    # Expected: a = (r/n) * tau_motor / I
    # ============================================================
    @testset "Acceleration vs inertia" begin
        winch.brake = false
        winch.f_coulomb = 0.0
        winch.c_vf = 0.0
        seg.unit_stiffness = 0.0
        seg.unit_damping = 0.0
        tau_motor = 1.0

        for I_test in [0.1, 0.5, 1.0]
            winch.inertia_total = I_test
            test_init!(sam; prn=false)

            dt = 0.001
            n_steps = 20
            for _ in 1:n_steps
                next_step!(sam; set_values=[tau_motor],
                           dt=dt, vsm_interval=0)
            end

            # v(t) = a*t, so a = v / t
            a_measured = winch.vel / (n_steps * dt)
            a_expected = (r / n) * tau_motor / I_test

            @test a_measured ≈ a_expected rtol=0.05
        end
    end

    # ============================================================
    # Test 3: Acceleration vs Coulomb friction
    # Zero stiffness, zero viscous friction, small epsilon.
    # After spin-up (so smooth_sign ≈ 1):
    #   a = (r/n)/I * (tau_motor - f_coulomb * r/n)
    # ============================================================
    @testset "Acceleration vs Coulomb friction" begin
        winch.brake = false
        winch.c_vf = 0.0
        winch.friction_epsilon = 0.01
        winch.inertia_total = 0.1
        seg.unit_stiffness = 0.0
        seg.unit_damping = 0.0
        tau_motor = 1.0
        I = winch.inertia_total

        for f_c in [0.5, 2.0, 5.0]
            winch.f_coulomb = f_c
            test_init!(sam; prn=false)

            dt = 0.001
            # Spin up so smooth_sign(ω, 0.01) ≈ 1
            for _ in 1:50
                next_step!(sam; set_values=[tau_motor],
                           dt=dt, vsm_interval=0)
            end
            v0 = winch.vel

            # Measure acceleration over next steps
            n_meas = 20
            for _ in 1:n_meas
                next_step!(sam; set_values=[tau_motor],
                           dt=dt, vsm_interval=0)
            end
            v1 = winch.vel

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
    # Expected: v_term = tau_motor * n / (c_vf * r)
    # ============================================================
    @testset "Terminal velocity vs viscous friction" begin
        winch.brake = false
        winch.f_coulomb = 0.0
        winch.friction_epsilon = 6.0
        winch.inertia_total = 0.01  # Small I for fast settling
        seg.unit_stiffness = 0.0
        seg.unit_damping = 0.0
        tau_motor = 1.0
        I = winch.inertia_total

        for c_vf_test in [50.0, 100.0, 200.0]
            winch.c_vf = c_vf_test
            test_init!(sam; prn=false)

            # Time constant: I / (c_vf * (r/n)^2)
            tau = I / (c_vf_test * (r / n)^2)
            t_settle = 10 * tau
            dt = 0.001
            n_steps = Int(ceil(t_settle / dt))

            for _ in 1:n_steps
                next_step!(sam; set_values=[tau_motor],
                           dt=dt, vsm_interval=0)
            end

            v_expected = tau_motor * n / (c_vf_test * r)
            @test winch.vel ≈ v_expected rtol=0.05
        end
    end

    # ============================================================
    # Test 5: Steady-state tether stretch
    # Brake on (l0 constant), weight sags under gravity.
    # At equilibrium: (unit_stiffness / len) * (len - l0) = m*g
    # Exact: extension = m*g*l0 / (unit_stiffness - m*g)
    # ============================================================
    @testset "Steady-state tether stretch" begin
        winch.brake = true
        winch.f_coulomb = 1.0
        winch.c_vf = 0.5
        winch.friction_epsilon = 6.0
        winch.inertia_total = 0.1

        unit_k = 50000.0  # Low but not too low [N]
        unit_d = unit_k * 0.01  # Damping for settling
        seg.unit_stiffness = unit_k
        seg.unit_damping = unit_d
        test_init!(sam; prn=false)

        # Compute total mass at weight point:
        # m = extra_mass + rho * pi * (d/2)^2 * l0 / 2
        l0 = 50.0
        d = seg.diameter  # meters
        rho = set.rho_tether
        extra_m = sam.sys_struct.points[:weight].extra_mass
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

        weight_z = sam.sys_struct.points[:weight].pos_w[3]
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
        winch.brake = false
        winch.f_coulomb = 0.0
        winch.friction_epsilon = 6.0
        winch.inertia_total = 0.01
        winch.c_vf = 100.0

        # Stiff tether to transmit gravity force
        seg.unit_stiffness = 120000.0
        seg.unit_damping = 500.0
        test_init!(sam; prn=false)

        tau_motor = 0.5

        # Time constant: I / (c_vf * (r/n)^2)
        tau_tc = winch.inertia_total /
            (winch.c_vf * (r / n)^2)
        t_settle = 100 * tau_tc
        dt = 0.1
        n_steps = Int(ceil(t_settle / dt))

        for _ in 1:n_steps
            next_step!(sam; set_values=[tau_motor],
                       dt=dt, vsm_interval=0)
        end

        steady = calc_steady_torque(sam)
        @test sam.sys_struct.winches[1].acc ≈ 0.0 atol=0.01
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
        using SymbolicAWEModels: Point, Segment, Tether, Winch,
            Transform, SystemStructure, SymbolicAWEModel,
            STATIC, DYNAMIC, init!

        points = [
            Point(:mass, [0.0, 0.0, -100.0], DYNAMIC;
                  extra_mass=5.0),
            Point(:anchor, [0.0, 0.0, 0.0], STATIC),
        ]
        transforms = [
            Transform(:tf, deg2rad(-80.0), 0.0, 0.0;
                base_pos=[0, 0, 0], base_point=:anchor,
                rot_point=:mass)
        ]

        # Route 2: auto-generate 4 segments between mass
        # and anchor
        tethers = [Tether(:line;
            start_point=:mass, end_point=:anchor,
            n_segments=4)]
        winches = [Winch(:winch, set, [:line];
            winch_point=:anchor)]

        sys2 = SystemStructure("route2_test", set;
            points, tethers, winches, transforms,
            prn=false)

        # 2 original + 3 intermediate = 5 points
        @test length(sys2.points) == 5
        # 4 auto-generated segments
        @test length(sys2.segments) == 4
        # Segment names follow convention
        @test sys2.segments[:line_seg_1].idx == 1
        @test sys2.segments[:line_seg_4].idx == 4

        # Tether start/end resolved correctly
        teth = sys2.tethers[:line]
        @test teth.start_point_idx ==
            sys2.points[:mass].idx
        @test teth.end_point_idx ==
            sys2.points[:anchor].idx
        @test length(teth.segment_idxs) == 4

        # Winch point is anchor
        @test sys2.winches[:winch].winch_point_idx ==
            sys2.points[:anchor].idx

        # Intermediate points evenly spaced
        for i in 1:3
            pt = sys2.points[Symbol("line_point_$i")]
            expected_z = -100.0 + i * 25.0
            @test pt.pos_w[3] ≈ expected_z atol=0.1
        end

        # Segments have correct l0 (100m / 4 = 25m each)
        for i in 1:4
            s = sys2.segments[Symbol("line_seg_$i")]
            @test s.l0 ≈ 25.0 atol=0.01
        end

        # Build and simulate: mass should settle under gravity
        sam2 = SymbolicAWEModel(set, sys2)
        test_init!(sam2)

        # Brake on, stiff tether
        w2 = sam2.sys_struct.winches[:winch]
        w2.brake = true
        for s in sam2.sys_struct.segments
            s.unit_stiffness = 50000.0
            s.unit_damping = 500.0
        end
        test_init!(sam2; prn=false)

        next_step!(sam2; dt=3.0, vsm_interval=0)

        # Exact equilibrium: each segment carries weight of
        # all mass below its upper end.
        # extension_j = T_j * l0 / (unit_k - T_j)
        l0 = 25.0
        g = set.g_earth
        unit_k = 50000.0
        d = sam2.sys_struct.segments[1].diameter
        rho = set.rho_tether
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

        mass_z = sam2.sys_struct.points[:mass].pos_w[3]
        @test mass_z ≈ expected_z rtol=0.01
    end

    # ============================================================
    # Test 8: Tether without winch (constant l0)
    # A tether with segments but no winch. Segments should use
    # constant rest length (from get_l0). Verify the system
    # builds and simulates without errors.
    # ============================================================
    @testset "Tether without winch" begin
        using SymbolicAWEModels: Point, Segment, Tether,
            Transform, SystemStructure, SymbolicAWEModel,
            STATIC, DYNAMIC

        points = [
            Point(:top, [0.0, 0.0, -50.0], DYNAMIC;
                  extra_mass=5.0),
            Point(:bot, [0.0, 0.0, 0.0], STATIC),
        ]
        segments = [
            Segment(:seg, :top, :bot, 50000.0, 500.0,
                    0.001; l0=50.0)
        ]
        tethers = [Tether(:free_tether, [:seg], 50.0)]
        transforms = [
            Transform(:tf, deg2rad(-80.0), 0.0, 0.0;
                base_pos=[0, 0, 0], base_point=:bot,
                rot_point=:top)
        ]

        # No winches — should build fine
        sys3 = SystemStructure("no_winch", set;
            points, segments, tethers, transforms,
            prn=false)

        @test length(sys3.winches) == 0
        @test length(sys3.tethers) == 1
        teth = sys3.tethers[:free_tether]
        @test teth.start_point_idx ==
            sys3.points[:top].idx
        @test teth.end_point_idx ==
            sys3.points[:bot].idx

        sam3 = SymbolicAWEModel(set, sys3)
        test_init!(sam3)

        # Simulate: mass settles under gravity with stiff
        # constant-l0 tether
        for _ in 1:2000
            next_step!(sam3; dt=0.001, vsm_interval=0)
        end

        # Should reach equilibrium (not diverge)
        top_z = sam3.sys_struct.points[:top].pos_w[3]
        @test isfinite(top_z)
        @test top_z < -49.0  # stretched slightly by gravity
    end

    rm(tmpdir; recursive=true)
end
nothing
