# SPDX-FileCopyrightText: 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

# test_steering_config.jl - Steering and depower configuration tests
#
# Tests SteeringConfig construction, reference resolution,
# base length capture, and set_steering!/set_depower! logic
# for both shared-line and dedicated depower modes.

using Test
using SymbolicAWEModels
using SymbolicAWEModels: apply_steering_config!,
    capture_steering_base_lengths!,
    resolve_steering_config!, load_sys_struct_from_yaml
using KiteUtils

# Minimal YAML with three steering/depower segments
const STEER_TEST_YAML = """
materials:
  headers: [name, youngs_modulus, density,
            damping_per_stiffness]
  data:
    - [test_mat, 1000.0, 724, 0.001]

points:
  headers: [name, pos_cad, type, wing_idx,
            transform_idx, extra_mass,
            body_frame_damping, world_frame_damping,
            area, drag_coeff]
  data:
    - [ground, [0, 0, 0], STATIC, nothing, nothing,
       0, 0, 0, 0, 0]
    - [kcu, [0, 0, -50], DYNAMIC, nothing, nothing,
       5, 0, 0, 0, 0]
    - [left, [0, -5, -55], DYNAMIC, nothing, nothing,
       1, 0, 0, 0, 0]
    - [right, [0, 5, -55], DYNAMIC, nothing, nothing,
       1, 0, 0, 0, 0]

segments:
  headers: [name, point_i, point_j, l0, diameter_mm,
            unit_stiffness, unit_damping, compression_frac]
  data:
    - [power_seg, ground, kcu, 50.0, 2.0,
       120000.0, 350.0, 0.0]
    - [steer_left, kcu, left, 10.0, 1.0,
       120000.0, 350.0, 0.0]
    - [steer_right, kcu, right, 12.0, 1.0,
       120000.0, 350.0, 0.0]

tethers:
  headers: [name, segment_idxs]
  data:
    - [main_tether, [power_seg]]
    - [left_tether, [steer_left]]
    - [right_tether, [steer_right]]

winches:
  headers: [name, tether_idxs, winch_point]
  data:
    - [main_winch, [main_tether, left_tether,
       right_tether], ground]

steering:
  steer_segments: [steer_left, steer_right]
  steer_gain: 2.5
"""

# Same geometry but with a dedicated depower segment
const STEER_DEPOWER_YAML = """
materials:
  headers: [name, youngs_modulus, density,
            damping_per_stiffness]
  data:
    - [test_mat, 1000.0, 724, 0.001]

points:
  headers: [name, pos_cad, type, wing_idx,
            transform_idx, extra_mass,
            body_frame_damping, world_frame_damping,
            area, drag_coeff]
  data:
    - [ground, [0, 0, 0], STATIC, nothing, nothing,
       0, 0, 0, 0, 0]
    - [kcu, [0, 0, -50], DYNAMIC, nothing, nothing,
       5, 0, 0, 0, 0]
    - [left, [0, -5, -55], DYNAMIC, nothing, nothing,
       1, 0, 0, 0, 0]
    - [right, [0, 5, -55], DYNAMIC, nothing, nothing,
       1, 0, 0, 0, 0]

segments:
  headers: [name, point_i, point_j, l0, diameter_mm,
            unit_stiffness, unit_damping, compression_frac]
  data:
    - [power_seg, ground, kcu, 50.0, 2.0,
       120000.0, 350.0, 0.0]
    - [steer_left, kcu, left, 10.0, 1.0,
       120000.0, 350.0, 0.0]
    - [steer_right, kcu, right, 12.0, 1.0,
       120000.0, 350.0, 0.0]

tethers:
  headers: [name, segment_idxs]
  data:
    - [main_tether, [power_seg]]
    - [left_tether, [steer_left]]
    - [right_tether, [steer_right]]

winches:
  headers: [name, tether_idxs, winch_point]
  data:
    - [main_winch, [main_tether, left_tether,
       right_tether], ground]

steering:
  steer_segments: [steer_left, steer_right]
  steer_gain: 2.5
  depower_segment: power_seg
  depower_gain: 3.0
"""

const STEER_TEST_SETTINGS = """
system:
    log_file: "data/steer_test"
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
    struc_geometry_path: "refine_struc_geometry.yaml"
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
    h_ref: 6.0
    profile_law: 0
"""

function build_steer_sys(yaml_string)
    tmpdir = mktempdir()
    yaml_path = joinpath(tmpdir, "geometry.yaml")
    write(yaml_path, yaml_string)
    settings_path = joinpath(tmpdir, "settings.yaml")
    write(settings_path, STEER_TEST_SETTINGS)
    system_yaml = "system:\n  sim_settings: settings.yaml\n"
    write(joinpath(tmpdir, "system.yaml"), system_yaml)
    set_data_path(tmpdir)
    set = Settings("system.yaml")
    return load_sys_struct_from_yaml(
        yaml_path; system_name="steer_test", set)
end

@testset "SteeringConfig" begin

    # ============================================================
    # Construction
    # ============================================================
    @testset "Constructor with symbol refs" begin
        cfg = SteeringConfig(
            steer_left=:seg_l, steer_right=:seg_r,
            steer_gain=2.0)
        @test cfg.steer_gain == 2.0
        @test cfg.depower_gain == 1.0
        @test cfg.depower_idx == 0
        @test cfg.steering == 0.0
        @test cfg.depower == 0.0
        @test cfg.steer_left_ref == :seg_l
        @test cfg.steer_right_ref == :seg_r
        @test isnothing(cfg.depower_ref)
    end

    @testset "Constructor with depower segment" begin
        cfg = SteeringConfig(
            steer_left=:sl, steer_right=:sr,
            steer_gain=1.5,
            depower_segment=:dp, depower_gain=3.0)
        @test cfg.steer_gain == 1.5
        @test cfg.depower_gain == 3.0
        @test cfg.depower_ref == :dp
    end

    @testset "Constructor with integer refs" begin
        cfg = SteeringConfig(
            steer_left=2, steer_right=3,
            steer_gain=1.0)
        @test cfg.steer_left_ref == 2
        @test cfg.steer_right_ref == 3
    end

    # ============================================================
    # Reference resolution
    # ============================================================
    @testset "resolve_steering_config!" begin
        names = Dict(:seg_l => Int64(2),
                     :seg_r => Int64(3),
                     :dp => Int64(1))

        cfg = SteeringConfig(
            steer_left=:seg_l, steer_right=:seg_r,
            steer_gain=1.0,
            depower_segment=:dp, depower_gain=2.0)
        resolve_steering_config!(cfg, names)

        @test cfg.steer_left_idx == 2
        @test cfg.steer_right_idx == 3
        @test cfg.depower_idx == 1
    end

    @testset "resolve with no depower" begin
        names = Dict(:sl => Int64(5), :sr => Int64(6))
        cfg = SteeringConfig(
            steer_left=:sl, steer_right=:sr,
            steer_gain=1.0)
        resolve_steering_config!(cfg, names)

        @test cfg.steer_left_idx == 5
        @test cfg.steer_right_idx == 6
        @test cfg.depower_idx == 0
    end

    # ============================================================
    # Shared-line mode (no dedicated depower segment)
    # ============================================================
    @testset "Shared-line steering" begin
        sys = build_steer_sys(STEER_TEST_YAML)
        cfg = sys.steering_config
        @test !isnothing(cfg)
        @test cfg.steer_gain == 2.5
        @test cfg.depower_idx == 0

        l0_left = cfg.l0_steer_left
        l0_right = cfg.l0_steer_right
        @test l0_left == 10.0
        @test l0_right == 12.0

        # Neutral steering: no change
        @test sys.segments[:steer_left].l0 == l0_left
        @test sys.segments[:steer_right].l0 == l0_right

        # Full right turn: left longer, right shorter
        set_steering!(sys, 1.0)
        @test cfg.steering == 1.0
        @test sys.segments[:steer_left].l0 ≈
            l0_left + 2.5
        @test sys.segments[:steer_right].l0 ≈
            l0_right - 2.5

        # Full left turn: left shorter, right longer
        set_steering!(sys, -1.0)
        @test sys.segments[:steer_left].l0 ≈
            l0_left - 2.5
        @test sys.segments[:steer_right].l0 ≈
            l0_right + 2.5

        # Half right steering
        set_steering!(sys, 0.5)
        @test sys.segments[:steer_left].l0 ≈
            l0_left + 1.25
        @test sys.segments[:steer_right].l0 ≈
            l0_right - 1.25

        # Back to neutral
        set_steering!(sys, 0.0)
        @test sys.segments[:steer_left].l0 ≈ l0_left
        @test sys.segments[:steer_right].l0 ≈ l0_right

        # Depower on shared lines
        set_depower!(sys, 1.0)
        @test cfg.depower == 1.0
        @test sys.segments[:steer_left].l0 ≈
            l0_left + cfg.depower_gain
        @test sys.segments[:steer_right].l0 ≈
            l0_right + cfg.depower_gain

        # Combined steering + depower
        set_steering!(sys, 0.5)
        set_depower!(sys, 0.4)
        dp_offset = 0.4 * cfg.depower_gain
        @test sys.segments[:steer_left].l0 ≈
            l0_left + 1.25 + dp_offset
        @test sys.segments[:steer_right].l0 ≈
            l0_right - 1.25 + dp_offset
    end

    # ============================================================
    # Dedicated depower segment mode
    # ============================================================
    @testset "Dedicated depower segment" begin
        sys = build_steer_sys(STEER_DEPOWER_YAML)
        cfg = sys.steering_config
        @test !isnothing(cfg)
        @test cfg.depower_gain == 3.0
        @test cfg.depower_idx > 0

        l0_left = cfg.l0_steer_left
        l0_right = cfg.l0_steer_right
        l0_depower = cfg.l0_depower
        @test l0_depower == 50.0

        # Steering only affects steer segments
        set_steering!(sys, 1.0)
        @test sys.segments[:steer_left].l0 ≈
            l0_left + 2.5
        @test sys.segments[:steer_right].l0 ≈
            l0_right - 2.5
        @test sys.segments[:power_seg].l0 ≈ l0_depower

        # Depower only affects dedicated segment
        set_steering!(sys, 0.0)
        set_depower!(sys, 0.5)
        @test sys.segments[:steer_left].l0 ≈ l0_left
        @test sys.segments[:steer_right].l0 ≈ l0_right
        @test sys.segments[:power_seg].l0 ≈
            l0_depower + 0.5 * 3.0

        # Full depower
        set_depower!(sys, 1.0)
        @test sys.segments[:power_seg].l0 ≈
            l0_depower + 3.0

        # Combined: independent
        set_steering!(sys, -0.5)
        set_depower!(sys, 0.8)
        @test sys.segments[:steer_left].l0 ≈
            l0_left - 1.25
        @test sys.segments[:steer_right].l0 ≈
            l0_right + 1.25
        @test sys.segments[:power_seg].l0 ≈
            l0_depower + 0.8 * 3.0
    end

    # ============================================================
    # Input clamping
    # ============================================================
    @testset "Input clamping" begin
        sys = build_steer_sys(STEER_TEST_YAML)
        cfg = sys.steering_config

        set_steering!(sys, 5.0)
        @test cfg.steering == 1.0

        set_steering!(sys, -10.0)
        @test cfg.steering == -1.0

        set_depower!(sys, 2.0)
        @test cfg.depower == 1.0

        set_depower!(sys, -0.5)
        @test cfg.depower == 0.0
    end

    # ============================================================
    # Error on missing config
    # ============================================================
    @testset "Error without config" begin
        yaml_no_steer = """
materials:
  headers: [name, youngs_modulus, density,
            damping_per_stiffness]
  data:
    - [test_mat, 1000.0, 724, 0.001]

points:
  headers: [name, pos_cad, type, wing_idx,
            transform_idx, extra_mass,
            body_frame_damping, world_frame_damping,
            area, drag_coeff]
  data:
    - [ground, [0, 0, 0], STATIC, nothing, nothing,
       0, 0, 0, 0, 0]
    - [mass, [0, 0, -50], DYNAMIC, nothing, nothing,
       5, 0, 0, 0, 0]

segments:
  headers: [name, point_i, point_j, l0, diameter_mm,
            unit_stiffness, unit_damping, compression_frac]
  data:
    - [seg, ground, mass, 50.0, 2.0,
       120000.0, 350.0, 0.0]

tethers:
  headers: [name, segment_idxs]
  data:
    - [tether, [seg]]

winches:
  headers: [name, tether_idxs, winch_point]
  data:
    - [winch, [tether], ground]
"""
        sys = build_steer_sys(yaml_no_steer)
        @test isnothing(sys.steering_config)
        @test_throws ErrorException set_steering!(sys, 0.5)
        @test_throws ErrorException set_depower!(sys, 0.5)
    end
end
