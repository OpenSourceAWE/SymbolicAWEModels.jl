# SPDX-FileCopyrightText: 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

# test_tether_init_len.jl - Tether initial length scaling tests
#
# Tests the Tether.init_len feature: scaling pos_w before transforms.
# Uses Route 2 (auto-generated) tethers and YAML-specified init_len.
# All tests use reinit! directly on a SystemStructure (no ODE compilation).

using Test
using SymbolicAWEModels
using SymbolicAWEModels: KVec3
using KiteUtils
using LinearAlgebra

# ================================================================
# Route 2 tether: ground→[auto mid]→top, 2 segments, init_len in YAML
# ================================================================
const INIT_LEN_YAML_ROUTE2 = """
materials:
  headers: [name, youngs_modulus, density, damping_per_stiffness]
  data:
    - [test_mat, 120000.0, 724, 0.001]

points:
  headers: [name, pos_cad, type, wing_idx, transform_idx,
            extra_mass, body_frame_damping, world_frame_damping,
            area, drag_coeff]
  data:
    - [ground, [0.0, 0.0, 0.0], STATIC, nothing, nothing,
       0.0, 0.0, 0.0, 0.0, 0.0]
    - [top, [0.0, 0.0, -100.0], DYNAMIC, nothing, nothing,
       1.0, 0.0, 0.0, 0.0, 0.0]

tethers:
  headers: [name, start_point, end_point, n_segments, material, init_len]
  data:
    - [main_tether, ground, top, 2, test_mat, 200.0]

winches:
  headers: [name, tether_idxs, winch_point]
  data:
    - [main_winch, [main_tether], ground]
"""

# Route 1 tether (explicit segments) with init_len in YAML
const INIT_LEN_YAML_ROUTE1 = """
materials:
  headers: [name, youngs_modulus, density, damping_per_stiffness]
  data:
    - [test_mat, 120000.0, 724, 0.001]

points:
  headers: [name, pos_cad, type, wing_idx, transform_idx,
            extra_mass, body_frame_damping, world_frame_damping,
            area, drag_coeff]
  data:
    - [ground, [0.0, 0.0, 0.0], STATIC, nothing, nothing,
       0.0, 0.0, 0.0, 0.0, 0.0]
    - [mid, [0.0, 0.0, -50.0], DYNAMIC, nothing, nothing,
       1.0, 0.0, 0.0, 0.0, 0.0]
    - [top, [0.0, 0.0, -100.0], DYNAMIC, nothing, nothing,
       1.0, 0.0, 0.0, 0.0, 0.0]

segments:
  headers: [name, point_i, point_j, l0, diameter_mm,
            unit_stiffness, unit_damping, compression_frac]
  data:
    - [s1, ground, mid, 50.0, 4.0, 120000.0, 350.0, 0.0]
    - [s2, mid, top, 50.0, 4.0, 120000.0, 350.0, 0.0]

tethers:
  headers: [name, segment_idxs, init_len]
  data:
    - [main_tether, [s1, s2], 200.0]

winches:
  headers: [name, tether_idxs, winch_point]
  data:
    - [main_winch, [main_tether], ground]
"""

# Route 2 with downstream point connected via bridle (non-tether segment)
const INIT_LEN_DOWNSTREAM_ROUTE2 = """
materials:
  headers: [name, youngs_modulus, density, damping_per_stiffness]
  data:
    - [test_mat, 120000.0, 724, 0.001]

points:
  headers: [name, pos_cad, type, wing_idx, transform_idx,
            extra_mass, body_frame_damping, world_frame_damping,
            area, drag_coeff]
  data:
    - [ground, [0.0, 0.0, 0.0], STATIC, nothing, nothing,
       0.0, 0.0, 0.0, 0.0, 0.0]
    - [top, [0.0, 0.0, -100.0], DYNAMIC, nothing, nothing,
       1.0, 0.0, 0.0, 0.0, 0.0]
    - [downstream, [10.0, 0.0, -100.0], DYNAMIC, nothing, nothing,
       1.0, 0.0, 0.0, 0.0, 0.0]

segments:
  headers: [name, point_i, point_j, l0, diameter_mm,
            unit_stiffness, unit_damping, compression_frac]
  data:
    - [bridle, top, downstream, 10.0, 4.0, 120000.0, 350.0, 0.0]

tethers:
  headers: [name, start_point, end_point, n_segments, material, init_len]
  data:
    - [main_tether, ground, top, 2, test_mat, 200.0]

winches:
  headers: [name, tether_idxs, winch_point]
  data:
    - [main_winch, [main_tether], ground]
"""

# Route 2 with loop (non-tether segment from top back to ground)
const INIT_LEN_LOOP_ROUTE2 = """
materials:
  headers: [name, youngs_modulus, density, damping_per_stiffness]
  data:
    - [test_mat, 120000.0, 724, 0.001]

points:
  headers: [name, pos_cad, type, wing_idx, transform_idx,
            extra_mass, body_frame_damping, world_frame_damping,
            area, drag_coeff]
  data:
    - [ground, [0.0, 0.0, 0.0], STATIC, nothing, nothing,
       0.0, 0.0, 0.0, 0.0, 0.0]
    - [top, [0.0, 0.0, -100.0], DYNAMIC, nothing, nothing,
       1.0, 0.0, 0.0, 0.0, 0.0]

segments:
  headers: [name, point_i, point_j, l0, diameter_mm,
            unit_stiffness, unit_damping, compression_frac]
  data:
    - [back_loop, top, ground, 100.0, 4.0, 120000.0, 350.0, 0.0]

tethers:
  headers: [name, start_point, end_point, n_segments, material]
  data:
    - [main_tether, ground, top, 2, test_mat]

winches:
  headers: [name, tether_idxs, winch_point]
  data:
    - [main_winch, [main_tether], ground]
"""

const INIT_LEN_SETTINGS = """
system:
    log_file: "data/init_len_test"
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
    e_tether: 120000.0
    rel_damping: 0.001
    d_tether: 4.0

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

@testset "Tether init_len Tests" begin
    tmpdir = mktempdir()
    write(joinpath(tmpdir, "settings.yaml"), INIT_LEN_SETTINGS)
    write(joinpath(tmpdir, "system.yaml"),
          "system:\n  sim_settings: settings.yaml\n")
    set_data_path(tmpdir)
    set = Settings("system.yaml")

    # ================================================================
    # Test 1: Route 2 auto-gen + YAML init_len → scale to 2×
    # ================================================================
    @testset "Route 2 YAML init_len: scale to 2x" begin
        yaml_path = joinpath(tmpdir, "r2_yaml.yaml")
        write(yaml_path, INIT_LEN_YAML_ROUTE2)
        sys = load_sys_struct_from_yaml(
            yaml_path; system_name="init_len_r2_yaml", set=set)

        # init_len=200 already set from YAML; no programmatic change needed
        SymbolicAWEModels.reinit!(sys, set)

        mid = sys.points[:main_tether_point_1]
        @test mid.pos_w ≈ KVec3(0, 0, -100)
        @test sys.points[:top].pos_w ≈ KVec3(0, 0, -200)
        @test sys.winches[:main_winch].tether_len ≈ 200.0
        # pos_cad unchanged
        @test mid.pos_cad ≈ KVec3(0, 0, -50)
        @test sys.points[:top].pos_cad ≈ KVec3(0, 0, -100)
    end

    # ================================================================
    # Test 2: Route 1 explicit segments + YAML init_len
    # ================================================================
    @testset "Route 1 YAML init_len: scale to 2x" begin
        yaml_path = joinpath(tmpdir, "r1_yaml.yaml")
        write(yaml_path, INIT_LEN_YAML_ROUTE1)
        sys = load_sys_struct_from_yaml(
            yaml_path; system_name="init_len_r1_yaml", set=set)

        SymbolicAWEModels.reinit!(sys, set)

        @test sys.points[:mid].pos_w ≈ KVec3(0, 0, -100)
        @test sys.points[:top].pos_w ≈ KVec3(0, 0, -200)
        @test sys.segments[:s1].l0 ≈ 100.0
        @test sys.segments[:s2].l0 ≈ 100.0
        @test sys.points[:top].pos_cad ≈ KVec3(0, 0, -100)
    end

    # ================================================================
    # Test 3: Route 2 + downstream point translation
    # ================================================================
    @testset "Route 2 downstream point translation" begin
        yaml_path = joinpath(tmpdir, "r2_downstream.yaml")
        write(yaml_path, INIT_LEN_DOWNSTREAM_ROUTE2)
        sys = load_sys_struct_from_yaml(
            yaml_path; system_name="init_len_r2_downstream", set=set)

        SymbolicAWEModels.reinit!(sys, set)

        @test sys.points[:top].pos_w ≈ KVec3(0, 0, -200)
        @test sys.points[:downstream].pos_w ≈ KVec3(10, 0, -200)
    end

    # ================================================================
    # Test 4: Error — downstream non-tether segment connects back to start
    # ================================================================
    @testset "Error on loop to start" begin
        yaml_path = joinpath(tmpdir, "r2_loop.yaml")
        write(yaml_path, INIT_LEN_LOOP_ROUTE2)
        sys = load_sys_struct_from_yaml(
            yaml_path; system_name="init_len_r2_loop", set=set)

        sys.tethers[:main_tether].init_len = 200.0
        @test_throws ErrorException SymbolicAWEModels.reinit!(sys, set)
    end

    # ================================================================
    # Test 5: Idempotency — repeated reinit! gives same result
    # ================================================================
    @testset "Idempotency" begin
        yaml_path = joinpath(tmpdir, "r2_yaml.yaml")
        sys = load_sys_struct_from_yaml(
            yaml_path; system_name="init_len_r2_idem", set=set)

        SymbolicAWEModels.reinit!(sys, set)
        mid_pos = copy(sys.points[:main_tether_point_1].pos_w)
        top_pos = copy(sys.points[:top].pos_w)

        SymbolicAWEModels.reinit!(sys, set)
        @test sys.points[:main_tether_point_1].pos_w ≈ mid_pos
        @test sys.points[:top].pos_w ≈ top_pos
    end

end
