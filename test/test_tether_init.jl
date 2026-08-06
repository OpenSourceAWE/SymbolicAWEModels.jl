# SPDX-FileCopyrightText: 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_tether_init.jl - Tether initial length scaling tests
#
# Tests the Tether.init_stretched_length feature: scaling pos_w before transforms.
# Uses Route 2 (auto-generated) tethers and YAML-specified init_stretched_length.
# All tests use reinit! directly on a SystemStructure (no ODE compilation).

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

using Test
using SymbolicAWEModels
using SymbolicAWEModels: KVec3, Point
using KiteUtils
using LinearAlgebra

# ================================================================
# Route 2 tether: ground→[auto mid]→top, 2 segments, init_stretched_length in YAML
# ================================================================
const INIT_LEN_YAML_ROUTE2 = """
variables:
  test_mat:
    unit_damping: 0.0
    density: 724.0

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
  headers: [name, start_point, end_point, n_segments, unit_damping,
            density, init_stretched_length]
  data:
    - [main_tether, ground, top, 2, test_mat, 200.0]

winches:
  headers: [name, tether_idxs, winch_point]
  data:
    - [main_winch, [main_tether], ground]
"""

# Route 1 tether (explicit segments) with init_stretched_length in YAML
const INIT_LEN_YAML_ROUTE1 = """
variables:
  test_mat:
    unit_damping: 0.0
    density: 724.0

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
  headers: [name, segment_idxs, init_stretched_length]
  data:
    - [main_tether, [s1, s2], 200.0]

winches:
  headers: [name, tether_idxs, winch_point]
  data:
    - [main_winch, [main_tether], ground]
"""

# Route 2 with downstream point connected via bridle (non-tether segment)
const INIT_LEN_DOWNSTREAM_ROUTE2 = """
variables:
  test_mat:
    unit_damping: 0.0
    density: 724.0

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
  headers: [name, start_point, end_point, n_segments, unit_damping,
            density, init_stretched_length]
  data:
    - [main_tether, ground, top, 2, test_mat, 200.0]

winches:
  headers: [name, tether_idxs, winch_point]
  data:
    - [main_winch, [main_tether], ground]
"""

# Route 2 with loop (non-tether segment from top back to ground)
const INIT_LEN_LOOP_ROUTE2 = """
variables:
  test_mat:
    unit_damping: 0.0
    density: 724.0

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
  headers: [name, start_point, end_point, n_segments, unit_damping,
            density, init_stretched_length]
  data:
    - [main_tether, ground, top, 2, test_mat, 100.0]

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
    struc_geometry_path: "particle_structural_geometry.yaml"
    aero_geometry_path: "aero_geometry.yaml"
    mass: 0.0

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
    upwind_elevation: 0.0
    wind_vec: [0.0, 0.0, 0.0]
    h_ref: 6.0
    profile_law: 0
"""

@testset "Tether init_stretched_length Tests" begin
    tmpdir = mktempdir()
    write(joinpath(tmpdir, "settings.yaml"), INIT_LEN_SETTINGS)
    write(joinpath(tmpdir, "system.yaml"),
          "system:\n  sim_settings: settings.yaml\n")
    set_data_path(tmpdir)
    set = Settings("system.yaml")

    # ================================================================
    # Test 1: Route 2 auto-gen + YAML init_stretched_length → scale to 2×
    # ================================================================
    @testset "Route 2 YAML init_stretched_length: scale to 2x" begin
        yaml_path = joinpath(tmpdir, "r2_yaml.yaml")
        write(yaml_path, INIT_LEN_YAML_ROUTE2)
        sys = load_sys_struct_from_yaml(
            yaml_path; system_name="init_stretched_length_r2_yaml", set=set)

        # init_stretched_length=200 already set from YAML; no programmatic change needed
        SymbolicAWEModels.reinit!(sys, set)

        mid = sys.points[:main_tether_point_1]
        @test mid.pos_w ≈ KVec3(0, 0, -100)
        @test sys.points[:top].pos_w ≈ KVec3(0, 0, -200)
        @test sys.tethers[:main_tether].len ≈ 200.0
        # pos_cad unchanged
        @test mid.pos_cad ≈ KVec3(0, 0, -50)
        @test sys.points[:top].pos_cad ≈ KVec3(0, 0, -100)
    end

    # ================================================================
    # Test 2: Route 1 explicit segments + YAML init_stretched_length
    # ================================================================
    @testset "Route 1 YAML init_stretched_length: scale to 2x" begin
        yaml_path = joinpath(tmpdir, "r1_yaml.yaml")
        write(yaml_path, INIT_LEN_YAML_ROUTE1)
        sys = load_sys_struct_from_yaml(
            yaml_path; system_name="init_stretched_length_r1_yaml", set=set)

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
            yaml_path; system_name="init_stretched_length_r2_downstream", set=set)

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
        @test_throws ErrorException load_sys_struct_from_yaml(
            yaml_path; system_name="init_stretched_length_r2_loop", set=set)
    end

    # ================================================================
    # Test 5: Idempotency — repeated reinit! gives same result
    # ================================================================
    @testset "Idempotency" begin
        yaml_path = joinpath(tmpdir, "r2_yaml.yaml")
        sys = load_sys_struct_from_yaml(
            yaml_path; system_name="init_stretched_length_r2_idem", set=set)

        SymbolicAWEModels.reinit!(sys, set)
        mid_pos = copy(sys.points[:main_tether_point_1].pos_w)
        top_pos = copy(sys.points[:top].pos_w)

        SymbolicAWEModels.reinit!(sys, set)
        @test sys.points[:main_tether_point_1].pos_w ≈ mid_pos
        @test sys.points[:top].pos_w ≈ top_pos
    end

    # ================================================================
    # Test 6: Multi-tether with a STATIC anchor and a winched DYNAMIC anchor
    # ================================================================
    @testset "Multi-tether: static + winch anchors stay fixed" begin
        multi_yaml = """
variables:
  test_mat:
    unit_damping: 0.0
    density: 724.0

points:
  headers: [name, pos_cad, type, wing_idx, transform_idx,
            extra_mass, body_frame_damping, world_frame_damping,
            area, drag_coeff]
  data:
    - [ground_static, [10.0, 0.0, 0.0], STATIC, nothing, nothing,
       0.0, 0.0, 0.0, 0.0, 0.0]
    - [ground_winch, [-10.0, 0.0, 0.0], DYNAMIC, nothing, nothing,
       1.0, 0.0, 0.0, 0.0, 0.0]
    - [top, [0.0, 0.0, -100.0], DYNAMIC, nothing, nothing,
       1.0, 0.0, 0.0, 0.0, 0.0]

tethers:
  headers: [name, start_point, end_point, n_segments, unit_damping, density]
  data:
    - [tether_static, ground_static, top, 2, test_mat]
    - [tether_winch, ground_winch, top, 2, test_mat]

winches:
  headers: [name, tether_idxs, winch_point]
  data:
    - [winch_b, [tether_winch], ground_winch]
"""
        yaml_path = joinpath(tmpdir, "multi.yaml")
        write(yaml_path, multi_yaml)
        sys = load_sys_struct_from_yaml(
            yaml_path; system_name="init_stretched_length_multi", set=set)

        sys.tethers[:tether_static].init_stretched_len = 200.0
        SymbolicAWEModels.reinit!(sys, set)

        ground_static = sys.points[:ground_static].pos_w
        @test norm(sys.points[:top].pos_w - ground_static) ≈ 200.0
        @test sys.points[:ground_static].pos_w ≈ KVec3(10, 0, 0)
        @test sys.points[:ground_winch].pos_w ≈ KVec3(-10, 0, 0)

        sys.tethers[:tether_winch].init_stretched_len = 100.0
        @test_logs (:info,) match_mode=:any SymbolicAWEModels.reinit!(sys, set)

        # Placed by the mean displacement of both roots: standoff is
        # ≈ the mean target (150), offset slightly because the two
        # tethers pull in different directions, and top is drawn off
        # the static-tether line toward that mean direction.
        ground_static = sys.points[:ground_static].pos_w
        @test isapprox(norm(sys.points[:top].pos_w - ground_static),
                       150.0; atol=0.05)
        @test sys.points[:top].pos_w[1] < -4.95
        @test sys.points[:ground_winch].pos_w ≈ KVec3(-10, 0, 0)
    end

    # ================================================================
    # Test 7: init_stretched_len on a non-root tether is an error
    # ================================================================
    @testset "Error on non-root init_stretched_len" begin
        stacked_yaml = """
variables:
  test_mat:
    unit_damping: 0.0
    density: 724.0

points:
  headers: [name, pos_cad, type, wing_idx, transform_idx,
            extra_mass, body_frame_damping, world_frame_damping,
            area, drag_coeff]
  data:
    - [ground, [0.0, 0.0, 0.0], STATIC, nothing, nothing,
       0.0, 0.0, 0.0, 0.0, 0.0]
    - [mid, [0.0, 0.0, -100.0], DYNAMIC, nothing, nothing,
       1.0, 0.0, 0.0, 0.0, 0.0]
    - [top, [0.0, 0.0, -200.0], DYNAMIC, nothing, nothing,
       1.0, 0.0, 0.0, 0.0, 0.0]

tethers:
  headers: [name, start_point, end_point, n_segments, unit_damping,
            density, init_stretched_length]
  data:
    - [lower, ground, mid, 2, test_mat, 100.0]
    - [upper, mid, top, 2, test_mat, 100.0]

winches:
  headers: [name, tether_idxs, winch_point]
  data:
    - [winch_a, [lower], ground]
"""
        yaml_path = joinpath(tmpdir, "stacked.yaml")
        write(yaml_path, stacked_yaml)
        @test_throws ErrorException load_sys_struct_from_yaml(
            yaml_path; system_name="init_stretched_length_nonroot", set=set)
    end

    # ================================================================
    # Test 8: init_tether_force / init_stretch_frac derive len
    # ================================================================
    @testset "init force and stretch_frac" begin
        yaml_path = joinpath(tmpdir, "r2_force.yaml")
        write(yaml_path, INIT_LEN_YAML_ROUTE2)
        sys = load_sys_struct_from_yaml(yaml_path;
            system_name="init_stretched_length_r2_force", set=set)
        SymbolicAWEModels.reinit!(sys, set)

        tether = sys.tethers[:main_tether]
        segs = sys.segments
        stretched = sum(segs[si].len for si in tether.segment_idxs)
        k = SymbolicAWEModels.tether_unit_stiffness(tether, segs)
        @test stretched ≈ 200.0

        # default: force 0, no frac → len == stretched
        @test tether.init_tether_force == 0.0
        @test isnothing(tether.init_stretch_frac)
        SymbolicAWEModels.apply_tether_init_forces!(sys)
        @test tether.len ≈ stretched

        # stretch_frac → len = frac * stretched; clears the force
        tether.init_stretch_frac = 0.8
        @test isnothing(tether.init_tether_force)
        SymbolicAWEModels.apply_tether_init_forces!(sys)
        @test tether.len ≈ 0.8 * stretched

        # stretch_frac > 1 → slack: len longer than stretched
        tether.init_stretch_frac = 1.1
        SymbolicAWEModels.apply_tether_init_forces!(sys)
        @test tether.len ≈ 1.1 * stretched

        # non-positive stretch_frac errors
        tether.init_stretch_frac = 0.0
        @test_throws ErrorException SymbolicAWEModels.apply_tether_init_forces!(sys)

        # force → len = stretched * (1 - force/k); clears the frac
        tether.init_tether_force = 0.1 * k
        @test isnothing(tether.init_stretch_frac)
        SymbolicAWEModels.apply_tether_init_forces!(sys)
        @test tether.len ≈ stretched * 0.9

        # only one of force / frac may have a value at a time
        @test_throws ErrorException Tether(:both, [:s1];
            tether_force=1.0, stretch_frac=0.9)
    end

    # ================================================================
    # Test 9: placement translates a rigid body anchored via a BODY_STATIC point
    # ================================================================
    @testset "Placement moves a rigid body (BODY_STATIC anchor)" begin
        seg_len = 1.0
        inertia = [0.01, 0.1, 0.1]
        bodies = [
            Body(:root; mass=1.0, inertia_principal=inertia,
                pos=[0.5seg_len, 0.0, 0.0], type=STATIC),
            Body(:tip; mass=1.0, inertia_principal=inertia,
                pos=[1.5seg_len, 0.0, 0.0], type=DYNAMIC),
        ]
        joints = [ElasticJoint(:j, :root, :tip; anchor_a=[seg_len/2, 0, 0],
            anchor_b=[-seg_len/2, 0, 0], stiffness_axial=1e4,
            stiffness_shear=1e4, stiffness_torsion=1e3, stiffness_bending=1e3,
            damping_trans=10.0, damping_rot=5.0)]
        # Ground 5 m below the tip; standoff 6 m forces the tip body up ~1 m.
        points = [
            Point(:ground, [2.0seg_len, 0.0, -5.0], STATIC),
            Point(:tip_anchor, [2.0seg_len, 0.0, 0.0], BODY_STATIC;
                  body=:tip, anchor_b=[seg_len/2, 0.0, 0.0]),
        ]
        segments = [Segment(:tether_seg, :ground, :tip_anchor,
            1e4, 10.0, 0.01; l0=5.0)]
        tethers = [Tether(:tether, [:tether_seg], 6.0)]
        winches = [Winch(:winch, set, [:tether]; winch_point=:ground)]
        sys = SystemStructure("init_stretched_length_body", set; points,
            segments, tethers, winches, bodies=bodies,
            elastic_joints=joints)

        SymbolicAWEModels.reinit!(sys, set)

        tip = sys.bodies[:tip]
        @test tip.pos_w[3] > 0.5
        R_b_to_w = SymbolicAWEModels.quaternion_to_rotation_matrix(tip.Q_b_to_w)
        anchor_w = tip.pos_w .+ R_b_to_w * sys.points[:tip_anchor].anchor_b
        @test sys.points[:tip_anchor].pos_w ≈ anchor_w
        @test norm(sys.points[:tip_anchor].pos_w -
                   sys.points[:ground].pos_w) ≈ 6.0
        # The clamped root body must not move.
        @test sys.bodies[:root].pos_w ≈ KVec3(0.5seg_len, 0, 0)
    end

    # ================================================================
    # Test 10: placement translates a RIGID_DYNAMICS wing whose free end
    # is a WING point. A WING point carries its body in wing_idx (its
    # body_idx is 0), so collecting moved bodies via body_idx alone leaves
    # the wing behind while its points move — the pos~anchor constraint
    # then snaps everything back at simulation time.
    # ================================================================
    @testset "Placement moves a rigid wing (WING-point free end)" begin
        wing_yaml = """
variables:
  test_mat:
    unit_damping: 0.0
    density: 724.0

points:
  headers: [name, pos_cad, type, wing_idx, transform_idx,
            extra_mass, body_frame_damping, world_frame_damping,
            area, drag_coeff, body_idx]
  data:
    - [le_left,   [-0.5, 1.0, 2.0], BODY_STATIC, main_wing, ~, 0.5, 0.0, 0.0, 0.0, 0.0, main_wing]
    - [te_left,   [0.5,  1.0, 2.2], BODY_STATIC, main_wing, ~, 0.5, 0.0, 0.0, 0.0, 0.0, main_wing]
    - [le_center, [-0.5, 0.0, 2.5], BODY_STATIC, main_wing, ~, 0.5, 0.0, 0.0, 0.0, 0.0, main_wing]
    - [te_center, [0.5,  0.0, 2.7], BODY_STATIC, main_wing, ~, 0.5, 0.0, 0.0, 0.0, 0.0, main_wing]
    - [le_right,  [-0.5,-1.0, 2.0], BODY_STATIC, main_wing, ~, 0.5, 0.0, 0.0, 0.0, 0.0, main_wing]
    - [te_right,  [0.5, -1.0, 2.2], BODY_STATIC, main_wing, ~, 0.5, 0.0, 0.0, 0.0, 0.0, main_wing]
    - [ground,    [0.0, 0.0, 0.0], STATIC, ~, ~, 0.0, 0.0, 0.0, 0.0, 0.0, ~]

segments:
  headers: [name, point_i, point_j, l0, diameter_mm,
            unit_stiffness, unit_damping, compression_frac]
  data:
    - [tether_seg, ground, le_center, 0.0, 4.0, 120000.0, 350.0, 0.0]

wings:
  data:
    - name: main_wing
      dynamics_type: RIGID_DYNAMICS
      aero_mode: AERO_NONE
      transform_idx: 0
      y_damping: 0.0

tethers:
  headers: [name, segment_idxs, init_stretched_length]
  data:
    - [main_tether, [tether_seg], 5.0990195135927845]

winches:
  headers: [name, tether_idxs, winch_point]
  data:
    - [main_winch, [main_tether], ground]
"""
        yaml_path = joinpath(tmpdir, "wing.yaml")
        write(yaml_path, wing_yaml)
        sys = load_sys_struct_from_yaml(yaml_path;
            system_name="init_stretched_length_wing", set=set,
            aero_mode=AeroNone())

        wing = sys.bodies[:main_wing]
        pos_cad = copy(wing.pos_cad)
        # Free end le_center at [-0.5, 0, 2.5]; standoff = 2× its distance to
        # ground, so the whole wing translates by delta = [-0.5, 0, 2.5].
        delta = KVec3(-0.5, 0.0, 2.5)

        SymbolicAWEModels.reinit!(sys, set)

        ground = sys.points[:ground].pos_w
        le_center = sys.points[:le_center].pos_w
        # Free end sits at the stretched standoff.
        @test norm(le_center - ground) ≈ 5.0990195135927845
        # The wing body translated with its points (the regression: the body
        # stayed at pos_cad while the points moved).
        @test wing.pos_w ≈ pos_cad .+ delta
        # Body and its WING points stay rigid: point→origin offset preserved.
        @test le_center .- wing.pos_w ≈
              sys.points[:le_center].pos_cad .- pos_cad
    end

end
nothing
