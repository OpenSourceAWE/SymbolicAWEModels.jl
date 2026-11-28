# SPDX-FileCopyrightText: 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_tether_init.jl - Tether initial length scaling tests
#
# Tests the Tether.init_stretched_length feature: scaling pos_w before transforms.
# Uses Route 2 (auto-generated) tethers and YAML-specified init_stretched_length.
# All tests use reinit! directly on a SystemDefinition (no ODE compilation).

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

using Test
using SymbolicAWEModels
using SymbolicAWEModels: KVec3, KiteGeometry
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

# The beam of `tube_rider_beam` as an authoring document.
const TUBE_RIDER_BEAM_YAML = """
bodies:
  headers: [name, extra_mass, inertia_principal, pos, type]
  data:
    - [node1, 1.0, [0.01, 0.01, 0.01], [0.0, 0.0, 0.0], DYNAMIC]
    - [node2, 1.0, [0.01, 0.01, 0.01], [1.0, 0.0, 0.0], DYNAMIC]
    - [node3, 1.0, [0.01, 0.01, 0.01], [2.0, 0.0, 0.0], DYNAMIC]
    - [node4, 1.0, [0.01, 0.01, 0.01], [3.0, 0.0, 0.0], DYNAMIC]

tubes:
  headers: [name, bodies, diameter, pressure, law, model, EA, GA, GJ, EIy, EIz,
            shear_coeff, damping]
  data:
    - [j1, [node1, node2], 0.1, 30000.0, breukels2011, timoshenko, 10000.0,
       1500.0, 50.0, 100.0, 100.0, 0.8333333333333334, 0.05]
    - [j2, [node2, node3], 0.1, 30000.0, breukels2011, timoshenko, 10000.0,
       1500.0, 50.0, 100.0, 100.0, 0.8333333333333334, 0.05]
    - [j3, [node3, node4], 0.1, 30000.0, breukels2011, timoshenko, 10000.0,
       1500.0, 50.0, 100.0, 100.0, 0.8333333333333334, 0.05]

points:
  headers: [name, pos_cad, type, tube, extra_mass]
  data:
    - [ground, [0.5, 0.0, -5.0], STATIC, nothing, 0.0]
    - [beam_anchor, [0.5, 0.0, 0.0], BODY_STATIC, j1, 0.0]
    - [tip_anchor, [1.5, 0.0, 0.0], BODY_STATIC, j2, 0.0]
    - [tail, [1.5, 0.0, -1.0], DYNAMIC, nothing, 1.0]

segments:
  headers: [name, point_i, point_j, l0, diameter_mm, unit_stiffness, unit_damping]
  data:
    - [tether_seg, ground, beam_anchor, 5.0, 10.0, 10000.0, 10.0]
    - [tail_seg, tip_anchor, tail, 1.0, 10.0, 10000.0, 10.0]

tethers:
  headers: [name, segment_idxs, init_stretched_length]
  data:
    - [tether, [tether_seg], 6.0]

winches:
  headers: [name, tether_idxs, winch_point]
  data:
    - [winch, [tether], ground]
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

"""The origin of body `body` (a name or an index) of `sys` in its initial pose."""
body_origin(sys, body) = SymbolicAWEModels.origin_ENU(SymbolicAWEModel(nothing, sys),
    body isa Symbol ? KiteGeometry.row_index(sys.bodies, body) : body)

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
        sys, _ = load_authoring(
            yaml_path; name="init_stretched_length_r2_yaml", set=set)

        # init_stretched_length=200 already set from YAML; no programmatic change needed
        place!(sys)

        mid = sys.points[:main_tether_point_1]
        @test mid.pos_ENU ≈ KVec3(0, 0, -100)
        @test sys.points[:top].pos_ENU ≈ KVec3(0, 0, -200)
        sam = SymbolicAWEModel(set, sys)
        @test run_value(sam, :tethers, :main_tether, :len) ≈ 200.0
        # The stretched placement is the initial pose.
        @test run_value(sam, :points, :main_tether_point_1, :pos_w) == mid.pos_ENU
        @test run_value(sam, :points, :top, :pos_w) == sys.points[:top].pos_ENU
    end

    # ================================================================
    # Test 2: Route 1 explicit segments + YAML init_stretched_length
    # ================================================================
    @testset "Route 1 YAML init_stretched_length: scale to 2x" begin
        yaml_path = joinpath(tmpdir, "r1_yaml.yaml")
        write(yaml_path, INIT_LEN_YAML_ROUTE1)
        sys, _ = load_authoring(
            yaml_path; name="init_stretched_length_r1_yaml", set=set)

        place!(sys)

        @test sys.points[:mid].pos_ENU ≈ KVec3(0, 0, -100)
        @test sys.points[:top].pos_ENU ≈ KVec3(0, 0, -200)
        @test sys.segments[:s1].rest_length ≈ 100.0
        @test sys.segments[:s2].rest_length ≈ 100.0
        @test run_value(SymbolicAWEModel(set, sys), :points, :top, :pos_w) ==
              sys.points[:top].pos_ENU
    end

    # ================================================================
    # Test 3: Route 2 + downstream point translation
    # ================================================================
    @testset "Route 2 downstream point translation" begin
        yaml_path = joinpath(tmpdir, "r2_downstream.yaml")
        write(yaml_path, INIT_LEN_DOWNSTREAM_ROUTE2)
        sys, _ = load_authoring(
            yaml_path; name="init_stretched_length_r2_downstream", set=set)

        place!(sys)

        @test sys.points[:top].pos_ENU ≈ KVec3(0, 0, -200)
        @test sys.points[:downstream].pos_ENU ≈ KVec3(10, 0, -200)
    end

    # ================================================================
    # Test 4: Error — downstream non-tether segment connects back to start
    # ================================================================
    @testset "Error on loop to start" begin
        yaml_path = joinpath(tmpdir, "r2_loop.yaml")
        write(yaml_path, INIT_LEN_LOOP_ROUTE2)
        @test_throws ErrorException load_authoring(
            yaml_path; name="init_stretched_length_r2_loop", set=set)
    end

    # ================================================================
    # Test 5: Idempotency — repeated reinit! gives same result
    # ================================================================
    @testset "Idempotency" begin
        yaml_path = joinpath(tmpdir, "r2_yaml.yaml")
        sys, _ = load_authoring(
            yaml_path; name="init_stretched_length_r2_idem", set=set)

        place!(sys)
        mid_pos = copy(sys.points[:main_tether_point_1].pos_ENU)
        top_pos = copy(sys.points[:top].pos_ENU)

        place!(sys)
        @test sys.points[:main_tether_point_1].pos_ENU ≈ mid_pos
        @test sys.points[:top].pos_ENU ≈ top_pos
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
        sys, _ = load_authoring(
            yaml_path; name="init_stretched_length_multi", set=set)

        sys.tethers[:tether_static].init_stretched_length = 200.0
        place!(sys)

        ground_static = sys.points[:ground_static].pos_ENU
        @test norm(sys.points[:top].pos_ENU - ground_static) ≈ 200.0
        @test sys.points[:ground_static].pos_ENU ≈ KVec3(10, 0, 0)
        @test sys.points[:ground_winch].pos_ENU ≈ KVec3(-10, 0, 0)

        sys.tethers[:tether_winch].init_stretched_length = 100.0
        top_before = copy(sys.points[:top].pos_ENU)
        place!(sys)

        # top moves along the line from the mean anchor (the origin) through where
        # it was, until the two tethers' mean length is the mean target, 150 m.
        placed = copy(sys.points[:top].pos_ENU)
        @test normalize(placed) ≈ normalize(top_before)
        @test (norm(placed - ground_static) + norm(placed - KVec3(-10, 0, 0))) / 2 ≈ 150.0
        @test sys.points[:ground_winch].pos_ENU ≈ KVec3(-10, 0, 0)

        place!(sys)
        @test sys.points[:top].pos_ENU ≈ placed

        # Tethers that already have their lengths stay where they are.
        sys.tethers[:tether_static].init_stretched_length = norm(placed - ground_static)
        sys.tethers[:tether_winch].init_stretched_length = norm(placed - KVec3(-10, 0, 0))
        place!(sys)
        @test sys.points[:top].pos_ENU ≈ placed
    end

    # ================================================================
    # Test 7: init_stretched_length on a non-root tether is an error
    # ================================================================
    @testset "Error on non-root init_stretched_length" begin
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
        @test_throws ErrorException load_authoring(
            yaml_path; name="init_stretched_length_nonroot", set=set)
    end

    # ================================================================
    # Test 8: init_tether_force / init_stretch_frac derive len
    # ================================================================
    @testset "init force and stretch_frac" begin
        yaml_path = joinpath(tmpdir, "r2_force.yaml")
        write(yaml_path, INIT_LEN_YAML_ROUTE2)
        sys, _ = load_authoring(yaml_path;
            name="init_stretched_length_r2_force", set=set)
        place!(sys)

        sam = SymbolicAWEModel(set, sys)
        segment_idxs = sys.tethers[:main_tether].segments
        stretched = sum(run_value(sam, :segments, si, :len) for si in segment_idxs)
        k = run_value(sam, :segments, first(segment_idxs), :unit_stiffness)
        tether = sys.tethers[:main_tether]
        unstretched() =
            (place!(sys); sum(sys.segments[si].rest_length for si in segment_idxs))
        @test stretched ≈ 200.0

        # default: no force, no frac → len == stretched
        @test isnothing(tether.init_tether_force)
        @test isnothing(tether.init_stretch_frac)
        @test unstretched() ≈ stretched

        # stretch_frac → len = frac * stretched
        tether.init_stretch_frac = 0.8
        @test unstretched() ≈ 0.8 * stretched

        # stretch_frac > 1 → slack: len longer than stretched
        tether.init_stretch_frac = 1.1
        @test unstretched() ≈ 1.1 * stretched

        # non-positive stretch_frac errors
        tether.init_stretch_frac = 0.0
        @test_throws ErrorException place!(sys)

        # force and frac both given: placing refuses the tether
        tether.init_tether_force = 0.1 * k
        @test_throws ErrorException place!(sys)

        # force → len = stretched * (1 - force/k)
        tether.init_stretch_frac = nothing
        @test unstretched() ≈ stretched * 0.9
    end

    # ================================================================
    # Test 9: placement translates a rigid body anchored via a BODY_STATIC point
    # ================================================================
    @testset "Placement moves a rigid body (BODY_STATIC anchor)" begin
        seg_len = 1.0
        inertia = [0.01, 0.1, 0.1]
        materials = KiteGeometry.MaterialTable(KiteGeometry.Material[],
                                               Dict{String, Any}())
        push!(materials.rows, Material(; name=:breukels2011, law="breukels2011",
            parameters=KiteGeometry.tube_law_parameters(:breukels2011)))
        bodies = [
            Body(; name=:root, pos_ENU=[0.5seg_len, 0.0, 0.0], extra_mass=1.0,
                extra_inertia_KA=Matrix(Diagonal(inertia)), fixed=true),
            Body(; name=:tip, pos_ENU=[1.5seg_len, 0.0, 0.0], extra_mass=1.0,
                extra_inertia_KA=Matrix(Diagonal(inertia))),
        ]
        tubes = [Tube(; name=:j, bodies=(:root, :tip), diameter=0.1, pressure=3e4,
            material=:breukels2011,
            model=ElasticTubeModel(; stiffness_axial=1e4, stiffness_shear=1e4,
                                   stiffness_torsion=1e3, stiffness_bending=1e3,
                                   damping=0.05, anchor_a=[seg_len/2, 0.0, 0.0],
                                   anchor_b=[-seg_len/2, 0.0, 0.0]))]
        # Ground 5 m below the tip; standoff 6 m forces the tip body up ~1 m.
        points = [
            Point(; name=:ground, pos_ENU=[2.0seg_len, 0.0, -5.0], fixed=true),
            Point(; name=:tip_anchor, pos_ENU=[2.0seg_len, 0.0, 0.0], body=:tip,
                  model=PointModel(; anchor_b=[seg_len/2, 0.0, 0.0])),
        ]
        segments = [KiteGeometry.line_segment(:tether_seg, (:ground, :tip_anchor), set,
            materials; rest_length=5.0, diameter=0.01, unit_stiffness=1e4,
            model=SegmentSpring(; unit_damping=10.0))]
        tethers = [Tether(; name=:tether, segments=[:tether_seg],
                          init_stretched_length=6.0)]
        winches = [Winch(; name=:winch, tethers=[:tether], winch_point=:ground)]
        sys = system_definition("init_stretched_length_body", set; points,
            segments, materials, tethers, winches, bodies=bodies, tubes)

        place!(sys)

        tip = sys.bodies[:tip]
        @test body_origin(sys, :tip)[3] > 0.5
        R_b_to_w = SymbolicAWEModels.quaternion_to_rotation_matrix(tip.Q_KA2ENU)
        anchor_w = body_origin(sys, :tip) .+
            R_b_to_w * sys.points[:tip_anchor].anchor_b
        @test sys.points[:tip_anchor].pos_ENU ≈ anchor_w
        @test norm(sys.points[:tip_anchor].pos_ENU -
                   sys.points[:ground].pos_ENU) ≈ 6.0
        # The clamped root body must not move.
        @test body_origin(sys, :root) ≈ KVec3(0.5seg_len, 0, 0)
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
      angular_damping: [0.0, 0.0, 0.0]

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
        sys, _ = load_authoring(yaml_path;
            name="init_stretched_length_wing", set=set)

        # Free end le_center at [-0.5, 0, 2.5]; standoff = 2× its distance to
        # ground, so the whole wing translates by delta = [-0.5, 0, 2.5]. The wing
        # origin is the mean of its six equal-mass points as authored.
        delta = KVec3(-0.5, 0.0, 2.5)
        authored_origin = KVec3(0.0, 0.0, 13.6 / 6)

        place!(sys)

        ground = sys.points[:ground].pos_ENU
        le_center = sys.points[:le_center].pos_ENU
        # Free end sits at the stretched standoff.
        @test norm(le_center - ground) ≈ 5.0990195135927845
        # The wing body translated with its points (the regression: the body
        # stayed where it was authored while the points moved).
        @test body_origin(sys, :main_wing) ≈ authored_origin .+ delta
        # Body and its WING points stay rigid: point→origin offset preserved.
        @test le_center .- body_origin(sys, :main_wing) ≈
              KVec3(-0.5, 0.0, 2.5) .- authored_origin
    end

    # ================================================================
    # Test 11: placement translates a free-floating beam anchored only through
    # tube-riding BODY_STATIC points. Such a point has body_idx 0 and carries
    # its association in tube_idx, so collecting moved bodies via body_idx
    # alone leaves the whole beam behind while the tether's free end moves —
    # the ride constraint then snaps the point back, silently.
    # ================================================================
    """A beam of four bodies on three tubes, hung from the ground by a tether to a point
    riding the first tube, with a tail off a point riding the second."""
    function tube_rider_beam(name)
        inertia = [0.01, 0.01, 0.01]
        materials = KiteGeometry.MaterialTable(KiteGeometry.Material[],
                                               Dict{String, Any}())
        push!(materials.rows, Material(; name=:breukels2011, law="breukels2011",
            parameters=KiteGeometry.tube_law_parameters(:breukels2011)))
        bodies = [Body(; name=Symbol(:node, i), pos_ENU=[Float64(i - 1), 0.0, 0.0],
                       extra_mass=1.0, extra_inertia_KA=Matrix(Diagonal(inertia)))
                  for i in 1:4]
        tubes = [Tube(; name=Symbol(:j, i), bodies=(Symbol(:node, i), Symbol(:node, i + 1)),
                      diameter=0.1, pressure=3e4, material=:breukels2011,
                      model=TimoshenkoTubeModel(; EA=1.0e4, GA=1500.0, GJ=50.0,
                                                EIy=100.0, EIz=100.0, shear_coeff=5/6,
                                                damping=0.05)) for i in 1:3]
        points = [
            Point(; name=:ground, pos_ENU=[0.5, 0.0, -5.0], fixed=true),
            Point(; name=:beam_anchor, pos_ENU=[0.5, 0.0, 0.0], tube=:j1),
            Point(; name=:tip_anchor, pos_ENU=[1.5, 0.0, 0.0], tube=:j2),
            Point(; name=:tail, pos_ENU=[1.5, 0.0, -1.0], extra_mass=1.0),
        ]
        line(name, ends, rest_length) = KiteGeometry.line_segment(name, ends, set,
            materials; rest_length, diameter=0.01, unit_stiffness=1e4,
            model=SegmentSpring(; unit_damping=10.0))
        segments = [line(:tether_seg, (:ground, :beam_anchor), 5.0),
                    line(:tail_seg, (:tip_anchor, :tail), 1.0)]
        tethers = [Tether(; name=:tether, segments=[:tether_seg],
                          init_stretched_length=6.0)]
        winches = [Winch(; name=:winch, tethers=[:tether], winch_point=:ground)]
        return system_definition(name, set; points, segments, materials, tethers,
                               winches, bodies, tubes)
    end

    @testset "Placement moves a beam (tube-anchored points)" begin
        sys = tube_rider_beam("init_stretched_length_beam")

        delta = KVec3(0.0, 0.0, 1.0)
        place!(sys)

        @test norm(sys.points[:beam_anchor].pos_ENU -
                   sys.points[:ground].pos_ENU) ≈ 6.0
        # Every beam body translates, including node4, which no point rides and
        # which is reachable only through the joint graph.
        for i in eachindex(sys.bodies)
            @test body_origin(sys, i) ≈
                  KVec3(Float64(i - 1), 0.0, 0.0) .+ delta
        end
        # tail hangs off the far end of the beam, reachable from the tether's
        # free end only through the rigid-sibling hop across the beam.
        @test sys.points[:tail].pos_ENU ≈ KVec3(1.5, 0.0, -1.0) .+ delta
        @test KiteGeometry.segment_length(sys.segments[:tail_seg],
                                          KiteGeometry.Positions(sys.points)) ≈ 1.0
    end

    # ================================================================
    # Test 11b: the beam of Test 11 lands on one pose and one set of rest lengths
    # whether it is loaded from YAML, built in Julia, or placed again after its
    # transform changes.
    # ================================================================
    @testset "YAML, Julia and re-placement agree (tube riders)" begin
        same_pose(a, b) = all(
            [all(p.pos_ENU ≈ q.pos_ENU for (p, q) in zip(a.points, b.points)),
             all(body_origin(a, i) ≈ body_origin(b, i) for i in eachindex(a.bodies)),
             all(p.Q_KA2ENU ≈ q.Q_KA2ENU for (p, q) in zip(a.bodies, b.bodies)),
             all(p.rest_length ≈ q.rest_length
                 for (p, q) in zip(a.segments, b.segments))])

        yaml_path = joinpath(tmpdir, "tube_rider_beam.yaml")
        write(yaml_path, TUBE_RIDER_BEAM_YAML)
        from_yaml, _ = load_authoring(
            yaml_path; name="tube_rider_beam_yaml", set=set)
        from_julia = tube_rider_beam("tube_rider_beam_julia")
        @test same_pose(from_yaml, from_julia)

        turn(elevation) = Transform(; name=:tf, elevation, azimuth=0.0, heading=0.0,
            base_pos=[0.5, 0.0, -5.0], base_point=1, rot_point=[2 => 1.0],
            points=[2, 3, 4], bodies=[1, 2, 3, 4])
        replaced = tube_rider_beam("tube_rider_beam_replaced")
        place!(replaced, [turn(deg2rad(30))])
        place!(replaced, [turn(deg2rad(60))])
        placed_once = tube_rider_beam("tube_rider_beam_once")
        place!(placed_once, [turn(deg2rad(60))])
        @test same_pose(replaced, placed_once)
        @test !same_pose(replaced, from_julia)
    end

    # ================================================================
    # Test 12: set_unstretched_length! writes len and l0 and moves nothing
    # ================================================================
    @testset "set_unstretched_length! sets len without re-placing" begin
        yaml_path = joinpath(tmpdir, "r1_yaml.yaml")
        sys, _ = load_authoring(
            yaml_path; name="init_stretched_length_r1_yaml", set=set)
        place!(sys)

        sam = SymbolicAWEModel(set, sys)
        tether_len() = run_value(sam, :tethers, :main_tether, :len)
        point_pos_w() = [run_value(sam, :points, i, :pos_w) for i in eachindex(sys.points)]
        @test tether_len() ≈ 200.0
        placed_w = point_pos_w()
        placed_initial = [copy(point.pos_ENU) for point in sys.points]

        set_unstretched_length!(sam, :main_tether, 150.0)

        @test tether_len() ≈ 150.0
        @test run_value(sam, :segments, :s1, :rest_length) ≈ 75.0
        @test run_value(sam, :segments, :s2, :rest_length) ≈ 75.0
        @test all(pos ≈ was for (pos, was) in zip(point_pos_w(), placed_w))
        @test all(point.pos_ENU ≈ was
                  for (point, was) in zip(sys.points, placed_initial))
        # The stretched geometry is untouched, so the rope is now taut.
        @test sum(run_value(sam, :segments, idx, :len)
                  for idx in sys.tethers[:main_tether].segments) ≈ 200.0

        # `reinit!` re-places the structure, so it derives `len` afresh.
        place!(sys)
        @test run_value(SymbolicAWEModel(set, sys), :tethers, :main_tether, :len) ≈ 200.0

        @test_throws ErrorException set_unstretched_length!(sam, :main_tether, 0.0)
        @test_throws ErrorException set_unstretched_length!(sam, :main_tether, -1.0)
    end

    @testset "a tether without a winch builds on the kernel backend" begin
        materials = KiteGeometry.MaterialTable(KiteGeometry.Material[],
                                               Dict{String, Any}())
        points = [Point(; name=:ground, pos_ENU=[0.0, 0.0, 0.0], fixed=true),
                  Point(; name=:middle, pos_ENU=[0.0, 0.0, -5.0], extra_mass=0.5),
                  Point(; name=:mass, pos_ENU=[0.0, 0.0, -10.0], extra_mass=1.0)]
        line(name, ends) = KiteGeometry.line_segment(name, ends, set, materials;
            rest_length=5.0, diameter=0.01, unit_stiffness=1e4,
            model=SegmentSpring(; unit_damping=10.0))
        segments = [line(:s1, (:ground, :middle)), line(:s2, (:middle, :mass))]
        tethers = [Tether(; name=:line, segments=[:s1, :s2])]
        sys = system_definition("unwinched_tether", set; points, segments, materials,
                                tethers)
        sam = SymbolicAWEModel(set, sys; backend=KernelBackend())
        init!(sam; prn=false)
        next_step!(sam; dt=0.01, vsm_interval=0)
        @test run_value(sam, :segments, :s1, :rest_length) ≈ 5.0
        @test run_value(sam, :tethers, :line, :len) ≈ 10.0
    end

end
nothing
