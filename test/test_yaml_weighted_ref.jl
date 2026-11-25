# SPDX-FileCopyrightText: 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_yaml_weighted_ref.jl - Test weighted z/y_ref_points
# parsing from YAML
#
# Verifies that nested weighted ref point specs like
#   z_ref_points: [1, [[12, 0.7], [11, 0.3]]]
# are correctly parsed into a wing frame's weighted reference points.

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

using Test
using LinearAlgebra
using SymbolicAWEModels
using SymbolicAWEModels: KVec3, VortexStepMethod
using KiteUtils

# Minimal YAML with a PARTICLE_DYNAMICS wing using weighted z_ref_points
# 3 LE/TE pairs to match the 3 aero sections in
# 2plate_kite/vsm_settings.yaml.
# Layout mirrors particle_structural_geometry.yaml.
const WEIGHTED_REF_YAML = """
points:
  headers: [idx, pos_cad, type, wing_idx, transform_idx,
            extra_mass, body_frame_damping,
            world_frame_damping, area, drag_coeff]
  data:
    - [1, [-0.5, 1.0, 2.0], DYNAMIC, 1, 1,
       0.1, 10.0, 0.0, 0.0, 0.0]
    - [2, [0.5, 1.0, 2.3], DYNAMIC, 1, 1,
       0.1, 10.0, 0.0, 0.0, 0.0]
    - [3, [-0.5, 0.0, 2.5], DYNAMIC, 1, 1,
       0.1, 10.0, 0.0, 0.0, 0.0]
    - [4, [0.5, 0.0, 2.8], DYNAMIC, 1, 1,
       0.1, 10.0, 0.0, 0.0, 0.0]
    - [5, [-0.5, -1.0, 2.0], DYNAMIC, 1, 1,
       0.1, 10.0, 0.0, 0.0, 0.0]
    - [6, [0.5, -1.0, 2.3], DYNAMIC, 1, 1,
       0.1, 10.0, 0.0, 0.0, 0.0]
    - [7, [0.0, 0.0, 0.0], DYNAMIC, 1, 1,
       1.0, 0.0, 0.0, 0.1, 1.0]
    - [8, [0.0, 0.0, -20.0], STATIC, 1, 1,
       0.0, 0.0, 0.0, 0.0, 0.0]
    - [9, [1.0, 0.0, 0.0], DYNAMIC, 1, 1,
       1.0, 0.0, 0.0, 0.1, 1.0]

segments:
  headers: [idx, point_i, point_j, l0, diameter_mm,
            unit_stiffness, unit_damping, compression_frac]
  data:
    - [1, 1, 2, 0, 1.0, 5000.0, 10.0, 1.0]
    - [2, 3, 4, 0, 1.0, 5000.0, 10.0, 1.0]
    - [3, 5, 6, 0, 1.0, 5000.0, 10.0, 1.0]
    - [4, 1, 7, 0, 1.0, 5000.0, 10.0, 0.01]
    - [5, 3, 7, 0, 1.0, 5000.0, 10.0, 0.01]
    - [6, 5, 7, 0, 1.0, 5000.0, 10.0, 0.01]
    - [7, 7, 9, 0, 1.0, 5000.0, 10.0, 0.01]

stations:
  headers: [name, point_idxs, type, moment_frac, damping]
  data:
    - [ts_left,   [1, 2], STATIC, 0.75, 50.0]
    - [ts_center, [3, 4], STATIC, 0.75, 50.0]
    - [ts_right,  [5, 6], STATIC, 0.75, 50.0]

wings:
  data:
    - idx: 1
      dynamics_type: PARTICLE_DYNAMICS
      aero_mode: AERO_NONE
      stations: [ts_left, ts_center, ts_right]
      origin_idx: [[7, 0.7], [9, 0.3]]
      z_ref_points: [7, [[3, 0.7], [5, 0.3]]]
      y_ref_points: [1, 5]

transforms:
  data:
    - idx: 1
      elevation: 50
      azimuth: 0.0
      heading: 0.0
      wing_idx: 1
      base_pos: [0.0, 0.0, 0.0]
      base_point_idx: 8
"""

@testset "Weighted z/y_ref_points from YAML" begin
    pkg_root = dirname(@__DIR__)
    src_data = joinpath(pkg_root, "data", "2plate_kite")

    tmpdir = mktempdir()
    data_path = joinpath(tmpdir, "2plate_kite")
    cp(src_data, data_path; force=true)

    yaml_path = joinpath(tmpdir, "geometry.yaml")
    write(yaml_path, WEIGHTED_REF_YAML)

    set_data_path(data_path)
    set = Settings("system.yaml")

    vsm_set = VortexStepMethod.VSMSettings(
        joinpath(data_path, "vsm_settings.yaml");
        data_prefix=false)

    # This should not throw — the bug causes
    # MethodError: no method matching Int64(::Vector{Real})
    sys, transforms = load_authoring(yaml_path;
        name="weighted_ref_test", set)

    wing = SymbolicAWEModels.SystemState(sys; set, vsm_set).wings[1]

    # z_ref_points[1] should be a single point (idx 7, kcu)
    z_p1, z_p2 = wing.wing_frame.z_ref_points
    @test first.(z_p1) == [7]
    @test last.(z_p1) == [1.0]

    # z_ref_points[2] should be weighted: point 3 @ 0.7,
    # point 5 @ 0.3
    @test first.(z_p2) == [3, 5]
    @test last.(z_p2) ≈ [0.7, 0.3]

    # y_ref_points should be simple single-point refs
    y_p1, y_p2 = wing.wing_frame.y_ref_points
    @test first.(y_p1) == [1]
    @test first.(y_p2) == [5]

    # Weighted origin_idx: point 7 @ 0.7, point 9 @ 0.3
    @test wing.origin == wing.wing_frame.origin
    @test first.(wing.origin) == [7, 9]
    @test last.(wing.origin) ≈ [0.7, 0.3]

    # The wing origin is the weighted centroid of points 7 and 9
    @test SymbolicAWEModels.origin_ENU(wing) ≈
          0.7 .* sys.points[7].pos_ENU .+ 0.3 .* sys.points[9].pos_ENU

    # N-point (4-point) weighted origin should also parse
    # and resolve correctly
    yaml_4pt = replace(WEIGHTED_REF_YAML,
        "origin_idx: [[7, 0.7], [9, 0.3]]" =>
        "origin_idx: [[1, 0.1], [3, 0.2], [5, 0.3], [9, 0.4]]")
    yaml_path_4pt = joinpath(tmpdir, "geometry_4pt.yaml")
    write(yaml_path_4pt, yaml_4pt)
    sys_4pt, _ = load_authoring(yaml_path_4pt;
        name="weighted_ref_4pt_test", set)
    wing_4pt = SymbolicAWEModels.SystemState(sys_4pt; set, vsm_set).wings[1]
    @test first.(wing_4pt.origin) == [1, 3, 5, 9]
    @test last.(wing_4pt.origin) ≈ [0.1, 0.2, 0.3, 0.4]

    # Compile and init — weighted refs must survive
    # symbolic equation generation and scalarization
    sam = SymbolicAWEModel(set, sys; transforms, vsm_set)
    integ = init!(sam; prn=false)
    @test integ !== nothing

    # The fitted pose read back after a step must use the
    # weights, not just the first point of each reference.
    next_step!(sam)
    points = sam.state.points
    wing = sam.state.bodies[1]
    blend(ref) = sum(w * points[i].pos_w for (i, w) in ref)
    @test wing.pos_w ≈ blend(wing.origin)
    @test wing.vel_w ≈ sum(w * points[i].vel_w for (i, w) in wing.origin)
    z_axis = normalize(blend(wing.wing_frame.z_ref_points[2]) -
                       blend(wing.wing_frame.z_ref_points[1]))
    @test wing.R_b_to_w[:, 3] ≈ z_axis
end
nothing
