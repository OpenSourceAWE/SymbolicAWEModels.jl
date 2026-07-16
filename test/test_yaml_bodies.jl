# SPDX-FileCopyrightText: 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_yaml_bodies.jl - Validate loading rigid bodies, a Timoshenko beam joint
# and a body-anchored point from a struc_geometry YAML (the `bodies`,
# `timoshenko_joints` blocks and the point `body_idx`/`anchor_b` fields added to
# load_sys_struct_from_yaml). A clamped root body + free tip body joined by one
# Timoshenko element must load with references resolved and take a successful
# next_step!, matching the programmatic path in test_timoshenko_joint.jl.

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

@isdefined(test_init!) || include(joinpath(@__DIR__, "util.jl"))

using Test
using SymbolicAWEModels
using KiteUtils
using LinearAlgebra

SETTINGS_YAML = """
system:
    log_file: "data/yaml_bodies_test"
    g_earth: 9.81
solver:
    solver: "FBDF"
    abs_tol: 1.0e-8
    rel_tol: 1.0e-8
kite:
    model: ""
    foil_file: "ram_air_kite/ram_air_kite_foil.dat"
    physical_model: "yaml_bodies_test"
    mass: 0.0
tether:
    cd_tether: 0.958
    unit_damping: 0.0
    unit_stiffness: 0.0
    rho_tether: 724.0
    e_tether: 5.5e10
winch:
    winch_model: "TorqueControlledMachine"
    drum_radius: 0.110
    gear_ratio: 1.0
    inertia_total: 0.024
    f_coulomb: 122.0
    c_vf: 30.6
environment:
    rho_0: 1.225
    v_wind: 0.0
    upwind_dir: -90.0
    upwind_elevation: 0.0
    wind_vec: [0.0, 0.0, 0.0]
    profile_law: 0
"""

STRUC_YAML = """
bodies:
  headers: [name, mass, inertia_principal, pos, type]
  data:
    - [nodeA, 1.0, [0.01, 0.01, 0.01], [0.0, 0.0, 0.0], STATIC]
    - [nodeB, 1.0, [0.01, 0.01, 0.01], [1.0, 0.0, 0.0], DYNAMIC]

timoshenko_joints:
  headers: [name, body_a, body_b, EA, GA, GJ, EIy, EIz, shear_coeff,
            damping_trans, damping_rot]
  data:
    - [joint, nodeA, nodeB, 10000.0, 1500.0, 50.0, 100.0, 100.0, 0.8333,
       200.0, 3.0]

points:
  headers: [name, pos_cad, type, body_idx, anchor_b]
  data:
    - [tip_anchor, [1.0, 0.0, 0.0], BODY_STATIC, nodeB, [0.0, 0.0, 0.0]]
"""

@testset "YAML body + Timoshenko joint loading" begin
    pkg_root = dirname(@__DIR__)
    tmpdir = mktempdir()
    data_path = joinpath(tmpdir, "2plate_kite")
    cp(joinpath(pkg_root, "data", "2plate_kite"), data_path; force=true)
    write(joinpath(data_path, "settings.yaml"), SETTINGS_YAML)
    write(joinpath(data_path, "system.yaml"),
        "system:\n  sim_settings: settings.yaml\n")
    struc_yaml = joinpath(data_path, "yaml_bodies_test.yaml")
    write(struc_yaml, STRUC_YAML)
    set_data_path(data_path)
    set = Settings("system.yaml")

    sys = load_sys_struct_from_yaml(struc_yaml;
        system_name="yaml_bodies_test", set)

    @testset "Structure wiring from YAML" begin
        @test length(sys.bodies) == 2
        @test sys.bodies[:nodeA].type == STATIC
        @test sys.bodies[:nodeB].type == DYNAMIC
        @test sys.bodies[:nodeA].mass == 1.0
        @test length(sys.timoshenko_joints) == 1
        joint = sys.timoshenko_joints[:joint]
        @test joint.body_a_idx == sys.bodies[:nodeA].idx
        @test joint.body_b_idx == sys.bodies[:nodeB].idx
        @test joint.rest_length ≈ 1.0          # taken from placed geometry
        @test joint.EA ≈ 10000.0
        anchor = sys.points[:tip_anchor]
        @test anchor.type == BODY_STATIC
        @test anchor.body_idx == sys.bodies[:nodeB].idx
    end

    @testset "next_step! runs on the loaded model" begin
        sam = SymbolicAWEModel(set, sys)
        test_init!(sam; prn=false)
        next_step!(sam; dt=0.001, vsm_interval=0)
        tip = sam.sys_struct.bodies[:nodeB].pos_w
        @test all(isfinite, tip)
        @test tip[1] ≈ 1.0 atol=1e-3           # axial position essentially held
    end

    rm(tmpdir; recursive=true)
end
nothing
