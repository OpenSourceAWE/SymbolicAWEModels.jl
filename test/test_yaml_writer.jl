# SPDX-FileCopyrightText: 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_yaml_writer.jl - save_sys_struct_to_yaml has to reproduce a structure well
# enough to build the same ODE from it. A structure is loaded, its parameters are
# perturbed away from the source file so that agreement cannot come from both
# sides re-reading that file, then it is saved, reloaded and rebuilt. The
# reloaded model must match in both the state vector `u` and the right-hand side
# `f(u, p, t)`: `u` alone would pass even if every parameter were dropped, since
# stiffnesses, rest lengths and masses live in `p`.
#
# KNOWN GAP: a rigid wing placed by a `transforms` block does not reload yet —
# `reinit!` reports "Wing/rot position and base position overlap", i.e. the wing
# body is not being carried by the transform the way it is when the same data is
# read from a hand-written file. The bodies/joints structure below has no
# transform and does round-trip. Extend this test to the 2plate rigid geometry
# once that is fixed.

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
    log_file: "data/yaml_writer_test"
    g_earth: 9.81
solver:
    solver: "FBDF"
    abs_tol: 1.0e-8
    rel_tol: 1.0e-8
kite:
    model: ""
    foil_file: "ram_air_kite/ram_air_kite_foil.dat"
    physical_model: "yaml_writer_test"
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
    - [nodeC, 1.0, [0.01, 0.01, 0.01], [2.0, 0.0, 0.0], DYNAMIC]

timoshenko_joints:
  headers: [name, body_a, body_b, EA, GA, GJ, EIy, EIz, shear_coeff,
            damping_trans, damping_rot]
  data:
    - [joint_ab, nodeA, nodeB, 10000.0, 1500.0, 50.0, 100.0, 100.0, 0.8333,
       200.0, 3.0]
    - [joint_bc, nodeB, nodeC, 10000.0, 1500.0, 50.0, 100.0, 100.0, 0.8333,
       200.0, 3.0]

points:
  headers: [name, pos_cad, type, body_idx, anchor_b, extra_mass]
  data:
    - [tip_anchor, [2.0, 0.0, 0.0], BODY_STATIC, nodeC, [0.0, 0.0, 0.0], 0.0]
    - [free_a, [2.0, 0.0, -1.0], DYNAMIC, nothing, nothing, 0.5]
    - [free_b, [2.0, 0.0, -2.0], DYNAMIC, nothing, nothing, 0.5]

segments:
  headers: [name, point_i, point_j, l0, diameter_mm, unit_stiffness,
            unit_damping, compression_frac, density]
  data:
    - [seg_1, tip_anchor, free_a, 0.9, 2.0, 1000.0, 10.0, 0.1, 724.0]
    - [seg_2, free_a, free_b, 0.9, 2.0, 1000.0, 10.0, 0.1, 724.0]
"""

"""Build a model on `sys`, initialise it and return `(u, du)`."""
function rhs_state(set, sys)
    sam = SymbolicAWEModel(set, sys)
    init!(sam; prn = false)
    integ = sam.integrator
    u = copy(integ.u)
    du = similar(u)
    integ.f(du, u, integ.p, integ.t)
    return u, du
end

@testset "save_sys_struct_to_yaml round trip" begin
    pkg_root = dirname(@__DIR__)
    tmpdir = mktempdir()
    data_path = joinpath(tmpdir, "2plate_kite")
    cp(joinpath(pkg_root, "data", "2plate_kite"), data_path; force = true)
    write(joinpath(data_path, "settings.yaml"), SETTINGS_YAML)
    write(joinpath(data_path, "system.yaml"),
        "system:\n  sim_settings: settings.yaml\n")
    source_yaml = joinpath(data_path, "yaml_writer_test.yaml")
    write(source_yaml, STRUC_YAML)
    set_data_path(data_path)
    set = Settings("system.yaml")

    sys = load_sys_struct_from_yaml(source_yaml;
        system_name = "yaml_writer_test", set, prn = false)

    for (k, segment) in enumerate(sys.segments)
        segment.unit_stiffness *= 1 + 0.01 * k
        segment.unit_damping *= 1.05
        segment.l0 *= 0.99
    end
    for point in sys.points
        point.extra_mass += 0.05
    end
    for joint in sys.timoshenko_joints
        joint.EA *= 1.1
        joint.damping_trans *= 1.2
    end

    saved_yaml = joinpath(tmpdir, "round_trip.yaml")
    save_sys_struct_to_yaml(sys, saved_yaml; prn = false)
    @test isfile(saved_yaml)

    reloaded = load_sys_struct_from_yaml(saved_yaml;
        system_name = "yaml_writer_test", set, prn = false)

    @testset "components survive the write" begin
        for field in (:points, :segments, :bodies, :timoshenko_joints)
            @test length(getproperty(reloaded, field)) ==
                  length(getproperty(sys, field))
        end
        @test all(a.unit_stiffness ≈ b.unit_stiffness
                  for (a, b) in zip(sys.segments, reloaded.segments))
        @test all(a.l0 ≈ b.l0
                  for (a, b) in zip(sys.segments, reloaded.segments))
        @test all(a.extra_mass ≈ b.extra_mass
                  for (a, b) in zip(sys.points, reloaded.points))
        @test all(a.EA ≈ b.EA for (a, b) in
                  zip(sys.timoshenko_joints, reloaded.timoshenko_joints))
        @test all(a.damping_trans ≈ b.damping_trans for (a, b) in
                  zip(sys.timoshenko_joints, reloaded.timoshenko_joints))
        @test all(maximum(abs, a.pos_cad .- b.pos_cad) < 1e-12
                  for (a, b) in zip(sys.points, reloaded.points))
    end

    u_before, du_before = rhs_state(set, sys)
    u_after, du_after = rhs_state(set, reloaded)

    @testset "state vector" begin
        @test length(u_after) == length(u_before)
        @test u_after ≈ u_before rtol = 1e-10
    end

    @testset "right-hand side" begin
        @test length(du_after) == length(du_before)
        scale = max(maximum(abs, du_before), 1.0)
        @test maximum(abs, du_after .- du_before) / scale < 1e-8
    end

    rm(tmpdir; recursive = true)
end
nothing
