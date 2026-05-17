# SPDX-FileCopyrightText: 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_deserialization.jl - Deserialization and settings sync tests
#
# Verifies that after loading a cached .bin model:
# 1. sam.set === sam.sys_struct.set (reference identity)
# 2. Modified settings propagate through deserialized models

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

using Test
using SymbolicAWEModels
using SymbolicAWEModels: KVec3, MVec3
using KiteUtils
using LinearAlgebra

const DESER_YAML = """
points:
  headers: [name, pos_cad, type, wing_idx, transform_idx,
            extra_mass, body_frame_damping,
            world_frame_damping, area, drag_coeff]
  data:
    - [test_point, [0.0, 0.0, 100.0], DYNAMIC, nothing,
       nothing, 1.0, 0.0, 0.0, 0.5, 1.0]
"""

@testset "Deserialization Tests" begin
    tmpdir = mktempdir()

    yaml_path = joinpath(tmpdir, "deser_geometry.yaml")
    write(yaml_path, DESER_YAML)

    settings_yaml = """
system:
    log_file: "data/2plate"
    g_earth: 9.81
solver:
    solver: "FBDF"
    abs_tol: 0.01
    rel_tol: 0.01
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
    settings_path = joinpath(tmpdir, "settings.yaml")
    write(settings_path, settings_yaml)

    system_yaml = """
system:
  sim_settings: settings.yaml
"""
    system_path = joinpath(tmpdir, "system.yaml")
    write(system_path, system_yaml)

    set_data_path(tmpdir)

    @testset "sam.set delegates to sys_struct.set" begin
        set = Settings("system.yaml")
        sys = load_sys_struct_from_yaml(
            yaml_path; system_name="deser_test_1", set)

        sam = SymbolicAWEModel(set, sys)
        init!(sam; remake=true, prn=false)

        # sam.set must be the exact same object as
        # sam.sys_struct.set
        @test sam.set === sam.sys_struct.set

        # Mutating via sam.set should be visible via
        # sys_struct.set
        sam.set.g_earth = 3.14
        @test sam.sys_struct.set.g_earth == 3.14

        # Assigning sam.set = ... should error
        @test_throws ErrorException (sam.set = set)
    end

    @testset "Deserialized model uses current settings" begin
        # Build and cache a model with zero wind
        set1 = Settings("system.yaml")
        set1.g_earth = 0.0
        set1.v_wind = 0.0
        sys1 = load_sys_struct_from_yaml(
            yaml_path; system_name="deser_test_2", set=set1)
        sam1 = SymbolicAWEModel(set1, sys1)
        init!(sam1; remake=true, prn=false)

        # Simulate - no gravity, no wind → stationary
        for _ in 1:10
            next_step!(sam1; dt=0.01, vsm_interval=0)
        end
        vel_no_grav = copy(
            sam1.sys_struct.points[:test_point].vel_w)
        @test norm(vel_no_grav) < 1e-6

        # Now create a new SAM with gravity enabled,
        # loading from the same .bin cache
        set2 = Settings("system.yaml")
        set2.g_earth = 9.81
        set2.v_wind = 0.0
        sys2 = load_sys_struct_from_yaml(
            yaml_path; system_name="deser_test_2", set=set2)
        sam2 = SymbolicAWEModel(set2, sys2)
        # remake=false → loads from cache
        init!(sam2; remake=false, prn=false)

        @test sam2.set === sam2.sys_struct.set
        @test sam2.set.g_earth == 9.81

        # Simulate - with gravity, point should accelerate
        for _ in 1:10
            next_step!(sam2; dt=0.01, vsm_interval=0)
        end
        vel_with_grav = copy(
            sam2.sys_struct.points[:test_point].vel_w)
        # Z velocity should be negative (falling in NED)
        @test vel_with_grav[3] < -0.5
    end
end
