# SPDX-FileCopyrightText: 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

# test_yaml_update.jl - Tests for update_sys_struct_from_yaml!
#
# Verifies:
# 1. Unchanged YAML round-trip: pos_cad and l0 unchanged
# 2. Modified pos_cad: only the changed point is updated
# 3. Modified segment l0: segment updated correctly
# 4. l0=nothing in YAML: auto-calc from pos_cad

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

using Test
using SymbolicAWEModels
using SymbolicAWEModels: KVec3
using LinearAlgebra: norm
using KiteUtils

const YAML_UPDATE_BASE = """
points:
  headers: [name, pos_cad, type, wing_idx, transform_idx,
            extra_mass, body_frame_damping,
            world_frame_damping, area, drag_coeff]
  data:
    - [pt_a, [1.0, 2.0, 3.0], DYNAMIC, nothing, nothing,
       1.0, 0.0, 0.0, 0.0, 0.0]
    - [pt_b, [4.0, 5.0, 6.0], DYNAMIC, nothing, nothing,
       1.0, 0.0, 0.0, 0.0, 0.0]

segments:
  headers: [name, point_i, point_j, l0,
            diameter_mm, unit_stiffness, unit_damping,
            compression_frac]
  data:
    - [seg_ab, pt_a, pt_b, 10.0,
       1.0, 5000.0, 10.0, 1.0]
"""

const SETTINGS_YAML_UPDATE = """
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
    profile_law: 0
"""

@testset "update_sys_struct_from_yaml!" begin
    tmpdir = mktempdir()
    yaml_path = joinpath(tmpdir, "geometry.yaml")
    write(yaml_path, YAML_UPDATE_BASE)

    settings_path = joinpath(tmpdir, "settings.yaml")
    write(settings_path, SETTINGS_YAML_UPDATE)
    system_path = joinpath(tmpdir, "system.yaml")
    write(system_path, """
system:
  sim_settings: settings.yaml
""")

    set_data_path(tmpdir)
    set = Settings("system.yaml")

    sys = load_sys_struct_from_yaml(yaml_path;
        system_name="yaml_update_test", set=set)

    # --------------------------------------------------------
    # Test 1: Unchanged YAML round-trip
    # --------------------------------------------------------
    @testset "Unchanged round-trip" begin
        orig_a = copy(sys.points[:pt_a].pos_cad)
        orig_b = copy(sys.points[:pt_b].pos_cad)
        orig_l0 = sys.segments[:seg_ab].l0

        update_sys_struct_from_yaml!(sys, yaml_path)

        @test sys.points[:pt_a].pos_cad == orig_a
        @test sys.points[:pt_b].pos_cad == orig_b
        @test sys.segments[:seg_ab].l0 == orig_l0
    end

    # --------------------------------------------------------
    # Test 2: Modified point position
    # --------------------------------------------------------
    @testset "Modified point position" begin
        modified_yaml = replace(
            YAML_UPDATE_BASE,
            "[1.0, 2.0, 3.0]" => "[10.0, 20.0, 30.0]")
        mod_path = joinpath(tmpdir, "modified_pos.yaml")
        write(mod_path, modified_yaml)

        orig_b = copy(sys.points[:pt_b].pos_cad)

        update_sys_struct_from_yaml!(sys, mod_path)

        @test sys.points[:pt_a].pos_cad ==
              KVec3(10.0, 20.0, 30.0)
        # pt_b should be unchanged
        @test sys.points[:pt_b].pos_cad == orig_b
    end

    # --------------------------------------------------------
    # Test 3: Modified segment l0
    # --------------------------------------------------------
    @testset "Modified segment l0" begin
        modified_yaml = replace(
            YAML_UPDATE_BASE,
            "10.0,\n       1.0, 5000.0" =>
            "42.5,\n       1.0, 5000.0")
        mod_path = joinpath(tmpdir, "modified_l0.yaml")
        write(mod_path, modified_yaml)

        update_sys_struct_from_yaml!(sys, mod_path)

        @test sys.segments[:seg_ab].l0 == 42.5
    end

    # --------------------------------------------------------
    # Test 4: l0=nothing auto-calculates from pos_cad
    # --------------------------------------------------------
    @testset "l0=nothing auto-calc from pos_cad" begin
        # First reset to base YAML so pos_cad is known
        write(yaml_path, YAML_UPDATE_BASE)
        update_sys_struct_from_yaml!(sys, yaml_path)

        expected_l0 = norm(
            KVec3(1.0, 2.0, 3.0) - KVec3(4.0, 5.0, 6.0))

        nothing_yaml = replace(
            YAML_UPDATE_BASE,
            "10.0,\n       1.0, 5000.0" =>
            "nothing,\n       1.0, 5000.0")
        nothing_path = joinpath(tmpdir, "nothing_l0.yaml")
        write(nothing_path, nothing_yaml)

        update_sys_struct_from_yaml!(sys, nothing_path)

        @test sys.segments[:seg_ab].l0 ≈ expected_l0
    end

    rm(tmpdir; recursive=true)
end
nothing
