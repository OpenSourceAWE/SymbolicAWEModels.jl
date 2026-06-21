# SPDX-FileCopyrightText: 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_getter_allocations.jl - state extraction allocation tests
#
# The monolithic in-place getter (`get_all_state`) writes the integrator state
# into the SystemStructure fields with zero allocations. This test asserts that
# `update_sys_struct!` (and the getter it calls) allocate nothing per step.

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

@isdefined(test_init!) || include(joinpath(@__DIR__, "util.jl"))

using Test
using SymbolicAWEModels
using KiteUtils

const GETTER_ALLOC_YAML = """
materials:
  headers: [name, youngs_modulus, density, damping_per_stiffness]
  data:
    - [test_material, 55000000000.0, 724, 0.00077]

points:
  headers: [name, pos_cad, type, wing_idx, transform_idx, extra_mass, body_frame_damping, world_frame_damping, area, drag_coeff]
  data:
    - [anchor, [0.0, 0.0, 0.0], STATIC, nothing, nothing, 0.0, 0.0, 0.0, 0.0, 0.0]
    - [mass_point, [0.0, 0.0, -10.0], DYNAMIC, nothing, nothing, 1.0, 0.0, 0.0, 0.0, 0.0]

segments:
  headers: [name, point_i, point_j, l0, diameter_mm, unit_stiffness, unit_damping, compression_frac]
  data:
    - [test_segment, anchor, mass_point, 10.0, 5.0, 1000.0, 10.0, 0.1]
"""

const GETTER_ALLOC_SETTINGS = """
system:
    log_file: "data/segment_test"
    g_earth:     9.81

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

"Measure allocations of one in-place getter call (typed function barrier)."
measure_getter_alloc(getter, integ, sys_struct) =
    @allocated getter(integ, sys_struct)

"Measure allocations of one `update_sys_struct!` call (typed function barrier)."
measure_update_alloc(prob, integ, sys_struct) =
    @allocated SymbolicAWEModels.update_sys_struct!(prob, integ, sys_struct)

@testset "Getter allocations" begin
    tmpdir = mktempdir()
    write(joinpath(tmpdir, "getter_alloc_geometry.yaml"), GETTER_ALLOC_YAML)
    write(joinpath(tmpdir, "settings.yaml"), GETTER_ALLOC_SETTINGS)
    write(joinpath(tmpdir, "system.yaml"), "system:\n  sim_settings: settings.yaml\n")

    set_data_path(tmpdir)
    set = Settings("system.yaml")
    yaml_path = joinpath(tmpdir, "getter_alloc_geometry.yaml")
    sys = load_sys_struct_from_yaml(yaml_path; system_name="segment_test", set=set)
    sam = SymbolicAWEModel(set, sys)
    test_init!(sam)
    next_step!(sam; dt=0.05, vsm_interval=0)

    prob = sam.prob
    integ = sam.integrator
    sys_struct = sam.sys_struct
    getter = prob.get_all_state

    # warm up, then assert zero allocation
    measure_getter_alloc(getter, integ, sys_struct)
    measure_update_alloc(prob, integ, sys_struct)

    @test measure_getter_alloc(getter, integ, sys_struct) == 0
    @test measure_update_alloc(prob, integ, sys_struct) == 0
end
