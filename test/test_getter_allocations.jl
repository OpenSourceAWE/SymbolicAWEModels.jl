# SPDX-FileCopyrightText: 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_getter_allocations.jl - state extraction allocation tests
#
# The monolithic in-place getter (`get_all_state`) writes the integrator state
# into the SystemStructure fields with zero allocations. This test asserts that
# `update_sys_struct!` (and the getter it calls) allocate nothing per step.
#
# The model below carries no wing, which for a long time made this file blind to
# the read-back path most kites actually take: reading `sys.wings` rebuilt its
# filtered view on every access, and the per-step allocation that cost went
# unnoticed here because the fixture had no wing to filter for. The second testset
# covers a model that has one.

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

@isdefined(test_init!) || include(joinpath(@__DIR__, "util.jl"))

using Test
using SymbolicAWEModels
using SymbolicAWEModels: VortexStepMethod
using KiteUtils

const GETTER_ALLOC_YAML = """
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

"""
Measure allocations of one `sync_params!` call (typed function barrier).

The counterpart of the getter: the getter carries the integrator's state out to
the struct, this carries the struct's parameters back in, and `next_step!` runs it
twice a step. Its readers are closures of many types, so without grouping them by
type it costs a dynamic dispatch and a boxed return per parameter.
"""
measure_sync_alloc(prob, integ, sys_struct) =
    @allocated SymbolicAWEModels.sync_params!(prob.param_sync, integ, sys_struct)

"""
    SYNC_ALLOC_BUDGET

What one `sync_params!` may allocate, per backend. The kernel backend writes its
own flat buffer and allocates nothing. The monolith finishes through an MTK `setp`
setter and keeps `Any` buffers for its array and callable kinds, which together
cost a fixed 256 B. Neither grows with the model: the scalar kind, the one with an
entry per parameter, is zero on both.
"""
SYNC_ALLOC_BUDGET = Dict(:KernelBackend => 0, :MonolithBackend => 256)

@testset "Getter allocations" begin
    tmpdir = mktempdir()
    write(joinpath(tmpdir, "getter_alloc_geometry.yaml"), GETTER_ALLOC_YAML)
    write(joinpath(tmpdir, "settings.yaml"), GETTER_ALLOC_SETTINGS)
    write(joinpath(tmpdir, "system.yaml"), "system:\n  sim_settings: settings.yaml\n")

    set_data_path(tmpdir)
    set = Settings("system.yaml")
    yaml_path = joinpath(tmpdir, "getter_alloc_geometry.yaml")
    sys = load_sys_struct_from_yaml(yaml_path; system_name="segment_test", set=set)
    for backend in (MonolithBackend(), KernelBackend())
        sys = load_sys_struct_from_yaml(yaml_path; system_name="segment_test", set=set)
        sam = SymbolicAWEModel(set, sys; backend)
        # init!, not test_init!: the RHS allocation check belongs to the tests
        # about the RHS, and this file is about the sync paths either side of it.
        init!(sam)
        next_step!(sam; dt=0.05, vsm_interval=0)

        prob = sam.prob
        integ = sam.integrator
        sys_struct = sam.sys_struct
        getter = prob.get_all_state

        # warm up, then assert zero allocation
        measure_getter_alloc(getter, integ, sys_struct)
        measure_update_alloc(prob, integ, sys_struct)
        measure_sync_alloc(prob, integ, sys_struct)

        @testset "$(nameof(typeof(backend)))" begin
            @test measure_getter_alloc(getter, integ, sys_struct) == 0
            @test measure_update_alloc(prob, integ, sys_struct) == 0
            allocated = measure_sync_alloc(prob, integ, sys_struct)
            @test allocated <= SYNC_ALLOC_BUDGET[nameof(typeof(backend))]
            allocated <= SYNC_ALLOC_BUDGET[nameof(typeof(backend))] ||
                println("sync_params! allocated $allocated B, over the " *
                        "$(SYNC_ALLOC_BUDGET[nameof(typeof(backend))]) B budget")
        end
    end
end

"""
    WING_ALLOC_BUDGET

What one `update_sys_struct!` may allocate on a model that carries a wing, by
backend. Zero is out of reach on that path for now: a fitted wing's pose is rebuilt
by `wing_kinematics_from_points!` and its scalars by `write_wing_scalars!`, both of
which work in freshly allocated vectors, and `quaternion_to_rotation_matrix` returning
a heap matrix is over a third of the kernel's total on its own. These are ratchets
against a new allocation inside a per-point loop, not targets — measured at 14.6 kB
and 2.3 kB on the particle wing, 10.0 kB and 1.3 kB on the rigid one.
"""
WING_ALLOC_BUDGET = Dict(:KernelBackend => 22_000, :MonolithBackend => 4_000)

@testset "Getter allocations with a wing" begin
    root = mktempdir()
    geometries = ("particle" => "particle_structural_geometry.yaml",
                  "rigid" => "rigid_structural_geometry.yaml")
    for backend in (MonolithBackend(), KernelBackend())
        backend_name = nameof(typeof(backend))
        data_path = joinpath(root, string(backend_name), "2plate_kite")
        mkpath(dirname(data_path))
        cp(joinpath(dirname(@__DIR__), "data", "2plate_kite"), data_path; force=true)
        for (name, geometry) in geometries
            set_data_path(data_path)
            set = Settings("system.yaml")
            vsm_set = VortexStepMethod.VSMSettings(
                joinpath(data_path, "vsm_settings.yaml"); data_prefix=false)
            sys = load_sys_struct_from_yaml(joinpath(data_path, geometry);
                system_name="wing_alloc_$name", set=set, vsm_set=vsm_set)
            sam = SymbolicAWEModel(set, sys; backend)
            init!(sam; prn=false, remake=true)
            next_step!(sam; dt=0.05, vsm_interval=0)

            prob, integ, sys_struct = sam.prob, sam.integrator, sam.sys_struct
            measure_update_alloc(prob, integ, sys_struct)

            @testset "$name wing / $backend_name" begin
                allocated = measure_update_alloc(prob, integ, sys_struct)
                @test allocated <= WING_ALLOC_BUDGET[backend_name]
                allocated <= WING_ALLOC_BUDGET[backend_name] ||
                    println("update_sys_struct! allocated $allocated B, " *
                            "over the $(WING_ALLOC_BUDGET[backend_name]) B budget")
            end
        end
    end
end
