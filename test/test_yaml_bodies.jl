# SPDX-FileCopyrightText: 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_yaml_bodies.jl - Validate loading rigid bodies, a tube between them and a
# body-anchored point from a struc_geometry YAML (the `bodies` and `tubes` blocks
# and the point `body_idx`/`anchor_b` fields of load_authoring). A
# clamped root body and a free tip body joined by one tube must load with
# references resolved and sag under gravity to where the same beam settled when
# it was written as a `timoshenko_joints` or `elastic_joints` row.

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
  headers: [name, extra_mass, inertia_principal, pos, type]
  data:
    - [nodeA, 1.0, [0.01, 0.01, 0.01], [0.0, 0.0, 0.0], STATIC]
    - [nodeB, 1.0, [0.01, 0.01, 0.01], [1.0, 0.0, 0.0], DYNAMIC]

TUBES

points:
  headers: [name, pos_cad, type, body_idx, anchor_b]
  data:
    - [tip_anchor, [1.0, 0.0, 0.0], BODY_STATIC, nodeB, [0.0, 0.0, 0.0]]
"""

TUBE_ROWS = Dict(
    "timoshenko" => """
tubes:
  headers: [name, bodies, diameter, pressure, law, model, EA, GA, GJ, EIy, EIz,
            shear_coeff, damping]
  data:
    - [tube, [nodeA, nodeB], 0.1, 30000.0, breukels2011, timoshenko, 10000.0,
       1500.0, 50.0, 100.0, 100.0, 0.8333, 0.05]
""",
    "elastic" => """
tubes:
  headers: [name, bodies, diameter, pressure, law, model, anchor_a, anchor_b,
            stiffness_axial, stiffness_shear, stiffness_torsion, stiffness_bending,
            damping]
  data:
    - [tube, [nodeA, nodeB], 0.1, 30000.0, breukels2011, elastic, [0.5, 0.0, 0.0],
       [-0.5, 0.0, 0.0], 10000.0, 1500.0, 50.0, 100.0, 0.05]
""")

# Where nodeB came to rest after 400 steps of 5 ms under v0.19.0, with the
# same beam written as a `timoshenko_joints` or `elastic_joints` row.
SETTLED_TIP = Dict(
    "timoshenko" => (pos_w = [0.9992186305467761, 0.0, -0.04051691939002452],
                     Q_b_to_w = [0.999699629386848, 0.0, 0.0245081823846388, 0.0]),
    "elastic" => (pos_w = [0.9993996077904469, 0.0, -0.031035544080983943],
                  Q_b_to_w = [0.999699758818729, 0.0, 0.024502902231676792, 0.0]))

function load_beam(model_name)
    pkg_root = dirname(@__DIR__)
    data_path = joinpath(mktempdir(), "2plate_kite")
    cp(joinpath(pkg_root, "data", "2plate_kite"), data_path; force=true)
    write(joinpath(data_path, "settings.yaml"), SETTINGS_YAML)
    write(joinpath(data_path, "system.yaml"),
        "system:\n  sim_settings: settings.yaml\n")
    struc_yaml = joinpath(data_path, "yaml_bodies_test.yaml")
    write(struc_yaml, replace(STRUC_YAML, "TUBES" => TUBE_ROWS[model_name]))
    set_data_path(data_path)
    set = Settings("system.yaml")
    sys, transforms = load_authoring(struc_yaml;
        name="yaml_bodies_$model_name", set)
    return set, sys, transforms
end

@testset "YAML body + tube loading" begin
    set, sys, transforms = load_beam("timoshenko")

    @testset "Structure wiring from YAML" begin
        @test length(sys.bodies) == 2
        @test sys.bodies[:nodeA].fixed
        @test !sys.bodies[:nodeB].fixed
        @test sys.bodies[:nodeA].extra_mass == 1.0
        @test length(sys.tubes) == 1
        tube = sys.tubes[:tube]
        model = run_value(SymbolicAWEModel(set, sys), :tubes, :tube, :model)
        @test tube.bodies[1] == KiteGeometry.row_index(sys.bodies, :nodeA)
        @test tube.bodies[2] == KiteGeometry.row_index(sys.bodies, :nodeB)
        @test (tube.diameter, tube.pressure, sys.materials[tube.material].name) ==
              (0.1, 30000.0, :breukels2011)
        @test model isa TimoshenkoTube
        @test model.rest_length ≈ 1.0     # taken from placed geometry
        @test model.EA ≈ 10000.0
        anchor = sys.points[:tip_anchor]
        @test !anchor.fixed
        @test anchor.body == KiteGeometry.row_index(sys.bodies, :nodeB)
    end

    @testset "a tube row without rigidities takes them from its law" begin
        # A row leaves its rigidities unset; the run's element derives them.
        material(law) = Material(; name=law, law=String(law),
                                 parameters=KiteGeometry.tube_law_parameters(law))
        function run_tube(tube, law=:breukels2011)
            bodies = [Body(; name, pos_ENU, extra_mass=1.0,
                           extra_inertia_KA=Matrix(Diagonal([0.01, 0.01, 0.01])))
                      for (name, pos_ENU) in ((:a, [0.0, 0.0, 0.0]), (:b, [1.0, 0.0, 0.0]))]
            strut = system_definition("yaml_bodies_strut", set;
                materials=[material(law)], bodies, tubes=[tube])
            return run_value(SymbolicAWEModel(set, strut), :tubes, 1, :model)
        end
        model = run_tube(Tube(; name=:strut, bodies=(1, 2), diameter=0.1, pressure=30000.0,
                             material=:breukels2011))
        EA, GA, EI, GJ = tube_linear_rigidities(0.05, 0.3)
        @test (model.EA, model.GA, model.GJ) == (EA, GA, GJ)
        @test model.EIy == model.EIz == EI
        partial = run_tube(Tube(; name=:strut, bodies=(1, 2), diameter=0.1,
                                pressure=30000.0, material=:breukels2011,
                                model=TimoshenkoTubeModel(; EA=1.0)))
        @test (partial.EA, partial.GA) == (1.0, GA)
        @test_throws "unknown" run_tube(Tube(; name=:strut, bodies=(1, 2), diameter=0.1,
                                             pressure=3e4, material=:no_such_law),
                                        :no_such_law)
    end

    for model_name in ("timoshenko", "elastic")
        @testset "$model_name tube settles where its joint row did" begin
            set, sys, transforms = model_name == "timoshenko" ?
                (set, sys, transforms) : load_beam(model_name)
            sam = SymbolicAWEModel(set, sys; transforms)
            test_init!(sam; prn=false)
            for _ in 1:400
                next_step!(sam; dt=0.005, vsm_interval=0)
            end
            tip(quantity) = run_value(sam, :bodies, :nodeB, quantity)
            @test tip(:pos_w) ≈ SETTLED_TIP[model_name].pos_w atol=1e-8
            @test tip(:Q_b_to_w) ≈ SETTLED_TIP[model_name].Q_b_to_w atol=1e-8
        end
    end
end
nothing
