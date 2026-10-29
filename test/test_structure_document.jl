# SPDX-FileCopyrightText: 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_structure_document.jl - SystemStructure → awesIO structure document →
# SystemStructure, in both encodings and through a run's log, against the
# structure_schema.yml KiteGeometry ships and the golden document in test/data.
# Regenerate the golden by pointing save_structure_document at
# data/2plate_kite_structure.yml.

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

using Test
using SymbolicAWEModels
using SymbolicAWEModels: KiteGeometry
using KiteUtils
using VortexStepMethod
using JSON
using JSONSchema
using YAML

"""
    two_plate_kite(set, vsm_set) -> SystemStructure

The rigid-dynamics 2-plate kite from `data/2plate_kite`, the system the golden
document describes.
"""
function two_plate_kite(set, vsm_set)
    return load_sys_struct_from_yaml(
        joinpath(get_data_path(), "rigid_structural_geometry.yaml");
        system_name="2plate_kite", set, vsm_set, prn=false)
end

"""
    bodies_and_tubes(set) -> SystemStructure

A clamped root body, a free tip body and a hub whose own mass sits off its origin,
linked by one Timoshenko tube and one elastic tube, with a segment between two
body-anchored points.
"""
function bodies_and_tubes(set)
    bodies = Body[
        Body(:root; extra_mass=1.0, inertia_principal=[0.01, 0.01, 0.01],
             pos=[0.0, 0.0, 0.0], type=STATIC),
        Body(:tip; extra_mass=1.0, inertia_principal=[0.01, 0.01, 0.01],
             pos=[1.0, 0.0, 0.0]),
        Body(:hub; extra_mass=0.5, inertia_principal=[0.02, 0.03, 0.04],
             pos=[2.0, 0.0, 0.0], com_offset_b=[0.01, 0.0, 0.0])]
    points = Point[
        Point(:tip_anchor, [1.0, 0.0, 0.0], BODY_STATIC; body=:tip),
        Point(:hub_anchor, [2.0, 0.0, 0.0], BODY_STATIC; body=:hub)]
    segments = Segment[
        Segment(:link, :tip_anchor, :hub_anchor, 1000.0, 10.0, 0.004;
                l0=1.0, density=724.0)]
    tubes = Tube[
        Tube(:beam, :root, :tip; diameter=0.12, pressure=3.0e4,
             model=TimoshenkoTube(EA=1.0e4, GA=1500.0, GJ=50.0, EIy=100.0,
                                  EIz=100.0, shear_coeff=0.8333, damping=0.05)),
        Tube(:spring, :tip, :hub; diameter=0.1, pressure=3.0e4,
             anchor_a=[0.1, 0.0, 0.0], anchor_b=[-0.1, 0.0, 0.0],
             model=ElasticTube(stiffness_axial=1.2e5, stiffness_shear=4.0e4,
                               stiffness_torsion=900.0, stiffness_bending=1500.0,
                               damping=0.002))]
    return SystemStructure("structure_document_bodies", set; points, segments,
        bodies, tubes, prn=false)
end

"""
    pretensioned_tether(set; winches=Winch[]) -> SystemStructure

A ground point and a top point on one two-segment tether whose rest length is 99 % of
its length, reeled by `winches`.
"""
function pretensioned_tether(set; winches=Winch[])
    points = Point[Point(:ground, [0.0, 0.0, 0.0], STATIC),
                   Point(:top, [0.0, 0.0, 10.0], DYNAMIC; extra_mass=1.0)]
    tethers = Tether[Tether(:line, 10.0; start_point=:ground, end_point=:top,
                            n_segments=2, stretch_frac=0.99)]
    return SystemStructure("structure_document_tether", set; points, tethers, winches,
        prn=false)
end

"""World position of `body`'s centre of mass."""
com_w(body) = body.pos_ENU .+ SymbolicAWEModels.initial_rotation(body) * body.com_offset_b

"""World position of `tube`'s end on body B in `sys`."""
function tube_end(sys, tube)
    body = sys.bodies[tube.body_b_idx]
    return body.pos_ENU .+ SymbolicAWEModels.initial_rotation(body) * tube.anchor_b_b
end

"""The `header` cells of the rows of `block` in `document`, by row name."""
function cells(document, block, header)
    table = document[block]
    column = findfirst(==(header), table["headers"])
    return Dict(row[1] => row[column] for row in table["data"])
end

@testset verbose = true "Structure document" begin
    pkg_root = dirname(@__DIR__)
    set_data_path(joinpath(pkg_root, "data", "2plate_kite"))
    set = Settings("system.yaml")
    vsm_set = VortexStepMethod.VSMSettings(
        joinpath(get_data_path(), "vsm_settings.yaml"); data_prefix=false)
    schema = Schema(YAML.load_file(joinpath(pkgdir(KiteGeometry), "src", "awesio",
                                            "structure_schema.yml")))
    golden = YAML.load_file(joinpath(@__DIR__, "data", "2plate_kite_structure.yml"))
    tmpdir = mktempdir()

    @testset "connectivity_sha hashes the preimage the schema documents" begin
        # The schema's worked example reads `3;1,2;2,3;2;1,2;1;1,2,3;`.
        @test connectivity_sha(
            (3, [(1, 2), (2, 3)]), (2, [(1, 2)]), (1, [(1, 2, 3)])) ==
            "1f9a31aca7b6d655aaae4f5de9872d3fd489f907f90316e8f33e3125b1fd18ea"
    end

    @testset "the 2plate kite writes the golden document" begin
        document = structure_document(two_plate_kite(set, vsm_set))
        @test isnothing(JSONSchema.validate(schema, document))
        @test golden == document
    end

    @testset "a document holds the initial pose, not where a run moved the system" begin
        sys = two_plate_kite(set, vsm_set)
        state = SystemState(sys)
        state.points[:kcu].pos_w .+= [1.0, 0.0, 0.0]
        kcu = cells(structure_document(sys), "points", "pos_ENU")["kcu"]
        @test kcu ≈ sys.points[:kcu].pos_ENU
        @test !(kcu ≈ state.points[:kcu].pos_w)
    end

    @testset "a point names its body, and carries its own mass and drag" begin
        document = structure_document(two_plate_kite(set, vsm_set))
        @test cells(document, "points", "body")["le_left"] == "main_wing"
        @test [cells(document, "points", header)["kcu"]
               for header in ("extra_mass", "area", "drag_coeff")] == [1.0, 0.1, 1.0]
    end

    @testset "a point is its definition row" for
            sys in (two_plate_kite(set, vsm_set), pretensioned_tether(set))
        point = last(sys.points)
        row = sys.definition.points[point.idx]
        point.extra_mass = 2.0
        point.area = 0.3
        point.fix_static = true
        @test (row.extra_mass, row.area, row.fix_static) == (2.0, 0.3, true)
        row.drag_coeff = 0.7
        @test point.drag_coeff == 0.7
        @test point === row
    end

    @testset "a segment is its definition row" for
            sys in (two_plate_kite(set, vsm_set), pretensioned_tether(set))
        segment = last(sys.segments)
        row = sys.definition.segments[segment.idx]
        segment.l0 = 1.5
        segment.diameter = 0.002
        segment.unit_damping = 3.0
        segment.compression_frac = 0.2
        segment.compression_damping_frac = 0.5
        @test (row.rest_length, row.diameter, row.unit_damping, row.compression_frac,
               row.compression_damping_frac) == (1.5, 0.002, 3.0, 0.2, 0.5)
        row.rest_length = 2.5
        @test segment.l0 == 2.5
        @test segment === row
    end

    @testset "a winch is its definition row" for sys in (two_plate_kite(set, vsm_set),
            pretensioned_tether(set;
                                winches=[Winch(:reel, set, [:line]; winch_point=:ground)]))
        winch = only(sys.winches)
        row = sys.definition.winches[winch.idx]
        winch.brake = 1.0
        winch.model.friction_epsilon = 3.0
        @test row.model.brake == 1.0
        document = structure_document(sys)
        @test only(values(cells(document, "winches", "friction_epsilon"))) == 3.0
        @test winch.model.dynamics.friction_epsilon == 3.0
        @test winch === row
    end

    @testset "a segment named by a number keeps that name as its row's" begin
        @test Segment(3, :ground, :top, 1.0, 0.0, 0.001).name === Symbol(3)
    end

    @testset "a station names the wing whose twist it carries" begin
        document = structure_document(two_plate_kite(set, vsm_set))
        @test all(==("main_wing"), values(cells(document, "spanwise_stations", "wing")))
    end

    @testset "a document with wings reads back in its initial pose" begin
        sys = two_plate_kite(set, vsm_set)
        reread = sys_struct_from_document(structure_document(sys); set, vsm_set,
                                          prn=false)
        @test [point.pos_ENU for point in reread.points] ≈
              [point.pos_ENU for point in sys.points]
        @test reread.wings[:main_wing].pos_ENU ≈ sys.wings[:main_wing].pos_ENU
        @test reread.wings[:main_wing].Q_KA_to_ENU ≈ sys.wings[:main_wing].Q_KA_to_ENU
    end

    @testset "bodies and tubes round-trip through $extension" for
            extension in (".yml", ".json")
        sys = bodies_and_tubes(set)
        document = structure_document(sys)
        @test isnothing(JSONSchema.validate(schema, document))
        path = joinpath(tmpdir, "bodies_and_tubes" * extension)
        save_structure_document(path, sys)
        reread = load_structure_document(path; set, prn=false)
        @test structure_document(reread) == document
        @test reread.bodies[:root].type == STATIC
        @test reread.tubes[:beam].body_a_idx == reread.bodies[:root].idx
        @test reread.tubes[:beam].model.EA ≈ 1.0e4
        @test reread.tubes[:spring].model.stiffness_axial ≈ 1.2e5
        @test reread.points[:hub_anchor].body_idx == reread.bodies[:hub].idx
        for name in (:tip, :hub)
            body, reread_body = sys.bodies[name], reread.bodies[name]
            @test reread_body.total_mass ≈ body.total_mass
            @test reread_body.inertia_principal ≈ body.inertia_principal
            @test com_w(reread_body) ≈ com_w(body)
        end
        spring, reread_spring = sys.tubes[:spring], reread.tubes[:spring]
        @test tube_end(reread, reread_spring) ≈ tube_end(sys, spring)
    end

    @testset "a document is read in its initial pose, its tethers' rest lengths kept" begin
        sys = pretensioned_tether(set)
        @test sum(segment.l0 for segment in sys.segments) ≈ 9.9
        document = structure_document(sys)
        reread = sys_struct_from_document(document; set, prn=false)
        @test structure_document(reread) == document
        @test [point.pos_w for point in SystemState(reread).points] ==
              [point.pos_ENU for point in sys.points]
    end

    @testset "a run's log carries the structure document of the system that ran" begin
        sys = two_plate_kite(set, vsm_set)
        sam = SymbolicAWEModel(set, sys)
        init!(sam; prn=false)
        logger = Logger(sam, 3)
        state = SysState(sam)
        for step in 0:2
            step > 0 && next_step!(sam)
            update_sys_state!(state, sam)
            log!(logger, state)
        end
        save_log(logger, sys, "2plate_kite_run"; path=tmpdir)
        run_log = load_log("2plate_kite_run"; path=tmpdir)
        topology = JSON.parse(run_log.metadata["topology"])
        @test isnothing(JSONSchema.validate(schema, topology))
        @test topology["metadata"]["connectivity_sha"] ==
              golden["metadata"]["connectivity_sha"]
        @test length(run_log.syslog) == 3
    end

    @testset "a document whose connectivity_sha does not describe it is refused" begin
        document = structure_document(bodies_and_tubes(set))
        document["metadata"]["connectivity_sha"] = repeat("0", 64)
        @test_throws "does not describe" sys_struct_from_document(document; set)
    end

    @testset "a document of another major version is refused" begin
        document = structure_document(bodies_and_tubes(set))
        document["metadata"]["awesIO_version"] = "1.0.0"
        @test_throws "major version" sys_struct_from_document(document; set)
    end

    @testset "a nonlinear stiffness law cannot be written yet" begin
        points = Point[Point(:ground, [0.0, 0.0, 0.0], STATIC),
                       Point(:top, [0.0, 0.0, 1.0], DYNAMIC; extra_mass=1.0)]
        segments = Segment[Segment(:link, :ground, :top, strain -> 1000.0 * strain,
                                   10.0, 0.004; l0=1.0, density=724.0)]
        sys = SystemStructure("nonlinear_link", set; points, segments, prn=false)
        @test_throws "awesIO#1" structure_document(sys)
    end
end
