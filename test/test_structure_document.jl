# SPDX-FileCopyrightText: 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_structure_document.jl - SystemStructure → structure document →
# SystemStructure, in both encodings, against a vendored copy of
# structure_schema.yml and the golden document in test/data. Regenerate the
# golden by pointing save_structure_document at data/2plate_kite_structure.yml
# with the description and note below.

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

using Test
using SymbolicAWEModels
using KiteUtils
using VortexStepMethod
using JSONSchema
using YAML

GOLDEN_DESCRIPTION = "Resolved structure of the 2-plate kite shipped in data/"
GOLDEN_NOTE = "Written by SymbolicAWEModels; the golden file of " *
              "test_structure_document.jl."

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
    bodies_and_joints(set) -> SystemStructure

A clamped root body, a free tip body and a hub, linked by one Timoshenko beam and
one elastic joint, with a segment between two body-anchored points.
"""
function bodies_and_joints(set)
    bodies = Body[
        Body(:root; mass=1.0, inertia_principal=[0.01, 0.01, 0.01],
             pos=[0.0, 0.0, 0.0], type=STATIC),
        Body(:tip; mass=1.0, inertia_principal=[0.01, 0.01, 0.01],
             pos=[1.0, 0.0, 0.0]),
        Body(:hub; mass=0.5, inertia_principal=[0.02, 0.03, 0.04],
             pos=[2.0, 0.0, 0.0], com_offset_b=[0.01, 0.0, 0.0])]
    points = Point[
        Point(:tip_anchor, [1.0, 0.0, 0.0], BODY_STATIC; body=:tip),
        Point(:hub_anchor, [2.0, 0.0, 0.0], BODY_STATIC; body=:hub)]
    segments = Segment[
        Segment(:link, :tip_anchor, :hub_anchor, 1000.0, 10.0, 0.004;
                l0=1.0, density=724.0)]
    timoshenko_joints = TimoshenkoJoint[
        TimoshenkoJoint(:beam, :root, :tip; EA=1.0e4, GA=1500.0, GJ=50.0,
                        EIy=100.0, EIz=100.0, shear_coeff=0.8333, damping=0.05,
                        rest_length=1.0, radius=0.06)]
    elastic_joints = ElasticJoint[
        ElasticJoint(:spring, :tip, :hub; anchor_a=[0.1, 0.0, 0.0],
                     anchor_b=[-0.1, 0.0, 0.0], stiffness_axial=1.2e5,
                     stiffness_shear=4.0e4, stiffness_torsion=900.0,
                     stiffness_bending=1500.0, damping=0.002)]
    return SystemStructure("structure_document_bodies", set; points, segments,
        bodies, elastic_joints, timoshenko_joints, prn=false)
end

@testset verbose = true "Structure document" begin
    pkg_root = dirname(@__DIR__)
    set_data_path(joinpath(pkg_root, "data", "2plate_kite"))
    set = Settings("system.yaml")
    vsm_set = VortexStepMethod.VSMSettings(
        joinpath(get_data_path(), "vsm_settings.yaml"); data_prefix=false)
    schema = Schema(YAML.load_file(joinpath(@__DIR__, "data",
                                            "structure_schema.yml")))
    tmpdir = mktempdir()

    @testset "connectivity_sha hashes the preimage the schema documents" begin
        # The schema's worked example reads `8;1,2;2,3;3,4;4,5;4,7;`.
        @test SymbolicAWEModels.connectivity_sha(
            8, [(1, 2), (2, 3), (3, 4), (4, 5), (4, 7)]) ==
            "d98529e23af6047ce9f49f172743d5dcef053f6a768a3b53de6d47260afd60b7"
    end

    @testset "the 2plate kite writes the golden document" begin
        document = structure_document(two_plate_kite(set, vsm_set);
            description=GOLDEN_DESCRIPTION, note=GOLDEN_NOTE)
        @test isnothing(JSONSchema.validate(schema, document))
        @test YAML.load_file(joinpath(@__DIR__, "data",
                                      "2plate_kite_structure.yml")) == document
    end

    @testset "the 2plate kite round-trips through $extension" for
            extension in (".yml", ".json")
        sys = two_plate_kite(set, vsm_set)
        document = structure_document(sys; description=GOLDEN_DESCRIPTION,
                                      note=GOLDEN_NOTE)
        path = joinpath(tmpdir, "2plate_kite_structure" * extension)
        save_structure_document(path, sys; description=GOLDEN_DESCRIPTION,
                                note=GOLDEN_NOTE)
        reread = load_structure_document(path; set, vsm_set, prn=false)
        @test structure_document(reread; description=GOLDEN_DESCRIPTION,
                                 note=GOLDEN_NOTE) == document
        @test length(reread.points) == length(sys.points)
        @test reread.wings[:main_wing].mass ≈ sys.wings[:main_wing].mass
        @test reread.stations[:left].point_idxs == sys.stations[:left].point_idxs
    end

    @testset "bodies and joints round-trip" begin
        sys = bodies_and_joints(set)
        document = structure_document(sys)
        @test isnothing(JSONSchema.validate(schema, document))
        path = joinpath(tmpdir, "bodies_and_joints.yml")
        save_structure_document(path, sys)
        reread = load_structure_document(path; set, prn=false)
        @test structure_document(reread) == document
        @test reread.bodies[:root].type == STATIC
        @test reread.timoshenko_joints[:beam].body_a_idx ==
              reread.bodies[:root].idx
        @test reread.elastic_joints[:spring].stiffness_axial ≈ 1.2e5
        @test reread.points[:hub_anchor].body_idx == reread.bodies[:hub].idx
    end

    @testset "a document whose connectivity_sha does not describe it is refused" begin
        document = structure_document(bodies_and_joints(set))
        document["metadata"]["connectivity_sha"] = repeat("0", 64)
        @test_throws "does not describe" sys_struct_from_document(document; set)
    end

    @testset "a nonlinear stiffness law cannot be written yet" begin
        sys = bodies_and_joints(set)
        sys.segments[:link].unit_stiffness = strain -> 1000.0 * strain
        @test_throws "awesIO#1" structure_document(sys)
    end
end
