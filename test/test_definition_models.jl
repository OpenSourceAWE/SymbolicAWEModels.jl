# SPDX-FileCopyrightText: 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_definition_models.jl - SymbolicAWEModels' own authoring YAMLs load through
# KiteGeometry, its columns read into the models it registers.

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

using Test
using SymbolicAWEModels
using KiteUtils

"""The component of `components` named `name`."""
named(components, name) = only(filter(row -> row.name == Symbol(name), components))

@testset verbose = true "Definition models" begin
    kite = joinpath(dirname(@__DIR__), "data", "2plate_kite")
    set_data_path(kite)
    set = Settings("system.yaml")
    load_kite(dynamics) =
        first(load_authoring(joinpath(kite, "$(dynamics)_structural_geometry.yaml"); set))

    @testset "the $dynamics 2-plate kite reads SAM's columns" for
            (dynamics, point_type, station_type, moment_frac, damping) in
            (("particle", "DYNAMIC", "STATIC", 0.75, 50.0),
             ("rigid", "BODY_STATIC", "DYNAMIC", 0.25, 100.0))
        system = load_kite(dynamics)
        le_left = named(system.points, "le_left")
        @test le_left.model isa SymbolicAWEModels.PointModel
        @test le_left.model.type == point_type
        @test le_left.model.body_frame_damping == 10.0
        @test le_left.model.world_frame_damping == 0.0
        @test le_left.area == 0.0
        @test le_left.drag_coeff == 0.0
        @test isnothing(le_left.tube)
        center = named(system.spanwise_stations, "center")
        @test center.type == station_type
        @test center.moment_frac == moment_frac
        @test center.damping == damping
        @test isnothing(center.stiffness)
        pulley = named(system.pulleys, "left")
        @test pulley.type == "DYNAMIC"
        @test isnothing(pulley.efficiency)
        @test only(system.wings).model isa SymbolicAWEModels.WingModel
        @test length(system.canopy_faces) == 2
    end

    @testset "a wing's aero_z_offset reads into its model" begin
        @test only(load_kite("rigid").wings).aero_z_offset == 0.0
        @test isnothing(only(load_kite("particle").wings).aero_z_offset)
    end

    @testset "a body's dynamics type reads into its model" begin
        authored = read(joinpath(kite, "rigid_structural_geometry.yaml"), String)
        file = joinpath(mktempdir(), "rigid_structural_geometry.yaml")
        write(file, authored * """
            bodies:
              headers: [name, extra_mass, pos, inertia_principal, type]
              data:
                - [hub, 1.0, [0.0, 0.0, 0.0], [1.0, 1.0, 1.0], STATIC]
            """)
        body = named(first(load_authoring(file; set)).bodies, "hub")
        @test body.model isa SymbolicAWEModels.BodyModel
        @test body.type == "STATIC"
    end

    @testset "a point's damping reads as one number or one per axis" begin
        authored = read(joinpath(kite, "rigid_structural_geometry.yaml"), String)
        authored_row = "main_wing, main_transform, 0.1, 10.0, 0.0,"
        file = joinpath(mktempdir(), "rigid_structural_geometry.yaml")
        damped_row = "main_wing, main_transform, 0.1, [1, 2, 3], 4,"
        write(file, replace(authored, authored_row => damped_row; count=1))
        point = first(first(load_authoring(file; set)).points)
        @test point.model.body_frame_damping == [1.0, 2.0, 3.0]
        @test point.model.world_frame_damping === 4.0
    end

    @testset "a tube reads into the model it names, a point into its anchor" begin
        beam(tubes) = """
            bodies:
              headers: [name, extra_mass, inertia_principal, pos, type]
              data:
                - [nodeA, 1.0, [0.01, 0.01, 0.01], [0.0, 0.0, 0.0], STATIC]
                - [nodeB, 1.0, [0.01, 0.01, 0.01], [1.0, 0.0, 0.0], DYNAMIC]
            points:
              headers: [name, pos_cad, type, body_idx, anchor_b]
              data:
                - [tip, [1.0, 0.0, 0.0], BODY_STATIC, nodeB, [0.0, 0.0, 0.0]]
            tubes:
            """ * tubes
        function load_beam(tubes)
            file = joinpath(mktempdir(), "beam.yaml")
            write(file, beam(tubes))
            return first(load_authoring(file; set))
        end
        system = load_beam("""
              headers: [name, bodies, diameter, pressure, law, model, anchor_a,
                        stiffness_axial, stiffness_shear, stiffness_torsion,
                        stiffness_bending]
              data:
                - [tube, [nodeA, nodeB], 0.1, 30000.0, breukels2011, elastic,
                   [0.5, 0.0, 0.0], 10000.0, 1500.0, 50.0, 100.0]
            """)
        tube = only(system.tubes)
        @test tube.model isa SymbolicAWEModels.ElasticTubeModel
        @test tube.stiffness_bending == 100.0
        @test tube.anchor_a == [0.5, 0.0, 0.0]
        @test only(system.points).anchor_b == zeros(3)
        @test named(system.bodies, "nodeA").type == "STATIC"
        tube = only(load_beam("""
              headers: [name, bodies, diameter, pressure, law, EA, damping]
              data:
                - [tube, [nodeA, nodeB], 0.1, 30000.0, breukels2011, 10000.0, 0.05]
            """).tubes)
        @test tube.model isa SymbolicAWEModels.TimoshenkoTubeModel
        @test tube.EA == 10000.0
        @test isnothing(tube.GJ)
    end

    @testset "a point on a wing naming a tube and no body rides the tube" begin
        authored = replace(read(joinpath(kite, "rigid_structural_geometry.yaml"), String),
                           "\r\n" => "\n")
        rider_row = "    - [rider, [0.0, 0.0, 2.5], BODY_STATIC, main_wing, " *
                    "main_transform, 0.1, 10.0, 0.0, 0.0, 0.0, nothing, strut]\n"
        authored = replace(authored,
            "area, drag_coeff, body_idx]\n  data:\n" =>
            "area, drag_coeff, body_idx, tube]\n  data:\n" * rider_row; count=1)
        file = joinpath(mktempdir(), "rigid_structural_geometry.yaml")
        write(file, authored * """
            bodies:
              headers: [name, extra_mass, inertia_principal, pos, type]
              data:
                - [nodeA, 1.0, [0.01, 0.01, 0.01], [0.0, 0.0, 2.5], DYNAMIC]
                - [nodeB, 1.0, [0.01, 0.01, 0.01], [0.0, 0.5, 2.5], DYNAMIC]
            tubes:
              headers: [name, bodies, diameter, pressure, law, EA]
              data:
                - [strut, [nodeA, nodeB], 0.1, 30000.0, breukels2011, 10000.0]
            """)
        sys = load_sys_struct_from_yaml(file; system_name="tube_rider", set,
                                        aero_mode=AeroNone(), prn=false)
        rider = sys.points[:rider]
        @test rider.tube_idx == sys.tubes[:strut].idx
    end
end
