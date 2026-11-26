# SPDX-FileCopyrightText: 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_membrane.jl - membrane triangles: forces from their energy, wrinkling, and a
# membrane canopy in the time simulation and in find_equilibrium!

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

@isdefined(test_init!) || include(joinpath(@__DIR__, "util.jl"))

using Test
using SymbolicAWEModels
using SymbolicAWEModels: KernelBackend, MonolithBackend, update_sys_struct!,
    membrane_strain_energy, membrane_dissipation, membrane_corner_forces
using SymbolicAWEModels.ForwardDiff
using SymbolicAWEModels.Symbolics
using KiteUtils
using LinearAlgebra

const FABRIC = (; youngs_modulus = 5e8, poisson_ratio = 0.3, thickness = 3e-4,
                damping = 2e3)

"""The energy [J] of a triangle of `FABRIC` with corners at the columns of `corners`,
the rest lengths `rest_lengths` and keeping `slack_frac` of its taut energy."""
triangle_energy(corners, rest_lengths; slack_frac = 1e-4) =
    membrane_strain_energy(corners[:, 1], corners[:, 2], corners[:, 3], rest_lengths,
        FABRIC.youngs_modulus, FABRIC.poisson_ratio, FABRIC.thickness, slack_frac)

"""The rest lengths of the edges 1→2, 2→3 and 3→1 of the triangle `corners`."""
edge_lengths(corners) = [norm(corners[:, mod1(edge + 1, 3)] - corners[:, edge])
                         for edge in 1:3]

"""The right triangle of unit legs, in the plane z = 0, one corner per column."""
const UNIT_TRIANGLE = [0.0 1.0 0.0; 0.0 0.0 1.0; 0.0 0.0 0.0]

"""`corners` stretched by `sx` along x and `sy` along y."""
stretched(corners, sx, sy) = Diagonal([sx, sy, 1.0]) * corners

@testset "membrane triangles" begin
    @testset "corner forces are minus the gradient of the energy and dissipation" begin
        rest_lengths = edge_lengths(UNIT_TRIANGLE)
        params = (; membranes = [(; FABRIC..., rest_lengths, slack_frac = 1e-4)])
        tilt = 0.01 .* [0.0 0.1 -0.1; 0.0 0.2 0.03; 0.0 0.05 0.02]
        velocities = [0.3 -0.2 0.1; 0.0 0.4 -0.1; 0.2 0.0 -0.3]
        @testset "$state" for (state, sx, sy) in (("taut", 1.05, 1.02),
                                                  ("wrinkled", 1.05, 0.9),
                                                  ("slack", 0.8, 0.8))
            corners = stretched(UNIT_TRIANGLE, sx, sy) .+ tilt
            symbolic = membrane_corner_forces(params, 1, collect.(eachcol(corners)),
                                              collect.(eachcol(velocities)))
            force = Symbolics.value.(substitute.(symbolic, (Dict(),); fold = Val(true)))
            energy(x) = triangle_energy(reshape(x, 3, 3), rest_lengths)
            dissipation(v) = membrane_dissipation(corners[:, 1], corners[:, 2],
                corners[:, 3], v[1:3], v[4:6], v[7:9], rest_lengths,
                FABRIC.thickness, FABRIC.damping)
            expected = .-(ForwardDiff.gradient(energy, vec(corners)) .+
                          ForwardDiff.gradient(dissipation, vec(velocities)))
            @test vec(force) ≈ expected rtol = 1e-10
        end
    end

    @testset "a wrinkled triangle stores only the energy of its tension" begin
        rest_lengths = edge_lengths(UNIT_TRIANGLE)
        area = 0.5
        e1 = (1.05^2 - 1) / 2
        uniaxial = FABRIC.youngs_modulus * e1^2 / 2 * area * FABRIC.thickness
        for sy in (0.9, 0.7)
            corners = stretched(UNIT_TRIANGLE, 1.05, sy)
            @test triangle_energy(corners, rest_lengths; slack_frac = 0) ≈ uniaxial rtol =
                1e-6
        end
        @test triangle_energy(stretched(UNIT_TRIANGLE, 0.8, 0.8), rest_lengths;
                              slack_frac = 0) == 0
        taut = stretched(UNIT_TRIANGLE, 1.05, 1.02)
        @test triangle_energy(taut, rest_lengths; slack_frac = 0) ≈
              triangle_energy(taut, rest_lengths; slack_frac = 1)
    end

    @testset "a flat sheet compressed in plane carries no compressive stress" begin
        nodes = [[i, j, 0.0] for j in 0:2 for i in 0:2]
        node(i, j) = 3j + i + 1
        triangles = [t for j in 0:1, i in 0:1
                     for t in ((node(i, j), node(i + 1, j), node(i + 1, j + 1)),
                               (node(i, j), node(i + 1, j + 1), node(i, j + 1)))]
        rest = reduce(hcat, nodes)
        rest_lengths = [edge_lengths(rest[:, collect(t)]) for t in triangles]
        sheet_energy(x, slack_frac) = sum(zip(triangles, rest_lengths)) do (t, lengths)
            triangle_energy(reshape(x, 3, :)[:, collect(t)], lengths; slack_frac)
        end
        deformed = vec(stretched(rest, 1.05, 0.9))
        top, right = node.(0:2, 2), node.(2, 0:2)
        # Force the sheet puts on its top and right edges, summed along each.
        function edge_loads(slack_frac)
            force = reshape(-ForwardDiff.gradient(x -> sheet_energy(x, slack_frac),
                                                  deformed), 3, :)
            return sum(force[2, top]), sum(force[1, right])
        end
        across, along = edge_loads(0.0)
        e1 = (1.05^2 - 1) / 2
        height = 2.0
        tension = 1.05 * FABRIC.youngs_modulus * e1 * FABRIC.thickness * height
        @test -along ≈ tension rtol = 1e-6
        @test abs(across) < 1e-6 * abs(along)
        svk_across, _ = edge_loads(1.0)
        @test svk_across > 1e-2 * abs(along)
    end
end

"""Canopy faces of the 2plate kite, each quad split into two membrane triangles."""
const MEMBRANE_CANOPY = """
canopy_faces:
  headers: [name, wing, points, model, youngs_modulus, poisson_ratio, thickness, damping]
  data:
    - [canopy_1a, main_wing, [le_left, te_left, te_center], sam_membrane, 5.0e8, 0.3, 3.0e-4, 1.0e5]
    - [canopy_1b, main_wing, [le_left, te_center, le_center], sam_membrane, 5.0e8, 0.3, 3.0e-4, 1.0e5]
    - [canopy_2a, main_wing, [le_center, te_center, te_right], sam_membrane, 5.0e8, 0.3, 3.0e-4, 1.0e5]
    - [canopy_2b, main_wing, [le_center, te_right, le_right], sam_membrane, 5.0e8, 0.3, 3.0e-4, 1.0e5]
"""

"""
    membrane_geometry() -> String

Path of a copy of the 2plate kite's particle geometry whose canopy is
[`MEMBRANE_CANOPY`](@ref).
"""
function membrane_geometry()
    source = read(joinpath(dirname(@__DIR__), "data", "2plate_kite",
                           "particle_structural_geometry.yaml"), String)
    path = joinpath(mktempdir(), "membrane_structural_geometry.yaml")
    write(path, replace(source, r"canopy_faces:\n(  .*\n)*" => MEMBRANE_CANOPY))
    return path
end

@testset "a membrane canopy" begin
    # `init!` writes the assembled model under the data path; an earlier test file may
    # have left it on a directory it has since removed.
    data_path = joinpath(mktempdir(), "2plate_kite")
    cp(joinpath(dirname(@__DIR__), "data", "2plate_kite"), data_path)
    set_data_path(data_path)

    @testset "a canopy face without the membrane model adds no membrane" begin
        sys, _, _ = load_2plate_sys("2plate_kite")
        @test isempty(sys.membranes)
    end

    sys, set, transforms = load_2plate_sys("2plate_membrane";
                                           geometry = membrane_geometry())
    @testset "loads from the canopy faces with its rest shape" begin
        @test length(sys.membranes) == 4
        membrane = sys.membranes[:canopy_1a]
        @test membrane.point_idxs == Tuple(sys.points[name].idx
                                           for name in (:le_left, :te_left, :te_center))
        @test membrane.rest_lengths ≈ [norm(sys.points[:te_left].pos_ENU -
                                            sys.points[:le_left].pos_ENU),
                                       norm(sys.points[:te_center].pos_ENU -
                                            sys.points[:te_left].pos_ENU),
                                       norm(sys.points[:le_left].pos_ENU -
                                            sys.points[:te_center].pos_ENU)]
    end

    for membrane in sys.membranes
        membrane.rest_lengths .*= 0.98
    end
    monolith = SymbolicAWEModel(set, sys; transforms, backend = MonolithBackend())
    kernel = SymbolicAWEModel(set, sys; transforms, backend = KernelBackend())
    @testset "both backends put the same membrane force on the canopy" begin
        forces = map((monolith, kernel)) do sam
            init!(sam; prn = false)
            update_sys_struct!(sam.prob, sam.integrator, sam.state)
            reduce(hcat, [point.force for point in sam.state.points])
        end
        @test forces[1] ≈ forces[2] rtol = 1e-8
    end

    @testset "find_equilibrium! finds where a long damped run of it comes to rest" begin
        start_at!(kernel; elevation = -90.0, v_wind = 0.0)
        start_state = copy(kernel.integrator.u)
        result = find_equilibrium!(kernel)
        @test result.converged
        equilibrium = point_positions(kernel)
        SymbolicAWEModels.settle!(kernel, start_state)
        set_world_frame_damping(kernel.state.points, 2.0)
        for _ in 1:6000
            next_step!(kernel; dt = 0.05)
        end
        @test largest_distance(point_positions(kernel), equilibrium) < 2e-3
        spring_net = kernel_2plate("particle_structural_geometry.yaml")
        start_at!(spring_net; elevation = -90.0, v_wind = 0.0)
        find_equilibrium!(spring_net)
        @test largest_distance(point_positions(spring_net), equilibrium) > 1e-3
    end
end
