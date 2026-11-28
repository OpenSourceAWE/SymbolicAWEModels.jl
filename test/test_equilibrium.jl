# SPDX-FileCopyrightText: 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_equilibrium.jl - find_equilibrium!, the structure's equilibrium at a held pose

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

@isdefined(test_init!) || include(joinpath(@__DIR__, "util.jl"))

using Test
using SymbolicAWEModels
using SymbolicAWEModels: KernelBackend, EquilibriumProblem, generalized_force, at_rest
using KiteUtils
using LinearAlgebra

"""
    stretch!(sam, factor)

Move every free point and body of `sam` radially away from the base by `factor`, a
start off the equilibrium that leaves the pose unturned.
"""
function stretch!(sam, factor)
    problem = EquilibriumProblem(sam)
    u = copy(sam.integrator.u)
    positions = reshape(view(u, problem.layout.position), 3, :)
    for position in eachcol(positions)
        position .= problem.base .+ factor .* (position .- problem.base)
    end
    SymbolicAWEModels.settle!(sam, u)
    return sam
end

# `init!` writes the assembled model under the data path; an earlier test file may
# have left it on a directory it has since removed.
data_path = joinpath(mktempdir(), "2plate_kite")
cp(joinpath(dirname(@__DIR__), "data", "2plate_kite"), data_path)
set_data_path(data_path)

particle = kernel_2plate("particle_structural_geometry.yaml")
rigid = kernel_2plate("rigid_structural_geometry.yaml")

@testset "find_equilibrium!" begin
    @testset "a hanging kite settles where a long damped run comes to rest" begin
        start_at!(particle; elevation = -90.0, v_wind = 0.0)
        start = point_positions(particle)
        start_state = copy(particle.integrator.u)
        result = find_equilibrium!(particle)
        @test result.converged
        equilibrium = point_positions(particle)
        SymbolicAWEModels.settle!(particle, start_state)
        set_world_frame_damping(particle.state.points, 2.0)
        for _ in 1:6000
            next_step!(particle; dt = 0.05)
        end
        set_world_frame_damping(particle.state.points, 0.0)
        @test largest_distance(point_positions(particle), equilibrium) < 2e-3
        @test largest_distance(equilibrium, start) > 0.1
    end

    @testset "a $name kite in wind returns from a stretched start" for (name, sam) in
            (("particle", particle), ("rigid", rigid))
        start_at!(sam; elevation = 50.0, v_wind = 15.51)
        start = copy(sam.integrator.u)
        @test find_equilibrium!(sam).converged
        plain = point_positions(sam)
        SymbolicAWEModels.settle!(sam, start)
        stretch!(sam, 1.01)
        result = find_equilibrium!(sam)
        @test result.converged
        @test largest_distance(point_positions(sam), plain) < 1e-3
        problem = EquilibriumProblem(sam)
        @test all(iszero, sam.integrator.u[problem.rows])
    end

    @testset "damping does not move the equilibrium" begin
        start_at!(particle; elevation = 50.0, v_wind = 15.51)
        start = copy(particle.integrator.u)
        find_equilibrium!(particle)
        undamped = point_positions(particle)
        SymbolicAWEModels.settle!(particle, start)
        set_world_frame_damping(particle.state.points, 5.0)
        find_equilibrium!(particle)
        set_world_frame_damping(particle.state.points, 0.0)
        @test largest_distance(point_positions(particle), undamped) < 1e-9
    end

    @testset "a pose off its equilibrium reports the moment holding it" begin
        start_at!(particle; elevation = -60.0, v_wind = 0.0)
        result = find_equilibrium!(particle)
        @test result.converged
        problem = EquilibriumProblem(particle)
        g = particle.set.g_earth
        weight_moment = sum(particle.state.points) do point
            cross(point.pos_w .- problem.base, [0.0, 0.0, -point.total_mass * g])
        end
        @test norm(weight_moment) > 10
        @test result.moment ≈ weight_moment rtol = 1e-5
    end

    @testset "a search cut short says so with what is left" begin
        start_at!(rigid; elevation = 50.0, v_wind = 15.51)
        stretch!(rigid, 1.01)
        result = find_equilibrium!(rigid; max_iterations = 1)
        @test !result.converged
        @test result.iterations == 1
        @test isfinite(result.residual) && result.residual > 0
    end

    @testset "a model on the monolith backend is refused by name" begin
        sys, set, transforms = load_2plate_sys("2plate_kite")
        monolith = SymbolicAWEModel(set, sys; transforms,
                                    backend = SymbolicAWEModels.MonolithBackend())
        @test_throws "KernelBackend" find_equilibrium!(monolith)
    end
end
