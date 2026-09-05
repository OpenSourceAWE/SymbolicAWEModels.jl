# Copyright (c) 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_vsm_solve_failure.jl — a VSM solve that does not converge, on the particle
# 2plate kite. Without `vsm_warn_on_fail` `next_step!` throws; with it, it warns,
# keeps the circulation and the point forces of the last converged solve, holds
# the `vsm_interval` schedule, and solves again as soon as the solve recovers.
#
# Non-convergence is forced through the solver rather than through a pose: the
# LOOP solver converges on `normalized_error < rtol`, so `rtol = 0` never
# converges and `max_iterations` bounds what that costs.

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

@isdefined(test_init!) || include(joinpath(@__DIR__, "util.jl"))

using Test
using SymbolicAWEModels
using SymbolicAWEModels: VortexStepMethod, VSMSolveFailure, warn_or_rethrow,
    wing_points
using KiteUtils
using Logging

"""
    break_solve!(solver)

Make every solve of `solver` fail: no relative tolerance is ever met, and the
iteration is cut short so the failure is cheap.
"""
function break_solve!(solver)
    solver.rtol = 0.0
    solver.max_iterations = 2
    return nothing
end

"""
    caught_by_refresh(failure, vsm_warn_on_fail)

Hand `failure` to [`warn_or_rethrow`](@ref) where `refresh_aero!` does: from
inside a `catch` block, the only place `rethrow` is allowed.
"""
function caught_by_refresh(failure, vsm_warn_on_fail)
    try
        throw(failure)
    catch caught
        warn_or_rethrow(caught, vsm_warn_on_fail)
    end
end

"""
    vsm_warnings(records)

The warnings a step emitted about a failed VSM solve.
"""
vsm_warnings(records) = filter(records) do record
    record.level == Logging.Warn && occursin("VSM solve failed", record.message)
end

@testset "VSM solve failure" begin
    @testset "only a failed solve is downgraded" begin
        @test_throws VSMSolveFailure caught_by_refresh(
            VSMSolveFailure("VSM solve failed."), false)
        @test_logs (:warn,) caught_by_refresh(
            VSMSolveFailure("VSM solve failed."), true)
        # A corrupted frozen state is fatal with the kwarg set as well.
        @test_throws AssertionError caught_by_refresh(
            AssertionError("non-finite point force"), true)
    end

    pkg_root = dirname(@__DIR__)
    tmpdir = mktempdir()
    data_path = joinpath(tmpdir, "2plate_kite")
    cp(joinpath(pkg_root, "data", "2plate_kite"), data_path; force=true)
    set_data_path(data_path)
    set = Settings("system.yaml")
    vsm_set = VortexStepMethod.VSMSettings(
        joinpath(data_path, "vsm_settings.yaml"); data_prefix=false)

    sys = load_sys_struct_from_yaml(
        joinpath(data_path, "particle_structural_geometry.yaml");
        system_name="wing_test_REFINE", set, vsm_set)
    sam = SymbolicAWEModel(set, sys)
    test_init!(sam)

    wing = sys.wings[1]
    solver = wing.vsm_solver
    rtol = solver.rtol
    max_iterations = solver.max_iterations
    dt = 0.05

    next_step!(sam; dt, vsm_interval=1)
    gamma = copy(solver.sol.gamma_distribution)
    point_forces = [copy(point.aero_force_b) for point in wing_points(sys, wing)]

    @testset "throws without the kwarg" begin
        break_solve!(solver)
        @test_throws VSMSolveFailure next_step!(sam; dt, vsm_interval=1)
        @test solver.sol.gamma_distribution == gamma
    end

    @testset "warns and reuses the last converged solve" begin
        records = Test.TestLogger()
        with_logger(records) do
            next_step!(sam; dt, vsm_interval=1, vsm_warn_on_fail=true)
        end
        @test length(vsm_warnings(records.logs)) == 1
        # The failed solve moved the circulation; what is stored is the old one.
        @test solver.lr.gamma_new != gamma
        @test solver.sol.gamma_distribution == gamma
        @test [copy(point.aero_force_b)
               for point in wing_points(sys, wing)] == point_forces
    end

    @testset "vsm_interval keeps its schedule while solves fail" begin
        records = Test.TestLogger()
        with_logger(records) do
            for _ in 1:4
                next_step!(sam; dt, vsm_interval=2, vsm_warn_on_fail=true)
            end
        end
        @test length(vsm_warnings(records.logs)) == 2
        @test solver.sol.gamma_distribution == gamma
    end

    @testset "solves again once the solver recovers" begin
        solver.rtol = rtol
        solver.max_iterations = max_iterations
        records = Test.TestLogger()
        with_logger(records) do
            next_step!(sam; dt, vsm_interval=1, vsm_warn_on_fail=true)
        end
        @test isempty(vsm_warnings(records.logs))
        @test solver.sol.gamma_distribution != gamma
        @test all(isfinite, solver.sol.gamma_distribution)
    end
end
