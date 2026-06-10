# SPDX-FileCopyrightText: 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_weighted_ref_points.jl - Unit tests for the
# WeightedRefPoints constructors and weight validation.
# These are pure (no model compilation) and exercise the
# input-coercion paths used by the YAML loader and the
# programmatic Wing constructors.

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

using Test
using SymbolicAWEModels
using SymbolicAWEModels: WeightedRefPoints, validate_weights!

@testset "WeightedRefPoints constructors" begin
    @testset "single Symbol ref" begin
        rp = WeightedRefPoints(:le)
        @test rp.refs == [:le]
        @test isempty(rp.ids)
        @test rp.weights == [1.0]
    end

    @testset "single String ref" begin
        rp = WeightedRefPoints("le")
        @test rp.refs == [:le]
        @test isempty(rp.ids)
        @test rp.weights == [1.0]
    end

    @testset "single resolved index" begin
        rp = WeightedRefPoints(7)
        @test isempty(rp.refs)
        @test rp.ids == [7]
        @test rp.weights == [1.0]
    end

    @testset "equal-weight average" begin
        rp = WeightedRefPoints([:le, :te])
        @test rp.refs == [:le, :te]
        @test rp.weights ≈ [0.5, 0.5]
    end

    @testset "explicit weighted tuples" begin
        rp = WeightedRefPoints([(:le, 0.7), (:te, 0.3)])
        @test rp.refs == [:le, :te]
        @test rp.weights ≈ [0.7, 0.3]
    end

    @testset "integer refs via tuples" begin
        rp = WeightedRefPoints([(2, 0.25), (4, 0.75)])
        @test rp.refs == [2, 4]
        @test rp.weights ≈ [0.25, 0.75]
    end

    @testset "weights normalized when sum != 1" begin
        rp = @test_logs (:warn,) WeightedRefPoints(
            [(:le, 2.0), (:te, 2.0)])
        @test rp.weights ≈ [0.5, 0.5]
    end

    @testset "non-positive weight sum errors" begin
        @test_throws ErrorException WeightedRefPoints(
            [(:le, 0.0), (:te, 0.0)])
    end

    @testset "empty vector errors" begin
        @test_throws ErrorException WeightedRefPoints(Symbol[])
    end

    @testset "deleted identity passthrough" begin
        rp = WeightedRefPoints(:le)
        @test_throws MethodError WeightedRefPoints(rp)
    end
end

@testset "validate_weights! normalization" begin
    w = [3.0, 1.0]
    @test_logs (:warn,) validate_weights!(w)
    @test w ≈ [0.75, 0.25]

    w_ok = [0.4, 0.6]
    @test_logs validate_weights!(w_ok)
    @test w_ok ≈ [0.4, 0.6]

    @test_throws ErrorException validate_weights!([0.0, 0.0])
end
nothing
