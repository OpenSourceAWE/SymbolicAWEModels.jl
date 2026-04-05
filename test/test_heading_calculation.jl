# Copyright (c) 2025 Bart van de Lint, Jelle Poland
# SPDX-License-Identifier: MPL-2.0

#!/usr/bin/env julia
"""
Test tangential sphere heading calculation.

Verifies that `calc_heading` correctly projects the body x-axis
onto the tangent plane of the tether sphere and returns the
angle from the elevation direction (x_t) toward azimuthal (y_t).

Tests run with both origin and non-origin base positions to
verify heading is computed relative to the sphere center.

Usage:
  julia --project=test test/test_heading_calculation.jl
  julia --project=test test/test_heading_calculation.jl special
"""

using Test
using LinearAlgebra
using SymbolicAWEModels: sym_normalize, sym_calc_R_t_to_w,
    calc_heading, calc_R_t_to_w

# Support selective test execution via command-line args
const test_patterns = isempty(ARGS) ? String[] : ARGS

println("Running heading calculation tests...")
if !isempty(test_patterns)
    println("Filtering tests matching: ", test_patterns)
end

function should_run_test(test_name::String)
    isempty(test_patterns) && return true
    return any(p -> occursin(lowercase(p),
        lowercase(test_name)), test_patterns)
end

function test_circular_path(base_pos, label)
    @testset "Circular Path ($label)" begin
        elevation = deg2rad(45)
        radius = 100.0
        circle_radius = 20.0
        center = base_pos + [radius * cos(elevation), 0.0,
                  radius * sin(elevation)]

        n_points = 36
        angles = range(0, 2π, length=n_points + 1)[1:end-1]

        for angle in angles
            wing_pos = center + circle_radius *
                [cos(angle), 0.0, sin(angle)]
            rel_pos = wing_pos - base_pos
            R_t_to_w = sym_calc_R_t_to_w(rel_pos)

            tangent = [-sin(angle), 0.0, cos(angle)]
            e_x = sym_normalize(tangent)

            e_x_t = R_t_to_w' * e_x
            expected = atan(e_x_t[2], e_x_t[1])

            R_b_to_w = hcat(e_x, zeros(3), zeros(3))
            @test calc_heading(R_b_to_w, rel_pos) ≈
                expected atol = 1e-10
        end
    end
end

function test_horizontal_circle(base_pos, label)
    @testset "Horizontal Circle ($label)" begin
        elevation = deg2rad(45)
        radius = 100.0
        n_points = 24
        azimuths = range(0, 2π,
            length=n_points + 1)[1:end-1]

        for azimuth in azimuths
            wing_pos = base_pos + radius * [
                cos(elevation) * cos(azimuth),
                cos(elevation) * sin(azimuth),
                sin(elevation)]
            rel_pos = wing_pos - base_pos

            tangent = [-sin(azimuth), cos(azimuth), 0.0]
            e_x = sym_normalize(tangent)

            R_b_to_w = hcat(e_x, zeros(3), zeros(3))
            heading_val = calc_heading(R_b_to_w, rel_pos)

            @test abs(heading_val) > π / 4
        end
    end
end

function test_special_cases(base_pos, label)
    @testset "Special Cases ($label)" begin
        rel_pos = [70.0, 70.0, 70.0]
        R_t_to_w = calc_R_t_to_w(rel_pos)

        # e_x along x_t (elevation dir) → heading = 0
        e_x = R_t_to_w[:, 1]
        R_b_to_w = hcat(e_x, zeros(3), zeros(3))
        @test calc_heading(R_b_to_w, rel_pos) ≈
            0.0 atol=1e-10

        # e_x along y_t (azimuthal) → heading = π/2
        e_x = R_t_to_w[:, 2]
        R_b_to_w = hcat(e_x, zeros(3), zeros(3))
        @test calc_heading(R_b_to_w, rel_pos) ≈
            π/2 atol=1e-10

        # Roundtrip: construct e_x at known heading angles
        for h in [0, π/6, π/4, π/3, π/2, 2π/3, π,
                  -π/6, -π/4]
            e_x_w = R_t_to_w * [cos(h), sin(h), 0.0]
            R_b_to_w = hcat(e_x_w, zeros(3), zeros(3))
            @test calc_heading(R_b_to_w, rel_pos) ≈
                h atol=1e-10
        end
    end
end

base_positions = [
    (zeros(3), "origin base"),
    ([10.0, -5.0, 3.0], "non-origin base"),
]

if should_run_test("circular")
@testset "Heading - Circular Path" begin
    for (bp, label) in base_positions
        test_circular_path(bp, label)
    end
end
end

if should_run_test("horizontal")
@testset "Heading - Horizontal Circle" begin
    for (bp, label) in base_positions
        test_horizontal_circle(bp, label)
    end
end
end

if should_run_test("special")
@testset "Heading - Special Cases" begin
    for (bp, label) in base_positions
        test_special_cases(bp, label)
    end
end
end

println("\n=== All Heading Tests Passed ===\n")

nothing
