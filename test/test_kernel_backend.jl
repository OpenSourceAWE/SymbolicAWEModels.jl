# SPDX-FileCopyrightText: 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Runs the monolith's own test files against the `KernelBackend`, with nothing
# swapped but the backend. `default_backend!` changes what `SymbolicAWEModel`
# builds when its constructor is not given one, so each file runs unmodified and
# the assertions stay the monolith's.

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

using Test
using SymbolicAWEModels

@isdefined(test_init!) || include(joinpath(@__DIR__, "util.jl"))

"""
    SCHEDULED_PARITY_FILES

The test files the `KernelBackend` currently covers. It grows as the backend
does; the goal is every file `runtests.jl` runs except the monolith-only ones
(linearization and control-function generation).
"""
SCHEDULED_PARITY_FILES = [
    "test_point.jl",
    "test_segment.jl",
    "test_segment_nonlinear.jl",
    "test_pulley.jl",
    "test_tether_init.jl",
    "test_tether_winch.jl",
    "test_rigid_body.jl",
    "test_joint.jl",
    "test_timoshenko_joint.jl",
    "test_flap_beam.jl",
    "test_wing.jl",
    "test_getter_allocations.jl",
    "test_heading_calculation.jl",
    "test_match_aero_sections.jl",
    "test_multi_section_group.jl",
    "test_principal_body_frame.jl",
    "test_quaternion_auto_groups.jl",
    "test_quaternion_conversions.jl",
    "test_section_alignment.jl",
    "test_static_twist.jl",
    "test_weighted_ref_points.jl",
    "test_wing_dynamics.jl",
    "test_yaml_bodies.jl",
    "test_transform.jl",
    "test_yaml_weighted_ref.jl",
]

@testset verbose = true "KernelBackend parity" begin
    previous = default_backend()
    default_backend!(KernelBackend())
    try
        for file in SCHEDULED_PARITY_FILES
            println("--> [kernel] $file")
            @testset "$file" begin
                include(joinpath(@__DIR__, file))
            end
        end
    finally
        default_backend!(previous)
    end
end
nothing
