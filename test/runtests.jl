# SPDX-FileCopyrightText: 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

ENV["MPLBACKEND"] = "Agg"
using KiteUtils
using Test
using SymbolicAWEModels

include("util.jl")

# Set up data path for 2plate_kite tests
pkg_root = dirname(@__DIR__)
src_data_path = joinpath(pkg_root, "data", "2plate_kite")
tmpdir = mktempdir()
data_path = joinpath(tmpdir, "2plate_kite")
cp(src_data_path, data_path; force=true)
@show data_path
set_data_path(data_path)

"""
    KERNEL_PARITY_FILES

The test files the [`KernelBackend`](@ref) currently covers, run instead of the
whole suite when `SYMAWE_TEST_BACKEND=kernel`. It grows as the backend does; the
files it leaves out are listed at startup.
"""
KERNEL_PARITY_FILES = [
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

backend_name = lowercase(strip(get(ENV, "SYMAWE_TEST_BACKEND", "")))
exclude = ["test_for_precompile.jl", "test_menu.jl"]
all_files = sort(filter(readdir(@__DIR__)) do f
    startswith(f, "test_") && endswith(f, ".jl") && f ∉ exclude
end)

if backend_name in ("", "monolith")
    test_files = all_files
elseif backend_name == "kernel"
    default_backend!(KernelBackend())
    test_files = sort(KERNEL_PARITY_FILES)
    unknown = filter(!in(all_files), test_files)
    isempty(unknown) ||
        error("KERNEL_PARITY_FILES names files the suite does not have: $unknown")
else
    error("Unknown SYMAWE_TEST_BACKEND=\"$backend_name\"; " *
          "expected \"monolith\" or \"kernel\".")
end
println("Running the test suite on the $(nameof(typeof(default_backend()))).")

skipped = filter(!in(test_files), all_files)
if !isempty(skipped)
    println("Skipping $(length(skipped)) file(s) this backend does not cover yet:")
    foreach(f -> println("    $f"), skipped)
end

# Filter by test_args if provided, e.g.:
#   Pkg.test(; test_args=["test_bench", "test_wing"])
if !isempty(ARGS)
    test_files = filter(test_files) do f
        name = replace(f, ".jl" => "")
        any(arg -> name == arg, ARGS)
    end
end

@testset verbose = true "Testing SymbolicAWEModels..." begin
    for f in test_files
        println("--> $f")
        include(f)
    end
end
nothing
