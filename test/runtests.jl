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
    KERNEL_SKIP_FILES

Why each test file is left out when `SYMAWE_TEST_BACKEND=kernel`; every other
file in the suite runs. Entries are printed with their reason at startup, so a
gap in the [`KernelBackend`](@ref) has to be named here to be tolerated rather
than going unnoticed.
"""
KERNEL_SKIP_FILES = Dict(
    "test_aqua.jl" => "package-quality checks, backend-independent",
    "test_helpers.jl" => "covers the test helpers themselves, builds no model",
    "test_yaml_variables.jl" => "YAML loader only, builds no model",
    "test_bench.jl" => "benchmarks; every file here already asserts a " *
                       "zero-allocation RHS through test_init!",
    "test_backend_parity.jl" => "builds both backends itself, so the " *
                                "monolith run already covers it",
)

backend_name = lowercase(strip(get(ENV, "SYMAWE_TEST_BACKEND", "")))
exclude = ["test_for_precompile.jl", "test_menu.jl"]
all_files = sort(filter(readdir(@__DIR__)) do f
    startswith(f, "test_") && endswith(f, ".jl") && f ∉ exclude
end)

skip_reasons = Dict{String, String}()
if backend_name in ("", "monolith")
elseif backend_name == "kernel"
    default_backend!(KernelBackend())
    skip_reasons = KERNEL_SKIP_FILES
    unknown = sort(filter(!in(all_files), collect(keys(skip_reasons))))
    isempty(unknown) ||
        error("KERNEL_SKIP_FILES names files the suite does not have: $unknown")
else
    error("Unknown SYMAWE_TEST_BACKEND=\"$backend_name\"; " *
          "expected \"monolith\" or \"kernel\".")
end
test_files = filter(!in(keys(skip_reasons)), all_files)

println("Running the test suite on the $(nameof(typeof(default_backend()))).")
if !isempty(skip_reasons)
    println("Skipping $(length(skip_reasons)) of $(length(all_files)) file(s):")
    foreach(sort(collect(keys(skip_reasons)))) do f
        println("    $f — $(skip_reasons[f])")
    end
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
