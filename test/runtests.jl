# SPDX-FileCopyrightText: 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

ENV["MPLBACKEND"] = "Agg"
using KiteUtils
using Test
using SymbolicAWEModels

include("util.jl")
include("suite_groups.jl")

# Set up data path for 2plate_kite tests
pkg_root = dirname(@__DIR__)
src_data_path = joinpath(pkg_root, "data", "2plate_kite")
tmpdir = mktempdir()
data_path = joinpath(tmpdir, "2plate_kite")
cp(src_data_path, data_path; force=true)
@show data_path
set_data_path(data_path)

backend_name = test_backend_name()
backend_name == "kernel" && default_backend!(KernelBackend())

all_files = all_test_files()
skip_reasons = backend_skip_reasons(backend_name, all_files)
test_files = backend_test_files(backend_name, all_files)

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
    requested = replace.(ARGS, ".jl" => "")
    unknown = filter(!in(replace.(all_files, ".jl" => "")), requested)
    isempty(unknown) || error("No such test file(s): $unknown")
    test_files = filter(f -> replace(f, ".jl" => "") in requested, test_files)
end

@testset verbose = true "Testing SymbolicAWEModels..." begin
    for f in test_files
        println("--> $f")
        include(f)
    end
end
nothing
