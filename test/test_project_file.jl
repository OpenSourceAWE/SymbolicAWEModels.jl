# SPDX-FileCopyrightText: 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_project_file.jl - Test `project_file`, which resolves the entries of a
# project file (`system.yaml`) against the data path.
#
# Verifies:
# 1. Every entry of a shipped model's project file resolves to a file that exists
# 2. An entry the project does not name, a data path with no project file and an
#    empty `system:` block all throw an error naming the missing entry, where
#    `optional_project_file` gives `nothing`
# 3. A project file in a subdirectory of the data path resolves against that
#    subdirectory

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

using Test
using SymbolicAWEModels
using KiteUtils

pkg_root = dirname(@__DIR__)
previous_data_path = get_data_path()

@testset "every entry of a model's project file resolves to a file" begin
    set_data_path(joinpath(pkg_root, "data", "2plate_kite"))
    for entry in ("sim_settings", "structural_geometry", "aero_geometry",
                  "vsm_settings")
        @test isfile(project_file(entry))
    end
end

@testset "an entry the project does not name throws, naming the entry" begin
    set_data_path(joinpath(pkg_root, "data", "beam"))
    @test isfile(project_file("sim_settings"))
    @test_throws r"names no system.aero_geometry" project_file("aero_geometry")
    @test isnothing(SymbolicAWEModels.optional_project_file("vsm_settings"))
end

@testset "a data path without a project file throws" begin
    set_data_path(mktempdir())
    @test_throws ArgumentError project_file("sim_settings")
    @test isnothing(SymbolicAWEModels.optional_project_file("sim_settings"))
end

@testset "a project file with an empty system block names nothing" begin
    set_data_path(mktempdir())
    write(joinpath(get_data_path(), "system.yaml"), "system:\n")
    @test isnothing(SymbolicAWEModels.optional_project_file("vsm_settings"))
    @test_throws ArgumentError project_file("vsm_settings")
end

@testset "a project in a subdirectory resolves against that subdirectory" begin
    set_data_path(joinpath(pkg_root, "data"))
    @test project_file("vsm_settings", "base/system.yaml") ==
          joinpath(pkg_root, "data", "base", "vsm_settings.yaml")
end

set_data_path(previous_data_path)
