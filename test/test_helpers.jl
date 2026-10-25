# SPDX-FileCopyrightText: 2022, 2024, 2025 Uwe Fechner
# SPDX-License-Identifier: MIT

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

using Test
using SymbolicAWEModels

@testset "copy_data and copy_examples honour force" begin
    path=pwd()
    tmpdir=mktempdir()
    mkpath(tmpdir)
    cd(tmpdir)
    SymbolicAWEModels.copy_data(; force=true)
    SymbolicAWEModels.copy_examples(; force=true)
    @test isfile(joinpath(tmpdir, "examples", "menu.jl"))
    if ! Sys.iswindows()
        rm(tmpdir, recursive=true)
    end
    cd(path)

    path=pwd()
    tmpdir=mktempdir()
    mkpath(tmpdir)
    menu_file = joinpath(tmpdir, "examples", "menu.jl")
    mkpath(joinpath(tmpdir, "examples"))
    touch(menu_file)
    cd(tmpdir)
    SymbolicAWEModels.copy_data(; force=false)
    SymbolicAWEModels.copy_examples(; force=false)
    @test isfile(joinpath(tmpdir, "examples", "menu.jl"))
    @test filesize(menu_file) == 0
    if ! Sys.iswindows()
        rm(tmpdir, recursive=true)
    end
    cd(path)
end

@testset "loading the package points the output path at ./output" begin
    output_path_before = get_output_path()
    set_output_path()
    cd(mktempdir()) do
        SymbolicAWEModels.__init__()
        @test get_output_path() == "output"
        @test isdir("output")
    end
    set_output_path(output_path_before)
end

@testset "test environment has no dev-only packages or Manifest.toml" begin
    @test ! ("TestEnv" ∈ keys(Pkg.project().dependencies))
    @test ! ("Revise" ∈ keys(Pkg.project().dependencies))
    @test ! ("Plots" ∈ keys(Pkg.project().dependencies))

    @test ! isfile(joinpath(dirname(@__DIR__), "Manifest.toml"))
end
nothing
