# SPDX-FileCopyrightText: 2022, 2024, 2025 Uwe Fechner
# SPDX-License-Identifier: MIT

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

using Test
using SymbolicAWEModels

@testset "Testing helper functions..." begin
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

    @test ! ("TestEnv" ∈ keys(Pkg.project().dependencies))
    @test ! ("Revise" ∈ keys(Pkg.project().dependencies))
    @test ! ("Plots" ∈ keys(Pkg.project().dependencies))

    root = dirname(@__DIR__)
    @test ! isfile(joinpath(root, "Manifest.toml"))
    project_mtime = mtime(joinpath(root, "Project.toml"))
    # Allow 1s tolerance: on CI, git checkout assigns timestamps in sequence order,
    # so Project.toml may be fractionally newer due to checkout ordering.
    @test mtime(joinpath(root, "Manifest-v1.11.toml.default")) >= project_mtime - 1.0
    @test mtime(joinpath(root, "Manifest-v1.12.toml.default")) >= project_mtime - 1.0
end

@testset "yaml_n_unrefined_sections" begin
    n_sections = SymbolicAWEModels.yaml_n_unrefined_sections
    @test n_sections((; n_unrefined_sections=4)) == 4
    @test n_sections((; twist_surfaces=[:a, :b, :c])) == 3
    # An explicit value wins over the twist-surface count.
    @test n_sections((; n_unrefined_sections=5, twist_surfaces=[:a])) == 5
    @test n_sections((; n_unrefined_sections=nothing, twist_surfaces=[:a])) == 1
    @test isnothing(n_sections((; name=:wing)))
    @test isnothing(n_sections((; twist_surfaces=nothing)))
end
nothing
