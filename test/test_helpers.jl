# SPDX-FileCopyrightText: 2022, 2024, 2025 Uwe Fechner
# SPDX-License-Identifier: MIT

using Test
using SymbolicAWEModels
using Pkg

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
end
