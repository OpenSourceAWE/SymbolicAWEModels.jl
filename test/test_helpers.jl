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
    SymbolicAWEModels.copy_examples()
    @test isfile(joinpath(tmpdir, "examples", "menu.jl"))
    if ! Sys.iswindows()
        rm(tmpdir, recursive=true)
    end
    cd(path)
    path=pwd()
    tmpdir=mktempdir()
    mkpath(tmpdir)
    cd(tmpdir)
    SymbolicAWEModels.install_examples(false)
    @test isfile(joinpath(tmpdir, "examples", "menu.jl"))
    if ! Sys.iswindows()
        rm(tmpdir, recursive=true)
    end
    cd(path)

    @test ! ("TestEnv" ∈ keys(Pkg.project().dependencies))
    @test ! ("Revise" ∈ keys(Pkg.project().dependencies))
    @test ! ("Plots" ∈ keys(Pkg.project().dependencies))
    # # ensure that BenchmarkTools is not in the main environment
    # oldprpath = Pkg.project().path
    # if ! Pkg.project().ispackage
    #     Pkg.activate(".")
    # end
    # hasbm = ("BenchmarkTools" ∈ keys(Pkg.project().dependencies))
    # Pkg.activate(oldprpath)
    # @test ! hasbm
end
