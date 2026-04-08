# SPDX-FileCopyrightText: 2022, 2024, 2025 Uwe Fechner
# SPDX-License-Identifier: MIT

using Pkg
if ! ("PackageCompiler" ∈ keys(Pkg.project().dependencies))
    @info "Installing PackageCompiler ..."
    Pkg.add("PackageCompiler")
end
if ! ("GLMakie" ∈ keys(Pkg.project().dependencies))
    @info "Installing GLMakie ..."
    Pkg.add("GLMakie")
end
@info "Loading packages ..."
using KiteUtils, KitePodModels, SymbolicAWEModels, GLMakie
using PackageCompiler

@info "Creating sysimage ..."
push!(LOAD_PATH,joinpath(pwd(),"src"))

PackageCompiler.create_sysimage(
    [:KiteUtils, :KitePodModels, :SymbolicAWEModels, :GLMakie];
    sysimage_path="kps-image_tmp.so",
    include_transitive_dependencies=true,
    precompile_execution_file=joinpath("test", "test_for_precompile.jl")
)
nothing
