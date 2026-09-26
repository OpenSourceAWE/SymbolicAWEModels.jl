# SPDX-FileCopyrightText: 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_makie_static_figure.jl - A SystemStructure drawn into a CairoMakie Axis3
# and saved as a vector PDF, the route docs/src/exported_functions.md documents.

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

@isdefined(test_init!) || include(joinpath(@__DIR__, "util.jl"))

using Test
using SymbolicAWEModels
import CairoMakie
import MakieControlPlots

@testset "plot! draws a system into a CairoMakie Axis3 that saves as a vector PDF" begin
    sys = load_2plate_particle_sys("2plate_segment_role")

    role_colors = Dict(:wing => :black, :bridle => :steelblue, :tether => :darkorange)
    fig = CairoMakie.Figure(size=(600, 500))
    ax = CairoMakie.Axis3(fig[1, 1]; aspect=:data,
        xlabel="x [m]", ylabel="y [m]", zlabel="z [m]")
    plots = CairoMakie.plot!(ax, sys;
        segment_color=segment -> role_colors[segment_role(sys, segment)],
        show_orient=false, show_wing_frame=false)
    CairoMakie.Legend(fig[1, 2],
        [CairoMakie.LineElement(color=color) for color in values(role_colors)],
        string.(keys(role_colors)))
    pdf_path = joinpath(mktempdir(), "kite.pdf")
    CairoMakie.save(pdf_path, fig; backend=CairoMakie)

    @test plots[:segments] in ax.scene.plots
    role_color(segment) = CairoMakie.to_color(role_colors[segment_role(sys, segment)])
    @test plots[:segment_colors_obs][] == role_color.(sys.segments)
    pdf = read(pdf_path, String)
    @test startswith(pdf, "%PDF")
    @test !occursin("/Image", pdf)
end
nothing
