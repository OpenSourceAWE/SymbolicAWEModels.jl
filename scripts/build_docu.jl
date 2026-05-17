# Copyright (c) 2025 Bart van de Lint, Jelle Poland
# SPDX-License-Identifier: LGPL-3.0-only

# build and display the html documentation locally
# run with: julia --project=docs scripts/build_docu.jl

using Pkg
docs_project = joinpath(dirname(@__DIR__), "docs")
if Pkg.project().path != joinpath(docs_project, "Project.toml")
    Pkg.activate(docs_project)
end

generated_figures = [
    "tether_sys_struct.png",
    "tether_sim.gif",
    "winch_sys_struct.png",
    "winch_sim.gif",
    "pulley_sys_struct.png",
    "pulley_sim.gif",
    "2plate_kite_structure.png",
]
assets_dir = joinpath(docs_project, "src", "assets")
missing_figures = filter(fig -> !isfile(joinpath(assets_dir, fig)), generated_figures)

if !isempty(missing_figures)
    println("Generating documentation figures: $(join(missing_figures, ", "))")
    include(joinpath(docs_project, "generate_figures.jl"))
end

using LiveServer; servedocs(launch_browser=true)
