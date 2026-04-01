# Copyright (c) 2025 Bart van de Lint, Jelle Poland
# SPDX-License-Identifier: MPL-2.0

# build and display the html documentation locally
# run with: julia --project=docs scripts/build_docu.jl

using Pkg
docs_project = joinpath(dirname(@__DIR__), "docs")
if Pkg.project().path != joinpath(docs_project, "Project.toml")
    Pkg.activate(docs_project)
end
using LiveServer; servedocs(launch_browser=true)
