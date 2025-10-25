# SPDX-FileCopyrightText: 2025 Uwe Fechner, Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

using Makie
using SymbolicAWEModels
using Documenter

DocMeta.setdocmeta!(SymbolicAWEModels, :DocTestSetup, :(using SymbolicAWEModels); recursive=true)

makedocs(;
    modules=[SymbolicAWEModels],
    authors="Uwe Fechner <fechner@aenarete.eu>, Bart van de Lint <bart@vandelint.net> and contributors",
    repo="https://github.com/OpenSourceAWE/SymbolicAWEModels.jl/blob/{commit}{path}#{line}",
    sitename="SymbolicAWEModels.jl",
    format=Documenter.HTML(;
        repolink = "https://github.com/OpenSourceAWE/SymbolicAWEModels.jl",
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://OpenSourceAWE.github.io/SymbolicAWEModels.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Getting Started" => "getting_started.md",
        "Examples" => "examples.md",
        "Custom Model" => "tutorial_system_structure.md",
        "API Reference" => [
            "SymbolicAWEModel" => "symbolic_awe_model.md",
            "System Structure" => "system_structure.md",
            "State Management" => "state_management.md",
            "Simulation Functions" => "exported_functions.md",
        ],
        "Parameters" => "parameters.md",
        "Private functions" => "private_functions.md",
        "Developers" => "developers.md",
    ],
)

deploydocs(;
    repo="github.com/OpenSourceAWE/SymbolicAWEModels.jl",
    devbranch="main",
)
