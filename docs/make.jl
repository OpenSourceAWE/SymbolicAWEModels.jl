# SPDX-FileCopyrightText: 2025 Uwe Fechner, Bart van de Lint
# SPDX-License-Identifier: MIT

using Pkg
if ("TestEnv" ∈ keys(Pkg.project().dependencies))
    if ! ("Documents" ∈ keys(Pkg.project().dependencies))
        using TestEnv; TestEnv.activate()
    end
end
using ControlPlots, VortexStepMethod
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
        "Types" => "types.md",
        "Functions" => "functions.md",
        "SymbolicAWEModel" => [
		    "Overview" => "symbolic_awe_model/ram_air_kite.md",
		    "Examples" => "symbolic_awe_model/examples_ram_air.md",
		    "Tutorial" => "symbolic_awe_model/tutorial_system_structure.md"
	    ],
        "Parameters" => "parameters.md",
        "Quickstart" => "quickstart.md",
        "Advanced usage" => "advanced.md",
    ],
)

deploydocs(;
    repo="github.com/OpenSourceAWE/SymbolicAWEModels.jl",
    devbranch="main",
)
