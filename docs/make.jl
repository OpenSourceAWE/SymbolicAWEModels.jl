# SPDX-FileCopyrightText: 2025 Uwe Fechner, Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

using Pkg
Pkg.activate(@__DIR__)

using Makie, VortexStepMethod
using SymbolicAWEModels
using Documenter
using Literate

# --- Convert Literate .jl sources to .md ---
# Write to a temp dir first, then copy only if content changed.
# This avoids retriggering LiveServer's file watcher on every build.
literate_dir = joinpath(@__DIR__, "src", "literate")
output_dir = joinpath(@__DIR__, "src")
tmp_dir = mktempdir()

for file in readdir(literate_dir)
    endswith(file, ".jl") || continue
    edit_url = joinpath("literate", file)
    Literate.markdown(joinpath(literate_dir, file), tmp_dir;
        execute=false, documenter=true,
        codefence="```julia" => "```",
        postprocess=content -> begin
            content = replace(content,
                r"^[^\n]*#hide *\n"m => "")
            replace(content,
                r"^EditURL = \"[^\"]*\"$"m =>
                "EditURL = \"$edit_url\"")
        end)
    md_name = replace(file, ".jl" => ".md")
    src = joinpath(tmp_dir, md_name)
    dst = joinpath(output_dir, md_name)
    if !isfile(dst) || read(src) != read(dst)
        cp(src, dst; force=true)
    end
end

DocMeta.setdocmeta!(SymbolicAWEModels, :DocTestSetup, :(using SymbolicAWEModels); recursive=true)

makedocs(;
    modules=[SymbolicAWEModels],
    authors="Uwe Fechner <fechner@aenarete.eu>, Bart van de Lint <bart@vandelint.net> and contributors",
    repo="https://github.com/OpenSourceAWE/SymbolicAWEModels.jl/blob/{commit}{path}#{line}",
    sitename="SymbolicAWEModels.jl",
    warnonly=[:cross_references],
    format=Documenter.HTML(;
        repolink = "https://github.com/OpenSourceAWE/SymbolicAWEModels.jl",
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://OpenSourceAWE.github.io/SymbolicAWEModels.jl",
        assets=String[],
        size_threshold=300 * 1024,
    ),
    pages=[
        "Home" => "index.md",
        "Getting started" => "getting_started.md",
        "Building a system using Julia" => "tutorial_julia.md",
        "Building a system using YAML" => "tutorial_yaml.md",
        "Compilation pipeline" => "pipeline.md",
        "Examples" => "examples.md",
        "VSM coupling" => "vsm_coupling.md",
        "Coordinate frames" => "coordinate_frames.md",
        "Types" => "exported_types.md",
        "Functions" => "exported_functions.md",
        "Parameters" => "parameters.md",
        "Developer guide" => "developers.md",
        "Private API" => "private_functions.md",
    ],
)

deploydocs(;
    repo="github.com/OpenSourceAWE/SymbolicAWEModels.jl",
    devbranch="main",
)
