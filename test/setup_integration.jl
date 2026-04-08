# SPDX-FileCopyrightText: 2026 Bart van de Lint
#
# SPDX-License-Identifier: MPL-2.0

# Integration test for end-user setup.
# Run from the repo root: julia test/setup_integration.jl
#
# This runs everything in a single Julia session to avoid repeated
# startup/compilation overhead.
#
# README examples are extracted and executed directly from README.md
# so the test always stays in sync with the documentation.

using Test

const REPO_ROOT = dirname(@__DIR__)

# Create a temporary directory simulating a fresh user project
const USER_DIR = mktempdir()
println("  tmpdir: $USER_DIR")
cd(USER_DIR)

# Activate the temp project and dev-install SymbolicAWEModels
using Pkg
Pkg.activate(".")
Pkg.develop(path=REPO_ROOT)

using SymbolicAWEModels

"""
    extract_readme_code(heading_pattern) -> String

Extract the first Julia code block following a heading that
matches `heading_pattern` in README.md. Lines containing
`Pkg`, `pkg"`, `plot(`, and `plot!(` lines are filtered
out since the project is already activated and plotting
requires a display.
"""
function extract_readme_code(heading_pattern)
    readme = read(joinpath(REPO_ROOT, "README.md"), String)
    lines = split(readme, '\n')
    found_heading = false
    in_block = false
    code_lines = String[]
    skip_patterns = [r"using Pkg", r"^pkg\"",
        r"using GLMakie", r"^plot[!\(]", r"^\s+plot[!\(]",
        r"sleep\("]
    for line in lines
        if !found_heading
            if occursin(heading_pattern, line)
                found_heading = true
            end
            continue
        end
        if !in_block
            if startswith(line, "```julia")
                in_block = true
            end
            continue
        end
        if startswith(line, "```")
            break
        end
        if !any(p -> occursin(p, line), skip_patterns)
            push!(code_lines, line)
        end
    end
    isempty(code_lines) && error(
        "No code block found after heading " *
        "matching: $heading_pattern")
    return join(code_lines, '\n')
end

@testset "End-user setup" begin
    @testset "copy_data()" begin
        SymbolicAWEModels.copy_data()
        for d in ["data/2plate_kite", "data/base",
                  "data/saddle_form"]
            @test isdir(d)
        end
    end

    @testset "copy_examples()" begin
        SymbolicAWEModels.copy_examples()
        for f in ["menu.jl", "hanging_mass.jl",
                   "catenary_line.jl", "pulley.jl",
                   "saddle_form.jl", "coupled_2plate_kite.jl",
                   "coupled_2plate_kite_linear_vsm.jl",
                   "coupled_tether_deflection.jl",
                   "coupled_linearize.jl",
                   "sam_tutorial.jl"]
            @test isfile(joinpath("examples", f))
        end
    end

    @testset "Examples use GLMakie" begin
        for f in ["menu.jl", "hanging_mass.jl",
                  "catenary_line.jl", "pulley.jl",
                  "saddle_form.jl", "coupled_2plate_kite.jl"]
            content = read(joinpath("examples", f), String)
            @test occursin("using GLMakie", content)
        end
    end

    @testset "Run examples" begin
        Pkg.activate(joinpath(USER_DIR, "examples"))
        # The examples Project.toml has [sources] path=".."
        # which works in the repo but not in this temp dir,
        # so we need to explicitly dev SymbolicAWEModels here.
        Pkg.develop(path=REPO_ROOT)
        Pkg.instantiate()
        for f in ["hanging_mass.jl", "catenary_line.jl",
                  "pulley.jl", "saddle_form.jl",
                  "sam_tutorial.jl",
                  "coupled_tether_deflection.jl",
                  "coupled_2plate_kite.jl",
                  "coupled_2plate_kite_linear_vsm.jl",
                  "coupled_linearize.jl"]
            @testset "run $f" begin
                println("  Running $f...")
                # Each example runs in its own module to avoid
                # top-level variable conflicts between examples.
                mod = Module(Symbol(f))
                Base.include(mod,
                    joinpath(USER_DIR, "examples", f))
            end
        end
    end

    @testset "README pendulum example" begin
        println("  Running README pendulum example...")
        code = extract_readme_code(r"^### .*pendulum"i)
        eval(Meta.parseall(code))
    end

    @testset "README 2plate kite example" begin
        println("  Running README 2plate kite example...")
        code = extract_readme_code(r"^### .*2-Plate Kite"i)
        eval(Meta.parseall(code))
    end
end
nothing
