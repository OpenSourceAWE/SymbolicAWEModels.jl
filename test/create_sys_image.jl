# SPDX-FileCopyrightText: 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

using PackageCompiler

# --- Julia Standard Library & General Utilities ---
using Pkg, TOML, DocStringExtensions, LinearAlgebra, Parameters, Printf, Serialization, SHA,
      CodecZlib, Tar, Statistics, Suppressor, Timers, GLMakie

# --- Numerical, Modeling & Scientific Computing ---
using ModelingToolkit, ControlSystemsBase, RecipesBase, StaticArrays, SymbolicIndexingInterface

# --- Solvers (Nonlinear, Differential Equations) ---
using NonlinearSolve, OrdinaryDiffEqBDF, OrdinaryDiffEqCore, OrdinaryDiffEqNonlinearSolve, SteadyStateDiffEq

# --- Open Source AWE Packages ---
using AtmosphericModels, KiteUtils, VortexStepMethod

GC.gc(true)
let mem = Sys.free_memory() / 1024^2
    @info "Free memory: $(round(mem; digits=1)) MB"
    swap_gb = 0.0
    if Sys.islinux()
        swap_info = read(`swapon --show --bytes --noheadings`, String)
        if !isempty(strip(swap_info))
            swap_size = sum(parse(Int, split(line)[3]) for line in split(strip(swap_info), '\n') if !isempty(line))
            swap_gb = swap_size / 1024^3
            @info "Swap size: $(round(swap_gb; digits=1)) GB"
        else
            @info "No swap configured"
        end
    elseif Sys.iswindows()
        pf_info = read(`wmic pagefile list /format:list`, String)
        m = match(r"AllocatedBaseSize=(\d+)", pf_info)
        if m !== nothing
            swap_gb = parse(Int, m.captures[1]) / 1024
            @info "Page file size: $(round(swap_gb; digits=1)) GB"
        else
            @info "Page file size: unknown"
        end
    end
    if haskey(ENV, "JULIA_IMAGE_THREADS")
        @info "JULIA_IMAGE_THREADS: $(ENV["JULIA_IMAGE_THREADS"])"
    else
        free_gb = mem / 1024
        if free_gb + swap_gb < 36.0
            @error "JULIA_IMAGE_THREADS is not defined and total available memory ($(round(free_gb + swap_gb; digits=1)) GB free RAM + swap) is less than 36 GB. System image creation may fail!"
        else
            @info "JULIA_IMAGE_THREADS not defined!"
        end
    end
end

@info "Creating sysimage ..."
PackageCompiler.create_sysimage(
    [:Pkg, :TOML, :DocStringExtensions, :LinearAlgebra, :Parameters, :Printf, :Serialization, :SHA, :CodecZlib, :Tar, :Statistics, :Suppressor, :Timers, :GLMakie, :ModelingToolkit, :ControlSystemsBase, :RecipesBase, :StaticArrays, :SymbolicIndexingInterface, :NonlinearSolve, :OrdinaryDiffEqBDF, :OrdinaryDiffEqCore, :OrdinaryDiffEqNonlinearSolve, :SteadyStateDiffEq, :AtmosphericModels, :KiteUtils, :VortexStepMethod];
    sysimage_path="kps-image_tmp.so",
    precompile_execution_file=joinpath("test", "test_for_precompile.jl")
)
