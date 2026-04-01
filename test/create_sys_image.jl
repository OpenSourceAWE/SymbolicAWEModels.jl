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

@info "Creating sysimage ..."
PackageCompiler.create_sysimage(
    [:Pkg, :TOML, :DocStringExtensions, :LinearAlgebra, :Parameters, :Printf, :Serialization, :SHA, :CodecZlib, :Tar, :Statistics, :Suppressor, :Timers, :GLMakie, :ModelingToolkit, :ControlSystemsBase, :RecipesBase, :StaticArrays, :SymbolicIndexingInterface, :NonlinearSolve, :OrdinaryDiffEqBDF, :OrdinaryDiffEqCore, :OrdinaryDiffEqNonlinearSolve, :SteadyStateDiffEq, :AtmosphericModels, :KiteUtils, :VortexStepMethod];
    sysimage_path="kps-image_tmp.so",
    precompile_execution_file=joinpath("test", "test_for_precompile.jl")
)
