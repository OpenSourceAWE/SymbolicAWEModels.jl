# Copyright (c) 2025 Bart van de Lint and Uwe Fechner
# SPDX-License-Identifier: LGPL-3.0-only

module SymbolicAWEModels

#======================================================================#
#                         DEPENDENCIES
#======================================================================#

# --- Julia Standard Library & General Utilities ---
using Pkg
using TOML
using DocStringExtensions
using LinearAlgebra
using Parameters
using Printf
using Serialization
using SHA
using CodecZlib
using Tar
using Statistics
using Suppressor
using Timers

# --- Numerical, Modeling & Scientific Computing ---
using ModelingToolkit
using ControlSystemsBase
using RecipesBase
using StaticArrays
using SymbolicIndexingInterface
using SymbolicIndexingInterface: AbstractIndexer

# --- Solvers (Nonlinear, Differential Equations) ---
using NonlinearSolve
using OrdinaryDiffEqBDF
using OrdinaryDiffEqCore
using OrdinaryDiffEqNonlinearSolve
using SteadyStateDiffEq

# --- Open Source AWE Packages ---
using AtmosphericModels
using KiteUtils
using VortexStepMethod
using ForwardDiff

#======================================================================#
#                  IMPORTS (for extending functions)
#======================================================================#

import KiteUtils: init!, next_step!, update_sys_state!, SysState
import ModelingToolkit: t_nounits as t, D_nounits as D
import ModelingToolkit.SciMLBase: successful_retcode, init

#======================================================================#
#                          EXPORTS
#                 (The Public API of this Module)
#======================================================================#

# --- KiteUtils ---
export update_from_sysstate!, get_data_path, set_data_path, se
export SysState, SysLog, Settings, AbstractKiteModel
export Logger, log!, save_log, load_log
export load_settings

# --- Types ---
# Core Model
export SymbolicAWEModel
# System Structure Components
export SystemStructure, Point, Group, Segment, Pulley, Tether, Winch, Wing, Transform
export AbstractWing, BaseWing, VSMWing
export NameRef, NamedCollection, WeightedRefPoints
# Enums
export DynamicsType, DYNAMIC, QUASI_STATIC, WING, STATIC
export SegmentType, POWER_LINE, STEERING_LINE, BRIDLE
export WingType, QUATERNION, REFINE
export AeroMode, AERO_NONE, AERO_DIRECT, AERO_LINEARIZED

# --- High-Level Simulation Functions (Workers) ---
export sim!, sim_reposition!

# --- Low-Level Simulation Functions ---
export find_steady_state!
export linearize!
export update_segment_forces!
export set_world_frame_damping
export set_body_frame_damping
export segment_stretch_stats
export calc_steady_torque

# --- Getter Functions ---
export winch_force
export unstretched_length
export tether_length

# --- Helper Functions ---
export init_module
export update_plot_observables!
export animate
export load_sys_struct_from_yaml
export update_yaml_from_sys_struct!
export update_sys_struct_from_yaml!
export replay
export record
export plot_sphere_trajectory
export plot_body_frame
export plot_aoa

set_zero_subnormals(true)       # required to avoid drastic slow down on Intel CPUs when numbers become very small

#======================================================================#
#                       TYPE DEFINITIONS
#======================================================================#

# Type definitions
"""
    const SimFloat = Float64

This type is used for all real variables, used in the Simulation. Possible alternatives: Float32, Double64, Dual
Other types than Float64 or Float32 do require support of Julia types by the solver. 
"""
const SimFloat = Float64

"""
   const KVec3    = MVector{3, SimFloat}

Basic 3-dimensional vector, stack allocated, mutable.
"""
const KVec3    = MVector{3, SimFloat}
const KVec4    = MVector{4, SimFloat}
const MVec3    = MVector{3, Float64}  # Used by VortexStepMethod functions

"""
   const SVec3    = SVector{3, SimFloat}

Basic 3-dimensional vector, stack allocated, immutable.
"""
const SVec3    = SVector{3, SimFloat}  

# Defined in ext/SymbolicAWEModelsMakieExt.jl
function plot end
# Defined in ext/SymbolicAWEModelsMakieExt.jl
function plot! end
function update_plot_observables! end
function animate end
function replay end
function record end
function plot_sphere_trajectory end
function plot_body_frame end
function plot_aoa end
function find_steady_state! end
function update_yaml_from_sys_struct! end
function make_lin_sys_state end
function create_model_archive end

function __init__()
    data_dir = joinpath(pwd(), "data")
    if isdir(data_dir) && isfile(joinpath(data_dir, "2plate_kite", "system.yaml"))
        set_data_path(data_dir)
    end
end

include("system_structure/system_structure.jl")
include("vsm_refine.jl")
include("symbolic_awe_model.jl")
include("model_management.jl")
include("yaml_loader.jl")
include("tether_properties.jl")
include("linearize.jl")
include("generate_system/generate_system.jl")
include("simulate.jl")

# rotate a 3d vector around the x axis in the yz plane - following the right hand rule
function rotate_around_x(vec, angle::T) where T
    result = zeros(T, 3)
    result[1] = vec[1]
    result[2] = cos(angle) * vec[2] - sin(angle) * vec[3]
    result[3] = sin(angle) * vec[2] + cos(angle) * vec[3]
    result
end

# rotate a 3d vector around the y axis in the xz plane - following the right hand rule
function rotate_around_y(vec, angle::T) where T
    result = zeros(T, 3)
    result[1] = cos(angle) * vec[1] + sin(angle) * vec[3]
    result[2] = vec[2]
    result[3] = -sin(angle) * vec[1] + cos(angle) * vec[3]
    result
end

# rotate a 3d vector around the z axis in the yx plane - following the right hand rule
function rotate_around_z(vec, angle::T) where T
    result = zeros(T, 3)
    result[1] = cos(angle) * vec[1] - sin(angle) * vec[2]
    result[2] = sin(angle) * vec[1] + cos(angle) * vec[2]
    result[3] = vec[3]
    result
end

"""
    copy_examples()

Copy all example scripts to the folder "examples"
(it will be created if it doesn't exist).
"""
function copy_examples(; force=false)
    src_data_path = joinpath(@__DIR__, "..", "examples")
    dst_data_path = abspath(joinpath(pwd(), "examples"))
    copy_dir(src_data_path, dst_data_path; force)
end

"""
    copy_bin()

Copy all example scripts to the folder "bin"
(it will be created if it doesn't exist).
"""
function copy_bin(; force=false)
    src_data_path = joinpath(@__DIR__, "..", "bin")
    dst_data_path = abspath(joinpath(pwd(), "bin"))
    copy_dir(src_data_path, dst_data_path; force)
end

"""
    copy_data()

Copy all data scripts to the folder "data"
(it will be created if it doesn't exist).
"""
function copy_data(; force=false)
    src_data_path = joinpath(@__DIR__, "..", "data")
    dst_data_path = abspath(joinpath(pwd(), "data"))
    copy_dir(src_data_path, dst_data_path; force)
end

"""
    copy_dir(src_dir, dst_dir)

Copies all files from `src_dir` to `dst_dir`.
Overwrites existing files if force=true.
Creates `dst_dir` if it does not exist.
"""
function copy_dir(src_dir::AbstractString, dst_dir::AbstractString; force=false)
    if !isdir(dst_dir)
        mkdir(dst_dir)
    end
    for file in readdir(src_dir)
        src_file = joinpath(src_dir, file)
        dst_file = joinpath(dst_dir, file)
        if force || (isfile(src_file) && !isfile(dst_file))
            cp(src_file, dst_file; force=true)
            chmod(dst_file, 0o774)
        elseif isdir(src_file)
            copy_dir(src_file, dst_file; force)
        end
    end
end

"""
    get_example_packages()

Get the list of packages from examples/Project.toml, excluding SymbolicAWEModels itself.
This ensures init_module installs the correct dependencies for running examples.
"""
function get_example_packages()
    examples_project_path = joinpath(@__DIR__, "..", "examples", "Project.toml")
    if !isfile(examples_project_path)
        @warn "examples/Project.toml not found, using default package list"
        return ["KiteUtils", "GLMakie"]
    end

    examples_project = TOML.parsefile(examples_project_path)
    deps = get(examples_project, "deps", Dict())
    # Exclude SymbolicAWEModels itself (it's already in the user's project)
    return sort([name for name in keys(deps) if name != "SymbolicAWEModels"])
end

"""
    init_module(; force=false, add_pkg=true)

Initialize the module in the current working directory.

This function performs the following actions:

- Copies all files from the module's `data` directory to the current working directory's `data` folder (`pwd()/data`). Existing files in the destination are NOT overwritten unless `force=true`.
- Copies all example scripts from the module to the current working directory's `examples` folder (`pwd()/examples`). The folder is created if it does not exist. Existing files are NOT overwritten unless `force=true`.
- Installs all required packages if they are not already installed. This occurs only if `add_pkg=true` (default). The packages are automatically determined from `examples/Project.toml`.

# Keyword Arguments
- `force::Bool=false`: If `true`, existing files in the destination directories will be overwritten. If `false` (default), existing files will be preserved.
- `add_pkg::Bool=true`: If `true` (default), installs required packages if they are not already present. If `false`, package installation is skipped.
"""
function init_module(; force=false, add_pkg=true)
    copy_data(; force)
    copy_examples(; force)

    if add_pkg
        # Install required packages if not already present
        pkgs = get_example_packages()
        println("Installing example dependencies: ", join(pkgs, ", "))
        for pkg in pkgs
            if !(pkg in keys(Pkg.project().dependencies))
                Pkg.add(pkg)
            end
        end
    end

    println("Initialization complete! Examples and data files are prepared in the current directory.")
end

include("precompile.jl")

end
