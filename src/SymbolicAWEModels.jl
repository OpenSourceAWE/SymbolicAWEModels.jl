# Copyright (c) 2020, 2021, 2022, 2024 Uwe Fechner, Bart van de Lint and Daan van Wolffelaar
# SPDX-License-Identifier: MIT

#= Models of a kite-power system in implicit form: residual = f(y, yd)

This model implements a 3D mass-spring system with reel-out. It uses six tether segments (the number can be
configured in the file data/settings.yaml). Two kite models are provided, the one point and the four point
kite model. The spring constant and the damping decrease with the segment length. The aerodynamic kite forces are
calculated, depending on reel-out speed, depower and steering settings. 

Scientific background: http://arxiv.org/abs/1406.6218 =#

module SymbolicAWEModels

using PrecompileTools: @setup_workload, @compile_workload 
using Serialization, StaticArrays, LinearAlgebra, Statistics, Parameters,
      DocStringExtensions, OrdinaryDiffEqCore, OrdinaryDiffEqBDF, OrdinaryDiffEqNonlinearSolve,
      NonlinearSolve, SHA
using KiteUtils, WinchModels, AtmosphericModels
using Pkg
using VortexStepMethod
import Base.zero
import KiteUtils.AbstractKiteModel
import KiteUtils.calc_elevation
import KiteUtils.calc_heading
import KiteUtils.calc_course
import KiteUtils.SysState
import OrdinaryDiffEqCore.init
import OrdinaryDiffEqCore.step!
using ModelingToolkit, SymbolicIndexingInterface
using ModelingToolkit: t_nounits as t, D_nounits as D
using ADTypes: AutoFiniteDiff
import ModelingToolkit.SciMLBase: successful_retcode
import OrdinaryDiffEqCore.ODEIntegrator

export SymbolicAWEModel                                 # constants and types
export copy_examples, copy_bin, update_sys_state!       # helper functions
export find_steady_state!                               # low level workers
export init_sim!, init!, reinit!, next_step!            # high level workers
export winch_force, unstretched_length, tether_length # getters
export create_ram_sys_struct, create_simple_ram_sys_struct
import LinearAlgebra: norm
export SystemStructure, Point, Group, Segment, Pulley, Tether, Winch, Wing, Transform
export DynamicsType, DYNAMIC, QUASI_STATIC, WING, STATIC
export SegmentType, POWER_LINE, STEERING_LINE, BRIDLE

set_zero_subnormals(true)       # required to avoid drastic slow down on Intel CPUs when numbers become very small

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

"""
   const SVec3    = SVector{3, SimFloat}

Basic 3-dimensional vector, stack allocated, immutable.
"""
const SVec3    = SVector{3, SimFloat}  

# Defined in ext/SymbolicAWEModelsControlPlotsExt.jl
function plot end

function __init__()
    if isdir(joinpath(pwd(), "data")) && isfile(joinpath(pwd(), "data", "system.yaml"))
        set_data_path(joinpath(pwd(), "data"))
    end
end

include("system_structure.jl")
include("symbolic_awe_model.jl") # include code, specific for the ram air kite model
include("linearize.jl")
include("mtk_model.jl")

function upwind_dir(v_wind_gnd)
    if v_wind_gnd[1] == 0.0 && v_wind_gnd[2] == 0.0
        return NaN
    end
    wind_dir = atan(v_wind_gnd[2], v_wind_gnd[1])
    -(wind_dir + π/2)
end

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
function copy_examples()
    PATH = "examples"
    if ! isdir(PATH) 
        mkdir(PATH)
    end
    src_path = joinpath(dirname(pathof(@__MODULE__)), "..", PATH)
    copy_files("examples", readdir(src_path))
end

function copy_model_settings()
    files = ["settings.yaml", "ram_air_kite_body.obj", "ram_air_kite_foil.dat", "system.yaml", "settings_ram.yaml", 
             "system_ram.yaml", "ram_air_kite_foil_cd_polar.csv", "ram_air_kite_foil_cl_polar.csv", "ram_air_kite_foil_cm_polar.csv"]
    dst_path = abspath(joinpath(pwd(), "data"))
    copy_files("data", files)
    set_data_path(joinpath(pwd(), "data"))
    println("Copied $(length(files)) files to $(dst_path) !")
end

function install_examples(add_packages=true)
    copy_examples()
    copy_settings()
    copy_bin()
    copy_model_settings()
    if add_packages
        Pkg.add(["KiteUtils", "KitePodModels", "WinchModels", "ControlPlots", 
                 "LaTeXStrings", "StatsBase", "Timers", "Rotations"])
    end
end

function copy_files(relpath, files)
    if ! isdir(relpath) 
        mkdir(relpath)
    end
    src_path = joinpath(dirname(pathof(@__MODULE__)), "..", relpath)
    for file in files
        cp(joinpath(src_path, file), joinpath(relpath, file), force=true)
        chmod(joinpath(relpath, file), 0o774)
    end
    files
end

"""
    copy_bin()

Copy the scripts create_sys_image and run_julia to the folder "bin"
(it will be created if it doesn't exist).
"""
function copy_bin()
    PATH = "bin"
    if ! isdir(PATH) 
        mkdir(PATH)
    end
    src_path = joinpath(dirname(pathof(@__MODULE__)), "..", PATH)
    cp(joinpath(src_path, "create_sys_image2"), joinpath(PATH, "create_sys_image"), force=true)
    cp(joinpath(src_path, "run_julia"), joinpath(PATH, "run_julia"), force=true)
    chmod(joinpath(PATH, "create_sys_image"), 0o774)
    chmod(joinpath(PATH, "run_julia"), 0o774)
    PATH = "test"
    if ! isdir(PATH) 
        mkdir(PATH)
    end
    src_path = joinpath(dirname(pathof(@__MODULE__)), "..", PATH)
    cp(joinpath(src_path, "create_sys_image2.jl"), joinpath(PATH, "create_sys_image.jl"), force=true)
    cp(joinpath(src_path, "test_for_precompile.jl"), joinpath(PATH, "test_for_precompile.jl"), force=true)
    cp(joinpath(src_path, "update_packages.jl"), joinpath(PATH, "update_packages.jl"), force=true)
    chmod(joinpath(PATH, "create_sys_image.jl"), 0o664)
    chmod(joinpath(PATH, "test_for_precompile.jl"), 0o664)
    chmod(joinpath(PATH, "update_packages.jl"), 0o664)
end

include("precompile.jl")

end
