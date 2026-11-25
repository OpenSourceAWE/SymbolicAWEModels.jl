# Copyright (c) 2025 Bart van de Lint, Jelle Poland
# SPDX-License-Identifier: LGPL-3.0-only

using YAML: YAML

"""
    project_file(entry, project=KiteUtils.PROJECT)

Path of the file the project file names under `system.<entry>`, resolved against
the data path. `entry` is `"sim_settings"`, `"vsm_settings"`, `"aero_geometry"`
or `"structural_geometry"`. Throws an `ArgumentError` naming the project file
where it does not exist or names no such file.
"""
function project_file(entry, project=KiteUtils.PROJECT)
    path = optional_project_file(entry, project)
    isnothing(path) && throw(ArgumentError(
        "$(joinpath(get_data_path(), project)) names no system.$entry"))
    return path
end

"""
    optional_project_file(entry, project=KiteUtils.PROJECT)

As [`project_file`](@ref), but `nothing` where the project file does not exist
or names no such file.
"""
function optional_project_file(entry, project=KiteUtils.PROJECT)
    data_path = get_data_path()
    project_path = joinpath(data_path, project)
    isfile(project_path) || return nothing
    entries = get(YAML.load_file(project_path), "system", nothing)
    entries isa AbstractDict && haskey(entries, entry) || return nothing
    return joinpath(data_path, dirname(project), entries[entry])
end

function parse_wing_type(text::String)
    text_upper = uppercase(text)
    text_upper == "PARTICLE_DYNAMICS" && return PARTICLE_DYNAMICS
    text_upper == "RIGID_DYNAMICS" && return RIGID_DYNAMICS
    error("Unknown WingType: $text")
end

function parse_principal_frame_method(text::String)
    text_upper = uppercase(text)
    text_upper == "EIGEN_DECOMP" && return EIGEN_DECOMP
    text_upper == "Y_ROTATION" && return Y_ROTATION
    error("Unknown PrincipalFrameMethod: $text")
end

function parse_aero_mode(text::String)
    key = lowercase(replace(text, "_" => ""))
    key in ("aeronone", "none") && return AeroNone()
    key in ("aerodirect", "direct") && return AeroDirect()
    key in ("aerolinearized", "linearized") && return AeroLinearized()
    key in ("aeroplate", "plate") && return AeroPlate()
    key in ("continuousaero", "continuous") && return ContinuousAero()
    key in ("aeropressure", "pressure") && return AeroPressure()
    error("Unknown aero model: $text")
end
