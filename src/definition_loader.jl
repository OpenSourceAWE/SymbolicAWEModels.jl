# Copyright (c) 2025 Bart van de Lint, Jelle Poland
# SPDX-License-Identifier: LGPL-3.0-only

using OrderedCollections: OrderedDict
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

"""
    load_sys_struct_from_yaml(yaml_path; system_name="from_yaml", set=nothing,
                              ignore_l0=false, dynamics_type=nothing,
                              aero_mode=nothing) -> (sys, transforms)

The structure of an authoring YAML and its KiteGeometry transforms, read, completed and
placed by KiteGeometry's `load_authoring`. `set` defaults to the `base` settings;
`dynamics_type` and `aero_mode` are written into every wing's row in place of its own.
"""
function load_sys_struct_from_yaml(yaml_path::AbstractString; system_name="from_yaml",
        set::Union{Nothing, Settings}=nothing, ignore_l0::Bool=false,
        dynamics_type::Union{Nothing, WingType}=nothing,
        aero_mode::Union{Nothing, AbstractAeroModel}=nothing)
    isnothing(set) && (set = load_settings("base"))
    authored = YAML.load_file(yaml_path; dicttype=OrderedDict{String, Any})
    isnothing(dynamics_type) ||
        set_wing_column!(authored, "dynamics_type", string(dynamics_type))
    isnothing(aero_mode) ||
        set_wing_column!(authored, "aero_mode", aero_mode_name(aero_mode))
    return load_authoring(authored; set, name=system_name, ignore_l0)
end

"""The name of an aero mode a wing row gives."""
aero_mode_name(mode) = string(nameof(typeof(mode)))
aero_mode_name(::AeroSpec{M}) where {M} = string(nameof(M))

"""Write `value` into column `column` of every row of the `wings` block of `authored`, an
authoring YAML, a table of `headers` and `data` or a list of mappings."""
function set_wing_column!(authored, column, value)
    wings = get(authored, "wings", nothing)
    isnothing(wings) && return authored
    rows = wings["data"]
    if haskey(wings, "headers")
        headers = wings["headers"]
        index = findfirst(==(column), headers)
        if isnothing(index)
            push!(headers, column)
            foreach(row -> push!(row, value), rows)
        end
        isnothing(index) || foreach(row -> row[index] = value, rows)
    else
        foreach(row -> row[column] = value, rows)
    end
    return authored
end

"""The tube row `tube` of `definition` as a run simulates it: its element, with the
rigidities it leaves unset from its material and its anchors measured from each body's
frame origin, the centres of the bodies' own masses lying `com_offsets` from them."""
function definition_tube(definition, tube, com_offsets)
    tube_model = tube.model
    columns = (; (field => getfield(tube_model, field)
                  for field in fieldnames(typeof(tube_model))
                  if !(field in (:anchor_a, :anchor_b)) &&
                     !isnothing(getfield(tube_model, field)))...)
    model = tube_model isa ElasticTubeModel ? ElasticTube(; columns...) :
        TimoshenkoTube(; columns...)
    model = with_law_rigidities(model, definition.materials[tube.material],
                                tube.diameter, tube.pressure)
    model.anchor_a_b .= something(tube_model.anchor_a, zeros(3)) + com_offsets[1]
    model.anchor_b_b .= something(tube_model.anchor_b, zeros(3)) + com_offsets[2]
    return KiteGeometry.Tube(; tube.name, tube.bodies, tube.diameter, tube.pressure,
                             tube.material, model)
end

# ==================== STRUCTURE DOCUMENTS ==================== #

"""
    sys_struct_from_document(document::AbstractDict)

The `SystemDefinition` of a parsed structure document, in the initial pose it gives.
"""
sys_struct_from_document(document::AbstractDict) = SystemDefinition(document)

"""
    save_structure_document(path, sys::SystemDefinition)

Write `structure_document``(sys)` to `path`, as YAML for a `.yml` or `.yaml`
extension and as JSON for `.json`.
"""
function save_structure_document(path::AbstractString, sys::SystemDefinition)
    for material in sys.materials
        material.law == "unit_stiffness" && error("material $(material.name) is a " *
            "nonlinear stiffness law, which a structure document cannot hold yet " *
            "(1-Bart-1/awesIO#1)")
    end
    document = structure_document(sys)
    write(path, document_format(path) == :json ? JSON.json(document, 2) :
                YAML.write(document))
    return path
end

"""
    load_structure_document(path)

The `SystemDefinition` of the structure document at `path`, YAML or JSON by extension.
"""
function load_structure_document(path::AbstractString)
    text = read(path, String)
    document = document_format(path) == :json ? JSON.parse(text) :
        YAML.load(text; dicttype=OrderedDict{String, Any})
    return sys_struct_from_document(document)
end

"""The encoding `path`'s extension asks for: `:yaml` or `:json`."""
function document_format(path::AbstractString)
    extension = lowercase(splitext(path)[2])
    extension in (".yml", ".yaml") && return :yaml
    extension == ".json" && return :json
    error("$path: a structure document is written as .yml, .yaml or .json.")
end

