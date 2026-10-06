# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

using JSON
using OrderedCollections: OrderedDict

"""Version of `structure_schema.yml` this writer emits and this reader accepts."""
const AWESIO_VERSION = "1.0.0"

"""The `metadata.schema` every conforming structure document carries."""
const STRUCTURE_SCHEMA = "structure_schema.yml"

"""
Columns of every block but `tubes`, each with its unit: first those
`structure_schema.yml` requires, in its order, then SymbolicAWEModels' own.
"""
const DOCUMENT_COLUMNS = OrderedDict(
    "points" => ["name" => "-", "type" => "-", "body" => "-", "pos_ENU" => "m",
                 "extra_mass" => "kg", "drag_area" => "m^2", "drag_coefficient" => "-",
                 "wing" => "-"],
    "segments" => ["name" => "-", "points" => "-", "l0" => "m", "diameter" => "m",
                   "density" => "kg/m^3", "unit_stiffness" => "N", "unit_damping" => "N*s",
                   "compression_frac" => "-", "compression_damping_frac" => "-"],
    "stations" => ["name" => "-", "wing" => "-", "type" => "-", "points" => "-",
                   "stiffness" => "N*m", "damping" => "N*m*s", "moment_frac" => "-"],
    "pulleys" => ["name" => "-", "segments" => "-", "type" => "-", "efficiency" => "-"],
    "tethers" => ["name" => "-", "start_point" => "-", "end_point" => "-",
                  "segments" => "-"],
    "winches" => ["name" => "-", "tethers" => "-", "winch_point" => "-",
                  "gear_ratio" => "-", "drum_radius" => "m", "model" => "-"],
    "wings" => ["name" => "-", "canopy_material" => "-", "aero" => "-"],
    "bodies" => ["name" => "-", "type" => "-", "pos_ENU" => "m", "Q_KA_to_ENU" => "-",
                 "extra_mass" => "kg", "extra_inertia_KA" => "kg*m^2", "wing" => "-",
                 "apparent_mass" => "kg"],
)

"""The columns every `tubes` row has, each with its unit."""
const TUBE_COLUMNS = ["name" => "-", "bodies" => "-", "diameter" => "m",
                      "pressure" => "Pa", "law" => "-", "model" => "-"]
"""Unit of each column a tube model ([`YAML_TUBE_MODELS`](@ref)) or a tube end adds."""
const TUBE_MODEL_UNITS = Dict(
    :EA => "N", :GA => "N", :GJ => "N*m^2", :EIy => "N*m^2", :EIz => "N*m^2",
    :shear_coeff => "-", :stiffness_axial => "N/m", :stiffness_shear => "N/m",
    :stiffness_torsion => "N*m", :stiffness_bending => "N*m", :damping => "s",
    :anchor_a => "m", :anchor_b => "m")

# ==================== SHARED HELPERS ==================== #

"""
    connectivity_sha(sections...) -> String

Lowercase hex SHA-256 of the connectivity preimage `structure_schema.yml`
documents. Each section is a count and its elements, each element a tuple of
one-based row numbers: the count, a semicolon, then every element comma separated
and semicolon terminated.
"""
function connectivity_sha(sections...)
    preimage = IOBuffer()
    for (count, elements) in sections
        print(preimage, count, ';')
        for element in elements
            join(preimage, element, ',')
            print(preimage, ';')
        end
    end
    return bytes2hex(sha256(take!(preimage)))
end

"""
    linear_rigidity(value, quantity) -> Float64

`value` as the number its stiffness column holds, erroring when it is a
nonlinear law that the document has no way to carry.
"""
function linear_rigidity(value, quantity::AbstractString)
    value isa Real && return Float64(value)
    error("$quantity is a nonlinear law ($(repr(value))); $STRUCTURE_SCHEMA " *
          "carries only a placeholder for one, so it does not round-trip " *
          "(1-Bart-1/awesIO#1).")
end

"""Three-component vector as the plain `Float64` list the document holds."""
vector3(value) = Float64[value[1], value[2], value[3]]

"""A 3×3 matrix as the three rows the document holds."""
matrix3(value) = [Float64[value[row, 1], value[row, 2], value[row, 3]] for row in 1:3]

"""A component name as a `Symbol`, passing an absent reference through."""
optional_symbol(::Nothing) = nothing
optional_symbol(name::AbstractString) = Symbol(name)

"""
    component_name(component) -> String

The name every other block refers `component` by, erroring when it has none.
"""
function component_name(component)
    isnothing(component.name) && error(
        "$(typeof(component)) $(component.idx) has no name; a structure " *
        "document refers to every component by name.")
    return string(component.name)
end

"""Name of `collection[idx]`, or `nothing` for the absent reference 0."""
ref_name(collection, idx::Integer) =
    idx == 0 ? nothing : component_name(collection[idx])

"""Names of `collection` at `idxs`, in order."""
ref_names(collection, idxs) = [component_name(collection[idx]) for idx in idxs]

"""
    document_connectivity(n_points, segment_points, n_bodies, tube_bodies) -> String

The `connectivity_sha` of a document with these points and segments, bodies and
tubes, and no canopy faces.
"""
document_connectivity(n_points, segment_points, n_bodies, tube_bodies) =
    connectivity_sha((n_points, segment_points), (n_bodies, tube_bodies), (0, ()))

# ==================== WRITING ==================== #

"""
    structure_document(sys::SystemStructure; name=sys.name, description="", note="")

Render `sys` in its initial pose as a document conforming to `structure_schema.yml`:
one `headers`/`units`/`data` table per block, every reference by name, every
position in the world frame. Returns nested `OrderedDict`s and `Vector`s, which
[`save_structure_document`](@ref) encodes as YAML or JSON.
"""
function structure_document(sys::SystemStructure; name::AbstractString=sys.name,
        description::AbstractString="", note::AbstractString="")
    tube_bodies = [(tube.body_a_idx, tube.body_b_idx) for tube in sys.tubes]
    sha = document_connectivity(length(sys.points),
        [segment.point_idxs for segment in sys.segments], length(sys.bodies), tube_bodies)
    document = OrderedDict{String, Any}(
        "metadata" => OrderedDict{String, Any}(
            "name" => name,
            "description" => description,
            "note" => note,
            "awesIO_version" => AWESIO_VERSION,
            "schema" => STRUCTURE_SCHEMA,
            "n_points" => length(sys.points),
            "connectivity_sha" => sha,
        ))
    document["points"] = document_table("points", point_rows(sys))
    document["segments"] = document_table("segments", segment_rows(sys))
    document["stations"] = document_table("stations", station_rows(sys))
    document["pulleys"] = document_table("pulleys", pulley_rows(sys))
    document["tethers"] = document_table("tethers", tether_rows(sys))
    document["winches"] = document_table("winches", winch_rows(sys))
    document["wings"] = document_table("wings", wing_rows(sys))
    document["bodies"] = document_table("bodies", body_rows(sys))
    document["tubes"] = tube_table(sys)
    return document
end

"""One `headers`/`units`/`data` block from `columns` (header => unit pairs)."""
table(columns, rows) = OrderedDict{String, Any}(
    "headers" => first.(columns), "units" => last.(columns), "data" => rows)

"""The `block` table of a structure document."""
document_table(block::AbstractString, rows) = table(DOCUMENT_COLUMNS[block], rows)

"""
Rows of the `points` block. The body a point is fixed to names the wing it is
part of, so `wing` carries only a point that has no body.
"""
function point_rows(sys::SystemStructure)
    return map(sys.points) do point
        body = ref_name(sys.bodies, point.body_idx)
        wing = isnothing(body) && point.is_wing_node ?
            component_name(sys.wings[point.wing_idx]) : nothing
        Any[component_name(point), string(point.type), body, vector3(point.pos_ENU),
            point.extra_mass, point.area, point.drag_coeff, wing]
    end
end

"""Rows of the `segments` block."""
segment_rows(sys::SystemStructure) =
    [Any[component_name(segment),
         ref_names(sys.points, segment.point_idxs),
         segment.l0, segment.diameter, segment.density,
         linear_rigidity(segment.unit_stiffness,
                         "segment $(segment.name) unit_stiffness"),
         segment.unit_damping,
         segment.compression_frac, segment.compression_damping_frac]
     for segment in sys.segments]

"""Rows of the `stations` block."""
station_rows(sys::SystemStructure) =
    [Any[component_name(station),
         station_wing(sys, station),
         string(station.type),
         ref_names(sys.points, station.point_idxs),
         station.stiffness, station.damping, station.moment_frac]
     for station in sys.stations]

"""
    station_wing(sys, station) -> String

Name of the wing whose twist `station` carries: its own reference where it has
one, and the wing its points sit on where it does not.
"""
function station_wing(sys::SystemStructure, station::Station)
    station.wing_idx == 0 || return component_name(sys.wings[station.wing_idx])
    wings = unique(point_wing(sys, sys.points[idx]) for idx in station.point_idxs)
    (length(wings) == 1 && !isnothing(only(wings))) || error(
        "Station $(station.name) names no wing and its points do not agree on " *
        "one ($(join(wings, ", "))); a structure document states it.")
    return only(wings)
end

"""Name of the wing a point is part of: its own reference, or its body's."""
function point_wing(sys::SystemStructure, point::Point)
    point.body_idx == 0 && return point.is_wing_node ?
        component_name(sys.wings[point.wing_idx]) : nothing
    body = sys.bodies[point.body_idx]
    return is_wing(body) ? component_name(body) : ref_name(sys.wings, body.wing_idx)
end

"""Rows of the `pulleys` block."""
pulley_rows(sys::SystemStructure) =
    [Any[component_name(pulley),
         ref_names(sys.segments, pulley.segment_idxs),
         string(pulley.type), pulley.efficiency]
     for pulley in sys.pulleys]

"""Rows of the `tethers` block."""
tether_rows(sys::SystemStructure) =
    [Any[component_name(tether),
         component_name(sys.points[tether.start_point_idx]),
         component_name(sys.points[tether.end_point_idx]),
         ref_names(sys.segments, tether.segment_idxs)]
     for tether in sys.tethers]

"""Rows of the `winches` block."""
winch_rows(sys::SystemStructure) =
    [Any[component_name(winch),
         ref_names(sys.tethers, winch.tether_idxs),
         component_name(sys.points[winch.winch_point_idx]),
         winch.gear_ratio, winch.drum_radius, string(nameof(typeof(winch.model)))]
     for winch in sys.winches]

"""Rows of the `wings` block: every wing is also the body of its name, with no canopy."""
wing_rows(sys::SystemStructure) =
    [Any[component_name(wing), nothing, string(nameof(typeof(wing.aero)))]
     for wing in sys.wings]

"""
Rows of the `bodies` block: each body's own mass properties, without the points it
carries, about the centre of that mass, which is the origin the row places.
"""
body_rows(sys::SystemStructure) =
    [Any[component_name(body), string(body.type), vector3(body_mass_centre(body)),
         Float64.(body.Q_KA_to_ENU), body.extra_mass, matrix3(body.extra_inertia_b),
         ref_name(sys.wings, body.wing_idx), body.apparent_mass]
     for body in sys.bodies]

"""World position [m] of the centre of `body`'s own mass in its initial pose."""
body_mass_centre(body) = body.pos_ENU + initial_rotation(body) * body.extra_com_offset_b

"""
    tube_table(sys) -> OrderedDict

The `tubes` block: the columns every row has, then the columns of each tube model
in `sys`, as the `tubes` table of the YAML loader reads them, then each end as an
offset from its body's mass centre in the body's frame.
"""
function tube_table(sys::SystemStructure)
    model_fields = unique(field for tube in sys.tubes
                          for field in last(YAML_TUBE_MODELS[tube_model_name(tube.model)]))
    columns = [TUBE_COLUMNS; [string(field) => TUBE_MODEL_UNITS[field]
                              for field in [model_fields; :anchor_a; :anchor_b]]]
    rows = map(sys.tubes) do tube
        body_a, body_b = sys.bodies[tube.body_a_idx], sys.bodies[tube.body_b_idx]
        Any[component_name(tube),
            ref_names(sys.bodies, (tube.body_a_idx, tube.body_b_idx)),
            tube.diameter, tube.pressure, string(tube.law), tube_model_name(tube.model),
            [tube_model_value(tube, field) for field in model_fields]...,
            vector3(tube.anchor_a_b - body_a.extra_com_offset_b),
            vector3(tube.anchor_b_b - body_b.extra_com_offset_b)]
    end
    return table(columns, rows)
end

"""The `model` column naming `model` in [`YAML_TUBE_MODELS`](@ref)."""
tube_model_name(model::AbstractTubeModel) =
    only(name for (name, (type, _)) in YAML_TUBE_MODELS if model isa type)

"""`tube.model`'s `field` as a document value, or `nothing` where it has none."""
function tube_model_value(tube::Tube, field::Symbol)
    hasfield(typeof(tube.model), field) || return nothing
    value = getfield(tube.model, field)
    field in rigidity_fields(tube.model) || return value
    return linear_rigidity(value, "tube $(tube.name) $field")
end

# ==================== READING ==================== #

"""
    sys_struct_from_document(doc::AbstractDict; set=nothing, vsm_set=nothing,
                             wind_mode=ProfileWind(), prn=true)

Build the `SystemStructure` a parsed structure document describes, starting in the
initial pose the document gives, without [`place!`](@ref). `set` supplies what the
schema has no column for — the winch friction and inertia, the tether material
defaults — and falls back to the `base` settings. A column the reader does not know is
ignored.

The document is the truth for what it carries: each body's own mass, inertia and
origin are taken from its row rather than re-derived from the points. A document
with wings is refused: it holds the wings where they are placed, but not the frame
their aerodynamic geometry is written in.
"""
function sys_struct_from_document(doc::AbstractDict; set=nothing, vsm_set=nothing,
        wind_mode::WindMode=ProfileWind(), prn::Bool=true)
    metadata = doc["metadata"]
    check_document_version(metadata)
    rows = Dict(block => document_rows(doc, block)
                for block in [collect(keys(DOCUMENT_COLUMNS)); "tubes"])
    check_connectivity(metadata, rows)
    isempty(rows["wings"]) || error(
        "The document has wings ($(join((row["name"] for row in rows["wings"]), ", "))); " *
        "reading a wing back is not supported yet (SymbolicAWEModels.jl#396).")
    resolved_set = isnothing(set) ? load_settings("base") : set

    sys_struct = SystemStructure(string(metadata["name"]), resolved_set;
        points = Point[read_point(row) for row in rows["points"]],
        segments = Segment[read_segment(row) for row in rows["segments"]],
        pulleys = Pulley[read_pulley(row) for row in rows["pulleys"]],
        tethers = Tether[read_tether(row) for row in rows["tethers"]],
        winches = Winch[read_winch(row, resolved_set) for row in rows["winches"]],
        bodies = Body[read_body(row) for row in rows["bodies"]],
        tubes = load_yaml_tubes(doc, Symbol),
        placed = true, vsm_set, wind_mode, prn)

    for row in rows["bodies"]
        apply_body_row!(sys_struct.bodies[Symbol(row["name"])], row)
    end
    update_mass_properties!(sys_struct; prn)
    return sys_struct
end

"""
    document_rows(doc, block) -> Vector{OrderedDict{String, Any}}

The rows of `block`, each a header → value mapping, so every column is addressed
by name. An absent block reads as an empty one.
"""
function document_rows(doc::AbstractDict, block::AbstractString)
    haskey(doc, block) || return OrderedDict{String, Any}[]
    table = doc[block]
    headers = string.(table["headers"])
    return [OrderedDict{String, Any}(zip(headers, row)) for row in table["data"]]
end

"""Refuse a document of an unknown schema or major version; warn on a minor one."""
function check_document_version(metadata::AbstractDict)
    metadata["schema"] == STRUCTURE_SCHEMA || error(
        "The document declares schema $(repr(metadata["schema"])); this " *
        "reader understands $STRUCTURE_SCHEMA.")
    version = string(metadata["awesIO_version"])
    version == AWESIO_VERSION && return nothing
    first(split(version, '.')) == first(split(AWESIO_VERSION, '.')) || error(
        "$STRUCTURE_SCHEMA $version is a major version this reader, written " *
        "against $AWESIO_VERSION, does not understand.")
    @warn "Reading a $STRUCTURE_SCHEMA $version document with a reader " *
          "written against $AWESIO_VERSION."
    return nothing
end

"""Check that `connectivity_sha` and `n_points` describe the document's own tables."""
function check_connectivity(metadata::AbstractDict, rows)
    points, bodies = rows["points"], rows["bodies"]
    length(points) == metadata["n_points"] || error(
        "The document holds $(length(points)) points but its metadata claims " *
        "$(metadata["n_points"]).")
    point_row = Dict(row["name"] => i for (i, row) in enumerate(points))
    body_row = Dict(row["name"] => i for (i, row) in enumerate(bodies))
    computed = document_connectivity(length(points),
        [Tuple(point_row[name] for name in row["points"]) for row in rows["segments"]],
        length(bodies),
        [Tuple(body_row[name] for name in row["bodies"]) for row in rows["tubes"]])
    computed == metadata["connectivity_sha"] || error(
        "connectivity_sha $(metadata["connectivity_sha"]) does not describe the " *
        "document's own points, segments, bodies and tubes (computed $computed).")
    return nothing
end

"""
    named_model(name, ::Type{T}) -> T

The `T` that `name` names, as an instance. The schema names a model without its
parameters, so the named type must take none.
"""
function named_model(name::AbstractString, ::Type{T}) where T
    symbol = Symbol(name)
    isdefined(@__MODULE__, symbol) || error(
        "The document names $name, which SymbolicAWEModels does not define.")
    model_type = getfield(@__MODULE__, symbol)
    model_type isa Type && model_type <: T || error(
        "The document names $name where a $T was expected.")
    try
        return model_type()
    catch err
        err isa MethodError || rethrow()
        error("$name takes parameters $STRUCTURE_SCHEMA does not carry, so it " *
              "cannot be built from its name alone.")
    end
end

"""
    optional_columns(row, columns) -> NamedTuple

The keyword arguments `row` carries among `columns` (column => keyword pairs), each
as a `Float64`, leaving out a column the row lacks so the constructor's default holds.
"""
optional_columns(row, columns) =
    (; (keyword => Float64(row[column]) for (column, keyword) in columns
        if !isnothing(get(row, column, nothing)))...)

read_point(row) = Point(Symbol(row["name"]), vector3(row["pos_ENU"]),
    parse_dynamics_type(row["type"]); body = optional_symbol(row["body"]),
    extra_mass = Float64(row["extra_mass"]), area = Float64(row["drag_area"]),
    drag_coeff = Float64(row["drag_coefficient"]))

read_segment(row) = Segment(Symbol(row["name"]), Symbol(row["points"][1]),
    Symbol(row["points"][2]),
    linear_rigidity(row["unit_stiffness"], "segment $(row["name"]) unit_stiffness"),
    Float64(get(row, "unit_damping", 0.0)), Float64(row["diameter"]);
    l0 = Float64(row["l0"]), density = Float64(row["density"]),
    optional_columns(row, ("compression_frac" => :compression_frac,
                           "compression_damping_frac" => :compression_damping_frac))...)

read_pulley(row) = Pulley(Symbol(row["name"]), Symbol(row["segments"][1]),
    Symbol(row["segments"][2]), parse_dynamics_type(row["type"]);
    efficiency = Float64(row["efficiency"]))

read_tether(row) = Tether(Symbol(row["name"]), Symbol.(row["segments"]);
    start_point = Symbol(row["start_point"]), end_point = Symbol(row["end_point"]))

function read_winch(row, set)
    model = get(row, "model", nothing)
    return Winch(Symbol(row["name"]), Symbol.(row["tethers"]),
        Float64(row["gear_ratio"]), Float64(row["drum_radius"]),
        set.f_coulomb, set.c_vf, set.inertia_total;
        winch_point = Symbol(row["winch_point"]),
        (isnothing(model) ? (;) :
         (; model = named_model(model, AbstractWinchModel)))...)
end

read_body(row) = Body(Symbol(row["name"]); extra_mass = Float64(row["extra_mass"]),
    inertia = yaml_matrix3(row["extra_inertia_KA"]), pos = vector3(row["pos_ENU"]),
    Q_b_to_w = Float64.(row["Q_KA_to_ENU"]), type = parse_dynamics_type(row["type"]))

"""
Overwrite the own mass properties `SystemStructure` derives for `body` with the
document's, about the origin its row places.
"""
function apply_body_row!(body::Body, row)
    body.extra_mass = Float64(row["extra_mass"])
    body.apparent_mass = Float64(get(row, "apparent_mass", 0.0))
    body.extra_inertia_b .= yaml_matrix3(row["extra_inertia_KA"])
    body.extra_com_offset_b .= 0.0
    return body
end

# ==================== YAML AND JSON ==================== #

"""
    save_structure_document(path, sys::SystemStructure; kwargs...)

Write [`structure_document`](@ref)`(sys; kwargs...)` to `path`, as YAML for a
`.yml` or `.yaml` extension and as JSON for `.json`.
"""
function save_structure_document(path::AbstractString, sys::SystemStructure;
                                 kwargs...)
    document = structure_document(sys; kwargs...)
    write(path, encode_document(document, document_format(path)))
    return path
end

"""
    save_log(logger::Logger, sys::SystemStructure, name="sim_log"; path="")

Save `logger` as an uncompressed .arrow file whose table metadata holds
[`structure_document`](@ref)`(sys)` as JSON under the key `topology`, read back into
`SysLog.metadata` by `load_log`.
"""
function KiteUtils.save_log(logger::Logger, sys::SystemStructure, name="sim_log";
                            path="")
    topology = encode_document(structure_document(sys), :json)
    return save_log(logger, name, false; path, metadata=Dict("topology" => topology))
end

"""
    load_structure_document(path; kwargs...)

Read the structure document at `path` — YAML or JSON by extension — and build the
`SystemStructure` it describes. `kwargs` reach
[`sys_struct_from_document`](@ref).
"""
function load_structure_document(path::AbstractString; kwargs...)
    text = read(path, String)
    document = document_format(path) == :json ? JSON.parse(text) : YAML.load(text)
    return sys_struct_from_document(document; kwargs...)
end

"""The encoding `path`'s extension asks for: `:yaml` or `:json`."""
function document_format(path::AbstractString)
    extension = lowercase(splitext(path)[2])
    extension in (".yml", ".yaml") && return :yaml
    extension == ".json" && return :json
    error("$path: a structure document is written as .yml, .yaml or .json.")
end

"""A structure document as text, in either encoding."""
encode_document(document, format::Symbol) =
    format == :json ? JSON.json(document, 2) : YAML.write(document)
