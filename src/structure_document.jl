# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

using JSON
using OrderedCollections: OrderedDict

"""Version of `structure_schema.yml` this writer emits and this reader accepts."""
const AWESIO_VERSION = "0.1.0"

"""The `metadata.schema` every conforming structure document carries."""
const STRUCTURE_SCHEMA = "structure_schema.yml"

"""Columns of every document block, in the order `structure_schema.yml` fixes them."""
const DOCUMENT_HEADERS = OrderedDict(
    "points" => ["name", "type", "body", "wing", "pos_cad"],
    "segments" => ["name", "point_a", "point_b", "l0", "diameter", "density",
                   "unit_stiffness", "unit_damping", "compression_frac",
                   "compression_damping_frac"],
    "stations" => ["name", "type", "points", "moment_frac", "stiffness", "damping"],
    "pulleys" => ["name", "segment_a", "segment_b", "type", "efficiency"],
    "tethers" => ["name", "start_point", "end_point", "segments"],
    "winches" => ["name", "tethers", "winch_point", "model", "gear_ratio",
                  "drum_radius"],
    "bodies" => ["name", "type", "aero", "wing", "mass", "apparent_mass",
                 "inertia_principal", "com_offset_b", "pos_cad"],
    "elastic_joints" => ["name", "body_a", "body_b", "anchor_a", "anchor_b",
                         "stiffness_axial", "stiffness_shear", "stiffness_torsion",
                         "stiffness_bending", "damping", "radius"],
    "timoshenko_joints" => ["name", "body_a", "body_b", "anchor_a", "anchor_b",
                            "EA", "GA", "GJ", "EIy", "EIz", "shear_coeff",
                            "damping", "rest_length", "radius"],
)

# ==================== SHARED HELPERS ==================== #

"""
    connectivity_sha(n_points, endpoints) -> String

Lowercase hex SHA-256 of the connectivity preimage `structure_schema.yml`
documents: the point count, a semicolon, then every segment's two endpoints as
one-based row numbers into the points block, comma separated and semicolon
terminated.
"""
function connectivity_sha(n_points::Integer, endpoints)
    preimage = IOBuffer()
    print(preimage, n_points, ';')
    for (point_a, point_b) in endpoints
        print(preimage, point_a, ',', point_b, ';')
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

"""A component name as a `Symbol`, passing an absent reference through."""
optional_symbol(::Nothing) = nothing
optional_symbol(name::AbstractString) = Symbol(name)

"""A length as a `Float64`, passing an absent quantity through."""
optional_length(::Nothing) = nothing
optional_length(value) = Float64(value)

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

# ==================== WRITING ==================== #

"""
    structure_document(sys::SystemStructure; name=sys.name, description="", note="")

Render `sys` as a document conforming to `structure_schema.yml`: one
`headers`/`data` table per block, every reference by name, every component
carrying its own geometry and material. Returns nested `OrderedDict`s and
`Vector`s, which [`save_structure_document`](@ref) encodes as YAML or JSON.

The schema describes structure only, so the transforms that place the system in
the world, the live state, and the point masses and drag properties it has no
column for are left out.
"""
function structure_document(sys::SystemStructure; name::AbstractString=sys.name,
        description::AbstractString="", note::AbstractString="")
    endpoints = [segment.point_idxs for segment in sys.segments]
    document = OrderedDict{String, Any}(
        "metadata" => OrderedDict{String, Any}(
            "name" => name,
            "description" => description,
            "note" => note,
            "awesIO_version" => AWESIO_VERSION,
            "schema" => STRUCTURE_SCHEMA,
            "n_points" => length(sys.points),
            "connectivity_sha" => connectivity_sha(length(sys.points), endpoints),
        ))
    document["points"] = document_table("points", point_rows(sys))
    document["segments"] = document_table("segments", segment_rows(sys))
    document["stations"] = document_table("stations", station_rows(sys))
    document["pulleys"] = document_table("pulleys", pulley_rows(sys))
    document["tethers"] = document_table("tethers", tether_rows(sys))
    document["winches"] = document_table("winches", winch_rows(sys))
    document["bodies"] = document_table("bodies", body_rows(sys))
    document["elastic_joints"] =
        document_table("elastic_joints", elastic_joint_rows(sys))
    document["timoshenko_joints"] =
        document_table("timoshenko_joints", timoshenko_joint_rows(sys))
    return document
end

"""One `headers`/`data` block of a structure document."""
document_table(block::AbstractString, rows) = OrderedDict{String, Any}(
    "headers" => DOCUMENT_HEADERS[block], "data" => rows)

"""Rows of the `points` block. A point names its wing when it is a wing node."""
point_rows(sys::SystemStructure) =
    [Any[component_name(point),
         string(point.type),
         ref_name(sys.bodies, point.body_idx),
         point.is_wing_node ? component_name(sys.wings[point.wing_idx]) : nothing,
         vector3(point.pos_cad)]
     for point in sys.points]

"""Rows of the `segments` block."""
segment_rows(sys::SystemStructure) =
    [Any[component_name(segment),
         component_name(sys.points[segment.point_idxs[1]]),
         component_name(sys.points[segment.point_idxs[2]]),
         segment.l0, segment.diameter, segment.density,
         linear_rigidity(segment.unit_stiffness,
                         "segment $(segment.name) unit_stiffness"),
         segment.unit_damping,
         segment.compression_frac, segment.compression_damping_frac]
     for segment in sys.segments]

"""Rows of the `stations` block."""
station_rows(sys::SystemStructure) =
    [Any[component_name(station),
         string(station.type),
         ref_names(sys.points, station.point_idxs),
         station.moment_frac, station.stiffness, station.damping]
     for station in sys.stations]

"""Rows of the `pulleys` block."""
pulley_rows(sys::SystemStructure) =
    [Any[component_name(pulley),
         component_name(sys.segments[pulley.segment_idxs[1]]),
         component_name(sys.segments[pulley.segment_idxs[2]]),
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
         string(nameof(typeof(winch.model))),
         winch.gear_ratio, winch.drum_radius]
     for winch in sys.winches]

"""Rows of the `bodies` block. A wing is a body whose `aero` is not null."""
body_rows(sys::SystemStructure) =
    [Any[component_name(body),
         string(body.type),
         body.aero isa AeroNone ? nothing : string(nameof(typeof(body.aero))),
         ref_name(sys.wings, body.wing_idx),
         body.mass, body.apparent_mass,
         vector3(body.inertia_principal), vector3(body.com_offset_b),
         vector3(body.pos_cad)]
     for body in sys.bodies]

"""Rows of the `elastic_joints` block."""
elastic_joint_rows(sys::SystemStructure) =
    [Any[component_name(joint),
         component_name(sys.bodies[joint.body_a_idx]),
         component_name(sys.bodies[joint.body_b_idx]),
         vector3(joint.anchor_a_b), vector3(joint.anchor_b_b),
         linear_rigidity(joint.stiffness_axial, "joint $(joint.name) axial"),
         linear_rigidity(joint.stiffness_shear, "joint $(joint.name) shear"),
         linear_rigidity(joint.stiffness_torsion, "joint $(joint.name) torsion"),
         linear_rigidity(joint.stiffness_bending, "joint $(joint.name) bending"),
         joint.damping, joint.radius]
     for joint in sys.elastic_joints]

"""Rows of the `timoshenko_joints` block."""
timoshenko_joint_rows(sys::SystemStructure) =
    [Any[component_name(joint),
         component_name(sys.bodies[joint.body_a_idx]),
         component_name(sys.bodies[joint.body_b_idx]),
         vector3(joint.anchor_a_b), vector3(joint.anchor_b_b),
         linear_rigidity(joint.EA, "joint $(joint.name) EA"),
         linear_rigidity(joint.GA, "joint $(joint.name) GA"),
         linear_rigidity(joint.GJ, "joint $(joint.name) GJ"),
         linear_rigidity(joint.EIy, "joint $(joint.name) EIy"),
         linear_rigidity(joint.EIz, "joint $(joint.name) EIz"),
         joint.shear_coeff, joint.damping, joint.rest_length, joint.radius]
     for joint in sys.timoshenko_joints]

# ==================== READING ==================== #

"""
    sys_struct_from_document(doc::AbstractDict; set=nothing, vsm_set=nothing,
                             wind_mode=ProfileWind(), prn=true)

Build the `SystemStructure` a parsed structure document describes. `set`
supplies what the schema has no column for — the winch friction and inertia, the
tether material defaults — and falls back to the `base` settings.

The document is the truth for what it carries: each body's mass, inertia, COM
offset and CAD origin are taken from its row rather than re-derived from the
points.
"""
function sys_struct_from_document(doc::AbstractDict; set=nothing, vsm_set=nothing,
        wind_mode::WindMode=ProfileWind(), prn::Bool=true)
    metadata = doc["metadata"]
    check_document_version(metadata)
    rows = Dict(block => document_rows(doc, block) for block in keys(DOCUMENT_HEADERS))
    check_connectivity(metadata, rows["points"], rows["segments"])
    resolved_set = isnothing(set) ? load_settings("base") : set

    wing_rows = filter(is_wing_row, rows["bodies"])
    stations = Station[read_station(row) for row in rows["stations"]]

    sys_struct = SystemStructure(string(metadata["name"]), resolved_set;
        points = Point[read_point(row) for row in rows["points"]],
        stations,
        segments = Segment[read_segment(row) for row in rows["segments"]],
        pulleys = Pulley[read_pulley(row) for row in rows["pulleys"]],
        tethers = Tether[read_tether(row) for row in rows["tethers"]],
        winches = Winch[read_winch(row, resolved_set) for row in rows["winches"]],
        wings = Body[read_wing(row, rows["stations"], rows["points"],
                               resolved_set, vsm_set) for row in wing_rows],
        bodies = Body[read_body(row) for row in rows["bodies"] if !is_wing_row(row)],
        elastic_joints = ElasticJoint[read_elastic_joint(row)
                                      for row in rows["elastic_joints"]],
        timoshenko_joints = TimoshenkoJoint[read_timoshenko_joint(row)
                                            for row in rows["timoshenko_joints"]],
        vsm_set, wind_mode, prn)

    for row in rows["bodies"]
        apply_body_row!(sys_struct.bodies[Symbol(row["name"])], row)
    end
    reinit!(sys_struct, resolved_set; prn)
    return sys_struct
end

"""A body row carries aero, which is what makes it a wing."""
is_wing_row(row) = !isnothing(row["aero"])

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
function check_connectivity(metadata::AbstractDict, points, segments)
    length(points) == metadata["n_points"] || error(
        "The document holds $(length(points)) points but its metadata claims " *
        "$(metadata["n_points"]).")
    row_number = Dict(row["name"] => i for (i, row) in enumerate(points))
    endpoints = [(row_number[row["point_a"]], row_number[row["point_b"]])
                 for row in segments]
    computed = connectivity_sha(length(points), endpoints)
    computed == metadata["connectivity_sha"] || error(
        "connectivity_sha $(metadata["connectivity_sha"]) does not describe the " *
        "document's own points and segments (computed $computed).")
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

read_point(row) = Point(Symbol(row["name"]), vector3(row["pos_cad"]),
    parse_dynamics_type(row["type"]); body = optional_symbol(row["body"]),
    wing = optional_symbol(row["wing"]))

read_segment(row) = Segment(Symbol(row["name"]), Symbol(row["point_a"]),
    Symbol(row["point_b"]),
    linear_rigidity(row["unit_stiffness"], "segment $(row["name"]) unit_stiffness"),
    Float64(row["unit_damping"]), Float64(row["diameter"]);
    l0 = Float64(row["l0"]), density = Float64(row["density"]),
    compression_frac = Float64(row["compression_frac"]),
    compression_damping_frac = Float64(row["compression_damping_frac"]))

read_station(row) = Station(Symbol(row["name"]), Symbol.(row["points"]),
    parse_dynamics_type(row["type"]), Float64(row["moment_frac"]);
    stiffness = Float64(row["stiffness"]), damping = Float64(row["damping"]))

read_pulley(row) = Pulley(Symbol(row["name"]), Symbol(row["segment_a"]),
    Symbol(row["segment_b"]), parse_dynamics_type(row["type"]);
    efficiency = Float64(row["efficiency"]))

read_tether(row) = Tether(Symbol(row["name"]), Symbol.(row["segments"]);
    start_point = Symbol(row["start_point"]), end_point = Symbol(row["end_point"]))

read_winch(row, set) = Winch(Symbol(row["name"]), Symbol.(row["tethers"]),
    Float64(row["gear_ratio"]), Float64(row["drum_radius"]),
    set.f_coulomb, set.c_vf, set.inertia_total;
    winch_point = Symbol(row["winch_point"]),
    model = named_model(row["model"], AbstractWinchModel))

read_body(row) = Body(Symbol(row["name"]); mass = Float64(row["mass"]),
    inertia_principal = vector3(row["inertia_principal"]),
    pos = vector3(row["pos_cad"]),
    com_offset_b = vector3(row["com_offset_b"]),
    type = parse_dynamics_type(row["type"]), wing = optional_symbol(row["wing"]))

"""
    read_wing(row, station_rows, point_rows, set, vsm_set) -> Body

The aero-carrying body `row` describes. The document has no wing column on a
station, so a wing's stations are those whose points name it.
"""
function read_wing(row, station_rows, point_rows, set, vsm_set)
    parse_dynamics_type(row["type"]) == KINEMATIC && error(
        "Wing $(row["name"]) is KINEMATIC, a PARTICLE_DYNAMICS wing whose body " *
        "frame is fitted to reference points $STRUCTURE_SCHEMA has no column for.")
    on_wing = Set(point["name"] for point in point_rows
                  if point["wing"] == row["name"])
    stations = [Symbol(station["name"]) for station in station_rows
                if all(in(on_wing), station["points"])]
    return VSMWing(Symbol(row["name"]), set, stations, vsm_set;
        transform = 0, dynamics_type = RIGID_DYNAMICS,
        aero = named_model(row["aero"], AbstractAeroModel),
        pos_cad = vector3(row["pos_cad"]), mass = Float64(row["mass"]),
        inertia_diag = KVec3(vector3(row["inertia_principal"])))
end

read_elastic_joint(row) = ElasticJoint(Symbol(row["name"]), Symbol(row["body_a"]),
    Symbol(row["body_b"]); anchor_a = vector3(row["anchor_a"]),
    anchor_b = vector3(row["anchor_b"]),
    stiffness_axial = linear_rigidity(row["stiffness_axial"],
                                      "joint $(row["name"]) axial"),
    stiffness_shear = linear_rigidity(row["stiffness_shear"],
                                      "joint $(row["name"]) shear"),
    stiffness_torsion = linear_rigidity(row["stiffness_torsion"],
                                        "joint $(row["name"]) torsion"),
    stiffness_bending = linear_rigidity(row["stiffness_bending"],
                                        "joint $(row["name"]) bending"),
    damping = Float64(row["damping"]), radius = optional_length(row["radius"]))

read_timoshenko_joint(row) = TimoshenkoJoint(Symbol(row["name"]),
    Symbol(row["body_a"]), Symbol(row["body_b"]);
    anchor_a = vector3(row["anchor_a"]), anchor_b = vector3(row["anchor_b"]),
    EA = linear_rigidity(row["EA"], "joint $(row["name"]) EA"),
    GA = linear_rigidity(row["GA"], "joint $(row["name"]) GA"),
    GJ = linear_rigidity(row["GJ"], "joint $(row["name"]) GJ"),
    EIy = linear_rigidity(row["EIy"], "joint $(row["name"]) EIy"),
    EIz = linear_rigidity(row["EIz"], "joint $(row["name"]) EIz"),
    shear_coeff = Float64(row["shear_coeff"]), damping = Float64(row["damping"]),
    rest_length = Float64(row["rest_length"]),
    radius = optional_length(row["radius"]))

"""Overwrite the mass properties `SystemStructure` derives with the document's."""
function apply_body_row!(body::Body, row)
    body.mass = Float64(row["mass"])
    body.apparent_mass = Float64(row["apparent_mass"])
    body.inertia_principal .= vector3(row["inertia_principal"])
    body.com_offset_b .= vector3(row["com_offset_b"])
    body.pos_cad .= vector3(row["pos_cad"])
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
