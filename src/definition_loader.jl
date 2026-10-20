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

function parse_dynamics_type(text::String)
    text_upper = uppercase(text)
    text_upper == "STATIC" && return STATIC
    text_upper == "DYNAMIC" && return DYNAMIC
    text_upper == "BODY_STATIC" && return BODY_STATIC
    text_upper == "KINEMATIC" && return KINEMATIC
    error("Unknown DynamicsType: $text")
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
                              ignore_l0=false, dynamics_type=nothing, aero_mode=nothing,
                              vsm_set=nothing, wind_mode=ProfileWind(), prn=true)

Build a `SystemStructure` from an authoring YAML, read by KiteGeometry's `load_authoring`
in its design pose and placed by its transforms ([`place!`](@ref)). `set` defaults to the
`base` settings; `dynamics_type` and `aero_mode` override every wing's.
`wind_mode=PerPointWind()` gives each point its own settable wind. `prn=false` silences the
messages loading a structure prints.
"""
function load_sys_struct_from_yaml(yaml_path::AbstractString; system_name="from_yaml",
        set::Union{Nothing, Settings}=nothing, ignore_l0::Bool=false,
        dynamics_type::Union{Nothing, WingType}=nothing,
        aero_mode::Union{Nothing, AbstractAeroModel}=nothing,
        vsm_set::Union{Nothing, VortexStepMethod.VSMSettings}=nothing,
        wind_mode::WindMode=ProfileWind(), prn::Bool=true)
    isnothing(set) && (set = load_settings("base"))
    definition, transforms = load_authoring(yaml_path; set, name=system_name,
                                            place=false)
    return SystemStructure(definition, transforms, set; placed=false, ignore_l0,
                           dynamics_type, aero_mode, vsm_set, wind_mode, prn)
end

"""
    SystemStructure(definition::SystemDefinition, transforms, set; placed=true,
                    ignore_l0=false, dynamics_type=nothing, aero_mode=nothing,
                    vsm_set=nothing, wind_mode=ProfileWind(), prn=true)

The structure of a KiteGeometry `definition` and its `transforms`, in the pose the
definition gives, or, with `placed=false`, placed from it ([`place!`](@ref)). Its
components are read from the definition's rows and the models SymbolicAWEModels registers
for them ([`register_definition_models`](@ref)).
"""
function SystemStructure(definition::SystemDefinition, transforms, set::Settings;
        placed=true, ignore_l0=false, dynamics_type=nothing, aero_mode=nothing,
        vsm_set=nothing, wind_mode::WindMode=ProfileWind(), prn=true)
    point_transforms = member_transforms(transforms, :points, length(definition.points))
    body_transforms = member_transforms(transforms, :bodies, length(definition.bodies))
    n_wings = length(definition.wings)
    stations = Station[definition_station(definition, station)
                       for station in definition.spanwise_stations]
    wings = Body[definition_wing(definition, w, set, stations, body_transforms[w];
                                 dynamics_type, aero_mode, vsm_set)
                 for w in 1:n_wings]
    return SystemStructure(definition.metadata.name, set;
        points=Point[definition_point(definition, point, point_transforms[i])
                     for (i, point) in enumerate(definition.points)],
        stations,
        segments=Segment[definition_segment(definition, segment, set)
                         for segment in definition.segments],
        pulleys=Pulley[definition_pulley(definition, pulley)
                       for pulley in definition.pulleys],
        tethers=Tether[definition_tether(definition, tether)
                       for tether in definition.tethers],
        winches=Winch[definition_winch(definition, winch, set)
                      for winch in definition.winches],
        wings,
        transforms=Transform[definition_transform(definition, transform, transforms)
                             for transform in transforms],
        bodies=Body[definition_body(definition.bodies[b], body_transforms[b])
                    for b in n_wings+1:length(definition.bodies)],
        tubes=Tube[definition_tube(definition, tube) for tube in definition.tubes],
        placed, ignore_l0, vsm_set, wind_mode, prn, definition)
end

"""The name of the transform moving each of `count` members of `block`, `nothing` for
none."""
function member_transforms(transforms, block, count)
    member = Vector{Union{Nothing, Int, Symbol}}(nothing, count)
    for transform in transforms
        member[getfield(transform, block)] .= Ref(row_name(transform))
    end
    return member
end

"""A row's name as a `Symbol`, or its index where the authoring YAML gave it none."""
row_name(row) = something(tryparse(Int, String(row.name)), row.name)

"""The name of row `i` of `rows`, as [`row_name`](@ref) gives it, `nothing` for
`nothing`."""
name_of(rows, i) = isnothing(i) ? nothing : row_name(rows[i])

"""A reference an authoring YAML cell gives: an index, or a name as a `Symbol`."""
cell_ref(value::Integer) = Int(value)
cell_ref(value::Union{AbstractString, Symbol}) = Symbol(value)

"""The weighted reference points a cell gives: one reference, a list averaged, or
`[reference, weight]` pairs."""
ref_points_spec(value) = value isa AbstractVector ?
    [entry isa AbstractVector ? (cell_ref(entry[1]), Float64(entry[2])) : cell_ref(entry)
     for entry in value] : cell_ref(value)

"""Keyword arguments from `pairs` whose value is set."""
set_kwargs(pairs...) = (; (key => value for (key, value) in pairs if !isnothing(value))...)

function definition_point(definition, point, transform)
    type = isnothing(point.type) ?
        (point.fixed ? STATIC : isnothing(point.body) ? DYNAMIC : BODY_STATIC) :
        parse_dynamics_type(point.type)
    return Point(row_name(point), point.pos_ENU, type;
        set_kwargs(:wing => maybe(cell_ref, point.wing_idx), :transform => transform,
                   :body => point_body(definition, point),
                   :anchor_b => point.anchor_b,
                   :tube => maybe(cell_ref, point.tube),
                   :extra_mass => point.extra_mass,
                   :body_frame_damping => point.body_frame_damping,
                   :world_frame_damping => point.world_frame_damping,
                   :area => point.area, :drag_coeff => point.drag_coeff,
                   :fix_sphere => point.fix_sphere, :fix_static => point.fix_static)...)
end

"""The body `point` rides: none for a tube rider whose body is the one KiteGeometry gives
a `BODY_STATIC` point naming its wing and no body, which is that wing's."""
function point_body(definition, point)
    isnothing(point.body) && return nothing
    wing = maybe(cell_ref, point.wing_idx)
    wing_body = point.body <= length(definition.wings) &&
        wing in (point.body, name_of(definition.wings, point.body))
    return !isnothing(point.tube) && wing_body ? nothing :
        name_of(definition.bodies, point.body)
end

"""`f(value)`, or `nothing` for `nothing`."""
maybe(f, value) = isnothing(value) ? nothing : f(value)

function definition_station(definition, station)
    isnothing(station.type) && error("station $(station.name) needs its `type`")
    return Station(row_name(station), name_of.((definition.points,), station.points),
        parse_dynamics_type(station.type),
        something(station.moment_frac, 0.0);
        set_kwargs(:damping => station.damping, :stiffness => station.stiffness,
                   :wing => name_of(definition.wings, station.wing),
                   :bodies => maybe(refs -> cell_ref.(refs), station.bodies),
                   :flap_bodies => maybe(refs -> cell_ref.(refs), station.flap_bodies),
                   :flap_points => maybe(refs -> cell_ref.(refs), station.flap_points),
                   :flap_axis => station.flap_axis)...)
end

function definition_segment(definition, segment, set)
    material = definition.materials[segment.material]
    material.law == "linear" || error("segment $(segment.name): material law " *
        "$(material.law) has no force law here; give its unit_stiffness in Julia")
    spring = segment.model
    points = name_of.((definition.points,), segment.points)
    return Segment(row_name(segment), set, points...;
        l0=segment.rest_length, diameter_mm=1000 * segment.diameter,
        youngs_modulus=material.parameters["youngs_modulus"],
        density=material.parameters["density"], unit_damping=spring.unit_damping,
        spring.compression_frac, spring.compression_damping_frac)
end

function definition_pulley(definition, pulley)
    isnothing(pulley.type) && error("pulley $(pulley.name) needs its `type`")
    return Pulley(row_name(pulley), name_of.((definition.segments,), pulley.segments)...,
        parse_dynamics_type(pulley.type);
        set_kwargs(:efficiency => pulley.efficiency, :damping => pulley.damping,
                   :brake => pulley.brake,
                   :friction_epsilon => pulley.friction_epsilon)...)
end

definition_tether(definition, tether) = Tether(row_name(tether),
    name_of.((definition.segments,), tether.segments), tether.init_stretched_length;
    start_point=name_of(definition.points, tether.start_point),
    end_point=name_of(definition.points, tether.end_point),
    tether_force=tether.init_tether_force, stretch_frac=tether.init_stretch_frac)

definition_winch(definition, winch, set) = Winch(row_name(winch), set,
    name_of.((definition.tethers,), winch.tethers);
    winch_point=name_of(definition.points, winch.winch_point),
    set_kwargs(:init_vel => winch.init_vel, :brake => winch.brake,
               :speed_controlled => winch.speed_controlled,
               :friction_epsilon => winch.friction_epsilon)...)

function definition_transform(definition, transform, transforms)
    body = transform.body
    wing = !isnothing(body) && body <= length(definition.wings) ?
        name_of(definition.wings, body) : nothing
    rot_point = isnothing(wing) ?
        name_of(definition.points, first(only(transform.rot_point))) : nothing
    return Transform(row_name(transform), transform.elevation, transform.azimuth,
        transform.heading; base_point=name_of(definition.points, transform.base_point),
        transform.base_pos, wing, rot_point,
        base_transform=name_of(transforms, transform.base_transform),
        transform.elevation_vel, transform.azimuth_vel, transform.turn_rate)
end

function definition_body(body, transform)
    type = maybe(parse_dynamics_type, body.type)
    return Body(row_name(body); extra_mass=something(body.extra_mass, 0.0),
        inertia=Matrix{SimFloat}(something(body.extra_inertia_KA, zeros(3, 3))),
        pos=body.pos_ENU, Q_b_to_w=Vector{SimFloat}(body.Q_KA2ENU),
        set_kwargs(:type => type, :transform => transform,
                   :wing => maybe(cell_ref, body.wing),
                   :angular_damping => body.angular_damping,
                   :world_frame_damping => body.world_frame_damping,
                   :body_frame_damping => body.body_frame_damping,
                   :fix_sphere => body.fix_sphere, :fix_static => body.fix_static,
                   :ext_force_w => body.ext_force_w, :ext_force_b => body.ext_force_b,
                   :ext_moment_b => body.ext_moment_b,
                   :principal_frame_method =>
                       maybe(parse_principal_frame_method, body.principal_frame_method))...)
end

function definition_tube(definition, tube)
    tube_model = tube.model
    columns = (; (field => getfield(tube_model, field)
                  for field in fieldnames(typeof(tube_model))
                  if !(field in (:anchor_a, :anchor_b)) &&
                     !isnothing(getfield(tube_model, field)))...)
    model = tube_model isa ElasticTubeModel ? ElasticTube(; columns...) :
        TimoshenkoTube(; columns...)
    return Tube(row_name(tube), name_of.((definition.bodies,), tube.bodies)...;
        tube.diameter, tube.pressure,
        law=definition.materials[tube.material].law, model,
        set_kwargs(:anchor_a => tube_model.anchor_a,
                   :anchor_b => tube_model.anchor_b)...)
end

"""
    definition_wing(definition, w, set, stations, transform; dynamics_type, aero_mode,
                    vsm_set)

Wing `w` of `definition`, its mesh authored in the wing's KA frame. `dynamics_type` and
`aero_mode` override the wing row's; `stations` collects the stations a flat-plate wing
adds.
"""
function definition_wing(definition, w, set, stations, transform; dynamics_type,
                         aero_mode, vsm_set)
    wing = definition.wings[w]
    dynamics = something(dynamics_type, maybe(parse_wing_type, wing.dynamics_type))
    mode = something(aero_mode, maybe(parse_aero_mode, wing.aero_mode),
        dynamics == RIGID_DYNAMICS ? AeroLinearized() : AeroDirect())
    station_refs = [cell_ref(station) for station in something(wing.stations, [])]
    built = load_wing(mode, definition, w, set, dynamics, vsm_set, station_refs,
                      stations; transform)
    built.mesh_frame_pos .= wing.pos_ENU
    built.R_mesh_to_w .= quaternion_to_rotation_matrix(wing.Q_KA2ENU)
    return built
end

"""
    load_wing(mode::AbstractAeroModel, definition, w, set, dynamics_type, vsm_set,
              station_refs, stations; transform)

Build wing `w` of `definition` for its aero `mode`: a [`VSMWing`](@ref) by default.
Add a method to build a wing for a custom aero mode.
"""
function load_wing(mode::AbstractAeroModel, definition, w, set, dynamics_type, vsm_set,
                   station_refs, stations; transform)
    wing = definition.wings[w]
    pos_ENU = dynamics_type == PARTICLE_DYNAMICS || isnothing(wing.pos_cad) ? nothing :
        wing.pos_ENU + quaternion_to_rotation_matrix(wing.Q_KA2ENU) * wing.pos_cad
    aero_z_offset = dynamics_type == RIGID_DYNAMICS ? wing.aero_z_offset : nothing
    return VSMWing(row_name(wing), set, station_refs, vsm_set;
        transform, dynamics_type, aero=mode, pos_ENU,
        extra_mass=definition.bodies[w].extra_mass,
        set_kwargs(:angular_damping => wing.angular_damping,
                   :z_ref_points => maybe(ref_pair_spec, wing.z_ref_points),
                   :y_ref_points => maybe(ref_pair_spec, wing.y_ref_points),
                   :origin => maybe(ref_points_spec, wing.origin_idx),
                   :aero_scale_chord => wing.aero_scale_chord,
                   :aero_z_offset => aero_z_offset,
                   :principal_frame_method =>
                       maybe(parse_principal_frame_method, wing.principal_frame_method),
                   :com => wing.com,
                   :unit_inertia => maybe(unit_inertia_spec, wing.unit_inertia))...)
end

"""The z or y reference points of a wing: two weighted point groups."""
ref_pair_spec(value) = (ref_points_spec(value[1]), ref_points_spec(value[2]))

"""A wing's per-unit-mass inertia as a 3×3 matrix or `[Ixx, Iyy, Izz, Ixy, Ixz, Iyz]`."""
unit_inertia_spec(value) = first(value) isa AbstractVector ?
    permutedims(reduce(hcat, Vector{Float64}.(value))) : Float64.(value)

# ==================== STRUCTURE DOCUMENTS ==================== #

"""
    structure_document(sys::SystemStructure)

The awesIO structure document of `sys`'s definition in its initial pose, refusing one
whose segments carry a nonlinear stiffness law, which the schema cannot hold yet
(1-Bart-1/awesIO#1).
"""
function KiteGeometry.structure_document(sys::SystemStructure)
    for material in sys.definition.materials
        material.law == "unit_stiffness" && error("material $(material.name) is a " *
            "nonlinear stiffness law, which a structure document cannot hold yet " *
            "(1-Bart-1/awesIO#1)")
    end
    return structure_document(sys.definition)
end

"""
    sys_struct_from_document(document::AbstractDict; set=nothing, kwargs...)

The `SystemStructure` of a parsed structure document, in the initial pose it gives. `set`
falls back to the `base` settings; `kwargs` reach
[`SystemStructure(::SystemDefinition, transforms, set)`](@ref).
"""
function sys_struct_from_document(document::AbstractDict; set=nothing, kwargs...)
    isnothing(set) && (set = load_settings("base"))
    return SystemStructure(SystemDefinition(document), KiteGeometry.Transform[], set;
                           kwargs...)
end

"""
    save_structure_document(path, sys::SystemStructure)

Write [`structure_document`](@ref)`(sys)` to `path`, as YAML for a `.yml` or `.yaml`
extension and as JSON for `.json`.
"""
function save_structure_document(path::AbstractString, sys::SystemStructure)
    document = structure_document(sys)
    write(path, document_format(path) == :json ? JSON.json(document, 2) :
                YAML.write(document))
    return path
end

"""
    load_structure_document(path; kwargs...)

The `SystemStructure` of the structure document at `path`, YAML or JSON by extension;
`kwargs` reach [`sys_struct_from_document`](@ref).
"""
function load_structure_document(path::AbstractString; kwargs...)
    text = read(path, String)
    document = document_format(path) == :json ? JSON.parse(text) :
        YAML.load(text; dicttype=OrderedDict{String, Any})
    return sys_struct_from_document(document; kwargs...)
end

"""The encoding `path`'s extension asks for: `:yaml` or `:json`."""
function document_format(path::AbstractString)
    extension = lowercase(splitext(path)[2])
    extension in (".yml", ".yaml") && return :yaml
    extension == ".json" && return :json
    error("$path: a structure document is written as .yml, .yaml or .json.")
end

"""
    save_log(logger::Logger, sys::SystemStructure, name="sim_log"; path="")

Save `logger` as an uncompressed .arrow file whose metadata holds the structure document
of `sys` under the key `topology`.
"""
KiteUtils.save_log(logger::Logger, sys::SystemStructure, name="sim_log"; path="") =
    save_log(logger, sys.definition, name, false; path)
