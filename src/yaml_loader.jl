# Copyright (c) 2025 Bart van de Lint, Jelle Poland
# SPDX-License-Identifier: LGPL-3.0-only

using YAML
using Logging
using LinearAlgebra
using VortexStepMethod, KiteUtils

# ------------------ helpers ------------------

"""
    get_field_or_nothing(::Type{T}, row::NamedTuple,
                         field::Symbol) where T

Convert field to type T if present, otherwise return nothing.

# Examples
```julia
get_field_or_nothing(Int64, row, :idx)  # -> Int64 or nothing
get_field_or_nothing(Tuple{Int64,Int64}, row, :pair)
    # -> (Int64, Int64) or nothing
```
"""
function get_field_or_nothing(::Type{T}, row::NamedTuple,
                               field::Symbol) where T
    if !haskey(row, field) || isnothing(row[field])
        return nothing
    end
    return convert_to_type(T, row[field])
end

"""
    convert_to_type(::Type{T}, value) where T

Convert value to type T. Handles special cases like Tuples.
"""
convert_to_type(::Type{T}, value) where T = T(value)

# Special handling for Tuple types
function convert_to_type(
        ::Type{Tuple{T,T}}, value) where T
    values = Vector{Int}(value)
    return (T(values[1]), T(values[2]))
end

"""
    resolve_references(row::NamedTuple, property_tables::Dict{String, Dict{String, NamedTuple}})

Resolve string references in a row by looking them up in property tables.
If a field value is a String, check if it exists as a key in any property table.
If found, merge those properties into the current row (current row takes precedence).
"""
function resolve_references(row::NamedTuple, property_tables::Dict{String, Dict{String, NamedTuple}})
    resolved = Dict{Symbol, Any}(pairs(row))

    for (_, value) in pairs(row)
        # Check if this is a string reference (unquoted strings in YAML)
        if value isa String
            # Try to find this reference in any property table
            for (_, lookup_dict) in property_tables
                if haskey(lookup_dict, value)
                    # Found! Merge these properties (current row takes precedence)
                    ref_props = lookup_dict[value]
                    for (key, field_value) in pairs(ref_props)
                        # Skip 'name': it identifies the referenced item, not a property.
                        key === :name && continue
                        if !haskey(resolved, key) || resolved[key] === nothing
                            resolved[key] = field_value
                        end
                    end
                    break
                end
            end
        end
    end

    return NamedTuple(resolved)
end

"""
    calculate_derived_properties!(props::Dict{Symbol, Any})

Calculate derived properties like `unit_stiffness` and `unit_damping` from material properties.
Modifies props in-place.
"""
function calculate_derived_properties!(props::Dict{Symbol, Any})
    # Store EA; the spring k is computed later as EA / len in generate_system.jl.
    if haskey(props, :youngs_modulus) && haskey(props, :diameter_mm)
        # Check if we need to calculate (missing, nothing, or is a string reference)
        need_calculation = !haskey(props, :unit_stiffness) ||
                          props[:unit_stiffness] === nothing ||
                          props[:unit_stiffness] isa String

        if need_calculation
            diameter_m = Float64(props[:diameter_mm]) / 1000.0  # mm to m
            area = π * (diameter_m / 2)^2
            youngs_modulus = Float64(props[:youngs_modulus])
            props[:unit_stiffness] = youngs_modulus * area
        end
    end

    # Calculate unit_damping from damping coefficient if missing or is a string
    if haskey(props, :damping_per_stiffness) && haskey(props, :unit_stiffness)
        # Only calculate if unit_stiffness is now a number
        if props[:unit_stiffness] isa Number
            need_damping_calc = !haskey(props, :unit_damping) ||
                               props[:unit_damping] === nothing ||
                               props[:unit_damping] isa String

            if need_damping_calc
                props[:unit_damping] = Float64(props[:damping_per_stiffness]) * Float64(props[:unit_stiffness])
            end
        end
    end

    # Set default unit_damping if still missing
    if !haskey(props, :unit_damping) || props[:unit_damping] === nothing || props[:unit_damping] isa String
        props[:unit_damping] = 0.0
    end
end

function parse_table(table)::Vector{NamedTuple}
    haskey(table, "data") || throw(ArgumentError("table is missing `data`"))

    rows = table["data"]
    (isnothing(rows) || isempty(rows)) && return NamedTuple[]

    # Dict first row => dict format; Array first row => header format.
    first_row = first(rows)

    if first_row isa AbstractDict
        out = NamedTuple[]
        for row in rows
            named_row = NamedTuple{Tuple(Symbol.(keys(row)))}(
                Tuple(values(row)))
            push!(out, named_row)
        end
        return out
    else
        # Array format: requires headers
        haskey(table, "headers") ||
            throw(ArgumentError(
                "table with array rows requires `headers`"))
        headers = String.(table["headers"])

        out = NamedTuple[]
        for (k, row) in enumerate(rows)
            # skip empty or comment rows
            if isempty(row) ||
               (isa(row[1], String) && startswith(row[1], "#"))
                continue
            end
            # allow missing trailing columns (fill with nothing)
            if length(row) < length(headers)
                row = vcat(row, fill(nothing,
                    length(headers) - length(row)))
            end
            if length(row) > length(headers)
                @warn "Skipping row $k: has $(length(row)) " *
                      "values, expected $(length(headers))."
                continue
            end
            named_row = NamedTuple{Tuple(Symbol.(headers))}(Tuple(row))
            push!(out, named_row)
        end
        return out
    end
end

"""
    extract_args(row, args_spec, mappings)

Extract positional constructor arguments from a YAML row.

For each name in `args_spec`, this helper first checks for a
mapping in `mappings`, then falls back to `row[arg_name]`.
Throws an error if a required argument is missing.
"""
function extract_args(row, args_spec, mappings)
    args = []
    for arg_name in args_spec
        if haskey(mappings, arg_name)
            push!(args, mappings[arg_name](row))
        elseif haskey(row, arg_name)
            push!(args, row[arg_name])
        else
            error("Missing required arg $arg_name")
        end
    end
    return args
end

"""
    call_yaml_constructor(Constructor, row::NamedTuple,
        args_spec, kwargs_spec; mappings=Dict())

Generic YAML-to-constructor caller. Extracts positional
args and kwargs from YAML row and calls constructor.

# Arguments
- `Constructor`: Constructor function to call
- `row::NamedTuple`: Parsed YAML row
- `args_spec::Vector{Symbol}`: Names for positional args
- `kwargs_spec::Vector{Symbol}`: Names for kwargs

# Keyword Arguments
- `mappings::Dict{Symbol, Function}`: Mapping functions
  that take the row and return the arg value

# Example
```julia
row = (idx=1, x=0.0, y=0.0, z=0.0, type="STATIC")
point = call_yaml_constructor(Point, row,
    [:idx, :pos_cad, :type],  # positional args
    [:extra_mass, :wing_idx];       # kwargs
    mappings=Dict(
        :pos_cad => r -> [Float64(r.x),
            Float64(r.y), Float64(r.z)],
        :type => r -> parse_dynamics_type(
            String(r.type))
    ))
```
"""
function call_yaml_constructor(
        Constructor,
        row::NamedTuple,
        args_spec::Vector{Symbol},
    ::Vector{Union{}};
        mappings::Dict{Symbol, <:Function}=
            Dict{Symbol, Function}())
    args = extract_args(row, args_spec, mappings)
    return Constructor(args...)
end

function call_yaml_constructor(
        Constructor,
        row::NamedTuple,
        args_spec::Vector{Symbol},
        kwargs_spec::Vector;
        mappings::Dict{Symbol, <:Function}=
            Dict{Symbol, Function}())
    args = extract_args(row, args_spec, mappings)

    # Extract keyword arguments (only if present)
    kwargs = Dict{Symbol, Any}()
    for kwarg_name in kwargs_spec
        if haskey(mappings, kwarg_name)
            kwargs[kwarg_name] = mappings[kwarg_name](row)
        elseif haskey(row, kwarg_name)
            val = row[kwarg_name]
            # Skip nothing values (Julia nothing or YAML "nothing" string)
            if !isnothing(val) && val != "nothing"
                kwargs[kwarg_name] = val
            end
        end
    end

    if isempty(kwargs)
        return Constructor(args...)
    else
        return Constructor(args...; kwargs...)
    end
end

function parse_dynamics_type(text::String)
    text_upper = uppercase(text)
    text_upper == "STATIC" && return STATIC
    text_upper == "DYNAMIC" && return DYNAMIC
    text_upper == "WING" && error(
        "DynamicsType `WING` was removed: use `BODY_STATIC` (rigid wing, with " *
        "`body_idx` = the wing) or `DYNAMIC` (particle wing), and make the point " *
        "a member of one of the wing's twist_surfaces.")
    text_upper == "BODY_STATIC" && return BODY_STATIC
    text_upper == "KINEMATIC" && return KINEMATIC
    error("Unknown DynamicsType: $text")
end

function parse_segment_type(text::String)
    text_upper = uppercase(text)
    text_upper == "POWER_LINE" && return POWER_LINE
    text_upper == "STEERING_LINE" && return STEERING_LINE
    text_upper == "BRIDLE" && return BRIDLE
    error("Unknown SegmentType: $text")
end

function parse_wing_type(text::String)
    text_upper = uppercase(text)
    text_upper == "PARTICLE_DYNAMICS" && return PARTICLE_DYNAMICS
    text_upper == "RIGID_DYNAMICS" && return RIGID_DYNAMICS
    if text_upper == "REFINE"
        @warn "WingType \"$text\" is deprecated; use \"PARTICLE_DYNAMICS\" instead."
        return PARTICLE_DYNAMICS
    end
    if text_upper == "QUATERNION"
        @warn "WingType \"$text\" is deprecated; use \"RIGID_DYNAMICS\" instead."
        return RIGID_DYNAMICS
    end
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
    load_wing(mode::AbstractAeroModel, row, idx, data, set, wing_type, vsm_set,
              yaml_to_ref, yaml_parse_ref_points, yaml_parse_origin, twist_surfaces)

Build a wing from a parsed YAML `row`, dispatched on its aero `mode`. The default
(VSM-backed modes) builds a [`VSMWing`](@ref); [`AeroPlate`](@ref) builds a
flat-plate wing via [`load_plate_wing`](@ref). Add a method to load a wing for a
custom aero mode.
"""
function load_wing(mode::AbstractAeroModel, row, idx, data, set, wing_type,
                   vsm_set, yaml_to_ref, yaml_parse_ref_points,
                   yaml_parse_origin, twist_surfaces)
    # PARTICLE derives pos_cad from the origin point during resolution; RIGID
    # reads it from the row when given.
    pos_cad = wing_type == PARTICLE_DYNAMICS ? (row -> nothing) :
        (row -> yaml_vec3(row, :pos_cad))
    # aero_z_offset only applies to RIGID wings.
    kwargs_spec = [:transform, :y_damping, :angular_damping, :dynamics_type,
        :aero, :z_ref_points, :y_ref_points, :origin, :pos_cad,
        :aero_scale_chord, :principal_frame_method,
        :mass, :com, :unit_inertia]
    wing_type == RIGID_DYNAMICS && push!(kwargs_spec, :aero_z_offset)
    return call_yaml_constructor(VSMWing, row,
        [:name, :set, :twist_surfaces, :vsm_set], kwargs_spec;
        mappings=Dict(
            :set => row -> set,
            :aero => row -> mode,
            :twist_surfaces => row -> hasfield(typeof(row), :twist_surfaces) &&
                !isnothing(row.twist_surfaces) ?
                [yaml_to_ref(ref) for ref in row.twist_surfaces] : [],
            :vsm_set => row -> vsm_set,
            :dynamics_type => row -> wing_type,
            :name => row -> yaml_row_name(row, idx),
            :transform => row -> yaml_ref_field(row, :transform_idx, yaml_to_ref),
            :pos_cad => pos_cad,
            :z_ref_points => row -> yaml_parse_ref_points(row, :z_ref_points),
            :y_ref_points => row -> yaml_parse_ref_points(row, :y_ref_points),
            :origin => row -> yaml_parse_origin(row, :origin_idx),
            :aero_scale_chord => row ->
                hasfield(typeof(row), :aero_scale_chord) &&
                !isnothing(row.aero_scale_chord) ?
                    float(row.aero_scale_chord) : 0.0,
            :principal_frame_method => row ->
                hasfield(typeof(row), :principal_frame_method) &&
                !isnothing(row.principal_frame_method) ?
                    parse_principal_frame_method(String(row.principal_frame_method)) :
                    EIGEN_DECOMP,
            :mass => row -> yaml_float(row, :mass),
            :com => row -> yaml_vec3(row, :com),
            :unit_inertia => row -> begin
                value = yaml_field(row, :unit_inertia)
                isnothing(value) ? nothing : Float64.(value)
            end
        ))
end

"""
    parse_tether_init(row, tether_name)
        -> (stretched_length, tether_force, stretch_frac)

Read `init_stretched_length`, `init_tether_force` and
`init_stretch_frac` from a tether YAML row (each `nothing` if
absent). Errors if the deprecated `init_unstretched_length` field is
present — the unstretched rest length is now derived from the placed
stretched length with `init_tether_force` or `init_stretch_frac`.
"""
function parse_tether_init(row, tether_name)
    if !isnothing(yaml_field(row, :init_unstretched_length))
        error("Tether $tether_name: init_unstretched_length is " *
              "deprecated; it is derived from the placed " *
              "init_stretched_length with init_tether_force or " *
              "init_stretch_frac.")
    end
    stretched_length = yaml_float(row, :init_stretched_length)
    tether_force = yaml_float(row, :init_tether_force)
    stretch_frac = yaml_float(row, :init_stretch_frac)
    return stretched_length, tether_force, stretch_frac
end

function yaml_field(row, field)
    hasfield(typeof(row), field) || return nothing
    getfield(row, field)
end

function yaml_float(row, field)
    value = yaml_field(row, field)
    isnothing(value) ? nothing : Float64(value)
end

"""
    yaml_row_name(row, i)

Name for a YAML row: its `name` field (as a `Symbol`) when present, else the
1-based row index `i`. Shared by the body/joint loaders so components can be
referenced either by name or by position.
"""
yaml_row_name(row, i) =
    (hasfield(typeof(row), :name) && !isnothing(row.name)) ? Symbol(row.name) : i

"""
    yaml_block_empty(data, key) -> Bool

`true` when the top-level YAML `key` is absent or holds no `data` rows.
"""
function yaml_block_empty(data, key)
    haskey(data, key) || return true
    table = data[key]
    (haskey(table, "data") && !isnothing(table["data"]) &&
     !isempty(table["data"])) || return true
    return false
end

"""
    yaml_vec3(row, field) -> KVec3 or nothing

Read a 3-element vector field, or `nothing` when absent.
"""
function yaml_vec3(row, field)
    value = yaml_field(row, field)
    isnothing(value) ? nothing : KVec3(value...)
end

"""
    yaml_matrix3(value) -> Matrix{SimFloat}

Convert a nested `[[…],[…],[…]]` YAML value into a 3×3 matrix (rows preserved).
"""
yaml_matrix3(value) =
    Matrix{SimFloat}(permutedims(reduce(hcat,
        (Vector{SimFloat}(row) for row in value))))

"""
    yaml_ref_field(row, field, to_ref) -> reference or nothing

Resolve an optional reference column `field` (e.g. `:transform_idx`) through
`to_ref`, or `nothing` when the field is absent or empty.
"""
yaml_ref_field(row, field, to_ref) =
    (hasfield(typeof(row), field) && !isnothing(getfield(row, field))) ?
        to_ref(getfield(row, field)) : nothing

"""
    load_property_table(data, key) -> Dict{String, NamedTuple}

Index a reference table (`materials`, `elements`, `segment_properties`) by its
`name` column, so [`resolve_references`](@ref) can merge shared properties.
"""
function load_property_table(data, key)
    table = Dict{String, NamedTuple}()
    for row in parse_table(data[key])
        table[String(row.name)] = row
    end
    return table
end

"""
    load_yaml_bodies(data, yaml_to_ref) -> Vector{Body}

Build the plain rigid [`Body`](@ref)s from a `bodies` YAML block (empty when the
block is absent). Field names mirror the [`Body`](@ref) constructor: required
`name`, `mass`, `pos`, and one of `inertia_principal` (3-vector) or `inertia`
(3×3); optional `type` (`DYNAMIC`/`STATIC`), `transform_idx`, `vel`, `Q_b_to_w`
(4-vector), `omega_b`, `com_offset_b`, `damping` (scalar or 3-vector),
`fix_sphere`, `ext_force_w`, `ext_force_b`, `ext_moment_b`,
`principal_frame_method`.
"""
function load_yaml_bodies(data, yaml_to_ref)
    bodies = Body[]
    yaml_block_empty(data, "bodies") && return bodies
    for (i, row) in enumerate(parse_table(data["bodies"]))
        name = yaml_row_name(row, i)
        mass = yaml_float(row, :mass)
        isnothing(mass) && error("Body $name: missing required `mass`.")
        pos = yaml_vec3(row, :pos)
        isnothing(pos) && error("Body $name: missing required `pos`.")
        kwargs = Dict{Symbol, Any}(:mass => mass, :pos => pos)
        inertia_principal = yaml_vec3(row, :inertia_principal)
        isnothing(inertia_principal) ||
            (kwargs[:inertia_principal] = inertia_principal)
        inertia = yaml_field(row, :inertia)
        isnothing(inertia) || (kwargs[:inertia] = yaml_matrix3(inertia))
        vel = yaml_vec3(row, :vel); isnothing(vel) || (kwargs[:vel] = vel)
        quat = yaml_field(row, :Q_b_to_w)
        isnothing(quat) || (kwargs[:Q_b_to_w] = Vector{SimFloat}(quat))
        omega = yaml_vec3(row, :omega_b)
        isnothing(omega) || (kwargs[:ω_b] = omega)
        com_offset = yaml_vec3(row, :com_offset_b)
        isnothing(com_offset) || (kwargs[:com_offset_b] = com_offset)
        damping = yaml_field(row, :damping)
        isnothing(damping) ||
            (kwargs[:damping] = damping isa Real ? SimFloat(damping) :
                KVec3(damping...))
        fix_sphere = yaml_field(row, :fix_sphere)
        isnothing(fix_sphere) || (kwargs[:fix_sphere] = Bool(fix_sphere))
        body_type = yaml_field(row, :type)
        isnothing(body_type) ||
            (kwargs[:type] = parse_dynamics_type(String(body_type)))
        transform = yaml_ref_field(row, :transform_idx, yaml_to_ref)
        isnothing(transform) || (kwargs[:transform] = transform)
        for (field, key) in ((:ext_force_w, :ext_force_w),
                             (:ext_force_b, :ext_force_b),
                             (:ext_moment_b, :ext_moment_b))
            vec = yaml_vec3(row, field)
            isnothing(vec) || (kwargs[key] = vec)
        end
        pfm = yaml_field(row, :principal_frame_method)
        isnothing(pfm) || (kwargs[:principal_frame_method] =
            parse_principal_frame_method(String(pfm)))
        push!(bodies, Body(name; kwargs...))
    end
    return bodies
end

"""
    load_yaml_joints(::Type{Joint}, data, key, required, optional; yaml_to_ref)

Build two-body joints of type `Joint` from the `key` YAML block (empty when the
block is absent). Every row needs `body_a`, `body_b` and each `required` scalar;
`anchor_a`/`anchor_b` (3-vectors) and each `optional` scalar are read when
present. Shared by the [`TimoshenkoJoint`](@ref) and [`ElasticJoint`](@ref)
loaders — scalar fields from YAML are linear; callable/nonlinear laws are
supplied programmatically.
"""
function load_yaml_joints(::Type{Joint}, data, key, required, optional;
                          yaml_to_ref) where Joint
    joints = Joint[]
    yaml_block_empty(data, key) && return joints
    for (i, row) in enumerate(parse_table(data[key]))
        name = yaml_row_name(row, i)
        body_a = yaml_to_ref(yaml_field(row, :body_a))
        body_b = yaml_to_ref(yaml_field(row, :body_b))
        (isnothing(body_a) || isnothing(body_b)) &&
            error("$(nameof(Joint)) $name: requires `body_a` and `body_b`.")
        kwargs = Dict{Symbol, Any}()
        for field in required
            value = yaml_float(row, field)
            isnothing(value) &&
                error("$(nameof(Joint)) $name: missing `$field`.")
            kwargs[field] = value
        end
        for field in (:anchor_a, :anchor_b)
            anchor = yaml_vec3(row, field)
            isnothing(anchor) || (kwargs[field] = anchor)
        end
        for field in optional
            value = yaml_float(row, field)
            isnothing(value) || (kwargs[field] = value)
        end
        push!(joints, Joint(name, body_a, body_b; kwargs...))
    end
    return joints
end

"""
    load_sys_struct_from_yaml(yaml_path; system_name, set, ...)

Build a `SystemStructure` from a component-based structural
YAML file. See source for full documentation of expected blocks.
"""
function load_sys_struct_from_yaml(yaml_path::AbstractString; system_name="from_yaml", set::Union{Nothing,Settings}=nothing, ignore_l0::Bool=false, dynamics_type::Union{Nothing,WingType}=nothing, aero_mode::Union{Nothing,AbstractAeroModel}=nothing, vsm_set::Union{Nothing,VortexStepMethod.VSMSettings}=nothing, wing_type::Union{Nothing,WingType}=nothing)
    if !isnothing(wing_type)
        if !isnothing(dynamics_type)
            error("Cannot specify both `wing_type` and `dynamics_type`; `wing_type` is deprecated, use `dynamics_type`.")
        end
        @warn "Keyword argument `wing_type` is deprecated; use `dynamics_type` instead."
        dynamics_type = wing_type
    end
    data = YAML.load_file(yaml_path)

    # Use provided settings or fall back to base settings
    local resolved_set = (set === nothing ? load_settings("base") : set)

    # Note: Name resolution is now handled by SystemStructure.assign_indices_and_resolve!

    # Helper to convert raw reference to proper type (Int or Symbol).
    yaml_to_ref = function (val)
        isnothing(val) && return nothing
        if val isa Integer
            return Int(val)
        elseif val isa String
            val == "nothing" && return nothing
            return Symbol(val)
        elseif val isa Symbol
            return val
        else
            return Int(val)
        end
    end

    # Ref spec to WeightedRefPoints: Symbol, Int, [a,b] (avg), or [[a,w],...] (weights).
    yaml_convert_ref = function (value)
        if value isa Vector && !isempty(value) && value[1] isa Vector
            return map(value) do entry
                if !(entry isa Vector) || length(entry) != 2
                    throw(ArgumentError("Invalid weighted reference point entry $(repr(entry)); expected format [[id, weight], ...]"))
                end
                if !(entry[2] isa Number)
                    throw(ArgumentError("Invalid weighted reference point weight $(repr(entry[2])) in entry $(repr(entry)); expected format [[id, weight], ...] with numeric weight"))
                end
                (yaml_to_ref(entry[1]), Float64(entry[2]))
            end
        elseif value isa Vector
            return [yaml_to_ref(entry) for entry in value]
        else
            return yaml_to_ref(value)
        end
    end

    # Parse a pair of reference-point refs (z/y axes) from a YAML row field.
    yaml_parse_ref_points = function (row, field)
        !hasfield(typeof(row), field) && return nothing
        val = getfield(row, field)
        val === nothing && return nothing

        if length(val) != 2
            throw(ArgumentError("ref_points must have 2 elements"))
        end

        point_1 = yaml_convert_ref(val[1])
        point_2 = yaml_convert_ref(val[2])
        return (point_1, point_2)
    end

    # Parse a single weighted-ref field (e.g. origin_idx).
    yaml_parse_origin = function (row, field)
        !hasfield(typeof(row), field) && return nothing
        val = getfield(row, field)
        val === nothing && return nothing
        return yaml_convert_ref(val)
    end

    # Load points - SystemStructure handles resolution
    points = Point[]

    if haskey(data, "points")
        point_rows = parse_table(data["points"])
        for (i, row) in enumerate(point_rows)
            # Raw references are passed; SystemStructure resolves them.
            point = call_yaml_constructor(Point, row,
                [:name, :pos_cad, :type],
                [:wing, :transform, :body, :joint, :vel_w, :extra_mass,
                 :body_frame_damping, :world_frame_damping,
                 :area, :drag_coeff, :fix_sphere, :fix_static];
                mappings=Dict(
                    :pos_cad => row -> KVec3(row.pos_cad...),
                    :type => row -> parse_dynamics_type(String(row.type)),
                    :name => row -> yaml_row_name(row, i),
                    # Pass raw references - constructor handles defaults
                    :wing => row -> yaml_ref_field(row, :wing_idx, yaml_to_ref),
                    :transform => row -> yaml_ref_field(row, :transform_idx, yaml_to_ref),
                    # A BODY_STATIC/wing node rides a Body; anchor_b is its body-frame offset.
                    :body => row -> haskey(row, :body_idx) ? yaml_to_ref(row.body_idx) :
                        yaml_ref_field(row, :body, yaml_to_ref),
                    :joint => row -> yaml_ref_field(row, :joint, yaml_to_ref),
                    :vel_w => row -> yaml_vec3(row, :vel_w)
                ))

            point.pos_w .= point.pos_cad
            point.vel_w .= 0.0
            push!(points, point)
        end
    end

    isempty(points) && !haskey(data, "bodies") &&
        error("No points or bodies found in YAML file $(yaml_path).")

    # Build property tables for reference resolution: materials, old-style
    # `elements`, and `segment_properties` (kept for backward compatibility).
    property_tables = Dict{String, Dict{String, NamedTuple}}()
    for key in ("materials", "elements", "segment_properties")
        haskey(data, key) || continue
        property_tables[key] = load_property_table(data, key)
    end

    # Load segments - SystemStructure handles resolution
    segments = Segment[]
    if haskey(data, "segments")
        segment_rows = parse_table(data["segments"])

        for (i, row) in enumerate(segment_rows)
            # Resolve material references and calculate derived properties
            resolved_row = resolve_references(row, property_tables)
            props = Dict{Symbol, Any}(pairs(resolved_row))
            calculate_derived_properties!(props)

            # Convert back to NamedTuple for constructor
            resolved_row = NamedTuple(props)

            # Deprecation check: error if old `type` column present
            if haskey(resolved_row, :type)
                error("Segment YAML `type` column " *
                    "(SegmentType) is removed. Delete " *
                    "the `type` header and column from " *
                    "your YAML segments block.")
            end

            # Create Segment (name, set, point_i, point_j; kwargs)
            segment = call_yaml_constructor(
                Segment, resolved_row,
                [:name, :set, :point_i, :point_j],
                [:l0, :diameter_mm, :unit_stiffness,
                 :unit_damping, :compression_frac, :density];
                mappings=Dict(
                    :set => row -> resolved_set,
                    :point_i => row -> yaml_to_ref(row.point_i),
                    :point_j => row -> yaml_to_ref(row.point_j),
                    :name => row -> yaml_row_name(row, i)
                ))

            push!(segments, segment)
        end
    end

    # Load pulleys - SystemStructure handles resolution
    pulleys = Pulley[]
    if haskey(data, "pulleys")
        pulley_rows = parse_table(data["pulleys"])
        for (i, row) in enumerate(pulley_rows)
            # Create Pulley using new constructor (name, segment_i, segment_j, type)
            pulley = call_yaml_constructor(Pulley, row,
                [:name, :segment_i, :segment_j, :type],
                Vector{Union{}}();
                mappings=Dict(
                    :segment_i => row -> yaml_to_ref(row.segment_i),
                    :segment_j => row -> yaml_to_ref(row.segment_j),
                    :name => row -> yaml_row_name(row, i),
                    :type => row -> parse_dynamics_type(String(row.type))
                ))
            push!(pulleys, pulley)
        end
    end

    # Load twist_surfaces (optional, for deformable wings) - SystemStructure handles resolution
    twist_surfaces = TwistSurface[]
    if !yaml_block_empty(data, "twist_surfaces")
        twist_surface_rows = parse_table(data["twist_surfaces"])

        for (i, row) in enumerate(twist_surface_rows)
            # Create TwistSurface using new constructor (name, points, type, moment_frac)
            twist_surface = call_yaml_constructor(TwistSurface, row,
                [:name, :points, :type, :moment_frac],
                [:damping, :stiffness, :wing, :bodies, :flap_bodies,
                 :flap_axis];
                mappings=Dict(
                    :points => row -> haskey(row, :points) ?
                        [yaml_to_ref(p) for p in row.points] :
                        (haskey(row, :point_idxs) ?
                         [yaml_to_ref(p) for p in row.point_idxs] : NameRef[]),
                    :name => row -> yaml_row_name(row, i),
                    :type => row -> parse_dynamics_type(String(row.type)),
                    :moment_frac => row -> haskey(row, :moment_frac) &&
                        !isnothing(row.moment_frac) ? float(row.moment_frac) : 0.0,
                    :wing => row -> haskey(row, :wing) && !isnothing(row.wing) ?
                        yaml_to_ref(row.wing) : 0,
                    :bodies => row -> haskey(row, :bodies) && !isnothing(row.bodies) ?
                        [yaml_to_ref(b) for b in row.bodies] : NameRef[],
                    :flap_bodies => row -> haskey(row, :flap_bodies) &&
                        !isnothing(row.flap_bodies) ?
                        [yaml_to_ref(b) for b in row.flap_bodies] : NameRef[]
                ))
            push!(twist_surfaces, twist_surface)
        end
    end

    # Load tethers (optional) - SystemStructure handles resolution
    tethers = Tether[]
    if !yaml_block_empty(data, "tethers")
        tether_rows = parse_table(data["tethers"])
        for (i, row) in enumerate(tether_rows)
            tether_name = yaml_row_name(row, i)
            # Detect Route 1 vs Route 2
            has_segments = hasfield(typeof(row),
                :segment_idxs) &&
                !isnothing(row.segment_idxs)
            if has_segments
                # Route 1: explicit segments
                segment_refs = [yaml_to_ref(segment_ref)
                    for segment_ref in row.segment_idxs]
                start_ref = yaml_ref_field(row, :start_point, yaml_to_ref)
                end_ref = yaml_ref_field(row, :end_point, yaml_to_ref)
                stretched_length, tether_force, stretch_frac =
                    parse_tether_init(row, tether_name)
                tether = Tether(tether_name, segment_refs, stretched_length;
                    start_point=start_ref, end_point=end_ref,
                    tether_force, stretch_frac)
            else
                # Route 2: auto-generation
                start_ref = yaml_to_ref(row.start_point)
                end_ref = yaml_to_ref(row.end_point)
                n_segments = Int(row.n_segments)
                stretched_length, tether_force, stretch_frac =
                    parse_tether_init(row, tether_name)
                # Resolve material reference if present
                resolved = resolve_references(
                    row, property_tables)
                props = Dict{Symbol, Any}(
                    pairs(resolved))
                calculate_derived_properties!(props)
                prop_float = key -> begin
                    value = get(props, key, NaN)
                    isnothing(value) ? NaN : Float64(value)
                end
                unit_stiffness = prop_float(:unit_stiffness)
                unit_damping = prop_float(:unit_damping)
                diameter_mm = prop_float(:diameter_mm)
                diameter = isnan(diameter_mm) ? NaN : diameter_mm * 0.001
                density = prop_float(:density)
                tether = Tether(tether_name, stretched_length;
                    start_point=start_ref, end_point=end_ref,
                    n_segments,
                    unit_stiffness, unit_damping,
                    diameter, density, tether_force, stretch_frac)
            end
            push!(tethers, tether)
        end
    end

    # Load winches (optional) - SystemStructure handles resolution
    winches = Winch[]
    if !yaml_block_empty(data, "winches")
        winch_rows = parse_table(data["winches"])
        for (i, row) in enumerate(winch_rows)
            # Create Winch using constructor (name, set, tethers; winch_point)
            winch = call_yaml_constructor(Winch, row,
                [:name, :set, :tethers],
                [:winch_point, :init_vel,
                 :brake, :speed_controlled, :friction_epsilon];
                mappings=Dict(
                    :set => row -> resolved_set,
                    :tethers => row -> [yaml_to_ref(tether_ref)
                        for tether_ref in row.tether_idxs],
                    :winch_point => row -> yaml_to_ref(row.winch_point),
                    :name => row -> yaml_row_name(row, i)
                ))
            push!(winches, winch)
        end
    end

    # Load wings (optional)
    wings = Body[]
    if !yaml_block_empty(data, "wings")
        wing_rows = parse_table(data["wings"])

        for (i, row) in enumerate(wing_rows)
            # Old `type` field is still parsed, with a deprecation warning.
            resolved_wing_type = if !isnothing(dynamics_type)
                dynamics_type
            else
                raw_type_field = if hasfield(typeof(row), :dynamics_type) && !isnothing(row.dynamics_type)
                    String(row.dynamics_type)
                elseif hasfield(typeof(row), :type) && !isnothing(row.type)
                    @warn "Wing YAML field `type` is deprecated; rename to `dynamics_type`."
                    String(row.type)
                else
                    error("Wing entry missing required `dynamics_type` field.")
                end
                parse_wing_type(raw_type_field)
            end

            # Determine aero_mode: kwarg > YAML > default.
            resolved_aero_mode = if !isnothing(aero_mode)
                aero_mode
            elseif hasfield(typeof(row), :aero_mode) &&
                    !isnothing(row.aero_mode)
                parse_aero_mode(String(row.aero_mode))
            else
                resolved_wing_type == RIGID_DYNAMICS ? AeroLinearized() :
                    AeroDirect()
            end

            wing = load_wing(resolved_aero_mode, row, i, data,
                resolved_set, resolved_wing_type, vsm_set, yaml_to_ref,
                yaml_parse_ref_points, yaml_parse_origin, twist_surfaces)
            push!(wings, wing)
        end
    end

    # Load transforms (optional) - SystemStructure handles resolution
    transforms = Transform[]
    if !yaml_block_empty(data, "transforms")
        transform_rows = parse_table(data["transforms"])

        for (i, row) in enumerate(transform_rows)
            # Create Transform using new constructor (name, elevation, azimuth, heading; kwargs)
            transform = call_yaml_constructor(Transform, row,
                [:name, :elevation, :azimuth, :heading],
                [:base_point, :base_pos, :base_transform,
                 :wing, :rot_point,
                 :elevation_vel, :azimuth_vel, :turn_rate];
                mappings=Dict(
                    :elevation => row -> deg2rad(row.elevation),
                    :azimuth => row -> deg2rad(row.azimuth),
                    :heading => row -> deg2rad(row.heading),
                    :elevation_vel => row -> hasfield(typeof(row), :elevation_vel) && !isnothing(row.elevation_vel) ?
                        deg2rad(row.elevation_vel) : 0.0,
                    :azimuth_vel => row -> hasfield(typeof(row), :azimuth_vel) && !isnothing(row.azimuth_vel) ?
                        deg2rad(row.azimuth_vel) : 0.0,
                    :turn_rate => row -> hasfield(typeof(row), :turn_rate) && !isnothing(row.turn_rate) ?
                        deg2rad(row.turn_rate) : 0.0,
                    :base_pos => row -> KVec3(row.base_pos...),
                    :base_point => row -> yaml_to_ref(row.base_point_idx),
                    :base_transform => row ->
                        yaml_ref_field(row, :base_transform_idx, yaml_to_ref),
                    :rot_point => row ->
                        yaml_ref_field(row, :rot_point_idx, yaml_to_ref),
                    :wing => row -> yaml_ref_field(row, :wing_idx, yaml_to_ref),
                    :name => row -> yaml_row_name(row, i)
                ))
            push!(transforms, transform)
            elev_deg = rad2deg(transform.elevation)
            azim_deg = rad2deg(transform.azimuth)
            head_deg = rad2deg(transform.heading)
            elev_vel_deg = rad2deg(transform.elevation_vel)
            azim_vel_deg = rad2deg(transform.azimuth_vel)
            turn_rate_deg = rad2deg(transform.turn_rate)
            info_msg = "  ✓ Transform $(i) created: " *
                       "elevation=$(elev_deg)°, azimuth=$(azim_deg)°, heading=$(head_deg)°"
            if elev_vel_deg != 0.0 || azim_vel_deg != 0.0
                info_msg *= ", elevation_vel=$(elev_vel_deg)°/s, azimuth_vel=$(azim_vel_deg)°/s"
            end
            if turn_rate_deg != 0.0
                info_msg *= ", turn_rate=$(turn_rate_deg)°/s"
            end
            @info info_msg
        end
    end

    # Plain rigid bodies + beam/elastic joints (empty when the blocks are absent).
    bodies = load_yaml_bodies(data, yaml_to_ref)
    timoshenko_joints = load_yaml_joints(TimoshenkoJoint, data,
        "timoshenko_joints", (:EA, :GA, :GJ, :EIy, :EIz),
        (:shear_coeff, :damping_trans, :damping_rot, :rest_length, :radius);
        yaml_to_ref)
    elastic_joints = load_yaml_joints(ElasticJoint, data, "elastic_joints",
        (:stiffness_axial, :stiffness_shear, :stiffness_torsion,
         :stiffness_bending), (:damping_trans, :damping_rot, :radius);
        yaml_to_ref)

    # The SystemStructure constructor handles WING->STATIC when no wings exist.
    return SystemStructure(system_name, resolved_set; points, twist_surfaces,
        segments, pulleys, tethers, winches, wings,
        transforms, bodies, elastic_joints, timoshenko_joints,
        ignore_l0, vsm_set)
end
