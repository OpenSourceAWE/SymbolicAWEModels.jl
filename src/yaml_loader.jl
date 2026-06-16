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
                        # Skip 'name' — it identifies the referenced
                        # item, not a property to inherit.
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
    # Calculate unit_stiffness from material properties if missing or if it's a string (material name)
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

    # Check format: if first row is a Dict,
    # use dict format; if Array, use header format
    first_row = first(rows)

    if first_row isa AbstractDict
        # Dict format: each row is already a dict with named keys
        # Convert each dict to a NamedTuple
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
    text_upper == "WING" && return WING
    text_upper == "QUASI_STATIC" && return QUASI_STATIC
    text_upper == "FIXED" && return FIXED
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

function parse_aero_mode(text::String)
    key = lowercase(replace(text, "_" => ""))
    key in ("aeronone", "none") && return AeroNone()
    key in ("aerodirect", "direct") && return AeroDirect()
    key in ("aerolinearized", "linearized") && return AeroLinearized()
    key in ("aeroplate", "plate") && return AeroPlate()
    error("Unknown aero model: $text")
end

"""
    load_wing(mode::AbstractAeroModel, row, idx, data, set, wing_type, vsm_set,
              yaml_to_ref, yaml_parse_ref_points, yaml_parse_origin, twist_surfaces)

Build a wing from a parsed YAML `row`, dispatched on its aero `mode`. The default
(VSM-backed modes) builds a [`VSMWing`](@ref); [`AeroPlate`](@ref) builds a
flat-plate wing via [`load_plate_wing`](@ref). Add a method to load a wing for a
custom aero mode.

Unused parameters (`data`, `twist_surfaces`) are accepted with `_` for dispatch
compatibility with the [`AeroPlate`](@ref) method which needs them.

"""
function load_wing(mode::AbstractAeroModel, row, idx, _data, set, wing_type,
                   vsm_set, yaml_to_ref, yaml_parse_ref_points,
                   yaml_parse_origin, _twist_surfaces)
    if wing_type == PARTICLE_DYNAMICS
        # PARTICLE_DYNAMICS wings need z_ref_points, y_ref_points, origin
        # Pass raw values - constructor handles defaults
        return call_yaml_constructor(VSMWing, row,
            [:name, :set, :twist_surfaces, :vsm_set],
            [:transform, :y_damping, :angular_damping,
             :dynamics_type, :aero,
             :z_ref_points, :y_ref_points, :origin, :pos_cad,
             :aero_scale_chord];
            mappings=Dict(
                :set => _row -> set,
                :twist_surfaces => row -> hasfield(typeof(row), :twist_surfaces) &&
                    !isnothing(row.twist_surfaces) ?
                    [yaml_to_ref(twist_surface_ref) for twist_surface_ref in row.twist_surfaces] : [],
                :vsm_set => _row -> vsm_set,
                :dynamics_type => _row -> wing_type,
                :aero => _row -> mode,
                :name => row -> begin
                    if haskey(row, :name) && !isnothing(row.name)
                        Symbol(row.name)
                    else
                        idx  # Use index as name if no name provided
                    end
                end,
                :transform => row -> begin
                    if hasfield(typeof(row), :transform_idx) && !isnothing(row.transform_idx)
                        yaml_to_ref(row.transform_idx)
                    else
                        nothing  # Constructor handles default
                    end
                end,
                :z_ref_points => row ->
                    yaml_parse_ref_points(row, :z_ref_points),
                :y_ref_points => row ->
                    yaml_parse_ref_points(row, :y_ref_points),
                :origin => row ->
                    yaml_parse_origin(row, :origin_idx),
                :aero_scale_chord => row ->
                    hasfield(typeof(row), :aero_scale_chord) && !isnothing(row.aero_scale_chord) ?
                        float(row.aero_scale_chord) : 0.0,
                :pos_cad => _row -> begin
                    # Note: pos_cad will be set from origin point position after resolution
                    # For now, return nothing - SystemStructure will handle this
                    nothing
                end
            ))
    else  # RIGID_DYNAMICS
        # Pass raw values - constructor handles defaults
        return call_yaml_constructor(VSMWing, row,
            [:name, :set, :twist_surfaces, :vsm_set],
            [:transform, :y_damping, :angular_damping,
             :dynamics_type, :aero, :aero_scale_chord,
             :aero_z_offset, :pos_cad,
             :z_ref_points, :y_ref_points, :origin];
            mappings=Dict(
                :set => _row -> set,
                :aero => _row -> mode,
                :twist_surfaces => row -> hasfield(typeof(row), :twist_surfaces) &&
                    !isnothing(row.twist_surfaces) ?
                    [yaml_to_ref(twist_surface_ref) for twist_surface_ref in row.twist_surfaces] : [],
                :vsm_set => _row -> vsm_set,
                :dynamics_type => _row -> wing_type,
                :name => row -> begin
                    if haskey(row, :name) && !isnothing(row.name)
                        Symbol(row.name)
                    else
                        idx
                    end
                end,
                :transform => row -> begin
                    if hasfield(typeof(row), :transform_idx) &&
                       !isnothing(row.transform_idx)
                        yaml_to_ref(row.transform_idx)
                    else
                        nothing
                    end
                end,
                :pos_cad => row -> begin
                    if !hasfield(typeof(row), :pos_cad) ||
                       row.pos_cad === nothing
                        return nothing
                    end
                    KVec3(row.pos_cad...)
                end,
                :aero_scale_chord => row ->
                    hasfield(typeof(row), :aero_scale_chord) &&
                    !isnothing(row.aero_scale_chord) ?
                        float(row.aero_scale_chord) : 0.0,
                :z_ref_points => row ->
                    yaml_parse_ref_points(row, :z_ref_points),
                :y_ref_points => row ->
                    yaml_parse_ref_points(row, :y_ref_points),
                :origin => row ->
                    yaml_parse_origin(row, :origin_idx)
            ))
    end
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

    # Convert one YAML ref spec to WeightedRefPoints input form.
    # Supports:
    #   :name         → Symbol
    #   7             → Int
    #   [a, b]        → equal-weight average
    #   [[a, w], ...] → explicit weights as (name, weight) tuples
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

    # Parse [a, b] or [a, [[id, w], ...]] style reference-point
    # fields from YAML rows (pair of refs for z/y axes).
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
    # Accepts the same shapes as one side of yaml_parse_ref_points.
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
            # Create Point using new constructor (name as first positional arg)
            # Raw references are passed - SystemStructure will resolve them
            point = call_yaml_constructor(Point, row,
                [:name, :pos_cad, :type],
                [:wing, :transform, :extra_mass,
                 :body_frame_damping, :world_frame_damping,
                 :area, :drag_coeff];
                mappings=Dict(
                    :pos_cad => row -> KVec3(row.pos_cad...),
                    :type => row -> parse_dynamics_type(String(row.type)),
                    :name => row -> haskey(row, :name) && !isnothing(row.name) ? Symbol(row.name) : i,
                    # Pass raw references - constructor handles defaults
                    :wing => row -> haskey(row, :wing_idx) ? yaml_to_ref(row.wing_idx) : nothing,
                    :transform => row -> haskey(row, :transform_idx) ? yaml_to_ref(row.transform_idx) : nothing,
                    :body_frame_damping => row -> haskey(row, :body_frame_damping) ? row.body_frame_damping : nothing,
                    :world_frame_damping => row -> haskey(row, :world_frame_damping) ? row.world_frame_damping : nothing
                ))

            point.pos_w .= point.pos_cad
            point.vel_w .= 0.0
            push!(points, point)
        end
    end

    isempty(points) &&
        error("No points found in YAML file $(yaml_path).")

    # Build property tables for reference resolution
    property_tables = Dict{String, Dict{String, NamedTuple}}()

    # Load materials table
    if haskey(data, "materials")
        materials_dict = Dict{String, NamedTuple}()
        material_rows = parse_table(data["materials"])
        for row in material_rows
            name = String(row.name)
            materials_dict[name] = row
        end
        property_tables["materials"] = materials_dict
    end

    # Load elements table (old-style segment properties)
    if haskey(data, "elements")
        elements_dict = Dict{String, NamedTuple}()
        element_rows = parse_table(data["elements"])
        for row in element_rows
            name = String(row.name)
            elements_dict[name] = row
        end
        property_tables["elements"] = elements_dict
    end

    # Load segment_properties table (for backward compatibility)
    if haskey(data, "segment_properties")
        segment_props_dict = Dict{String, NamedTuple}()
        prop_rows = parse_table(data["segment_properties"])
        for row in prop_rows
            name = String(row.name)
            segment_props_dict[name] = row
        end
        property_tables["segment_properties"] = segment_props_dict
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
                 :unit_damping, :compression_frac];
                mappings=Dict(
                    :set => _row -> resolved_set,
                    :point_i => row -> yaml_to_ref(row.point_i),
                    :point_j => row -> yaml_to_ref(row.point_j),
                    :name => row -> begin
                        if haskey(row, :name) &&
                                !isnothing(row.name)
                            Symbol(row.name)
                        else
                            i
                        end
                    end
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
                    :name => row -> begin
                        if haskey(row, :name) && !isnothing(row.name)
                            Symbol(row.name)
                        else
                            i  # Use index as name if no name provided
                        end
                    end,
                    :type => row -> parse_dynamics_type(String(row.type))
                ))
            push!(pulleys, pulley)
        end
    end

    # Load twist_surfaces (optional, for deformable wings) - SystemStructure handles resolution
    twist_surfaces = TwistSurface[]
    if haskey(data, "twist_surfaces") &&
       haskey(data["twist_surfaces"], "data") &&
       data["twist_surfaces"]["data"] !== nothing &&
       !isempty(data["twist_surfaces"]["data"])
        twist_surface_rows = parse_table(data["twist_surfaces"])

        for (i, row) in enumerate(twist_surface_rows)
            # Create TwistSurface using new constructor (name, points, type, moment_frac)
            twist_surface = call_yaml_constructor(TwistSurface, row,
                [:name, :points, :type, :moment_frac],
                [:damping];
                mappings=Dict(
                    :points => row -> [yaml_to_ref(point) for point in row.point_idxs],
                    :name => row -> begin
                        if haskey(row, :name) && !isnothing(row.name)
                            Symbol(row.name)
                        else
                            i  # Use index as name if no name provided
                        end
                    end,
                    :type => row -> parse_dynamics_type(
                        String(row.type))
                ))
            push!(twist_surfaces, twist_surface)
        end
    end

    # Load tethers (optional) - SystemStructure handles resolution
    tethers = Tether[]
    if haskey(data, "tethers") &&
       haskey(data["tethers"], "data") &&
       data["tethers"]["data"] !== nothing &&
       !isempty(data["tethers"]["data"])
        tether_rows = parse_table(data["tethers"])
        for (i, row) in enumerate(tether_rows)
            tether_name = if haskey(row, :name) &&
                             !isnothing(row.name)
                Symbol(row.name)
            else
                i
            end
            # Detect Route 1 vs Route 2
            has_segments = hasfield(typeof(row),
                :segment_idxs) &&
                !isnothing(row.segment_idxs)
            if has_segments
                # Route 1: explicit segments
                segment_refs = [yaml_to_ref(segment_ref)
                    for segment_ref in row.segment_idxs]
                start_ref = if hasfield(typeof(row),
                        :start_point) &&
                        !isnothing(row.start_point)
                    yaml_to_ref(row.start_point)
                else
                    nothing
                end
                end_ref = if hasfield(typeof(row),
                        :end_point) &&
                        !isnothing(row.end_point)
                    yaml_to_ref(row.end_point)
                else
                    nothing
                end
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
                unit_stiffness = get(props, :unit_stiffness,
                    NaN)
                unit_stiffness = isnothing(unit_stiffness) ? NaN :
                    Float64(unit_stiffness)
                unit_damping = get(props, :unit_damping,
                    NaN)
                unit_damping = isnothing(unit_damping) ? NaN :
                    Float64(unit_damping)
                diameter_mm = get(props, :diameter_mm,
                    NaN)
                diameter_mm = isnothing(diameter_mm) ? NaN :
                    Float64(diameter_mm)
                diameter = isnan(diameter_mm) ? NaN :
                    diameter_mm * 0.001
                tether = Tether(tether_name, stretched_length;
                    start_point=start_ref, end_point=end_ref,
                    n_segments,
                    unit_stiffness, unit_damping,
                    diameter, tether_force, stretch_frac)
            end
            push!(tethers, tether)
        end
    end

    # Load winches (optional) - SystemStructure handles resolution
    winches = Winch[]
    if haskey(data, "winches") &&
       haskey(data["winches"], "data") &&
       data["winches"]["data"] !== nothing &&
       !isempty(data["winches"]["data"])
        winch_rows = parse_table(data["winches"])
        for (i, row) in enumerate(winch_rows)
            # Create Winch using constructor (name, set, tethers; winch_point)
            winch = call_yaml_constructor(Winch, row,
                [:name, :set, :tethers],
                [:winch_point, :init_vel,
                 :brake, :speed_controlled, :friction_epsilon];
                mappings=Dict(
                    :set => _row -> resolved_set,
                    :tethers => row -> [yaml_to_ref(tether_ref)
                        for tether_ref in row.tether_idxs],
                    :winch_point => row -> begin
                        yaml_to_ref(row.winch_point)
                    end,
                    :name => row -> begin
                        if haskey(row, :name) &&
                           !isnothing(row.name)
                            Symbol(row.name)
                        else
                            i
                        end
                    end
                ))
            push!(winches, winch)
        end
    end

    # Load wings (optional)
    wings = AbstractWing[]
    if haskey(data, "wings") &&
       haskey(data["wings"], "data") &&
       data["wings"]["data"] !== nothing &&
       !isempty(data["wings"]["data"])
        wing_rows = parse_table(data["wings"])

        for (i, row) in enumerate(wing_rows)
            # Use provided dynamics_type parameter or parse from YAML
            # Support old `type` field with deprecation warning
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

            # Build kwargs based on wing type - SystemStructure handles resolution
            # Determine aero_mode: kwarg > YAML > default
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
    if haskey(data, "transforms") &&
       haskey(data["transforms"], "data") &&
       data["transforms"]["data"] !== nothing &&
       !isempty(data["transforms"]["data"])
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
                    :base_transform => row -> begin
                        if hasfield(typeof(row), :base_transform_idx) && !isnothing(row.base_transform_idx)
                            yaml_to_ref(row.base_transform_idx)
                        else
                            nothing
                        end
                    end,
                    :rot_point => row -> begin
                        if hasfield(typeof(row), :rot_point_idx) && !isnothing(row.rot_point_idx)
                            yaml_to_ref(row.rot_point_idx)
                        else
                            nothing
                        end
                    end,
                    :wing => row -> begin
                        if !hasfield(typeof(row), :wing_idx) || row.wing_idx === nothing
                            return nothing
                        end
                        yaml_to_ref(row.wing_idx)
                    end,
                    :name => row -> begin
                        if haskey(row, :name) && !isnothing(row.name)
                            Symbol(row.name)
                        else
                            i  # Use index as name if no name provided
                        end
                    end
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

    # SystemStructure constructor now handles WING→STATIC
    # conversion when no wings are defined
    return SystemStructure(system_name, resolved_set; points, twist_surfaces,
        segments, pulleys, tethers, winches, wings,
        transforms, ignore_l0, vsm_set)
end

# ==================== EXPORT (YAML) ==================== #

"""
    dynamics_type_to_str(dt::DynamicsType) -> String

Convert a `DynamicsType` enum to its YAML string representation.
"""
function dynamics_type_to_str(dt::DynamicsType)
    dt == DYNAMIC && return "DYNAMIC"
    dt == QUASI_STATIC && return "QUASI_STATIC"
    dt == WING && return "WING"
    dt == STATIC && return "STATIC"
    dt == FIXED && return "FIXED"
    return "DYNAMIC"
end

"""
    wing_type_to_str(wt::WingType) -> String

Convert a `WingType` enum to its YAML string representation.
"""
function wing_type_to_str(wt::WingType)
    wt == RIGID_DYNAMICS && return "RIGID_DYNAMICS"
    wt == PARTICLE_DYNAMICS && return "PARTICLE_DYNAMICS"
    return "RIGID_DYNAMICS"
end

"""
    aero_mode_to_str(aero::AbstractAeroModel) -> String

Convert an aero mode to YAML string representation.
"""
function aero_mode_to_str(aero::AbstractAeroModel)
    aero isa AeroNone && return "AeroNone"
    aero isa AeroDirect && return "AeroDirect"
    aero isa AeroLinearized && return "AeroLinearized"
    aero isa AeroPlate && return "AeroPlate"
    error("Unknown aero mode: $(typeof(aero)). " *
          "Add an overload of aero_mode_to_str for this type.")
end

"""
    name_or_nothing(name) -> Any

Return the name as-is (Symbol) if non-nothing, otherwise `nothing`.
Useful for fields that may have `nothing` name.
"""
name_or_nothing(name) = name === nothing ? nothing : name

"""
    ref_to_yaml(ref) -> Any

Convert a `NameRef` back to the appropriate YAML representation.
Symbols become strings, integers stay as integers.
"""
ref_to_yaml(ref::Symbol) = String(ref)
ref_to_yaml(ref::Int) = ref
ref_to_yaml(::Nothing) = nothing

"""
    ref_point_to_yaml(wrp::WeightedRefPoints) -> Any

Convert a `WeightedRefPoints` to YAML-friendly representation:
- Single point → name string or integer
- Equal-weight average → array of refs
- Weighted → array of [ref, weight] pairs
"""
function ref_point_to_yaml(wrp::WeightedRefPoints)
    if length(wrp.refs) == 1 && length(wrp.ids) == 0
        return ref_to_yaml(wrp.refs[1])
    elseif length(wrp.ids) == 1 && isempty(wrp.refs)
        return wrp.ids[1]
    end
    # Multi-point case
    has_explicit_refs = !isempty(wrp.refs)
    _has_ids = !isempty(wrp.ids)
    n = length(has_explicit_refs ? wrp.refs : wrp.ids)
    if all(w -> abs(w - 1.0/n) < 1e-10, wrp.weights)
        # Equal weights → simple list
        if has_explicit_refs
            return [ref_to_yaml(r) for r in wrp.refs]
        else
            return collect(wrp.ids)
        end
    else
        # Weighted → list of [ref, weight] pairs
        result = []
        for i in eachindex(wrp.weights)
            ref = has_explicit_refs ? ref_to_yaml(wrp.refs[i]) : wrp.ids[i]
            push!(result, [ref, wrp.weights[i]])
        end
        return result
    end
end

"""
    ref_point_tuple_to_yaml(tup) -> Vector

Convert a tuple of two `WeightedRefPoints` to YAML-friendly `[a, b]` format.
"""
function ref_point_tuple_to_yaml(tup)
    return [ref_point_to_yaml(tup[1]), ref_point_to_yaml(tup[2])]
end

"""
    vec3_to_yaml(v) -> Vector{Float64}

Convert a 3-element vector to a plain Float64 array for YAML output.
"""
vec3_to_yaml(v) = [Float64(v[1]), Float64(v[2]), Float64(v[3])]

"""
    to_yaml_flow(value) -> String

Convert a value to a YAML flow-style representation string.

Always returns valid YAML, unlike `repr()` which can emit
Julia type annotations like `Any[...]` for vectors with
non-concrete element types.

# Examples
```julia
to_yaml_flow("foo")            # -> "\"foo\""
to_yaml_flow(42)               # -> "42"
to_yaml_flow(1.5)              # -> "1.5"
to_yaml_flow(nothing)          # -> "nothing"
to_yaml_flow(["a", "b"])      # -> "[\\"a\\", \\"b\\"]"
to_yaml_flow([["a", 0.5]])    # -> "[[\"a\", 0.5]]"
```
"""
function to_yaml_flow(value::AbstractString)
    return repr(value)
end

function to_yaml_flow(value::Symbol)
    return repr(String(value))
end

function to_yaml_flow(value::Integer)
    return string(value)
end

function to_yaml_flow(value::AbstractFloat)
    return string(value)
end

function to_yaml_flow(::Nothing)
    return "nothing"
end

function to_yaml_flow(value::AbstractVector)
    parts = join((to_yaml_flow(v) for v in value), ", ")
    return "[$parts]"
end

function to_yaml_flow(value::Tuple)
    return to_yaml_flow(collect(value))
end

"""
    save_sys_struct_to_yaml(sys::SystemStructure, yaml_path::AbstractString;
                           materials::Bool=true)

Export a `SystemStructure` to a YAML file in the same component-based format
that `load_sys_struct_from_yaml` can read back.

# Arguments
- `sys::SystemStructure`: The system structure to export.
- `yaml_path::AbstractString`: Output YAML file path.

# Keyword Arguments
- `materials::Bool=true`: Whether to include a materials section.
"""
function save_sys_struct_to_yaml(sys::SystemStructure, yaml_path::AbstractString;
                                 materials::Bool=true)
    open(yaml_path, "w") do io
        # Write header comment
        write(io, "##############################\n")
        write(io, "## System Structure ##########\n")
        write(io, "##############################\n")
        write(io, "# This file was generated by save_sys_struct_to_yaml\n\n")

        # ==================== Materials ==================== #
        if materials
            write(io, "materials:\n")
            write(io, "  headers: [name, youngs_modulus, density, damping_per_stiffness]\n")
            write(io, "  data:\n")
            write(io, "    - [dyneema, 55000000000.0, 724, 0.00077]\n")
            write(io, "\n")
        end

        # ==================== Points ==================== #
        write(io, "points:\n")
        write(io, "  headers: [name, pos_cad, type, wing_idx, transform_idx, extra_mass, body_frame_damping, world_frame_damping, area, drag_coeff]\n")
        write(io, "  data:\n")
        for point in sys.points
            point.name === nothing && continue
            name_str = point.name isa Symbol ? String(point.name) : string(point.name)
            pos = vec3_to_yaml(point.pos_cad)
            dtype = dynamics_type_to_str(point.type)
            widx = wing_type_to_yaml_ref(point.wing_idx, sys.wings)
            tidx = transform_to_yaml_ref(point.transform_idx, sys.transforms)
            bfd = vec3_to_yaml(point.body_frame_damping)
            wfd = vec3_to_yaml(point.world_frame_damping)
            write(io, "    - [$(to_yaml_flow(name_str)), $pos, $(to_yaml_flow(dtype)), $(to_yaml_flow(widx)), $(to_yaml_flow(tidx)), $(Float64(point.extra_mass)), $bfd, $wfd, $(Float64(point.area)), $(Float64(point.drag_coeff))]\n")
        end
        write(io, "\n")

        # ==================== Segments ==================== #
        write(io, "segments:\n")
        write(io, "  headers: [name, point_i, point_j, l0, diameter_mm, unit_stiffness, unit_damping, compression_frac]\n")
        write(io, "  data:\n")
        for seg in sys.segments
            seg.name === nothing && continue
            name_str = seg.name isa Symbol ? String(seg.name) : string(seg.name)
            p1 = idx_to_name_or_number(seg.point_idxs[1], sys.points)
            p2 = idx_to_name_or_number(seg.point_idxs[2], sys.points)
            l0_str = seg.l0 == 0.0 ? "nothing" : string(Float64(seg.l0))
            write(io, "    - [$(to_yaml_flow(name_str)), $(to_yaml_flow(p1)), $(to_yaml_flow(p2)), $l0_str, $(Float64(seg.diameter * 1000.0)), $(Float64(seg.unit_stiffness)), $(Float64(seg.unit_damping)), $(Float64(seg.compression_frac))]\n")
        end
        write(io, "\n")

        # ==================== Pulleys ==================== #
        if !isempty(sys.pulleys)
            write(io, "pulleys:\n")
            write(io, "  headers: [name, segment_i, segment_j, type]\n")
            write(io, "  data:\n")
            for pulley in sys.pulleys
                pulley.name === nothing && continue
                name_str = pulley.name isa Symbol ? String(pulley.name) : string(pulley.name)
                s1 = idx_to_name_or_number(pulley.segment_idxs[1], sys.segments)
                s2 = idx_to_name_or_number(pulley.segment_idxs[2], sys.segments)
                write(io, "    - [$(to_yaml_flow(name_str)), $(to_yaml_flow(s1)), $(to_yaml_flow(s2)), $(to_yaml_flow(dynamics_type_to_str(pulley.type)))]\n")
            end
            write(io, "\n")
        end

        # ==================== Twist Surfaces ==================== #
        if !isempty(sys.twist_surfaces)
            write(io, "twist_surfaces:\n")
            write(io, "  headers: [name, point_idxs, type, moment_frac, damping]\n")
            write(io, "  data:\n")
            for ts in sys.twist_surfaces
                ts.name === nothing && continue
                name_str = ts.name isa Symbol ? String(ts.name) : string(ts.name)
                pts = "[$(join([to_yaml_flow(idx_to_name_or_number(idx, sys.points)) for idx in ts.point_idxs], ", "))]"
                write(io, "    - [$(to_yaml_flow(name_str)), $pts, $(to_yaml_flow(dynamics_type_to_str(ts.type))), $(Float64(ts.moment_frac)), $(Float64(ts.damping))]\n")
            end
            write(io, "\n")
        end

        # ==================== Tethers ==================== #
        if !isempty(sys.tethers)
            write(io, "tethers:\n")
            write(io, "  headers: [name, start_point, end_point, n_segments, material, diameter_mm, init_stretched_length]\n")
            write(io, "  data:\n")
            for tether in sys.tethers
                tether.name === nothing && continue
                name_str = tether.name isa Symbol ? String(tether.name) : string(tether.name)
                sp = idx_to_name_or_number(tether.start_point_idx, sys.points)
                ep = idx_to_name_or_number(tether.end_point_idx, sys.points)
                nseg = tether.n_segments > 0 ? tether.n_segments : length(tether.segment_idxs)
                stretched = something(tether.init_stretched_len, tether.stretched_len)
                stretched_str = stretched isa Number ? string(Float64(stretched)) : string(Float64(tether.stretched_len))
                write(io, "    - [$(to_yaml_flow(name_str)), $(to_yaml_flow(sp)), $(to_yaml_flow(ep)), $nseg, $(to_yaml_flow("dyneema")), $(Float64(tether.diameter * 1000.0)), $stretched_str]\n")
            end
            write(io, "\n")
        end

        # ==================== Winches ==================== #
        if !isempty(sys.winches)
            write(io, "winches:\n")
            write(io, "  headers: [name, tether_idxs, winch_point]\n")
            write(io, "  data:\n")
            for winch in sys.winches
                winch.name === nothing && continue
                name_str = winch.name isa Symbol ? String(winch.name) : string(winch.name)
                tether_refs = "[$(join([to_yaml_flow(idx_to_name_or_number(idx, sys.tethers)) for idx in winch.tether_idxs], ", "))]"
                wp = idx_to_name_or_number(winch.winch_point_idx, sys.points)
                write(io, "    - [$(to_yaml_flow(name_str)), $tether_refs, $(to_yaml_flow(wp))]\n")
            end
            write(io, "\n")
        end

        # ==================== Wings ==================== #
        if !isempty(sys.wings)
            write(io, "wings:\n")
            write(io, "  data:\n")
            for wing in sys.wings
                wing.name === nothing && continue
                write(io, "    - name: $(string(wing.name))\n")
                write(io, "      dynamics_type: $(wing_type_to_str(wing.dynamics_type))\n")
                write(io, "      aero_mode: $(aero_mode_to_str(wing.aero))\n")

                # Point references
                wing_point_names = String[]
                for point in sys.points
                    if point.wing_idx == wing.idx && point.name !== nothing
                        push!(wing_point_names, string(point.name))
                    end
                end
                write(io, "      point_idxs: $(to_yaml_flow(wing_point_names))\n")

                # Twist surface references
                ts_names = [string(sys.twist_surfaces[idx].name) for idx in wing.twist_surface_idxs
                           if sys.twist_surfaces[idx].name !== nothing]
                if !isempty(ts_names)
                    write(io, "      twist_surfaces: $(to_yaml_flow(ts_names))\n")
                end

                # Origin and ref points
                if wing.origin !== nothing
                    write(io, "      origin_idx: $(to_yaml_flow(ref_point_to_yaml(wing.origin)))\n")
                end
                if wing.z_ref_points !== nothing
                    write(io, "      z_ref_points: $(to_yaml_flow(ref_point_tuple_to_yaml(wing.z_ref_points)))\n")
                end
                if wing.y_ref_points !== nothing
                    write(io, "      y_ref_points: $(to_yaml_flow(ref_point_tuple_to_yaml(wing.y_ref_points)))\n")
                end

                # Transform reference
                if wing.transform_idx > 0 && wing.transform_idx <= length(sys.transforms)
                    tf = sys.transforms[wing.transform_idx]
                    if tf.name !== nothing
                        write(io, "      transform_idx: $(string(tf.name))\n")
                    end
                end

                # Aero offset
                if hasproperty(wing.aero, :aero_z_offset)
                    write(io, "      aero_z_offset: $(Float64(wing.aero.aero_z_offset))\n")
                end
            end
            write(io, "\n")
        end

        # ==================== Transforms ==================== #
        if !isempty(sys.transforms)
            write(io, "transforms:\n")
            write(io, "  data:\n")
            for tf in sys.transforms
                tf.name === nothing && continue
                write(io, "    - name: $(string(tf.name))\n")
                write(io, "      elevation: $(Float64(rad2deg(tf.elevation)))\n")
                write(io, "      azimuth: $(Float64(rad2deg(tf.azimuth)))\n")
                write(io, "      heading: $(Float64(rad2deg(tf.heading)))\n")

                if tf.base_point_idx !== nothing && tf.base_point_idx > 0
                    bp = sys.points[tf.base_point_idx]
                    if bp.name !== nothing
                        write(io, "      base_point_idx: $(string(bp.name))\n")
                    end
                end

                if tf.base_pos !== nothing
                    write(io, "      base_pos: $(to_yaml_flow(vec3_to_yaml(tf.base_pos)))\n")
                end

                if tf.wing_idx !== nothing && tf.wing_idx > 0
                    w = sys.wings[tf.wing_idx]
                    if w.name !== nothing
                        write(io, "      wing_idx: $(string(w.name))\n")
                    end
                end

                if tf.rot_point_idx !== nothing && tf.rot_point_idx > 0
                    rp = sys.points[tf.rot_point_idx]
                    if rp.name !== nothing
                        write(io, "      rot_point_idx: $(string(rp.name))\n")
                    end
                end

                if tf.base_transform_idx !== nothing && tf.base_transform_idx > 0
                    bt = sys.transforms[tf.base_transform_idx]
                    if bt.name !== nothing
                        write(io, "      base_transform_idx: $(string(bt.name))\n")
                    end
                end

                if tf.elevation_vel != 0.0
                    write(io, "      elevation_vel: $(Float64(rad2deg(tf.elevation_vel)))\n")
                end
                if tf.azimuth_vel != 0.0
                    write(io, "      azimuth_vel: $(Float64(rad2deg(tf.azimuth_vel)))\n")
                end
                if tf.turn_rate != 0.0
                    write(io, "      turn_rate: $(Float64(rad2deg(tf.turn_rate)))\n")
                end
            end
        end
    end

    @info "System structure saved to $yaml_path"
    return nothing
end

# ==================== Lookup helpers for export ==================== #

"""
    idx_to_name_or_number(idx, collection) -> Any

Convert a resolved index back to its symbolic name (if available) or keep
the numeric index. Used when exporting YAML to produce human-readable references.
"""
function idx_to_name_or_number(idx, collection)
    if idx < 1 || idx > length(collection)
        return idx
    end
    item = collection[idx]
    if item.name !== nothing
        name = item.name
        return name isa Symbol ? String(name) : name
    end
    return idx
end

"""
    wing_type_to_yaml_ref(wing_idx, wings) -> Any

Convert a wing index to its name reference for YAML export.
Returns the wing name as a string, or the integer index if nameless.
"""
function wing_type_to_yaml_ref(wing_idx, wings)
    if wing_idx < 1 || wing_idx > length(wings)
        return 0
    end
    wing = wings[wing_idx]
    if wing.name !== nothing
        name = wing.name
        return name isa Symbol ? String(name) : name
    end
    return wing_idx
end

"""
    transform_to_yaml_ref(tf_idx, transforms) -> Any

Convert a transform index to its name reference for YAML export.
"""
function transform_to_yaml_ref(tf_idx, transforms)
    if tf_idx < 1 || tf_idx > length(transforms)
        return 0
    end
    tf = transforms[tf_idx]
    if tf.name !== nothing
        name = tf.name
        return name isa Symbol ? String(name) : name
    end
    return tf_idx
end
