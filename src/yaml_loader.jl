# Copyright (c) 2025 Bart van de Lint, Jelle Poland
# SPDX-License-Identifier: MPL-2.0

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
    vec = Vector{Int}(value)
    return (T(vec[1]), T(vec[2]))
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
                    for (k, v) in pairs(ref_props)
                        # Skip 'name' — it identifies the referenced
                        # item, not a property to inherit.
                        k === :name && continue
                        if !haskey(resolved, k) || resolved[k] === nothing
                            resolved[k] = v
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
            d_m = Float64(props[:diameter_mm]) / 1000.0  # mm to m
            A = π * (d_m / 2)^2
            E = Float64(props[:youngs_modulus])
            props[:unit_stiffness] = E * A
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

function parse_table(tbl)::Vector{NamedTuple}
    haskey(tbl, "data") || throw(ArgumentError("table is missing `data`"))

    rows = tbl["data"]
    (isnothing(rows) || isempty(rows)) && return NamedTuple[]

    # Check format: if first row is a Dict,
    # use dict format; if Array, use header format
    first_row = first(rows)

    if first_row isa AbstractDict
        # Dict format: each row is already a dict with named keys
        # Convert each dict to a NamedTuple
        out = NamedTuple[]
        for row in rows
            nt = NamedTuple{Tuple(Symbol.(keys(row)))}(
                Tuple(values(row)))
            push!(out, nt)
        end
        return out
    else
        # Array format: requires headers
        haskey(tbl, "headers") ||
            throw(ArgumentError(
                "table with array rows requires `headers`"))
        headers = String.(tbl["headers"])

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
            nt = NamedTuple{Tuple(Symbol.(headers))}(Tuple(row))
            push!(out, nt)
        end
        return out
    end
end

"""
    _extract_args(row, args_spec, mappings)

Extract positional constructor arguments from a YAML row.

For each name in `args_spec`, this helper first checks for a
mapping in `mappings`, then falls back to `row[arg_name]`.
Throws an error if a required argument is missing.
"""
function _extract_args(row, args_spec, mappings)
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
    args = _extract_args(row, args_spec, mappings)
    return Constructor(args...)
end

function call_yaml_constructor(
        Constructor,
        row::NamedTuple,
        args_spec::Vector{Symbol},
        kwargs_spec::Vector;
        mappings::Dict{Symbol, <:Function}=
            Dict{Symbol, Function}())
    args = _extract_args(row, args_spec, mappings)

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

function parse_dynamics_type(s::String)
    s_upper = uppercase(s)
    s_upper == "STATIC" && return STATIC
    s_upper == "DYNAMIC" && return DYNAMIC
    s_upper == "WING" && return WING
    s_upper == "QUASI_STATIC" && return QUASI_STATIC
    error("Unknown DynamicsType: $s")
end

function parse_segment_type(s::String)
    s_upper = uppercase(s)
    s_upper == "POWER_LINE" && return POWER_LINE
    s_upper == "STEERING_LINE" && return STEERING_LINE
    s_upper == "BRIDLE" && return BRIDLE
    error("Unknown SegmentType: $s")
end

function parse_wing_type(s::String)
    s_upper = uppercase(s)
    s_upper == "REFINE" && return REFINE
    s_upper == "QUATERNION" && return QUATERNION
    error("Unknown WingType: $s")
end

function parse_aero_mode(s::String)
    s_upper = uppercase(s)
    s_upper == "AERO_NONE" && return AERO_NONE
    s_upper == "AERO_DIRECT" && return AERO_DIRECT
    s_upper == "AERO_LINEARIZED" && return AERO_LINEARIZED
    error("Unknown AeroMode: $s")
end

"""
        load_sys_struct_from_yaml(yaml_path::AbstractString; system_name="from_yaml", set=nothing)

Build a `SystemStructure` from a component-based structural YAML file.

**IMPORTANT**: All indices (points, segments, etc.) must be sequential
starting from 1 with no gaps.

# Expected top-level blocks
- `points`: table with headers `[id,x,y,z,type,mass,body_damping,world_damping]`
  - `type`: STATIC, DYNAMIC, WING, or QUASI_STATIC
  - `id` must be sequential: 1, 2, 3, ...

- `segments`: table with one of two formats:
  - Direct format: `[id,point_i,point_j,type,l0,diameter_mm,unit_stiffness,unit_damping,compression_frac]`
  - Named format: `[name,point_i,point_j]` (requires `segment_properties` block)

- `segment_properties`: (optional) table with headers `[name,type,l0,diameter_mm,unit_stiffness,unit_damping,compression_frac]`
  - Used with named segment format for shared properties across symmetric segments

- `pulleys`: table with headers `[id,segment_i,segment_j,type]`
  - `type`: DYNAMIC or QUASI_STATIC

- `groups`: (optional) table with headers `[id,point_idxs,gamma,type,reference_chord_frac]`
- `tethers`: (optional) table with headers
  - Route 1 (explicit segments): `[name, segment_idxs]`, optional lengths
  - Route 2 (auto-generated): `[name, start_point, end_point, n_segments, material]`, optional lengths
  - `init_stretched_length`: scales `pos_w` before transforms [m]
  - `init_unstretched_length`: rope length for segment l0 [m]
- `winches`: (optional) table with headers `[name, tether_idxs, winch_point]`
- `wings`: (optional, typically from VSM configuration)
- `transforms`: (optional, typically from settings)
"""
function load_sys_struct_from_yaml(yaml_path::AbstractString; system_name="from_yaml", set::Union{Nothing,Settings}=nothing, ignore_l0::Bool=false, wing_type::Union{Nothing,WingType}=nothing, aero_mode::Union{Nothing,AeroMode}=nothing, vsm_set::Union{Nothing,VortexStepMethod.VSMSettings}=nothing)
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

    # Parse [a, b] or [a, [[id, w], ...]] style reference-point
    # fields from YAML rows. Weighted specs like
    #   [[2, 0.7], [4, 0.3]]
    # are converted to tuples for WeightedRefPoints.
    yaml_parse_ref_points = function (row, field)
        !hasfield(typeof(row), field) && return nothing
        val = getfield(row, field)
        val === nothing && return nothing

        if length(val) != 2
            throw(ArgumentError("ref_points must have 2 elements"))
        end

        convert_ref = function (v)
            if v isa Vector && !isempty(v) && v[1] isa Vector
                # Weighted: [[id, weight], ...] → tuples
                return map(v) do x
                    if !(x isa Vector) || length(x) != 2
                        throw(ArgumentError("Invalid weighted reference point entry $(repr(x)); expected format [[id, weight], ...]"))
                    end
                    if !(x[2] isa Number)
                        throw(ArgumentError("Invalid weighted reference point weight $(repr(x[2])) in entry $(repr(x)); expected format [[id, weight], ...] with numeric weight"))
                    end
                    (yaml_to_ref(x[1]), Float64(x[2]))
                end
            elseif v isa Vector
                # Multiple equal-weight refs: [a, b, ...]
                return [yaml_to_ref(x) for x in v]
            else
                return yaml_to_ref(v)
            end
        end
        p1 = convert_ref(val[1])
        p2 = convert_ref(val[2])
        return (p1, p2)
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
                    :pos_cad => r -> KVec3(r.pos_cad...),
                    :type => r -> parse_dynamics_type(String(r.type)),
                    :name => r -> haskey(r, :name) && !isnothing(r.name) ? Symbol(r.name) : i,
                    # Pass raw references - constructor handles defaults
                    :wing => r -> haskey(r, :wing_idx) ? yaml_to_ref(r.wing_idx) : nothing,
                    :transform => r -> haskey(r, :transform_idx) ? yaml_to_ref(r.transform_idx) : nothing,
                    :body_frame_damping => r -> haskey(r, :body_frame_damping) ? r.body_frame_damping : nothing,
                    :world_frame_damping => r -> haskey(r, :world_frame_damping) ? r.world_frame_damping : nothing
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
                    :set => r -> resolved_set,
                    :point_i => r -> yaml_to_ref(r.point_i),
                    :point_j => r -> yaml_to_ref(r.point_j),
                    :name => r -> begin
                        if haskey(r, :name) &&
                                !isnothing(r.name)
                            Symbol(r.name)
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
                    :segment_i => r -> yaml_to_ref(r.segment_i),
                    :segment_j => r -> yaml_to_ref(r.segment_j),
                    :name => r -> begin
                        if haskey(r, :name) && !isnothing(r.name)
                            Symbol(r.name)
                        else
                            i  # Use index as name if no name provided
                        end
                    end,
                    :type => r -> parse_dynamics_type(String(r.type))
                ))
            push!(pulleys, pulley)
        end
    end

    # Load groups (optional, for deformable wings) - SystemStructure handles resolution
    groups = Group[]
    if haskey(data, "groups") &&
       haskey(data["groups"], "data") &&
       data["groups"]["data"] !== nothing &&
       !isempty(data["groups"]["data"])
        group_rows = parse_table(data["groups"])

        for (i, row) in enumerate(group_rows)
            # Create Group using new constructor (name, points, type, moment_frac)
            group = call_yaml_constructor(Group, row,
                [:name, :points, :type, :moment_frac],
                [:damping];
                mappings=Dict(
                    :points => r -> [yaml_to_ref(p) for p in r.point_idxs],
                    :name => r -> begin
                        if haskey(r, :name) && !isnothing(r.name)
                            Symbol(r.name)
                        else
                            i  # Use index as name if no name provided
                        end
                    end,
                    :type => r -> parse_dynamics_type(
                        String(r.type))
                ))
            push!(groups, group)
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
                segs = [yaml_to_ref(s)
                    for s in row.segment_idxs]
                sp = if hasfield(typeof(row),
                        :start_point) &&
                        !isnothing(row.start_point)
                    yaml_to_ref(row.start_point)
                else
                    nothing
                end
                ep = if hasfield(typeof(row),
                        :end_point) &&
                        !isnothing(row.end_point)
                    yaml_to_ref(row.end_point)
                else
                    nothing
                end
                il = if hasfield(typeof(row),
                        :init_stretched_length) &&
                        !isnothing(
                            row.init_stretched_length)
                    Float64(row.init_stretched_length)
                else
                    nothing
                end
                tl = if hasfield(typeof(row),
                        :init_unstretched_length) &&
                        !isnothing(
                            row.init_unstretched_length)
                    Float64(
                        row.init_unstretched_length)
                else
                    nothing
                end
                # Default: unstretched = stretched
                ul = !isnothing(tl) ? tl :
                    !isnothing(il) ? il :
                    error("Tether $tether_name: " *
                        "init_unstretched_length " *
                        "or init_stretched_length " *
                        "is required")
                tether = Tether(tether_name, segs, ul;
                    start_point=sp, end_point=ep,
                    stretched_length=il)
            else
                # Route 2: auto-generation
                sp = yaml_to_ref(row.start_point)
                ep = yaml_to_ref(row.end_point)
                n_seg = Int(row.n_segments)
                il = if hasfield(typeof(row),
                        :init_stretched_length) &&
                        !isnothing(
                            row.init_stretched_length)
                    Float64(row.init_stretched_length)
                else
                    nothing
                end
                tl = if hasfield(typeof(row),
                        :init_unstretched_length) &&
                        !isnothing(
                            row.init_unstretched_length)
                    Float64(
                        row.init_unstretched_length)
                else
                    nothing
                end
                # Resolve material reference if present
                resolved = resolve_references(
                    row, property_tables)
                props = Dict{Symbol, Any}(
                    pairs(resolved))
                calculate_derived_properties!(props)
                us = get(props, :unit_stiffness,
                    NaN)
                us = isnothing(us) ? NaN :
                    Float64(us)
                ud = get(props, :unit_damping,
                    NaN)
                ud = isnothing(ud) ? NaN :
                    Float64(ud)
                d_mm = get(props, :diameter_mm,
                    NaN)
                d_mm = isnothing(d_mm) ? NaN :
                    Float64(d_mm)
                d = isnan(d_mm) ? NaN :
                    d_mm * 0.001
                # Default: unstretched = stretched
                ul = !isnothing(tl) ? tl :
                    !isnothing(il) ? il :
                    error("Tether $tether_name: " *
                        "init_unstretched_length " *
                        "or init_stretched_length " *
                        "is required")
                tether = Tether(tether_name, ul;
                    start_point=sp, end_point=ep,
                    n_segments=n_seg,
                    unit_stiffness=us, unit_damping=ud,
                    diameter=d,
                    stretched_length=il)
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
                 :brake, :speed_controlled,
                 :friction_epsilon];
                mappings=Dict(
                    :set => r -> resolved_set,
                    :tethers => r -> [yaml_to_ref(t)
                        for t in r.tether_idxs],
                    :winch_point => r -> begin
                        yaml_to_ref(r.winch_point)
                    end,
                    :name => r -> begin
                        if haskey(r, :name) &&
                           !isnothing(r.name)
                            Symbol(r.name)
                        else
                            i
                        end
                    end
                ))
            push!(winches, winch)
        end
    end

    # Load wings (optional)
    wings = VSMWing[]
    if haskey(data, "wings") &&
       haskey(data["wings"], "data") &&
       data["wings"]["data"] !== nothing &&
       !isempty(data["wings"]["data"])
        wing_rows = parse_table(data["wings"])

        # Validate vsm_set is provided when wings are defined
        if isnothing(vsm_set)
            error("Wings are defined in YAML but vsm_set was not provided to load_sys_struct_from_yaml. " *
                  "Please pass a VortexStepMethod.VSMSettings object via the vsm_set keyword argument.")
        end

        for (i, row) in enumerate(wing_rows)
            # Use provided wing_type parameter or parse from YAML
            wt = isnothing(wing_type) ? parse_wing_type(String(row.type)) : wing_type

            # Build kwargs based on wing type - SystemStructure handles resolution
            # Determine aero_mode: kwarg > YAML > default
            am = if !isnothing(aero_mode)
                aero_mode
            elseif hasfield(typeof(row), :aero_mode) &&
                    !isnothing(row.aero_mode)
                parse_aero_mode(String(row.aero_mode))
            else
                wt == QUATERNION ? AERO_LINEARIZED :
                    AERO_DIRECT
            end

            if wt == REFINE
                # REFINE wings need z_ref_points, y_ref_points, origin
                # Pass raw values - constructor handles defaults
                wing = call_yaml_constructor(VSMWing, row,
                    [:name, :set, :groups, :vsm_set],
                    [:transform, :y_damping, :angular_damping,
                     :wing_type, :aero_mode,
                     :z_ref_points, :y_ref_points, :origin, :pos_cad,
                     :aero_scale_chord];
                    mappings=Dict(
                        :set => r -> resolved_set,
                        :groups => r -> [],  # REFINE wings don't have groups
                        :vsm_set => r -> vsm_set,
                        :wing_type => r -> wt,
                        :aero_mode => r -> am,
                        :name => r -> begin
                            if haskey(r, :name) && !isnothing(r.name)
                                Symbol(r.name)
                            else
                                i  # Use index as name if no name provided
                            end
                        end,
                        :transform => r -> begin
                            if hasfield(typeof(r), :transform_idx) && !isnothing(r.transform_idx)
                                yaml_to_ref(r.transform_idx)
                            else
                                nothing  # Constructor handles default
                            end
                        end,
                        :z_ref_points => r ->
                            yaml_parse_ref_points(r, :z_ref_points),
                        :y_ref_points => r ->
                            yaml_parse_ref_points(r, :y_ref_points),
                        :origin => r -> begin
                            if !hasfield(typeof(r), :origin_idx) || r.origin_idx === nothing
                                return nothing
                            end
                            yaml_to_ref(r.origin_idx)
                        end,
                        :aero_scale_chord => r ->
                            hasfield(typeof(r), :aero_scale_chord) && !isnothing(r.aero_scale_chord) ?
                                float(r.aero_scale_chord) : 0.0,
                        :pos_cad => r -> begin
                            # Note: pos_cad will be set from origin point position after resolution
                            # For now, return nothing - SystemStructure will handle this
                            nothing
                        end
                    ))
            else  # QUATERNION
                # Pass raw values - constructor handles defaults
                wing = call_yaml_constructor(VSMWing, row,
                    [:name, :set, :groups, :vsm_set],
                    [:transform, :y_damping, :angular_damping,
                     :wing_type, :aero_mode, :aero_scale_chord,
                     :aero_z_offset, :pos_cad,
                     :z_ref_points, :y_ref_points, :origin];
                    mappings=Dict(
                        :set => r -> resolved_set,
                        :aero_mode => r -> am,
                        :groups => r -> hasfield(typeof(r), :groups) &&
                            !isnothing(r.groups) ?
                            [yaml_to_ref(g) for g in r.groups] : [],
                        :vsm_set => r -> vsm_set,
                        :wing_type => r -> wt,
                        :name => r -> begin
                            if haskey(r, :name) && !isnothing(r.name)
                                Symbol(r.name)
                            else
                                i
                            end
                        end,
                        :transform => r -> begin
                            if hasfield(typeof(r), :transform_idx) &&
                               !isnothing(r.transform_idx)
                                yaml_to_ref(r.transform_idx)
                            else
                                nothing
                            end
                        end,
                        :pos_cad => r -> begin
                            if !hasfield(typeof(r), :pos_cad) ||
                               r.pos_cad === nothing
                                return nothing
                            end
                            KVec3(r.pos_cad...)
                        end,
                        :aero_scale_chord => r ->
                            hasfield(typeof(r), :aero_scale_chord) &&
                            !isnothing(r.aero_scale_chord) ?
                                float(r.aero_scale_chord) : 0.0,
                        :z_ref_points => r ->
                            yaml_parse_ref_points(r, :z_ref_points),
                        :y_ref_points => r ->
                            yaml_parse_ref_points(r, :y_ref_points),
                        :origin => r -> begin
                            if !hasfield(typeof(r), :origin_idx) ||
                               r.origin_idx === nothing
                                return nothing
                            end
                            yaml_to_ref(r.origin_idx)
                        end
                    ))
            end
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
                    :elevation => r -> deg2rad(r.elevation),
                    :azimuth => r -> deg2rad(r.azimuth),
                    :heading => r -> deg2rad(r.heading),
                    :elevation_vel => r -> hasfield(typeof(r), :elevation_vel) && !isnothing(r.elevation_vel) ?
                        deg2rad(r.elevation_vel) : 0.0,
                    :azimuth_vel => r -> hasfield(typeof(r), :azimuth_vel) && !isnothing(r.azimuth_vel) ?
                        deg2rad(r.azimuth_vel) : 0.0,
                    :turn_rate => r -> hasfield(typeof(r), :turn_rate) && !isnothing(r.turn_rate) ?
                        deg2rad(r.turn_rate) : 0.0,
                    :base_pos => r -> KVec3(r.base_pos...),
                    :base_point => r -> yaml_to_ref(r.base_point_idx),
                    :base_transform => r -> begin
                        if hasfield(typeof(r), :base_transform_idx) && !isnothing(r.base_transform_idx)
                            yaml_to_ref(r.base_transform_idx)
                        else
                            nothing
                        end
                    end,
                    :rot_point => r -> begin
                        if hasfield(typeof(r), :rot_point_idx) && !isnothing(r.rot_point_idx)
                            yaml_to_ref(r.rot_point_idx)
                        else
                            nothing
                        end
                    end,
                    :wing => r -> begin
                        if !hasfield(typeof(r), :wing_idx) || r.wing_idx === nothing
                            return nothing
                        end
                        yaml_to_ref(r.wing_idx)
                    end,
                    :name => r -> begin
                        if haskey(r, :name) && !isnothing(r.name)
                            Symbol(r.name)
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
    return SystemStructure(system_name, resolved_set; points, groups,
        segments, pulleys, tethers, winches, wings,
        transforms, ignore_l0, vsm_set)
end

"""
    update_yaml_from_sys_struct!(sys_struct::SystemStructure,
                                  source_struc_yaml::AbstractString,
                                  dest_struc_yaml::AbstractString,
                                  source_aero_yaml::AbstractString,
                                  dest_aero_yaml::AbstractString)

Update point positions in structural and aerodynamic YAML files from
the current state of a SystemStructure.

# Arguments
- `sys_struct`: SystemStructure with current point positions
- `source_struc_yaml`: Path to source structural geometry YAML file
- `dest_struc_yaml`: Path to destination structural YAML file
- `source_aero_yaml`: Path to source aero geometry YAML file
- `dest_aero_yaml`: Path to destination aero YAML file

Source and destination paths must be different for each pair.

# Example
```julia
sys = load_sys_struct_from_yaml("refine_struc_geometry.yaml"; ...)
sam = SymbolicAWEModel(set, sys)
# ... run simulation ...
update_yaml_from_sys_struct!(sys,
    "refine_struc_geometry.yaml",
    "refine_struc_geometry_stable.yaml",
    "aero_geometry.yaml",
    "aero_geometry_stable.yaml")
```
"""
function update_yaml_from_sys_struct!(sys_struct::SystemStructure,
                                      source_struc_yaml::AbstractString,
                                      dest_struc_yaml::AbstractString,
                                      source_aero_yaml::AbstractString,
                                      dest_aero_yaml::AbstractString)
    # Helper to format coordinate with rounding
    function format_coord(val::Float64)
        # Round to 4 decimals, set small values to 0
        rounded = abs(val) < 1e-4 ? 0.0 : round(val, digits=4)
        return rounded
    end

    # Update pos_b for REFINE wing points based on current wing orientation
    for wing in sys_struct.wings
        if wing.wing_type == REFINE
            R_w_to_b = wing.R_b_to_w'  # transpose to get world-to-body
            for point in sys_struct.points
                if point.wing_idx == wing.idx
                    point.pos_b .= R_w_to_b * (point.pos_w - wing.pos_w)
                end
            end
        end
    end

    # Build position dictionary from system structure (body-frame positions)
    positions = Dict{Int, Vector{Float64}}()
    for point in sys_struct.points
        positions[point.idx] = copy(point.pos_b)
    end

    # Build segment l0 dictionary from system structure
    segment_l0s = Dict{Int, Float64}()
    for seg in sys_struct.segments
        segment_l0s[seg.idx] = seg.l0
    end

    # Update structural geometry YAML
    struc_full_path = isabspath(source_struc_yaml) ? source_struc_yaml :
                      joinpath(pwd(), source_struc_yaml)

    if !isfile(struc_full_path)
        error("Source structural YAML file not found: $struc_full_path")
    end

    dest_struc_full_path = isabspath(dest_struc_yaml) ? dest_struc_yaml :
                          joinpath(pwd(), dest_struc_yaml)

    lines = readlines(struc_full_path)
    n_points_updated = 0
    n_segments_updated = 0
    in_points_section = false
    in_segments_section = false

    for (i, line) in enumerate(lines)
        # Track which section we're in
        if occursin(r"^points:", line)
            in_points_section = true
            in_segments_section = false
        elseif occursin(r"^segments:", line)
            in_points_section = false
            in_segments_section = true
        elseif occursin(r"^\w+:", line)  # New section starts
            in_points_section = false
            in_segments_section = false
        end

        # Update lines in the points section
        if in_points_section
            # Match: "- [idx, [x, y, z], ..." where coordinates are floats
            m = match(r"^(\s*-\s*\[)(\d+)(,\s*\[)([-+]?\d+\.?\d*(?:[eE][-+]?\d+)?,\s*[-+]?\d+\.?\d*(?:[eE][-+]?\d+)?,\s*[-+]?\d+\.?\d*(?:[eE][-+]?\d+)?)(\].*)", line)
            if m !== nothing
                point_idx = parse(Int, something(m.captures[2]))
                if haskey(positions, point_idx)
                    new_pos = positions[point_idx]
                    x = format_coord(new_pos[1])
                    y = format_coord(new_pos[2])
                    z = format_coord(new_pos[3])
                    new_coords = "$x, $y, $z"
                    lines[i] = something(m.captures[1]) * something(m.captures[2]) *
                              something(m.captures[3]) * new_coords * something(m.captures[5])
                    n_points_updated += 1
                end
            end
        end

        # Update lines in the segments section
        if in_segments_section
            # Match: "- [idx, point_i, point_j, l0, ...]"
            # Format: [idx, point_i, point_j, l0, diameter_mm, ...]
            # We want to update the l0 field (4th field)
            m = match(r"^(\s*-\s*\[)(\d+)(,\s*\d+,\s*\d+,\s*)([-+]?\d+\.?\d*(?:[eE][-+]?\d+)?)(.*)", line)
            if m !== nothing
                seg_idx = parse(Int, something(m.captures[2]))

                if haskey(segment_l0s, seg_idx)
                    new_l0 = format_coord(segment_l0s[seg_idx])
                    # Reconstruct line with updated l0
                    lines[i] = something(m.captures[1]) * something(m.captures[2]) * something(m.captures[3]) *
                              string(new_l0) * something(m.captures[5])
                    n_segments_updated += 1
                end
            end
        end
    end

    # Build tether init_len dictionary
    tether_init_lens = Dict{String, Float64}()
    for tether in sys_struct.tethers
        name = string(tether.name)
        if !isnothing(tether.init_stretched_len)
            tether_init_lens[name] =
                tether.init_stretched_len
        end
    end

    # Update tether init_len values in tethers section
    n_tethers_updated = 0
    in_tethers_section = false
    tether_init_len_col = 0
    for (i, line) in enumerate(lines)
        if occursin(r"^tethers:", line)
            in_tethers_section = true
        elseif occursin(r"^\w+:", line) &&
               !occursin(r"^tethers:", line)
            if in_tethers_section
                in_tethers_section = false
                tether_init_len_col = 0
            end
        end

        if in_tethers_section
            # Find init_len column index from headers
            hm = match(r"headers:\s*\[(.+)\]", line)
            if hm !== nothing
                cols = split(something(hm.captures[1]), r",\s*")
                for (ci, c) in enumerate(cols)
                    if strip(c) == "init_stretched_length"
                        tether_init_len_col = ci
                    end
                end
            end

            # Update data rows if we know the column
            if tether_init_len_col > 0
                dm = match(
                    r"^(\s*-\s*\[)([\w-]+)(.*)", line)
                if dm !== nothing
                    name = String(something(dm.captures[2]))
                    if haskey(tether_init_lens, name)
                        # Parse fields, update init_len
                        rest = something(dm.captures[3])
                        fields = split(
                            strip(rest, [',', ']']),
                            r",\s*")
                        col = tether_init_len_col - 1
                        if col <= length(fields)
                            fields[col] = " " * string(
                                format_coord(
                                    tether_init_lens[
                                        name]))
                            lines[i] = dm.captures[1] *
                                name * "," *
                                join(fields, ",") * "]"
                            n_tethers_updated += 1
                        end
                    end
                end
            end
        end
    end

    @info "Updated structural positions and segments" n_points=n_points_updated n_segments=n_segments_updated n_tethers=n_tethers_updated

    # Write updated structural YAML
    @info "Writing updated structural YAML" source=struc_full_path dest=dest_struc_full_path
    open(dest_struc_full_path, "w") do io
        for line in lines
            println(io, line)
        end
    end

    # Update aerodynamic geometry YAML
    aero_full_path = isabspath(source_aero_yaml) ? source_aero_yaml :
                    joinpath(pwd(), source_aero_yaml)

    if !isfile(aero_full_path)
        error("Source aero YAML file not found: $aero_full_path")
    end

    dest_aero_full_path = isabspath(dest_aero_yaml) ? dest_aero_yaml :
                         joinpath(pwd(), dest_aero_yaml)

    aero_lines = readlines(aero_full_path)
    n_aero_updated = 0

    for (i, line) in enumerate(aero_lines)
        # Match wing section data lines with point references in comments
        # Format: "- [airfoil_id, LE_x, LE_y, LE_z, TE_x, TE_y, TE_z]"
        # Look for comment indicating point mapping
        comment_match = match(r"#.*points?\s+(\d+).*\(LE\).*and\s+(\d+).*\(TE\)", line)

        if comment_match !== nothing
            le_idx = parse(Int, something(comment_match.captures[1]))
            te_idx = parse(Int, something(comment_match.captures[2]))

            # Check next line for the actual data
            if i < length(aero_lines)
                data_line = aero_lines[i+1]
                m = match(r"^(\s*-\s*\[)(\d+)(,\s*)([-+]?\d+\.?\d*(?:[eE][-+]?\d+)?),\s*([-+]?\d+\.?\d*(?:[eE][-+]?\d+)?),\s*([-+]?\d+\.?\d*(?:[eE][-+]?\d+)?),\s*([-+]?\d+\.?\d*(?:[eE][-+]?\d+)?),\s*([-+]?\d+\.?\d*(?:[eE][-+]?\d+)?),\s*([-+]?\d+\.?\d*(?:[eE][-+]?\d+)?)(\])", data_line)

                if m !== nothing && haskey(positions, le_idx) && haskey(positions, te_idx)
                    le_pos = positions[le_idx]
                    te_pos = positions[te_idx]

                    airfoil_id = m.captures[2]
                    new_data = "$(m.captures[1])$airfoil_id$(m.captures[3])" *
                              "$(format_coord(le_pos[1])), $(format_coord(le_pos[2])), $(format_coord(le_pos[3])), " *
                              "$(format_coord(te_pos[1])), $(format_coord(te_pos[2])), $(format_coord(te_pos[3]))$(m.captures[10])"

                    aero_lines[i+1] = new_data
                    n_aero_updated += 1
                end
            end
        end
    end

    @info "Updated aerodynamic positions" n_sections=n_aero_updated

    # Write updated aero YAML
    @info "Writing updated aero YAML" source=aero_full_path dest=dest_aero_full_path
    open(dest_aero_full_path, "w") do io
        for line in aero_lines
            println(io, line)
        end
    end

    @info "Done!"
    return nothing
end

"""
    update_sys_struct_from_yaml!(sys_struct::SystemStructure,
                                  struc_yaml::AbstractString)

Update an existing `SystemStructure` in-place from a (possibly
modified) structural geometry YAML file. Inverse of
`update_yaml_from_sys_struct!`.

Updates `pos_cad` for points, `l0` for segments, and
`init_stretched_len`/`init_unstretched_len` for tethers,
matched by symbolic name. When `l0` is `nothing` in the
YAML, it is auto-calculated from endpoint `pos_cad`.

Only raw geometry is updated. Call `reinit!(sys_struct, set)` afterward
to recompute derived quantities (`pos_b`, `pos_w`, wing frames, etc.).

Unmatched names are silently skipped (the YAML may contain a subset of
components).

# Arguments
- `sys_struct`: The SystemStructure to update in-place.
- `struc_yaml`: Path to the structural geometry YAML file.

# Example
```julia
sys = load_sys_struct_from_yaml("refine_struc_geometry.yaml"; ...)
# ... edit YAML externally ...
update_sys_struct_from_yaml!(sys, "refine_struc_geometry.yaml")
```
"""
function update_sys_struct_from_yaml!(
        sys_struct::SystemStructure,
        struc_yaml::AbstractString)
    yaml_path = isabspath(struc_yaml) ? struc_yaml :
                joinpath(pwd(), struc_yaml)
    isfile(yaml_path) ||
        error("YAML file not found: $yaml_path")

    data = YAML.load_file(yaml_path)

    # --- Update points ---
    n_points = 0
    if haskey(data, "points")
        point_rows = parse_table(data["points"])
        for row in point_rows
            haskey(row, :name) || continue
            name = Symbol(row.name)
            haskey(sys_struct.points, name) || continue

            point = sys_struct.points[name]
            point.pos_cad .= KVec3(row.pos_cad...)
            n_points += 1
        end
    end

    # --- Update segment l0 ---
    n_segments = 0
    if haskey(data, "segments")
        segment_rows = parse_table(data["segments"])
        for row in segment_rows
            haskey(row, :name) || continue
            name = Symbol(row.name)
            haskey(sys_struct.segments, name) || continue

            seg = sys_struct.segments[name]

            # l0: use YAML value, or auto-calc from pos_cad
            l0_val = haskey(row, :l0) ? row.l0 : nothing
            if !isnothing(l0_val) && l0_val != "nothing"
                seg.l0 = Float64(l0_val)
            else
                seg.l0 = segment_cad_length(
                    seg, sys_struct.points)
            end

            n_segments += 1
        end
    end

    # --- Update tether init lengths ---
    n_tethers = 0
    if haskey(data, "tethers")
        tether_rows = parse_table(data["tethers"])
        for row in tether_rows
            haskey(row, :name) || continue
            name = Symbol(row.name)
            haskey(sys_struct.tethers, name) || continue

            tether = sys_struct.tethers[name]

            if hasfield(typeof(row),
                    :init_stretched_length) &&
               !isnothing(row.init_stretched_length)
                tether.init_stretched_len =
                    Float64(row.init_stretched_length)
            end
            if hasfield(typeof(row),
                    :init_unstretched_length) &&
               !isnothing(
                    row.init_unstretched_length)
                tether.init_unstretched_len =
                    Float64(
                        row.init_unstretched_length)
            end

            n_tethers += 1
        end
    end

    @info "update_sys_struct_from_yaml!" n_points n_segments n_tethers
    return nothing
end
