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
    s_upper == "PARTICLE_DYNAMICS" && return PARTICLE_DYNAMICS
    s_upper == "RIGID_DYNAMICS" && return RIGID_DYNAMICS
    if s_upper == "REFINE"
        @warn "WingType \"$s\" is deprecated; use \"PARTICLE_DYNAMICS\" instead."
        return PARTICLE_DYNAMICS
    end
    if s_upper == "QUATERNION"
        @warn "WingType \"$s\" is deprecated; use \"RIGID_DYNAMICS\" instead."
        return RIGID_DYNAMICS
    end
    error("Unknown WingType: $s")
end

function parse_aero_mode(s::String)
    s_upper = uppercase(s)
    s_upper == "AERO_NONE" && return AERO_NONE
    s_upper == "AERO_DIRECT" && return AERO_DIRECT
    s_upper == "AERO_LINEARIZED" && return AERO_LINEARIZED
    s_upper == "AERO_PLATE" && return AERO_PLATE
    error("Unknown AeroMode: $s")
end

"""
    _load_plate_wing(row, idx, data, set, wt, am,
                     yaml_to_ref, yaml_parse_ref_points)

Load a PlateWing from YAML wing row + surfaces block.
CL/CD interpolations are created from Settings polar data.
"""
function _load_plate_wing(row, idx, data, set, wt, am,
                          yaml_to_ref, yaml_parse_ref_points,
                          yaml_parse_origin)
    name = if haskey(row, :name) && !isnothing(row.name)
        Symbol(row.name)
    else
        idx
    end

    # CL/CD from settings polar data
    cl_interp, cd_interp = create_plate_interpolations(
        set.alpha_cl, set.cl_list, set.cd_list;
        alpha_cd=set.alpha_cd)

    # Parse wing-level parameters
    drag_corr = hasfield(typeof(row), :drag_corr) &&
        !isnothing(row.drag_corr) ? float(row.drag_corr) : 0.93
    cmq = hasfield(typeof(row), :cmq) &&
        !isnothing(row.cmq) ? float(row.cmq) : 0.0
    smc = hasfield(typeof(row), :smc) &&
        !isnothing(row.smc) ? float(row.smc) : 0.0
    cord_length = hasfield(typeof(row), :cord_length) &&
        !isnothing(row.cord_length) ?
        float(row.cord_length) : 1.0
    y_damping = hasfield(typeof(row), :y_damping) &&
        !isnothing(row.y_damping) ?
        float(row.y_damping) : 150.0

    # Parse reference points
    z_ref = yaml_parse_ref_points(row, :z_ref_points)
    y_ref = yaml_parse_ref_points(row, :y_ref_points)
    origin = yaml_parse_origin(row, :origin_idx)
    transform = if hasfield(typeof(row), :transform_idx) &&
                   !isnothing(row.transform_idx)
        yaml_to_ref(row.transform_idx)
    else
        nothing
    end

    # Load surfaces from YAML
    surfaces = PlateSurface[]
    if haskey(data, "surfaces") &&
       haskey(data["surfaces"], "data") &&
       data["surfaces"]["data"] !== nothing
        surf_rows = parse_table(data["surfaces"])
        for (si, sr) in enumerate(surf_rows)
            sname = haskey(sr, :name) && !isnothing(sr.name) ?
                Symbol(sr.name) : nothing
            x_airf = KVec3(sr.x_airf...)
            y_airf = KVec3(sr.y_airf...)
            area = float(sr.area)
            point = yaml_to_ref(sr.point_idx)
            twist = hasfield(typeof(sr), :twist) &&
                !isnothing(sr.twist) ?
                float(sr.twist) : 0.0
            push!(surfaces, PlateSurface(
                sname, x_airf, y_airf, area, point;
                twist))
        end
    end

    PlateWing(name, surfaces, cl_interp, cd_interp;
              dynamics_type=wt, transform, y_damping,
              drag_corr, cmq, smc, cord_length,
              z_ref_points=z_ref, y_ref_points=y_ref,
              origin)
end

"""
    parse_tether_init(row, tether_name) -> (stretched_len, force)

Read `init_stretched_length` and `init_tether_force` from a tether
YAML row (each `nothing` if absent). Errors if the deprecated
`init_unstretched_length` field is present — the unstretched rest
length is now derived from the placed stretched length and
`init_tether_force`.
"""
function parse_tether_init(row, tether_name)
    if hasfield(typeof(row), :init_unstretched_length) &&
       !isnothing(row.init_unstretched_length)
        error("Tether $tether_name: init_unstretched_length is " *
              "deprecated; it is derived from the placed " *
              "init_stretched_length and init_tether_force.")
    end
    isl = hasfield(typeof(row), :init_stretched_length) &&
        !isnothing(row.init_stretched_length) ?
        Float64(row.init_stretched_length) : nothing
    itf = hasfield(typeof(row), :init_tether_force) &&
        !isnothing(row.init_tether_force) ?
        Float64(row.init_tether_force) : nothing
    return isl, itf
end

"""
    load_sys_struct_from_yaml(yaml_path; system_name, set, ...)

Build a `SystemStructure` from a component-based structural
YAML file. See source for full documentation of expected blocks.
"""
function load_sys_struct_from_yaml(yaml_path::AbstractString; system_name="from_yaml", set::Union{Nothing,Settings}=nothing, ignore_l0::Bool=false, dynamics_type::Union{Nothing,WingType}=nothing, aero_mode::Union{Nothing,AeroMode}=nothing, vsm_set::Union{Nothing,VortexStepMethod.VSMSettings}=nothing, wing_type::Union{Nothing,WingType}=nothing)
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
    yaml_convert_ref = function (v)
        if v isa Vector && !isempty(v) && v[1] isa Vector
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
            return [yaml_to_ref(x) for x in v]
        else
            return yaml_to_ref(v)
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

        p1 = yaml_convert_ref(val[1])
        p2 = yaml_convert_ref(val[2])
        return (p1, p2)
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
                isl, itf = parse_tether_init(
                    row, tether_name)
                tether = Tether(tether_name, segs, isl;
                    start_point=sp, end_point=ep,
                    tether_force=itf)
            else
                # Route 2: auto-generation
                sp = yaml_to_ref(row.start_point)
                ep = yaml_to_ref(row.end_point)
                n_seg = Int(row.n_segments)
                isl, itf = parse_tether_init(
                    row, tether_name)
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
                tether = Tether(tether_name, isl;
                    start_point=sp, end_point=ep,
                    n_segments=n_seg,
                    unit_stiffness=us, unit_damping=ud,
                    diameter=d, tether_force=itf)
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
    wings = AbstractWing[]
    if haskey(data, "wings") &&
       haskey(data["wings"], "data") &&
       data["wings"]["data"] !== nothing &&
       !isempty(data["wings"]["data"])
        wing_rows = parse_table(data["wings"])

        for (i, row) in enumerate(wing_rows)
            # Use provided dynamics_type parameter or parse from YAML
            # Support old `type` field with deprecation warning
            wt = if !isnothing(dynamics_type)
                dynamics_type
            else
                raw_wt_field = if hasfield(typeof(row), :dynamics_type) && !isnothing(row.dynamics_type)
                    String(row.dynamics_type)
                elseif hasfield(typeof(row), :type) && !isnothing(row.type)
                    @warn "Wing YAML field `type` is deprecated; rename to `dynamics_type`."
                    String(row.type)
                else
                    error("Wing entry missing required `dynamics_type` field.")
                end
                parse_wing_type(raw_wt_field)
            end

            # Build kwargs based on wing type - SystemStructure handles resolution
            # Determine aero_mode: kwarg > YAML > default
            am = if !isnothing(aero_mode)
                aero_mode
            elseif hasfield(typeof(row), :aero_mode) &&
                    !isnothing(row.aero_mode)
                parse_aero_mode(String(row.aero_mode))
            else
                wt == RIGID_DYNAMICS ? AERO_LINEARIZED :
                    AERO_DIRECT
            end

            if am == AERO_PLATE
                # PlateWing — load surfaces and CL/CD from settings
                wing = _load_plate_wing(row, i, data,
                    resolved_set, wt, am, yaml_to_ref,
                    yaml_parse_ref_points,
                    yaml_parse_origin)
                push!(wings, wing)
                continue
            end

            # VSMWing — validate vsm_set
            if isnothing(vsm_set)
                error("VSMWing defined in YAML but vsm_set " *
                      "was not provided.")
            end

            if wt == PARTICLE_DYNAMICS
                # PARTICLE_DYNAMICS wings need z_ref_points, y_ref_points, origin
                # Pass raw values - constructor handles defaults
                wing = call_yaml_constructor(VSMWing, row,
                    [:name, :set, :groups, :vsm_set],
                    [:transform, :y_damping, :angular_damping,
                     :dynamics_type, :aero_mode,
                     :z_ref_points, :y_ref_points, :origin, :pos_cad,
                     :aero_scale_chord];
                    mappings=Dict(
                        :set => r -> resolved_set,
                        :groups => r -> hasfield(typeof(r), :groups) &&
                            !isnothing(r.groups) ?
                            [yaml_to_ref(g) for g in r.groups] : [],
                        :vsm_set => r -> vsm_set,
                        :dynamics_type => r -> wt,
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
                        :origin => r ->
                            yaml_parse_origin(r, :origin_idx),
                        :aero_scale_chord => r ->
                            hasfield(typeof(r), :aero_scale_chord) && !isnothing(r.aero_scale_chord) ?
                                float(r.aero_scale_chord) : 0.0,
                        :pos_cad => r -> begin
                            # Note: pos_cad will be set from origin point position after resolution
                            # For now, return nothing - SystemStructure will handle this
                            nothing
                        end
                    ))
            else  # RIGID_DYNAMICS
                # Pass raw values - constructor handles defaults
                wing = call_yaml_constructor(VSMWing, row,
                    [:name, :set, :groups, :vsm_set],
                    [:transform, :y_damping, :angular_damping,
                     :dynamics_type, :aero_mode, :aero_scale_chord,
                     :aero_z_offset, :pos_cad,
                     :z_ref_points, :y_ref_points, :origin];
                    mappings=Dict(
                        :set => r -> resolved_set,
                        :aero_mode => r -> am,
                        :groups => r -> hasfield(typeof(r), :groups) &&
                            !isnothing(r.groups) ?
                            [yaml_to_ref(g) for g in r.groups] : [],
                        :vsm_set => r -> vsm_set,
                        :dynamics_type => r -> wt,
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
                        :origin => r ->
                            yaml_parse_origin(r, :origin_idx)
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
