# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

"""
Logging functionality for SystemStructure objects.

This module provides efficient logging capabilities for recording SystemStructure
states over time, with automatic field detection, multiple output formats, and
persistence features.
"""

using Serialization
using Dates

# ==================== LOGGING DATA STRUCTURES ==================== #

"""
    LoggingConfig

Configuration for controlling which fields get logged.

$(TYPEDFIELDS)
"""
struct LoggingConfig
    include_fields::Union{Vector{String}, Nothing}     # Specific fields to include
    exclude_fields::Union{Vector{String}, Nothing}     # Specific fields to exclude
    include_patterns::Union{Vector{Regex}, Nothing}    # Regex patterns to include
    exclude_patterns::Union{Vector{Regex}, Nothing}    # Regex patterns to exclude
end

"""
    LoggingConfig(; include_fields=nothing, exclude_fields=nothing,
                    include_patterns=nothing, exclude_patterns=nothing)

Create a logging configuration. If all parameters are `nothing`, all fields are logged.

# Examples
```julia
# Log only position and velocity
config = LoggingConfig(include_patterns=[r"\\.pos_w", r"\\.vel_w"])

# Exclude internal fields
config = LoggingConfig(exclude_patterns=[r"\\.internal_", r"\\.cache_"])

# Specific field selection
config = LoggingConfig(include_fields=["points.pos_w[1]", "points.mass"])
```
"""
LoggingConfig(; include_fields=nothing, exclude_fields=nothing,
                include_patterns=nothing, exclude_patterns=nothing) =
    LoggingConfig(include_fields, exclude_fields, include_patterns, exclude_patterns)

"""
    SystemLogger

A logging system for efficiently recording SystemStructure states over time.

The logger uses pre-allocated memory for efficient data collection and leverages
the existing reflection-based field collection system for automatic adaptation
to structural changes.

$(TYPEDFIELDS)
"""
mutable struct SystemLogger
    const data::Matrix{SimFloat}        # Pre-allocated matrix: [time_steps × variables]
    const time_stamps::Vector{SimFloat} # Time stamps for each logged step
    const metadata::Dict{String, Any}   # Field names, types, sizes for reconstruction
    current_step::Int                   # Current logging position
    capacity::Int                      # Maximum number of time steps
end

# ==================== FIELD FILTERING ==================== #

"""
    _collect_logged_fields(sys::SystemStructure, fields_metadata::Vector)

Collect only the specific fields listed in the metadata, in the same order.
This ensures the collected data matches exactly what the logger expects.
"""
function _collect_logged_fields(sys::SystemStructure, fields_metadata::Vector)
    # Pre-allocate result vector
    total_size = sum(meta["size"] for meta in fields_metadata)
    result = Vector{SimFloat}(undef, total_size)

    # Create a map from component type to objects for efficient lookup
    component_map = Dict(
        "points" => sys.points,
        "groups" => sys.groups,
        "segments" => sys.segments,
        "pulleys" => sys.pulleys,
        "tethers" => sys.tethers,
        "winches" => sys.winches,
        "wings" => sys.wings,
        "transforms" => sys.transforms,
        "system" => [sys]
    )

    current_idx = 1

    for field_meta in fields_metadata
        field_name = field_meta["name"]
        field_size = field_meta["size"]

        # Parse field name like "points.pos_w[1]" or "points.mass"
        parts = split(field_name, ".")
        component_type = parts[1]
        field_part = parts[2]

        objects = component_map[component_type]

        # Check if this is a vector component or scalar field
        if occursin("[", field_part)
            # Vector component: "pos_w[1]"
            bracket_match = match(r"(\w+)\[(\d+)\]", field_part)
            if bracket_match !== nothing
                field_name_base = bracket_match.captures[1]
                component_idx = parse(Int, bracket_match.captures[2])

                # Collect this component from all objects
                for obj in objects
                    value = getfield(obj, Symbol(field_name_base))
                    result[current_idx] = value[component_idx]
                    current_idx += 1
                end
            end
        else
            # Scalar field: "mass"
            field_symbol = Symbol(field_part)
            for obj in objects
                result[current_idx] = getfield(obj, field_symbol)
                current_idx += 1
            end
        end
    end

    return result
end

"""
    should_log_field(field_name::String, config::LoggingConfig)

Determine if a field should be logged based on the configuration.

Returns `true` if the field should be logged, `false` otherwise.
"""
function should_log_field(field_name::String, config::LoggingConfig)
    # If no configuration provided, log everything
    if config.include_fields === nothing && config.exclude_fields === nothing &&
       config.include_patterns === nothing && config.exclude_patterns === nothing
        return true
    end

    # Check explicit exclude list first
    if config.exclude_fields !== nothing && field_name in config.exclude_fields
        return false
    end

    # Check exclude patterns
    if config.exclude_patterns !== nothing
        for pattern in config.exclude_patterns
            if occursin(pattern, field_name)
                return false
            end
        end
    end

    # If include lists are specified, field must match at least one
    has_include_criteria = config.include_fields !== nothing || config.include_patterns !== nothing
    if has_include_criteria
        # Check explicit include list
        if config.include_fields !== nothing && field_name in config.include_fields
            return true
        end

        # Check include patterns
        if config.include_patterns !== nothing
            for pattern in config.include_patterns
                if occursin(pattern, field_name)
                    return true
                end
            end
        end

        # No include criteria matched
        return false
    end

    # No exclude criteria matched and no include criteria specified
    return true
end

"""
    list_available_fields(sys::SystemStructure)

List all fields that could potentially be logged from a SystemStructure.

Useful for configuring selective logging.
"""
function list_available_fields(sys::SystemStructure)
    all_fields = String[]

    # Helper function to collect field names
    function collect_field_names(objects, type_name)
        if isempty(objects)
            return
        end
        T = typeof(objects[1])
        field_names = fieldnames(T)
        field_types = T.types

        for (name, field_type) in zip(field_names, field_types)
            if field_type <: SimFloat
                push!(all_fields, "$(type_name).$(name)")
            elseif field_type <: AbstractVector{SimFloat}
                first_value = getfield(objects[1], name)
                n_components = length(first_value)
                for component_idx in 1:n_components
                    push!(all_fields, "$(type_name).$(name)[$(component_idx)]")
                end
            elseif field_type <: AbstractMatrix{SimFloat}
                first_value = getfield(objects[1], name)
                n_components = length(first_value)
                for component_idx in 1:n_components
                    push!(all_fields, "$(type_name).$(name)[$(component_idx)]")
                end
            end
        end
    end

    # Collect from all component types
    for (type_name, objects) in [("points", sys.points), ("groups", sys.groups),
                                  ("segments", sys.segments), ("pulleys", sys.pulleys),
                                  ("tethers", sys.tethers), ("winches", sys.winches),
                                  ("wings", sys.wings), ("transforms", sys.transforms)]
        collect_field_names(objects, type_name)
    end

    # System structure itself
    collect_field_names([sys], "system")

    return sort(all_fields)
end

# ==================== METADATA GENERATION ==================== #

"""
    _collect_field_metadata(obj, prefix="")

Collect metadata about all SimFloat fields and vectors in an object.
Returns a vector of dictionaries containing field information for reconstruction.
"""
function _collect_field_metadata(obj, prefix="")
    metadata = Dict{String, Any}[]
    T = typeof(obj)
    field_names = fieldnames(T)
    field_types = T.types

    for (name, field_type) in zip(field_names, field_types)
        value = getfield(obj, name)
        field_path = isempty(prefix) ? string(name) : "$(prefix).$(name)"

        if field_type <: SimFloat
            push!(metadata, Dict(
                "name" => field_path,
                "type" => "SimFloat",
                "size" => 1,
                "start_idx" => nothing  # Will be filled during collection
            ))
        elseif field_type <: AbstractVector{SimFloat}
            push!(metadata, Dict(
                "name" => field_path,
                "type" => "Vector{SimFloat}",
                "size" => length(value),
                "shape" => size(value),
                "start_idx" => nothing  # Will be filled during collection
            ))
        elseif field_type <: AbstractMatrix{SimFloat}
            push!(metadata, Dict(
                "name" => field_path,
                "type" => "Matrix{SimFloat}",
                "size" => length(value),
                "shape" => size(value),
                "start_idx" => nothing  # Will be filled during collection
            ))
        end
    end

    return metadata
end

"""
    _collect_grouped_scalar_metadata(objects, type_name, config::LoggingConfig)

Collect metadata for all scalar fields across all objects of a given type.
Returns grouped metadata where each scalar field gets one entry per component type.
Only includes fields that pass the configuration filter.
"""
function _collect_grouped_scalar_metadata(objects, type_name, config::LoggingConfig)
    metadata = Dict{String, Any}[]
    if isempty(objects)
        return metadata
    end

    # Get the field structure from the first object
    T = typeof(objects[1])
    field_names = fieldnames(T)
    field_types = T.types

    # Collect scalar fields
    for (name, field_type) in zip(field_names, field_types)
        if field_type <: SimFloat
            field_name = "$(type_name).$(name)"
            if should_log_field(field_name, config)
                push!(metadata, Dict(
                    "name" => field_name,
                    "type" => "Vector{SimFloat}",  # Grouped scalars become vectors
                    "size" => length(objects),
                    "shape" => (length(objects),),
                    "start_idx" => nothing  # Will be filled during collection
                ))
            end
        end
    end

    return metadata
end

"""
    _collect_grouped_vector_metadata(objects, type_name, config::LoggingConfig)

Collect metadata for all vector/matrix fields across all objects of a given type.
Returns grouped metadata where each vector component gets one entry per component type.
Only includes fields that pass the configuration filter.
"""
function _collect_grouped_vector_metadata(objects, type_name, config::LoggingConfig)
    metadata = Dict{String, Any}[]
    if isempty(objects)
        return metadata
    end

    # Get the field structure from the first object
    T = typeof(objects[1])
    field_names = fieldnames(T)
    field_types = T.types

    # Collect vector/matrix fields, component by component
    for (name, field_type) in zip(field_names, field_types)
        if field_type <: AbstractVector{SimFloat}
            # Get the size from the first object
            first_value = getfield(objects[1], name)
            n_components = length(first_value)
            # Create metadata for each component
            for component_idx in 1:n_components
                field_name = "$(type_name).$(name)[$(component_idx)]"
                if should_log_field(field_name, config)
                    push!(metadata, Dict(
                        "name" => field_name,
                        "type" => "Vector{SimFloat}",  # Grouped components become vectors
                        "size" => length(objects),
                        "shape" => (length(objects),),
                        "start_idx" => nothing  # Will be filled during collection
                    ))
                end
            end
        elseif field_type <: AbstractMatrix{SimFloat}
            # Get the size from the first object
            first_value = getfield(objects[1], name)
            n_components = length(first_value)
            # Create metadata for each component (flattened)
            for component_idx in 1:n_components
                field_name = "$(type_name).$(name)[$(component_idx)]"
                if should_log_field(field_name, config)
                    push!(metadata, Dict(
                        "name" => field_name,
                        "type" => "Vector{SimFloat}",  # Grouped components become vectors
                        "size" => length(objects),
                        "shape" => (length(objects),),
                        "start_idx" => nothing  # Will be filled during collection
                    ))
                end
            end
        end
    end

    return metadata
end

"""
    _generate_system_metadata(sys::SystemStructure, config::LoggingConfig)

Generate complete metadata for all loggable fields in a SystemStructure using grouped approach.
Only includes fields that pass the configuration filter.
Returns a dictionary with field information and total variable count.
"""
function _generate_system_metadata(sys::SystemStructure, config::LoggingConfig)
    all_metadata = Dict{String, Any}[]
    current_idx = 1

    # Process each component type using grouped approach
    # First all scalars, then all vectors

    # Collect scalar metadata for each component type
    for (type_name, objects) in [("points", sys.points), ("groups", sys.groups),
                                  ("segments", sys.segments), ("pulleys", sys.pulleys),
                                  ("tethers", sys.tethers), ("winches", sys.winches),
                                  ("wings", sys.wings), ("transforms", sys.transforms)]
        scalar_meta = _collect_grouped_scalar_metadata(objects, type_name, config)
        for meta in scalar_meta
            meta["start_idx"] = current_idx
            current_idx += meta["size"]
        end
        append!(all_metadata, scalar_meta)
    end

    # Collect vector metadata for each component type
    for (type_name, objects) in [("points", sys.points), ("groups", sys.groups),
                                  ("segments", sys.segments), ("pulleys", sys.pulleys),
                                  ("tethers", sys.tethers), ("winches", sys.winches),
                                  ("wings", sys.wings), ("transforms", sys.transforms)]
        vector_meta = _collect_grouped_vector_metadata(objects, type_name, config)
        for meta in vector_meta
            meta["start_idx"] = current_idx
            current_idx += meta["size"]
        end
        append!(all_metadata, vector_meta)
    end

    # SystemStructure itself
    sys_scalar_meta = _collect_grouped_scalar_metadata([sys], "system", config)
    for meta in sys_scalar_meta
        meta["start_idx"] = current_idx
        current_idx += meta["size"]
    end
    append!(all_metadata, sys_scalar_meta)

    sys_vector_meta = _collect_grouped_vector_metadata([sys], "system", config)
    for meta in sys_vector_meta
        meta["start_idx"] = current_idx
        current_idx += meta["size"]
    end
    append!(all_metadata, sys_vector_meta)

    return Dict(
        "fields" => all_metadata,
        "total_variables" => current_idx - 1,
        "system_info" => Dict(
            "name" => sys.name,
            "n_points" => length(sys.points),
            "n_groups" => length(sys.groups),
            "n_segments" => length(sys.segments),
            "n_pulleys" => length(sys.pulleys),
            "n_tethers" => length(sys.tethers),
            "n_winches" => length(sys.winches),
            "n_wings" => length(sys.wings),
            "n_transforms" => length(sys.transforms)
        )
    )
end

# ==================== LOGGER INITIALIZATION ==================== #

"""
    create_logger(sys::SystemStructure, capacity::Int=10000;
                  include_fields=nothing, exclude_fields=nothing,
                  include_patterns=nothing, exclude_patterns=nothing,
                  config::Union{LoggingConfig, Nothing}=nothing)

Create a SystemLogger for efficiently recording system states over time.

Automatically generates metadata about all loggable fields in the system and
pre-allocates memory for efficient data collection. Supports selective logging
through field filtering.

# Arguments
- `sys::SystemStructure`: The system to log
- `capacity::Int=10000`: Maximum number of time steps to log

# Keyword Arguments
- `include_fields`: Field names to include (String, Vector{String}, Regex, or Vector{Regex})
- `exclude_fields`: Field names to exclude (String, Vector{String}, Regex, or Vector{Regex})
- `include_patterns`: Regex patterns for fields to include (Regex or Vector{Regex})
- `exclude_patterns`: Regex patterns for fields to exclude (Regex or Vector{Regex})
- `config::LoggingConfig`: Pre-configured logging configuration (overrides other parameters)

If no filtering parameters are provided, all compatible fields are logged.
Regex patterns passed to `*_fields` parameters are automatically moved to `*_patterns`.

# Returns
- `SystemLogger`: Initialized logger ready for data collection

# Examples
```julia
# Log everything (default behavior)
logger = create_logger(sys, 5000)

# Log only position fields (multiple equivalent syntaxes)
logger = create_logger(sys, 5000, include_fields=r"\\.pos_w")
logger = create_logger(sys, 5000, include_fields=[r"\\.pos_w"])
logger = create_logger(sys, 5000, include_patterns=r"\\.pos_w")
logger = create_logger(sys, 5000, include_patterns=[r"\\.pos_w"])

# Log position and velocity fields
logger = create_logger(sys, 5000, include_fields=[r"\\.pos_w", r"\\.vel_w"])
logger = create_logger(sys, 5000, include_patterns=[r"\\.pos_w", r"\\.vel_w"])

# Exclude internal computation fields
logger = create_logger(sys, 5000, exclude_fields=r"\\.internal_")

# Use specific field names (strings)
logger = create_logger(sys, 5000, include_fields=["points.pos_w[1]", "points.mass"])

# Use pre-configured settings
config = LoggingConfig(include_patterns=[r"\\.pos_w"])
logger = create_logger(sys, 5000, config=config)
```
"""
function create_logger(sys::SystemStructure, capacity::Int=10000;
                       include_fields=nothing, exclude_fields=nothing,
                       include_patterns=nothing, exclude_patterns=nothing,
                       config::Union{LoggingConfig, Nothing}=nothing)
    # Use provided config or create one from parameters
    if config === nothing
        # Handle common user mistakes: convert regex to pattern parameters
        if include_fields isa Regex
            include_patterns = [include_fields]
            include_fields = nothing
        elseif include_fields isa Vector && !isempty(include_fields) && include_fields[1] isa Regex
            include_patterns = include_fields
            include_fields = nothing
        end

        if exclude_fields isa Regex
            exclude_patterns = [exclude_fields]
            exclude_fields = nothing
        elseif exclude_fields isa Vector && !isempty(exclude_fields) && exclude_fields[1] isa Regex
            exclude_patterns = exclude_fields
            exclude_fields = nothing
        end

        # Convert single patterns to vectors
        if include_patterns isa Regex
            include_patterns = [include_patterns]
        end
        if exclude_patterns isa Regex
            exclude_patterns = [exclude_patterns]
        end

        config = LoggingConfig(include_fields, exclude_fields, include_patterns, exclude_patterns)
    end

    metadata = _generate_system_metadata(sys, config)
    n_vars = metadata["total_variables"]

    # Provide helpful feedback if no fields were selected
    if n_vars == 0
        all_fields = list_available_fields(sys)
        if !isempty(all_fields)
            @warn """No fields matched the filtering criteria.

                     Available fields include: $(join(first(all_fields, 5), ", "))$(length(all_fields) > 5 ? " ... (+$(length(all_fields)-5) more)" : "")

                     Common regex patterns:
                     - r"\\\\.pos_w" matches position fields
                     - r"\\\\.vel_w" matches velocity fields
                     - r"\\\\.mass" matches mass fields

                     Use list_available_fields(sys) to see all available fields."""
        end
    end

    # Pre-allocate memory
    data = zeros(SimFloat, capacity, n_vars)
    time_stamps = zeros(SimFloat, capacity)

    return SystemLogger(data, time_stamps, metadata, 0, capacity)
end

"""
    resize_logger!(logger::SystemLogger, new_capacity::Int)

Resize a logger's capacity, preserving existing data.
"""
function resize_logger!(logger::SystemLogger, new_capacity::Int)
    old_data = logger.data
    old_times = logger.time_stamps
    n_vars = size(old_data, 2)

    # Create new arrays
    new_data = zeros(SimFloat, new_capacity, n_vars)
    new_time_stamps = zeros(SimFloat, new_capacity)

    # Copy existing data
    copy_steps = min(logger.current_step, new_capacity)
    if copy_steps > 0
        new_data[1:copy_steps, :] .= old_data[1:copy_steps, :]
        new_time_stamps[1:copy_steps] .= old_times[1:copy_steps]
    end

    # Update logger fields
    logger.data .= new_data
    logger.time_stamps .= new_time_stamps
    logger.capacity = new_capacity

    return logger
end

# ==================== LOGGING FUNCTIONS ==================== #

"""
    log!(logger::SystemLogger, sys::SystemStructure, time::SimFloat)

Efficiently log the current state of a SystemStructure.

Uses the existing `sys.vec` property for zero-copy data extraction and stores
the flattened state vector along with the timestamp.

# Arguments
- `logger::SystemLogger`: The logger to record to
- `sys::SystemStructure`: The system to log
- `time::SimFloat`: Current simulation time

# Returns
- `Bool`: `true` if logged successfully, `false` if logger is full

# Example
```julia
logger = create_logger(sys, 1000)
for t in 0.0:0.1:10.0
    # ... simulation step ...
    log!(logger, sys, t)
end
```
"""
function log!(logger::SystemLogger, sys::SystemStructure, time::SimFloat)
    # Check if logger has capacity
    if logger.current_step >= logger.capacity
        # Automatically resize by doubling capacity
        new_capacity = logger.capacity * 2
        @warn "Logger capacity exceeded, resizing from $(logger.capacity) to $new_capacity"
        resize_logger!(logger, new_capacity)
    end

    # Increment step counter
    logger.current_step += 1
    step = logger.current_step

    # Record timestamp
    logger.time_stamps[step] = time

    # Record state vector using filtered collection that matches logger metadata
    state_vec = _collect_logged_fields(sys, logger.metadata["fields"])
    logger.data[step, :] .= state_vec

    return true
end

"""
    log!(logger::SystemLogger, sys::SystemStructure)

Log the current state without a timestamp (uses step number as time).
"""
function log!(logger::SystemLogger, sys::SystemStructure)
    return log!(logger, sys, SimFloat(logger.current_step + 1))
end

"""
    is_full(logger::SystemLogger)

Check if the logger has reached its capacity.
"""
is_full(logger::SystemLogger) = logger.current_step >= logger.capacity

"""
    n_logged(logger::SystemLogger)

Get the number of states currently logged.
"""
n_logged(logger::SystemLogger) = logger.current_step

# ==================== DATA EXTRACTION ==================== #

"""
    get_log_data(logger::SystemLogger; format=:matrix)

Extract logged data in various formats.

# Arguments
- `logger::SystemLogger`: The logger containing the data
- `format::Symbol`: Output format (`:matrix`, `:named_tuple`, `:dict`)

# Returns
- `:matrix`: Returns `(times, data)` where `data` is a matrix [time_steps × variables]
- `:named_tuple`: Returns a named tuple with time and individual field vectors
- `:dict`: Returns a dictionary with field names as keys

# Example
```julia
# Get raw matrix data
times, data = get_log_data(logger, format=:matrix)

# Get structured data with field names
result = get_log_data(logger, format=:named_tuple)
plot(result.time, result.points[1].pos_w[:, 1])  # Plot x-position of first point
```
"""
function get_log_data(logger::SystemLogger; format=:matrix)
    n_steps = logger.current_step
    if n_steps == 0
        if format == :matrix
            return (SimFloat[], zeros(SimFloat, 0, 0))
        elseif format == :named_tuple
            return (time=SimFloat[],)
        elseif format == :dict
            return Dict("time" => SimFloat[])
        end
    end

    times = logger.time_stamps[1:n_steps]
    data = logger.data[1:n_steps, :]

    if format == :matrix
        return (times, data)
    elseif format == :named_tuple || format == :dict
        return _extract_structured_data(logger, times, data, format)
    else
        error("Unknown format: $format. Use :matrix, :named_tuple, or :dict")
    end
end

"""
    _extract_structured_data(logger, times, data, format)

Extract logged data into structured format using metadata.
"""
function _extract_structured_data(logger::SystemLogger, times, data, format)
    fields_meta = logger.metadata["fields"]
    result_data = Dict{String, Any}()
    result_data["time"] = times

    for field_meta in fields_meta
        name = field_meta["name"]
        start_idx = field_meta["start_idx"]
        field_size = field_meta["size"]
        field_type = field_meta["type"]

        # Extract the data slice
        field_data = data[:, start_idx:start_idx+field_size-1]

        if field_type == "SimFloat"
            # Scalar field - extract as vector over time
            result_data[name] = field_data[:, 1]
        elseif field_type == "Vector{SimFloat}"
            # Vector field - keep as matrix [time × vector_components]
            result_data[name] = field_data
        elseif field_type == "Matrix{SimFloat}"
            # Matrix field - reshape back to original dimensions
            original_shape = field_meta["shape"]
            n_steps = size(field_data, 1)
            # Reshape to [time × matrix_rows × matrix_cols]
            # Ensure original_shape is properly handled for splatting
            if original_shape isa Tuple
                reshaped = reshape(field_data, n_steps, original_shape...)
            else
                # Convert to tuple if it's an array or other iterable
                reshaped = reshape(field_data, n_steps, Tuple(original_shape)...)
            end
            result_data[name] = reshaped
        end
    end

    if format == :dict
        return result_data
    elseif format == :named_tuple
        # Convert nested dictionary to nested named tuple
        return _dict_to_nested_named_tuple(result_data)
    end
end

"""
    _dict_to_nested_named_tuple(d::Dict)

Convert a nested dictionary to a nested named tuple with proper structure.
"""
function _dict_to_nested_named_tuple(d::Dict)
    # Group fields by component type and index
    grouped = Dict{String, Dict{String, Any}}()

    for (key, value) in d
        if key == "time"
            grouped["time"] = value
            continue
        end

        # Parse field names like "points.pos_w[1]" or "system.wind_elevation"
        if occursin(".", key)
            parts = split(key, ".")
            comp_type = parts[1]      # "points" or "system"
            field_part = parts[2]     # "pos_w[1]" or "wind_elevation"

            if !haskey(grouped, comp_type)
                grouped[comp_type] = Dict{String, Any}()
            end

            # The field is already grouped, so we can directly store it
            grouped[comp_type][field_part] = value
        else
            # Direct field
            grouped[key] = value
        end
    end

    # Convert to named tuple structure
    result = Dict{Symbol, Any}()

    for (key, value) in grouped
        if key == "time"
            result[:time] = value
        elseif isa(value, Dict)
            # Convert grouped field dictionary to named tuple
            # Fields like "pos_w[1]", "pos_w[2]", "scalar_field" become symbols
            result[Symbol(key)] = NamedTuple{Tuple(Symbol.(keys(value)))}(values(value))
        else
            result[Symbol(key)] = value
        end
    end

    return NamedTuple{Tuple(Symbol.(keys(result)))}(values(result))
end

"""
    get_field_data(logger::SystemLogger, field_name::String)

Extract data for a specific field by name.

# Example
```julia
pos_data = get_field_data(logger, "points[1].pos_w")  # Returns [time_steps × 3] matrix
```
"""
function get_field_data(logger::SystemLogger, field_name::String)
    fields_meta = logger.metadata["fields"]

    for field_meta in fields_meta
        if field_meta["name"] == field_name
            start_idx = field_meta["start_idx"]
            field_size = field_meta["size"]
            n_steps = logger.current_step

            data = logger.data[1:n_steps, start_idx:start_idx+field_size-1]

            if field_meta["type"] == "SimFloat"
                return data[:, 1]  # Return as vector
            elseif field_meta["type"] == "Matrix{SimFloat}"
                # Reshape back to original matrix structure
                original_shape = field_meta["shape"]
                # Ensure original_shape is properly handled for splatting
                if original_shape isa Tuple
                    return reshape(data, n_steps, original_shape...)
                else
                    # Convert to tuple if it's an array or other iterable
                    return reshape(data, n_steps, Tuple(original_shape)...)
                end
            else
                return data  # Return as matrix for vectors
            end
        end
    end

    error("Field '$field_name' not found in logger metadata")
end

"""
    list_fields(logger::SystemLogger)

List all available field names in the logger.
"""
function list_fields(logger::SystemLogger)
    return [field["name"] for field in logger.metadata["fields"]]
end

# ==================== PERSISTENCE ==================== #

"""
    save_log(logger::SystemLogger, filename::String)

Save a SystemLogger to disk for later analysis.

The logger is saved in a Julia-native format that preserves all metadata
and can be loaded back with `load_log()`.

# Arguments
- `logger::SystemLogger`: The logger to save
- `filename::String`: Path to save file (will add .jls extension if not present)

# Example
```julia
save_log(logger, "simulation_results")  # Saves to "simulation_results.jls"
```
"""
function save_log(logger::SystemLogger, filename::String)
    # Add extension if not present
    filepath = endswith(filename, ".jls") ? filename : filename * ".jls"

    # Create a lightweight save structure
    save_data = Dict(
        "data" => logger.data[1:logger.current_step, :],
        "time_stamps" => logger.time_stamps[1:logger.current_step],
        "metadata" => logger.metadata,
        "current_step" => logger.current_step,
        "saved_at" => string(now()),
        "version" => "1.0"
    )

    # Use Julia's built-in serialization
    open(filepath, "w") do io
        serialize(io, save_data)
    end

    @info "Logger saved to $filepath ($(logger.current_step) steps, $(size(logger.data, 2)) variables)"
    return filepath
end

"""
    load_log(filename::String)

Load a previously saved SystemLogger from disk.

# Arguments
- `filename::String`: Path to the saved logger file

# Returns
- `SystemLogger`: The loaded logger

# Example
```julia
logger = load_log("simulation_results.jls")
times, data = get_log_data(logger)
```
"""
function load_log(filename::String)
    # Add extension if not present
    filepath = endswith(filename, ".jls") ? filename : filename * ".jls"

    if !isfile(filepath)
        error("File not found: $filepath")
    end

    # Load the data
    save_data = open(filepath, "r") do io
        deserialize(io)
    end

    # Validate version compatibility
    if !haskey(save_data, "version")
        @warn "Loading legacy logger format - some features may not work correctly"
    end

    # Reconstruct the logger
    n_steps, n_vars = size(save_data["data"])
    capacity = n_steps  # Set capacity to actual data size

    # Create new arrays with the exact size needed
    data = zeros(SimFloat, capacity, n_vars)
    time_stamps = zeros(SimFloat, capacity)

    # Copy the loaded data
    data[1:n_steps, :] .= save_data["data"]
    time_stamps[1:n_steps] .= save_data["time_stamps"]

    logger = SystemLogger(
        data,
        time_stamps,
        save_data["metadata"],
        save_data["current_step"],
        capacity
    )

    @info "Logger loaded from $filepath ($(logger.current_step) steps, $(size(logger.data, 2)) variables)"
    return logger
end

"""
    export_csv(logger::SystemLogger, filename::String; selected_fields=nothing)

Export logged data to CSV format for analysis in external tools.

# Arguments
- `logger::SystemLogger`: The logger containing the data
- `filename::String`: Path to CSV file (will add .csv extension if not present)
- `selected_fields`: Vector of field names to export (exports all if `nothing`)

# Example
```julia
# Export all fields
export_csv(logger, "results")

# Export specific fields
export_csv(logger, "positions", selected_fields=["points[1].pos_w", "points[2].pos_w"])
```
"""
function export_csv(logger::SystemLogger, filename::String; selected_fields=nothing)
    # Add extension if not present
    filepath = endswith(filename, ".csv") ? filename : filename * ".csv"

    # Get the data in dictionary format
    data_dict = get_log_data(logger, format=:dict)
    times = data_dict["time"]

    # Determine which fields to export
    available_fields = list_fields(logger)
    if selected_fields === nothing
        fields_to_export = available_fields
    else
        # Validate that selected fields exist
        missing_fields = setdiff(selected_fields, available_fields)
        if !isempty(missing_fields)
            error("Fields not found: $(missing_fields)")
        end
        fields_to_export = selected_fields
    end

    # Prepare header and data for CSV with vector-valued cells
    headers = ["time"]
    csv_rows = []

    # Initialize rows with time data as strings to allow mixed content
    n_rows = length(times)
    for i in 1:n_rows
        push!(csv_rows, [string(times[i])])
    end

    for field_name in fields_to_export
        field_data = data_dict[field_name]
        push!(headers, field_name)

        # With the new grouped approach, all fields are matrices (vectors over time)
        # Each row contains a vector of values from all objects of that type
        for i in 1:n_rows
            row_vector = field_data[i, :]
            # Format as vector string like "[1.0, 2.0, 3.0]"
            vector_str = "[" * join(row_vector, ", ") * "]"
            push!(csv_rows[i], vector_str)
        end
    end

    # Write to CSV
    open(filepath, "w") do io
        # Write header
        println(io, join(headers, ","))

        # Write data rows
        for row in csv_rows
            println(io, join(row, ","))
        end
    end

    @info "Data exported to $filepath ($(length(headers)) columns, $(n_rows) rows)"
    return filepath
end

# ==================== CONVENIENCE METHODS ==================== #

"""
    log_state!(sys::SystemStructure, logger::SystemLogger, time::SimFloat)

Convenience method to log the current state. Same as `log!(logger, sys, time)`.
"""
log_state!(sys::SystemStructure, logger::SystemLogger, time::SimFloat) = log!(logger, sys, time)

"""
    log_state!(sys::SystemStructure, logger::SystemLogger)

Convenience method to log the current state without timestamp.
"""
log_state!(sys::SystemStructure, logger::SystemLogger) = log!(logger, sys)

"""
    Base.summary(logger::SystemLogger)

Provide a summary of the logger's current state.
"""
function Base.summary(logger::SystemLogger)
    return "SystemLogger($(logger.current_step)/$(logger.capacity) steps, $(size(logger.data, 2)) variables)"
end

"""
    Base.show(io::IO, logger::SystemLogger)

Custom display for SystemLogger objects.
"""
function Base.show(io::IO, logger::SystemLogger)
    n_vars = size(logger.data, 2)
    n_logged = logger.current_step
    capacity = logger.capacity
    usage_pct = round(100 * n_logged / capacity, digits=1)

    print(io, "SystemLogger:\n")
    print(io, "  Steps: $n_logged/$capacity ($(usage_pct)%)\n")
    print(io, "  Variables: $n_vars\n")

    if n_logged > 0
        time_range = (logger.time_stamps[1], logger.time_stamps[n_logged])
        print(io, "  Time range: $(time_range[1]) → $(time_range[2])\n")
    end

    # Show system info if available
    if haskey(logger.metadata, "system_info")
        sys_info = logger.metadata["system_info"]
        print(io, "  System: $(sys_info["name"]) ($(sys_info["n_points"]) points, $(sys_info["n_wings"]) wings)\n")
    end

    # Show a few example field names
    fields = list_fields(logger)
    if length(fields) > 0
        n_show = min(3, length(fields))
        sample_fields = fields[1:n_show]
        if length(fields) > n_show
            print(io, "  Sample fields: $(join(sample_fields, ", ")) ... (+$(length(fields) - n_show) more)")
        else
            print(io, "  Fields: $(join(sample_fields, ", "))")
        end
    end
end

"""
    callback_logger(logger::SystemLogger)

Create a callback function for automatic logging during simulation.

Returns a function that can be used as a callback in ODE solvers or
custom simulation loops to automatically log system states.

# Example
```julia
logger = create_logger(sys, 10000)
callback = callback_logger(logger)

# In simulation loop:
for (i, t) in enumerate(time_steps)
    # ... simulation step ...
    callback(sys, t, i)  # Automatically logs every step
end
```
"""
function callback_logger(logger::SystemLogger)
    return function(sys::SystemStructure, time::SimFloat, step::Int=0)
        log!(logger, sys, time)
    end
end

"""
    selective_logger(logger::SystemLogger, condition_fn)

Create a selective logging callback that only logs when a condition is met.

# Arguments
- `logger::SystemLogger`: The logger to use
- `condition_fn`: Function that takes `(sys, time, step)` and returns `Bool`

# Example
```julia
# Log every 10th step
selective_cb = selective_logger(logger, (sys, t, step) -> step % 10 == 0)

# Log when wing elevation changes significantly
last_elevation = Ref(0.0)
elevation_cb = selective_logger(logger, (sys, t, step) -> begin
    current_elev = sys.wings[1].elevation
    if abs(current_elev - last_elevation[]) > 0.1
        last_elevation[] = current_elev
        return true
    end
    return false
end)
```
"""
function selective_logger(logger::SystemLogger, condition_fn)
    return function(sys::SystemStructure, time::SimFloat, step::Int=0)
        if condition_fn(sys, time, step)
            log!(logger, sys, time)
        end
    end
end
