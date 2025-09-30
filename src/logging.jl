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
    _generate_system_metadata(sys::SystemStructure)

Generate complete metadata for all loggable fields in a SystemStructure.
Returns a dictionary with field information and total variable count.
"""
function _generate_system_metadata(sys::SystemStructure)
    all_metadata = Dict{String, Any}[]
    current_idx = 1

    # Process each component type
    for (i, point) in enumerate(sys.points)
        point_meta = _collect_field_metadata(point, "points[$i]")
        for meta in point_meta
            meta["start_idx"] = current_idx
            current_idx += meta["size"]
        end
        append!(all_metadata, point_meta)
    end

    for (i, group) in enumerate(sys.groups)
        group_meta = _collect_field_metadata(group, "groups[$i]")
        for meta in group_meta
            meta["start_idx"] = current_idx
            current_idx += meta["size"]
        end
        append!(all_metadata, group_meta)
    end

    for (i, segment) in enumerate(sys.segments)
        segment_meta = _collect_field_metadata(segment, "segments[$i]")
        for meta in segment_meta
            meta["start_idx"] = current_idx
            current_idx += meta["size"]
        end
        append!(all_metadata, segment_meta)
    end

    for (i, pulley) in enumerate(sys.pulleys)
        pulley_meta = _collect_field_metadata(pulley, "pulleys[$i]")
        for meta in pulley_meta
            meta["start_idx"] = current_idx
            current_idx += meta["size"]
        end
        append!(all_metadata, pulley_meta)
    end

    for (i, tether) in enumerate(sys.tethers)
        tether_meta = _collect_field_metadata(tether, "tethers[$i]")
        for meta in tether_meta
            meta["start_idx"] = current_idx
            current_idx += meta["size"]
        end
        append!(all_metadata, tether_meta)
    end

    for (i, winch) in enumerate(sys.winches)
        winch_meta = _collect_field_metadata(winch, "winches[$i]")
        for meta in winch_meta
            meta["start_idx"] = current_idx
            current_idx += meta["size"]
        end
        append!(all_metadata, winch_meta)
    end

    for (i, wing) in enumerate(sys.wings)
        wing_meta = _collect_field_metadata(wing, "wings[$i]")
        for meta in wing_meta
            meta["start_idx"] = current_idx
            current_idx += meta["size"]
        end
        append!(all_metadata, wing_meta)
    end

    for (i, transform) in enumerate(sys.transforms)
        transform_meta = _collect_field_metadata(transform, "transforms[$i]")
        for meta in transform_meta
            meta["start_idx"] = current_idx
            current_idx += meta["size"]
        end
        append!(all_metadata, transform_meta)
    end

    # SystemStructure itself
    sys_meta = _collect_field_metadata(sys, "system")
    for meta in sys_meta
        meta["start_idx"] = current_idx
        current_idx += meta["size"]
    end
    append!(all_metadata, sys_meta)

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
    create_logger(sys::SystemStructure, capacity::Int=10000)

Create a SystemLogger for efficiently recording system states over time.

Automatically generates metadata about all loggable fields in the system and
pre-allocates memory for efficient data collection.

# Arguments
- `sys::SystemStructure`: The system to log
- `capacity::Int=10000`: Maximum number of time steps to log

# Returns
- `SystemLogger`: Initialized logger ready for data collection
"""
function create_logger(sys::SystemStructure, capacity::Int=10000)
    metadata = _generate_system_metadata(sys)
    n_vars = metadata["total_variables"]

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

    # Record state vector using the existing vec property (zero-copy)
    state_vec = vec(sys.vec)  # Flatten the 2D column vector to 1D
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

        # Parse field names like "points[1].pos_w" or "system.wind_elevation"
        if occursin("[", key)
            # Component with index: "points[1].pos_w"
            parts = split(key, ".")
            component_part = parts[1]  # "points[1]"
            field_name = parts[2]      # "pos_w"

            # Extract component type and index
            comp_match = match(r"(\w+)\[(\d+)\]", component_part)
            if comp_match !== nothing
                comp_type = comp_match.captures[1]
                comp_idx = parse(Int, comp_match.captures[2])

                if !haskey(grouped, comp_type)
                    grouped[comp_type] = Dict{String, Any}()
                end
                if !haskey(grouped[comp_type], string(comp_idx))
                    grouped[comp_type][string(comp_idx)] = Dict{String, Any}()
                end

                grouped[comp_type][string(comp_idx)][field_name] = value
            end
        else
            # System-level field: "system.wind_elevation"
            parts = split(key, ".")
            if length(parts) == 2 && parts[1] == "system"
                if !haskey(grouped, "system")
                    grouped["system"] = Dict{String, Any}()
                end
                grouped["system"][parts[2]] = value
            else
                # Direct field
                grouped[key] = value
            end
        end
    end

    # Convert to named tuple structure
    result = Dict{Symbol, Any}()

    for (key, value) in grouped
        if key == "time"
            result[:time] = value
        elseif isa(value, Dict) && any(k -> tryparse(Int, k) !== nothing, keys(value))
            # Component array (points, groups, etc.)
            indices = sort([parse(Int, k) for k in keys(value) if tryparse(Int, k) !== nothing])
            component_array = []
            for idx in indices
                idx_str = string(idx)
                if haskey(value, idx_str)
                    component_data = value[idx_str]
                    # Convert to named tuple
                    push!(component_array, NamedTuple{Tuple(Symbol.(keys(component_data)))}(values(component_data)))
                end
            end
            result[Symbol(key)] = component_array
        elseif isa(value, Dict)
            # System or other grouped data
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

    # Prepare header and data matrix
    headers = ["time"]
    data_matrix = reshape(times, :, 1)

    for field_name in fields_to_export
        field_data = data_dict[field_name]

        if field_data isa Vector
            # Scalar field over time
            push!(headers, field_name)
            data_matrix = hcat(data_matrix, field_data)
        elseif field_data isa Matrix
            # Vector field over time
            n_components = size(field_data, 2)
            for i in 1:n_components
                push!(headers, "$(field_name)[$i]")
                data_matrix = hcat(data_matrix, field_data[:, i])
            end
        end
    end

    # Write to CSV
    open(filepath, "w") do io
        # Write header
        println(io, join(headers, ","))

        # Write data rows
        for i in 1:size(data_matrix, 1)
            println(io, join(data_matrix[i, :], ","))
        end
    end

    @info "Data exported to $filepath ($(length(headers)) columns, $(size(data_matrix, 1)) rows)"
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
