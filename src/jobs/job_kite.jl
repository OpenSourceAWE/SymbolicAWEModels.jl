# Copyright (c) 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# A kite as the job runner serves it, its parameter menus, and the job key
# (`api/openapi.yaml`).

"""
    JobParameter(name, unit, min, max, default; step=nothing)

One entry of a kite's job menu: a parameter a person may set, in `unit`, from `min`
to `max`, taking `default` when left out. `step` is the increment a request form
offers.
"""
struct JobParameter
    name::String
    unit::String
    min::Float64
    max::Float64
    default::Float64
    step::Union{Nothing, Float64}
    function JobParameter(name, unit, min, max, default; step=nothing)
        min <= default <= max ||
            throw(ArgumentError("$name: default $default outside [$min, $max]"))
        return new(name, unit, min, max, default, step)
    end
end

"""
    JobKite(name; steady_menu, sim_menu, sweep_menu, max_axes=1,
            steady, simulate, sweep_point)

A kite that [`serve_jobs`](@ref) runs jobs for: a menu of [`JobParameter`](@ref)s per
job kind, and one function per kind. Each takes the job's parameters as a
`Dict{String, Float64}` holding every parameter of its menu, and returns a `SysLog`
whose metadata carries the structure document (see [`sys_log`](@ref)):

- `steady(params)`: a one-frame log of the quasi-steady state;
- `simulate(params; progress)`: the log of a simulation, calling
  `progress(sim_time)` as it goes; `progress` throws when the job is cancelled;
- `sweep_point(params)`: the log of one operating point of a sweep, which varies at
  most `max_axes` of the `sweep_menu` parameters at once.
"""
struct JobKite
    name::String
    steady_menu::Vector{JobParameter}
    sim_menu::Vector{JobParameter}
    sweep_menu::Vector{JobParameter}
    max_axes::Int
    steady::Function
    simulate::Function
    sweep_point::Function
    function JobKite(name; steady_menu, sim_menu, sweep_menu, max_axes=1,
                     steady, simulate, sweep_point)
        1 <= max_axes <= 4 ||
            throw(ArgumentError("$name: max_axes $max_axes outside [1, 4]"))
        return new(name, steady_menu, sim_menu, sweep_menu, max_axes,
                   steady, simulate, sweep_point)
    end
end

"""A job the runner refuses before its kite runs."""
struct JobRejected <: Exception
    msg::String
end

Base.showerror(io::IO, err::JobRejected) = print(io, "job rejected: ", err.msg)

"""The menu of `kite` as the heartbeat sends it: its parameters per job kind."""
function menu_document(kite::JobKite)
    return Dict("steady" => kind_menu(kite.steady_menu),
                "sim" => kind_menu(kite.sim_menu),
                "sweep" => merge(kind_menu(kite.sweep_menu),
                                 Dict("max_axes" => kite.max_axes)))
end

kind_menu(menu::Vector{JobParameter}) = Dict("parameters" => menu_entry.(menu))

function menu_entry(parameter::JobParameter)
    entry = Dict{String, Any}("name" => parameter.name, "unit" => parameter.unit,
                              "min" => parameter.min, "max" => parameter.max,
                              "default" => parameter.default)
    isnothing(parameter.step) || (entry["step"] = parameter.step)
    return entry
end

"""
Every parameter of `menu` but the `swept` ones, from `requested` or else its default.
Throws [`JobRejected`](@ref) for a parameter not on the menu, swept as well, or out
of range.
"""
function menu_values(menu::Vector{JobParameter}, requested::AbstractDict;
                     swept=())
    names = [parameter.name for parameter in menu]
    for name in keys(requested)
        name in names || throw(JobRejected("$name is not on the menu"))
        name in swept && throw(JobRejected("$name is both fixed and swept"))
    end
    chosen = Dict{String, Float64}()
    for parameter in menu
        parameter.name in swept && continue
        value = get(requested, parameter.name, parameter.default)
        chosen[parameter.name] = checked_value(parameter, value)
    end
    return chosen
end

function checked_value(parameter::JobParameter, value)
    (value isa Real && !(value isa Bool) && parameter.min <= value <= parameter.max) ||
        throw(JobRejected("$(parameter.name) = $value outside " *
                          "[$(parameter.min), $(parameter.max)] $(parameter.unit)"))
    return Float64(value)
end

"""
The axes of a sweep job, in the requested order, each value checked against the
kite's sweep menu.
"""
function sweep_axes(kite::JobKite, requested::AbstractDict)
    isempty(requested) && throw(JobRejected("a sweep needs at least one axis"))
    length(requested) <= kite.max_axes ||
        throw(JobRejected("$(length(requested)) axes, $(kite.name) sweeps at most " *
                          "$(kite.max_axes)"))
    axes = OrderedDict{String, Vector{Float64}}()
    for (name, axis_values) in requested
        index = findfirst(parameter -> parameter.name == name, kite.sweep_menu)
        isnothing(index) && throw(JobRejected("$name is not on the sweep menu"))
        (axis_values isa AbstractVector && !isempty(axis_values)) ||
            throw(JobRejected("axis $name has no values"))
        axes[name] = [checked_value(kite.sweep_menu[index], value)
                      for value in axis_values]
    end
    return axes
end

"""
    job_key(sandbox, kite, kind, parameters, version; axes=nothing)

The key the site recognises a duplicate job by: the hex SHA-256 of the RFC 8785
canonical JSON of the job's sandbox, kite name, kind, parameters (every menu
parameter, swept ones in `axes` instead), and the runner `version`.
"""
function job_key(sandbox, kite, kind, parameters, version; axes=nothing)
    job = Dict("sandbox" => sandbox, "kite" => kite, "kind" => kind,
               "parameters" => parameters, "version" => version)
    isnothing(axes) || (job["axes"] = axes)
    return bytes2hex(sha256(canonical_json(job)))
end

"""`value` as RFC 8785 (JCS) canonical JSON: sorted keys, no whitespace."""
function canonical_json(value::AbstractDict)
    members = [JSON.json(string(key)) * ":" * canonical_json(value[key])
               for key in sort!(collect(keys(value)); by=string)]
    return "{" * join(members, ",") * "}"
end
canonical_json(value::AbstractVector) = "[" * join(canonical_json.(value), ",") * "]"
canonical_json(value::AbstractString) = JSON.json(value)
canonical_json(value::Bool) = string(value)
canonical_json(::Nothing) = "null"
canonical_json(value::Real) = canonical_number(Float64(value))

"""A number as ECMAScript prints it, which RFC 8785 prescribes."""
function canonical_number(value::Float64)
    isfinite(value) || throw(ArgumentError("JSON has no number $value"))
    iszero(value) && return "0"
    text = string(abs(value))
    mantissa, exponent = occursin('e', text) ? split(text, 'e') : (text, "0")
    whole, fraction = split(mantissa, '.')
    digits = whole * fraction
    leading_zeros = findfirst(!=('0'), digits) - 1
    digits = rstrip(digits[leading_zeros+1:end], '0')
    k = length(digits)
    n = length(whole) - leading_zeros + parse(Int, exponent)  # value = 0.digits × 10^n
    body = if k <= n <= 21
        digits * "0"^(n - k)
    elseif 0 < n <= 21
        digits[1:n] * "." * digits[n+1:end]
    elseif -6 < n <= 0
        "0." * "0"^(-n) * digits
    else
        significand = k == 1 ? digits : digits[1] * "." * digits[2:end]
        significand * "e" * (n > 0 ? "+" : "-") * string(abs(n - 1))
    end
    return value < 0 ? "-" * body : body
end
