# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Flat MTK parameters replacing per-step `@register_symbolic` struct reads.

"""
    read_path(obj, path)

Walk `path` (a tuple of `Symbol` fields and `Int` indices) from `obj`, e.g.
`read_path(sys_struct, (:wings, 1, :aero, :engine, :aero_jac))`.
"""
read_path(obj, ::Tuple{}) = obj
read_path(obj, path::Tuple) = read_path(path_step(obj, first(path)), Base.tail(path))
path_step(obj, key::Symbol) = getproperty(obj, key)
path_step(obj, key::Int) = obj[key]

"""
    PathReader(path)

Serialisable closure reading a fixed `path` from a `sys_struct` at sync time
(holds only the path tuple, never the struct).
"""
struct PathReader{P<:Tuple}
    path::P
end
(reader::PathReader)(sys_struct) = read_path(sys_struct, reader.path)

"""
    ParamEntry

One flattened parameter: the symbolic `param`, a `read(sys_struct)` callable that
returns its live value, and a `kind` (`:scalar`, `:array`, or `:callable`).
"""
struct ParamEntry
    param::Any
    read::Any
    kind::Symbol
end

"""
    ParamRegistry

Single source of truth for the flattened parameters. Built during equation
generation; the `params` view records one [`ParamEntry`](@ref) per distinct field
it reads (memoised by cache key). After compilation it drives `sync_params!`.
"""
mutable struct ParamRegistry
    sys_struct::SystemStructure
    entries::Vector{ParamEntry}
    cache::Dict{Any, Any}
end
ParamRegistry(sys_struct::SystemStructure) =
    ParamRegistry(sys_struct, ParamEntry[], Dict{Any, Any}())

# ---- parameter constructors (value baked in as default) ----

"""Numeric scalar parameter named `name` (default `value`)."""
make_param(name::Symbol, value::Real) = only(@parameters $name = value)

"""Numeric vector parameter `name[1:n]` (default `value`)."""
make_array_param(name::Symbol, value::AbstractVector) =
    only(@parameters $name[1:length(value)] = collect(value))

"""Numeric matrix parameter `name[1:m,1:n]` (default `value`)."""
make_array_param(name::Symbol, value::AbstractMatrix) =
    only(@parameters $name[1:size(value, 1), 1:size(value, 2)] = collect(value))

"""
Callable parameter `name` (invoked symbolically as `name(x)`; default `value`).
For a leaf that is a function/interpolation/polar — MTK codegens the call and
ForwardDiff differentiates through it, so no `@register_symbolic` is needed.
"""
function make_callable_param(name::Symbol, value)
    T = typeof(value)
    return only(@parameters ($name::T)(..) = value)
end

"""
    leaf_param!(reg, key, name, reader, value)

Create (once, memoised on `key`) and record the flat parameter for a leaf
`value`. Numeric scalars/arrays become data params; any other (callable) leaf —
an interpolation or polar — becomes a callable param applied as `name(x)`.
`reader` reads the live value from a `sys_struct` at sync time.
"""
function leaf_param!(reg::ParamRegistry, key, name::Symbol, reader, value)
    cached = get(reg.cache, key, nothing)
    cached === nothing || return cached
    if value isa Real
        param, kind = make_param(name, value), :scalar
    elseif value isa AbstractArray{<:Real}
        param, kind = make_array_param(name, value), :array
    else
        param, kind = make_callable_param(name, value), :callable
    end
    push!(reg.entries, ParamEntry(param, reader, kind))
    reg.cache[key] = param
    return param
end

"""
    param_computed!(reg, name, reader)

Escape hatch for a value that is not a plain field read — `reader(sys_struct)`
computes it (e.g. a [`WindFactorReader`](@ref) building a callable wind-factor
from the atmospheric model). `reader` must be a named struct (serialisable), not
a closure over `sys_struct`.
"""
param_computed!(reg::ParamRegistry, name::Symbol, reader) =
    leaf_param!(reg, name, name, reader, reader(reg.sys_struct))

# ==================== BUILD-TIME VIEW ==================== #

"""Types the view descends *through* (everything else is a leaf)."""
param_descend(x) = x isa NamedCollection || x isa AbstractAeroModel ||
                   x isa AbstractWinchModel || x isa VSMEngine || x isa Settings

param_name(path::Tuple) = Symbol("p_", join(path, "_"))

"""
Top-level `params` view wrapping a [`ParamRegistry`](@ref).
`params.segments[i].l0` mirrors `sys_struct.segments[i].l0` (build-time only).
"""
struct ParamView
    reg::ParamRegistry
end

"""A partial path into `sys_struct` being resolved to a parameter."""
struct PathView
    reg::ParamRegistry
    path::Tuple
end

Base.getproperty(view::ParamView, sym::Symbol) =
    sym === :reg ? getfield(view, :reg) : PathView(getfield(view, :reg), (sym,))

Base.getindex(view::PathView, idx::Integer) =
    PathView(getfield(view, :reg), (getfield(view, :path)..., Int(idx)))

function Base.getproperty(view::PathView, sym::Symbol)
    (sym === :reg || sym === :path) && return getfield(view, sym)
    reg = getfield(view, :reg)
    path = (getfield(view, :path)..., sym)
    value = read_path(reg.sys_struct, path)
    param_descend(value) && return PathView(reg, path)
    return leaf_param!(reg, path, param_name(path), PathReader(path), value)
end

# ==================== SYNC ==================== #

"""
    ParamGroup

A `setp` setter plus the readers and preallocated value buffer for one parameter
kind. `eltype` is `SimFloat` for numeric scalars, `Any` for arrays/callables.
"""
struct ParamGroup{Setter, Buf}
    setter::Setter
    readers::Vector{Any}
    buffer::Buf
end

"""Bundle of the per-kind sync groups (each may be `nothing`)."""
struct ParamSync{S, A, C}
    scalar::S
    array::A
    callable::C
end

"""
    survivor_index(sys) -> Dict{String, param}

Map each parameter surviving `mtkcompile` to its name, keyed by both the full name
and the leaf name (after the last `₊` namespace separator) so a registry's bare
param matches its namespaced counterpart from a subsystem.
"""
function survivor_index(sys)
    index = Dict{String, Any}()
    for p in parameters(sys)
        name = string(ModelingToolkit.getname(ModelingToolkit.unwrap(p)))
        index[name] = p
        sep = findlast('₊', name)
        sep === nothing || (index[name[nextind(name, sep):end]] = p)
    end
    return index
end

entry_name(param) = string(ModelingToolkit.getname(ModelingToolkit.unwrap(param)))

"""
    build_param_sync(sys, registry) -> ParamSync | Nothing

Build the per-kind sync groups from the compiled system and the registry. Pruned
parameters (no surviving equation references them) are dropped, so a setter never
touches a parameter absent from the buffer.
"""
function build_param_sync(sys, registry::ParamRegistry)
    isempty(registry.entries) && return nothing
    index = survivor_index(sys)
    by_kind(k) = filter(e -> e.kind === k, registry.entries)
    grp(entries, ::Type{Buf}) where {Buf} = begin
        survivors, readers = Any[], Any[]
        for entry in entries
            survivor = get(index, entry_name(entry.param), nothing)
            survivor === nothing && continue
            push!(survivors, survivor); push!(readers, entry.read)
        end
        isempty(survivors) ? nothing :
            ParamGroup(setp(sys, survivors), readers, Buf(undef, length(survivors)))
    end
    scalar = grp(by_kind(:scalar), Vector{SimFloat})
    array = grp(by_kind(:array), Vector{Any})
    callable = grp(by_kind(:callable), Vector{Any})
    (scalar === nothing && array === nothing && callable === nothing) && return nothing
    return ParamSync(scalar, array, callable)
end

"""
    sync_params!(sync, target, sys_struct)

Copy every flattened field from the live `sys_struct` into `target`'s parameter
buffers (`target` is an `ODEProblem` or an `ODEIntegrator`). A no-op when there
are no flattened parameters.
"""
sync_params!(::Nothing, target, sys_struct) = nothing
function sync_params!(sync::ParamSync, target, sys_struct::SystemStructure)
    sync_group!(sync.scalar, target, sys_struct)
    sync_group!(sync.array, target, sys_struct)
    sync_group!(sync.callable, target, sys_struct)
    return nothing
end

sync_group!(::Nothing, target, sys_struct) = nothing
function sync_group!(group::ParamGroup, target, sys_struct::SystemStructure)
    readers = group.readers
    buffer = group.buffer
    @inbounds for k in eachindex(readers)
        buffer[k] = readers[k](sys_struct)
    end
    group.setter(target, buffer)
    return nothing
end
