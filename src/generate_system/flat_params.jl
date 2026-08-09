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
returns its live value, a `kind` (`:scalar`, `:array`, or `:callable`), and the
`path` it was read from (a `(name,)` tuple for computed leaves). The `path` lets a
backend resolve a runtime address per instance (see `ParamView`).
"""
struct ParamEntry
    param::Any
    read::Any
    kind::Symbol
    path::Any
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
    register_leaf!(reg, key, name, reader, value, path)

Create (once, memoised on `key`) and record the flat parameter for a leaf
`value` under symbol `name`. Numeric scalars/arrays become data params; any other
(callable) leaf — an interpolation or polar — becomes a callable param applied as
`name(x)`. `reader` reads the live value from a `sys_struct` at sync time; `path`
is stored on the entry for per-instance address resolution.
"""
function register_leaf!(reg::ParamRegistry, key, name::Symbol, reader, value, path)
    cached = get(reg.cache, key, nothing)
    cached === nothing || return cached
    if value isa Real
        param, kind = make_param(name, value), :scalar
    elseif value isa AbstractArray{<:Real}
        param, kind = make_array_param(name, value), :array
    else
        param, kind = make_callable_param(name, value), :callable
    end
    push!(reg.entries, ParamEntry(param, reader, kind, path))
    reg.cache[key] = param
    return param
end

"""Last `Symbol` in `path` (the struct field name), ignoring trailing indices."""
function leaf_symbol(path::Tuple)
    for key in Iterators.reverse(path)
        key isa Symbol && return key
    end
    error("param path $path has no symbol leaf")
end

"""
Symbolic name for a leaf `path` under backend `B`. The monolith bakes the full
per-instance path into the name (distinct symbol per instance); the network uses
the bare field name (one generic symbol per component *type*).
"""
param_symbol_name(::Type{<:ModelBackend}, path::Tuple) = param_name(path)
param_symbol_name(::Type{NetworkBackend}, path::Tuple) = leaf_symbol(path)

"""
Memoisation key for a leaf `path` under backend `B`: the full path for the
monolith, the bare field name for the network (so one generic symbol is reused
across a kernel's equations).
"""
param_cache_key(::Type{<:ModelBackend}, path::Tuple) = path
param_cache_key(::Type{NetworkBackend}, path::Tuple) = leaf_symbol(path)

"""
    leaf_param!(B, reg, path, reader, value)

Record the flat parameter for a leaf read at `path`, naming and memoising it per
the backend `B` policy ([`param_symbol_name`](@ref), [`param_cache_key`](@ref)).
"""
function leaf_param!(::Type{B}, reg::ParamRegistry, path::Tuple, reader, value,
                     suffix::Symbol = Symbol("")) where {B}
    key = param_cache_key(B, path)
    name = param_symbol_name(B, path)
    if suffix !== Symbol("")
        key = Symbol(key, suffix)
        name = Symbol(name, suffix)
    end
    return register_leaf!(reg, key, name, reader, value, path)
end

"""
    param_computed!(reg, name, reader)

Escape hatch for a value that is not a plain field read — `reader(sys_struct)`
computes it (e.g. a [`WindFactorReader`](@ref) building a callable wind-factor
from the atmospheric model). `reader` must be a named struct (serialisable), not
a closure over `sys_struct`.
"""
param_computed!(reg::ParamRegistry, name::Symbol, reader) =
    register_leaf!(reg, name, name, reader, reader(reg.sys_struct), (name,))

# ==================== BUILD-TIME VIEW ==================== #

"""Types the view descends *through* (everything else is a leaf)."""
param_descend(x) = x isa NamedCollection || x isa AbstractAeroModel ||
                   x isa AbstractWinchModel || x isa VSMEngine || x isa Settings

param_name(path::Tuple) = Symbol("p_", join(path, "_"))

"""
Top-level `params` view wrapping a [`ParamRegistry`](@ref), tagged with the
[`ModelBackend`](@ref) type `B` so leaf resolution dispatches on it.
`params.segments[i].l0` mirrors `sys_struct.segments[i].l0` (build-time only). The
default `ParamView(reg)` is a [`MonolithBackend`](@ref) view.
"""
struct ParamView{B<:ModelBackend}
    reg::ParamRegistry
    suffix::Symbol
end
ParamView{B}(reg::ParamRegistry) where {B<:ModelBackend} =
    ParamView{B}(reg, Symbol(""))
ParamView(reg::ParamRegistry) = ParamView{MonolithBackend}(reg, Symbol(""))

"""
    suffixed(view::ParamView{B}, suffix) -> ParamView{B}

A view sharing `view`'s registry but tagging every leaf it resolves with `suffix`, so
otherwise-aliasing bare field names (the `NetworkBackend` naming) become distinct params.
Used to disambiguate repeated reads of the same field across members of one combined edge.
"""
suffixed(view::ParamView{B}, suffix::Symbol) where {B} =
    ParamView{B}(getfield(view, :reg), suffix)

"""A partial path into `sys_struct` being resolved to a parameter under backend `B`."""
struct PathView{B<:ModelBackend}
    reg::ParamRegistry
    path::Tuple
    suffix::Symbol
end

Base.getproperty(view::ParamView{B}, sym::Symbol) where {B} =
    sym === :reg ? getfield(view, :reg) :
    sym === :suffix ? getfield(view, :suffix) :
    PathView{B}(getfield(view, :reg), (sym,), getfield(view, :suffix))

Base.getindex(view::PathView{B}, idx::Integer) where {B} =
    PathView{B}(getfield(view, :reg), (getfield(view, :path)..., Int(idx)),
                getfield(view, :suffix))

"""
    param_unknowns(params)

The symbolic parameters a `params` view recorded so far, in insertion order — passed
as the parameter list of an `@named` component `System` so every `params.…` read it
made is declared. Used by both backends' component assembly.
"""
param_unknowns(params::ParamView) = Any[entry.param for entry in params.reg.entries]

function Base.getproperty(view::PathView{B}, sym::Symbol) where {B}
    (sym === :reg || sym === :path || sym === :suffix) && return getfield(view, sym)
    reg = getfield(view, :reg)
    suffix = getfield(view, :suffix)
    path = (getfield(view, :path)..., sym)
    value = read_path(reg.sys_struct, path)
    param_descend(value) && return PathView{B}(reg, path, suffix)
    return leaf_param!(B, reg, path, PathReader(path), value, suffix)
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
