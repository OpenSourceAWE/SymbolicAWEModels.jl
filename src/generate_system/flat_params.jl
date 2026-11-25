# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Flat MTK parameters holding the values of a run (`run_value`).

"""
    read_path(obj, path)

Walk `path` (a tuple of `Symbol` fields and `Int` indices) from `obj`, e.g.
`read_path(sam, (:wings, 1, :aero, :aero_jac))`. From a run, a leading
`(collection, index, quantity)` is that item's [`run_value`](@ref).
"""
read_path(obj, ::Tuple{}) = obj
read_path(obj, path::Tuple) = read_path(path_step(obj, first(path)), Base.tail(path))
function read_path(sam::SymbolicAWEModel, path::Tuple{Any, Vararg})
    is_item_path(path) || return read_path(path_step(sam, first(path)), Base.tail(path))
    return read_path(run_value(sam, Val(path[1]), path[2], Val(path[3])), path[4:end])
end
@inline path_step(obj, key::Symbol) = getproperty(obj, key)
path_step(obj, key::Int) = obj[key]

"""Whether `path` starts at a quantity of an item of a run: `(collection, index,
quantity, …)`."""
is_item_path(path) = length(path) >= 3 && is_collection(path[1]) && path[2] isa Int &&
    path[3] isa Symbol

"""
    canonical_path(sam, path)

`path` as the parameter it reads is named: a quantity of wing `w` of `sam` that is its
body's, `(:wings, w, quantity)`, reads `(:bodies, wing_body(sam, w), quantity)`.
"""
function canonical_path(sam, path)
    length(path) == 3 && path[1] === :wings && path[2] isa Int &&
        !has_quantity(:wings, path[3]) || return path
    return (:bodies, wing_body(sam, path[2]), path[3])
end

"""
    PathReader(path)

Serialisable closure reading a fixed `path` from a run `sam` at sync time (holds only
the path, never the run).

The field names live in the type and the container indices stay values. A name has
to be a compile-time constant for the field read to resolve to an offset —
otherwise every step of the walk is a dynamic call that boxes both its argument
and its result — while keeping the indices out of the type lets all components of
one kind share a single compiled reader instead of one per index. A leading
`(collection, index, quantity)` reads that item's [`run_value`](@ref).
"""
struct PathReader{Fields, I <: Tuple}
    indices::I
end

function PathReader(path::Tuple)
    fields = map(step -> step isa Symbol ? step : nothing, path)
    indices = Tuple(step for step in path if !(step isa Symbol))
    return PathReader{fields, typeof(indices)}(indices)
end

"""The `(name, index, …)` tuple `reader` was built from."""
function path_of(reader::PathReader{Fields}) where {Fields}
    indices = getfield(reader, :indices)
    taken = 0
    steps = Vector{Any}(undef, length(Fields))
    for (k, field) in enumerate(Fields)
        steps[k] = field === nothing ? indices[taken += 1] : field
    end
    return Tuple(steps)
end

Base.getproperty(reader::PathReader, name::Symbol) =
    name === :path ? path_of(reader) : getfield(reader, name)

@generated function (reader::PathReader{Fields})(sam) where {Fields}
    walk, taken, start = :sam, 0, 1
    if length(Fields) >= 3 && Fields[1] isa Symbol &&
            hasfield(typeof(QUANTITIES), Fields[1]) && Fields[2] === nothing &&
            Fields[3] isa Symbol
        walk = :(run_value(sam, Val($(QuoteNode(Fields[1]))),
                           getfield(reader, :indices)[1], Val($(QuoteNode(Fields[3])))))
        taken, start = 1, 4
    end
    for field in Fields[start:end]
        if field === nothing
            taken += 1
            walk = :(path_step($walk, getfield(reader, :indices)[$taken]))
        else
            walk = :(path_step($walk, $(QuoteNode(field))))
        end
    end
    return walk
end

"""
    ParamEntry

One flattened parameter: the symbolic `param`, a `read(sam)` callable that
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
    "the run the parameters are read from"
    sam::SymbolicAWEModel
    entries::Vector{ParamEntry}
    cache::Dict{Any, Any}
end
ParamRegistry(sam::SymbolicAWEModel) = ParamRegistry(sam, ParamEntry[], Dict{Any, Any}())

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
`name(x)`. `reader` reads the live value from the run at sync time; `path`
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

"""
    leaf_param!(reg, path, reader, value)

Record the flat parameter for a leaf read at `path`, named and memoised by the full
path. Two reads of the same field on different components are therefore different
symbols, which is what lets one kernel read both endpoints of a joint without them
collapsing onto one parameter.
"""
leaf_param!(reg::ParamRegistry, path::Tuple, reader, value) =
    register_leaf!(reg, path, param_name(path), reader, value, path)

"""
    param_computed!(reg, name, reader)

Escape hatch for a value that is not a plain field read — `reader(sam)` computes it
(e.g. a [`WindFactorReader`](@ref) building a callable wind-factor from the atmospheric
model). `reader` must be a named struct (serialisable), not a closure over the run.
"""
param_computed!(reg::ParamRegistry, name::Symbol, reader) =
    register_leaf!(reg, name, name, reader, reader(reg.sam), (name,))

# ==================== BUILD-TIME VIEW ==================== #

"""Types the view descends *through* (everything else is a leaf)."""
param_descend(x) = x isa Rows || x isa AbstractAeroModel ||
                   x isa AbstractWinchModel || x isa WinchModel ||
                   x isa VSMEngine || x isa Settings || x isa UnsteadyAero

param_name(path::Tuple) = Symbol("p_", join(path, "_"))

"""
Top-level `params` view wrapping a [`ParamRegistry`](@ref): `params.segments[i].rest_length`
reads `run_value(sam, :segments, i, :rest_length)` and records the parameter that stands
for it (build-time only).
"""
struct ParamView
    reg::ParamRegistry
end

"""A partial path into the run on its way to a parameter."""
struct PathView
    reg::ParamRegistry
    path::Tuple
end

Base.getproperty(view::ParamView, sym::Symbol) =
    sym === :reg ? getfield(view, :reg) : PathView(getfield(view, :reg), (sym,))

Base.getindex(view::PathView, idx::Integer) =
    PathView(getfield(view, :reg), (getfield(view, :path)..., Int(idx)))

"""
    param_unknowns(params)

The symbolic parameters a `params` view recorded so far, in insertion order — passed
as the parameter list of a component `System` so every `params.…` read it made is
declared.
"""
param_unknowns(params::ParamView) = Any[entry.param for entry in params.reg.entries]

function Base.getproperty(view::PathView, sym::Symbol)
    (sym === :reg || sym === :path) && return getfield(view, sym)
    reg = getfield(view, :reg)
    path = canonical_path(reg.sam, (getfield(view, :path)..., sym))
    cached = get(reg.cache, path, nothing)
    cached === nothing || return cached
    value = read_path(reg.sam, path)
    param_descend(value) && return PathView(reg, path)
    return leaf_param!(reg, path, leaf_reader(path, value), value)
end

"""The reader of the leaf `value` at `path`: a [`PathReader`](@ref), copying into a
buffer of its own a vector built on read, and a [`StiffnessReader`](@ref) for a
segment's numeric `unit_stiffness`."""
leaf_reader(path, value) = PathReader(path)
leaf_reader(path, value::SVector) = CopyReader(PathReader(path), MVector(Tuple(value)))
leaf_reader(path, value::Real) = is_stiffness_path(path) ?
    StiffnessReader(PathReader(path)) : PathReader(path)

"""Whether `path` reads a segment's `unit_stiffness`."""
is_stiffness_path(path) = length(path) == 3 && path[1] === :segments &&
    path[3] === :unit_stiffness

"""
    CopyReader(reader, buffer)

Reads `reader` into `buffer` and returns it, so a vector of a run's state, read from
its columns, reaches a parameter buffer as the same array every sync rather than a
new one.
"""
struct CopyReader{R, B}
    reader::R
    buffer::B
end

(reader::CopyReader)(sam) = (reader.buffer .= reader.reader(sam))

"""
    StiffnessReader(reader)

Reads the numeric stiffness per unit length of a segment `reader` reads as a `Float64`.
"""
struct StiffnessReader{R}
    reader::R
end

(reader::StiffnessReader)(sam) = to_float64(reader.reader(sam))

"""The [`PathReader`](@ref) `reader` reads through, or `nothing` for a computed one."""
path_reader(reader) = nothing
path_reader(reader::PathReader) = reader
path_reader(reader::CopyReader) = reader.reader
path_reader(reader::StiffnessReader) = reader.reader

"""A reader of the same kind as `reader` reading `path` instead."""
read_at(::PathReader, path) = PathReader(path)
read_at(reader::CopyReader, path) = CopyReader(PathReader(path), similar(reader.buffer))
read_at(::StiffnessReader, path) = StiffnessReader(PathReader(path))

# ==================== SYNC ==================== #

"""
    ReaderGroup(targets, readers)

The readers sharing one concrete type, with the index each fills in the vector it
writes to.

Readers are closures of many types, so a `Vector{Any}` of them costs a dynamic
dispatch and a boxed return per element. Grouping by type lets
[`sync_readers!`](@ref) take one dispatch per group and then run a fully inferred
loop: at SK100 scale that is a dozen dispatches instead of 22756, twice a step.
"""
struct ReaderGroup{R}
    targets::Vector{Int}
    readers::Vector{R}
end

"""
    group_readers(targets, readers) -> Vector{ReaderGroup}

Split the parallel `targets`/`readers` vectors into one [`ReaderGroup`](@ref) per
concrete reader type, keeping the original order within each group.
"""
function group_readers(targets, readers)
    position = IdDict{Any, Int}()
    groups = ReaderGroup[]
    for (target, reader) in zip(targets, readers)
        key = typeof(reader)
        index = get(position, key, 0)
        if index == 0
            push!(groups, ReaderGroup(Int[], Vector{key}()))
            position[key] = index = length(groups)
        end
        push!(groups[index].targets, target)
        push!(groups[index].readers, reader)
    end
    return groups
end

"""
    sync_readers!(group::ReaderGroup, destination, sam)

Write every reader of `group` into `destination` at its target index. Called with
a concretely typed `group`, which is the function barrier the grouping exists for.
"""
function sync_readers!(group::ReaderGroup, destination, sam)
    targets, readers = group.targets, group.readers
    @inbounds for k in eachindex(targets)
        destination[targets[k]] = readers[k](sam)
    end
    return nothing
end

"""
    ParamGroup

A `setp` setter plus the type-grouped readers and preallocated value buffer for
one parameter kind. `eltype` is `SimFloat` for numeric scalars, `Any` for
arrays/callables.
"""
struct ParamGroup{Setter, Buf}
    setter::Setter
    groups::Vector{ReaderGroup}
    buffer::Buf
end

"""The per-kind sync groups of one sync class, each `nothing` where it holds none."""
struct ParamKinds{S, A, C}
    scalar::S
    array::A
    callable::C
end

"""
    ParamSync(live, aero, run)

The [`ParamKinds`](@ref) of each sync class: the parameters a run's `SysState` gives
(`live`, written every step), those an aerodynamic refresh changes (`aero`), and the
rest, which a run writes once when it starts and then only where
[`set_run_value!`](@ref) writes one (`run`).
"""
struct ParamSync{L, A, R}
    live::L
    aero::A
    run::R
end

"""Whether `path` reads a point's frozen aerodynamic force, which only an aerodynamic
refresh changes: a sync at another attitude would turn it with the wing."""
refreshed_path(path) = length(path) >= 3 && path[1] === :points && path[3] === :aero_force_b

"""The `(collection, field)` pairs whose parameters read a column of a run's `SysState`,
or what follows from them."""
const LIVE_FIELDS = Set([(:points, :pos_w), (:points, :wind_vec), (:points, :total_mass),
    (:bodies, :wind_vec), (:wings, :aero_force_b), (:wings, :aero_moment_b),
    (:bodies, :pos_w), (:bodies, :Q_b_to_w), (:bodies, :com_w),
    (:spanwise_stations, :twist), (:tethers, :len)])

"""The fields of a wing an aerodynamic refresh changes: its aero mode and its engine's."""
const AERO_FIELDS = Set([:aero, :aero_y, :aero_x, :aero_jac, :unsteady])

"""
    sync_class(path) -> Symbol

When the parameter read at `path` changes: `:live` every step, `:aero` with each
aerodynamic refresh, `:run` only where a run writes it.
"""
function sync_class(path)
    length(path) >= 3 || return :run
    (path[1], path[3]) in LIVE_FIELDS && return :live
    (refreshed_path(path) || path[3] in AERO_FIELDS) && return :aero
    return :run
end

"""
    survivor_index(sys) -> Dict{String, Vector{param}}

Map each name to *every* parameter surviving `mtkcompile` under it, keyed by both the
full name and the leaf name (after the last `₊` namespace separator) so a registry's
bare param matches its namespaced counterparts from subsystems.

A field read by both the parent equations and a subsystem is memoised to one registry
entry but compiles to one parameter per reader, all sharing a leaf name. Keeping only
the last would leave the others frozen at their build-time value, and whichever
observed equation reads a frozen one then writes it back over the run in
[`update_run_state!`](@ref).
"""
function survivor_index(sys)
    index = Dict{String, Vector{Any}}()
    record(key, param) = push!(get!(() -> Any[], index, key), param)
    for p in parameters(sys)
        name = string(ModelingToolkit.getname(ModelingToolkit.unwrap(p)))
        record(name, p)
        sep = findlast('₊', name)
        sep === nothing || record(name[nextind(name, sep):end], p)
    end
    return index
end

entry_name(param) = string(ModelingToolkit.getname(ModelingToolkit.unwrap(param)))

"""
    build_param_sync(sys, registry) -> ParamSync | Nothing

Build the per-class, per-kind sync groups from the compiled system and the registry.
Pruned parameters (no surviving equation references them) are dropped, so a setter
never touches a parameter absent from the buffer.
"""
function build_param_sync(sys, registry::ParamRegistry)
    isempty(registry.entries) && return nothing
    index = survivor_index(sys)
    function group(entries, ::Type{Buf}) where {Buf}
        survivors, readers = Any[], Any[]
        for entry in entries
            for survivor in get(index, entry_name(entry.param), ())
                push!(survivors, survivor); push!(readers, entry.read)
            end
        end
        isempty(survivors) && return nothing
        return ParamGroup(setp(sys, survivors), group_readers(eachindex(readers), readers),
                          Buf(undef, length(survivors)))
    end
    of(class, kind) = filter(e -> e.kind === kind && sync_class(e.path) === class,
                             registry.entries)
    kinds(class) = ParamKinds(group(of(class, :scalar), Vector{SimFloat}),
                              group(of(class, :array), Vector{Any}),
                              group(of(class, :callable), Vector{Any}))
    return ParamSync(kinds(:live), kinds(:aero), kinds(:run))
end

"""
    sync_params!(sync, target, sam)

Copy the parameters a run's `SysState` gives from the live `sam` into
`target`'s parameter buffers (`target` is an `ODEProblem` or an `ODEIntegrator`). A
no-op when there are no flattened parameters.
"""
sync_params!(::Nothing, target, sam) = nothing
sync_params!(sync::ParamSync, target, sam) =
    sync_kinds!(sync.live, target, sam)

"""
    sync_aero_params!(sync, target, sam)

[`sync_params!`](@ref), and the parameters an aerodynamic refresh changes: to call right
after a refresh.
"""
sync_aero_params!(::Nothing, target, sam) = nothing
function sync_aero_params!(sync::ParamSync, target, sam)
    sync_kinds!(sync.aero, target, sam)
    return sync_params!(sync, target, sam)
end

"""
    sync_run_params!(sync, target, sam)

Every parameter, from `sam`: where a run starts, before it holds any of its own.
"""
sync_run_params!(::Nothing, target, sam) = nothing
function sync_run_params!(sync::ParamSync, target, sam)
    sync_kinds!(sync.run, target, sam)
    return sync_aero_params!(sync, target, sam)
end

function sync_kinds!(kinds::ParamKinds, target, sam)
    sync_group!(kinds.scalar, target, sam)
    sync_group!(kinds.array, target, sam)
    sync_group!(kinds.callable, target, sam)
    return nothing
end

sync_group!(::Nothing, target, sam) = nothing
function sync_group!(group::ParamGroup, target, sam)
    buffer = group.buffer
    for reader_group in group.groups
        sync_readers!(reader_group, buffer, sam)
    end
    group.setter(target, buffer)
    return nothing
end

# ==================== RUN PARAMETERS ==================== #

"""
    build_parameter_access(sync, prob, path) -> ParameterAccess

The access to the parameter read from `path` in an integrator of `prob`, a
[`ProbWithAttributes`](@ref), whose parameters `sync` syncs: `nothing` where no
parameter is read from `path`.
"""
build_parameter_access(::Nothing, prob, path) = nothing
function build_parameter_access(::ParamSync, prob, path)
    survivors = get(survivor_index(prob.sys), string(param_name(path)), nothing)
    return isnothing(survivors) ? nothing : MonolithParameter(survivors)
end

"""The value of the parameter `access` stands for in `integrator`."""
read_parameter(access::MonolithParameter, integrator) =
    getp(integrator, first(access.parameters))(integrator)

"""Write `value` to every parameter `access` stands for in `integrator`."""
function write_parameter!(access::MonolithParameter, integrator, value)
    setp(integrator, access.parameters)(integrator,
                                        fill(value, length(access.parameters)))
    return nothing
end
