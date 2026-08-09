# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

"""
    SlotMap(symbols)

Name → position map for one of a kernel's buffers. Array-valued symbolics are
scalarized before codegen, so `pos(t)[2]` is stored under `:pos_2`; the *group*
`:pos` additionally resolves to all of `pos_1 … pos_n` in order, which is how a
whole vector is wired in one call. A name that carries no trailing index is its
own one-element group. A group is ordered by the element index its name carries,
never by buffer position, because simplification reorders variables freely.
"""
struct SlotMap
    names::Vector{Symbol}
    index::Dict{Symbol, Int}
    groups::Dict{Symbol, Vector{Int}}
end

function SlotMap(symbols)
    names = Symbol[KernelCodegen.scalar_name(sym) for sym in symbols]
    index = Dict(name => k for (k, name) in enumerate(names))
    members = Dict{Symbol, Vector{Tuple{Int, Int}}}()
    for (k, name) in enumerate(names)
        base, component = group_name(name)
        base === name || push!(get!(members, base, Tuple{Int, Int}[]), (component, k))
    end
    groups = Dict(base => [slot for (_, slot) in sort!(entries)]
                  for (base, entries) in members)
    return SlotMap(names, index, groups)
end

"""
    group_name(name) -> (base, component)

Split a scalarized slot name into the vector it belongs to and its element index:
`:pos_2` → `(:pos, 2)`. One trailing `_<digits>` is stripped, which is exactly what
scalarizing a vector variable appends. A name without such a suffix is returned
unchanged with element index `0`.
"""
function group_name(name::Symbol)
    text = string(name)
    cut = findlast('_', text)
    (cut === nothing || cut == length(text)) && return name, 0
    tail = text[(cut + 1):end]
    all(isdigit, tail) || return name, 0
    return Symbol(text[1:(cut - 1)]), parse(Int, tail)
end

Base.length(map::SlotMap) = length(map.names)

"""
    slots(map, name) -> Vector{Int}

The buffer positions `name` occupies: the whole vector when `name` is a group
(`:pos` → the slots of `pos_1 … pos_n`), otherwise the single named slot.
"""
function slots(map::SlotMap, name::Symbol)
    group = get(map.groups, name, nothing)
    group === nothing || return group
    position = get(map.index, name, nothing)
    position === nothing && error("no slot named $name; available: $(map.names)")
    return [position]
end

"""Whether `map` has any slot under `name`, as a slot or as a group."""
has_slot(map::SlotMap, name::Symbol) =
    haskey(map.index, name) || haskey(map.groups, name)

"""
    KernelReads(output_input, output_state, dstate_input, dstate_state)

Which of a kernel's inputs and states each of its outputs, and each of its state
derivatives, actually reads, as slot lists. This is the component's own Jacobian
sparsity; walked out through the wiring it gives the model's
([`state_dependencies`](@ref)), and collapsed to one bool per input it gives the
[`ComponentKernel`](@ref)'s `input_feeds_output`.
"""
struct KernelReads
    output_input::Vector{Vector{Int}}
    output_state::Vector{Vector{Int}}
    dstate_input::Vector{Vector{Int}}
    dstate_state::Vector{Vector{Int}}
end

"""Whether each of `count` inputs is read by any output of a kernel with `reads`."""
function feeds_output(reads::KernelReads, count::Int)
    feeds = zeros(Bool, count)
    for slots in reads.output_input
        feeds[slots] .= true
    end
    return feeds
end

"""
    ComponentKernel

One compiled component *type*. All three maps take
`(target, u, input, numeric, callables, instances, batch, t)` and run over every
instance named in `batch`: `rhs!` integrates the component's own states and is
`nothing` for a stateless component, `out!` writes its declared outputs and `obs!`
its remaining observed variables. The `SlotMap`s name each buffer's positions,
`input_feeds_output` marks the inputs `out!` reads — all the schedule needs to know
about the component's internal dependencies — and `reads` resolves that per slot.
"""
struct ComponentKernel{F, G, O, M}
    name::Symbol
    rhs!::F
    out!::G
    obs!::O
    mass_matrix::M
    states::SlotMap
    inputs::SlotMap
    outputs::SlotMap
    observables::SlotMap
    params::SlotMap
    callables::SlotMap
    param_defaults::Vector{SimFloat}
    callable_defaults::Vector{Any}
    input_feeds_output::Vector{Bool}
    reads::KernelReads
end

"""
    required_default(kernel, param, value)

`value` as a `SimFloat`, or an error naming the parameter that has none. Every
parameter needs a build-time default: the kernel's buffer is seeded from them and
only then overwritten per instance, so a missing one would leave a silent zero.
"""
function required_default(kernel, param, value)
    value === nothing && error("kernel $kernel: parameter $param has no default")
    return SimFloat(value)
end

"""
    compile_kernel(system, inputs, outputs; name, verbose=false)

Compile one component `System` into a [`ComponentKernel`](@ref). `inputs` and
`outputs` name variables of `system`; array-valued ones are scalarized, so an
output named `:pos` becomes the slots `pos_1 … pos_3`. Any parameter without a
build-time default is an error — the kernel's parameter buffer is seeded from
those defaults and only then overwritten per instance.
"""
function compile_kernel(system, inputs, outputs; name = nameof(system),
                        verbose = false)
    gen = KernelCodegen.generate_io_function(system, inputs, outputs; verbose)
    defaults = SimFloat[required_default(name, param, value)
                        for (param, value) in zip(gen.params, gen.param_defaults)]
    reads = KernelReads(gen.reads.output_input, gen.reads.output_state,
                        gen.reads.dstate_input, gen.reads.dstate_state)
    return ComponentKernel(name, gen.f, gen.g, gen.obs, gen.mass_matrix,
        SlotMap(gen.states), SlotMap(gen.inputs), SlotMap(gen.outputs),
        SlotMap(gen.obsstates), SlotMap(gen.params), SlotMap(gen.callable_params),
        defaults, Vector{Any}(gen.callable_defaults),
        feeds_output(reads, length(gen.inputs)), reads)
end
