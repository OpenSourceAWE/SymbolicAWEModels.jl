# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# The parameter object the runtime evaluates against, and the binding of a kernel's
# parameters to one instance. A kernel is compiled from one representative
# component, so its parameters carry that component's index in their recorded path
# (`(:points, 3, :extra_mass)`). Every other instance of the same kernel reads the
# same fields of *its* component, which is the recorded path with the container
# index swapped. The runtime address is `(instance, slot)` into our own buffer, so
# two reads of the same field name can never collide.
#
# This is the lowest layer that mentions the `SystemStructure`: the runtime below it
# sees only integers and buffers, and only ever reads `p.numeric` / `p.callables`.

"""
    ScheduledParams(numeric, callables)

The parameter object the assembled `ODEProblem` carries. `numeric` is the flat
buffer the generated kernels index — its layout is ours, `(instance, slot) → flat
index`, so the same field read on two instances can never collide. `callables`
holds the polars and rigidity laws, which cannot live in a numeric buffer: one
entry per kernel, each a vector of that kernel's instances' callables, so the
element type stays concrete and calling a polar dispatches statically.
"""
struct ScheduledParams{C <: Tuple}
    numeric::Vector{SimFloat}
    callables::C
end

"""
    ScheduledParamSync(slots, readers, callable_targets, callable_slots,
                       callable_readers)

Copies live `SystemStructure` fields into a [`ScheduledParams`](@ref):
`readers[k](sys_struct)` supplies `numeric[slots[k]]`, and a callable is written
into `callable_targets[k][callable_slots[k]]` — the instance's own callable vector,
held directly so the write needs no index arithmetic. This is the scheduled
backend's parameter sync, applied every step through the same
`ProbWithAttributes` machinery as the monolith's.
"""
struct ScheduledParamSync
    slots::Vector{Int}
    readers::Vector{Any}
    callable_targets::Vector{Any}
    callable_slots::Vector{Int}
    callable_readers::Vector{Any}
end

function sync_params!(sync::ScheduledParamSync, target, sys_struct::SystemStructure)
    numeric = target.p.numeric
    @inbounds for k in eachindex(sync.slots)
        numeric[sync.slots[k]] = sync.readers[k](sys_struct)
    end
    for k in eachindex(sync.callable_slots)
        sync.callable_targets[k][sync.callable_slots[k]] =
            sync.callable_readers[k](sys_struct)
    end
    return nothing
end

"""
    ScaledReader(reader, factor)

Reads `reader` and multiplies by `factor`. Used where a parameter is a fixed share
of a struct field — an unwinched tether segment's rest length is its tether's
length over the segment count.
"""
struct ScaledReader{R}
    reader::R
    factor::SimFloat
end

(reader::ScaledReader)(sys_struct) = reader.reader(sys_struct) * reader.factor

"""
    kernel_param_slots(kernel) -> Dict{Any, Int}

Map each of `kernel`'s numeric parameter symbolics to its slot in the kernel's
parameter buffer, so a registry entry can be matched to the slot it survived as.
"""
function kernel_param_slots(kernel::ComponentKernel)
    return Dict{Any, Int}(ModelingToolkit.unwrap(sym) => k
                          for (k, sym) in enumerate(kernel.param_syms))
end

"""
    kernel_callable_slots(kernel) -> Dict{Any, Int}

As [`kernel_param_slots`](@ref), for the kernel's callable parameters.
"""
function kernel_callable_slots(kernel::ComponentKernel)
    return Dict{Any, Int}(ModelingToolkit.unwrap(sym) => k
                          for (k, sym) in enumerate(kernel.callable_syms))
end

"""
    instance_readers(kernel, registry, index_map)

The live readers for one instance of `kernel`. Returns `(numeric, callable)`, each a
vector of `(slot, reader)` pairs: `reader(sys_struct)` is the current value of the
parameter at that slot of the instance's slice. `index_map` gives this instance's
index in each container the kernel reads (`Dict(:points => 7)`); parameters the
kernel minted itself (a pulley's damping, say) are not in the registry and keep
their build-time default.
"""
function instance_readers(kernel::ComponentKernel, registry::ParamRegistry, index_map)
    numeric_slots = kernel_param_slots(kernel)
    callable_slots = kernel_callable_slots(kernel)
    numeric = Tuple{Int, Any}[]
    callable = Tuple{Int, Any}[]
    built_from = Dict{Symbol, Int}()
    for entry in registry.entries
        record_build_index!(built_from, entry, index_map)
        if entry.kind === :array
            for (k, element) in enumerate(vec(collect(Symbolics.scalarize(entry.param))))
                slot = get(numeric_slots, ModelingToolkit.unwrap(element), nothing)
                slot === nothing && continue
                push!(numeric, (slot, entry_reader(entry, index_map, k)))
            end
        else
            target = entry.kind === :callable ? callable_slots : numeric_slots
            slot = get(target, ModelingToolkit.unwrap(entry.param), nothing)
            slot === nothing && continue
            push!(entry.kind === :callable ? callable : numeric,
                  (slot, entry_reader(entry, index_map, nothing)))
        end
    end
    return numeric, callable
end

"""
    entry_reader(entry, index_map, element)

The reader for one registry `entry` on the instance described by `index_map`: the
recorded path with its container index swapped, plus `element` when the entry is an
array whose components occupy separate slots. A computed entry (not a plain field
read) is instance-independent and is reused as it stands.
"""
function entry_reader(entry::ParamEntry, index_map, element)
    if !(entry.read isa PathReader)
        return element === nothing ? entry.read : ElementReader(entry.read, element)
    end
    path = entry.read.path
    index = get(index_map, first(path), nothing)
    remapped = index === nothing ? path : (first(path), index, path[3:end]...)
    return PathReader(element === nothing ? remapped : (remapped..., element))
end

"""
    record_build_index!(built_from, entry, index_map)

Record which component of each remapped container the kernel was built from, and
error on a second one. Two components of one container in a single kernel would
both be swapped to the same instance index, silently reading the wrong one; the
components are written so this cannot happen, and this says so if that changes.
"""
function record_build_index!(built_from, entry::ParamEntry, index_map)
    entry.read isa PathReader || return nothing
    path = entry.read.path
    (length(path) >= 2 && path[2] isa Int && haskey(index_map, first(path))) ||
        return nothing
    previous = get!(built_from, first(path), path[2])
    previous == path[2] || error(
        "kernel reads $(first(path))[$(path[2])] and $(first(path))[$previous]; a " *
        "kernel may read only one component of a container it is instanced over")
    return nothing
end
