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
    kernel_param_slots(kernel) -> Dict{Symbol, Vector{Int}}

Map each of `kernel`'s numeric parameter *names* to every slot it occupies, so a
registry entry can be matched to the slots it survived as. Names, not symbolics: a
parameter a nested component declared comes back namespaced
(`aero₊p_wings_1_aero_v_ind_1_2`) while the registry holds the bare symbol it was
created under, so the two are different objects for the same field.
"""
kernel_param_slots(kernel::ComponentKernel) = name_slots(kernel.params.names)

"""
    name_slots(names) -> Dict{Symbol, Vector{Int}}

`name => slots` for one buffer's names, where a namespaced name is registered under
its leaf (the part after the last `₊`) as well as in full. One field can occupy
several slots: a component that both reads a field itself and passes it to a nested
subsystem declares it twice, once bare and once namespaced, and only the namespaced
copy is the one the subsystem's equations read. Writing every match is right rather
than lucky — the names agree because the *field* is the same, so the value is too.
"""
function name_slots(names)
    slots = Dict{Symbol, Vector{Int}}()
    for (k, name) in enumerate(names)
        push!(get!(Vector{Int}, slots, name), k)
        leaf = leaf_name(name)
        leaf === name || push!(get!(Vector{Int}, slots, leaf), k)
    end
    return slots
end

"""The slots a parameter symbolic occupies in a [`name_slots`](@ref) map, by the
name it scalarizes to; empty when the kernel does not carry it."""
matched_slots(slots, param) =
    get(slots, KernelCodegen.scalar_name(param), EMPTY_SLOTS)

"""No slots — the shared empty result of [`matched_slots`](@ref)."""
const EMPTY_SLOTS = Int[]

"""The part of a namespaced name after the last `₊`, or the name itself."""
function leaf_name(name::Symbol)
    text = string(name)
    cut = findlast('₊', text)
    return cut === nothing ? name : Symbol(text[nextind(text, cut):end])
end

"""
    kernel_callable_slots(kernel) -> Dict{Symbol, Vector{Int}}

As [`kernel_param_slots`](@ref), for the kernel's callable parameters.
"""
kernel_callable_slots(kernel::ComponentKernel) = name_slots(kernel.callables.names)

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
                reader = entry_reader(entry, index_map, k)
                for slot in matched_slots(numeric_slots, element)
                    push!(numeric, (slot, reader))
                end
            end
        else
            target = entry.kind === :callable ? callable_slots : numeric_slots
            bound = entry.kind === :callable ? callable : numeric
            reader = entry_reader(entry, index_map, nothing)
            for slot in matched_slots(target, entry.param)
                push!(bound, (slot, reader))
            end
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
    remapped = remap_path(entry.read.path, index_map)
    return PathReader(element === nothing ? remapped : (remapped..., element))
end

"""
    remap_path(path, index_map) -> Tuple

`path` with the index following each container named in `index_map` swapped for that
instance's. The container may sit at any depth, so a panel addressed as
`wings[1].aero.panels[3].v_ind` is remapped on `panels` while `wings` stays put.
"""
function remap_path(path::Tuple, index_map)
    isempty(index_map) && return path
    remapped = collect(Any, path)
    for k in 1:(length(remapped) - 1)
        index = get(index_map, remapped[k], nothing)
        index === nothing && continue
        remapped[k + 1] isa Int || continue
        remapped[k + 1] = index
    end
    return Tuple(remapped)
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
    for k in 1:(length(path) - 1)
        haskey(index_map, path[k]) && path[k + 1] isa Int || continue
        previous = get!(built_from, path[k], path[k + 1])
        previous == path[k + 1] || error(
            "kernel reads $(path[k])[$(path[k + 1])] and $(path[k])[$previous]; a " *
            "kernel may read only one component of a container it is instanced over")
    end
    return nothing
end
