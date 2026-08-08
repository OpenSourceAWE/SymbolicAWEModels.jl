# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# The runtime: four concepts and nothing else — a *kernel* (one compiled
# component type, kernel.jl), an *instance* (one occurrence of a kernel, with its
# slice of every buffer), the *wiring* (which output slots sum into which input
# slot) and the *schedule* (the layered order the output maps run in). Nothing
# here knows about points, segments or bodies; the assembler translates topology
# into instances and wiring, and the runtime sees integers and buffers.

"""
    ComponentInstance

One occurrence of a kernel. `position` counts the instance among its own kernel's,
which is how it finds its callable parameters; every other field is this instance's
contiguous slice of the corresponding global buffer — `states` of `u`/`du`, `params`
of the numeric parameter vector, and so on.
"""
struct ComponentInstance
    kernel::Int
    position::Int
    states::UnitRange{Int}
    inputs::UnitRange{Int}
    outputs::UnitRange{Int}
    observables::UnitRange{Int}
    params::UnitRange{Int}
end

"""
    Wiring(sources, weights, offsets)

Flattened compressed-row map from output slots to input slots: input slot `k` is
the weighted sum of the output-buffer slots `sources[offsets[k]:offsets[k + 1] - 1]`.
Most weights are `1`; a weighted reference point — a wing frame fitted from a
blend of structural points — is the reason they exist at all, and having them here
is why such a point needs no component of its own. An input with no sources reads
zero, which is how an unconnected slot gets its value without anyone writing it.
"""
struct Wiring
    sources::Vector{Int}
    weights::Vector{SimFloat}
    offsets::Vector{Int}
end

"""
    GatherPlan(copy_targets, copy_sources, general)

One gather pass, split by how much work each slot needs. A slot fed by a single
unit-weighted source is a plain copy, which nearly all of them are — on a replicated
four-kite model the wiring averages 1.45 sources per slot, so the loop setup guarding
the sum costs more than the sum. The rest keep the weighted form in `general`. A slot
with no source appears in neither list: only a gather writes an input, so it keeps the
zero its buffer was built with and never needs touching.
"""
struct GatherPlan
    copy_targets::Vector{Int}
    copy_sources::Vector{Int}
    general::Vector{Int}
end

"""
    gather_plan(wiring, slots) -> GatherPlan

Sort `slots` into the copies, the weighted sums and the untouched, once, so that
[`gather!`](@ref) does not rediscover the shape of the wiring on every call.
"""
function gather_plan(wiring::Wiring, slots)
    copy_targets = Int[]
    copy_sources = Int[]
    general = Int[]
    for k in slots
        first_edge = wiring.offsets[k]
        edges = wiring.offsets[k + 1] - first_edge
        if edges == 0
            continue
        elseif edges == 1 && isone(wiring.weights[first_edge])
            push!(copy_targets, k)
            push!(copy_sources, wiring.sources[first_edge])
        else
            push!(general, k)
        end
    end
    return GatherPlan(copy_targets, copy_sources, general)
end

"""
    gather!(input, output, wiring, plan)

Sum every wired output slot into the input slots `plan` covers. This is the whole
coupling model: a segment writing force into a point, a ride point forwarding a
moment into its body and a winch driving a tether length all go through it, and
because it indexes by slot there is no uniform I/O width to pad to. A plan holds what
the next pass is about to read rather than the whole buffer, and each slot belongs to
exactly one plan — the schedule orders every producer before the earliest consumer, so
a slot that two layers read is already final when the first of them runs.
"""
function gather!(input, output, wiring::Wiring, plan::GatherPlan)
    targets = plan.copy_targets
    from = plan.copy_sources
    @inbounds for i in eachindex(targets, from)
        input[targets[i]] = output[from[i]]
    end
    sources = wiring.sources
    weights = wiring.weights
    offsets = wiring.offsets
    @inbounds for k in plan.general
        total = zero(eltype(input))
        for j in offsets[k]:(offsets[k + 1] - 1)
            total += weights[j] * output[sources[j]]
        end
        input[k] = total
    end
    return nothing
end

"""
    KernelSystem

An assembled model. `kernels` is a tuple so the evaluation loop unrolls over it and
every kernel call is statically dispatched. A *batch* is the instance list of one
kernel: `layers[l][k]` are the instances of kernel `k` whose outputs run in layer
`l`, and `state_batches[k]` the stateful instances of kernel `k`. `layer_inputs[l]` is
the [`GatherPlan`](@ref) that layer needs and `final_inputs` what the derivative and
observable passes still need on top, so each [`gather!`](@ref) touches only the slots
about to be used and no slot is gathered twice in a call. The `_active` lists hold the positions of the non-empty
batches of each pass, which is what [`run_batches!`](@ref) walks — most kernels are
idle in any one layer and a pass over the whole tuple charged for them.
"""
struct KernelSystem{K <: Tuple, N, M}
    kernels::K
    instances::Vector{ComponentInstance}
    layers::Vector{NTuple{N, Vector{Int}}}
    layer_active::Vector{Vector{Int}}
    layer_inputs::Vector{GatherPlan}
    final_inputs::GatherPlan
    state_batches::NTuple{N, Vector{Int}}
    state_active::Vector{Int}
    observable_batches::NTuple{N, Vector{Int}}
    observable_active::Vector{Int}
    wiring::Wiring
    n_states::Int
    n_inputs::Int
    n_outputs::Int
    n_observables::Int
    n_params::Int
    mass_matrix::M
    sparsity::SparseMatrixCSC{Bool, Int}
    compile_seconds::Vector{Float64}
end

# ======================= assembly ======================= #

"""
    SystemBuilder()

Accumulates kernels, instances and connections, then [`build_system`](@ref)s them
into a [`KernelSystem`](@ref). Buffer slices are handed out as instances are
added, so a connection can be recorded as soon as both endpoints exist.
"""
mutable struct SystemBuilder
    kernels::Vector{ComponentKernel}
    counts::Vector{Int}
    instances::Vector{ComponentInstance}
    targets::Vector{Int}
    sources::Vector{Int}
    weights::Vector{SimFloat}
    n_states::Int
    n_inputs::Int
    n_outputs::Int
    n_observables::Int
    n_params::Int
    compile_seconds::Vector{Float64}
    verbose::Bool
end

SystemBuilder(; verbose = false) =
    SystemBuilder(ComponentKernel[], Int[], ComponentInstance[], Int[], Int[],
                  SimFloat[], 0, 0, 0, 0, 0, Float64[], verbose)

"""
    add_kernel!(builder, kernel, seconds) -> Int

Register a compiled component type and return its index in the kernel table.
`seconds` is what compiling it cost, kept for the build-time breakdown.
"""
function add_kernel!(builder::SystemBuilder, kernel::ComponentKernel, seconds)
    push!(builder.kernels, kernel)
    push!(builder.counts, 0)
    push!(builder.compile_seconds, seconds)
    return length(builder.kernels)
end

"""
    add_instance!(builder, kernel_idx) -> Int

Add one occurrence of a registered kernel, reserving its slice of every buffer.
Returns the instance index.
"""
function add_instance!(builder::SystemBuilder, kernel_idx::Int)
    kernel = builder.kernels[kernel_idx]
    position = (builder.counts[kernel_idx] += 1)
    states = reserve!(builder, :n_states, length(kernel.states))
    inputs = reserve!(builder, :n_inputs, length(kernel.inputs))
    outputs = reserve!(builder, :n_outputs, length(kernel.outputs))
    observables = reserve!(builder, :n_observables, length(kernel.observables))
    params = reserve!(builder, :n_params, length(kernel.params))
    push!(builder.instances, ComponentInstance(kernel_idx, position, states, inputs,
                                               outputs, observables, params))
    return length(builder.instances)
end

"""Reserve `count` slots at the end of the buffer counted by `field`."""
function reserve!(builder::SystemBuilder, field::Symbol, count::Int)
    base = getfield(builder, field)
    setfield!(builder, field, base + count)
    return (base + 1):(base + count)
end

"""
    connect!(builder, source, out_name, target, in_name; weight=1)

Add `weight` times instance `source`'s output `out_name` to instance `target`'s
input `in_name`, slot by slot. Both names may be whole vectors (`:pos` covering
`pos_1 … pos_3`) and must then have the same width. Several sources may feed one
input; they sum.
"""
function connect!(builder::SystemBuilder, source::Int, out_name::Symbol,
                  target::Int, in_name::Symbol; weight = 1.0)
    from = builder.instances[source]
    into = builder.instances[target]
    out_slots = slots(builder.kernels[from.kernel].outputs, out_name)
    in_slots = slots(builder.kernels[into.kernel].inputs, in_name)
    length(out_slots) == length(in_slots) || error(
        "cannot connect $out_name ($(length(out_slots)) slots) to $in_name " *
        "($(length(in_slots)) slots): widths differ")
    for (out_slot, in_slot) in zip(out_slots, in_slots)
        push!(builder.sources, first(from.outputs) - 1 + out_slot)
        push!(builder.targets, first(into.inputs) - 1 + in_slot)
        push!(builder.weights, weight)
    end
    return nothing
end

"""
    build_system(builder) -> KernelSystem

Turn the recorded connections into a [`Wiring`](@ref) and a layered schedule.

A slot a layer gathers is left out of `final_inputs`: the schedule puts every
producer of a feeding input before the layer that reads it, so such a slot is
already final when its layer gathered it and the derivative pass can read it as is.
What remains for the final gather is the inputs no output map reads — a point's
aggregated force and its like — which is where the coupling comes back down.
"""
function build_system(builder::SystemBuilder)
    wiring = build_wiring(builder.targets, builder.sources, builder.weights,
                          builder.n_inputs)
    layers = build_schedule(builder, wiring)
    stateful = [i for (i, inst) in enumerate(builder.instances)
                if builder.kernels[inst.kernel].rhs! !== nothing]
    observed = [i for (i, inst) in enumerate(builder.instances)
                if builder.kernels[inst.kernel].obs! !== nothing]
    layer_slots = layer_gather_slots(builder, layers)
    layer_inputs = [gather_plan(wiring, slots) for slots in layer_slots]
    settled = settled_slots(builder, wiring, layers, layer_slots)
    final_inputs = gather_plan(wiring,
        [slot for slot in read_slots(builder,
                                     Iterators.flatten((stateful, observed)))
         if !(slot in settled)])
    state_batches = batch_by_kernel(builder, stateful)
    observable_batches = batch_by_kernel(builder, observed)
    return KernelSystem(Tuple(builder.kernels), builder.instances, layers,
        [active_batches(layer) for layer in layers], layer_inputs, final_inputs,
        state_batches, active_batches(state_batches),
        observable_batches, active_batches(observable_batches),
        wiring, builder.n_states, builder.n_inputs, builder.n_outputs,
        builder.n_observables, builder.n_params, global_mass_matrix(builder),
        state_sparsity(builder, wiring, layers), builder.compile_seconds)
end

"""
    active_batches(batches) -> Vector{Int}

The positions of the batches that hold an instance, in order, so a pass skips the
kernels with no work in it instead of calling into each of them to find none.
"""
active_batches(batches) = [k for k in eachindex(batches) if !isempty(batches[k])]

"""
    settled_slots(builder, wiring, layers, layer_slots) -> Set{Int}

The slots a layer gathered whose value cannot change afterwards, so the final gather
can leave them alone. Almost every feeding input qualifies — the schedule orders a
producer before the layer that reads its output. The exception is an instance whose
own output feeds its own input, which [`build_schedule`](@ref) deliberately does not
order against itself; that slot is only complete once its layer has run, so it stays
in `final_inputs`.
"""
function settled_slots(builder::SystemBuilder, wiring::Wiring, layers, layer_slots)
    owner = Vector{Int}(undef, builder.n_outputs)
    for (i, inst) in enumerate(builder.instances)
        owner[inst.outputs] .= i
    end
    layer_of = zeros(Int, length(builder.instances))
    for (index, layer) in enumerate(layers), batch in layer, i in batch
        layer_of[i] = index
    end
    settled = Set{Int}()
    for (index, slots) in enumerate(layer_slots), slot in slots
        complete = true
        for j in wiring.offsets[slot]:(wiring.offsets[slot + 1] - 1)
            layer_of[owner[wiring.sources[j]]] < index || (complete = false; break)
        end
        complete && push!(settled, slot)
    end
    return settled
end

"""
    layer_gather_slots(builder, layers) -> Vector{Vector{Int}}

The input slots each layer has to gather before its output maps run: those the maps
actually read, minus the ones an earlier layer already gathered. Both restrictions
are what the schedule earns. An input its instance's outputs ignore — a point's
aggregated force, say — is no use to the output pass and is left to `final_inputs`;
and since every producer is ordered before the earliest consumer, a slot two layers
read is already final when the first of them runs, so gathering it twice would
recompute the same sum.
"""
function layer_gather_slots(builder::SystemBuilder, layers)
    gathered = Set{Int}()
    per_layer = Vector{Int}[]
    for layer in layers
        fresh = Int[]
        for slot in output_read_slots(builder, Iterators.flatten(layer))
            slot in gathered && continue
            push!(gathered, slot)
            push!(fresh, slot)
        end
        push!(per_layer, fresh)
    end
    return per_layer
end

"""
    output_read_slots(builder, instances) -> Vector{Int}

The input slots `instances` read *in their output maps*, ascending. `read_slots` is
the same question for the instance as a whole.
"""
function output_read_slots(builder::SystemBuilder, instances)
    slots = Int[]
    for i in instances
        inst = builder.instances[i]
        feeds = builder.kernels[inst.kernel].input_feeds_output
        for (local_slot, input_slot) in enumerate(inst.inputs)
            feeds[local_slot] && push!(slots, input_slot)
        end
    end
    return sort!(unique!(slots))
end

"""
    read_slots(builder, instances) -> Vector{Int}

The input slots `instances` read, in ascending order so the gather walks the buffer
and the wiring forwards.
"""
function read_slots(builder::SystemBuilder, instances)
    slots = Int[]
    for i in instances
        append!(slots, builder.instances[i].inputs)
    end
    return sort!(unique!(slots))
end

"""
    build_wiring(targets, sources, weights, n_inputs) -> Wiring

Invert the recorded `(target input slot, source output slot, weight)` triples into
the compressed-row form [`gather!`](@ref) walks.
"""
function build_wiring(targets, sources, weights, n_inputs)
    counts = zeros(Int, n_inputs + 1)
    for target in targets
        counts[target + 1] += 1
    end
    offsets = cumsum(counts) .+ 1
    cursor = copy(offsets)
    ordered = Vector{Int}(undef, length(sources))
    scale = Vector{SimFloat}(undef, length(sources))
    for (target, source, weight) in zip(targets, sources, weights)
        ordered[cursor[target]] = source
        scale[cursor[target]] = weight
        cursor[target] += 1
    end
    return Wiring(ordered, scale, offsets)
end

"""
    build_schedule(builder, wiring) -> Vector{NTuple{N, Vector{Int}}}

Layer the instances so that every output map runs after the maps it reads. An
instance depends on another when one of the inputs its own outputs read is fed by
that instance's outputs; inputs the outputs ignore (a point's aggregated force,
say) impose no order, which is what lets the force flow back down the same chain.
Errors on a cycle, naming the kernels involved.
"""
function build_schedule(builder::SystemBuilder, wiring::Wiring)
    owner = Vector{Int}(undef, builder.n_outputs)
    for (i, inst) in enumerate(builder.instances)
        owner[inst.outputs] .= i
    end
    count = length(builder.instances)
    dependents = [Int[] for _ in 1:count]
    remaining = zeros(Int, count)
    for (i, inst) in enumerate(builder.instances)
        feeds = builder.kernels[inst.kernel].input_feeds_output
        for (local_slot, input_slot) in enumerate(inst.inputs)
            feeds[local_slot] || continue
            for j in wiring.offsets[input_slot]:(wiring.offsets[input_slot + 1] - 1)
                producer = owner[wiring.sources[j]]
                producer == i && continue
                push!(dependents[producer], i)
                remaining[i] += 1
            end
        end
    end
    layers = NTuple{length(builder.kernels), Vector{Int}}[]
    ready = [i for i in 1:count if remaining[i] == 0]
    placed = 0
    while !isempty(ready)
        push!(layers, batch_by_kernel(builder, ready))
        placed += length(ready)
        next = Int[]
        for i in ready, j in dependents[i]
            remaining[j] -= 1
            remaining[j] == 0 && push!(next, j)
        end
        ready = next
    end
    placed == count || throw(ScheduleCycleError(
        unique([builder.kernels[builder.instances[i].kernel].name
                for i in 1:count if remaining[i] > 0])))
    return layers
end

"""
    ScheduleCycleError(kernels)

Thrown when the output-dependency graph has a cycle — a genuine algebraic loop
between components. `kernels` names the component types still unscheduled.
"""
struct ScheduleCycleError <: Exception
    kernels::Vector{Symbol}
end

function Base.showerror(io::IO, e::ScheduleCycleError)
    print(io, "the component outputs form a cycle; involved kernels: ",
          join(e.kernels, ", "))
end

"""Split `instances` into one list per kernel, in kernel-table order."""
function batch_by_kernel(builder::SystemBuilder, instances)
    grouped = [Int[] for _ in builder.kernels]
    for i in instances
        push!(grouped[builder.instances[i].kernel], i)
    end
    return Tuple(grouped)
end

"""
    state_sparsity(builder, wiring, layers) -> SparseMatrixCSC{Bool, Int}

The Jacobian's nonzero pattern: entry `(i, j)` is set when `du[i]` reads `u[j]`.
Handing it to the integrator turns a dense finite-difference Jacobian — O(n) right-
hand side evaluations and a dense factorization — into a coloured sparse one. The
diagonal is always stored, read or not, because the solver's `W = M/γ - J` writes
there and can only write where the pattern has room.
"""
function state_sparsity(builder::SystemBuilder, wiring::Wiring, layers)
    rows = state_dependencies(builder, wiring, layers)
    reading = collect(1:(builder.n_states))
    read = copy(reading)
    for (state, sources) in enumerate(rows), source in sources
        source == state && continue
        push!(reading, state)
        push!(read, source)
    end
    return sparse(reading, read, true, builder.n_states, builder.n_states)
end

"""
    state_dependencies(builder, wiring, layers) -> Vector{Vector{Int}}

For each state, the states its derivative reads. An output reaches its instance's
own states plus everything reaching the outputs wired into the inputs it reads, so
one sweep in schedule order resolves the model; the sweep repeats until nothing
grows, which costs one confirming pass and converges an instance that feeds itself.
"""
function state_dependencies(builder::SystemBuilder, wiring::Wiring, layers)
    reach = [BitSet() for _ in 1:builder.n_outputs]
    grew = true
    while grew
        grew = false
        for layer in layers, batch in layer, instance in batch
            grew |= expand_reach!(reach, builder, wiring, instance)
        end
    end
    rows = [Int[] for _ in 1:builder.n_states]
    for inst in builder.instances
        reads = builder.kernels[inst.kernel].reads
        for (slot, state) in enumerate(inst.states)
            row = BitSet()
            collect_reach!(row, reach, wiring, inst, reads.dstate_state[slot],
                           reads.dstate_input[slot])
            rows[state] = collect(row)
        end
    end
    return rows
end

"""Grow one instance's outputs' reachable states; `true` if any of them grew."""
function expand_reach!(reach, builder::SystemBuilder, wiring::Wiring, instance::Int)
    inst = builder.instances[instance]
    reads = builder.kernels[inst.kernel].reads
    grew = false
    for (slot, output) in enumerate(inst.outputs)
        before = length(reach[output])
        collect_reach!(reach[output], reach, wiring, inst, reads.output_state[slot],
                       reads.output_input[slot])
        grew |= length(reach[output]) != before
    end
    return grew
end

"""Union into `into` the states `inst` reads directly, plus those its inputs reach."""
function collect_reach!(into::BitSet, reach, wiring::Wiring, inst::ComponentInstance,
                        state_slots, input_slots)
    for slot in state_slots
        push!(into, inst.states[slot])
    end
    for slot in input_slots
        input = inst.inputs[slot]
        for j in wiring.offsets[input]:(wiring.offsets[input + 1] - 1)
            union!(into, reach[wiring.sources[j]])
        end
    end
    return into
end

"""
    global_mass_matrix(builder)

Assemble the state mass matrix from the kernels' own. Returns `I` when every
component is a plain ODE, which is the case for all of them today.
"""
function global_mass_matrix(builder::SystemBuilder)
    diagonal = ones(builder.n_states)
    for inst in builder.instances
        mass = builder.kernels[inst.kernel].mass_matrix
        mass === I && continue
        diagonal[inst.states] .= diag(mass)
    end
    return all(isone, diagonal) ? I : Diagonal(diagonal)
end

# ======================= evaluation ======================= #

"""
    KernelBuffers(input, output, observable)

The scratch buffers one evaluation needs, at one element type. Cached per
`eltype(u)` so a ForwardDiff dual pass gets its own set.
"""
struct KernelBuffers{T}
    input::Vector{T}
    output::Vector{T}
    observable::Vector{T}
end

KernelBuffers{T}(system::KernelSystem) where {T} = KernelBuffers(
    zeros(T, system.n_inputs), zeros(T, system.n_outputs),
    zeros(T, system.n_observables))

"""
    KernelRHS(system)

The callable `(du, u, p, t)` an `ODEProblem` integrates. Each schedule layer
gathers the inputs and runs its batches' output maps; once every layer has run the
inputs are complete and the stateful kernels write their derivatives.
"""
struct KernelRHS{S}
    system::S
    scratch::KernelBuffers{SimFloat}
    others::Dict{DataType, Any}
end

KernelRHS(system::KernelSystem) = KernelRHS(system,
    KernelBuffers{SimFloat}(system), Dict{DataType, Any}())

"""Scratch buffers for element type `T`, allocated on first use. The solver's own
element type is a field rather than a lookup, so the hot path allocates nothing."""
buffers(rhs::KernelRHS, ::Type{SimFloat}) = rhs.scratch

function buffers(rhs::KernelRHS, ::Type{T}) where {T}
    cached = get(rhs.others, T, nothing)
    cached === nothing || return cached::KernelBuffers{T}
    fresh = KernelBuffers{T}(rhs.system)
    rhs.others[T] = fresh
    return fresh
end

function (rhs::KernelRHS)(du, u, p, t)
    system = rhs.system
    scratch = buffers(rhs, eltype(u))
    run_layers!(system, scratch, u, p, t)
    run_batches!(derivative_call, system.kernels, system.state_active,
                 system.state_batches, p.callables, system, scratch, du, u, p, t)
    return nothing
end

"""Gather and run every schedule layer, leaving the inputs the derivative and
observable passes read complete."""
function run_layers!(system::KernelSystem, scratch, u, p, t)
    for (index, layer) in enumerate(system.layers)
        gather!(scratch.input, scratch.output, system.wiring,
                system.layer_inputs[index])
        run_batches!(output_call, system.kernels, system.layer_active[index], layer,
                     p.callables, system, scratch, scratch.output, u, p, t)
    end
    gather!(scratch.input, scratch.output, system.wiring, system.final_inputs)
    return nothing
end

"""
    run_batches!(call, kernels, active, batches, callables, system, scratch, target,
                 u, p, t)

Apply `call` to every instance of every kernel that has work, one batch at a time.
`active` holds the positions of the non-empty batches, so a kernel idle in this
pass costs nothing; walking the kernel tuple instead cost a call frame per kernel
per pass, which on a real model was most of the right-hand side.
"""
function run_batches!(call::C, kernels::Tuple, active::Vector{Int}, batches,
                      callables::Tuple, system, scratch, target, u, p,
                      t) where {C}
    for k in active
        dispatch_batch!(call, kernels, callables, batches, k, system, scratch,
                        target, u, p, t)
    end
    return nothing
end

"""
    dispatch_batch!(call, kernels, callables, batches, k, system, scratch, target,
                    u, p, t)

Run batch `k` with the types its entries have. The kernels differ in type, so
indexing the tuple at a position only known at run time would erase that type and
copy the entry into a box; a ladder of comparisons against literal positions keeps
every call to [`run_batch!`](@ref) concrete instead, at a handful of integer
compares per batch.
"""
@generated function dispatch_batch!(call, kernels::NTuple{N, Any}, callables,
                                    batches, k::Int, system, scratch, target, u, p,
                                    t) where {N}
    body = :(return nothing)
    for i in N:-1:1
        body = Expr(:if, :(k == $i),
                    :(return run_batch!(call, kernels[$i], callables[$i],
                                        batches[$i], system, scratch, target,
                                        u, p, t)),
                    body)
    end
    return body
end

"""
    run_batch!(call, kernel, values, instances, system, scratch, target, u, p, t)

Apply `call` to one kernel's `instances`, with the per-instance calls statically
dispatched. `@noinline` so that one body is shared by every batch of this kernel
rather than copied into each caller.
"""
@noinline function run_batch!(call::C, kernel::ComponentKernel, values, instances,
                              system, scratch, target, u, p, t) where {C}
    for i in instances
        call(kernel, values, system.instances[i], scratch, target, u, p, t)
    end
    return nothing
end

"""Write one instance's declared outputs into `target`."""
function output_call(kernel::ComponentKernel, callables, inst, scratch, target,
                     u, p, t)
    kernel.out!(view(target, inst.outputs), view(u, inst.states),
                view(scratch.input, inst.inputs), view(p.numeric, inst.params),
                callables[inst.position], t)
    return nothing
end

"""Write one instance's state derivatives into `target`."""
function derivative_call(kernel::ComponentKernel, callables, inst, scratch, target,
                         u, p, t)
    kernel.rhs!(view(target, inst.states), view(u, inst.states),
                view(scratch.input, inst.inputs), view(p.numeric, inst.params),
                callables[inst.position], t)
    return nothing
end

"""Write one instance's remaining observed variables into `target`."""
function observable_call(kernel::ComponentKernel, callables, inst, scratch, target,
                         u, p, t)
    kernel.obs!(view(target, inst.observables), view(u, inst.states),
                view(scratch.input, inst.inputs), view(p.numeric, inst.params),
                callables[inst.position], t)
    return nothing
end

"""
    refresh_outputs!(rhs, u, p, t) -> KernelBuffers

Re-run every output map and every observed map for the state `u`, and return the
buffers. The state getter reads component results out of these after a step, where
the integrator's last RHS evaluation need not correspond to `u`.
"""
function refresh_outputs!(rhs::KernelRHS, u, p, t)
    system = rhs.system
    scratch = buffers(rhs, eltype(u))
    run_layers!(system, scratch, u, p, t)
    run_batches!(observable_call, system.kernels, system.observable_active,
                 system.observable_batches, p.callables, system, scratch,
                 scratch.observable, u, p, t)
    return scratch
end
