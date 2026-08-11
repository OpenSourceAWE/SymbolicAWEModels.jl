# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# The analytical Jacobian. The right-hand side is a layered composition — an output
# map writes `out = g(state, input)` per instance, a constant gather sets
# `input = G · out`, and a derivative map writes `du = f(state, input)` — so with
# `A = ∂out/∂state`, `B = ∂out/∂input`, `C = ∂du/∂state` and `D = ∂du/∂input`,
#
#     Y := ∂out/∂u = A + B·G·Y
#     J  = ∂du/∂u  = C + D·G·Y
#
# `A`, `B`, `C` and `D` are block diagonal over instances, and the schedule is
# acyclic, so `Y` needs no solve: one sweep in schedule order settles it. Global
# forward mode instead pays `n_states / chunk` passes over the whole model.

"""
    KernelDuals{W, T}

One kernel's forward-mode workspace, at that kernel's own width
`W = states + inputs`. Every instance of a kernel owns a disjoint slice of `u`, of
`input` and of `output`, so seeding partial `j` on local slot `j` of *every* instance
at once returns all of their blocks from a single pass without the seeds mixing — a
compressed differentiation whose colouring is the local slot index.

The buffers are packed, one instance after another, and `instances` is the kernel's
instance list rebased onto them; the global buffers are scattered and 25× larger, and
a set of them per kernel would be hundreds of megabytes of duals to hold and to walk.
`blocks` receives the local Jacobians as `(row, slot, position)`, the rows being the
kernel's outputs followed by its state derivatives.
"""
struct KernelDuals{W, T}
    kernel::Int
    batch::Vector{Int}
    instances::Vector{ComponentInstance}
    positions::Vector{Int}
    n_state::Int
    n_input::Int
    n_output::Int
    state::Vector{T}
    input::Vector{T}
    output::Vector{T}
    derivative::Vector{T}
    seeds::NTuple{W, ForwardDiff.Partials{W, SimFloat}}
    blocks::Array{SimFloat, 3}
end

"""The slice packed buffers give the `position`-th instance of a `width`-slot buffer."""
packed_slice(position::Int, width::Int) = ((position - 1) * width + 1):(position * width)

function KernelDuals(system::KernelSystem, kernel_idx::Int, batch::Vector{Int})
    kernel = system.kernels[kernel_idx]
    n_state, n_input = length(kernel.states), length(kernel.inputs)
    n_output = length(kernel.outputs)
    width = n_state + n_input
    dual = ForwardDiff.Dual{Nothing, SimFloat, width}
    count = length(batch)
    rebased = [ComponentInstance(kernel_idx, system.instances[index].position,
                                 packed_slice(position, n_state),
                                 packed_slice(position, n_input),
                                 packed_slice(position, n_output),
                                 system.instances[index].observables,
                                 system.instances[index].params)
               for (position, index) in enumerate(batch)]
    return KernelDuals{width, dual}(kernel_idx, batch, rebased, collect(1:count),
        n_state, n_input, n_output, zeros(dual, count * n_state),
        zeros(dual, count * n_input), zeros(dual, count * n_output),
        zeros(dual, count * n_state),
        ForwardDiff.construct_seeds(ForwardDiff.Partials{width, SimFloat}),
        zeros(SimFloat, n_output + n_state, width, count))
end

"""
    fill_blocks!(duals, kernel, system, u, input, numeric, callables, t)

Run one kernel's output and derivative maps over duals seeded on its own slots,
leaving every instance's local Jacobian in `duals.blocks`. One pass per kernel
*type*, at that kernel's width, and nothing per model state.
"""
function fill_blocks!(duals::KernelDuals{W, T}, kernel::ComponentKernel,
                      system::KernelSystem, u, input, numeric, callables,
                      t) where {W, T}
    seeds = duals.seeds
    n_state, n_input, n_output = duals.n_state, duals.n_input, duals.n_output
    @inbounds for (position, index) in enumerate(duals.batch)
        instance = system.instances[index]
        base = (position - 1) * n_state
        for j in 1:n_state
            duals.state[base + j] = T(u[instance.states[j]], seeds[j])
        end
        base = (position - 1) * n_input
        for j in 1:n_input
            duals.input[base + j] = T(input[instance.inputs[j]], seeds[n_state + j])
        end
    end

    n_output == 0 || kernel.out!(duals.output, duals.state, duals.input, numeric,
                                 callables, duals.instances, duals.positions, t)
    kernel.rhs! === nothing || kernel.rhs!(duals.derivative, duals.state, duals.input,
                                           numeric, callables, duals.instances,
                                           duals.positions, t)

    blocks = duals.blocks
    @inbounds for position in eachindex(duals.batch)
        base = (position - 1) * n_output
        for row in 1:n_output
            partials = ForwardDiff.partials(duals.output[base + row])
            for slot in 1:W
                blocks[row, slot, position] = partials[slot]
            end
        end
        kernel.rhs! === nothing && continue
        base = (position - 1) * n_state
        for row in 1:n_state
            partials = ForwardDiff.partials(duals.derivative[base + row])
            for slot in 1:W
                blocks[n_output + row, slot, position] = partials[slot]
            end
        end
    end
    return nothing
end

"""
    InputGather(slots, offsets, producer, row, weight)

Where one instance's input rows of `G·Y` come from, in compressed-row form: the `k`th
gathered row is local input slot `slots[k]`, and sums `weight` times row `row` of the
`Y` block of instance `producer` over `offsets[k]:offsets[k + 1] - 1`. The same wiring
[`gather!`](@ref) walks for values, resolved to producing instances so the
substitution never searches for one.

`slots` holds only the inputs the pass in hand reads — those the output maps read for
the sweep, those the derivatives read for the composition — because both the gather
and the product that follows it are linear in the number of rows, and an input the
pass ignores contributes a row of zeros times a column of zeros.
"""
struct InputGather
    slots::Vector{Int}
    offsets::Vector{Int}
    producer::Vector{Int}
    row::Vector{Int}
    weight::Vector{SimFloat}
end

"""
    KernelJacobian

`jac(J, u, p, t)` for a [`KernelSystem`](@ref), composed from per-kernel local
Jacobians and the constant wiring instead of differentiated globally.

Per instance, `rows` holds its outputs' `Y` block over the state columns `cols`, and
`nzindex` maps its derivative block over `dstate_cols` onto the nonzeros of `matrix`.
Both column sets are the instance's own reach, so every block is small and dense and
the substitution is a sweep of little matrix products rather than sparse algebra over
the whole model.

`blocks` and `shape` restate what the [`KernelDuals`](@ref) hold, at concrete types:
a workspace's element type carries its dual width, so a vector of them is abstract
and reaching through one inside the sweep left every operation on it dynamically
dispatched. The workspaces are touched once per kernel and the blocks once per
instance, so only the latter has to be concrete.

Each component is differentiated *numerically*, so a registered leaf without a
symbolic derivative — the wind profile, an aerodynamic polar — is differentiated by
running it on duals, which is exactly where the monolith's symbolic Jacobian fails.
"""
struct KernelJacobian{R}
    rhs::R
    duals::Vector{KernelDuals}
    blocks::Vector{Array{SimFloat, 3}}
    shape::Vector{NTuple{3, Int}}
    dual_of::Vector{Int}
    position::Vector{Int}
    order::Vector{Int}
    stateful::Vector{Int}
    gather::Vector{InputGather}
    dstate_gather::Vector{InputGather}
    cols::Vector{Vector{Int}}
    rows::Vector{Matrix{SimFloat}}
    dstate_cols::Vector{Vector{Int}}
    nzindex::Vector{Matrix{Int}}
    column_of::Vector{Int}
    input_scratch::Matrix{SimFloat}
    block_scratch::Matrix{SimFloat}
    matrix::SparseMatrixCSC{SimFloat, Int}
    scratch::Vector{SimFloat}
end

"""
    build_jacobian(rhs; prn=true) -> KernelJacobian or nothing

Plan the analytical Jacobian of `rhs`: the per-kernel dual workspaces, each
instance's column support and input gather, and the map from its derivative block
onto the nonzeros of the assembled matrix.

Returns `nothing`, with a warning, when an instance's own outputs feed an input its
outputs read. The schedule deliberately does not order such an instance against
itself, so its block is an algebraic loop that one sweep would not settle; those
models keep the solver's own Jacobian.
"""
function build_jacobian(rhs::KernelRHS; prn = true)
    system = rhs.system
    offenders = self_feeding_instances(system)
    if !isempty(offenders)
        prn && @warn "No analytical Jacobian: $(length(offenders)) instances feed " *
            "their own output map, which the layered substitution cannot settle " *
            "in one sweep. Involved kernels: " *
            join(unique(system.kernels[system.instances[i].kernel].name
                        for i in offenders), ", ")
        return nothing
    end

    reach = output_reach(system, system.wiring, system.layers)
    owner, local_row = output_owners(system)
    count = length(system.instances)

    duals = KernelDuals[]
    dual_of = zeros(Int, count)
    position = zeros(Int, count)
    for (kernel_idx, kernel) in enumerate(system.kernels)
        batch = [i for i in 1:count if system.instances[i].kernel == kernel_idx]
        isempty(batch) && continue
        length(kernel.states) + length(kernel.inputs) == 0 && continue
        push!(duals, KernelDuals(system, kernel_idx, batch))
        for (slot, index) in enumerate(batch)
            dual_of[index] = length(duals)
            position[index] = slot
        end
    end

    gather = [instance_gather(system, owner, local_row, i, :output) for i in 1:count]
    dstate_gather = [instance_gather(system, owner, local_row, i, :dstate)
                     for i in 1:count]
    cols = [instance_columns(system, reach, i, :output) for i in 1:count]
    dstate_cols = [instance_columns(system, reach, i, :dstate) for i in 1:count]
    rows = [zeros(SimFloat, length(system.kernels[system.instances[i].kernel].outputs),
                  length(cols[i])) for i in 1:count]

    matrix = SimFloat.(system.sparsity)
    nzindex = [nonzero_indices(matrix, system.instances[i].states, dstate_cols[i])
               for i in 1:count]

    order = [i for layer in system.layers for batch in layer for i in batch]
    stateful = [i for i in 1:count
                if system.kernels[system.instances[i].kernel].rhs! !== nothing]
    widest_input = maximum(length(kernel.inputs) for kernel in system.kernels)
    widest_state = maximum(length(kernel.states) for kernel in system.kernels)
    widest_columns = maximum(max(length(cols[i]), length(dstate_cols[i]))
                             for i in 1:count)

    return KernelJacobian(rhs, duals, [entry.blocks for entry in duals],
        [(entry.n_state, entry.n_input, entry.n_output) for entry in duals],
        dual_of, position, order, stateful, gather, dstate_gather,
        cols, rows, dstate_cols, nzindex, zeros(Int, system.n_states),
        zeros(SimFloat, widest_input, widest_columns),
        zeros(SimFloat, widest_state, widest_columns), matrix,
        zeros(SimFloat, system.n_states))
end

"""
    self_feeding_instances(system) -> Vector{Int}

The instances whose own outputs feed an input their output map reads.
[`build_schedule`](@ref) skips such an edge, so the schedule does not order the
instance against itself and one sweep of `Y = A + B·G·Y` would leave that block short
of its own contribution.
"""
function self_feeding_instances(system::KernelSystem)
    owner, _ = output_owners(system)
    wiring = system.wiring
    offenders = Int[]
    for (index, instance) in enumerate(system.instances)
        feeds = system.kernels[instance.kernel].input_feeds_output
        for (slot, input) in enumerate(instance.inputs)
            feeds[slot] || continue
            any(owner[wiring.sources[j]] == index
                for j in wiring.offsets[input]:(wiring.offsets[input + 1] - 1)) &&
                (push!(offenders, index); break)
        end
    end
    return offenders
end

"""
    output_owners(system) -> (owner, local_row)

Per output slot, the instance that writes it and the slot's position within that
instance's outputs — which is the row of its `Y` block.
"""
function output_owners(system::KernelSystem)
    owner = zeros(Int, system.n_outputs)
    local_row = zeros(Int, system.n_outputs)
    for (index, instance) in enumerate(system.instances)
        for (row, slot) in enumerate(instance.outputs)
            owner[slot] = index
            local_row[slot] = row
        end
    end
    return owner, local_row
end

"""
    instance_gather(system, owner, local_row, index, which) -> InputGather

Resolve one instance's wiring into the producing instances and rows the substitution
reads, one compressed row per input slot its output maps (`which = :output`) or its
state derivatives (`which = :dstate`) read.
"""
function instance_gather(system::KernelSystem, owner, local_row, index::Int,
                         which::Symbol)
    instance = system.instances[index]
    reads = system.kernels[instance.kernel].reads
    wanted = which === :output ? reads.output_input : reads.dstate_input
    slots = sort!(unique!(reduce(vcat, wanted; init = Int[])))
    wiring = system.wiring
    offsets = Vector{Int}(undef, length(slots) + 1)
    producer, row, weight = Int[], Int[], SimFloat[]
    for (position, slot) in enumerate(slots)
        offsets[position] = length(producer) + 1
        input = instance.inputs[slot]
        for j in wiring.offsets[input]:(wiring.offsets[input + 1] - 1)
            push!(producer, owner[wiring.sources[j]])
            push!(row, local_row[wiring.sources[j]])
            push!(weight, wiring.weights[j])
        end
    end
    offsets[end] = length(producer) + 1
    return InputGather(slots, offsets, producer, row, weight)
end

"""
    instance_columns(system, reach, index, which) -> Vector{Int}

The state columns one instance's blocks span: its own states plus everything reaching
the inputs its outputs read (`which = :output`) or its state derivatives read
(`which = :dstate`). One column set per instance rather than per row keeps each block
dense, and the rows that span less of it hold structural zeros.
"""
function instance_columns(system::KernelSystem, reach, index::Int, which::Symbol)
    instance = system.instances[index]
    reads = system.kernels[instance.kernel].reads
    wanted = which === :output ? reads.output_input : reads.dstate_input
    columns = BitSet(instance.states)
    for slots in wanted, slot in slots
        input = instance.inputs[slot]
        for j in system.wiring.offsets[input]:(system.wiring.offsets[input + 1] - 1)
            union!(columns, reach[system.wiring.sources[j]])
        end
    end
    return collect(columns)
end

"""
    nonzero_indices(matrix, rows, columns) -> Matrix{Int}

Where each entry of one instance's derivative block lands in `matrix`'s value array,
`0` for an entry the pattern does not hold. Those entries are structurally zero: a
derivative reads a column only through an input whose reach the pattern already
covers, so skipping them writes nothing that is not zero.
"""
function nonzero_indices(matrix::SparseMatrixCSC, rows, columns)
    index = zeros(Int, length(rows), length(columns))
    for (column_position, column) in enumerate(columns)
        span = nzrange(matrix, column)
        for (row_position, row) in enumerate(rows)
            found = searchsortedfirst(view(rowvals(matrix), span), row)
            found > length(span) && continue
            rowvals(matrix)[span[found]] == row || continue
            index[row_position, column_position] = span[found]
        end
    end
    return index
end

"""
    (jac::KernelJacobian)(J, u, p, t)

Fill `J` with `∂du/∂u` at `(u, p, t)`. One right-hand-side call settles the input
buffer, one dual pass per kernel gives every instance's local Jacobian, a sweep in
schedule order settles `Y`, and the stateful instances then compose their rows.
"""
function (jac::KernelJacobian)(J, u, p, t)
    rhs = jac.rhs
    system = rhs.system
    rhs(jac.scratch, u, p, t)
    input = buffers(rhs, SimFloat).input

    for duals in jac.duals
        fill_blocks!(duals, system.kernels[duals.kernel], system, u, input,
                     p.numeric, p.callables[duals.kernel], t)
    end

    for index in jac.order
        settle_outputs!(jac, system, index)
    end
    fill!(nonzeros(jac.matrix), zero(SimFloat))
    for index in jac.stateful
        compose_derivatives!(jac, system, index)
    end
    store_jacobian!(J, jac.matrix)
    return nothing
end

"""
    settle_outputs!(jac, system, index)

Fill one instance's `Y` block: gather `G·Y` over the rows its outputs read, multiply
by the instance's own `B`, and add `A` on its own state columns. Called in schedule
order, so every producer it reads is already settled.
"""
function settle_outputs!(jac::KernelJacobian, system::KernelSystem, index::Int)
    slot = jac.dual_of[index]
    slot == 0 && return nothing
    n_state, _, n_output = jac.shape[slot]
    n_output == 0 && return nothing
    columns = jac.cols[index]
    target = jac.rows[index]
    position = jac.position[index]
    entries = jac.gather[index]
    gathered = gather_blocks!(jac, entries, columns)

    blocks = jac.blocks[slot]
    block_product!(target, blocks, 0, n_state, position, n_output, entries.slots,
                   gathered, length(columns))
    @inbounds for (state, column) in enumerate(system.instances[index].states)
        place = jac.column_of[column]
        for row in 1:n_output
            target[row, place] += blocks[row, state, position]
        end
    end
    clear_columns!(jac, columns)
    return nothing
end

"""
    block_product!(target, blocks, row_offset, column_offset, position, n_row, slots,
                   gathered, n_column)

`target = B · gathered`, where row `k` of `gathered` is input slot `slots[k]` and `B`
is the corresponding window of one instance's local Jacobian at
`(row_offset, column_offset)`. Written out rather than left to `mul!` because the
blocks are a handful of rows wide: a BLAS call costs more than the product, and the
loop can skip the zeros `gathered` is full of.
"""
function block_product!(target, blocks::Array{SimFloat, 3}, row_offset::Int,
                        column_offset::Int, position::Int, n_row::Int,
                        slots::Vector{Int}, gathered, n_column::Int)
    @inbounds for column in 1:n_column
        for row in 1:n_row
            target[row, column] = zero(SimFloat)
        end
        for (inner, slot) in enumerate(slots)
            weight = gathered[inner, column]
            iszero(weight) && continue
            for row in 1:n_row
                target[row, column] += blocks[row_offset + row,
                                              column_offset + slot, position] * weight
            end
        end
    end
    return nothing
end

"""
    compose_derivatives!(jac, system, index)

Add one instance's rows of `J`: `C` on its own state columns plus `D` times the
gathered `G·Y` over the inputs its derivatives read, scattered onto the assembled
matrix's nonzeros.
"""
function compose_derivatives!(jac::KernelJacobian, system::KernelSystem, index::Int)
    slot = jac.dual_of[index]
    slot == 0 && return nothing
    n_state, _, n_output = jac.shape[slot]
    columns = jac.dstate_cols[index]
    position = jac.position[index]
    entries = jac.dstate_gather[index]
    gathered = gather_blocks!(jac, entries, columns)

    blocks = jac.blocks[slot]
    target = jac.block_scratch
    block_product!(target, blocks, n_output, n_state, position, n_state,
                   entries.slots, gathered, length(columns))
    @inbounds for (state, column) in enumerate(system.instances[index].states)
        place = jac.column_of[column]
        for row in 1:n_state
            target[row, place] += blocks[n_output + row, state, position]
        end
    end

    values = nonzeros(jac.matrix)
    index_map = jac.nzindex[index]
    @inbounds for place in eachindex(columns), row in 1:n_state
        entry = index_map[row, place]
        entry == 0 || (values[entry] += target[row, place])
    end
    clear_columns!(jac, columns)
    return nothing
end

"""
    gather_blocks!(jac, entries, columns) -> Matrix

`G·Y` for one instance, in the leading `rows × columns` corner of the shared scratch,
one row per entry of `entries.slots`. Marks `columns` in `column_of` so a producer's
own column set can be mapped onto this one by lookup; a producer column this instance
does not span carries a structural zero in the row being read, so dropping it loses
nothing. The caller clears the marks.
"""
function gather_blocks!(jac::KernelJacobian, entries::InputGather,
                        columns::Vector{Int})
    gathered = jac.input_scratch
    n_row = length(entries.slots)
    @inbounds for (place, column) in enumerate(columns)
        jac.column_of[column] = place
        for slot in 1:n_row
            gathered[slot, place] = zero(SimFloat)
        end
    end
    @inbounds for slot in 1:n_row
        for edge in entries.offsets[slot]:(entries.offsets[slot + 1] - 1)
            producer = entries.producer[edge]
            source = jac.rows[producer]
            weight = entries.weight[edge]
            row = entries.row[edge]
            for (place, column) in enumerate(jac.cols[producer])
                target = jac.column_of[column]
                target == 0 && continue
                gathered[slot, target] += weight * source[row, place]
            end
        end
    end
    return gathered
end

"""Undo [`gather_blocks!`](@ref)'s marks, so `column_of` stays all zero between uses."""
function clear_columns!(jac::KernelJacobian, columns::Vector{Int})
    @inbounds for column in columns
        jac.column_of[column] = 0
    end
    return nothing
end

"""
    store_jacobian!(target, assembled)

Copy the assembled sparse Jacobian into whatever matrix the solver handed us: its
value array directly when it carries the same pattern, entry by entry otherwise.
"""
function store_jacobian!(target::SparseMatrixCSC, assembled::SparseMatrixCSC)
    if target.colptr == assembled.colptr && target.rowval == assembled.rowval
        copyto!(nonzeros(target), nonzeros(assembled))
        return nothing
    end
    fill!(nonzeros(target), zero(eltype(target)))
    scatter_jacobian!(target, assembled)
    return nothing
end

function store_jacobian!(target::AbstractMatrix, assembled::SparseMatrixCSC)
    fill!(target, zero(eltype(target)))
    scatter_jacobian!(target, assembled)
    return nothing
end

"""Write every stored entry of `assembled` into `target` by index."""
function scatter_jacobian!(target, assembled::SparseMatrixCSC)
    values = nonzeros(assembled)
    rows = rowvals(assembled)
    for column in axes(assembled, 2), entry in nzrange(assembled, column)
        target[rows[entry], column] = values[entry]
    end
    return nothing
end
