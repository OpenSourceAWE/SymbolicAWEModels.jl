# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Initial conditions as MTK `Initial(x)` parameters (mirror of flat_params.jl).

"""
    InitialEntry

One initial-condition binding: the scalar state-variable terms `vars` (e.g.
`[pos[1,i], pos[2,i], pos[3,i]]`) and a `read(sys_struct)` returning their live
value (scalar or vector, element-aligned with `vars`).
"""
struct InitialEntry
    vars::Vector{Any}
    read::Any
end

"""
    InitialRegistry

Build-time record of every `bind_initial!` call. Transient: it lives only during
equation generation and drives `build_initial_sync` after compilation.
"""
mutable struct InitialRegistry
    sys_struct::SystemStructure
    entries::Vector{InitialEntry}
end
InitialRegistry(sys_struct::SystemStructure) =
    InitialRegistry(sys_struct, InitialEntry[])

# ==================== BUILD-TIME VIEW ==================== #

"""
Top-level `initial` view wrapping an [`InitialRegistry`](@ref).
`initial.points[i].pos_w` mirrors `sys_struct.points[i].pos_w`.
"""
struct InitialView
    reg::InitialRegistry
end

"""A partial path into `sys_struct` identifying an initial-condition field."""
struct InitialPath
    reg::InitialRegistry
    path::Tuple
end

Base.getproperty(view::InitialView, sym::Symbol) =
    sym === :reg ? getfield(view, :reg) :
    InitialPath(getfield(view, :reg), (sym,))

Base.getindex(view::InitialPath, idx::Integer) =
    InitialPath(getfield(view, :reg),
                (getfield(view, :path)..., Int(idx)))

Base.getproperty(view::InitialPath, sym::Symbol) =
    (sym === :reg || sym === :path) ? getfield(view, sym) :
    InitialPath(getfield(view, :reg),
                (getfield(view, :path)..., sym))

"""
    bind_initial!(initial_path, state_var) -> Vector{Pair}

Record that struct field `initial_path` (e.g. `initial.points[i].pos_w`) provides
the initial condition for `state_var` (a scalar state term or a collected vector
of them, e.g. `pos[:, i]`). Returns the constant default pair(s) (build-time
numeric value) to splice into the system `defaults`, which makes MTK expose a
settable `Initial(state_var)`.
"""
function bind_initial!(path_view::InitialPath, state_var)
    reg = getfield(path_view, :reg)
    path = getfield(path_view, :path)
    value = read_path(reg.sys_struct, path)
    vars = state_var isa AbstractVector ? collect(state_var) : Any[state_var]
    vals = value isa AbstractVector ? collect(value) : Any[value]
    length(vars) == length(vals) || error(
        "bind_initial! length mismatch at $path: " *
        "$(length(vars)) state vars vs $(length(vals)) values.")
    push!(reg.entries, InitialEntry(Vector{Any}(vars), PathReader(path)))
    return [vars[k] => vals[k] for k in eachindex(vars)]
end

# ==================== SYNC ==================== #

"""
    ElementReader(base, index)

Serialisable reader returning the `index`-th element of `base(sys_struct)`. Used
to map an array-valued struct field onto per-element `Initial` parameters
(scalars are indexed at 1, which is a no-op).
"""
struct ElementReader{R}
    base::R
    index::Int
end
(reader::ElementReader)(sys_struct) = reader.base(sys_struct)[reader.index]

"""
    InitialSync

A `setp` setter for a list of `Initial(state)` parameters, their per-element
readers, and a preallocated value buffer.
"""
struct InitialSync{Setter}
    setter::Setter
    readers::Vector{Any}
    buffer::Vector{SimFloat}
end

"""
    build_initial_sync(sys, registry) -> InitialSync | Nothing

Build the `Initial`-parameter sync from the compiled system and the registry. A
bound variable is kept when `Initial(var)` is a real parameter of `sys` — this
holds for surviving unknowns *and* for observed variables `mtkcompile` solves
during initialization (their `Initial` is an initialization constraint); only
variables removed entirely are dropped, so the setter never touches an absent
`Initial` parameter.
"""
function build_initial_sync(sys, registry::InitialRegistry)
    isempty(registry.entries) && return nothing
    init_params, readers = Any[], Any[]
    for entry in registry.entries
        for k in eachindex(entry.vars)
            init = ModelingToolkit.Initial(entry.vars[k])
            is_parameter(sys, init) || continue
            push!(init_params, init)
            push!(readers, ElementReader(entry.read, k))
        end
    end
    isempty(init_params) && return nothing
    return InitialSync(setp(sys, init_params), readers,
                       Vector{SimFloat}(undef, length(init_params)))
end

"""
    sync_initial!(sync, prob, sys_struct)

Copy every bound initial condition from the live `sys_struct` onto `prob`'s
`Initial` parameters. Must run before a fresh `init`/`solve` (a `reinit!` with
`reinit_dae` does not re-read `Initial`). A no-op when there are none.
"""
sync_initial!(::Nothing, prob, sys_struct) = nothing
function sync_initial!(sync::InitialSync, prob, sys_struct::SystemStructure)
    readers = sync.readers
    buffer = sync.buffer
    @inbounds for k in eachindex(readers)
        buffer[k] = readers[k](sys_struct)
    end
    sync.setter(prob, buffer)
    return nothing
end
