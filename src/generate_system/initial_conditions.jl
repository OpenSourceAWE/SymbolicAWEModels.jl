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
the initial condition for `state_var` (a scalar state term, a collected vector,
e.g. `pos[:, i]`, or a matrix, e.g. `R_b_to_w[:, :, i]` — arrays are flattened
column-major and must align element-wise with the read value). Returns the constant
default pair(s) (build-time numeric value) to splice into the system `defaults`,
which makes MTK expose a settable `Initial(state_var)`.
"""
function bind_initial!(path_view::InitialPath, state_var)
    reg = getfield(path_view, :reg)
    path = getfield(path_view, :path)
    value = read_path(reg.sys_struct, path)
    vars = state_var isa AbstractArray ? vec(collect(state_var)) : Any[state_var]
    vals = value isa AbstractArray ? vec(collect(value)) : Any[value]
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

Push bound initial conditions from the struct onto a problem. A bound variable
reaches the problem one of two ways: through its `Initial(state)` parameter when
that parameter exists (an initialization constraint mtkcompile kept), or — when
there is no such parameter (`build_initializeprob=false`) — by writing the
surviving unknown's `u0` directly. Holds a `setp`/`setu` setter, per-element
readers, and a preallocated value buffer for each path.
"""
struct InitialSync{PSet, USet}
    param_setter::PSet
    param_readers::Vector{Any}
    param_buffer::Vector{SimFloat}
    state_setter::USet
    state_readers::Vector{Any}
    state_buffer::Vector{SimFloat}
end

"""
    build_initial_sync(sys, registry) -> InitialSync | Nothing

Build the initial-condition sync from the compiled system and the registry. Each
bound variable is routed by how it survived `mtkcompile`: to its `Initial(var)`
parameter when that is a real parameter of `sys` (a surviving observed variable
or an init constraint), otherwise — when the variable survives as an unknown but
carries no `Initial` parameter — to a direct `u0` write. Variables removed
entirely are dropped.
"""
function build_initial_sync(sys, registry::InitialRegistry)
    isempty(registry.entries) && return nothing
    init_params, param_readers = Any[], Any[]
    state_vars, state_readers = Any[], Any[]
    for entry in registry.entries
        for k in eachindex(entry.vars)
            var = entry.vars[k]
            init = ModelingToolkit.Initial(var)
            # A surviving unknown's u0 must be written directly: with
            # build_initializeprob=false there is no init solve to apply its
            # Initial parameter, so setting that parameter alone leaves u0 stale.
            if is_variable(sys, var)
                push!(state_vars, var)
                push!(state_readers, ElementReader(entry.read, k))
            elseif is_parameter(sys, init)
                push!(init_params, init)
                push!(param_readers, ElementReader(entry.read, k))
            end
        end
    end
    (isempty(init_params) && isempty(state_vars)) && return nothing
    return InitialSync(
        isempty(init_params) ? nothing : setp(sys, init_params),
        param_readers, Vector{SimFloat}(undef, length(init_params)),
        isempty(state_vars) ? nothing : setu(sys, state_vars),
        state_readers, Vector{SimFloat}(undef, length(state_vars)))
end

"""
    sync_initial!(sync, prob, sys_struct)

Copy every bound initial condition from the live `sys_struct` onto `prob` — both
the `Initial` parameters and the directly-set unknown `u0` values. Must run
before a fresh `init`/`solve` (a `reinit!` with `reinit_dae` does not re-read
them). A no-op when there are none.
"""
sync_initial!(::Nothing, prob, sys_struct) = nothing
function sync_initial!(sync::InitialSync, prob, sys_struct::SystemStructure)
    if sync.param_setter !== nothing
        @inbounds for k in eachindex(sync.param_readers)
            sync.param_buffer[k] = sync.param_readers[k](sys_struct)
        end
        sync.param_setter(prob, sync.param_buffer)
    end
    if sync.state_setter !== nothing
        @inbounds for k in eachindex(sync.state_readers)
            sync.state_buffer[k] = sync.state_readers[k](sys_struct)
        end
        sync.state_setter(prob, sync.state_buffer)
    end
    return nothing
end
