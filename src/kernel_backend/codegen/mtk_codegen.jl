# SPDX-FileCopyrightText: 2019-2025 Frank Hellmann, Michael Lindner and Hans Würfel and Contributors
# SPDX-License-Identifier: MIT
#
# Vendored from NetworkDynamics.jl, branch `nonnumeric-params`, commit
# 68756bcebb0d98d591f2e0965cb84e08b898c716 (2026-08-06). See LICENSES/MIT.txt.
#
# Upstream sources, per section below:
#   ext/MTKExt_utils.jl          eq_type, rhs_differentials, getproperty_symbolic,
#                                generate_massmatrix, check_metadata, fix_metadata!,
#                                match_diff_states, get_variables_deriv
#   ext/MTKExt_simplification.jl _scalarize_eqs, _scalarize_system, _insert_sorted!,
#                                get_alias, _alias_connected_components,
#                                pick_best_alias_names, simplify_with_mtkcompile
#   ext/NetworkDynamicsMTKExt.jl _scalarize_io_syms, is_nonnumeric_param,
#                                _get_formulas, _collect_deps_on_obs,
#                                _all_dependencies, generate_io_function
#
# Removed while vendoring — all of it exists only to serve NetworkDynamics' fixed
# depth-2 coreloop, which this backend does not have:
#   * `ff_to_constraint` promotion of feedforward outputs to algebraic states
#   * `fftype` as a scheduling concept; `g` here always takes the same arguments
#   * `assume_io_coupling` and its `implicit_output` machinery
#   * `simplify_without_mtkcompile` and its equation reducer (Hungarian matching)
#   * the ComponentModel / Network / VIndex surface, component metadata,
#     init/guess formulas, and the model cache
#
# Deliberate deviations from upstream, to avoid two dependencies for two call
# sites: `@argcheck` is `@assert`, and Moshi's `@match` in `get_variables_deriv`
# is an `iscall`/`operation` test. `cse` defaults to `true` here (measured 167x
# compile and 3.8-9.8x runtime on these kernels, agreeing to <=8e-15).
# `fix_metadata!` reports its residue only under `verbose`: a scalarized array
# element never resolves by name, so upstream warns on every kernel we build.
# `drop_unused_inputs` replaces upstream's element-wise drop of inputs that no
# equation reads, and `mtkcompile_inputs` passes a partially-read array whole,
# because `mtkcompile` accepts neither part of an array nor an absent input.

"""
    KernelCodegen

Compile one ModelingToolkit `System` with declared inputs and outputs into plain
callable functions. `generate_io_function` is the only entry point; everything
else is upstream NetworkDynamics support code kept close to its original form so
it stays diffable against a fresh checkout.
"""
module KernelCodegen

using ModelingToolkit: ModelingToolkitBase
using ModelingToolkit.ModelingToolkitBase: Equation, System, Differential,
    equations, full_equations, get_variables, mtkcompile, getname, unwrap,
    parameters, unknowns, independent_variables, observed, iscall, operation,
    arguments, build_function
using Symbolics: Symbolics, fixpoint_sub
using SymbolicUtils: SymbolicUtils
using SymbolicUtils.Code: Let, Assignment, unwrap_const
using OrderedCollections: OrderedDict
using LinearAlgebra: Diagonal, I

export generate_io_function, scalar_name

"The concrete `BasicSymbolic` type every symbolic in these pipelines has."
const ST = SymbolicUtils.BasicSymbolicImpl.var"typeof(BasicSymbolicImpl)"{SymbolicUtils.SymReal}

# ======================= MTKExt_utils.jl ======================= #

"""
    eq_type(eq::Equation)

Checks the type of the equation. Returns:
- `(:explicit_diffeq, lhs_variable)` for explicit differential equations
- `(:implicit_diffeq, nothing)` for implicit differential equations
- `(:explicit_algebraic, lhs_variable)` for explicit algebraic equations
- `(:implicit_algebraic, nothing)` for implicit algebraic equations
"""
function eq_type(eq::Equation)
    rhs_differentials = _collect_differentials(eq.rhs)
    if !isempty(rhs_differentials)
        return (:implicit_diffeq, nothing)
    end

    if isdifferential(eq.lhs)
        return (:explicit_diffeq, only(eq.lhs.args))
    end

    lhsvars = get_variables(eq.lhs)
    rhsvars = get_variables(eq.rhs)
    commonvars = intersect(lhsvars, rhsvars)
    if isempty(commonvars) && length(lhsvars) == 1 && isequal(eq.lhs, only(lhsvars))
        return (:explicit_algebraic, eq.lhs)
    else
        return (:implicit_algebraic, nothing)
    end
end

isdifferential(s) = iscall(unwrap(s)) && operation(unwrap(s)) isa Differential

"The set of `D(x)` terms appearing on the right-hand side of `eqs`."
function rhs_differentials(eqs::Vector{Equation})
    diffs = Set{ST}()
    for eq in eqs
        _collect_differentials!(diffs, eq.rhs)
    end
    return diffs
end

_collect_differentials(ex) = _collect_differentials!(Set{ST}(), ex)

function _collect_differentials!(found, ex)
    if iscall(ex)
        if operation(ex) isa Differential
            push!(found, ex)
        else
            for arg in arguments(ex)
                _collect_differentials!(found, arg)
            end
        end
    end
    return found
end

"""
    getproperty_symbolic(sys, var; might_contain_toplevel_ns=true)

Like `getproperty` but works on a greater variety of "var"
- var can be Num or Symbolic (resolved using getname)
- strip namespace of sys if present (don't strip if `might_contain_top_level_ns=false`)
- for nested variables (foo₊bar₊baz) resolve them one by one
"""
function getproperty_symbolic(sys, var; might_contain_toplevel_ns=true)
    ns = string(getname(sys))
    varname = string(getname(var))
    if might_contain_toplevel_ns && startswith(varname, ns*"₊")
        if getname(sys) ∈ getname.(ModelingToolkitBase.get_systems(sys))
            @warn "Namespace :$ns appears multiple times, this might lead to unexpected, since it is not clear whether the namespace should be stripped or not."
        end
        varname = replace(varname, r"^"*ns*"₊" => "")
    end
    parts = split(varname, "₊")
    # Descend a `₊`-segment ONLY when it names an actual child subsystem of the
    # current node; the first segment that is not a subsystem starts a flat variable
    # name, so the remaining segments (re-joined by `₊`) are resolved as one variable.
    r = sys
    i = 1
    while i <= length(parts)
        first = (i == 1)
        if r isa System && Symbol(parts[i]) in getname.(ModelingToolkitBase.get_systems(r))
            r = getproperty(r, Symbol(parts[i]); namespace = !first)
            i += 1
        else
            r = getproperty(r, Symbol(join(parts[i:end], "₊")); namespace = !first)
            break
        end
    end
    unwrap(r)
end

"Diagonal mass matrix: 1 for each explicit differential equation, 0 for each
implicit algebraic one. Returns `I` when every equation is differential."
function generate_massmatrix(eqs::AbstractVector{Equation})
    V = map(eqs) do eq
        type = eq_type(eq)[1]
        if type === :explicit_diffeq
            1
        elseif type === :implicit_algebraic
            0
        else
            error("Cant build mass matrix entry for $(eq) of type $type")
        end
    end
    M = Diagonal(V)
    return M==I ? I : M
end

"The variables in `exprs` that lost their symbolic metadata during simplification."
function check_metadata(exprs)
    nometadata = []
    for ex in exprs
        if ex isa Equation
            _check_metadata!(nometadata, ex.rhs)
            _check_metadata!(nometadata, ex.lhs)
        else
            _check_metadata!(nometadata, ex)
        end
    end
    return unique!(nometadata)
end

function _check_metadata!(nometadata, expr)
    vars = get_variables_deriv(expr)
    for v in vars
        isnothing(Symbolics.metadata(v)) && push!(nometadata, v)
    end
end

"""
    _match_scalar_element(invalids, valid)

When `invalids` is a scalar array element (`pos[1]`) whose metadata-carrying
resolution `valid` came back as the whole array, return the matching scalarized
element of `valid` (`Symbolics.scalarize(valid)[i]`); otherwise return `valid`.
"""
function _match_scalar_element(invalids, valid)
    (valid === nothing) && return valid
    (Symbolics.symtype(valid) <: AbstractArray &&
     !(Symbolics.symtype(Symbolics.unwrap(invalids)) <: AbstractArray)) || return valid
    u = Symbolics.unwrap(invalids)
    (Symbolics.iscall(u) && Symbolics.operation(u) === getindex) || return valid
    idx = map(i -> i isa Integer ? Int(i) : Int(Symbolics.value(i)),
              Symbolics.arguments(u)[2:end])
    return collect(Symbolics.scalarize(valid))[idx...]
end

"Substitute metadata-carrying symbolics back into `invalid_eqs` in place, resolving
each stripped variable by name against `sys`. Scalarized array elements never
resolve, so the residual report is behind `verbose` rather than a warning."
function fix_metadata!(invalid_eqs, sys; verbose = false)
    missingmetadata = check_metadata(invalid_eqs)
    if isempty(missingmetadata)
        return invalid_eqs
    end

    metadatasubs = Dict()
    allsyms = ModelingToolkitBase.all_symbols(sys)
    filter!(s->!contains(repr(s), "Initial"), allsyms)
    allnames = string.(ModelingToolkitBase.getname.(allsyms))

    for invalids in missingmetadata
        invalidname = getname(invalids)
        valid = if hasproperty(sys, getname(invalidname))
            getproperty_symbolic(sys, invalids)
        else
            idxs = findall(contains(string(invalidname)), allnames)
            if length(idxs) == 1
                allsyms[only(idxs)]
            else
                verbose && @warn "Could not resolve invalid symbol $invalidname, options are $(allsyms[idxs])"
            end
        end
        metadatasubs[invalids] = _match_scalar_element(invalids, valid)
    end

    fixedeqs = [Symbolics.substitute(eq, metadatasubs) for eq in invalid_eqs]
    if verbose && !isempty(check_metadata(fixedeqs))
        @warn "Some transformation dropped metadata ($missingmetadata)! Could NOT be fixed. $(check_metadata(fixedeqs))"
    end
    invalid_eqs .= fixedeqs
end

"""
    match_diff_states(eqs, states)

Match a set of unknowns to a vector of equations. Returns a vector of states with
diff states at the correct places and non-diff states filled in arbitrarily.
"""
function match_diff_states(eqs, states)
    @assert length(eqs) == length(states) "Number of equations and states must match to be able to match them, got $(length(eqs)) equations and $(length(states)) states"
    states_new = Vector{eltype(states)}(undef, length(eqs))
    diff_state_set = Set(states)
    for (i, eq) in pairs(eqs)
        type, var = eq_type(eq)
        if type == :explicit_diffeq
            states_new[i] = pop!(diff_state_set, var)
        elseif type == :implicit_algebraic
            # do nothing
        else
            error("Got equation $eq of type $type, which is not expected at that stage (only explicit diff eqs and implicit alg. eqs). Something went wrong.")
        end
    end
    # Fill algebraic slots in the order states appear in the input vector (not hash order).
    alg_states = [s for s in states if s ∈ diff_state_set]
    ai = 0
    for i in eachindex(states_new)
        if !isassigned(states_new, i)
            ai += 1
            states_new[i] = alg_states[ai]
        end
    end
    @assert ai == length(alg_states) "Not all states could be matched, something went wrong"
    @assert length(eqs) == length(states_new)
    states_new
end

"Like `get_variables` but replaces any `D(x)` term with the inner variable `x`,
which `get_variables` intentionally treats as one atomic symbol."
function get_variables_deriv(ex)
    set = get_variables(ex)
    for s in set
        u = unwrap(s)
        if iscall(u) && operation(u) isa Differential
            pop!(set, s)
            push!(set, only(arguments(u)))
        end
    end
    set
end

# ======================= MTKExt_simplification.jl ======================= #

"""
    _scalarize_eqs(eqs) -> Vector{Equation}

Expand array-valued equations (`lhs::Array ~ rhs::Array`) into their scalar
component equations; scalar equations pass through unchanged. Used so array I/O
variables can be scalarized without vector equations leaking into the scalar-only
codegen pipeline.
"""
function _scalarize_eqs(eqs)
    out = Equation[]
    for eq in eqs
        sc = Symbolics.scalarize(eq)
        for e in (sc isa AbstractArray ? sc : (sc,))
            # drop trivial `x ~ x` identities (e.g. from scalarizing an array
            # reconstruction `v ~ [v[1],…]`), which would self-loop the obs walk
            isequal(e.lhs, e.rhs) && continue
            push!(out, e)
        end
    end
    return out
end

_isarrsym(v) = Symbolics.symtype(Symbolics.unwrap(v)) <: AbstractArray

# `vec` flattens a scalarized *matrix* symbol (e.g. a 3×N parameter) to a column-major
# vector; a scalarized vector is already 1-D, so `vec` is a no-op there.
_scal_flat(vs) = reduce(vcat,
    (_isarrsym(v) ? vec(collect(Symbolics.scalarize(v))) : [v] for v in vs); init=Any[])

"""
    _scalarize_system(sys) -> System

Rebuild `sys` with every array-valued unknown/parameter/equation expanded to its
scalar components, so the scalar-only codegen sees one consistent set of scalar
symbolics (avoids array/scalar identity mismatches between our scalarization and
mtkcompile's). A no-op when the system already has no array symbolics.
"""
function _scalarize_system(sys)
    any(_isarrsym, vcat(unknowns(sys), parameters(sys))) || return sys
    eqs = _scalarize_eqs(equations(sys))
    obs = _scalarize_eqs(observed(sys))
    unk = identity.(_scal_flat(unknowns(sys)))
    ps = identity.(_scal_flat(parameters(sys)))
    # Keep array-symbol default keys atomic: `AtomicArrayDict` (the System's
    # initial-conditions store) rejects indexed keys like `world_damping[3]`, so a
    # scalarized array parameter/unknown keeps its default on the whole-array key.
    defs = Dict{Any, Any}()
    for (k, v) in ModelingToolkitBase.get_initial_conditions(sys)
        defs[k] = v
    end
    iv = only(independent_variables(sys))
    return System(eqs, iv, unk, ps; observed=obs, initial_conditions=defs,
                  name=nameof(sys))
end

"Insert `newobs` (already in topological order) into `obseqs`, each just after the
last equation it depends on."
function _insert_sorted!(obseqs, newobs)
    for eq in newobs
        _insert_sorted!(obseqs, eq)
    end
    nothing
end

function _insert_sorted!(obseqs, eq::Equation)
    obssym = [e.lhs for e in obseqs]
    vars = get_variables_deriv(eq.rhs)
    idx = findlast(sym -> sym ∈ vars, obssym)
    last_dependency = isnothing(idx) ? 0 : idx
    insert!(obseqs, last_dependency+1, eq)
end

"The `(a, b)` pair of an equation that is a pure alias `a ~ b`, else `nothing`."
function get_alias(eq)
    vars = get_variables_deriv(eq)
    length(vars) == 2 || return nothing
    a, b = vars
    # parameters are not aliases — an equation like `x ~ V` (V a parameter) defines x,
    # it does not alias two unknowns
    any(ModelingToolkitBase.isparameter, (a, b)) && return nothing
    match = if isequal(unwrap_const(eq.lhs), 0)
        isequal(eq.rhs, a-b) || isequal(eq.rhs, b-a)
    else
        isequal(eq, a ~ b) || isequal(eq, b ~ a)
    end
    if match
        return (a, b)
    else
        return nothing
    end
end

"""
    _alias_connected_components(pairs)

Given alias pairs `[(a,b), (c,d), ...]`, returns groups of transitively connected
variables as a `Vector{Vector}`. E.g. pairs `(a,b), (b,c)` → `[[a, b, c]]`.
"""
function _alias_connected_components(pairs)
    adj = Dict{ST, Vector{ST}}()
    for (a, b) in pairs
        push!(get!(adj, a, ST[]), b)
        push!(get!(adj, b, ST[]), a)
    end
    visited = Set{ST}()
    groups = Vector{Vector{ST}}()
    for v in keys(adj)
        v ∈ visited && continue
        group = ST[]
        stack = ST[v]
        while !isempty(stack)
            u = pop!(stack)
            u ∈ visited && continue
            push!(visited, u)
            push!(group, u)
            for w in adj[u]
                w ∈ visited || push!(stack, w)
            end
        end
        push!(groups, group)
    end
    return groups
end

"""
    pick_best_alias_names(eqs, obseqs, states, outputs, inputs; verbose)

Consolidate the alias chains simplification leaves behind. Pure alias equations
`a ~ b` in `obseqs` are grouped into transitively connected clusters; one *main*
representative is picked per cluster (differential states first, then inputs, then
outputs, then the least deeply nested name); the substitution `nonmain → main` is
applied throughout `eqs`, `obseqs` and `states`, and canonical `nonmain ~ main`
observations are reinserted.
"""
function pick_best_alias_names(eqs, obseqs, states, outputs, inputs; verbose)
    alias_pairs = Tuple{ST,ST}[]
    alias_idx = Int[]
    for (i, eq) in enumerate(obseqs)
        alias = get_alias(eq)
        isnothing(alias) && continue
        push!(alias_pairs, alias)
        push!(alias_idx, i)
    end

    groups = _alias_connected_components(alias_pairs)

    diffstateset = Set{ST}()
    for eq in eqs
        type, var = eq_type(eq)
        type == :explicit_diffeq && push!(diffstateset, var)
    end
    inset = Set{ST}(inputs)
    outset = Set{ST}(outputs)
    sortf = function(s)
        prio = if s ∈ diffstateset
            0
        elseif s ∈ inset
            1
        elseif s ∈ outset
            2
        else
            3
        end
        # name as final tiebreak: when a class is structurally symmetric under renaming
        # nothing above distinguishes the members, so fall back to the name to keep the
        # representative deterministic across Julia versions.
        (prio, count('₊', string(s)), string(getname(s)))
    end
    alias_subs = Dict()
    new_alias_obs = Equation[]
    for group in groups
        sorted = sort!(collect(group), by=sortf)
        main = first(sorted)
        for r in @view(sorted[2:end])
            alias_subs[r] = main
            push!(new_alias_obs, r ~ main)
        end
    end
    verbose && @info "Found $(length(groups)) alias groups"

    eqs_new = Symbolics.substitute.(eqs, Ref(alias_subs))

    obseqs_new = let
        obseqs_new = copy(obseqs)
        obseqs_new = deleteat!(obseqs_new, alias_idx)
        obseqs_new = Symbolics.substitute.(obseqs_new, Ref(alias_subs))
        _insert_sorted!(obseqs_new, new_alias_obs)
        obseqs_new
    end

    states_new = Symbolics.substitute.(states, Ref(alias_subs))
    allunique(states_new) || error("Alias elimination resulted in duplicate state names: $(repr(states_new))! This should never happen.")

    eqs_new, obseqs_new, states_new
end

"""
    io_array_base(sym)

The array a scalarized element belongs to (`pos[2]` → `pos`), or the symbol itself
when it is not an array element.
"""
function io_array_base(sym)
    u = unwrap(sym)
    (iscall(u) && operation(u) === getindex) || return u
    return arguments(u)[1]
end

"""
    drop_unused_inputs(inputs, used) -> Set

The inputs to hide from `mtkcompile` because no equation reads them. An array is
all-in or all-out here: an array *some* of whose components are read stays, and
[`mtkcompile_inputs`](@ref) then passes it whole.
"""
function drop_unused_inputs(inputs, used)
    unused = setdiff(inputs, used)
    kept = Set(io_array_base(v) for v in setdiff(inputs, unused))
    return Set(v for v in unused if io_array_base(v) ∉ kept)
end

"""
    mtkcompile_inputs(inputs, used) -> Vector

The input list `mtkcompile` accepts, given the inputs left after
[`drop_unused_inputs`](@ref) and the set of variables the equations actually read.
`mtkcompile` refuses *part* of a declared array — "the entire array must be an
input" — and equally refuses an input that appears in no equation, so an array only
some of whose components are read is passed as the whole array symbolic and
everything else element by element. Our components read parts of vectors routinely:
a wrench reads `pos[3]` for the air density at its height, and a flap hinge about
`[0, 1, 0]` leaves six of its nine frame entries unread.
"""
function mtkcompile_inputs(inputs, used)
    split_bases = Set(io_array_base(sym) for sym in inputs if sym ∉ used)
    result = Any[]
    for sym in inputs
        base = io_array_base(sym)
        if !isequal(base, sym) && base ∈ split_bases
            any(isequal(base), result) || push!(result, base)
        else
            push!(result, sym)
        end
    end
    return result
end

"""
    simplify_with_mtkcompile(sys, allinputs, alloutputs; verbose)

Run `mtkcompile` with the declared inputs left unbound and the declared outputs
kept, then scalarize what it left as array symbolics and repair the metadata it
stripped. Returns `(sys, eqs, obseqs_sorted, states, params)` with `states`
ordered to match `eqs`.
"""
function simplify_with_mtkcompile(_sys, allinputs, alloutputs; verbose)
    missing_inputs = Set{ST}()
    sys = if ModelingToolkitBase.iscomplete(_sys)
        _sys
    else
        _openinputs = setdiff(allinputs, Set(parameters(_sys)))
        all_eq_vars = mapreduce(get_variables_deriv, union, full_equations(_sys), init=Set{ST}())
        if !(_openinputs ⊆ all_eq_vars)
            missing_inputs = drop_unused_inputs(_openinputs, all_eq_vars)
            verbose && @warn "The specified inputs ($missing_inputs) do not appear in the equations of the system!"
            _openinputs = setdiff(_openinputs, missing_inputs)
        end

        implicit_outputs = setdiff(alloutputs, all_eq_vars)
        if !isempty(implicit_outputs)
            throw(ArgumentError("The outputs $(getname.(implicit_outputs)) do not \
                appear in the equations of the system!"))
        end

        verbose && @info "Simplifying system with $(length(_openinputs)) inputs and $(length(alloutputs)) outputs"
        mtkcompile(_sys; inputs=mtkcompile_inputs(_openinputs, all_eq_vars),
                   outputs=alloutputs, simplify=false)
    end

    # extract the main equations and observed equations, scalarizing any array
    # equations mtkcompile left behind (e.g. observed `pos ~ [pos[1],…]`) so the
    # scalar-only downstream (fix_metadata!, codegen) never sees vector symbolics.
    eqs::Vector{Equation} = _scalarize_eqs(equations(sys))
    obseqs_sorted::Vector{Equation} = _scalarize_eqs(observed(sys))
    fix_metadata!(eqs, sys; verbose)
    fix_metadata!(obseqs_sorted, sys; verbose)

    # mtkcompile represents array-origin params/unknowns as the whole array;
    # scalarize them so they line up with the scalarized inputs/outputs/eqs.
    allparams = _scal_flat(parameters(sys)) # contains inputs!
    # mtkcompile/complete calls remove_bound_parameters_from_ps which removes params
    # lets push those bindings to the obseqs
    bps = ModelingToolkitBase.bound_parameters(sys)
    if !isempty(bps)
        bindings = ModelingToolkitBase.bindings(sys)
        newobs = [bp ~ bindings[bp] for bp in bps]
        _insert_sorted!(obseqs_sorted, newobs)
    end

    @assert allinputs ⊆ Set(allparams) ∪ missing_inputs
    params = setdiff(allparams, Set(allinputs))

    states = match_diff_states(eqs, _scal_flat(unknowns(sys)))

    sys, eqs, obseqs_sorted, states, params
end

# ======================= NetworkDynamicsMTKExt.jl ======================= #

"""
    _scalarize_io_syms(syms) -> Vector

Expand any array-valued symbolics in `syms` into their scalar components, leaving
scalars untouched, so array I/O variables (`pos(t)[1:3]`) become `pos[1], pos[2],
pos[3]` for the otherwise scalar-only codegen pipeline.
"""
function _scalarize_io_syms(syms)
    out = Any[]
    for s in syms
        u = Symbolics.unwrap(s)
        if Symbolics.symtype(u) <: AbstractArray
            append!(out, Symbolics.unwrap.(collect(Symbolics.scalarize(u))))
        else
            push!(out, u)
        end
    end
    return out
end

"""
    is_nonnumeric_param(p) -> Bool

Whether a parameter symbolic is nonnumeric, i.e. a callable operator (an
interpolant/polar declared as `(name::T)(..)`). Its `symtype` is a
`SymbolicUtils.FnType`. Numeric scalars and scalarized array elements have a
`Number` symtype and stay in the flat parameter buffer.
"""
is_nonnumeric_param(p) = Symbolics.symtype(p) <: SymbolicUtils.FnType

_all_rhs_symbols(term) = get_variables_deriv(term)
_all_rhs_symbols(eq::Equation) = get_variables_deriv(eq.rhs)
_all_rhs_symbols(eqs::Union{AbstractVector,AbstractDict}) =
    mapreduce(eq->get_variables_deriv(eq isa Pair ? eq.second : eq.rhs), ∪, eqs,
              init=Set{ST}())

"The symbols `term` depends on once every substitution in `dict` is resolved
recursively."
function _all_dependencies(term, dict)
    deps = Set{ST}()
    _recursive_collect_dependencies!(deps, term, dict)
    deps
end

function _recursive_collect_dependencies!(deps, term, dict)
    termdeps = _all_rhs_symbols(term)
    for sym in termdeps
        if haskey(dict, sym)
            _recursive_collect_dependencies!(deps, dict[sym], dict)
        else
            push!(deps, sym)
        end
    end
    deps
end

"""
    dependency_indices(eqs, obs_subs, symbols) -> Vector{Vector{Int}}

For each equation, the positions in `symbols` its right-hand side reads once every
observed substitution is resolved. `obs_subs` is in definition-before-use order, so
each observed variable is resolved once and reused, where `_all_dependencies`
re-descends the whole chain per call.
"""
function dependency_indices(eqs, obs_subs, symbols)
    position = Dict{ST, Int}(sym => k for (k, sym) in enumerate(symbols))
    resolved = Dict{ST, Set{Int}}()
    for (lhs, rhs) in obs_subs
        resolved[lhs] = _resolved_reads(rhs, position, resolved)
    end
    return [sort!(collect(_resolved_reads(eq.rhs, position, resolved))) for eq in eqs]
end

"The positions `term` reads, taking an already-resolved observed variable's own."
function _resolved_reads(term, position, resolved)
    reads = Set{Int}()
    for sym in get_variables_deriv(term)
        known = get(resolved, sym, nothing)
        if known === nothing
            index = get(position, sym, nothing)
            index === nothing || push!(reads, index)
        else
            union!(reads, known)
        end
    end
    return reads
end

"Build the `build_function` argument list for `eqs`: one `Let` block holding every
observed assignment `eqs` needs plus the equation assignments, returning the first
output, then the remaining outputs."
function _get_formulas(eqs, obs_subs)
    isempty(eqs) && return []
    obsdeps = _collect_deps_on_obs([eq.rhs for eq in eqs], obs_subs)
    obs_assignments = [Assignment(k, v) for (k,v) in obs_subs if k ∈ obsdeps]

    # implicit equations do not become assignments (lhs is 0), so we use the rhs directly
    eqs_assignments = [Assignment(eq.lhs, eq.rhs) for eq in eqs
                          if !isequal(eq.lhs, eq.rhs) && !isequal(unwrap_const(eq.lhs), 0)]
    out = [isequal(unwrap_const(eq.lhs), 0) ? eq.rhs : eq.lhs for eq in eqs]

    [Let(vcat(obs_assignments, eqs_assignments), out[1], false), out[2:end]...]
end

function _collect_deps_on_obs(terms, obs_subs)
    deps = Set{ST}()
    for term in terms
        _collect_deps_on_obs!(deps, obs_subs, term)
    end
    deps
end

function _collect_deps_on_obs!(deps, obs_subs, term)
    termdeps = get_variables_deriv(term)
    for sym in termdeps
        if haskey(obs_subs, sym)
            _collect_deps_on_obs!(deps, obs_subs, obs_subs[sym])
            push!(deps, sym)
        end
    end
    deps
end

"""
    scalar_name(sym) -> Symbol

Unique label for a symbolic. A scalarized array element `pos(t)[2]`, whose `getname`
collapses to `:pos`, becomes `:pos_2`, so every slot of a buffer has its own name.
"""
function scalar_name(sym)
    u = unwrap(sym)
    if iscall(u) && operation(u) === getindex
        args = arguments(u)
        return Symbol(getname(args[1]), "_", join(string.(args[2:end]), "_"))
    end
    return getname(sym)
end

"""
    param_default_value(defaults, p)

The build-time default of parameter `p` in `defaults` (a system's initial
conditions). A scalarized array element `p[k]` falls back to element `k` of its
whole-array default, which is where [`_scalarize_system`](@ref) keeps it. Returns
`nothing` when the parameter has no default.
"""
function param_default_value(defaults, p)
    haskey(defaults, p) && return unwrap_const(defaults[p])
    u = unwrap(p)
    if iscall(u) && operation(u) === getindex
        base = arguments(u)[1]
        haskey(defaults, base) || return nothing
        idx = map(i -> Int(unwrap_const(i)), arguments(u)[2:end])
        return collect(unwrap_const(defaults[base]))[idx...]
    end
    return nothing
end

"""
    compile_batched(func_expr, target_field) -> RuntimeGeneratedFunction

Wrap a `build_function` expression in a loop over component instances, giving

    (target, u, input, numeric, callables, instances, batch, t)

which writes `target[inst.target_field]` for every instance named in `batch`.
`target_field` is `:states` for a state-derivative map, `:outputs` for an output map
and `:observables` for an observed map.

The body appears once, so what is compiled grows with the number of component
*types*, not with the number of instances; calling the scalar body per instance
instead leaves a real call boundary, since kernel bodies are too large to inline.

The loop's own names are fixed rather than `gensym`ed, so two compiles of the same
component share a type: a `RuntimeGeneratedFunction`'s type is a hash of its argument
names and body. The `#` in each name keeps it from colliding with the body's own
variables.
"""
function compile_batched(func_expr, target_field::Symbol)
    out, state, input, param, callable, iv = func_expr.args[1].args
    body = func_expr.args[2]
    target, u, buffer, numeric, callables, instances, batch, index, inst =
        map(name -> Symbol("#batched#", name),
            (:target, :u, :input, :numeric, :callables, :instances, :batch,
             :index, :instance))
    loop = quote
        Base.@inbounds for $index in $batch
            $inst = $instances[$index]
            $out = Base.view($target, $inst.$target_field)
            $state = Base.view($u, $inst.states)
            $input = Base.view($buffer, $inst.inputs)
            $param = Base.view($numeric, $inst.params)
            $callable = $callables[$inst.position]
            $body
        end
        nothing
    end
    signature = Expr(:tuple, target, u, buffer, numeric, callables, instances,
                     batch, iv)
    return Symbolics._build_and_inject_function(Symbolics,
                                                Expr(:function, signature, loop))
end

"""
    generate_io_function(sys, inputs, outputs; verbose=false, cse=true)

Compile `sys` into the callables one component type needs. `inputs` and `outputs`
name variables of `sys` (array-valued ones are scalarized). Returns a named tuple:

- `f` — state derivatives, `nothing` when the component has no state
- `g` — the declared outputs
- `obs` — every remaining observed variable, `nothing` when there are none
- `mass_matrix`, and the symbolic `states`, `inputs`, `outputs`, `obsstates`,
  `params`, `callable_params` in the order the buffers use
- `reads`, which inputs and states each output and each state derivative depends
  on, as [`dependency_indices`](@ref) lists

All three callables take the same argument list,
`(target, u, input, numeric, callables, instances, batch, t)`, and run over every
instance named in `batch`, as [`compile_batched`](@ref) wraps them.
"""
function generate_io_function(sys, inputs, outputs; verbose=false, cse=true)
    allinputs = convert(Vector{ST},
                        _scalarize_io_syms(getproperty_symbolic.(Ref(sys), inputs)))
    alloutputs = convert(Vector{ST},
                         _scalarize_io_syms(getproperty_symbolic.(Ref(sys), outputs)))

    # Rebuild the system fully scalarized so mtkcompile and codegen see the same
    # scalar symbolics as the (already scalarized) I/O lists above.
    # Read the defaults off the hierarchical system: a nested subsystem's are
    # namespaced here, and flattening drops them from the rebuilt system's store.
    defaults = ModelingToolkitBase.initial_conditions(sys)
    sys = ModelingToolkitBase.expand_connections(_scalarize_system(sys))
    iv = only(independent_variables(sys))

    simplified, eqs, obseqs_sorted, states, params =
        simplify_with_mtkcompile(sys, allinputs, alloutputs; verbose)
    eqs, obseqs_sorted, states =
        pick_best_alias_names(eqs, obseqs_sorted, states, alloutputs, allinputs; verbose)

    diffs = rhs_differentials(vcat(eqs, obseqs_sorted))
    isempty(diffs) || error("Equations contain derivatives on the right-hand side: $diffs")

    # obs can only depend on parameters (including allinputs) or states
    obs_subs = OrderedDict(eq.lhs => eq.rhs for eq in obseqs_sorted)
    obs_deps = setdiff(_all_rhs_symbols(obs_subs), Set(keys(obs_subs)))
    if !(obs_deps ⊆ Set(params) ∪ Set(allinputs) ∪ Set(states) ∪ iv)
        @warn "obs_deps !⊆ params ∪ inputs ∪ unknowns. Difference: $(setdiff(obs_deps, Set(params) ∪ Set(allinputs) ∪ Set(states)))"
    end

    outeqs = Equation[]
    for out in alloutputs
        if out ∈ Set(states)
            push!(outeqs, out ~ out)
        elseif out ∈ keys(obs_subs)
            push!(outeqs, out ~ obs_subs[out])
            # delete from the obs equations but *not* from obs_subs, which other
            # equations may still reference
            deleteat!(obseqs_sorted, findfirst(eq -> isequal(eq.lhs, out), obseqs_sorted))
        else
            throw(ArgumentError("Output $out was neither found in states nor in observed equations."))
        end
    end

    mass_matrix = generate_massmatrix(eqs)
    verbose && @info "Generated mass matrix" mass_matrix

    reads = (output_input = dependency_indices(outeqs, obs_subs, allinputs),
             output_state = dependency_indices(outeqs, obs_subs, states),
             dstate_input = dependency_indices(eqs, obs_subs, allinputs),
             dstate_state = dependency_indices(eqs, obs_subs, states))

    # Nonnumeric (callable) params can't live in the flat Float64 buffer; they are
    # a separate argument after `p`, always present so the signature is uniform.
    callable_params = filter(is_nonnumeric_param, params)
    numeric_params = filter(!is_nonnumeric_param, params)
    args = (states, allinputs, numeric_params, callable_params, iv)

    compile(formulas, target_field) =
        compile_batched(last(build_function(formulas, args...;
                                            cse, expression = Val{true})),
                        target_field)
    f = isempty(eqs) ? nothing : compile(_get_formulas(eqs, obs_subs), :states)
    g = compile(_get_formulas(outeqs, obs_subs), :outputs)
    obsstates = [eq.lhs for eq in obseqs_sorted]
    obs = isempty(obsstates) ? nothing :
        compile(_get_formulas([s ~ s for s in obsstates], obs_subs), :observables)

    return (; f, g, obs, mass_matrix, states, inputs = allinputs,
            outputs = alloutputs, obsstates, reads,
            params = numeric_params,
            callable_params, simplified,
            param_defaults = [param_default_value(defaults, p) for p in numeric_params],
            callable_defaults = [param_default_value(defaults, p)
                                 for p in callable_params],
            equations = eqs, observed = obseqs_sorted, outputeqs = outeqs)
end

end # module KernelCodegen
