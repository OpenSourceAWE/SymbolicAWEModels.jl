# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

"""
    ModelBackend

Abstract supertype selecting how a [`SymbolicAWEModel`](@ref) is assembled and
which features it supports. Concrete backends: [`MonolithBackend`](@ref) (the
default) and [`KernelBackend`](@ref). The backend is a field of the model and is
the dispatch axis for all backend-varying behaviour (problem assembly,
linearization, control-function generation).

Backends differ at assembly only. Every equation is written once, in
`components.jl`, and both backends build from that one definition — the monolith
through the generators in `generate_system/`, the kernel by wrapping it in a
component. A quantity one backend reports and the other does not is a missing
assembly or readout, never a licence to restate the math on the other side.
"""
abstract type ModelBackend end

"""
    MonolithBackend()

Default backend. `create_sys!` builds one flattened ModelingToolkit `System` that
`mtkcompile` turns into a single RHS. Supports every feature (linearization,
control functions). Compile time grows with the total node count.
"""
struct MonolithBackend <: ModelBackend end

"""
    KernelBackend()

Backend that compiles one kernel per component *type* and runs them from a
build-time schedule over gather/scatter buffers. Compile time is flat in node
count, because refining a model adds instances of kernels that already exist.
Features without an implementation throw [`BackendUnsupportedError`](@ref).
"""
struct KernelBackend <: ModelBackend end

"""
    backend_tag(backend) -> String

The backend's mark in a serialized model's filename. Two backends assemble the
same `SystemStructure` into different artefacts, so they need separate cache
entries; the [`MonolithBackend`](@ref) tag is empty so its existing bins keep
loading.
"""
backend_tag(::ModelBackend) = ""
backend_tag(::KernelBackend) = "_kernel"

"""
    default_autodiff(backend)

How [`init!`](@ref) differentiates the right-hand side for the solver's Jacobian.

The [`MonolithBackend`](@ref) takes `AutoFiniteDiff()`: forward mode would compile its
single enormous right-hand side a second time, at `ForwardDiff.Dual`, which costs more
at the first `init!` than it saves.

The [`KernelBackend`](@ref) takes `AutoForwardDiff()`. Its right-hand side is small
per-kernel functions and `buffers` already keeps a scratch set per element type, so
the `Dual` specialization is one more compilation of each kernel rather than of one
enormous function. It is not free, but the first solve repays it: on a large kite a dense
finite-difference Jacobian over 1305 states is 1306 evaluations against forward
mode's ~126, which measured 87.8 s against 212.4 s for `init!` and 1.42× on the wall
clock of a step.
"""
default_autodiff(::ModelBackend) = AutoFiniteDiff()
default_autodiff(::KernelBackend) = AutoForwardDiff()

"""
    default_analytic_jacobian(backend) -> Bool

Whether [`init!`](@ref) gives the solver an analytical Jacobian rather than letting
it differentiate the right-hand side numerically.

The [`KernelBackend`](@ref) does, and it is the whole reason its structure is worth
keeping: the model is a layered composition of small components, so
[`build_jacobian`](@ref) differentiates each component at its own width and composes
the result through the constant wiring. On a large kite that is one pass per kernel against
109 chunk-12 passes over 1305 states.

The [`MonolithBackend`](@ref) does not. Its route is MTK's `jac=true`, which
differentiates the flattened right-hand side *symbolically*: on the smallest model in
the repo that costs 25.5 s against 3.3 s to build, and it does not run, because a
registered numerical leaf — the wind profile, an aerodynamic polar — has no symbolic
derivative and leaves a `Differential` in the generated matrix.
"""
default_analytic_jacobian(::ModelBackend) = false
default_analytic_jacobian(::KernelBackend) = true

"""
    default_sparse(backend) -> Bool

Whether [`init!`](@ref) hands the solver the Jacobian's sparsity pattern rather
than letting it factorize a dense matrix.

The [`KernelBackend`](@ref) does. It knows the pattern already — [`state_sparsity`](@ref)
walks it out of the wiring — and a structure is mostly zeros: on a 392-point beam the
Jacobian is 1301 states square and 5.15% dense, where factorizing it dense is most of a
step. That measured 35.5 ms a step against 50.2 ms. It also threads where the dense one
does not: eight models stepping at once manage 48.4 steps/s against 4.8, because a dense
factorization opens a BLAS thread pool per calling worker over the same cores.

The [`MonolithBackend`](@ref) does not, so the bins it has already written keep loading.
"""
default_sparse(::ModelBackend) = false
default_sparse(::KernelBackend) = true

"""
    default_linsolve(backend)

Which factorization [`init!`](@ref)'s solver uses, or `nothing` to leave the choice to
LinearSolve.

The [`KernelBackend`](@ref) takes `KLUFactorization`. [`default_sparse`](@ref) gives it a
sparse Jacobian, for which LinearSolve would otherwise pick UMFPACK; KLU refactorizes one
pattern over and over, which is what a BDF integrator does with it, and measured 1.3x
faster at every worker count on a large kite.

The [`MonolithBackend`](@ref) takes `nothing`: its Jacobian is dense, so neither sparse
factorization applies.
"""
default_linsolve(::ModelBackend) = nothing
default_linsolve(::KernelBackend) = KLUFactorization()

const DEFAULT_BACKEND = Ref{ModelBackend}(MonolithBackend())

"""
    default_backend()

The [`ModelBackend`](@ref) a [`SymbolicAWEModel`](@ref) is built with when its
constructor is not given one; [`MonolithBackend`](@ref) unless changed with
[`default_backend!`](@ref).
"""
default_backend() = DEFAULT_BACKEND[]

"""
    default_backend!(backend)

Set the [`ModelBackend`](@ref) new [`SymbolicAWEModel`](@ref)s default to, and
return it. Runs a whole model — or a whole test suite — on another backend without
touching the call sites.
"""
default_backend!(backend::ModelBackend) = (DEFAULT_BACKEND[] = backend; backend)

"""
    BackendUnsupportedError(feature, backend)

Thrown when `feature` (a name string) has no implementation for `backend` — for
example control-function generation on a [`KernelBackend`](@ref). The message
points at the backend that does support the feature.
"""
struct BackendUnsupportedError <: Exception
    feature::String
    backend::ModelBackend
end

function Base.showerror(io::IO, e::BackendUnsupportedError)
    print(io, "`", e.feature, "` is not available on the ",
          nameof(typeof(e.backend)),
          ". Rebuild the model with `MonolithBackend()` to use it.")
end
