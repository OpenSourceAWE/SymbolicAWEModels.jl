# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

"""
    ModelBackend

Abstract supertype selecting how a [`SymbolicAWEModel`](@ref) is assembled and
which features it supports. Concrete backends: [`MonolithBackend`](@ref) (the
default) and [`KernelBackend`](@ref). The backend is a field of the model and is
the dispatch axis for all backend-varying behaviour (problem assembly,
linearization, control-function generation).
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
enormous function. It is not free, but the first solve repays it: on SK100 a dense
finite-difference Jacobian over 1305 states is 1306 evaluations against forward
mode's ~126, which measured 87.8 s against 212.4 s for `init!` and 1.42× on the wall
clock of a step.
"""
default_autodiff(::ModelBackend) = AutoFiniteDiff()
default_autodiff(::KernelBackend) = AutoForwardDiff()

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
