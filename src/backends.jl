# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

"""
    ModelBackend

Abstract supertype selecting how a [`SymbolicAWEModel`](@ref) is assembled and
which features it supports. Concrete backends: [`MonolithBackend`](@ref) (the
default) and [`NetworkBackend`](@ref). The backend is a field of the model and is
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
    NetworkBackend()

Optional backend that assembles the model as a NetworkDynamics `Network`,
compiling one kernel per component *type* (flat in node count). The assembly
methods live in the NetworkDynamics package extension, so `using NetworkDynamics`
is required before building a model with this backend. Features without a network
implementation throw [`BackendUnsupportedError`](@ref).
"""
struct NetworkBackend <: ModelBackend end

"""
    BackendUnsupportedError(feature, backend)

Thrown when `feature` (a name string) has no implementation for `backend` — for
example control-function generation on a [`NetworkBackend`](@ref). The message
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
