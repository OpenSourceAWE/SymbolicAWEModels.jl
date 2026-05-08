# SPDX-FileCopyrightText: 2026 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

# Shared helpers for the test suite. Not picked up automatically
# by runtests.jl (filename must start with `test_`); included
# explicitly from runtests.jl.

using Test
using SymbolicAWEModels

"""
    validate_rhs_allocs(sam; max_bytes=0)

Run the ODE RHS three times (two warmups, one measured) and
`@test` that allocations are within `max_bytes`. A nonzero
count usually means a `Vector{Num}` intermediate (e.g.
`vec/scalar` or an unscalarised cross product) survived MTK
codegen.
"""
function validate_rhs_allocs(sam; max_bytes::Integer = 0)
    isnothing(sam.integrator) && error(
        "validate_rhs_allocs: integrator not initialised; " *
        "call init!(sam) first.")
    f = sam.integrator.f
    u = copy(sam.integrator.u)
    p = sam.integrator.p
    t = sam.integrator.t
    du = similar(u)
    f(du, u, p, t)
    f(du, u, p, t)
    bytes = @allocated f(du, u, p, t)
    @test bytes <= max_bytes
    return bytes
end

"""
    test_init!(sam; max_bytes=0, kwargs...)

Wrapper around `init!` for the test suite. Forwards `kwargs`
to `init!`, then runs `validate_rhs_allocs(sam; max_bytes)`
to ensure the generated ODE RHS is allocation-clean. Returns
the integrator (same as `init!`).
"""
function test_init!(sam; max_bytes::Integer = 0, kwargs...)
    integ = init!(sam; kwargs...)
    validate_rhs_allocs(sam; max_bytes)
    return integ
end
