# SPDX-FileCopyrightText: 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Shared helpers for the test suite. Not picked up automatically
# by runtests.jl (filename must start with `test_`); included
# explicitly from runtests.jl.

using Test
using SymbolicAWEModels
using Profile

"""
    validate_rhs_allocs(sam; max_bytes=0, diagnose=true)

Run the ODE RHS three times (two warmups, one measured) and
`@test` that allocations are within `max_bytes`. A nonzero
count usually means a `Vector{Num}` intermediate (e.g.
`vec/scalar` or an unscalarised cross product) survived MTK
codegen as a runtime broadcast. When `diagnose=true` and
allocations exceed `max_bytes`, samples `Profile.Allocs` and
prints the top allocation sites with stack traces.
"""
function validate_rhs_allocs(
    sam; max_bytes::Integer = 0, diagnose::Bool = true)
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
    if bytes > max_bytes && diagnose
        diagnose_rhs(f, du, u, p, t)
    end
    @test bytes <= max_bytes
    return bytes
end

"""
    validate_sysstate_roundtrip(sam; rtol=1e-10)

`@test` that this model survives a round trip through a
`SysLog` on disk: log the state, scramble the structure, load
it back, re-seed the integrator, and check `integrator.u` is
unchanged. That is the proof that everything reaching `u0` is
covered, without needing to know which fields those are — and
which they are depends on the model and on what the compiler
tears, so running this for every model the suite builds is what
covers the combinations.

`reinit_sys=false` keeps `init!` from calling
`reinit!(sys_struct, set)`, which would reset positions from
the CAD frame and discard the loaded state.

Scrambling matters: without it an `update_from_sysstate!` that
restored nothing would still pass, since the structure already
holds the right values. `sys.state_vars` scrambles every field a
component can carry state in, with no `DynamicsType` filter — a
field that never reaches `u0` is invisible to the assertion
anyway, so filtering can only lose coverage.

`rtol` is not zero because bodies are logged in the body frame
and rebuilt into the principal frame, which costs a few ULP in
the quaternion/matrix conversions.

The temporary directory is left to Julia's exit-time cleanup:
`load_log` memory-maps the Arrow file, and on Windows a mapped
file cannot be deleted while the mapping is alive.
"""
function validate_sysstate_roundtrip(sam; rtol = 1e-10)
    sys = sam.sys_struct
    expected = copy(sam.integrator.u)
    isempty(expected) && return nothing
    path = mktempdir()
    logger = Logger(sam, 1; precision = Float64)
    log!(logger, SysState(sam; precision = Float64))
    save_log(logger, "roundtrip", false; path)
    reloaded = load_log("roundtrip"; path)
    sys.state_vars = 1.5 .* vec(sys.state_vars) .+ 0.25
    update_from_sysstate!(sys, reloaded.syslog[1])
    init!(sam; remake = false, reinit_sys = false, prn = false)
    @test length(sam.integrator.u) == length(expected)
    @test isapprox(sam.integrator.u, expected; rtol)
    return nothing
end

"""
    test_init!(sam; max_bytes=0, diagnose=true, roundtrip=true, kwargs...)

Wrapper around `init!` for the test suite. Forwards `kwargs`
to `init!`, then runs `validate_rhs_allocs(sam; max_bytes,
diagnose)` to ensure the generated ODE RHS is allocation-
clean, and `validate_sysstate_roundtrip(sam)` unless
`roundtrip=false`. Returns the integrator (same as `init!`).
"""
function test_init!(
    sam;
    max_bytes::Integer = 0,
    diagnose::Bool = true,
    roundtrip::Bool = true,
    kwargs...)
    integ = init!(sam; kwargs...)
    validate_rhs_allocs(sam; max_bytes, diagnose)
    roundtrip && validate_sysstate_roundtrip(sam)
    return integ
end

function diagnose_rhs(f, du, u, p, t)
    f(du, u, p, t)
    f(du, u, p, t)
    Profile.Allocs.clear()
    GC.enable(false)
    try
        Profile.Allocs.@profile sample_rate=1.0 f(
            du, u, p, t)
    finally
        GC.enable(true)
    end
    results = Profile.Allocs.fetch()
    n = length(results.allocs)
    if n == 0
        println(stderr,
            "[validate_rhs_allocs] Profile.Allocs sampled " *
            "no allocations on the next call — they may be " *
            "GC-triggered or already amortised. Try " *
            "increasing run count or inspect generated eqs.")
        return nothing
    end
    println(stderr,
        "\n[validate_rhs_allocs] RHS allocated. ",
        "Profile.Allocs captured ", n,
        " sample(s); top sites:")
    sorted = sort(results.allocs; by = a -> -a.size)
    seen = Set{Tuple{Any, Symbol, Int}}()
    shown = 0
    for a in sorted
        st = a.stacktrace
        top_func = isempty(st) ? :_ : st[1].func
        top_line = isempty(st) ? 0 : st[1].line
        sig = (a.type, top_func, top_line)
        sig in seen && continue
        push!(seen, sig)
        shown += 1
        shown > 8 && break
        println(stderr,
            "  [", shown, "] type=", a.type,
            "  size=", a.size, " B")
        for fr in first(st, 6)
            println(stderr, "      ",
                fr.file, ":", fr.line, "  ", fr.func)
        end
    end
    return nothing
end
