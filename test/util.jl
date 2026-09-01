# SPDX-FileCopyrightText: 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Shared helpers for the test suite. Not picked up automatically
# by runtests.jl (filename must start with `test_`); included
# explicitly from runtests.jl.

using Test
using SymbolicAWEModels
using SymbolicAWEModels: panel_force_slots, panel_force_eqs
using SymbolicAWEModels.ModelingToolkit: Symbolics
using Profile
using LinearAlgebra

"""
    wind_frame_force(force_b, va_b; span_b=[0,1,0]) -> (; drag, lift, side)

Signed projections of a body-frame aero force onto the wind-axis triad, the same
basis the force reconstruction uses: `drag` along the apparent wind, `lift` along
`va × span`, `side` completing the triad. Comparing these separately is what keeps
a drag error visible — `|F|` is lift-dominated, so a norm budget tight enough to
look strict still hides several percent of drag.
"""
function wind_frame_force(force_b, va_b; span_b = [0.0, 1.0, 0.0])
    drag_dir = normalize(collect(va_b))
    lift_dir = normalize(cross(drag_dir, collect(span_b)))
    side_dir = cross(lift_dir, drag_dir)
    return (; drag = dot(force_b, drag_dir),
            lift = dot(force_b, lift_dir),
            side = dot(force_b, side_dir))
end

"""
    pin_wing_nodes!(sys, wing, pinned)

Freeze (or release) every node of `wing` via `fix_static`, so a deformed shape can
be held against the aero without the structure relaxing into it. Read live from
the struct, so it needs no rebuild.
"""
function pin_wing_nodes!(sys, wing, pinned::Bool = true)
    for point in sys.points
        point.is_wing_node && point.wing_idx == wing.idx &&
            (point.fix_static = pinned)
    end
    return nothing
end

"""
    wing_node_positions(sys, wing)

`point idx => copy(pos_w)` for every node of `wing`, the baseline a deformation is
applied from.
"""
wing_node_positions(sys, wing) = Dict(point.idx => copy(point.pos_w)
    for point in sys.points
    if point.is_wing_node && point.wing_idx == wing.idx)

"""
    deform_trailing_edge!(sys, wing, base_pos, delta)

Mirror-antisymmetric deformation from `base_pos`: each station's trailing
edge moves `±delta` along the wing's body z, signed by its span station. Prescribed
twist cannot do this on a `PARTICLE_DYNAMICS` wing — there a multi-point `STATIC`
surface is an inert aero-section marker (`station_eqs.jl`) and the free
points carry the deformation.
"""
function deform_trailing_edge!(sys, wing, base_pos, delta)
    normal_w = (wing.R_b_to_w::Matrix{Float64})[:, 3]
    for point in sys.points
        haskey(base_pos, point.idx) && (point.pos_w .= base_pos[point.idx])
    end
    for surface in sys.stations
        length(surface.point_idxs) >= 2 || continue
        trailing = sys.points[surface.point_idxs[2]]
        side = sign(trailing.pos_cad[2])
        side == 0 && continue
        trailing.pos_w .+= (side * delta) .* normal_w
    end
    return nothing
end

"""
    set_twist!(sys, twists)

Prescribe `twist` per twist-surface name, e.g. `(; left=0.06, right=-0.06)`.
Surfaces not named keep their value. Only meaningful for `STATIC` stations,
where the angle is held rather than integrated.
"""
function set_twist!(sys, twists)
    for surface in sys.stations
        haskey(twists, surface.name) && (surface.twist = twists[surface.name])
    end
    return nothing
end

"""
    asymmetric_twist(delta=0.06)

Mirror-antisymmetric twist: `+delta` on the left surface, `-delta` on the right,
centre flat. Passing `-delta` gives the mirrored configuration.
"""
asymmetric_twist(delta = 0.06) =
    Dict(:left => delta, :center => 0.0, :right => -delta)

"""
    settle_aero!(sam; steps=25, dt=0.02)

Run the model forward, then take one `dt=1e-5` refresh step so the model and a VSM
reference see the same geometry. Deliberately does **not** call `init!`: the whole
point is to measure aero in a moved, deformed state, and re-initialising would
discard exactly the motion the contract is about.
"""
function settle_aero!(sam; steps::Integer = 25, dt = 0.02)
    for _ in 1:steps
        next_step!(sam; dt, vsm_interval = 1)
    end
    next_step!(sam; dt = 1e-5, vsm_interval = 1)
    return nothing
end

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

"""
    scalar_equation_pairs(eq) -> Vector{Tuple}

`eq` as scalar `(left, right)` pairs, splitting an array-valued equation into one
pair per component.
"""
function scalar_equation_pairs(eq)
    left = eq.lhs
    is_array = left isa AbstractArray ||
        Symbolics.symbolic_type(left) === Symbolics.ArraySymbolic()
    is_array || return [(left, eq.rhs)]
    return collect(zip(collect(Symbolics.scalarize(left)),
                       collect(Symbolics.scalarize(eq.rhs))))
end

"""
    evaluate_panel_equations(sections, flow, coefficients, spanwise, scale,
                             orient, chord_weight) -> NamedTuple

One panel of `panel_force_eqs` evaluated numerically. The equations are built on
numeric `sections`/`flow` and substituted forward in the order they are emitted,
so every returned quantity is the expression the model compiles rather than a
reimplementation of it. `coefficients` maps the resulting angle of attack to
`(cl, cd, cm)`, which are folded into `force` and `couple`.
"""
function evaluate_panel_equations(sections, flow, coefficients, spanwise, scale,
                                  orient, chord_weight)
    slots = panel_force_slots(1)
    polar_symbols = [Symbolics.variable(name) for name in (:cl, :cd, :cm)]
    polars = Tuple(_ -> symbol for symbol in polar_symbols)
    values = Dict{Any, Any}()
    for equation in panel_force_eqs(slots, 1, sections, flow, polars, spanwise,
                                    scale, orient, chord_weight, nothing),
        (left, right) in scalar_equation_pairs(equation)
        values[Symbolics.value(left)] =
            Symbolics.substitute(right, values; fold = Val(true))
    end
    readout(slot) = Symbolics.value(values[Symbolics.value(slot)])
    angle_of_attack = Float64(readout(slots.alpha[1]))
    polar_values = Dict(Symbolics.value(symbol) => value for (symbol, value) in
                        zip(polar_symbols, coefficients(angle_of_attack)))
    fold(slot) = Float64(Symbolics.value(Symbolics.substitute(
        readout(slot), polar_values; fold = Val(true))))
    triple(slot) = [Float64(readout(slot[k, 1])) for k in 1:3]
    return (; chord = Float64(readout(slots.chord[1])),
              width = Float64(readout(slots.width[1])),
              q_dyn = Float64(readout(slots.q_dyn[1])),
              alpha = angle_of_attack,
              x_airf = triple(slots.x_airf), y_airf = triple(slots.y_airf),
              z_airf = triple(slots.z_airf), v_eff = triple(slots.v_eff),
              dir_lift = triple(slots.dir_lift),
              dir_drag = triple(slots.dir_drag),
              force = [fold(slots.panel_force[k, 1]) for k in 1:3],
              couple = [fold(slots.panel_couple[k, 1]) for k in 1:3])
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
