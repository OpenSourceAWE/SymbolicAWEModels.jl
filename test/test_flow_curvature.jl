# Copyright (c) 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_flow_curvature.jl
# The thin-airfoil pitch-rate moment increment in the continuous VSM modes. VSM
# applies it inside solve!, which a continuous mode does not run every step, so
# the term is re-derived symbolically from the sections' own trailing minus
# leading edge apparent wind. Every mode is built twice, with the solver flag and
# without, so the increment is read as the difference between two models instead
# of inferred from a total that also contains the polar cm.
#
# What a mode owes:
#   - a pitch rate that is blind to translation and reproduces omega·y_airf
#     under rigid rotation
#   - a rate no rigid omega can express under a sign-alternating twist
#   - a moment increment that damps for either sign of the rate
#   - nothing at all when the flag is off
#
# The two modes place the moment differently — ContinuousAero with the strut
# couple, AeroPressure with a pure-moment chordwise shape over its surface nodes
# (couple_shape) — so running both is what covers the placement, not just the rate.

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

@isdefined(test_init!) || include(joinpath(@__DIR__, "util.jl"))
@isdefined(write_pressure_fixture) ||
    include(joinpath(@__DIR__, "pressure_fixture.jl"))

using Test
using SymbolicAWEModels
using SymbolicAWEModels: VortexStepMethod, wing_points, aero_section_columns
using KiteUtils
using LinearAlgebra

"""
    build_curvature_case(case, root, flow_curvature)

Fresh data directory, `Settings` and `SystemStructure` for one mode of the
particle-dynamics 2plate kite, with the solver's `flow_curvature` flag set as
given. Each build gets its own data copy and `system_name`: the flag changes the
emitted equations, so the two must not share a compiled model.
"""
function build_curvature_case(case, root, flow_curvature)
    system_name = case.system_name * (flow_curvature ? "_on" : "_off")
    data_path = joinpath(root, system_name, "2plate_kite")
    mkpath(dirname(data_path))
    cp(joinpath(dirname(@__DIR__), "data", "2plate_kite"), data_path; force=true)
    case.surface_fixture && write_pressure_fixture(data_path)
    set_data_path(data_path)
    set = Settings("system.yaml")
    vsm_set = VortexStepMethod.VSMSettings(
        joinpath(data_path, "vsm_settings.yaml"); data_prefix=false)
    vsm_set.solver_settings.flow_curvature = flow_curvature
    if case.billowing
        for vsm_wing_settings in vsm_set.wings
            vsm_wing_settings.spanwise_panel_distribution =
                VortexStepMethod.BILLOWING
            vsm_wing_settings.billowing_percentage = 8.0
        end
    end
    sys = load_sys_struct_from_yaml(
        joinpath(data_path, "particle_structural_geometry.yaml");
        system_name, set, vsm_set, aero_mode=case.make())
    return set, sys
end

"""
    aero_subsystem(sam, wing)

The wing's aero subsystem in the compiled model, whose `pitch_rate` and `y_airf`
columns are read rather than recomputed, so the assertions cannot drift from the
equations they check.
"""
aero_subsystem(sam, wing) = getproperty(sam.prob.sys, Symbol("aero_$(wing.idx)"))

"""
    has_panel_readout(sam) -> Bool

Whether the backend exposes the per-panel intermediates as observed variables of one
compiled system. The [`KernelBackend`](@ref) assembles many small kernels and has no
such system, so the checks that read `pitch_rate`/`y_airf` are monolith-only. What it
places on the wing — the moment increment itself — is checked on both.
"""
has_panel_readout(sam) = !isnothing(sam.prob.sys)

"""
    panel_pitch_rates(sam, wing)

Every refined panel's rotation rate about its own spanwise axis, as the compiled
model computes it.
"""
panel_pitch_rates(sam, wing) =
    sam.integrator[collect(aero_subsystem(sam, wing).pitch_rate)]

"""
    panel_spanwise_axes(sam, wing)

Every refined panel's body-frame `y_airf`, the axis its pitch rate is about.
"""
function panel_spanwise_axes(sam, wing)
    y_airf = sam.integrator[collect(aero_subsystem(sam, wing).y_airf)]
    return [Vector(y_airf[:, i]) for i in axes(y_airf, 2)]
end

"""
    settle!(sam)

Sync the symbolic aero into the struct without re-solving the circulation, so a
comparison sees the injected state and the frozen induced velocity of the last
refresh.
"""
settle!(sam) = next_step!(sam; dt=1e-6, vsm_interval=0)

"""
    inject_velocity!(sam, nodes, velocities)

Add a world-frame velocity to each wing node and carry it into the integrator
state. Returns the largest injection error relative to the intended velocity,
which must be small or the test is asserting about a state the model never saw.
Measured before [`settle!`](@ref), since that step advances the state itself.
"""
function inject_velocity!(sam, nodes, velocities)
    intended = [Vector(node.vel_w) .+ velocities[k] for (k, node) in enumerate(nodes)]
    for (node, velocity) in zip(nodes, velocities)
        node.vel_w .+= velocity
    end
    init!(sam; reinit_sys=false, lin_vsm=false, prn=false)
    error_norm = maximum(norm(Vector(nodes[k].vel_w) .- intended[k])
                         for k in eachindex(nodes))
    settle!(sam)
    return error_norm / maximum(norm, intended)
end

"""
    rigid_rotation_field(nodes, rate)

World-frame velocity each node picks up from a rigid rotation `rate` about their
centroid.
"""
function rigid_rotation_field(nodes, rate)
    centroid = sum(node.pos_w for node in nodes) / length(nodes)
    return [cross(rate, node.pos_w - centroid) for node in nodes]
end

"""
    twist_rate_field(wing, nodes, speed)

World-frame velocity field that twists each strut about its own chord, the rate
running as a full cosine along the span: each strut's trailing edge moves along the
wing's body `z` and its leading edge against it, scaled by `cos(2π·span)`.

The rate therefore changes sign **twice** across the span. A rigid `omega` gives
`omega · y_airf`, which is monotone in the panel axis and can change sign once, so
no rigid motion reproduces this however the `y_airf` orientation convention runs —
which a simple left/right antisymmetric field does not guarantee, since `y_airf`
already flips at mid span.
"""
function twist_rate_field(wing, nodes, speed)
    column = aero_section_columns(wing, nodes)
    n_struts = Int(wing.vsm_wing.n_unrefined_sections)
    z_world = wing.R_b_to_w * [0.0, 0.0, 1.0]
    velocities = [zeros(3) for _ in nodes]
    for strut in 1:n_struts
        weight = cos(2pi * (strut - 1) / max(n_struts - 1, 1))
        velocities[column[(strut, :TE)]] .+= weight * speed .* z_world
        velocities[column[(strut, :LE)]] .-= weight * speed .* z_world
    end
    return velocities
end

"""
    rigid_rate_residual(rates, axes)

Fraction of `rates` no rigid body rate explains: the least-squares `omega`
minimising `sum((rate - omega·y)^2)` is removed and the remaining norm reported
relative to the original. Near zero means a rigid `omega` describes the motion.
"""
function rigid_rate_residual(rates, axes)
    gram = sum(axis * axis' for axis in axes)
    moments = sum(rates[i] * axes[i] for i in eachindex(rates))
    omega = gram \ moments
    residual = [rates[i] - dot(omega, axes[i]) for i in eachindex(rates)]
    return norm(residual) / norm(rates)
end

@testset "Flow curvature" begin
    root = mktempdir()
    # About the wing nodes' centroid: a few m/s against the ~15 m/s inflow.
    pitch_rate = 3.0
    twist_speed = 1.5
    cases = [
        (name="ContinuousAero", system_name="flow_curvature_cont",
            make=() -> ContinuousAero(), billowing=true, surface_fixture=false),
        (name="AeroPressure", system_name="flow_curvature_press",
            make=() -> AeroPressure(), billowing=false, surface_fixture=true),
    ]

    for case in cases
        @testset "$(case.name)" begin
            # Each model is built while its own case is the active data path, so
            # the two compiled models cache separately.
            set_on, sys_on = build_curvature_case(case, root, true)
            sam_on = SymbolicAWEModel(set_on, sys_on)
            set_off, sys_off = build_curvature_case(case, root, false)
            sam_off = SymbolicAWEModel(set_off, sys_off)
            wing_on, wing_off = sys_on.wings[1], sys_off.wings[1]

            @testset "the flag reaches the wing" begin
                @test SymbolicAWEModels.flow_curvature_enabled(wing_on)
                @test !SymbolicAWEModels.flow_curvature_enabled(wing_off)
                # The term is in the equations, so no shared cached model.
                @test SymbolicAWEModels.get_sys_struct_hash(sys_on) !=
                      SymbolicAWEModels.get_sys_struct_hash(sys_off)
            end

            test_init!(sam_on)
            test_init!(sam_off)
            nodes_on = wing_points(sys_on, wing_on)
            nodes_off = wing_points(sys_off, wing_off)
            settle!(sam_on)
            settle!(sam_off)

            @testset "no rate, no increment" begin
                # Both models sit on the same state, so any moment difference is
                # the increment alone. The kite is not rigid and is still settling
                # at init, so its stations carry a small genuine twist rate; what
                # must hold is that it is negligible beside an injected rate, and
                # that the moment it produces is negligible beside the total.
                moment_scale = norm(wing_on.aero_moment_b)
                @test moment_scale > 1e-3
                rel_moment = norm(wing_on.aero_moment_b - wing_off.aero_moment_b) /
                    moment_scale
                rate_note = ""
                if has_panel_readout(sam_on)
                    rates = panel_pitch_rates(sam_on, wing_on)
                    @test maximum(abs, rates) < 0.01 * pitch_rate
                    rate_note = "max|q|=$(round(maximum(abs, rates); sigdigits=3)) rad/s, "
                end
                println("  [$(case.name)] at rest: ", rate_note,
                    "rel_M=$(round(rel_moment; sigdigits=3))")
                @test rel_moment < 1e-3
            end

            @testset "the flag off leaves the rate bound to zero" begin
                has_panel_readout(sam_off) || @test_skip false
                has_panel_readout(sam_off) &&
                    @test all(iszero, panel_pitch_rates(sam_off, wing_off))
            end

            @testset "blind to translation" begin
                has_panel_readout(sam_on) || @test_skip false
                if has_panel_readout(sam_on)
                # The gather weights sum to zero, so a uniform velocity field
                # contributes nothing to the rate however large it is. Not exact
                # here only because re-initialising moves the structure a little,
                # which is the same drift the at-rest rate already shows; a rigid
                # rotation of the same magnitude moves the rate by O(pitch_rate).
                before = panel_pitch_rates(sam_on, wing_on)
                uniform = [[4.0, -2.0, 1.0] for _ in nodes_on]
                injected = inject_velocity!(sam_on, nodes_on, uniform)
                @test injected < 1e-3
                drift = maximum(abs, panel_pitch_rates(sam_on, wing_on) - before)
                println("  [$(case.name)] translation drift: ",
                    "$(round(drift; sigdigits=3)) rad/s")
                @test drift < 0.01 * pitch_rate
                inject_velocity!(sam_on, nodes_on, [-v for v in uniform])
                end
            end

            @testset "rigid rotation gives omega dot y_airf" begin
                has_panel_readout(sam_on) || @test_skip false
                if has_panel_readout(sam_on)
                before = panel_pitch_rates(sam_on, wing_on)
                axes = panel_spanwise_axes(sam_on, wing_on)
                # About the wing's own spanwise axis, so every panel sees it.
                omega_body = [0.0, pitch_rate, 0.0]
                omega_world = wing_on.R_b_to_w * omega_body
                field = rigid_rotation_field(nodes_on, omega_world)
                injected = inject_velocity!(sam_on, nodes_on, field)
                @test injected < 1e-3
                measured = panel_pitch_rates(sam_on, wing_on) - before
                expected = [dot(omega_body, axis) for axis in axes]
                @test norm(expected) > 0.5 * pitch_rate * sqrt(length(expected))
                rel_error = norm(measured - expected) / norm(expected)
                residual = rigid_rate_residual(measured, axes)
                println("  [$(case.name)] rigid rotation: rel_error=",
                    "$(round(rel_error; sigdigits=3)), rigid residual=",
                    "$(round(residual; sigdigits=3))")
                # The rate is gathered at the struts and interpolated, so it
                # agrees with the panel's own axis only to the mesh resolution.
                @test rel_error < 0.15
                @test residual < 0.15
                inject_velocity!(sam_on, nodes_on, [-v for v in field])
                end
            end

            @testset "opposite twist rates no rigid omega reproduces" begin
                has_panel_readout(sam_on) || @test_skip false
                if has_panel_readout(sam_on)
                before = panel_pitch_rates(sam_on, wing_on)
                axes = panel_spanwise_axes(sam_on, wing_on)
                field = twist_rate_field(wing_on, nodes_on, twist_speed)
                injected = inject_velocity!(sam_on, nodes_on, field)
                @test injected < 1e-3
                measured = panel_pitch_rates(sam_on, wing_on) - before
                @test maximum(abs, measured) > 0.1
                # The rate changes sign twice along the span, which no rigid
                # omega can do.
                @test minimum(measured) < 0.0 < maximum(measured)
                residual = rigid_rate_residual(measured, axes)
                println("  [$(case.name)] antisymmetric twist: max|q|=",
                    "$(round(maximum(abs, measured); sigdigits=3)) rad/s, ",
                    "rigid residual=$(round(residual; sigdigits=3))")
                @test residual > 0.5
                inject_velocity!(sam_on, nodes_on, [-v for v in field])
                end
            end

            @testset "the increment damps for either sign" begin
                # Both models are stepped through the same injected state, so the
                # moment difference is the increment and nothing else. A sign
                # error turns pitch damping into divergence, which no finiteness
                # check would catch.
                powers = Float64[]
                for sign in (1.0, -1.0)
                    omega_body = [0.0, sign * pitch_rate, 0.0]
                    field_on = rigid_rotation_field(nodes_on,
                                                    wing_on.R_b_to_w * omega_body)
                    field_off = rigid_rotation_field(nodes_off,
                                                     wing_off.R_b_to_w * omega_body)
                    inject_velocity!(sam_on, nodes_on, field_on)
                    inject_velocity!(sam_off, nodes_off, field_off)
                    increment = Vector(wing_on.aero_moment_b) -
                        Vector(wing_off.aero_moment_b)
                    push!(powers, dot(increment, omega_body))
                    println("  [$(case.name)] omega_y=$(sign * pitch_rate): ",
                        "dM=$(round(norm(increment); sigdigits=3)) N·m, ",
                        "dM·omega=$(round(powers[end]; sigdigits=3)) W")
                    @test norm(increment) > 1e-6 * norm(wing_off.aero_moment_b)
                    inject_velocity!(sam_on, nodes_on, [-v for v in field_on])
                    inject_velocity!(sam_off, nodes_off, [-v for v in field_off])
                end
                @test all(<(0.0), powers)
                # Linear in the rate, so reversing it reverses the increment.
                @test isapprox(powers[1], powers[2]; rtol=0.25)
            end

            @testset "stepping with the term stays finite" begin
                for _ in 1:5
                    next_step!(sam_on; dt=0.01, vsm_interval=1)
                end
                @test all(isfinite, wing_on.aero_force_b)
                @test all(isfinite, wing_on.aero_moment_b)
                has_panel_readout(sam_on) &&
                    @test all(isfinite, panel_pitch_rates(sam_on, wing_on))
            end
        end
    end

    rm(root; recursive=true, force=true)
end
nothing
