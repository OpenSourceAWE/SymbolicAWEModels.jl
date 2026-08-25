# Copyright (c) 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_live_polar.jl
# AeroPressure with live polars: every panel's polar regenerated from its deformed
# shape each solve instead of read off a flap angle.
# - construction fits a base airfoil per panel and records the control-point frame
# - the flap axis is gone: no panel carries a twist_surface, so the RHS has no δ
# - a refresh leaves every panel on a TAYLOR polar fitted about its own α
# - a chordwise deformation moves the polar, and the frame maths that produces it
#   agrees with a hand-computed chord frame
#
# What AeroPressure's scatter does with the resulting forces is in
# test_pressure_aero.jl and is unchanged by the polar source.

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

@isdefined(test_init!) || include(joinpath(@__DIR__, "util.jl"))

using Test
using SymbolicAWEModels
using SymbolicAWEModels: VortexStepMethod, refresh_particle_aero!, SimFloat,
                         chord_frame_coordinates, update_live_deflection!,
                         LivePolarState
using VortexStepMethod: AirfoilAero, TAYLOR, calculate_cl
using KiteUtils
using LinearAlgebra

@isdefined(write_pressure_fixture) ||
    include(joinpath(@__DIR__, "pressure_fixture.jl"))

@testset "AeroPressure live polars" begin
    pkg_root = dirname(@__DIR__)
    tmpdir = mktempdir()
    data_path = joinpath(tmpdir, "2plate_kite")
    cp(joinpath(pkg_root, "data", "2plate_kite"), data_path; force=true)
    write_pressure_fixture(data_path)

    set, sys = load_pressure_sys(data_path, AeroPressure(; live_polars=true))
    wing = sys.wings[1]
    mode = wing.aero
    panels = wing.vsm_aero.panels

    @testset "build" begin
        @test mode.live_polars
        @test mode.live isa LivePolarState
        state = mode.live
        @test length(state.control_point) == length(panels)
        @test all(length(pts) >= 2 for pts in state.control_point)
        @test length(state.source.base) == length(panels)
        @test all(length(d) == length(state.source.basis.x) for d in state.deflection)
        # The δ axis is gone, so no panel is wired to a flap surface.
        @test all(iszero, mode.panel_twist_surface)
    end

    wing.va_b .= SimFloat[15.0, 0.0, 3.5]
    wing.ω_b .= SimFloat[0.0, 0.0, 0.0]
    wing.vsm_solver.density = 1.225
    va_vals = zeros(SimFloat, 3, length(sys.points))
    for p in sys.points
        @views va_vals[:, p.idx] .= wing.va_b
    end
    refresh_particle_aero!(mode, wing, sys.points, va_vals)

    @testset "refresh fits every panel" begin
        @test all(p -> p.aero_model == TAYLOR, panels)
        @test all(p -> all(isfinite, p.cl_coeffs), panels)
        @test all(isfinite, mode.traction)
        @test norm(mode.traction) > 0.0
        # Every panel converged inside the window its own fit was built on.
        @test AirfoilAero.polar_drift(mode.live.source,
            collect(wing.vsm_solver.lr.alpha_dist)) <= 1.0
        # The fit is a local expansion: at its own centre it is the network's answer.
        for (i, panel) in enumerate(panels)
            @test calculate_cl(panel, mode.live.source.alpha_ref[i]) ≈ panel.cl_coeffs[1]
        end
    end

    @testset "chord frame coordinates" begin
        panel = panels[1]
        le_mid = 0.5 .* (Vector(panel.LE_point_1) .+ Vector(panel.LE_point_2))
        probe = le_mid .+ 0.4 * panel.chord .* Vector(panel.x_airf) .+
                0.03 * panel.chord .* Vector(panel.z_airf)
        fraction, offset = chord_frame_coordinates(panel, probe)
        @test fraction ≈ 0.4
        @test offset ≈ 0.03
        @test all(iszero, chord_frame_coordinates(panel, le_mid))
    end

    @testset "a chordwise deformation moves the polar" begin
        state = mode.live
        alpha = collect(SimFloat, wing.vsm_solver.lr.alpha_dist)
        # Fit about the angles the comparison reads, so the two answers differ only
        # in the shape and not in where the expansion was centred.
        SymbolicAWEModels.refit_live_polars!(mode, wing, alpha)
        flat = [calculate_cl(panels[i], alpha[i]) for i in eachindex(panels)]

        # Bend the mean line up by 2% of chord at mid chord, the deformation a chord
        # beam would report; the fixture's wing carries only its chord ends, so the
        # deflection is injected where update_live_deflection! would write it.
        for i in eachindex(panels)
            state.deflection[i] .= @. 0.08 * state.source.basis.x *
                                      (1 - state.source.basis.x)
        end
        SymbolicAWEModels.refit_live_polars!(mode, wing, alpha)
        cambered = [calculate_cl(panels[i], alpha[i]) for i in eachindex(panels)]
        @test all(cambered .> flat)

        # And the structure's own (undeformed) shape puts it back.
        update_live_deflection!(mode, wing, sys.points)
        @test all(d -> maximum(abs, d) < 1e-9, state.deflection)
        SymbolicAWEModels.refit_live_polars!(mode, wing, alpha)
        @test [calculate_cl(panels[i], alpha[i]) for i in eachindex(panels)] ≈ flat
    end
end
nothing
