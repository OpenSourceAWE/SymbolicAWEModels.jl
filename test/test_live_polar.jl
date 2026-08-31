# Copyright (c) 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_live_polar.jl
# AeroPressure with live polars: every panel's polar regenerated from its deformed
# shape each solve instead of read off a flap angle.
# - construction fits a base airfoil per panel and records the control-point frame
# - the flap axis is gone: no panel carries a twist_surface, so the RHS has no δ
# - a refresh leaves every panel on a rewritten polar centred on its own α
# - a chordwise deformation moves the polar, and the frame maths that produces it
#   agrees with a hand-computed chord frame
# - against the tabulated mode: the scatter is bit-identical, only the polar source
#   differs, and an undeformed panel's polar is the network's own answer
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
using VortexStepMethod: AirfoilAero, POLAR_VECTORS, calculate_cl, calculate_cm
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
        stations = SymbolicAWEModels.control_point_stations(
            sys.twist_surfaces, sys.points, wing)
        # Stations are read off the twist surfaces, not inferred from geometry.
        @test length(state.control_point) == length(stations)
        @test all(sort(a) == sort(b) for (a, b) in zip(state.control_point, stations))
        @test all(length(pts) >= 2 for pts in state.control_point)
        @test length(state.station_deflection) == length(state.control_point)
        @test length(state.panel_blend) == length(panels)
        # No panel may draw its load from more than one station: a panel straddling
        # two spanwise planes gets a chordwise profile with a step no basis can hold.
        for i in eachindex(panels)
            used = unique(mode.station_point[i])
            @test count(group -> any(in(group), used), stations) == 1
        end
        # A station may not name the same node twice: a twist surface lists its
        # chord ends alongside control points that often sit on them, and a fit
        # handed one chord fraction twice has no answer.
        for group in stations
            @test length(unique(group)) == length(group)
            @test length(unique(round.(sys.points[i].pos_cad; digits=9)
                                for i in group)) == length(group)
        end
        # Every panel blends between two neighbouring stations, in range.
        @test all(1 <= lo <= hi <= length(stations) && hi - lo <= 1 && 0 <= w <= 1
                  for (lo, hi, w) in state.panel_blend)
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

    @testset "refresh samples every panel" begin
        # A rewritten table keeps the shape the panel was built with; this wing's
        # polars carry no delta axis, and a live window is what marks it as live.
        @test all(p -> p.aero_model == POLAR_VECTORS, panels)
        @test all(p -> p.alpha_window > 0, panels)
        @test all(p -> all(isfinite, p.cl_coeffs), panels)
        @test all(isfinite, mode.traction)
        @test norm(mode.traction) > 0.0
        # Every panel converged inside the range its own polar was sampled over.
        @test AirfoilAero.polar_drift(mode.live.source,
            collect(wing.vsm_solver.lr.alpha_dist)) <= 1.0
        # A sample is the polar: at the reference angle it is the middle sample.
        offsets = mode.live.source.settings.offsets
        middle = (length(offsets) + 1) ÷ 2
        for (i, panel) in enumerate(panels)
            @test calculate_cl(panel, mode.live.source.alpha_ref[i]) ≈
                  panel.cl_coeffs[middle]
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
        # Sample about the angles the comparison reads, so the two answers differ
        # only in the shape and not in where the grid was centred.
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

    # The values the pipeline installs are the network's own answers at the panel's
    # own sample angles, so nothing between deforming the shape and writing the polar
    # has altered them. There is no fit residual left to bound separately: between
    # two samples the polar is the straight line between two exact points.
    @testset "an undeformed polar is the network's own answer" begin
        state = mode.live
        settings = state.source.settings
        solver = wing.vsm_solver
        alpha = collect(SimFloat, solver.lr.alpha_dist)
        offsets = settings.offsets
        for (i, panel) in enumerate(panels)
            reynolds = solver.density * norm(panel.va) * panel.chord / solver.mu
            samples = AirfoilAero.neuralfoil_aero(state.source.base[i],
                rad2deg.(alpha[i] .+ offsets), reynolds;
                model_size=settings.model_size, n_crit=settings.n_crit)
            @test panel.cl_coeffs ≈ samples.CL rtol = 1e-8
            @test panel.cd_coeffs ≈ samples.CD rtol = 1e-8
            @test panel.cm_coeffs ≈ samples.CM rtol = 1e-8
            @test panel.alpha_knots ≈ alpha[i] .+ offsets
            @test panel.alpha_ref ≈ alpha[i]
            @test panel.alpha_window ≈ maximum(abs, offsets)
        end
    end

    # The forces follow the deformed shape; so must the pattern that places them.
    # Before this, a deformation changed how hard a panel pulled but not where.
    @testset "the traction pattern follows the deformed shape" begin
        state = mode.live
        @test all(p -> all(isfinite, p), state.contour_cp)
        @test all(p -> all(>(0.0), p), state.contour_cf)
        # The contour is the one the surface→point map was built on, held at δ = 0.
        for (i, panel) in enumerate(panels)
            @test length(state.contour_cp[i]) == length(mode.station_point[i])
            @test state.contour_x[i] ≈
                  VortexStepMethod.section_surface(panel.section_aero, 0.0, 0.0)[1]
        end

        # Undeformed, the contour is the reference one and nothing has moved.
        @test all(i -> state.contour_y[i] ≈ state.contour_y_ref[i],
                  eachindex(panels))

        flat = copy(mode.traction)
        alpha = collect(SimFloat, wing.vsm_solver.lr.alpha_dist)
        for i in eachindex(panels)
            state.deflection[i] .= @. 0.08 * state.source.basis.x *
                                      (1 - state.source.basis.x)
        end
        SymbolicAWEModels.refit_live_polars!(mode, wing, alpha)
        SymbolicAWEModels.refresh_live_pressure!(mode, wing)
        SymbolicAWEModels.freeze_traction_pattern!(mode, wing)
        @test all(isfinite, mode.traction)
        @test !(mode.traction ≈ flat)
        # The contour moved with the deformation, so the normals and the segment
        # areas the traction is built on are the deformed section's.
        @test !(state.contour_y[1] ≈ state.contour_y_ref[1])
        @test state.contour_x[1] ≈
              VortexStepMethod.section_surface(panels[1].section_aero, 0.0, 0.0)[1]
        # And it is a different distribution, not the same one scaled up.
        scale = norm(mode.traction) / norm(flat)
        @test !(mode.traction ≈ scale .* flat)

        update_live_deflection!(mode, wing, sys.points)
        SymbolicAWEModels.refit_live_polars!(mode, wing, alpha)
        SymbolicAWEModels.refresh_live_pressure!(mode, wing)
        SymbolicAWEModels.freeze_traction_pattern!(mode, wing)
        @test mode.traction ≈ flat
        @test state.contour_y[1] ≈ state.contour_y_ref[1]
    end

    # Live polars change where cl/cd/cm come from and nothing else. The surface→point
    # map and the shares that scatter a panel's force over it are geometry, so the two
    # modes must agree on them exactly; if they ever drift, the live mode has started
    # to be a second scatter rather than a second polar source.
    @testset "against the tabulated mode" begin
        set_t, sys_t = load_pressure_sys(data_path, AeroPressure())
        wing_t = sys_t.wings[1]
        mode_t = wing_t.aero
        @test mode.station_point == mode_t.station_point
        @test mode.node_residual_share == mode_t.node_residual_share
        @test mode.node_couple_shape == mode_t.node_couple_shape
        @test mode.section_left_strut == mode_t.section_left_strut
        @test mode.section_left_weight ≈ mode_t.section_left_weight
        # The two never share a serialized model: one wires δ into the RHS, one does not.
        @test SymbolicAWEModels.aero_hash_id(mode) !=
              SymbolicAWEModels.aero_hash_id(mode_t)

        wing_t.va_b .= wing.va_b
        wing_t.ω_b .= wing.ω_b
        wing_t.vsm_solver.density = wing.vsm_solver.density
        va_t = zeros(SimFloat, 3, length(sys_t.points))
        for p in sys_t.points
            @views va_t[:, p.idx] .= wing_t.va_b
        end
        refresh_particle_aero!(mode_t, wing_t, sys_t.points, va_t)

        force = sum(eachcol(wing.vsm_solver.sol.f_body_3D))
        force_t = sum(eachcol(wing_t.vsm_solver.sol.f_body_3D))
        # Loose on magnitude on purpose: this fixture's table is not NeuralFoil's,
        # so the two are different data for the same section and only the direction
        # and the order of magnitude are a shared claim.
        @test 0.5 < norm(force) / norm(force_t) < 2.0
        @test dot(force, force_t) / (norm(force) * norm(force_t)) > 0.99
    end
end
nothing
