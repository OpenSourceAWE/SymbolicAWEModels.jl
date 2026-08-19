# Copyright (c) 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_pressure_aero.jl
# What is AeroPressure's alone: distributing VSM's per-section load onto structural
# points using the airfoil surface traction (−Cp·n̂ + cf·ŝ) as the pattern, anchored
# to VSM's f_body_3D. Uses the particle 2plate kite with a synthetic section_aero
# fixture (the human-readable .dat/Cp/cf "files-provided" route, no NeuralFoil):
# - construction builds the static surface→point map and passes the frame guard
# - a too-tight frame tolerance errors (misalignment guard)
# - the distributed point forces sum exactly to VSM's total (anchoring)
#
# The contract AeroPressure shares with ContinuousAero is in
# test_continuous_modes.jl, and the one it shares with every aero mode (including
# that it compiles, steps and stays bounded) in test_aero_modes.jl.

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

@isdefined(test_init!) || include(joinpath(@__DIR__, "util.jl"))

using Test
using SymbolicAWEModels
using SymbolicAWEModels: VortexStepMethod, refresh_particle_aero!, SimFloat,
                         loft_contour_node, aero_inflow_groups, wing_points,
                         aero_scatter_entries, surface_node_forces
using KiteUtils
using LinearAlgebra

@isdefined(write_pressure_fixture) ||
    include(joinpath(@__DIR__, "pressure_fixture.jl"))

"""
    load_pressure_sys(data_path, aero_mode)

`Settings` and `SystemStructure` for the particle 2plate kite carrying
`aero_mode`, reading the fixture-patched geometry under `data_path`.
"""
function load_pressure_sys(data_path, aero_mode)
    set_data_path(data_path)
    set = Settings("system.yaml")
    vsm_set = VortexStepMethod.VSMSettings(
        joinpath(data_path, "vsm_settings.yaml"); data_prefix=false)
    particle_yaml = joinpath(data_path, "particle_structural_geometry.yaml")
    sys = load_sys_struct_from_yaml(particle_yaml; system_name="pressure_test",
        set, vsm_set, aero_mode)
    return set, sys
end

@testset "AeroPressure" begin
    pkg_root = dirname(@__DIR__)
    tmpdir = mktempdir()
    data_path = joinpath(tmpdir, "2plate_kite")
    cp(joinpath(pkg_root, "data", "2plate_kite"), data_path; force=true)
    n_node = write_pressure_fixture(data_path)

    set, sys = load_pressure_sys(data_path, AeroPressure())
    wing = sys.wings[1]
    mode = wing.aero
    wing_pts = [p for p in sys.points if p.is_wing_node && p.wing_idx == wing.idx]

    @testset "surface→point map" begin
        @test mode isa AeroPressure
        @test !isempty(mode.station_point)
        @test length(mode.station_point) == length(wing.vsm_aero.panels)
        @test all(length(sp) == n_node for sp in mode.station_point)
        @test all(idx -> idx in [p.idx for p in wing_pts],
                  reduce(vcat, mode.station_point))
        @test all(p.section_aero !== nothing for p in wing.vsm_aero.panels)
    end

    @testset "per-section inflow groups" begin
        groups, section_group = aero_inflow_groups(mode, wing, wing_points(sys, wing))
        @test length(groups) == length(mode.section_left_strut)
        @test length(unique(section_group)) == length(groups)
        @test length(unique(groups)) > 1
        @test all(g -> isapprox(sum(last, g), 1.0; atol=1e-10), groups)
    end

    @testset "frame guard rejects misalignment" begin
        @test_throws ErrorException load_pressure_sys(
            data_path, AeroPressure(frame_tol_frac=0.01))
    end

    # Exercise the refresh directly (no ODE compile): set the operating point the
    # way update_sys_struct!/sync_aero_density! would, then distribute. Alpha is
    # 13°, clear of the α ≈ 8° sign change of the body-origin moment.
    wing.va_b .= SimFloat[15.0, 0.0, 3.5]
    wing.ω_b .= SimFloat[0.0, 0.0, 0.0]
    wing.vsm_solver.density = 1.225
    va_vals = zeros(SimFloat, 3, length(sys.points))
    for p in sys.points
        @views va_vals[:, p.idx] .= wing.va_b
    end
    refresh_particle_aero!(mode, wing, sys.points, va_vals)

    @testset "frozen traction pattern" begin
        n_panels = length(wing.vsm_aero.panels)
        @test size(mode.traction, 2) == sum(length, mode.station_point)
        @test size(mode.traction_net) == (3, n_panels)
        @test all(isfinite, mode.traction)
        @test all(isfinite, mode.traction_net)
        @test norm(mode.traction) > 0.0
        @test all(offset -> all(isfinite, offset), values(mode.point_offset))
    end

    # The distributed point forces sum to VSM's total, per panel and overall.
    # Built from `aero_scatter_entries` and `mode.point_offset` — the same map and
    # the same constants `aero_component` wires up — so this cannot drift from the
    # scatter it checks the way a re-derivation of the residual algebra would.
    @testset "correct total force (per VSM panel and total)" begin
        sol = wing.vsm_solver.sol
        nodes = wing_points(sys, wing)
        entries = aero_scatter_entries(mode, wing, nodes)
        point_force = Dict(p.idx => zeros(SimFloat, 3) for p in nodes)
        panel_sum = [zeros(SimFloat, 3) for _ in eachindex(wing.vsm_aero.panels)]
        for (panel, column, force_weight, _) in entries
            share = force_weight .* Vector(sol.f_body_3D[:, panel])
            point_force[nodes[column].idx] .+= share
            panel_sum[panel] .+= share
        end
        for (point_idx, offset) in mode.point_offset
            point_force[point_idx] .+= offset
        end
        # Per panel the scattered share is the whole panel force: the frozen
        # offsets cancel across a panel by construction.
        for (i, total) in enumerate(panel_sum)
            f_vsm = Vector(sol.f_body_3D[:, i])
            @test norm(total - f_vsm) <= 1e-9 * (norm(f_vsm) + 1)
        end
        total_pts = sum(values(point_force))
        total_vsm = vec(sum(sol.f_body_3D, dims=2))
        @test norm(total_vsm) > 1.0
        @test norm(total_pts - total_vsm) / norm(total_vsm) < 1e-10
    end

    # Force sums survive a placement error untouched, moments do not. Comparing
    # each panel's loads about the body origin where they are applied against the
    # same loads at the surface node they came from separates nearest-node lumping
    # from a wrong distribution. Both residuals are load-centre offsets in chords,
    # ‖ΔM‖/(‖F‖·c): ‖M‖ about the body origin passes through zero inside the flight
    # range, so normalising by it would divide by a vanishing quantity. Both sides
    # use src (`aero_scatter_entries`, `surface_node_forces`), not a re-derivation.
    @testset "moment placement (per VSM panel and total)" begin
        sol = wing.vsm_solver.sol
        panels = wing.vsm_aero.panels
        nodes = wing_points(sys, wing)
        rot_cad_to_body = wing.R_b_to_c'
        point_pos_b = Dict(p.idx => rot_cad_to_body * (p.pos_cad - wing.pos_cad)
                           for p in nodes)

        applied_total = zeros(SimFloat, 3)
        for (panel, column, force_weight, _) in aero_scatter_entries(mode, wing, nodes)
            applied_total .+= cross(point_pos_b[nodes[column].idx],
                force_weight .* Vector(sol.f_body_3D[:, panel]))
        end
        for (point_idx, offset) in mode.point_offset
            applied_total .+= cross(point_pos_b[point_idx], offset)
        end

        # The same loads at the contour node each one was sampled at, which is the
        # placement the frozen pattern intends.
        pattern_total = zeros(SimFloat, 3)
        worst_arm = 0.0
        column = 0
        for (i, panel) in enumerate(panels)
            xc, yc, _, _ = VortexStepMethod.section_surface(
                panel.section_aero, 0.0, panel.delta)
            assigned = mode.station_point[i]
            f_vsm = Vector(sol.f_body_3D[:, i])
            node_forces = surface_node_forces(mode.traction, column + 1,
                length(assigned), f_vsm, mode.traction_net[:, i])
            applied = zeros(SimFloat, 3)
            pattern = zeros(SimFloat, 3)
            for node in eachindex(assigned)
                column += 1
                applied .+= cross(point_pos_b[assigned[node]], node_forces[node])
                pattern .+= cross(loft_contour_node(panel, xc, yc, node),
                                  node_forces[node])
            end
            pattern_total .+= pattern
            worst_arm = max(worst_arm,
                norm(applied - pattern) / (norm(f_vsm) * panel.chord))
        end

        moment_vsm = Vector(sol.moment)
        mean_chord = sum(p.chord for p in panels) / length(panels)
        lever_scale = norm(vec(sum(sol.f_body_3D, dims=2))) * mean_chord
        @test lever_scale > 0.0
        lumping = norm(applied_total - pattern_total) / lever_scale
        pattern_error = norm(pattern_total - moment_vsm) / lever_scale
        println("  [AeroPressure] moment placement: worst per-panel arm ",
            "$(round(worst_arm; sigdigits=3)) chord; traction pattern vs VSM ",
            "$(round(pattern_error; sigdigits=3)) chord; lumping ",
            "$(round(lumping; sigdigits=3)) chord")
        # The scatter may not move a load farther than the frame guard admits.
        @test worst_arm <= mode.frame_tol_frac
        # Where the traction pattern puts the load centre against where VSM's own
        # Cm does. The fixture's Cp and polars/1.csv are unrelated datasets, so
        # 0.13 chord is their disagreement; negating the traction while keeping
        # the anchored force scores 0.27, which is the margin this bound has.
        @test pattern_error < 0.20
        # Nearest-node quantization, the only part the scatter itself controls.
        @test lumping < 0.05
        # Mz on its own, so a yaw-placement error cannot hide inside the norm.
        @test abs(applied_total[3] - pattern_total[3]) / lever_scale < 0.05
    end

    # The scatter's identities are arithmetic on the frozen params, so they hold
    # exactly rather than approximately: each panel's nodes sum to its net, the
    # pattern pushes the same way as the panel's VSM force (a flipped Cp/cf sign
    # would invert it while leaving the anchored total untouched), and the whole
    # pattern scales with dynamic pressure and nothing else.
    @testset "traction pattern algebra" begin
        sol = wing.vsm_solver.sol
        n_panels = length(wing.vsm_aero.panels)
        column = 0
        for i in 1:n_panels
            node_sum = zeros(SimFloat, 3)
            for _ in eachindex(mode.station_point[i])
                column += 1
                node_sum .+= mode.traction[:, column]
            end
            @test norm(node_sum .- mode.traction_net[:, i]) <=
                1e-12 * (norm(mode.traction_net[:, i]) + 1)
        end

        # Each panel's frozen pattern pushes the same way as its VSM force.
        for i in 1:n_panels
            @test dot(mode.traction_net[:, i], Vector(sol.f_body_3D[:, i])) > 0
        end

        # The pattern scales with dynamic pressure and nothing else.
        traction_before = copy(mode.traction)
        wing.va_b .*= 2.0
        refresh_particle_aero!(mode, wing, sys.points, 2.0 .* va_vals)
        @test isapprox(mode.traction, 4.0 .* traction_before; rtol=1e-6)

        # Returning to the original inflow returns to the original pattern. Only
        # bounded, not exact: the circulation solve is warm-started from the
        # previous gamma (`use_gamma_prev`), so the converged state carries a
        # little history. A mode that latched state would drift by O(1).
        wing.va_b ./= 2.0
        refresh_particle_aero!(mode, wing, sys.points, va_vals)
        @test isapprox(mode.traction, traction_before; rtol=1e-2)
    end

    # Below the cutoff the aero vanishes completely, the constant per-point
    # offsets included — those are what the symbolic scatter adds on top of the
    # live panel force, so leaving them stale keeps applying a load.
    @testset "no stale load below vsm_min_wind" begin
        wing.va_b .= SimFloat[0.1, 0.0, 0.0]
        refresh_particle_aero!(mode, wing, sys.points, 0.0 .* va_vals)
        @test all(iszero, mode.v_ind)
        @test all(iszero, mode.traction)
        @test all(iszero, mode.traction_net)
        @test all(offset -> all(iszero, offset), values(mode.point_offset))
    end
end

# A cambered section is what a canopy is, and what the symmetric fixture above
# cannot exercise: with enough camber both surfaces lie on one side of the chord
# line, so orienting a normal "away from the chord" turns half of them inward.
# Uniform pressure over a closed contour must produce no net force, whatever the
# shape, so the net is a shape-independent check on the normals alone.
@testset "AeroPressure cambered contour closure" begin
    pkg_root = dirname(@__DIR__)
    for camber in (0.0, 0.10)
        tmpdir = mktempdir()
        data_path = joinpath(tmpdir, "2plate_kite")
        cp(joinpath(pkg_root, "data", "2plate_kite"), data_path; force=true)
        write_pressure_fixture(data_path; camber, uniform_cp=1.0)

        set, sys = load_pressure_sys(data_path, AeroPressure())
        wing = sys.wings[1]
        mode = wing.aero
        wing.va_b .= SimFloat[15.0, 0.0, 3.5]
        wing.ω_b .= SimFloat[0.0, 0.0, 0.0]
        wing.vsm_solver.density = 1.225
        va_vals = zeros(SimFloat, 3, length(sys.points))
        for p in sys.points
            @views va_vals[:, p.idx] .= wing.va_b
        end
        refresh_particle_aero!(mode, wing, sys.points, va_vals)

        # Scale by what one surface alone would contribute, so the tolerance is a
        # fraction of the force a flipped half would leave behind.
        dynamic_pressure = 0.5 * wing.vsm_solver.density * dot(wing.va_b, wing.va_b)
        for (i, panel) in enumerate(wing.vsm_aero.panels)
            scale = dynamic_pressure * panel.chord * panel.width
            @test norm(mode.traction_net[:, i]) / scale < 0.05
        end
    end
end
