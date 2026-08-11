# Copyright (c) 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_pressure_aero.jl
# AeroPressure: distribute VSM's per-section load onto structural points using the
# airfoil surface traction (−Cp·n̂ + cf·ŝ) as the pattern, anchored to VSM's
# f_body_3D. Uses the particle 2plate kite with a synthetic section_aero fixture
# (the human-readable .dat/Cp/cf "files-provided" route, no NeuralFoil):
# - construction builds the static surface→point map and passes the frame guard
# - a too-tight frame tolerance errors (misalignment guard)
# - the distributed point forces sum exactly to VSM's total (anchoring)
# - the mode integrates with the compiled RHS (init! + step stay finite)

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

@isdefined(test_init!) || include(joinpath(@__DIR__, "util.jl"))

using Test
using SymbolicAWEModels
using SymbolicAWEModels: VortexStepMethod, refresh_particle_aero!, SimFloat
using KiteUtils
using LinearAlgebra

include(joinpath(@__DIR__, "pressure_fixture.jl"))

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

    @testset "frame guard rejects misalignment" begin
        @test_throws ErrorException load_pressure_sys(
            data_path, AeroPressure(frame_tol_frac=0.01))
    end

    # Exercise the refresh directly (no ODE compile): set the operating point the
    # way update_sys_struct!/sync_aero_density! would, then distribute.
    wing.va_b .= SimFloat[15.0, 0.0, 1.5]
    wing.ω_b .= SimFloat[0.0, 0.0, 0.0]
    wing.vsm_solver.density = 1.225
    va_vals = zeros(SimFloat, 3, length(sys.points))
    for p in sys.points
        @views va_vals[:, p.idx] .= wing.va_b
    end
    refresh_particle_aero!(mode, wing, sys.points, va_vals)

    # Reconstruct the symbolic scatter numerically from the frozen params (same
    # (panel, node) column order as `aero_component`), using VSM `f_body_3D` as the
    # per-panel total (= the live `panel_force` at the solve operating point).
    @testset "correct total force (per VSM panel and total)" begin
        sol = wing.vsm_solver.sol
        n_panels = length(wing.vsm_aero.panels)
        point_force = Dict(p.idx => zeros(SimFloat, 3) for p in wing_pts)
        panel_sum = [zeros(SimFloat, 3) for _ in 1:n_panels]
        column = 0
        for i in 1:n_panels
            assigned = mode.station_point[i]
            n_nodes = length(assigned)
            residual = (Vector(sol.f_body_3D[:, i]) .-
                        mode.traction_net[:, i]) ./ n_nodes
            for node in 1:n_nodes
                column += 1
                node_force = mode.traction[:, column] .+ residual
                point_force[assigned[node]] .+= node_force
                panel_sum[i] .+= node_force
            end
        end
        # Per VSM panel: scattered node forces sum to that panel's VSM force.
        for i in 1:n_panels
            f_vsm = Vector(sol.f_body_3D[:, i])
            @test norm(panel_sum[i] - f_vsm) <= 1e-9 * (norm(f_vsm) + 1)
        end
        # Total: distributed point forces sum to the VSM total force.
        total_pts = sum(values(point_force))
        total_vsm = vec(sum(sol.f_body_3D, dims=2))
        @test norm(total_vsm) > 1.0
        @test norm(total_pts - total_vsm) / norm(total_vsm) < 1e-10
    end

    @testset "compiled RHS stays finite" begin
        sam = SymbolicAWEModel(set, sys)
        test_init!(sam)
        for _ in 1:5
            next_step!(sam; dt=0.01, vsm_interval=1)
        end
        @test all(isfinite, wing.aero_force_b)
        @test all(isfinite, wing.va_b)
        @test all(isfinite, mode.traction)
    end
end
