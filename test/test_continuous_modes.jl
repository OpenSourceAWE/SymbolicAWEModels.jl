# Copyright (c) 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_continuous_modes.jl
# The contract every *continuous* VSM mode owes — currently ContinuousAero and
# AeroPressure — on the particle-dynamics 2plate kite. A continuous mode freezes
# the circulation at each VSM refresh and re-derives the per-panel force
# symbolically every RHS step from the live inflow, so all of them owe:
#   - built mesh maps and a finite, nonzero frozen induced velocity
#   - per-section inflow gathered from that section's own bounding struts, never
#     a wing-wide mean (a wing-wide mean annihilates a rigid rotation exactly and
#     leaves the panels without rate damping)
#   - solve-point parity with a full VSM solve! on the same frozen state
#   - live response to the state between refreshes (aerodynamic damping)
#   - a surviving VSM refresh
#
# A new continuous mode inherits all of it by being added to `cases`. Mode-only
# behaviour lives in test_continuous_aero.jl / test_pressure_aero.jl, and the
# contract shared with the non-continuous modes in test_aero_modes.jl.

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

@isdefined(test_init!) || include(joinpath(@__DIR__, "util.jl"))
@isdefined(write_pressure_fixture) ||
    include(joinpath(@__DIR__, "pressure_fixture.jl"))

using Test
using SymbolicAWEModels
using SymbolicAWEModels: VortexStepMethod, aero_inflow_groups, wing_points
using KiteUtils
using LinearAlgebra

"""
    build_continuous_case(case, root)

Fresh data directory, `Settings` and `SystemStructure` for one continuous-mode
case. The mode's spanwise distribution and its own fixture (AeroPressure's
synthetic surface aero) are applied to that copy alone, so no case can disturb
another's geometry or model cache.
"""
function build_continuous_case(case, root)
    data_path = joinpath(root, case.system_name, "2plate_kite")
    mkpath(dirname(data_path))
    cp(joinpath(dirname(@__DIR__), "data", "2plate_kite"), data_path; force=true)
    case.surface_fixture && write_pressure_fixture(data_path)
    set_data_path(data_path)
    set = Settings("system.yaml")
    vsm_set = VortexStepMethod.VSMSettings(
        joinpath(data_path, "vsm_settings.yaml"); data_prefix=false)
    if case.billowing
        for vsm_wing_settings in vsm_set.wings
            vsm_wing_settings.spanwise_panel_distribution =
                VortexStepMethod.BILLOWING
            vsm_wing_settings.billowing_percentage = 8.0
        end
    end
    sys = load_sys_struct_from_yaml(
        joinpath(data_path, "particle_structural_geometry.yaml");
        system_name=case.system_name, set, vsm_set, aero_mode=case.make())
    return set, sys
end

"""
    rigid_rotation_field(points, rate)

World-frame velocity every point of `points` picks up from a rigid rotation
`rate` about their centroid. Its mean over `points` is exactly zero, so an inflow
that averages the whole wing cannot see it while a per-section inflow can.
"""
function rigid_rotation_field(points, rate)
    centroid = sum(point.pos_w for point in points) / length(points)
    return [cross(rate, point.pos_w - centroid) for point in points]
end

"""
    group_means(groups, values)

The weighted mean of `values` each group of [`aero_inflow_groups`](@ref) gathers,
`values` being indexed by the same wing-node column the groups reference.
"""
group_means(groups, values) =
    [sum(weight .* values[column] for (column, weight) in group)
     for group in groups]

@testset "Continuous aero modes" begin
    root = mktempdir()
    # About the wing nodes' centroid: 3 rad/s over the 2plate kite's ~1 m half
    # span is a few m/s against a ~15 m/s inflow.
    roll_rate = [3.0, 0.0, 0.0]
    cases = [
        # Tolerances sit a small factor above the measured values, which the
        # testsets print: parity rel_F 0.024/0.013 and 0.27°/0.11°, roll-rate
        # moment response 0.29/0.22 against 0 for a wing-wide mean inflow.
        (name="ContinuousAero", system_name="continuous_test",
            make=() -> ContinuousAero(), billowing=true, surface_fixture=false,
            parity_rtol=0.04, parity_deg=1.0, roll_response=0.10),
        (name="AeroPressure", system_name="pressure_test",
            make=() -> AeroPressure(), billowing=false, surface_fixture=true,
            parity_rtol=0.04, parity_deg=1.0, roll_response=0.10),
    ]

    for case in cases
        @testset "$(case.name)" begin
            set, sys = build_continuous_case(case, root)
            wing = sys.wings[1]
            mode = wing.aero
            nodes = wing_points(sys, wing)
            n_panels = length(wing.vsm_aero.panels)

            @testset "mesh maps" begin
                @test size(mode.v_ind) == (3, n_panels)
                @test length(mode.section_left_strut) == n_panels + 1
                @test length(mode.section_left_weight) == n_panels + 1
                @test all(0.0 .<= mode.section_left_weight .<= 1.0)
                n_struts = length(wing.vsm_wing.unrefined_sections)
                @test all(1 .<= mode.section_left_strut .<= n_struts - 1)
            end

            @testset "per-section inflow grouping" begin
                groups, section_group = aero_inflow_groups(mode, wing, nodes)
                @test length(groups) == length(mode.section_left_strut)
                @test length(unique(section_group)) == length(groups)
                @test length(unique(groups)) > 1
                @test all(g -> isapprox(sum(last, g), 1.0; atol=1e-10), groups)

                velocity = rigid_rotation_field(nodes, roll_rate)
                largest = maximum(norm, velocity)
                @test largest > 1.0
                @test norm(sum(velocity) / length(velocity)) < 1e-10 * largest
                @test maximum(norm, group_means(groups, velocity)) > 0.3 * largest
            end

            sam = SymbolicAWEModel(set, sys)
            test_init!(sam)

            @testset "frozen induced velocity" begin
                @test all(isfinite, mode.v_ind)
                @test norm(mode.v_ind) > 0.0
            end

            # Sync the symbolic forces (computed with the frozen v_ind from the
            # init refresh) into the struct with a near-zero step.
            next_step!(sam; dt=1e-4, vsm_interval=0)
            force_symbolic = copy(wing.aero_force_b)

            @testset "solve-point parity with full VSM" begin
                # Reference: full nonlinear solve + calc_forces! on the same panel
                # apparent wind the refresh used. Remaining differences: the
                # per-panel va assignment (nearest strut vs interpolated), the
                # corrected-AoA direction triad, and VSM's spanwise aero-center
                # weighting.
                VortexStepMethod.solve!(wing.vsm_solver, wing.vsm_aero)
                force_reference = vec(sum(wing.vsm_solver.sol.f_body_3D, dims=2))
                @test all(isfinite, force_symbolic)
                @test norm(force_reference) > 1.0
                rel_force = norm(force_symbolic - force_reference) /
                    norm(force_reference)
                cos_angle = dot(force_symbolic, force_reference) /
                    (norm(force_symbolic) * norm(force_reference))
                angle_deg = rad2deg(acos(clamp(cos_angle, -1.0, 1.0)))
                println("  [$(case.name)] solve-point parity: ",
                    "rel_F=$(round(rel_force; sigdigits=3)) ",
                    "(tol $(case.parity_rtol)), ",
                    "dir=$(round(angle_deg; digits=3))° ",
                    "(tol $(case.parity_deg)°)")
                @test rel_force < case.parity_rtol
                @test angle_deg < case.parity_deg
            end

            @testset "per-point apparent wind reaches the panels" begin
                # A rigid rotation about the wing nodes' centroid leaves their mean
                # velocity exactly unchanged, so a wing-wide mean inflow would
                # reproduce the previous moment to the last bit. Per-section inflow
                # turns it into a rate-damping moment.
                v_ind_frozen = copy(mode.v_ind)
                moment_before = copy(wing.aero_moment_b)
                @test norm(moment_before) > 1.0
                rotation = rigid_rotation_field(nodes, roll_rate)
                intended = [Vector(point.vel_w) .+ rotation[k]
                            for (k, point) in enumerate(nodes)]
                for (point, velocity) in zip(nodes, rotation)
                    point.vel_w .+= velocity
                end
                init!(sam; reinit_sys=false, lin_vsm=false, prn=false)

                # The rotation only reaches the panels if the re-init actually
                # carried it into the integrator state.
                injected = maximum(norm(Vector(nodes[k].vel_w) .- intended[k])
                                   for k in eachindex(nodes))
                @test injected < 1e-6 * maximum(norm, intended)

                next_step!(sam; dt=1e-6, vsm_interval=0)
                @test mode.v_ind == v_ind_frozen
                @test all(isfinite, wing.aero_moment_b)
                rel_moment = norm(wing.aero_moment_b - moment_before) /
                    norm(moment_before)
                println("  [$(case.name)] roll-rate moment response: ",
                    "rel_M=$(round(rel_moment; sigdigits=3)) ",
                    "(threshold $(case.roll_response)), injection error ",
                    "$(round(injected; sigdigits=3)) m/s")
                @test rel_moment > case.roll_response
            end

            @testset "aerodynamic damping (live forces, frozen circulation)" begin
                # vsm_interval=0: the circulation is never re-solved. AeroDirect
                # would hold the point forces exactly constant; a continuous mode's
                # forces are live functions of the state and must respond.
                force_before = copy(wing.aero_force_b)
                v_ind_before = copy(mode.v_ind)
                for _ in 1:10
                    next_step!(sam; dt=0.01, vsm_interval=0)
                end
                @test mode.v_ind == v_ind_before
                @test norm(wing.aero_force_b - force_before) > 1e-6
                @test all(isfinite, wing.aero_force_b)
            end

            @testset "stepping with VSM refresh" begin
                v_ind_before = copy(mode.v_ind)
                for _ in 1:5
                    next_step!(sam; dt=0.01, vsm_interval=1)
                end
                @test mode.v_ind != v_ind_before
                @test all(isfinite, mode.v_ind)
                @test all(isfinite, wing.aero_force_b)
            end
        end
    end

    rm(root; recursive=true, force=true)
end
nothing
