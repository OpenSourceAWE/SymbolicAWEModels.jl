# Copyright (c) 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_continuous_aero.jl
# ContinuousAero (frozen circulation, live symbolic forces) on the
# particle-dynamics 2plate kite:
# - mesh maps and frozen induced-velocity buffer are built and sane
# - solve-point parity: the symbolic per-refined-panel force sum matches a
#   full VSM solve!+calc_forces! on the same frozen state
# - aerodynamic damping: with vsm_interval=0 (circulation frozen forever) the
#   forces still respond to the moving state, unlike AeroDirect's frozen forces
# - the generated RHS stays allocation-free (test_init!)

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

@isdefined(test_init!) || include(joinpath(@__DIR__, "util.jl"))

using Test
using SymbolicAWEModels
using SymbolicAWEModels: VortexStepMethod, WING
using KiteUtils
using LinearAlgebra

@testset "ContinuousAero" begin
    pkg_root = dirname(@__DIR__)
    src_data_path = joinpath(pkg_root, "data", "2plate_kite")
    tmpdir = mktempdir()
    data_path = joinpath(tmpdir, "2plate_kite")
    cp(src_data_path, data_path; force=true)
    set_data_path(data_path)
    set = Settings("system.yaml")

    vsm_settings_path = joinpath(data_path, "vsm_settings.yaml")
    vsm_set = VortexStepMethod.VSMSettings(
        vsm_settings_path; data_prefix=false)
    # ContinuousAero requires the BILLOWING spanwise distribution.
    for vsm_wing_settings in vsm_set.wings
        vsm_wing_settings.spanwise_panel_distribution =
            VortexStepMethod.BILLOWING
        vsm_wing_settings.billowing_percentage = 8.0
    end
    particle_yaml = joinpath(data_path,
        "particle_structural_geometry.yaml")

    sys = load_sys_struct_from_yaml(particle_yaml;
        system_name="continuous_test", set, vsm_set,
        aero_mode=ContinuousAero())
    wing = sys.wings[1]
    mode = wing.aero
    n_panels = length(wing.vsm_aero.panels)

    @testset "mesh maps" begin
        @test mode isa ContinuousAero
        @test size(mode.v_ind) == (3, n_panels)
        @test length(mode.section_left_strut) == n_panels + 1
        @test length(mode.section_left_weight) == n_panels + 1
        @test all(0.0 .<= mode.section_left_weight .<= 1.0)
        n_struts = length(wing.vsm_wing.unrefined_sections)
        @test all(1 .<= mode.section_left_strut .<= n_struts - 1)
    end

    sam = SymbolicAWEModel(set, sys)
    test_init!(sam)

    @testset "frozen induced velocity" begin
        @test all(isfinite, mode.v_ind)
        @test norm(mode.v_ind) > 0.0
    end

    # Sync the symbolic forces (computed with the frozen v_ind from the init
    # refresh) into the struct with a near-zero step.
    next_step!(sam; dt=1e-4, vsm_interval=0)
    force_symbolic = copy(wing.aero_force_b)

    @testset "solve-point parity with full VSM" begin
        # Reference: full nonlinear solve + calc_forces! on the same panel
        # apparent wind the refresh used. Remaining differences: the per-panel
        # va assignment (nearest strut vs interpolated), the corrected-AoA
        # direction triad, and VSM's spanwise aero-center weighting.
        VortexStepMethod.solve!(wing.vsm_solver, wing.vsm_aero)
        force_reference = vec(sum(wing.vsm_solver.sol.f_body_3D, dims=2))
        @test all(isfinite, force_symbolic)
        @test norm(force_reference) > 1.0
        @test norm(force_symbolic - force_reference) /
              norm(force_reference) < 0.07
        cos_angle = dot(force_symbolic, force_reference) /
            (norm(force_symbolic) * norm(force_reference))
        @test cos_angle > cos(deg2rad(3))
    end

    @testset "aerodynamic damping (live forces, frozen circulation)" begin
        # vsm_interval=0: the circulation is never re-solved. AeroDirect would
        # hold the point forces exactly constant; ContinuousAero's forces are
        # live functions of the state and must respond to the motion.
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
