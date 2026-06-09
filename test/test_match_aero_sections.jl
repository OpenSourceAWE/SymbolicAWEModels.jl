# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_match_aero_sections.jl
# Tests for match_aero_sections_to_structure!:
# verifies geometry alignment and that use_prior_polar
# preserves refined panel polars across section count
# changes (no re-interpolation).

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

using Test
using SymbolicAWEModels
using SymbolicAWEModels: KVec3, VortexStepMethod, WING,
    PARTICLE_DYNAMICS, RIGID_DYNAMICS, SimFloat,
    match_aero_sections_to_structure!
using KiteUtils
using LinearAlgebra

# ── helpers ──────────────────────────────────────────────
pkg_root = dirname(@__DIR__)
src_data = joinpath(pkg_root, "data", "2plate_kite")
tmpdir   = mktempdir()
data_path = joinpath(tmpdir, "2plate_kite")
cp(src_data, data_path; force=true)
set_data_path(data_path)

struc_yaml = joinpath(data_path,
    "rigid_structural_geometry.yaml")
refine_yaml = joinpath(data_path,
    "particle_structural_geometry.yaml")

set = Settings("system.yaml")
set.g_earth = 0.0
vsm_set_path = joinpath(data_path, "vsm_settings.yaml")
vsm_set = VortexStepMethod.VSMSettings(
    vsm_set_path; data_prefix=false)

# ─────────────────────────────────────────────────────────
@testset "match_aero_sections — PARTICLE_DYNAMICS" begin

    @testset "geometry: LE/TE match structural points" begin
        sys = SymbolicAWEModels.load_sys_struct_from_yaml(
            refine_yaml;
            system_name="refine_geom", set, vsm_set)
        wing = sys.wings[1]
        points = sys.points
        vsm_w = wing.vsm_wing

        n_struct = count(
            p -> p.type == WING &&
                 p.wing_idx == wing.idx,
            points) ÷ 2

        # Add extra section to create mismatch
        extra = deepcopy(vsm_w.unrefined_sections[1])
        push!(vsm_w.unrefined_sections, extra)
        vsm_w.n_unrefined_sections =
            Int16(length(vsm_w.unrefined_sections))
        @test vsm_w.n_unrefined_sections != n_struct

        vsm_w.use_prior_polar = true
        wing.wing_segments = nothing

        match_aero_sections_to_structure!(
            wing, points)

        @test vsm_w.n_unrefined_sections == n_struct
        @test !isnothing(wing.wing_segments)
        @test length(wing.wing_segments) == n_struct

        R = wing.R_b_to_c
        origin = wing.pos_cad
        for (i, (le_idx, te_idx)) in
                enumerate(wing.wing_segments)
            sec = vsm_w.unrefined_sections[i]
            le_body = R' * (points[le_idx].pos_cad -
                            origin)
            te_body = R' * (points[te_idx].pos_cad -
                            origin)
            @test isapprox(Vector(sec.LE_point),
                           le_body; atol=1e-10)
            @test isapprox(Vector(sec.TE_point),
                           te_body; atol=1e-10)
        end
    end

    @testset "use_prior_polar preserves refined polars" begin
        sys = SymbolicAWEModels.load_sys_struct_from_yaml(
            refine_yaml;
            system_name="refine_polar", set, vsm_set)
        wing = sys.wings[1]
        points = sys.points
        vsm_w = wing.vsm_wing

        # Baseline refined polars from constructor
        n_refined = length(vsm_w.refined_sections)
        @test n_refined > 0
        @test all(
            s -> !isnothing(s.aero_data),
            vsm_w.refined_sections)
        baseline = [
            deepcopy(sec.aero_data)
            for sec in vsm_w.refined_sections
        ]

        n_struct = count(
            p -> p.type == WING &&
                 p.wing_idx == wing.idx,
            points) ÷ 2

        # Set all unrefined polars to a scaled value so
        # re-interpolation would produce uniform output
        # that differs from baseline
        orig_ad = vsm_w.unrefined_sections[1].aero_data
        @test !isnothing(orig_ad)
        uniform = Tuple(v .* 2.0 for v in orig_ad)
        for sec in vsm_w.unrefined_sections
            sec.aero_data = deepcopy(uniform)
        end

        # Verify non-trivial: some baseline refined
        # polars differ from the uniform unrefined data
        @test any(
            any(!isapprox(baseline[i][k], uniform[k])
                for k in eachindex(uniform))
            for i in eachindex(baseline))

        # Create mismatch: add 2 extra unrefined sections
        for _ in 1:2
            push!(vsm_w.unrefined_sections,
                  deepcopy(vsm_w.unrefined_sections[1]))
        end
        vsm_w.n_unrefined_sections =
            Int16(length(vsm_w.unrefined_sections))
        @test vsm_w.n_unrefined_sections != n_struct

        vsm_w.use_prior_polar = true
        wing.wing_segments = nothing

        match_aero_sections_to_structure!(
            wing, points)

        # Refined count unchanged (n_panels unchanged)
        @test length(vsm_w.refined_sections) == n_refined

        # Refined polars preserved — NOT re-interpolated
        for (i, sec) in
                enumerate(vsm_w.refined_sections)
            @test !isnothing(sec.aero_data)
            for k in eachindex(baseline[i])
                @test sec.aero_data[k] ≈ baseline[i][k]
            end
        end

        # Preserved polars differ from rebuilt unrefined
        # (re-interpolation would have made them equal)
        unrefined_ad = vsm_w.unrefined_sections[1].aero_data
        @test any(
            any(!isapprox(sec.aero_data[k],
                          unrefined_ad[k])
                for k in eachindex(sec.aero_data))
            for sec in vsm_w.refined_sections)

        # Unrefined sections rebuilt to match structure
        @test vsm_w.n_unrefined_sections == n_struct
    end

    @testset "errors when use_prior_polar=false" begin
        sys = SymbolicAWEModels.load_sys_struct_from_yaml(
            refine_yaml;
            system_name="refine_err", set, vsm_set)
        wing = sys.wings[1]
        points = sys.points
        vsm_w = wing.vsm_wing

        # Create mismatch
        push!(vsm_w.unrefined_sections,
              deepcopy(vsm_w.unrefined_sections[1]))
        vsm_w.n_unrefined_sections =
            Int16(length(vsm_w.unrefined_sections))

        vsm_w.use_prior_polar = false
        wing.wing_segments = nothing

        @test_throws ErrorException match_aero_sections_to_structure!(
            wing, points)
    end
end

# ─────────────────────────────────────────────────────────
@testset "match_aero_sections — RIGID_DYNAMICS" begin

    @testset "geometry: LE/TE match structural points" begin
        sys_q = SymbolicAWEModels.load_sys_struct_from_yaml(
            struc_yaml;
            system_name="quat_geom", set, vsm_set,
            dynamics_type=RIGID_DYNAMICS)
        wing = sys_q.wings[1]
        points = sys_q.points
        vsm_w = wing.vsm_wing

        # wing_segments populated by constructor
        @test !isnothing(wing.wing_segments)

        n_struct = length(wing.twist_surface_idxs)
        @test length(wing.wing_segments) == n_struct
        @test vsm_w.n_unrefined_sections == n_struct

        # Verify LE/TE positions match
        R = wing.R_b_to_c
        origin = wing.pos_cad
        for (i, (le_idx, te_idx)) in
                enumerate(wing.wing_segments)
            sec = vsm_w.unrefined_sections[i]
            le_body = R' * (points[le_idx].pos_cad -
                            origin)
            te_body = R' * (points[te_idx].pos_cad -
                            origin)
            @test isapprox(Vector(sec.LE_point),
                           le_body; atol=1e-10)
            @test isapprox(Vector(sec.TE_point),
                           te_body; atol=1e-10)
        end
    end

    @testset "preserve aero when n_twist_surfaces < n_aero" begin
        # When a RIGID_DYNAMICS wing has fewer twist DOFs
        # (twist_surfaces) than aero sections, the OBJ/VSM
        # geometry must stay intact — only the Voronoi
        # partition assigns sections to twist_surfaces.
        sys_q = SymbolicAWEModels.load_sys_struct_from_yaml(
            struc_yaml;
            system_name="quat_preserve", set, vsm_set,
            dynamics_type=RIGID_DYNAMICS)
        wing = sys_q.wings[1]
        points = sys_q.points
        vsm_w = wing.vsm_wing

        n_struct = length(wing.twist_surface_idxs)
        n_refined = length(vsm_w.refined_sections)

        # Add extra section so n_aero > n_twist_surfaces
        extra = deepcopy(vsm_w.unrefined_sections[1])
        push!(vsm_w.unrefined_sections, extra)
        vsm_w.n_unrefined_sections =
            Int16(length(vsm_w.unrefined_sections))
        n_aero_before = vsm_w.n_unrefined_sections
        @test n_aero_before > n_struct

        baseline_unrefined = [
            deepcopy(sec)
            for sec in vsm_w.unrefined_sections
        ]
        baseline_refined = [
            deepcopy(sec.aero_data)
            for sec in vsm_w.refined_sections
        ]

        wing.wing_segments = nothing
        match_aero_sections_to_structure!(
            wing, points;
            twist_surfaces=collect(sys_q.twist_surfaces))

        # Aero geometry untouched
        @test vsm_w.n_unrefined_sections == n_aero_before
        @test length(vsm_w.unrefined_sections) ==
              n_aero_before
        for (i, sec) in
                enumerate(vsm_w.unrefined_sections)
            @test sec.LE_point ≈
                  baseline_unrefined[i].LE_point
            @test sec.TE_point ≈
                  baseline_unrefined[i].TE_point
        end

        # Refined polars untouched (no rebuild ran)
        @test length(vsm_w.refined_sections) == n_refined
        for (i, sec) in
                enumerate(vsm_w.refined_sections)
            for k in eachindex(baseline_refined[i])
                @test sec.aero_data[k] ≈
                      baseline_refined[i][k]
            end
        end

        # wing_segments still populated — one per twist_surface
        @test !isnothing(wing.wing_segments)
        @test length(wing.wing_segments) == n_struct

        # Aero arrays remain twist_surface-count-sized
        n_twist_surfaces = length(wing.twist_surface_idxs)
        @test length(wing.aero_y) == 5 + n_twist_surfaces
        @test length(wing.aero_x) == 6 + n_twist_surfaces
        @test size(wing.aero_jac) ==
              (6 + n_twist_surfaces, 5 + n_twist_surfaces)
    end
end
nothing
