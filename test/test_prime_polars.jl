# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

# test_prime_polars.jl
# Tests for prime_polars_then_lock_unrefined_to_structure!:
# verifies geometry, unrefined polar propagation, and
# refined polar interpolation for both REFINE and QUATERNION
# wing types.

using Test
using SymbolicAWEModels
using SymbolicAWEModels: KVec3, VortexStepMethod, WING,
    REFINE, QUATERNION, SimFloat,
    prime_polars_then_lock_unrefined_to_structure!
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
    "quat_struc_geometry.yaml")
refine_yaml = joinpath(data_path,
    "refine_struc_geometry.yaml")
aero_yaml = joinpath(data_path, "aero_geometry.yaml")

SymbolicAWEModels.update_aero_yaml_from_struc_yaml!(
    struc_yaml, aero_yaml)

set = Settings("system.yaml")
set.g_earth = 0.0
vsm_set_path = joinpath(data_path, "vsm_settings.yaml")
vsm_set = VortexStepMethod.VSMSettings(
    vsm_set_path; data_prefix=false)

# ── mock polar helpers ───────────────────────────────────
alpha_mock = collect(-5.0:1.0:5.0)  # 11 points
n_alpha = length(alpha_mock)

"""Create recognizable mock polar for section `i`:
cl=i, cd=10i, cm=100i (constant across alpha)."""
function mock_polar(i)
    cl = fill(Float64(i), n_alpha)
    cd = fill(Float64(10i), n_alpha)
    cm = fill(Float64(100i), n_alpha)
    return (copy(alpha_mock), cl, cd, cm)
end

# ─────────────────────────────────────────────────────────
@testset "prime_polars — REFINE" begin

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

        prime_polars_then_lock_unrefined_to_structure!(
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

    # Shared setup for polar tests: 5 aero sections with
    # mock polars, 3 structural sections, cleared refined
    # aero_data so refine! freshly interpolates.
    sys = SymbolicAWEModels.load_sys_struct_from_yaml(
        refine_yaml;
        system_name="refine_polar", set, vsm_set)
    wing = sys.wings[1]
    points = sys.points
    vsm_w = wing.vsm_wing

    n_struct = count(
        p -> p.type == WING &&
             p.wing_idx == wing.idx,
        points) ÷ 2

    # Add 2 extra sections: 5 aero vs 3 struct
    for _ in 1:2
        push!(vsm_w.unrefined_sections,
              deepcopy(vsm_w.unrefined_sections[1]))
    end
    n_templates = length(vsm_w.unrefined_sections)
    vsm_w.n_unrefined_sections = Int16(n_templates)

    # Stamp recognizable mock polars on all 5 sections
    for i in 1:n_templates
        sec = vsm_w.unrefined_sections[i]
        sec.aero_model = VortexStepMethod.POLAR_VECTORS
        sec.aero_data = mock_polar(i)
    end

    # Clear refined section aero_data so refine! freshly
    # interpolates (makes _can_reuse_prior_refined_polar_data
    # return false).
    for sec in vsm_w.refined_sections
        sec.aero_data = nothing
        sec.aero_model = VortexStepMethod.POLAR_VECTORS
    end

    vsm_w.use_prior_polar = true
    wing.wing_segments = nothing

    prime_polars_then_lock_unrefined_to_structure!(
        wing, points)

    @testset "polars correct on unrefined sections" begin
        @test vsm_w.n_unrefined_sections == n_struct

        for i in 1:n_struct
            sec = vsm_w.unrefined_sections[i]
            @test sec.aero_model ==
                  VortexStepMethod.POLAR_VECTORS

            # Same index mapping as vsm_refine.jl
            tidx = n_struct == 1 ? 1 :
                round(Int, 1 + (i - 1) *
                    (n_templates - 1) /
                    (n_struct - 1))
            expected = mock_polar(tidx)

            @test sec.aero_data[1] ≈ expected[1]
            @test sec.aero_data[2] ≈ expected[2]
            @test sec.aero_data[3] ≈ expected[3]
            @test sec.aero_data[4] ≈ expected[4]
        end
    end

    @testset "polars correct on refined sections" begin
        refined = vsm_w.refined_sections
        n_panels = Int(vsm_w.n_panels)
        @test length(refined) == n_panels + 1

        # Unrefined cl values (constant per section)
        unrefined_cls = [
            Float64(
                vsm_w.unrefined_sections[i].aero_data[2][1]
            )
            for i in 1:n_struct
        ]
        cl_lo = minimum(unrefined_cls)
        cl_hi = maximum(unrefined_cls)

        for sec in refined
            @test sec.aero_model ==
                  VortexStepMethod.POLAR_VECTORS
            @test !isnothing(sec.aero_data)

            # Every refined cl value should be a weighted
            # interpolation of neighboring unrefined cls,
            # so it must lie within [cl_lo, cl_hi].
            cl_val = sec.aero_data[2][1]
            @test cl_lo - 1e-10 ≤ cl_val ≤ cl_hi + 1e-10
        end
    end
end

# ─────────────────────────────────────────────────────────
@testset "prime_polars — QUATERNION" begin

    @testset "geometry: LE/TE match structural points" begin
        sys_q = SymbolicAWEModels.load_sys_struct_from_yaml(
            struc_yaml;
            system_name="quat_geom", set, vsm_set,
            wing_type=QUATERNION)
        wing = sys_q.wings[1]
        points = sys_q.points
        vsm_w = wing.vsm_wing

        # wing_segments populated by constructor
        @test !isnothing(wing.wing_segments)

        n_struct = length(wing.group_idxs)
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

    @testset "vsm resize after count mismatch" begin
        sys_q = SymbolicAWEModels.load_sys_struct_from_yaml(
            struc_yaml;
            system_name="quat_resize", set, vsm_set,
            wing_type=QUATERNION)
        wing = sys_q.wings[1]
        points = sys_q.points
        vsm_w = wing.vsm_wing

        n_struct = length(wing.group_idxs)

        # Add extra section to create mismatch
        extra = deepcopy(vsm_w.unrefined_sections[1])
        push!(vsm_w.unrefined_sections, extra)
        vsm_w.n_unrefined_sections =
            Int16(length(vsm_w.unrefined_sections))
        @test vsm_w.n_unrefined_sections != n_struct

        vsm_w.use_prior_polar = true
        wing.wing_segments = nothing

        prime_polars_then_lock_unrefined_to_structure!(
            wing, points;
            groups=collect(sys_q.groups))

        @test vsm_w.n_unrefined_sections == n_struct

        # Verify linearization vectors resized
        ny = 3 + n_struct + 3
        nx = 3 + 3 + n_struct
        @test length(wing.vsm_y) == ny
        @test length(wing.vsm_x) == nx
        @test size(wing.vsm_jac) == (nx, ny)
    end
end
