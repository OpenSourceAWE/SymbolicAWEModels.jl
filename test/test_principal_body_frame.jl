# SPDX-FileCopyrightText: 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

# test_principal_body_frame.jl
#
# Tests that the principal frame (dynamics) and body
# frame (output) are correctly separated for QUATERNION
# wings with structural reference points.
#
# Key invariants:
# 1. pos_w = com_w + R_b_to_w * pos_b  (rigid body)
# 2. Group le_pos, chord, y_airf in body frame
#    relative to COM (same frame as pos_b)
# 3. R_b_to_p = R_p_to_c' * R_b_to_c (constant body→principal)

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

using Test
using SymbolicAWEModels
using SymbolicAWEModels: KVec3, QUATERNION, WING,
    VortexStepMethod
using LinearAlgebra

@testset "Principal-Body Frame Separation" begin
    pkg_root = dirname(@__DIR__)
    src_data = joinpath(
        pkg_root, "data", "2plate_kite")

    tmpdir = mktempdir()
    data_path = joinpath(tmpdir, "2plate_kite")
    cp(src_data, data_path; force=true)

    set_data_path(data_path)
    set = Settings("system.yaml")

    vsm_set = VortexStepMethod.VSMSettings(
        joinpath(data_path, "vsm_settings.yaml");
        data_prefix=false)

    yaml_path = joinpath(
        data_path, "quat_struc_geometry.yaml")
    sys = load_sys_struct_from_yaml(
        yaml_path;
        system_name="frame_test", set, vsm_set)

    wing = sys.wings[:main_wing]
    @test wing.wing_type == QUATERNION
    @test !isnothing(wing.origin_idx)
    @test !isnothing(wing.z_ref_points)

    # Compute R_p_to_w from Q_p_to_w (same as runtime)
    R_p_to_w = wing.R_p_to_c' * wing.R_b_to_c * wing.R_b_to_w'
    # More directly: wing.R_p_to_w uses Q_p_to_w
    R_p_to_w_from_q = wing.R_p_to_w

    # ---- Test 1: R_p_b is consistent ---- #
    @testset "R_p_b consistency" begin
        # Constant rotation from body to principal
        R_p_b_const = wing.R_p_to_c' * wing.R_b_to_c
        # Same thing from world-frame rotations
        R_p_b_world = R_p_to_w_from_q' * wing.R_b_to_w
        @test R_p_b_const ≈ R_p_b_world atol=1e-10
        # Verify R_p_b is not identity (body != principal)
        # This ensures the test is meaningful
        @test !isapprox(
            R_p_b_const, I(3); atol=0.01)
        println("  R_p_b angle: ",
            round(rad2deg(acos(clamp(
                (tr(R_p_b_const) - 1) / 2,
                -1, 1))); digits=2), "°")
    end

    # ---- Test 2: Rigid body constraint ---- #
    @testset "pos_w = com_w + R_b_to_w * pos_b" begin
        wing_pts = [p for p in sys.points
            if p.type == WING &&
               p.wing_idx == wing.idx]
        @test !isempty(wing_pts)

        for pt in wing_pts
            pos_reconstructed = wing.com_w .+
                wing.R_b_to_w * pt.pos_b
            err = norm(pt.pos_w - pos_reconstructed)
            @test err < 1e-10
            if err > 1e-10
                println("  FAIL point $(pt.name): " *
                    "err = $err")
            end
        end
        println("  All $(length(wing_pts)) WING " *
            "points satisfy rigid body constraint")
    end

    # ---- Test 3: Group geometry in principal frame ---- #
    @testset "Group geometry frame" begin
        wing_pts = [p for p in sys.points
            if p.type == WING &&
               p.wing_idx == wing.idx]

        for group in sys.groups
            # Find LE and TE points in this group
            le_pt = nothing
            te_pt = nothing
            for pt_idx in group.point_idxs
                pt = sys.points[pt_idx]
                name_str = string(pt.name)
                if occursin("le", name_str)
                    le_pt = pt
                elseif occursin("te", name_str)
                    te_pt = pt
                end
            end
            isnothing(le_pt) && continue
            isnothing(te_pt) && continue

            # chord_b = pos_b - le_pos (from point_eqs)
            chord_from_pos = te_pt.pos_b -
                group.le_pos
            # chord_from_pos should be parallel to
            # group.chord (same direction, maybe
            # different magnitude)
            chord_dir = normalize(group.chord)
            pos_dir = normalize(chord_from_pos)
            dot_val = abs(dot(chord_dir, pos_dir))
            @test dot_val > 0.99
            if dot_val < 0.99
                println("  FAIL group $(group.idx): " *
                    "chord misaligned, dot=$dot_val")
            end

            # le_pos should be close to the LE point's
            # pos_b (since LE is at the leading edge)
            le_err = norm(le_pt.pos_b - group.le_pos)
            # LE point may not exactly match le_pos
            # (le_pos comes from panel center, pos_b
            # from point mass position) but should be
            # in the same ballpark
            @test le_err < 1.0
        end
        println("  All $(length(sys.groups)) groups " *
            "have consistent geometry")
    end

    # ---- Test 4: com_w from body origin ---- #
    @testset "com_w = pos_w + R_b_to_w * com_offset_b" begin
        com_reconstructed = wing.pos_w .+
            wing.R_b_to_w * wing.com_offset_b
        @test wing.com_w ≈ com_reconstructed atol=1e-10
    end

    # ---- Test 5: Body ≠ principal origin ---- #
    @testset "COM offset is nonzero" begin
        # For the 2plate_kite with kcu at origin and
        # COM at wing centroid, offset should be nonzero
        @test norm(wing.com_offset_b) > 0.1
        println("  com_offset_b = ",
            round.(wing.com_offset_b; digits=3))
    end

    rm(tmpdir; recursive=true)
end
nothing
