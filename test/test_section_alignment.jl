# Copyright (c) 2025 Bart van de Lint, Jelle Poland
# SPDX-License-Identifier: MPL-2.0

# test_section_alignment.jl
# Check that VSM unrefined section LE/TE points align with
# structural WING-type points in both body and world frames.

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

using Test
using SymbolicAWEModels
using SymbolicAWEModels: KVec3, VortexStepMethod, WING,
    QUATERNION
using KiteUtils
using LinearAlgebra

# Use the 2plate_kite quat geometry
pkg_root = dirname(@__DIR__)
set_data_path(joinpath(pkg_root, "data", "2plate_kite"))

struc_yaml = joinpath(
    get_data_path(), "quat_struc_geometry.yaml")

set = Settings("system.yaml")
set.g_earth = 0.0
vsm_set_path = joinpath(
    get_data_path(), "vsm_settings.yaml")
vsm_set = VortexStepMethod.VSMSettings(vsm_set_path; data_prefix=false)

sys = SymbolicAWEModels.load_sys_struct_from_yaml(
    struc_yaml; system_name="2plate_kite",
    set=set, vsm_set=vsm_set)

wing = sys.wings[1]
vsm_wing = wing.vsm_wing
points = sys.points

# Collect WING-type points for this wing, sorted by y
# Exclude origin point (e.g. KCU) which is not an LE/TE
wing_pts = filter(
    p -> p.type == WING && p.wing_idx == wing.idx &&
         p.idx != wing.origin_idx,
    points)
sort!(wing_pts; by=p -> p.pos_cad[2], rev=true)

# Identify LE/TE pairs (LE has smaller x in CAD frame)
n_sections = length(wing_pts) ÷ 2
pairs = Vector{Tuple{typeof(wing_pts[1]),
                      typeof(wing_pts[1])}}()
for i in 1:2:length(wing_pts)
    p1, p2 = wing_pts[i], wing_pts[i+1]
    if p1.pos_cad[1] < p2.pos_cad[1]
        push!(pairs, (p1, p2))  # (LE, TE)
    else
        push!(pairs, (p2, p1))
    end
end

com = wing.pos_cad
R = wing.R_b_to_c

println("== Wing info ==")
println("  COM = ", round.(com; digits=4))
println("  R_b_to_c = ", round.(R; digits=4))
println("  n_unrefined = ",
    vsm_wing.n_unrefined_sections)
println("  n_wing_pts = ", length(wing_pts))
println()

@testset "Body-frame section alignment" begin
    n = vsm_wing.n_unrefined_sections
    @test n == n_sections

    for i in 1:n
        sec = vsm_wing.unrefined_sections[i]
        le_pt, te_pt = pairs[i]

        # Body-frame positions of structural points
        le_body = R' * (le_pt.pos_cad - com)
        te_body = R' * (te_pt.pos_cad - com)

        sec_le = Vector(sec.LE_point)
        sec_te = Vector(sec.TE_point)

        println("Section $i (body frame):")
        println("  struct LE = ",
            round.(le_body; digits=4))
        println("  vsm    LE = ",
            round.(sec_le; digits=4))
        println("  LE err = ",
            round.(sec_le - le_body; digits=6))
        println("  TE err = ",
            round.(sec_te - te_body; digits=6))

        @test isapprox(sec_le, le_body; atol=1e-10)
        @test isapprox(sec_te, te_body; atol=1e-10)
    end
end

# Apply transforms to get world-frame positions
println("\n== Applying transforms ==")
SymbolicAWEModels.reinit!(sys, set)

println("  R_b_to_w = ", round.(wing.R_b_to_w; digits=4))
println("  pos_w = ", round.(wing.pos_w; digits=4))
println()

@testset "World-frame panel alignment" begin
    R_bw = wing.R_b_to_w
    T_bw = wing.pos_w
    panels = wing.vsm_aero.panels

    for i in 1:length(pairs)
        le_pt, te_pt = pairs[i]

        # Structural points in world frame (set by
        # reinit! via transforms)
        le_w = Vector(le_pt.pos_w)
        te_w = Vector(te_pt.pos_w)

        # Find the panel whose center is closest to
        # the midpoint of this LE-TE pair in world
        mid_w = (le_w + te_w) / 2
        best_dist = Inf
        best_panel = panels[1]
        for panel in panels
            pc_body = (Vector(panel.LE_point_1) +
                       Vector(panel.LE_point_2)) / 2
            pc_w = R_bw * pc_body + T_bw
            d = norm(pc_w - mid_w)
            if d < best_dist
                best_dist = d
                best_panel = panel
            end
        end

        # Panel corners → world frame
        p_le1_w = R_bw * Vector(best_panel.LE_point_1) +
                  T_bw
        p_le2_w = R_bw * Vector(best_panel.LE_point_2) +
                  T_bw
        p_te1_w = R_bw * Vector(best_panel.TE_point_1) +
                  T_bw
        p_te2_w = R_bw * Vector(best_panel.TE_point_2) +
                  T_bw

        # For unrefined panels, LE_point_1 and
        # LE_point_2 are from adjacent sections.
        # Check that one of them matches LE and one
        # matches TE (or that the average matches).
        le_panel_w = (p_le1_w + p_le2_w) / 2
        te_panel_w = (p_te1_w + p_te2_w) / 2

        # Check individual corner alignment with
        # structural points
        println("Section $i (world frame):")
        println("  struct LE = ",
            round.(le_w; digits=4))
        println("  struct TE = ",
            round.(te_w; digits=4))
        println("  panel LE1 = ",
            round.(p_le1_w; digits=4))
        println("  panel LE2 = ",
            round.(p_le2_w; digits=4))
        println("  panel TE1 = ",
            round.(p_te1_w; digits=4))
        println("  panel TE2 = ",
            round.(p_te2_w; digits=4))

        # At least one panel corner should match each
        # structural point (for unrefined, section i
        # becomes panel boundary)
        le_match = min(
            norm(p_le1_w - le_w),
            norm(p_le2_w - le_w))
        te_match = min(
            norm(p_te1_w - te_w),
            norm(p_te2_w - te_w))
        println("  LE closest dist = ",
            round(le_match; digits=6))
        println("  TE closest dist = ",
            round(te_match; digits=6))
        println()

        # Panel corner should coincide with a
        # structural point
        @test le_match < 0.01
        @test te_match < 0.01
    end
end

# Also check unrefined sections → world directly
@testset "World-frame unrefined section alignment" begin
    R_bw = wing.R_b_to_w
    T_bw = wing.pos_w
    n = vsm_wing.n_unrefined_sections

    for i in 1:n
        sec = vsm_wing.unrefined_sections[i]
        le_pt, te_pt = pairs[i]

        # Section LE/TE → world
        sec_le_w = R_bw * Vector(sec.LE_point) + T_bw
        sec_te_w = R_bw * Vector(sec.TE_point) + T_bw

        # Structural points in world frame
        le_w = Vector(le_pt.pos_w)
        te_w = Vector(te_pt.pos_w)

        le_err = norm(sec_le_w - le_w)
        te_err = norm(sec_te_w - te_w)

        println("Section $i (world via R_b_to_w):")
        println("  struct LE_w = ",
            round.(le_w; digits=4))
        println("  sec    LE_w = ",
            round.(sec_le_w; digits=4))
        println("  LE err = ", round(le_err; digits=6))
        println("  TE err = ", round(te_err; digits=6))

        @test le_err < 1e-10
        @test te_err < 1e-10
    end
end
nothing
