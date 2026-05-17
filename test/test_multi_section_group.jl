# Copyright (c) 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Verifies that a QUATERNION wing can have FEWER groups
# than unrefined VSM sections — i.e. one twist DOF drives
# multiple aero sections via the Voronoi partition in
# compute_spatial_group_mapping!.

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

using Test
using SymbolicAWEModels
using SymbolicAWEModels: VortexStepMethod, QUATERNION,
    compute_spatial_group_mapping!,
    match_aero_sections_to_structure!
using KiteUtils
using LinearAlgebra

pkg_root = dirname(@__DIR__)
src_data = joinpath(pkg_root, "data", "2plate_kite")
tmpdir = mktempdir()
data_path = joinpath(tmpdir, "2plate_kite")
cp(src_data, data_path; force=true)
set_data_path(data_path)

struc_yaml = joinpath(data_path,
    "quat_struc_geometry.yaml")

set = Settings("system.yaml")
vsm_set_path = joinpath(data_path, "vsm_settings.yaml")
vsm_set = VortexStepMethod.VSMSettings(
    vsm_set_path; data_prefix=false)

@testset "Multi-section group partition" begin
    sys = SymbolicAWEModels.load_sys_struct_from_yaml(
        struc_yaml; system_name="multi_section",
        set, vsm_set, wing_type=QUATERNION)
    wing = sys.wings[1]
    vsm_w = wing.vsm_wing

    # Baseline: 3 groups, 3 unrefined sections (1:1)
    @test length(sys.groups) == 3
    @test vsm_w.n_unrefined_sections == 3
    for group in sys.groups
        @test length(group.unrefined_section_idxs) == 1
    end

    # Inject a 4th unrefined section by duplicating an
    # existing one. Now n_groups (3) < n_unrefined (4).
    extra = deepcopy(vsm_w.unrefined_sections[2])
    push!(vsm_w.unrefined_sections, extra)
    vsm_w.n_unrefined_sections = Int16(4)
    wing.wing_segments = nothing

    # match_aero_sections_to_structure! should NOT
    # collapse aero back to 3 sections in this case
    match_aero_sections_to_structure!(
        wing, sys.points; groups=sys.groups)
    @test vsm_w.n_unrefined_sections == 4
    @test length(vsm_w.unrefined_sections) == 4
    @test !isnothing(wing.wing_segments)
    @test length(wing.wing_segments) == 3

    # Re-run partition: every section assigned, no
    # overlaps, every group claims ≥ 1 section.
    compute_spatial_group_mapping!(
        wing, sys.groups, sys.points)
    assigned = Int64[]
    for group in sys.groups
        @test !isempty(group.unrefined_section_idxs)
        append!(assigned, group.unrefined_section_idxs)
    end
    @test sort(assigned) == [1, 2, 3, 4]
    @test length(unique(assigned)) == 4

    # Wing aero arrays sized by n_groups, not n_unrefined.
    @test length(wing.aero_y) == 5 + 3
    @test length(wing.aero_x) == 6 + 3
    @test size(wing.aero_jac) == (6 + 3, 5 + 3)
end

@testset "n_groups > n_unrefined errors" begin
    sys = SymbolicAWEModels.load_sys_struct_from_yaml(
        struc_yaml; system_name="too_many_groups",
        set, vsm_set, wing_type=QUATERNION)
    wing = sys.wings[1]
    vsm_w = wing.vsm_wing

    # Drop down to 2 unrefined sections while keeping
    # the 3 groups → must error.
    pop!(vsm_w.unrefined_sections)
    vsm_w.n_unrefined_sections = Int16(2)

    @test_throws ErrorException compute_spatial_group_mapping!(
        wing, sys.groups, sys.points)
end
nothing
