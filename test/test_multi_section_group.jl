# Copyright (c) 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Verifies that a RIGID_DYNAMICS wing can have FEWER twist_surfaces
# than unrefined VSM sections — i.e. one twist DOF drives
# multiple aero sections via the Voronoi partition in
# compute_spatial_twist_surface_mapping!.

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

using Test
using SymbolicAWEModels
using SymbolicAWEModels: VortexStepMethod, RIGID_DYNAMICS,
    compute_spatial_twist_surface_mapping!,
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
    "rigid_structural_geometry.yaml")

set = Settings("system.yaml")
vsm_set_path = joinpath(data_path, "vsm_settings.yaml")
vsm_set = VortexStepMethod.VSMSettings(
    vsm_set_path; data_prefix=false)

@testset "Multi-section twist_surface partition" begin
    sys = SymbolicAWEModels.load_sys_struct_from_yaml(
        struc_yaml; system_name="multi_section",
        set, vsm_set, dynamics_type=RIGID_DYNAMICS)
    wing = sys.wings[1]
    vsm_w = wing.vsm_wing

    # Baseline: 3 twist_surfaces, 3 unrefined sections (1:1)
    @test length(sys.twist_surfaces) == 3
    @test vsm_w.n_unrefined_sections == 3
    for twist_surface in sys.twist_surfaces
        @test length(twist_surface.unrefined_section_idxs) == 1
    end

    # Inject a 4th unrefined section by duplicating an
    # existing one. Now n_twist_surfaces (3) < n_unrefined (4).
    extra = deepcopy(vsm_w.unrefined_sections[2])
    push!(vsm_w.unrefined_sections, extra)
    vsm_w.n_unrefined_sections = Int16(4)
    wing.wing_segments = nothing

    # match_aero_sections_to_structure! should NOT
    # collapse aero back to 3 sections in this case
    match_aero_sections_to_structure!(
        wing, sys.points; twist_surfaces=sys.twist_surfaces)
    @test vsm_w.n_unrefined_sections == 4
    @test length(vsm_w.unrefined_sections) == 4
    @test !isnothing(wing.wing_segments)
    @test length(wing.wing_segments) == 3

    # Re-run partition: every section assigned, no
    # overlaps, every twist_surface claims ≥ 1 section.
    compute_spatial_twist_surface_mapping!(
        wing, sys.twist_surfaces, sys.points)
    assigned = Int64[]
    for twist_surface in sys.twist_surfaces
        @test !isempty(twist_surface.unrefined_section_idxs)
        append!(assigned, twist_surface.unrefined_section_idxs)
    end
    @test sort(assigned) == [1, 2, 3, 4]
    @test length(unique(assigned)) == 4

    # Wing aero arrays sized by n_twist_surfaces, not n_unrefined.
    @test length(wing.aero_y) == 5 + 3
    @test length(wing.aero_x) == 6 + 3
    @test size(wing.aero_jac) == (6 + 3, 5 + 3)
end

@testset "n_twist_surfaces > n_unrefined errors" begin
    sys = SymbolicAWEModels.load_sys_struct_from_yaml(
        struc_yaml; system_name="too_many_twist_surfaces",
        set, vsm_set, dynamics_type=RIGID_DYNAMICS)
    wing = sys.wings[1]
    vsm_w = wing.vsm_wing

    # Drop down to 2 unrefined sections while keeping
    # the 3 twist_surfaces → must error.
    pop!(vsm_w.unrefined_sections)
    vsm_w.n_unrefined_sections = Int16(2)

    @test_throws ErrorException compute_spatial_twist_surface_mapping!(
        wing, sys.twist_surfaces, sys.points)
end
nothing
