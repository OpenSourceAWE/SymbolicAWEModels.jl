# Copyright (c) 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Verifies that a RIGID_DYNAMICS wing can have FEWER stations
# than unrefined VSM sections — i.e. one twist DOF drives
# multiple aero sections via the Voronoi partition in
# compute_spatial_station_mapping!.

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

using Test
using SymbolicAWEModels
using SymbolicAWEModels: VortexStepMethod, RIGID_DYNAMICS,
    compute_spatial_station_mapping!,
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

@testset "Multi-section station partition" begin
    sys = SymbolicAWEModels.load_sys_struct_from_yaml(
        struc_yaml; system_name="multi_section",
        set, vsm_set, dynamics_type=RIGID_DYNAMICS)
    wing = sys.wings[1]
    vsm_w = wing.vsm_wing

    # Baseline: 3 stations, 3 unrefined sections (1:1)
    @test length(sys.stations) == 3
    @test vsm_w.n_unrefined_sections == 3
    for station in sys.stations
        @test length(station.unrefined_section_idxs) == 1
    end

    # Every panel belongs to exactly one station, the one nearest its centre, so the
    # mirror-symmetric kite gives its left and right stations the same number of panels.
    n_panels = length(wing.vsm_aero.panels)
    owned = reduce(vcat, [station.panel_idxs for station in sys.stations])
    @test sort(owned) == collect(1:n_panels)
    station_center(station) = sum(wing.R_b_to_c' * (sys.points[i].pos_cad - wing.pos_cad)
                                  for i in station.point_idxs) / length(station.point_idxs)
    offset = [0.0, 0.0, wing.aero_z_offset]
    for station in sys.stations, panel_idx in station.panel_idxs
        corners = wing.vsm_aero.panels[panel_idx].corner_points
        center = vec(sum(corners; dims=2)) / 4 - offset
        @test all(norm(center - station_center(station)) <= norm(center - station_center(other))
                  for other in sys.stations)
    end
    @test length(sys.stations[:left].panel_idxs) == length(sys.stations[:right].panel_idxs)

    # The wing mass sits on its points, so no station takes a share of it. Moved onto
    # the wing body, it is shared by panel area: all of it, the same left and right.
    @test all(station.body_mass == 0 for station in sys.stations)
    for point in sys.points
        SymbolicAWEModels.wing_frame_member(point, wing.idx) && (point.extra_mass = 0.0)
    end
    compute_spatial_station_mapping!(wing, sys.stations, sys.points)
    @test sum(station.body_mass for station in sys.stations) ≈ wing.mass
    @test sys.stations[:left].body_mass ≈ sys.stations[:right].body_mass
    @test sys.stations[:center].body_mass > 0

    # Inject a 4th unrefined section by duplicating an
    # existing one. Now n_stations (3) < n_unrefined (4).
    extra = deepcopy(vsm_w.unrefined_sections[2])
    push!(vsm_w.unrefined_sections, extra)
    vsm_w.n_unrefined_sections = Int16(4)
    wing.wing_segments = nothing

    # match_aero_sections_to_structure! should NOT
    # collapse aero back to 3 sections in this case
    match_aero_sections_to_structure!(
        wing, sys.points; stations=sys.stations)
    @test vsm_w.n_unrefined_sections == 4
    @test length(vsm_w.unrefined_sections) == 4
    @test !isnothing(wing.wing_segments)
    @test length(wing.wing_segments) == 3

    # Re-run partition: every section assigned, no
    # overlaps, every station claims ≥ 1 section.
    compute_spatial_station_mapping!(
        wing, sys.stations, sys.points)
    assigned = Int64[]
    for station in sys.stations
        @test !isempty(station.unrefined_section_idxs)
        append!(assigned, station.unrefined_section_idxs)
    end
    @test sort(assigned) == [1, 2, 3, 4]
    @test length(unique(assigned)) == 4
    @test sort(reduce(vcat, [station.panel_idxs for station in sys.stations])) ==
        collect(1:length(wing.vsm_aero.panels))

    # Wing aero arrays sized by n_stations, not n_unrefined.
    @test length(wing.aero_y) == 5 + 3
    @test length(wing.aero_x) == 6 + 3
    @test size(wing.aero_jac) == (6 + 3, 5 + 3)
end

@testset "Rigid wing mass in total_mass" begin
    sys = SymbolicAWEModels.load_sys_struct_from_yaml(
        struc_yaml; system_name="rigid_total_mass",
        set, vsm_set, dynamics_type=RIGID_DYNAMICS)
    wing = sys.wings[1]

    # The six wing nodes and the kcu ride the wing body, so their 1.6 kg is the
    # wing's mass and counts once, next to the two 0.1 kg bridle points.
    @test wing.mass ≈ 1.6
    @test sys.total_mass ≈ 1.8
    wing.mass = 2.0
    @test sys.total_mass ≈ 2.2
end

@testset "n_stations > n_unrefined errors" begin
    sys = SymbolicAWEModels.load_sys_struct_from_yaml(
        struc_yaml; system_name="too_many_stations",
        set, vsm_set, dynamics_type=RIGID_DYNAMICS)
    wing = sys.wings[1]
    vsm_w = wing.vsm_wing

    # Drop down to 2 unrefined sections while keeping
    # the 3 stations → must error.
    pop!(vsm_w.unrefined_sections)
    vsm_w.n_unrefined_sections = Int16(2)

    @test_throws ErrorException compute_spatial_station_mapping!(
        wing, sys.stations, sys.points)
end
nothing
