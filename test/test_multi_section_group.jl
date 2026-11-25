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
    match_aero_sections_to_structure!, wing_body, initial_body_position
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

"""The run value `quantity` of station `id`, an index or a name, of `sam`."""
station_value(sam, id, quantity) = run_value(sam, :spanwise_stations, id, quantity)

@testset "Multi-section station partition" begin
    sys, _ = load_authoring(
        struc_yaml; name="multi_section", set)
    sam = SymbolicAWEModel(set, sys; vsm_set)
    mode, stations = sam.aero[1], eachindex(sam.definition.spanwise_stations)
    body_idx = wing_body(sam, 1)
    vsm_w = mode.vsm_wing
    sections(s) = station_value(sam, s, :unrefined_section_idxs)
    panels(s) = station_value(sam, s, :panel_idxs)
    body_mass(s) = station_value(sam, s, :body_mass)

    # Baseline: 3 stations, 3 unrefined sections (1:1)
    @test length(stations) == 3
    @test vsm_w.n_unrefined_sections == 3
    for s in stations
        @test length(sections(s)) == 1
    end

    # Every panel belongs to exactly one station, the one nearest its centre, a tie within
    # 1e-9 relative going outboard, so the mirror-symmetric kite gives its left and right
    # stations the same number of panels.
    n_panels = length(mode.vsm_aero.panels)
    owned = reduce(vcat, [panels(s) for s in stations])
    @test sort(owned) == collect(1:n_panels)
    function station_center(s)
        points = sam.definition.spanwise_stations[s].points
        return sum(initial_body_position(sam, body_idx, sam.definition.points[i].pos_ENU)
                   for i in points) / length(points)
    end
    offset = [0.0, 0.0, mode.aero_z_offset]
    for s in stations, panel_idx in panels(s)
        corners = mode.vsm_aero.panels[panel_idx].corner_points
        center = vec(sum(corners; dims=2)) / 4 - offset
        distance(other) = norm(center - station_center(other))
        @test all(distance(s) <= distance(other) * (1 + 1e-9) for other in stations)
    end
    @test length(panels(:left)) == length(panels(:right))

    # The wing mass sits on its points, so no station takes a share of it. Given to
    # the wing body as its own mass, it is shared by panel area: all of it, the same
    # left and right.
    @test all(body_mass(s) == 0 for s in stations)
    sam.definition.bodies[body_idx].extra_mass = 1.6
    compute_spatial_station_mapping!(sam, 1)
    @test sum(body_mass(s) for s in stations) ≈
        run_value(sam, :bodies, body_idx, :extra_mass)
    @test body_mass(:left) ≈ body_mass(:right)
    @test body_mass(:center) > 0

    # Inject a 4th unrefined section by duplicating an
    # existing one. Now n_stations (3) < n_unrefined (4).
    extra = deepcopy(vsm_w.unrefined_sections[2])
    push!(vsm_w.unrefined_sections, extra)
    vsm_w.n_unrefined_sections = Int16(4)
    mode.wing_segments = nothing

    # match_aero_sections_to_structure! should NOT
    # collapse aero back to 3 sections in this case
    match_aero_sections_to_structure!(mode, sam, 1)
    @test vsm_w.n_unrefined_sections == 4
    @test length(vsm_w.unrefined_sections) == 4
    @test !isnothing(mode.wing_segments)
    @test length(mode.wing_segments) == 3

    # Re-run partition: every section assigned, no
    # overlaps, every station claims ≥ 1 section.
    compute_spatial_station_mapping!(sam, 1)
    assigned = Int64[]
    for s in stations
        @test !isempty(sections(s))
        append!(assigned, sections(s))
    end
    @test sort(assigned) == [1, 2, 3, 4]
    @test length(unique(assigned)) == 4
    @test sort(reduce(vcat, [panels(s) for s in stations])) ==
        collect(1:length(mode.vsm_aero.panels))

    # Wing aero arrays sized by n_stations, not n_unrefined.
    @test length(mode.aero_y) == 5 + 3
    @test length(mode.aero_x) == 6 + 3
    @test size(mode.aero_jac) == (6 + 3, 5 + 3)
end

"""The total mass of the points riding the body of wing 1 of `sam`: extra_mass plus their
segment halves."""
function riding_mass(sam)
    body_idx = wing_body(sam, 1)
    riders = [i for (i, point) in enumerate(sam.definition.points)
              if !isnothing(point.body) &&
                 SymbolicAWEModels.item_index(sam, :bodies, point.body) == body_idx]
    mass = sum(something(sam.definition.points[i].extra_mass, 0.0) for i in riders)
    segment(k, quantity) = run_value(sam, :segments, k, quantity)
    for (k, row) in enumerate(sam.definition.segments), i in row.points
        i in riders && (mass += SymbolicAWEModels.segment_half_mass(
            segment(k, :rest_length), segment(k, :diameter), segment(k, :density)))
    end
    return mass
end

@testset "Rigid wing total_mass is its extra_mass plus its riders" begin
    sys, _ = load_authoring(
        struc_yaml; name="rigid_total_mass", set)
    sam = SymbolicAWEModel(set, sys; vsm_set)

    # The six wing nodes and the kcu carry 1.6 kg of extra_mass; it stays on them
    # and reaches the wing once, with the halves of the bridles and tether they hold.
    @test run_value(sam, :wings, 1, :extra_mass) == 0
    @test riding_mass(sam) > 1.6
    @test run_value(sam, :wings, 1, :total_mass) ≈ riding_mass(sam)

    # Mass on the wing and on its points both count, with no warning.
    wing_yaml = replace(read(struc_yaml, String),
        "      aero_z_offset: 0.0" => "      aero_z_offset: 0.0\n      extra_mass: 2.0")
    both_yaml = joinpath(data_path, "rigid_both_masses.yaml")
    write(both_yaml, wing_yaml)
    load_both() = load_authoring(
        both_yaml; name="rigid_total_mass", set)
    both, _ = @test_logs min_level=Base.CoreLogging.Warn load_both()
    both_sam = SymbolicAWEModel(set, both; vsm_set)
    @test run_value(both_sam, :wings, 1, :extra_mass) == 2.0
    @test run_value(both_sam, :wings, 1, :total_mass) ≈ 2.0 + riding_mass(both_sam)
end

@testset "n_stations > n_unrefined errors" begin
    sys, _ = load_authoring(
        struc_yaml; name="too_many_stations", set)
    sam = SymbolicAWEModel(set, sys; vsm_set)
    vsm_w = sam.aero[1].vsm_wing

    # Drop down to 2 unrefined sections while keeping
    # the 3 stations → must error.
    pop!(vsm_w.unrefined_sections)
    vsm_w.n_unrefined_sections = Int16(2)

    @test_throws ErrorException compute_spatial_station_mapping!(sam, 1)
end
nothing
