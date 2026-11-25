# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# ==================== STATION MAPPING ==================== #

"""
    compute_spatial_station_mapping!(sam, w)

Give each unrefined VSM section of wing `w` of `sam`, and each of its panels, to the
station closest in the body frame, into its engine's `station_sections` and
`station_panels`, and split the wing's own mass by panel ([`share_body_mass!`](@ref)).
A station's twist moment sums over its panels, and one owning several sections drives
them as a rigid unit. Errors if there are more stations than sections or a station
claims none.
"""
function compute_spatial_station_mapping!(sam, w)
    engine = sam.aero[w].engine
    station_idxs = wing_station_idxs(sam, w)
    the_vsm_wing = engine.vsm_wing
    n_unrefined = the_vsm_wing.n_unrefined_sections
    n_stations = length(station_idxs)

    n_stations <= n_unrefined || error(
        "Wing $w: n_stations ($n_stations) > n_unrefined sections ($n_unrefined). " *
        "Reduce stations or increase aero resolution.")

    station_centers = [station_center_b(sam, w, s) for s in station_idxs]
    offset_vec = [0.0, 0.0, engine.aero_z_offset]
    unrefined_centers = [(section.LE_point + section.TE_point) / 2 .- offset_vec
                         for section in the_vsm_wing.unrefined_sections[1:n_unrefined]]

    engine.station_sections = [Int[] for _ in station_idxs]
    engine.station_panels = [Int[] for _ in station_idxs]
    for section_idx in 1:n_unrefined
        closest = nearest_station(unrefined_centers[section_idx], station_centers)
        push!(engine.station_sections[closest], section_idx)
    end
    for (panel_idx, panel) in enumerate(engine.vsm_aero.panels)
        panel_center = vec(sum(panel.corner_points; dims=2)) / 4 .- offset_vec
        push!(engine.station_panels[nearest_station(panel_center, station_centers)],
              panel_idx)
    end

    share_body_mass!(sam, w)

    for (k, s) in enumerate(station_idxs)
        isempty(engine.station_sections[k]) && error(
            "Wing $w: station $(sam.definition.spanwise_stations[s].name) claims no " *
            "unrefined sections (likely coincident station centres).")
    end
    return nothing
end

"""The centre of the points of station `s` of `sam` in the frame of the body of wing
`w`, from its origin, in the initial pose [m]."""
function station_center_b(sam, w, s)
    b = wing_body(sam, w)
    points = sam.definition.spanwise_stations[s].points
    positions = initial_positions(sam)
    return sum(initial_body_position(sam, b, positions[i]) for i in points) /
        length(points)
end

"""
    nearest_station(position, station_centers) -> Int

Index of the body-frame station centre nearest `position`. Centres within `1e-9` of
the nearest distance, relative, tie, and a tie goes to the outboard one (largest
`|y|`), so a mirror-symmetric wing splits symmetrically.
"""
function nearest_station(position, station_centers)
    dists = [norm(position - center) for center in station_centers]
    tied = findall(<=(minimum(dists) * (1 + 1e-9)), dists)
    return argmax(local_idx -> abs(station_centers[local_idx][2]), tied)
end

"""
    share_body_mass!(sam, w)

Set each station's `station_body_mass` in the engine of wing `w` of `sam`: the
`extra_mass` of the wing's body, which sits on no point, times the station's share of
the wing's panel area.
"""
function share_body_mass!(sam, w)
    engine = sam.aero[w].engine
    panels = engine.vsm_aero.panels
    wing_area = sum(panel.chord * panel.width for panel in panels)
    extra_mass = something(sam.definition.bodies[wing_body(sam, w)].extra_mass, 0.0)
    engine.station_body_mass = [
        extra_mass * sum(panels[i].chord * panels[i].width for i in station_panels;
                         init=0.0) / wing_area
        for station_panels in engine.station_panels]
    return nothing
end

# ==================== CONSTRUCTOR ==================== #

"""
    system_definition(name, set; points, segments, materials, spanwise_stations, pulleys,
                      tethers, winches, wings, bodies, tubes, ignore_l0=false)

The `SystemDefinition` named `name` of KiteGeometry rows, completed from `set` and placed:
a winch without a `gear_ratio` or `drum_radius` takes `set`'s, every wing's body is put in
its design pose and given `set.mass` where it carries none (KiteGeometry's
`complete_wings!`), and the whole is placed by KiteGeometry's `place!`, with every
segment's `rest_length` set to its placed length where `ignore_l0`. `materials` holds the
`Material` rows the segments and tubes name, as a vector or as the `MaterialTable`
KiteGeometry's `line_segment` and `expand_tether!` fill. A wing's body is one of `bodies`.

`name` is also the `system_name` that keys the compiled-model cache, so systems of the
same structural shape should share one.
"""
function system_definition(name::AbstractString, set;
        points=KiteGeometry.Point[], segments=KiteGeometry.Segment[],
        materials=KiteGeometry.Material[],
        spanwise_stations=KiteGeometry.SpanwiseStation[], pulleys=KiteGeometry.Pulley[],
        tethers=KiteGeometry.Tether[], winches=KiteGeometry.Winch[],
        wings=KiteGeometry.Wing[], bodies=KiteGeometry.Body[], tubes=KiteGeometry.Tube[],
        ignore_l0::Bool=false)
    for winch in winches
        winch.gear_ratio = something(winch.gear_ratio, set.gear_ratio)
        winch.drum_radius = something(winch.drum_radius, set.drum_radius)
    end
    materials = materials isa KiteGeometry.MaterialTable ? materials.rows : materials
    definition = SystemDefinition(;
        metadata=Metadata(; name, n_points=length(points)), materials, points, segments,
        spanwise_stations, pulleys, tethers, winches, wings,
        bodies, tubes)
    complete_wings!(definition, set)
    return place!(definition; ignore_l0)
end
