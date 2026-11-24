# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# ==================== SYSTEM STRUCTURE ==================== #

"""The bodies of the wings of `sys`, the first of its `bodies`."""
wing_bodies(sys) = sys.wings
wing_bodies(sys::KiteGeometry.SystemDefinition) = view(sys.bodies, 1:length(sys.wings))

"""
    calc_heading(sys::SystemState)

Calculate heading angles for all wings using the tangential
sphere frame method. Returns a vector of heading angles, one
per wing.
"""
function calc_heading(sys::SystemState)
    return [calc_heading(wing.R_b_to_w, wing.pos_w) for wing in sys.wings]
end

"""
    wing_frame_member(point, wing_idx) -> Bool

Whether `point` belongs to wing `wing_idx`'s frame: it places the wing body and,
without a mesh, shapes its own inertia. Its aerodynamic-surface nodes
(`is_wing_node`) plus any `BODY_STATIC` structural point sharing its `wing_idx`.
"""
wing_frame_member(point, wing_idx) =
    (point.is_wing_node || point_type(point) == BODY_STATIC) && point.wing_idx == wing_idx

"""
    carrier_body_idx(point, bodies) -> Int

Index of the body whose `total_mass` holds `point`'s mass: the DYNAMIC or STATIC
body a `BODY_STATIC` point rides, or the `RIGID_DYNAMICS` wing a wing node belongs
to. 0 when no integrated body carries it, and the point then carries its own weight.
"""
function carrier_body_idx(point, bodies)
    body = point.body
    !isnothing(body) && body_type(bodies[body]) != KINEMATIC && return body
    wing = point.wing_idx
    point.is_wing_node && !isnothing(wing) && body_type(bodies[wing]) != KINEMATIC &&
        return wing
    return 0
end

"""
    particle_wing_parts(wing, stations, points, bodies, wings, root)
        -> (point_idxs, body_idxs)

The member points of a PARTICLE_DYNAMICS `wing` and the bodies other than `wings` (the
first of `bodies`) sharing a beam component (`root`, from KiteGeometry's
`connected_body_groups`) with its structure, which carry its mass. The components are
seeded from its stations' member/flap bodies and the bodies its member points ride.
Members are the union of the wing's station points, falling back to `wing_frame_member`
points when it has no stations.
"""
function particle_wing_parts(wing, stations, points, bodies, wings, root)
    seed_roots = Set{Int64}()
    point_idxs = Set{Int64}()
    for station_idx in wing.station_idxs
        station = stations[station_idx]
        union!(point_idxs, station.points)
        for body_idx in Iterators.flatten(
                (something(station.bodies, ()), something(station.flap_bodies, ())))
            body_idx != 0 && body_idx != wing.idx && push!(seed_roots, root[body_idx])
        end
    end
    if isempty(point_idxs)
        for point in points
            wing_frame_member(point, wing.idx) && push!(point_idxs, point.idx)
        end
    end
    for idx in point_idxs
        body_idx = points[idx].body
        !isnothing(body_idx) && body_idx != wing.idx && push!(seed_roots, root[body_idx])
    end
    body_idxs = [body_idx for body_idx in eachindex(bodies)
                 if body_idx ∉ eachindex(wings) && root[body_idx] in seed_roots]
    return point_idxs, body_idxs
end

"""
    compute_spatial_station_mapping!(the_wing, stations, points)

Give each unrefined VSM section, and into `panel_idxs` each panel, to the station closest
in the body frame, and split the wing's own mass by panel ([`share_body_mass!`](@ref)). A
station's twist moment sums over `panel_idxs`, and one owning several sections drives them
as a rigid unit. Errors if there are more stations than sections or a station claims none.
"""
function compute_spatial_station_mapping!(the_wing, stations::AbstractVector,
                                          points::AbstractVector)
    the_vsm_wing = the_wing.vsm_wing
    n_unrefined = the_vsm_wing.n_unrefined_sections
    n_stations = length(the_wing.station_idxs)

    n_stations <= n_unrefined || error(
        "Wing $(the_wing.idx): n_stations " *
        "($n_stations) > n_unrefined sections " *
        "($n_unrefined). Reduce stations or increase " *
        "aero resolution.")

    # Compute station centers in body frame
    station_centers = Vector{MVec3}(undef, n_stations)
    for (local_idx, station_idx) in
            enumerate(the_wing.station_idxs)
        station = stations[station_idx]
        center = zeros(3)
        for pt_idx in station.points
            center += initial_body_position(the_wing, points[pt_idx])
        end
        station_centers[local_idx] =
            center / length(station.points)
    end

    offset_vec = [0.0, 0.0, the_wing.aero_z_offset]
    unrefined_centers = Vector{MVec3}(
        undef, n_unrefined)
    for i in 1:n_unrefined
        le_point =
            the_vsm_wing.unrefined_sections[i].LE_point
        te_point =
            the_vsm_wing.unrefined_sections[i].TE_point
        unrefined_centers[i] =
            (le_point + te_point) / 2 .- offset_vec
    end

    # Reset section and panel lists (we rebuild the partition)
    for station_idx in the_wing.station_idxs
        empty!(stations[station_idx].unrefined_section_idxs)
        empty!(stations[station_idx].panel_idxs)
    end

    for section_idx in 1:n_unrefined
        closest_local = nearest_station(unrefined_centers[section_idx], station_centers)
        push!(stations[the_wing.station_idxs[closest_local]].unrefined_section_idxs,
              Int64(section_idx))
    end
    for (panel_idx, panel) in enumerate(the_wing.vsm_aero.panels)
        panel_center = vec(sum(panel.corner_points; dims=2)) / 4 .- offset_vec
        closest_local = nearest_station(panel_center, station_centers)
        push!(stations[the_wing.station_idxs[closest_local]].panel_idxs,
              Int64(panel_idx))
    end

    share_body_mass!(the_wing, stations)

    # Every station must claim at least one section
    for station_idx in the_wing.station_idxs
        station = stations[station_idx]
        isempty(station.unrefined_section_idxs) && error(
            "Wing $(the_wing.idx): station " *
            "$(station.name) claims no unrefined " *
            "sections (likely coincident station centres).")
    end
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
    share_body_mass!(the_wing, stations)

Set each station's `body_mass`: the wing's `extra_mass`, which sits on no point, times
the station's share of the wing's panel area.
"""
function share_body_mass!(the_wing, stations::AbstractVector)
    panels = the_wing.vsm_aero.panels
    wing_area = sum(panel.chord * panel.width for panel in panels)
    for station_idx in the_wing.station_idxs
        station = stations[station_idx]
        station_area = sum(panels[i].chord * panels[i].width for i in station.panel_idxs;
                           init=0.0)
        station.body_mass = the_wing.extra_mass * station_area / wing_area
    end
    return nothing
end

# ==================== CONSTRUCTOR ==================== #

"""
    update_mass_properties!(sys_struct; prn=true)

Set every DYNAMIC and STATIC body's `total_mass`, `com_offset_b` and principal
inertia (`inertia_principal`, `R_b_to_p`): its own mass properties plus
each point it carries ([`carrier_body_idx`](@ref)) as a point mass of the point's
`total_mass`, at its anchor. A particle wing's `total_mass` sums its free member points
and its section bodies. Point masses are taken at the current `rest_length`.
"""
function update_mass_properties!(sys_struct; prn=true)
    (; points, bodies) = sys_struct
    stations = sys_struct.spanwise_stations
    for body in bodies
        body_type(body) == KINEMATIC && continue
        old_com_offset = copy(body.com_offset_b)
        combine_carried_points!(body, points, bodies)
        # A rigid wing's station leading edges are measured from its COM.
        for station_idx in body.station_idxs
            stations[station_idx].le_pos .+= old_com_offset .- body.com_offset_b
        end
        prn && body.idx in eachindex(sys_struct.wings) &&
            @info "RIGID_DYNAMICS wing $(body.idx):" *
            " total_mass=$(round(body.total_mass; digits=3))" *
            ", I=$(round.(body.inertia_principal; digits=4))" *
            ", com_offset_b=$(round.(body.com_offset_b; digits=4))"
    end
    root = KiteGeometry.connected_body_groups(length(bodies), structure_tubes(sys_struct))
    for wing in sys_struct.wings
        body_dynamics_type(wing) == PARTICLE_DYNAMICS || continue
        point_idxs, body_idxs = particle_wing_parts(wing, stations, points, bodies,
                                                    sys_struct.wings, root)
        wing.total_mass = sum(points[idx].total_mass for idx in point_idxs
                              if carrier_body_idx(points[idx], bodies) == 0; init=0.0) +
            sum(bodies[idx].total_mass for idx in body_idxs; init=0.0)
    end
    return nothing
end

"""
    combine_carried_points!(body, points, bodies)

Add the points `body` carries to its own mass properties: `total_mass`, the COM
`com_offset_b`, and the inertia about it diagonalised by `principal_frame_method`.
"""
function combine_carried_points!(body, points, bodies)
    carried = [(point.total_mass, carried_position_b(point, body)) for point in points
               if carrier_body_idx(point, bodies) == body.idx]
    total_mass = body.extra_mass + sum(first, carried; init=0.0)
    first_moment = body.extra_mass .* body.extra_com_offset_b
    for (mass, position) in carried
        first_moment .+= mass .* position
    end
    com = total_mass > 0 ? first_moment ./ total_mass : body.extra_com_offset_b
    inertia = something(body.extra_inertia_KA, zeros(3, 3)) .+
        point_mass_inertia(body.extra_mass, body.extra_com_offset_b .- com)
    for (mass, position) in carried
        inertia .+= point_mass_inertia(mass, position .- com)
    end
    inertia_principal, R_b_to_p = principal_axes(inertia, body.principal_frame_method)
    body.total_mass = total_mass
    body.com_offset_b .= com
    body.inertia_principal .= inertia_principal
    body.R_b_to_p .= R_b_to_p
    return nothing
end

"""
    carried_position_b(point, body) -> Vector

Body-frame position of a point `body` carries, from its body origin: a rider's
`anchor_b`, or where its initial position sits on the body in its initial pose.
"""
carried_position_b(point, body) = point.body == body.idx ? Vector(point.anchor_b) :
    initial_body_position(body, point)

"""
    system_definition(name, set; points, segments, materials, spanwise_stations, pulleys,
                      tethers, winches, wings, bodies, tubes, ignore_l0=false)

The `SystemDefinition` named `name` of KiteGeometry rows, completed from `set` and placed:
a winch without a `gear_ratio` or `drum_radius` takes `set`'s, every wing's body is put in
its design pose and given `set.mass` where it carries none (KiteGeometry's
`complete_wings!`), and the whole is placed by KiteGeometry's `place!`, with every
segment's `rest_length` set to its placed length where `ignore_l0`. `materials` holds the
`Material` rows the segments and tubes name, as a vector or as the `MaterialTable`
KiteGeometry's `line_segment` and `expand_tether!` fill. A wing's body is one of `bodies`,
which the wings' bodies lead.

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
        bodies=wings_first(wings, bodies), tubes)
    complete_wings!(definition, set)
    return place!(definition; ignore_l0)
end

"""`bodies` with the body of each of `wings` first, in the wings' order."""
function wings_first(wings, bodies)
    leading = map(wings) do wing
        index = findfirst(body -> body.name == wing.body, bodies)
        isnothing(index) && error("Wing $(wing.name): its body $(wing.body) is none of " *
                                  "`bodies`")
        index
    end
    return [bodies[leading]; bodies[setdiff(eachindex(bodies), leading)]]
end

"""The rows of the tubes of `sys`, each naming its two bodies by index."""
structure_tubes(sys) = sys.tubes
structure_tubes(state::SystemState) = getfield(state, :structure).tubes
