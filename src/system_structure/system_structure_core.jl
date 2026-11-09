# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# ==================== SYSTEM STRUCTURE ==================== #

"""
    SystemStructure

The old name of KiteGeometry's `SystemDefinition`: a kite system's points, segments,
tethers, winches, wings, bodies and tubes as rows of the same type a structure document
holds. Its `wings` are KiteGeometry's rows, each a wing's frame; a wing's body is the first
of `bodies` under the same name. `sys.stations` lists the `spanwise_stations` and the
stations a body that is no wing owns, kept on that body.
"""
const SystemStructure = KiteGeometry.SystemDefinition

"""The bodies of the wings of `sys`, the first of its `bodies`."""
wing_bodies(sys) = sys.wings
wing_bodies(sys::KiteGeometry.SystemDefinition) = view(sys.bodies, 1:length(sys.wings))

"""The stations of `sys`: its `spanwise_stations`, then the ones no wing owns."""
function structure_stations(sys)
    stations = collect(Station, getfield(sys, :spanwise_stations))
    for body in getfield(sys, :bodies)
        append!(stations, body_model(body).stations_without_wing)
    end
    return NamedCollection(stations)
end

function Base.getproperty(sys::KiteGeometry.SystemDefinition, name::Symbol)
    name === :stations && return structure_stations(sys)
    name === :name && return getfield(getfield(sys, :metadata), :name)
    return getfield(sys, name)
end


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
    build_name_dict(items::Vector) -> Dict{Symbol, Int64}

Build a name→index dictionary from a vector of items with optional `name` fields.
Items with `name=nothing` are skipped. Integer names are converted to Symbols.
"""
function build_name_dict(items::Vector)
    names = Dict{Symbol, Int64}()
    for (i, item) in enumerate(items)
        isnothing(item.name) && continue
        name = Symbol(item.name)
        haskey(names, name) &&
            error("Duplicate name '$name' found at indices $(names[name]) and $i")
        names[name] = i
    end
    return names
end


"""
    resolve_ref(ref::NameRef, name_dict::Dict{Symbol, Int64}, component_type::String) -> Int64

Resolve a reference (name or index) to an index using the name dictionary.
If ref is an integer, returns it directly. If ref is a symbol, looks up in dictionary.
"""
resolve_ref(index::Int, ::Dict{Symbol, Int64}, ::String) = Int64(index)
resolve_ref(::Nothing, ::Dict{Symbol, Int64}, ::String) = Int64(0)
function resolve_ref(name::Symbol, name_dict::Dict{Symbol, Int64}, component_type::String)
    haskey(name_dict, name) || error("Unknown $component_type name: $name")
    return name_dict[name]
end

"""
    resolve!(ref_pt::WeightedRefPoints, name_dict, type)

Resolve symbolic refs to integer indices, filling
`ref_pt.ids`. No-op if `refs` is empty (already resolved).
"""
function resolve!(
    ref_pt::WeightedRefPoints,
    name_dict::Dict{Symbol, Int64},
    component_type::String
)
    isempty(ref_pt.refs) && return
    ref_pt.ids = Int64[
        resolve_ref(ref, name_dict, component_type)
        for ref in ref_pt.refs]
end

"""
    expand_tethers!(points, segments, tethers, materials, set)

Expand each segment-count tether through KiteGeometry's `expand_tether!`, adding its inner
points and its segments, which carry its material and spring. Must be called before
`assign_indices_and_resolve!`.
"""
function expand_tethers!(points::Vector{Point}, segments::AbstractVector,
                         tethers::Vector{Tether}, materials, set::Settings)
    point_transforms = fill(nothing, length(points))
    table = KiteGeometry.MaterialTable(materials, Dict{String, Any}())
    seg_names = Set{Symbol}(Symbol(segment.name) for segment in segments)
    for (i, tether) in enumerate(tethers)
        isnothing(tether.start_point_ref) && continue
        first_segment = tether.segment_refs[1]
        (first_segment isa Int || first_segment in seg_names) && continue
        columns = tether.segment_columns
        spring = SegmentSpring(;
            unit_damping=nan_to_nothing(get(columns, :unit_damping, NaN)),
            damping_per_stiffness=nan_to_nothing(get(columns, :damping_per_stiffness, NaN)),
            compression_frac=nan_to_nothing(get(columns, :compression_frac, NaN)),
            compression_damping_frac=nan_to_nothing(
                get(columns, :compression_damping_frac, NaN)))
        line = Base.structdiff(columns, NamedTuple{fieldnames(SegmentSpring)})
        expanded = KiteGeometry.expand_tether!(points, segments, point_transforms,
            tether.name, tether.start_point_ref, tether.end_point_ref, tether.n_segments,
            (; stretched_length=tether.init_stretched_len), set, table;
            model=tether.model, segment_model=spring, line...)
        tethers[i] = expanded
        expanded.model.segment_columns = (;)
    end
    return nothing
end

"""
    mark_wing_nodes!(points, stations)

Set `point.is_wing_node` for every point that is a member of some station.
Requires `station.point_idxs` to be resolved first.
"""
function mark_wing_nodes!(points, stations)
    members = Set{Int64}()
    for station in stations, point_idx in station.point_idxs
        push!(members, point_idx)
    end
    for point in points
        point.is_wing_node = point.idx in members
    end
    return nothing
end

"""
    wing_frame_member(point, wing_idx) -> Bool

Whether `point` belongs to wing `wing_idx`'s frame: it places the wing body and,
without a mesh, shapes its own inertia. Its aerodynamic-surface nodes
(`is_wing_node`) plus any `BODY_STATIC` structural point sharing its `wing_idx`.
"""
wing_frame_member(point, wing_idx) =
    (point.is_wing_node || point.type == BODY_STATIC) && point.wing_idx == wing_idx

"""
    carrier_body_idx(point, bodies) -> Int

Index of the body whose `total_mass` holds `point`'s mass: the DYNAMIC or STATIC
body a `BODY_STATIC` point rides, or the `RIGID_DYNAMICS` wing a wing node belongs
to. 0 when no integrated body carries it, and the point then carries its own weight.
"""
function carrier_body_idx(point, bodies)
    point.type == BODY_STATIC && point.body_idx > 0 &&
        bodies[point.body_idx].type != KINEMATIC && return point.body_idx
    point.is_wing_node && point.wing_idx > 0 &&
        bodies[point.wing_idx].type != KINEMATIC && return point.wing_idx
    return 0
end

"""
    distribute_mass_over_points!(points, point_idxs, mass)

Split `mass` equally across `point_idxs`, writing each point's `extra_mass`. Used
for a wing given a lumped `set.mass` rather than per-point masses.
"""
function distribute_mass_over_points!(points, point_idxs, mass)
    isempty(point_idxs) && return nothing
    for idx in point_idxs
        points[idx].extra_mass = mass / length(point_idxs)
    end
    return nothing
end

"""
    particle_wing_parts(wing, stations, points, bodies, root)
        -> (point_idxs, body_idxs)

The member points of a PARTICLE_DYNAMICS `wing` and the non-wing bodies sharing a
beam component (`root`, from KiteGeometry's `connected_body_groups`) with its structure,
which carry its mass. The components are seeded from its stations' member/flap
bodies and the bodies its member points ride. Members are the union of the wing's
station points, falling back to `wing_frame_member` points when it has no stations.
"""
function particle_wing_parts(wing, stations, points, bodies, root)
    seed_roots = Set{Int64}()
    point_idxs = Set{Int64}()
    for station_idx in wing.station_idxs
        station = stations[station_idx]
        union!(point_idxs, station.point_idxs)
        for body_idx in Iterators.flatten(
                (station.body_idxs, station.flap_body_idxs))
            body_idx != 0 && body_idx != wing.idx && push!(seed_roots, root[body_idx])
        end
    end
    if isempty(point_idxs)
        for point in points
            wing_frame_member(point, wing.idx) && push!(point_idxs, point.idx)
        end
    end
    for idx in point_idxs
        body_idx = points[idx].body_idx
        body_idx != 0 && body_idx != wing.idx && push!(seed_roots, root[body_idx])
    end
    body_idxs = [body_idx for (body_idx, body) in enumerate(bodies)
                 if !is_wing(body) && root[body_idx] in seed_roots]
    return point_idxs, body_idxs
end

"""
    finalize_particle_wing_mass!(wing, stations, points, bodies, set, root)

Distribute `set.mass` over a PARTICLE_DYNAMICS `wing`'s member points when neither
they nor its section bodies (see [`particle_wing_parts`](@ref)) carry `extra_mass`,
warning when that leaves the wing massless. Deferred until point/station→body and
tube→body references resolve.
"""
function finalize_particle_wing_mass!(wing, stations, points, bodies, set, root)
    point_idxs, body_idxs = particle_wing_parts(wing, stations, points, bodies, root)
    own_mass = sum(points[idx].extra_mass for idx in point_idxs; init=0.0) +
        sum(bodies[idx].extra_mass for idx in body_idxs; init=0.0)
    own_mass > 0 && return nothing
    set_mass = hasproperty(set, :mass) ? set.mass : 0.0
    if set_mass > 0
        wing_point_idxs = [point.idx for point in points
            if wing_frame_member(point, wing.idx)]
        distribute_mass_over_points!(points, wing_point_idxs, set_mass)
    else
        @warn "Wing $(wing.idx) (PARTICLE_DYNAMICS) has zero mass — no member " *
            "point extra_mass and no connected body mass."
    end
    return nothing
end

"""
    assign_indices_and_resolve!(components, name_dicts)

Assign indices to all components based on their position in the vectors,
and resolve all references to indices.
"""
function assign_indices_and_resolve!(
    points::Vector{Point},
    stations::Vector{Station},
    segments::AbstractVector{<:Segment},
    pulleys::Vector{Pulley},
    tethers::Vector{Tether},
    winches::AbstractVector{<:Winch},
    wings::AbstractVector{<:Body}
)
    # First pass: assign indices based on position.
    for components in (points, stations, pulleys, tethers, winches, wings),
            (i, component) in enumerate(components)
        component.idx = i
    end

    # Build name lookup dictionaries
    point_names = build_name_dict(points)
    station_names = build_name_dict(stations)
    segment_names = build_name_dict(segments)
    pulley_names = build_name_dict(pulleys)
    tether_names = build_name_dict(tethers)
    winch_names = build_name_dict(winches)
    wing_names = build_name_dict(wings)

    # Points: resolve wing refs (body_ref resolved later, no bodies here).
    for point in points
        point.wing_idx = resolve_ref(point.wing_ref, wing_names, "wing")
        point.wing_idx in 0:length(wings) || error(
            "Point $(point.name): wing $(point.wing_ref) does not exist; the " *
            "system has $(length(wings)) wings.")
    end

    # Stations: resolve point_refs
    for station in stations
        station.point_idxs = Int64[resolve_ref(ref, point_names, "point") for ref in station.point_refs]
    end
    mark_wing_nodes!(points, stations)

    # Segments: resolve point_refs
    for segment in segments
        point1 = resolve_ref(segment.point_refs[1], point_names, "point")
        point2 = resolve_ref(segment.point_refs[2], point_names, "point")
        segment.point_idxs = (point1, point2)
    end

    # Pulleys: resolve segment_refs
    for pulley in pulleys
        segment1 = resolve_ref(pulley.segment_refs[1], segment_names, "segment")
        segment2 = resolve_ref(pulley.segment_refs[2], segment_names, "segment")
        pulley.segment_idxs = (segment1, segment2)
    end

    # Tethers: resolve segment_refs and start/end point refs
    for tether in tethers
        tether.segment_idxs = Int64[
            resolve_ref(ref, segment_names, "segment")
            for ref in tether.segment_refs]
        if !isnothing(tether.start_point_ref)
            tether.start_point_idx = resolve_ref(
                tether.start_point_ref, point_names, "point")
        end
        if !isnothing(tether.end_point_ref)
            tether.end_point_idx = resolve_ref(
                tether.end_point_ref, point_names, "point")
        end
    end

    # Winches: resolve tether_refs and winch_point_ref
    for winch in winches
        winch.tether_idxs = Int64[
            resolve_ref(ref, tether_names, "tether")
            for ref in winch.tether_refs]
        winch.winch_point_idx = resolve_ref(
            winch.winch_point_ref, point_names, "point")
    end

    # Wings: resolve station_refs and PARTICLE_DYNAMICS-specific refs
    for wing in wings
        # BaseWing fields
        wing.station_idxs = Int64[resolve_ref(ref, station_names, "station") for ref in wing.station_refs]

        # Body-frame reference points (any wing may carry them)
        if !isnothing(wing.origin)
            resolve!(wing.origin, point_names, "point")
        end
        if !isnothing(wing.z_ref_points)
            resolve!(something(wing.z_ref_points)[1],
                point_names, "point")
            resolve!(something(wing.z_ref_points)[2],
                point_names, "point")
        end
        if !isnothing(wing.y_ref_points)
            resolve!(something(wing.y_ref_points)[1],
                point_names, "point")
            resolve!(something(wing.y_ref_points)[2],
                point_names, "point")
        end

        # Resize per-mode aero state now that station_idxs are resolved.
        resize_aero_state!(wing.aero, wing)
    end

    return (point_names, station_names, segment_names, pulley_names,
            tether_names, winch_names, wing_names)
end

"""
    init_body_frame_from_ref_points!(wing, points; prn=true)

Initialize wing body frame (`Q_KA_to_ENU`, `pos_ENU`) from z/y
reference points. Shared by VSMWing PARTICLE_DYNAMICS and PlateWing.
"""
function init_body_frame_from_ref_points!(
    wing, points; prn=true
)
    isnothing(wing.origin) && return
    isnothing(wing.z_ref_points) && return
    isnothing(wing.y_ref_points) && return

    origin_pos = get_ref_position_from_points(
        points, wing.origin; field=:pos_ENU)
    wing.pos_ENU .= origin_pos

    R_b_to_w, _ = calc_particle_dynamics_wing_frame(
        points, wing.z_ref_points,
        wing.y_ref_points, wing.origin; field=:pos_ENU)
    wing.Q_KA_to_ENU .= rotation_matrix_to_quaternion(R_b_to_w)
    wing.R_b_to_p .= Matrix{SimFloat}(I, 3, 3)

    if prn
        origin_rounded = round.(origin_pos; digits=3)
        @info "Wing $(wing.name) ($(typeof(wing).name.name))" *
              ": origin=[$origin_rounded]"
    end
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
        for pt_idx in station.point_idxs
            center += initial_body_position(the_wing, points[pt_idx])
        end
        station_centers[local_idx] =
            center / length(station.point_idxs)
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
    setup_wing_frame!(wing, points; prn=true)

Compute a wing's body frame (`Q_KA_to_ENU`, `pos_ENU`) and, for `RIGID_DYNAMICS`, its own
COM offset and inertia (`extra_com_offset_b`, `extra_inertia_b`) from its mesh or points and
ref points; without ref points the frame keeps the world axes, origin at the wing's COM.
Without a mesh the inertia comes from the frame points' `extra_mass`; with neither a mesh
nor point masses, the constructor's inertia stays, or none for a massless wing.
"""
function setup_wing_frame!(wing, points; prn=true)
    if wing.dynamics_type == RIGID_DYNAMICS
        any(wing_frame_member(point, wing.idx)
            for point in points) || return nothing

        com, inertia_normalized = normalized_inertia(wing.aero, wing, points)

        # Body frame from ref points (else body = world orientation, origin = its own COM)
        origin = wing.origin
        z_ref = wing.z_ref_points
        y_ref = wing.y_ref_points
        if !isnothing(origin) && !isnothing(z_ref) && !isnothing(y_ref)
            origin_pos = get_ref_position_from_points(
                points, origin; field=:pos_ENU)
            wing.pos_ENU .= origin_pos
            R_b_to_w, _ = calc_particle_dynamics_wing_frame(
                points, z_ref, y_ref, origin; field=:pos_ENU)
            wing.Q_KA_to_ENU .= rotation_matrix_to_quaternion(R_b_to_w)
            wing.extra_com_offset_b .= R_b_to_w' * (com - origin_pos)
        else
            wing.pos_ENU .= com
            wing.Q_KA_to_ENU .= SimFloat[1, 0, 0, 0]
            wing.extra_com_offset_b .= 0.0
        end
        wing.com_offset_b .= wing.extra_com_offset_b

        # The hook returns per-unit-mass inertia [m²] in the world axes of `pos_ENU`.
        if !isnothing(inertia_normalized)
            R_initial = initial_rotation(wing)
            wing.extra_inertia_b .=
                R_initial' * (wing.extra_mass .* inertia_normalized) * R_initial
        elseif iszero(wing.extra_mass)
            wing.extra_inertia_b .= 0.0
        end
    else  # PARTICLE_DYNAMICS (VSM or flat-plate)
        init_body_frame_from_ref_points!(wing, points; prn)
    end
    return nothing
end

"""
    update_mass_properties!(sys_struct; prn=true)

Set every DYNAMIC and STATIC body's `total_mass`, `com_offset_b` and principal
inertia (`inertia_principal`, `R_b_to_p`): its own mass properties plus
each point it carries ([`carrier_body_idx`](@ref)) as a point mass of the point's
`total_mass`, at its anchor. A particle wing's `total_mass` sums its free member points
and its section bodies. Point masses are taken at the current `l0`, and the rigid wing
nodes' offsets from the new COM re-derived ([`init_undeformed_offsets!`](@ref)).
"""
function update_mass_properties!(sys_struct; prn=true)
    (; points, bodies, stations) = sys_struct
    write_total_mass!(sys_struct)
    for body in bodies
        body.type == KINEMATIC && continue
        old_com_offset = copy(body.com_offset_b)
        combine_carried_points!(body, points, bodies)
        # A rigid wing's station leading edges are measured from its COM.
        for station_idx in body.station_idxs
            stations[station_idx].le_pos .+= old_com_offset .- body.com_offset_b
        end
        prn && is_wing(body) && @info "RIGID_DYNAMICS wing $(body.idx):" *
            " total_mass=$(round(body.total_mass; digits=3))" *
            ", I=$(round.(body.inertia_principal; digits=4))" *
            ", com_offset_b=$(round.(body.com_offset_b; digits=4))"
    end
    root = KiteGeometry.connected_body_groups(length(bodies),
                                              placement_tubes(sys_struct.tubes))
    for wing in sys_struct.wings
        wing.dynamics_type == PARTICLE_DYNAMICS || continue
        point_idxs, body_idxs = particle_wing_parts(wing, stations, points, bodies, root)
        wing.total_mass = sum(points[idx].total_mass for idx in point_idxs
                              if carrier_body_idx(points[idx], bodies) == 0; init=0.0) +
            sum(bodies[idx].total_mass for idx in body_idxs; init=0.0)
    end
    init_undeformed_offsets!(bodies, points)
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
    inertia = body.extra_inertia_b .+
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
carried_position_b(point, body) = point.body_idx == body.idx ? Vector(point.anchor_b) :
    initial_body_position(body, point)

"""
    SystemStructure(name, set; points, stations, segments, pulleys, tethers,
                    winches, wings, bodies, tubes)

Constructs a `SystemStructure` object representing a complete mechanical system.
Resolves every symbolic reference to a numeric index, fills in derived properties
and validates the result.

`name` is also the `system_name` that keys the compiled-model cache, so systems
of the same structural shape should share one.

# Arguments
- `name::String`: Model identifier, free-form.
- `set::Settings`: Configuration parameters from `KiteUtils.jl`.

# Keyword Arguments
- `points`, `stations`, `segments`, etc.: Vectors of the system components.
- `vsm_set`: `VSMSettings` for VSM wings; read from the project file's
  `vsm_settings:` when omitted.
- `ignore_l0::Bool=false`: Set every segment `l0` to its placed length
  ([`relax_segments!`](@ref)).
- `placed::Bool=false`: The positions given are the initial pose already: keep them
  and every `l0` instead of running [`place!`](@ref).
- `prn::Bool=true`: If true, print info messages about auto-generated components.
- `definition`: The KiteGeometry definition the components were read from; written from
  them ([`system_definition`](@ref)) when omitted.

# Returns
- `SystemStructure`: A complete system ready for building a `SymbolicAWEModel`.
"""
function SystemStructure(name::AbstractString, set;
        points=Point[],
        stations=Station[],
        segments=Segment[],
        pulleys=Pulley[],
        tethers=Tether[],
        winches=Winch[],
        wings=Body[],
        bodies=Body[],
        tubes=Tube[],
        ignore_l0::Bool=false,
        placed::Bool=false,
        vsm_set=nothing,
        prn::Bool=true,
        definition::Union{Nothing, SystemDefinition}=nothing,
    )
    # Load VSMSettings if not provided and VSM wings exist
    has_vsm_wings = any(has_vsm_engine(wing.aero) for wing in wings)
    if isnothing(vsm_set) && has_vsm_wings
        vsm_set_path = optional_project_file("vsm_settings")
        if !isnothing(vsm_set_path)
            vsm_set = VortexStepMethod.VSMSettings(
                vsm_set_path; data_prefix=false)
        end
    end

    # Narrow to a concrete element type only when wings are homogeneous; mixed dynamics keep the broader type.
    if length(wings) > 0
        W = typeof(wings[1])
        if all(wing -> typeof(wing) === W, wings) && eltype(wings) !== W
            wings = convert(Vector{W}, wings)
        end
    end

    # Expand Route 2 tethers (auto-generate points/segments)
    materials = isnothing(definition) ? KiteGeometry.Material[] : definition.materials
    segments = convert(Vector{Segment}, segments)
    expand_tethers!(points, segments, tethers, materials, set)
    segments = SegmentRow[structure_segment!(materials, segment, set)
                          for segment in segments]

    # Assign indices and resolve all symbolic references to numeric indices.
    (point_names_dict, station_names_dict,
     segment_names_dict, pulley_names_dict,
     tether_names_dict, winch_names_dict, _) =
        assign_indices_and_resolve!(
            points, stations, segments, pulleys, tethers, winches, wings)

    # If no wings defined, convert wing nodes to STATIC
    if isempty(wings)
        wing_point_idxs = findall(point -> point.is_wing_node, points)
        if !isempty(wing_point_idxs)
            @warn "No wings provided but " *
                  "$(length(wing_point_idxs)) wing-node " *
                  "points found. Converting to STATIC."
            for point in points[wing_point_idxs]
                point.type, point.wing_idx, point.body_idx, point.tube_idx = STATIC, 0, 0, 0
                point.anchor_b, point.apparent_mass, point.is_wing_node = nothing, 0.0, false
            end
        end
    end

    # A surviving wing node is a wing's own structural node, so it must name one.
    for point in points
        point.is_wing_node && point.wing_idx == 0 && error(
            "Point $(point.name) is a station member, so it is a structural node " *
            "of a wing, but it names no wing. Give it a `wing`.")
    end

    for segment in segments
        (segment.l0 ≈ 0) && (segment.l0 = segment_initial_length(segment, points))
    end
    for tether in tethers
        # Route 1: auto-detect start/end from segment chain
        if isnothing(tether.start_point_ref) &&
           !isempty(tether.segment_idxs)
            seg_first = segments[tether.segment_idxs[1]]
            seg_last = segments[tether.segment_idxs[end]]
            tether.start_point_idx =
                seg_first.point_idxs[1]
            tether.end_point_idx =
                seg_last.point_idxs[2]
        end
    end
    for wing in wings
        wing_point_idxs = [point.idx for point in points
            if wing_frame_member(point, wing.idx)]
        set_mass = hasproperty(set, :mass) ? set.mass : 0.0
        if wing.dynamics_type == PARTICLE_DYNAMICS
            wing.extra_mass > 0 && @warn "Wing $(wing.idx) (PARTICLE_DYNAMICS): " *
                "`extra_mass=$(wing.extra_mass)` is ignored — particle wings carry " *
                "mass on their points and section bodies, not on the wing body."
            # set.mass is spread by finalize_particle_wing_mass! once body refs resolve.
        elseif iszero(wing.extra_mass) && set_mass > 0 &&
               all(iszero(points[idx].extra_mass) for idx in wing_point_idxs)
            distribute_mass_over_points!(points, wing_point_idxs, set_mass)
        end
    end

    # Body frame + COM/principal inertia (dynamics, independent of aero mode).
    for wing in wings
        setup_wing_frame!(wing, points; prn)
    end

    # Per-mode aero construction (dispatched; no-op for modes without an engine).
    for wing in wings
        setup_aero!(wing.aero, wing, points, stations; prn, vsm_set)
    end

    # Wings are bodies, placed first so idx stays 1..n_wings (wing arrays index by idx).
    prepend!(bodies, wings)
    wing_names_dict = build_name_dict(wings)
    for (i, body) in enumerate(bodies)
        body.idx = i
        body.wing_idx = resolve_ref(body.wing_ref, wing_names_dict, "wing")
    end
    rigid_body_names_dict = build_name_dict(bodies)

    # Body-anchored points: resolve the body, deriving anchor_b from pos_ENU when unset.
    for point in points
        point.body_idx = resolve_ref(
            point.body_ref, rigid_body_names_dict, "rigid_body")
        # A BODY_STATIC point given only a `wing` rides that wing's own body.
        if point.type == BODY_STATIC && point.body_idx == 0 && point.wing_idx > 0
            point.body_idx = point.wing_idx
        elseif point.body_idx > 0 && point.wing_idx > 0 &&
               is_wing(bodies[point.body_idx]) && point.body_idx != point.wing_idx
            error("Point $(point.name): `body` and `wing` name different wings " *
                "(body $(point.body_idx) vs wing $(point.wing_idx)); a wing is a " *
                "body, so they must reference the same one.")
        end
        if point.body_idx > 0 && iszero(point.anchor_b)
            body = bodies[point.body_idx]
            point.resolved_anchor_b = initial_body_position(body, point)
        end
    end

    # Tubes: assign indices, resolve their body references.
    for (i, tube) in enumerate(tubes)
        tube.idx = i
        tube.body_a_idx = resolve_ref(
            tube.body_a_ref, rigid_body_names_dict, "rigid_body")
        tube.body_b_idx = resolve_ref(
            tube.body_b_ref, rigid_body_names_dict, "rigid_body")
    end
    tube_names_dict = build_name_dict(tubes)

    # Stations: resolve owning-wing and (member + flap) body references.
    for station in stations
        station.wing_idx = resolve_ref(
            station.wing_ref, rigid_body_names_dict, "rigid_body")
        station.body_idxs = Int64[resolve_ref(ref, rigid_body_names_dict,
            "rigid_body") for ref in station.body_refs]
        station.flap_body_idxs = Int64[resolve_ref(ref, rigid_body_names_dict,
            "rigid_body") for ref in station.flap_body_refs]
        station.flap_point_idxs = Int64[resolve_ref(ref, point_names_dict,
            "point") for ref in station.flap_point_refs]
    end

    # Beam-anchored points: resolve the tube ref, derive beam_frac + offset.
    for point in points
        point.tube_idx = resolve_ref(point.tube_ref, tube_names_dict, "tube")
        point.tube_idx == 0 && continue
        tube = tubes[point.tube_idx]
        tube.model isa TimoshenkoTube || error("Point $(point.name) rides tube " *
            "$(tube.name), which is not a TimoshenkoTube beam element.")
        derive_point_beam_anchor!(point, tube, bodies)
    end

    # Particle wings carry mass on their (tube-connected) section bodies.
    body_groups = KiteGeometry.connected_body_groups(length(bodies), placement_tubes(tubes))
    for wing in wings
        wing.dynamics_type == PARTICLE_DYNAMICS &&
            finalize_particle_wing_mass!(
                wing, stations, points, bodies, set, body_groups)
    end

    wing_bodies = filter(is_wing, bodies)
    sys_struct = (; name,
        points=NamedCollection{Point}(points, point_names_dict),
        stations=NamedCollection{Station}(stations, station_names_dict),
        segments=NamedCollection{SegmentRow}(segments, segment_names_dict),
        pulleys=NamedCollection{Pulley}(pulleys, pulley_names_dict),
        tethers=NamedCollection{Tether}(tethers, tether_names_dict),
        winches=NamedCollection{Winch}(winches, winch_names_dict),
        bodies=NamedCollection{Body}(bodies, rigid_body_names_dict),
        wings=NamedCollection{Body}(wing_bodies, build_name_dict(wing_bodies)),
        tubes=NamedCollection{Tube}(tubes, tube_names_dict),
        materials)
    placed || place_structure!(sys_struct; ignore_l0, prn)
    init_derived_properties!(sys_struct; prn)

    # Panel→flap station map (structural; needs placed bodies + built panels).
    for wing in sys_struct.wings
        build_panel_station_map!(wing.aero, wing, sys_struct)
    end
    system = isnothing(definition) ? system_definition(sys_struct) : definition
    foreach(wing -> detach_vsm_engine!(wing, vsm_set), sys_struct.wings)
    return system
end
