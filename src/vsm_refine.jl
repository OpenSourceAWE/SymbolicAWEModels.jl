# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

"""
    const AERO_SCALE_CHORD = 0.0

Baseline chord-based aero scaling for PARTICLE_DYNAMICS wings; effective
multiplier is `1 + aero_scale_chord`, the engine's or else this default.
"""
const AERO_SCALE_CHORD = 0.0

"""
    identify_wing_segments(sam, wing_idx) -> Vector{Tuple{Int64, Int64}}

The `(le_point_idx, te_point_idx)` pair of each section of wing `wing_idx` of `sam`:
of each of its stations' points, the one with the smallest and the one with the largest
body-frame x in the wing's initial pose ([`initial_chord_position`](@ref)). A wing
without stations pairs its nodes ([`wing_points`](@ref)) two by two in point order, the
LE again the one with the smaller x.
"""
function identify_wing_segments(sam, wing_idx)
    body_idx = wing_body(sam, wing_idx)
    points = sam.definition.points
    chord_position(i) = initial_chord_position(sam, body_idx, points[i].pos_ENU)
    station_idxs = wing_station_idxs(sam, wing_idx)
    if !isempty(station_idxs)
        segments = Tuple{Int64, Int64}[]
        for station_idx in station_idxs
            station = sam.definition.spanwise_stations[station_idx]
            length(station.points) >= 2 || error(
                "Station $(station.name): need at least 2 point_idxs (LE/TE), got " *
                "$(length(station.points))")
            push!(segments, (argmin(chord_position, station.points),
                             argmax(chord_position, station.points)))
        end
        return segments
    end

    point_idxs = wing_points(sam, wing_idx)
    @assert iseven(length(point_idxs)) "Wing must have even number of points (LE/TE pairs)"
    segments = Tuple{Int64, Int64}[]
    for k in 1:2:length(point_idxs)
        le_idx, te_idx = point_idxs[k], point_idxs[k + 1]
        chord_position(le_idx) < chord_position(te_idx) ?
            push!(segments, (le_idx, te_idx)) : push!(segments, (te_idx, le_idx))
    end
    return segments
end

"""
    check_span_order(mode, sam, wing_idx)

Error unless the unrefined VSM sections of the aero `mode` of wing `wing_idx` of `sam`
run from +y to -y: the first section's leading edge lies further along
`spanwise_direction` than the last one's.
"""
function check_span_order(mode, sam, wing_idx)
    sections = mode.vsm_wing.unrefined_sections
    span = first(sections).LE_point - last(sections).LE_point
    dot(span, mode.vsm_wing.spanwise_direction) >= 0 || error(
        "Wing $(sam.definition.wings[wing_idx].name): its aerodynamic sections run " *
        "from -y to +y. Order the sections and the structural stations from +y to -y.")
    return nothing
end

"""
    match_aero_sections_to_structure!(mode, sam, wing_idx)

Reconcile the aerodynamic sections of the aero `mode` of wing `wing_idx` of `sam` with
its structural geometry.

RIGID_DYNAMICS wings own their aero panel geometry (mesh- or YAML-defined) and keep it;
only the station→section mapping (`wing_segments`) is recorded. PARTICLE_DYNAMICS wings
deform with their structural points, so each unrefined section is rebuilt onto its
structural LE/TE pair ([`identify_wing_segments`](@ref)): a 1:1 copy when counts match,
otherwise `use_prior_polar` and existing `refined_sections` are required to preserve
polars.
"""
function match_aero_sections_to_structure!(mode, sam, wing_idx)
    body_idx = wing_body(sam, wing_idx)
    if body_dynamics_type(sam, body_idx) == RIGID_DYNAMICS
        mode.wing_segments = identify_wing_segments(sam, wing_idx)
        return nothing
    end

    station_idxs = wing_station_idxs(sam, wing_idx)
    if isempty(station_idxs)
        n_points = length(wing_points(sam, wing_idx))
        n_points % 2 == 0 || error(
            "Wing $wing_idx: no stations and odd number of wing nodes ($n_points). " *
            "Define stations to specify LE/TE pairs.")
        n_struct_sections = n_points ÷ 2
    else
        n_struct_sections = length(station_idxs)
        for station_idx in station_idxs
            station = sam.definition.spanwise_stations[station_idx]
            length(station.points) >= 2 || error(
                "PARTICLE_DYNAMICS wing $wing_idx: station $(station.name) needs at " *
                "least 2 points (LE = point_idxs[1], TE = point_idxs[end]), got " *
                "$(length(station.points))")
        end
    end

    vsm_wing = mode.vsm_wing
    n_aero_sections = length(vsm_wing.unrefined_sections)
    counts_differ = n_struct_sections != n_aero_sections
    if counts_differ
        vsm_wing.use_prior_polar || error(
            "Wing $wing_idx: structural sections ($n_struct_sections) do not match " *
            "aerodynamic sections ($n_aero_sections). Set use_prior_polar=true to " *
            "rebuild unrefined sections from structural geometry.")
        isempty(vsm_wing.refined_sections) && error(
            "Wing $wing_idx: cannot rebuild unrefined sections because no refined " *
            "sections exist to preserve polars from.")
    end

    wing_segments = identify_wing_segments(sam, wing_idx)
    mode.wing_segments = wing_segments
    length(wing_segments) == n_struct_sections || error(
        "Wing $wing_idx: failed to identify structural LE/TE pairs.")

    original_sections = vsm_wing.unrefined_sections
    n_original = length(original_sections)
    n_original > 0 || error("Wing $wing_idx: aerodynamic geometry has zero unrefined " *
                            "sections.")
    new_sections = Vector{VortexStepMethod.Section}(undef, n_struct_sections)
    points = sam.definition.points
    for (i, (le_idx, te_idx)) in enumerate(wing_segments)
        source_idx = if counts_differ
            n_struct_sections == 1 ? 1 :
                round(Int, 1 + (i - 1) * (n_original - 1) / (n_struct_sections - 1))
        else
            i
        end
        source_section = original_sections[source_idx]
        le_body = initial_body_position(sam, body_idx, points[le_idx].pos_ENU)
        te_body = initial_body_position(sam, body_idx, points[te_idx].pos_ENU)
        section = VortexStepMethod.Section()
        VortexStepMethod.reinit!(section, le_body, te_body, source_section.aero_model,
                                 source_section.aero_data, source_section.section_aero)
        new_sections[i] = section
    end

    vsm_wing.unrefined_sections = new_sections
    vsm_wing.n_unrefined_sections = Int16(n_struct_sections)
    check_span_order(mode, sam, wing_idx)
    refine!(vsm_wing; recompute_mapping=true, sort_sections=false)
    VortexStepMethod.reinit!(mode.vsm_aero)
    return nothing
end

"""
    build_point_to_vsm_point_mapping(wing_segments)

Invert the station-derived `wing_segments` — `(le_point_idx, te_point_idx)`
per unrefined section — into the `structural point → (section_idx, :LE/:TE)` map.

Each unrefined section contributes exactly its LE and TE structural points. Any
interior points on the same station (chord control points on the beam) are not
aero station points and are simply absent from the map; there is no constraint on
the total wing-node count.
"""
function build_point_to_vsm_point_mapping(
    wing_segments::AbstractVector{Tuple{Int64, Int64}},
)
    point_to_vsm_point = Dict{Int64, Tuple{Int64, Symbol}}()
    for (section_idx, (le_idx, te_idx)) in enumerate(wing_segments)
        point_to_vsm_point[le_idx] = (Int64(section_idx), :LE)
        point_to_vsm_point[te_idx] = (Int64(section_idx), :TE)
    end
    return point_to_vsm_point
end

# Distribute a panel force/moment to four corner nodes while preserving force and moment
function compute_aerostruc_loads(panel, F_panel::SVector{3}, M_panel::SVector{3};
    reference_point::SVector{3}=SVector(0.0, 0.0, 0.0))

    # Approximate quad corners from aero center, chord, width, and local axes
    c_vec = panel.x_airf
    s_vec = panel.y_airf
    chord = panel.chord
    width = panel.width
    r_ac = SVector{3}(panel.aero_center)
    # Assume aero center at quarter-chord midspan
    r_le_mid = r_ac - 0.25 * chord * c_vec
    r_te_mid = r_ac + 0.75 * chord * c_vec
    half_span = 0.5 * width * s_vec
    le_left = r_le_mid - half_span
    le_right = r_le_mid + half_span
    te_right = r_te_mid + half_span
    te_left = r_te_mid - half_span

    nodes = (le_left, le_right, te_right, te_left)
    r_ref = reference_point

    # Midpoints of LE/TE edges
    r_le_mid = 0.5 * (nodes[1] + nodes[2])
    r_te_mid = 0.5 * (nodes[3] + nodes[4])

    # Relative positions
    r_le_rel = r_le_mid - r_ac
    r_te_rel = r_te_mid - r_ac
    chord_dir = r_le_rel - r_te_rel  # chord direction (LE→TE)
    chord_norm_sq = dot(chord_dir, chord_dir)
    if chord_norm_sq < 1e-12
        # Degenerate chord: just split equally and spanwise-preserve the torque
        F_le = 0.5 * F_panel
        F_te = 0.5 * F_panel
    else
        # Minimum-norm split that preserves moment about r_ref
        w_le = clamp(-dot(r_te_rel, chord_dir) / chord_norm_sq, 0.0, 1.0)
        r_weighted = w_le * r_le_rel + (1 - w_le) * r_te_rel
        M_cp = M_panel - cross(r_ac - r_ref, F_panel)
        M_target = M_cp - cross(r_weighted, F_panel)
        ΔF = cross(M_target, chord_dir) / chord_norm_sq
        F_le = w_le * F_panel + ΔF
        F_te = (1 - w_le) * F_panel - ΔF
    end

    # Spanwise split preserving moment about r_ref
    span_split = function (F::SVector{3}, r_mid::SVector{3}, r_left::SVector{3}, r_right::SVector{3})
        left_moment = cross(r_left - r_ref, F)
        right_moment = cross(r_right - r_ref, F)
        m_target = cross(r_mid - r_ref, F)
        moment_diff = left_moment - right_moment
        denom = dot(moment_diff, moment_diff)
        split_weight = denom < 1e-14 ? 0.5 :
            clamp(dot(moment_diff, m_target - right_moment) / denom, 0.0, 1.0)
        return (split_weight * F, (1 - split_weight) * F)
    end

    F_le_left, F_le_right = span_split(F_le, r_le_mid, nodes[1], nodes[2])
    # nodes[3] = TE right, nodes[4] = TE left → pass left/right accordingly
    F_te_left, F_te_right = span_split(F_te, r_te_mid, nodes[4], nodes[3])

    return (F_le_left, F_le_right, F_te_right, F_te_left, nodes)
end

"""
    vsm_point_to_struct(mode, wing_idx) -> Dict{Tuple{Int64, Symbol}, Int64}

The structural point of each `(section_idx, :LE/:TE)` of the VSM sections of the aero
`mode` of `PARTICLE_DYNAMICS` wing `wing_idx`: the inverse of its `point_to_vsm_point`.
"""
function vsm_point_to_struct(mode, wing_idx)
    point_to_vsm_point =
        mode.point_to_vsm_point::Union{Nothing, Dict{Int64, Tuple{Int64, Symbol}}}
    isnothing(point_to_vsm_point) && error(
        "PARTICLE_DYNAMICS wing $wing_idx missing point_to_vsm_point mapping")
    return Dict{Tuple{Int64, Symbol}, Int64}(
        key => point_idx for (point_idx, key) in point_to_vsm_point)
end

"""
    distribute_panel_forces_to_points!(mode, sam, wing_idx)

Distribute the VSM forces of `PARTICLE_DYNAMICS` wing `wing_idx` of `sam` to its
structural points' `aero_force`, from the refined panel forces of its solved aero
`mode`. Each refined panel's body-frame force/moment is split into corner-node forces
([`compute_aerostruc_loads`](@ref), moment-preserving about the chosen reference) and
accumulated on the structural LE/TE points of its parent section, which
`refined_panel_mapping` names (1:1); the wing's other nodes get none.
"""
function distribute_panel_forces_to_points!(mode, sam, wing_idx)
    body_idx = wing_body(sam, wing_idx)
    @assert(body_dynamics_type(sam, body_idx) == PARTICLE_DYNAMICS,
            "Can only distribute forces for PARTICLE_DYNAMICS wings")
    sol = mode.vsm_solver.sol
    panels = mode.vsm_aero.panels
    panel_to_section = mode.vsm_wing.refined_panel_mapping
    to_struct = vsm_point_to_struct(mode, wing_idx)
    R_b_to_w = body_rotation(sam, body_idx)
    scale = 1.0 + (isfinite(mode.aero_scale_chord) ? mode.aero_scale_chord :
                   AERO_SCALE_CHORD)

    forces = Dict(i => zero(SVec3) for i in wing_points(sam, wing_idx))
    for (panel_idx, panel) in enumerate(panels)
        section_idx = Int64(panel_to_section[panel_idx])
        le_key, te_key = (section_idx, :LE), (section_idx, :TE)
        (haskey(to_struct, le_key) && haskey(to_struct, te_key)) || continue
        panel_force = scale .* SVector{3}(sol.f_body_3D[:, panel_idx])
        panel_moment = scale .* SVector{3}(sol.m_body_3D[:, panel_idx])
        F_le_left, F_le_right, F_te_right, F_te_left, _ =
            compute_aerostruc_loads(panel, panel_force, panel_moment)
        le_idx, te_idx = to_struct[le_key], to_struct[te_key]
        forces[le_idx] = get(forces, le_idx, zero(SVec3)) +
            R_b_to_w * (F_le_left + F_le_right)
        forces[te_idx] = get(forces, te_idx, zero(SVec3)) +
            R_b_to_w * (F_te_left + F_te_right)
    end
    for (i, force) in forces
        set_run_value!(sam, Val(:points), i, Val(:aero_force), force)
    end
    return nothing
end

"""
    update_vsm_wing_from_structure!(mode, sam, wing_idx)

Update the VSM section points (LE/TE) of the aero `mode` of `PARTICLE_DYNAMICS` wing
`wing_idx` of `sam` from the positions of its structural points, closing the two-way
coupling structural deformation → VSM sections → aero forces. Each structural wing node
maps 1:1 onto a VSM section point through the engine's `point_to_vsm_point`, and writes
it its position in the wing's body frame, `R_b_to_w' * (pos_w - origin)`.
"""
function update_vsm_wing_from_structure!(mode, sam, wing_idx)
    body_idx = wing_body(sam, wing_idx)
    @assert(body_dynamics_type(sam, body_idx) == PARTICLE_DYNAMICS,
            "Can only update wing geometry for PARTICLE_DYNAMICS wings")
    point_to_vsm_point =
        mode.point_to_vsm_point::Union{Nothing, Dict{Int64, Tuple{Int64, Symbol}}}
    isnothing(point_to_vsm_point) && error(
        "PARTICLE_DYNAMICS wing $wing_idx missing point_to_vsm_point mapping")
    for (point_idx, (section_idx, le_or_te)) in point_to_vsm_point
        pos_b = point_position_b(sam, body_idx, point_idx)
        section = mode.vsm_wing.unrefined_sections[section_idx]
        le_or_te == :LE ? (section.LE_point .= pos_b) : (section.TE_point .= pos_b)
    end
    check_span_order(mode, sam, wing_idx)
    refine!(mode.vsm_wing; recompute_mapping=false, sort_sections=false)
    VortexStepMethod.reinit!(mode.vsm_aero)
    return nothing
end

"""
    compute_station_geometry!(sam, wing_idx)

For each station of wing `wing_idx` of `sam` whose row gives no chord, take its leading edge
from the body origin, chord vector and spanwise axis from the nearest VSM refined
section, body frame, into the wing's engine (`station_le_b`, `station_chord`,
`station_y_airf`); zero for the others.
"""
function compute_station_geometry!(sam, wing_idx)
    engine = sam.aero[wing_idx].engine
    station_idxs = wing_station_idxs(sam, wing_idx)
    engine.station_le_b = [zero(KVec3) for _ in station_idxs]
    engine.station_chord = [zero(KVec3) for _ in station_idxs]
    engine.station_y_airf = [zero(KVec3) for _ in station_idxs]
    sections = engine.vsm_wing.refined_sections
    n_sec = length(sections)
    offset_vec = [0.0, 0.0, engine.aero_z_offset]
    for (k, station_idx) in enumerate(station_idxs)
        isnothing(sam.definition.spanwise_stations[station_idx].chord) || continue
        center = station_center_b(sam, wing_idx, station_idx)
        ksec = argmin([norm(center - ((Vector(section.LE_point) +
                                       Vector(section.TE_point)) / 2 .- offset_vec))
                       for section in sections])
        le_sec = Vector(sections[ksec].LE_point)
        te_sec = Vector(sections[ksec].TE_point)
        span_dir = zeros(3)
        ksec > 1 && (span_dir += normalize(Vector(sections[ksec - 1].LE_point) - le_sec))
        ksec < n_sec &&
            (span_dir += normalize(le_sec - Vector(sections[ksec + 1].LE_point)))
        engine.station_le_b[k] .= le_sec
        engine.station_chord[k] .= te_sec - le_sec
        engine.station_y_airf[k] .= normalize(span_dir)
    end
    return nothing
end

"""
    setup_particle_point_mapping!(mode, sam, wing_idx)

Build the structural↔panel point mapping and LE/TE `wing_segments` of the VSM aero
`mode` of `PARTICLE_DYNAMICS` wing `wing_idx` of `sam` where not already set. Errors if
the wing's frame reference points are missing.
"""
function setup_particle_point_mapping!(mode, sam, wing_idx)
    isnothing(mode.wing_segments) &&
        (mode.wing_segments = identify_wing_segments(sam, wing_idx))
    isnothing(mode.point_to_vsm_point) &&
        (mode.point_to_vsm_point = build_point_to_vsm_point_mapping(mode.wing_segments))
    isnothing(wing_frame(sam, wing_idx)) && error("PARTICLE_DYNAMICS wing " *
        "'$(sam.definition.wings[wing_idx].name)': z_ref_points and y_ref_points must " *
        "be specified")
    return nothing
end
