# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

"""
Helper functions for VSM wing types.

PARTICLE_DYNAMICS-specific functions (panel force distribution, structural geometry
updates) are at the bottom of this file.  The shared
`match_aero_sections_to_structure!` works for all VSMWing types.
"""

"""
    const AERO_SCALE_CHORD = 0.0

Baseline chord-based aero scaling for PARTICLE_DYNAMICS wings; effective
multiplier is `1 + (wing.aero_scale_chord or this default)`.
"""
const AERO_SCALE_CHORD = 0.0

"""
    identify_wing_segments(wing_points; twist_surfaces=nothing, wing_twist_surface_idxs=nothing)

Identify wing segments (LE/TE pairs) from wing nodes.

When `twist_surfaces` and `wing_twist_surface_idxs` are provided, uses twist_surface `point_idxs`
to determine LE (`point_idxs[1]`) and TE (`point_idxs[end]`) for each
section. Falls back to a consecutive-pair heuristic (sorted by point index)
when twist_surfaces are unavailable.

In both paths an x-coordinate check swaps LE/TE if needed (LE has
smaller `pos_cad[1]`).

# Arguments
- `wing_points::AbstractVector{Point}`: wing nodes for a wing.

# Keyword Arguments
- `twist_surfaces::Union{Nothing, AbstractVector{TwistSurface}}`: All twist_surfaces in the
  system (indexed by `wing_twist_surface_idxs`).
- `wing_twist_surface_idxs::Union{Nothing, AbstractVector{<:Integer}}`:
  Indices into `twist_surfaces` belonging to this wing.

# Returns
- `Vector{Tuple{Int64, Int64}}`: (le_point_idx, te_point_idx) pairs.
"""
function identify_wing_segments(
    wing_points::AbstractVector{Point};
    twist_surfaces::AbstractVector{TwistSurface}=TwistSurface[],
    wing_twist_surface_idxs::AbstractVector{<:Integer}=Int[]
)
    use_twist_surfaces = !isempty(twist_surfaces) &&
        !isempty(wing_twist_surface_idxs)

    if use_twist_surfaces
        segments = Tuple{Int64, Int64}[]
        for g_idx in wing_twist_surface_idxs
            twist_surface = twist_surfaces[g_idx]
            length(twist_surface.point_idxs) >= 2 || error(
                "TwistSurface $(twist_surface.name): need at least " *
                "2 point_idxs (LE/TE), got " *
                "$(length(twist_surface.point_idxs))")
            station = [wing_points[findfirst(p -> p.idx == idx, wing_points)]
                       for idx in twist_surface.point_idxs]
            le_point = argmin(p -> p.pos_cad[1], station)
            te_point = argmax(p -> p.pos_cad[1], station)
            push!(segments, (le_point.idx, te_point.idx))
        end
        return segments
    end

    # Fallback: consecutive-pair heuristic
    sorted_points = sort(wing_points, by=p->p.idx)

    n_points = length(sorted_points)
    @assert n_points % 2 == 0 (
        "Wing must have even number of points " *
        "(LE/TE pairs)")

    n_segments = n_points ÷ 2
    segments = Tuple{Int64, Int64}[]

    for i in 1:n_segments
        le_idx = 2*i - 1
        te_idx = 2*i
        le_point = sorted_points[le_idx]
        te_point = sorted_points[te_idx]

        if le_point.pos_cad[1] < te_point.pos_cad[1]
            push!(segments, (le_point.idx, te_point.idx))
        else
            push!(segments, (te_point.idx, le_point.idx))
        end
    end

    return segments
end

"""
    match_aero_sections_to_structure!(wing, points; twist_surfaces)

Reconcile a wing's aerodynamic sections with its structural geometry.

RIGID_DYNAMICS wings own their aero panel geometry (mesh- or
YAML-defined) and keep it; only the twist_surface→section mapping
(`wing.wing_segments`) is recorded. PARTICLE_DYNAMICS wings deform with
their structural points, so each unrefined section is rebuilt onto its
structural LE/TE pair: a 1:1 copy when counts match, otherwise
`use_prior_polar` and existing `refined_sections` are required to
preserve polars.

# Keyword Arguments
- `twist_surfaces::AbstractVector{TwistSurface}`: TwistSurfaces used for LE/TE identification
  via [`identify_wing_segments`](@ref).
"""
function match_aero_sections_to_structure!(
    wing::Body,
    points::AbstractVector{Point};
    twist_surfaces::AbstractVector{TwistSurface}=TwistSurface[]
)
    wing_points = [
        p for p in points if
        p.is_wing_node && p.wing_idx == wing.idx
    ]

    if wing.dynamics_type == RIGID_DYNAMICS
        wing.wing_segments = identify_wing_segments(
            wing_points; twist_surfaces=twist_surfaces,
            wing_twist_surface_idxs=wing.twist_surface_idxs)
        return nothing
    end

    wing_twist_surface_idxs = wing.twist_surface_idxs
    has_twist_surfaces = !isempty(twist_surfaces) &&
        !isempty(wing_twist_surface_idxs)

    if has_twist_surfaces
        n_struct_sections = length(wing_twist_surface_idxs)
        for g_idx in wing_twist_surface_idxs
            twist_surface = twist_surfaces[g_idx]
            length(twist_surface.point_idxs) >= 2 || error(
                "PARTICLE_DYNAMICS wing $(wing.idx): twist_surface " *
                "$(twist_surface.name) needs at least 2 points " *
                "(LE = point_idxs[1], TE = point_idxs[end]), got " *
                "$(length(twist_surface.point_idxs))")
        end
    else
        n_points = length(wing_points)
        n_points % 2 == 0 || error(
            "Wing $(wing.idx): no twist_surfaces and odd " *
            "number of wing nodes " *
            "($(n_points)). Define twist_surfaces to " *
            "specify LE/TE pairs.")
        n_struct_sections = n_points ÷ 2
    end

    n_aero_sections =
        length(wing.vsm_wing.unrefined_sections)
    counts_differ = n_struct_sections != n_aero_sections

    if counts_differ
        wing.vsm_wing.use_prior_polar || error(
            "Wing $(wing.idx): structural sections " *
            "($(n_struct_sections)) do not match " *
            "aerodynamic sections " *
            "($(n_aero_sections)). Set " *
            "use_prior_polar=true to rebuild " *
            "unrefined sections from structural " *
            "geometry."
        )

        isempty(wing.vsm_wing.refined_sections) &&
            error(
                "Wing $(wing.idx): cannot rebuild " *
                "unrefined sections because no " *
                "refined sections exist to " *
                "preserve polars from."
            )
    end

    wing_segments = identify_wing_segments(
        wing_points; twist_surfaces=twist_surfaces,
        wing_twist_surface_idxs=wing_twist_surface_idxs)
    wing.wing_segments = wing_segments
    length(wing_segments) == n_struct_sections || error(
        "Wing $(wing.idx): failed to identify " *
        "structural LE/TE pairs."
    )

    original_sections = wing.vsm_wing.unrefined_sections
    n_original = length(original_sections)
    n_original > 0 || error(
        "Wing $(wing.idx): aerodynamic geometry " *
        "has zero unrefined sections."
    )
    R_b_to_c = wing.R_b_to_c
    origin_cad = wing.pos_cad
    new_sections = Vector{VortexStepMethod.Section}(
        undef, n_struct_sections)

    for (i, (le_idx, te_idx)) in enumerate(wing_segments)
        source_idx = if counts_differ
            n_struct_sections == 1 ? 1 :
                round(Int,
                    1 + (i - 1) * (n_original - 1) /
                    (n_struct_sections - 1))
        else
            i  # 1:1 copy when counts match
        end
        source_section = original_sections[source_idx]

        le_body = R_b_to_c' *
            (points[le_idx].pos_cad - origin_cad)
        te_body = R_b_to_c' *
            (points[te_idx].pos_cad - origin_cad)

        section = VortexStepMethod.Section()
        VortexStepMethod.reinit!(
            section, le_body, te_body,
            source_section.aero_model,
            source_section.aero_data,
            source_section.section_aero
        )
        new_sections[i] = section
    end

    wing.vsm_wing.unrefined_sections = new_sections
    wing.vsm_wing.n_unrefined_sections =
        Int16(n_struct_sections)

    refine!(wing.vsm_wing;
        recompute_mapping=true, sort_sections=false)
    VortexStepMethod.reinit!(wing.vsm_aero)

    return nothing
end

"""
    build_point_to_vsm_point_mapping(wing_segments)

Invert the twist-surface-derived `wing_segments` — `(le_point_idx, te_point_idx)`
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
    distribute_panel_forces_to_points!(wing::Body, points::AbstractVector{Point})

Distribute VSM forces to structural points using refined panel forces.

After the VSM solve, each refined panel's body-frame force/moment is split into
corner-node forces ([`compute_aerostruc_loads`](@ref), moment-preserving about the
chosen reference) and accumulated on the structural LE/TE points of its parent
section, which `refined_panel_mapping` names (1:1).

# Arguments
- `wing::Body`: Wing with PARTICLE_DYNAMICS type and solved VSM state
- `points::AbstractVector{Point}`: All structural points (filtered to this wing's nodes)
"""
function distribute_panel_forces_to_points!(wing::Body, points::AbstractVector{Point})
    @assert wing.dynamics_type == PARTICLE_DYNAMICS "Can only distribute forces for PARTICLE_DYNAMICS wings"

    sol = wing.vsm_solver.sol
    panels = wing.vsm_aero.panels
    panel_to_section = wing.vsm_wing.refined_panel_mapping
    f_body = sol.f_body_3D
    m_body = sol.m_body_3D

    # Initialize all wing-node forces to zero
    for point in points
        if point.is_wing_node && point.wing_idx == wing.idx
            point.aero_force_b .= 0.0
        end
    end

    # Build inverse mapping: (section_idx, :LE/:TE) -> point_idx
    point_to_vsm_point =
        wing.point_to_vsm_point::Union{Nothing,
            Dict{Int64, Tuple{Int64, Symbol}}}
    isnothing(point_to_vsm_point) && error(
        "PARTICLE_DYNAMICS wing $(wing.idx) missing point_to_vsm_point mapping")
    point_to_vsm_point =
        point_to_vsm_point::Dict{Int64,
            Tuple{Int64, Symbol}}
    vsm_point_to_struct = Dict{Tuple{Int64, Symbol}, Int64}()
    for (point_idx, (section_idx, le_or_te)) in point_to_vsm_point
        vsm_point_to_struct[(section_idx, le_or_te)] = point_idx
    end

    # Determine offset of this wing's panels in the solver arrays
    start_idx = 1
    if hasproperty(wing.vsm_solver, :body_aero)
        for other_wing in wing.vsm_solver.body_aero.wings
            other_wing === wing && break
            start_idx += length(other_wing.vsm_aero.panels)
        end
    end

    # For each refined panel, split force/moment to LE/TE of its parent section
    n_panels_wing = length(panels)
    for local_panel_idx in 1:n_panels_wing
        panel_idx = start_idx + local_panel_idx - 1
        panel = panels[local_panel_idx]
        scale = 1.0 + (isfinite(wing.aero_scale_chord) ? wing.aero_scale_chord : AERO_SCALE_CHORD)
        panel_force = scale .* SVector{3}(f_body[:, panel_idx])
        panel_moment = scale .* SVector{3}(m_body[:, panel_idx])

        section_idx = panel_to_section[local_panel_idx]

        le_key = (Int64(section_idx), :LE)
        te_key = (Int64(section_idx), :TE)
        haskey(vsm_point_to_struct, le_key) || continue
        haskey(vsm_point_to_struct, te_key) || continue

        F_le_left, F_le_right, F_te_right, F_te_left, _ =
            compute_aerostruc_loads(panel, panel_force, panel_moment)
        F_le = F_le_left + F_le_right
        F_te = F_te_left + F_te_right

        points[vsm_point_to_struct[le_key]].aero_force_b .+= F_le
        points[vsm_point_to_struct[te_key]].aero_force_b .+= F_te
    end

    return nothing
end

"""
    update_vsm_wing_from_structure!(wing::Body, points::AbstractVector{Point})

Update VSM section points (LE/TE) directly from structural point positions, closing
the two-way coupling structural deformation → VSM sections → aero forces.

Each structural wing node maps 1:1 onto a VSM section point through
`wing.point_to_vsm_point`, and writes it its body-frame position
`pos_b = R_b_to_w' * (pos_w - origin)`.

# Notes
- Section points are stored in body frame coordinates
- `wing.R_b_to_w` and `wing.pos_w` are updated each timestep from structural geometry (symbolic equations)
- To get world coordinates: `world_pos = wing.R_b_to_w * section.LE_point + wing.pos_w`

# Arguments
- `wing::Body`: Wing with PARTICLE_DYNAMICS type
- `points::AbstractVector{Point}`: All structural points (filtered to this wing's nodes)
"""
function update_vsm_wing_from_structure!(wing::Body, points::AbstractVector{Point})
    @assert wing.dynamics_type == PARTICLE_DYNAMICS "Can only update wing geometry for PARTICLE_DYNAMICS wings"

    # R_b_to_w and origin are updated during simulation from structural geometry.
    R_b_to_w = wing.R_b_to_w::Matrix{SimFloat}
    origin = wing.pos_w::KVec3

    # Update each VSM section point directly from its corresponding structural point
    point_to_vsm_point =
        wing.point_to_vsm_point::Union{Nothing,
            Dict{Int64, Tuple{Int64, Symbol}}}
    isnothing(point_to_vsm_point) && error(
        "PARTICLE_DYNAMICS wing $(wing.idx) missing point_to_vsm_point mapping")
    point_to_vsm_point =
        point_to_vsm_point::Dict{Int64,
            Tuple{Int64, Symbol}}
    for (point_idx, (section_idx, le_or_te)) in point_to_vsm_point
        point = points[point_idx]

        # Calculate current position in body frame
        pos_b = R_b_to_w' * (point.pos_w - origin)

        # Get the section
        section = wing.vsm_wing.unrefined_sections[section_idx]

        # Set section point directly to body frame position
        if le_or_te == :LE
            section.LE_point .= pos_b
        else  # :TE
            section.TE_point .= pos_b
        end
    end

    refine!(wing.vsm_wing; recompute_mapping=false, sort_sections=false)
    VortexStepMethod.reinit!(wing.vsm_aero)
    # Do NOT reinit! the wing; body_aero reinit! updates panels in refresh_aero!.
    return nothing
end

"""
    compute_twist_surface_geometry!(wing, twist_surfaces, points)

For each of `wing`'s twist surfaces with an unset chord, derive its leading-edge
position, chord vector, and spanwise airfoil axis from the nearest VSM refined
section (body frame). Fills the chord for any twist surface left unset.
"""
function compute_twist_surface_geometry!(wing, twist_surfaces, points)
    for twist_surface_idx in wing.twist_surface_idxs
        twist_surface = twist_surfaces[twist_surface_idx]
        iszero(twist_surface.chord) || continue
        center = zeros(3)
        for pt_idx in twist_surface.point_idxs
            center .+= wing.R_b_to_c' *
                (points[pt_idx].pos_cad - wing.pos_cad)
        end
        center ./= length(twist_surface.point_idxs)

        sections = wing.vsm_wing.refined_sections
        n_sec = length(sections)
        offset_vec = [0.0, 0.0, wing.aero_z_offset]
        ksec = argmin([
            norm(center -
                ((Vector(section.LE_point) +
                  Vector(section.TE_point)) / 2 .- offset_vec))
            for section in sections])
        le_sec = Vector(sections[ksec].LE_point)
        te_sec = Vector(sections[ksec].TE_point)
        span_dir = zeros(3)
        ksec > 1 && (span_dir += normalize(
            Vector(sections[ksec - 1].LE_point) - le_sec))
        ksec < n_sec && (span_dir += normalize(
            le_sec - Vector(sections[ksec + 1].LE_point)))

        twist_surface.le_pos .= le_sec
        twist_surface.chord .= te_sec - le_sec
        twist_surface.y_airf .= normalize(span_dir)
    end
    return nothing
end

"""
    setup_particle_point_mapping!(wing, points, twist_surfaces)

For a VSM `PARTICLE_DYNAMICS` `wing`, build the structural↔panel point mapping and
LE/TE `wing_segments` if not already set. Errors if the required body-frame
`z_ref_points`/`y_ref_points` are missing.
"""
function setup_particle_point_mapping!(wing, points, twist_surfaces)
    wing_point_idxs = findall(
        point -> point.is_wing_node && point.wing_idx == wing.idx, points)
    wing_pts = [points[idx] for idx in wing_point_idxs]
    if isnothing(wing.wing_segments)
        wing.wing_segments = identify_wing_segments(wing_pts;
            twist_surfaces=twist_surfaces,
            wing_twist_surface_idxs=wing.twist_surface_idxs)
    end
    if isnothing(wing.point_to_vsm_point)
        wing.point_to_vsm_point =
            build_point_to_vsm_point_mapping(wing.wing_segments)
    end
    isnothing(wing.z_ref_points) && error(
        "PARTICLE_DYNAMICS wing '$(wing.name)': z_ref_points must be specified")
    isnothing(wing.y_ref_points) && error(
        "PARTICLE_DYNAMICS wing '$(wing.name)': y_ref_points must be specified")
    return nothing
end
