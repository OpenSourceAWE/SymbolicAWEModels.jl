# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

"""
Helper functions for REFINE wing type that applies VSM panel forces directly to
structural points.
"""

# Baseline chord-based aero scaling for REFINE wings.
# Effective multiplier = 1 + (wing.aero_scale_chord or default below).
const AERO_SCALE_CHORD = 0.0

"""
    identify_wing_segments(wing_points::AbstractVector{Point})

Identify wing segments (LE/TE pairs) from WING-type points.

Assumes points are organized in pairs along the span, with even-numbered points
being leading edge and odd-numbered being trailing edge (or vice versa).

# Arguments
- `wing_points::AbstractVector{Point}`: All WING-type points for a wing (sorted by index)

# Returns
- `Vector{Tuple{Int64, Int64}}`: Vector of (le_point_idx, te_point_idx) pairs defining segments
"""
function identify_wing_segments(wing_points::AbstractVector{Point})
    # Sort points by index to ensure consistent ordering
    sorted_points = sort(wing_points, by=p->p.idx)

    n_points = length(sorted_points)
    @assert n_points % 2 == 0 "Wing must have even number of points (LE/TE pairs)"

    n_segments = n_points ÷ 2
    segments = Tuple{Int64, Int64}[]

    # Group consecutive pairs: (point[1], point[2]), (point[3], point[4]), ...
    for i in 1:n_segments
        le_idx = 2*i - 1
        te_idx = 2*i
        le_point = sorted_points[le_idx]
        te_point = sorted_points[te_idx]

        # Determine which is LE and which is TE by x-coordinate (LE has smaller x)
        if le_point.pos_cad[1] < te_point.pos_cad[1]
            push!(segments, (le_point.idx, te_point.idx))
        else
            push!(segments, (te_point.idx, le_point.idx))
        end
    end

    return segments
end

"""
    build_point_to_vsm_point_mapping(wing_points::AbstractVector{Point}, vsm_wing::VortexStepMethod.AbstractWing)

Build 1:1 mapping from structural WING points to VSM wing section points (LE/TE) using closest-point distance.

For each VSM section point (LE/TE), finds the closest structural point in CAD frame.

# Constraint
Requires: `length(wing_points) == 2 * length(vsm_wing.unrefined_sections)`

# Arguments
- `wing_points::AbstractVector{Point}`: Structural WING-type points
- `vsm_wing::VortexStepMethod.AbstractWing`: VSM wing with sections

# Returns
- `Dict{Int64, Tuple{Int64, Symbol}}`: Mapping structural_point_idx -> (section_idx, :LE or :TE)

# Algorithm
1. For each section in vsm_wing.sections:
   - Find closest unused structural point to section.LE_point → assign to (section_idx, :LE)
   - Find closest unused structural point to section.TE_point → assign to (section_idx, :TE)
2. Distance measured in CAD/body frame using norm(point.pos_cad - section_point)
"""
function build_point_to_vsm_point_mapping(
    wing_points::AbstractVector{Point},
    vsm_wing::VortexStepMethod.AbstractWing
)
    n_points = length(wing_points)
    n_sections = length(vsm_wing.unrefined_sections)

    # Validate 1:1 correspondence constraint
    if n_points != 2 * n_sections
        error("REFINE wing requires n_structural_points ($(n_points)) == " *
              "2 * n_vsm_sections ($(n_sections))")
    end

    point_to_vsm_point = Dict{Int64, Tuple{Int64, Symbol}}()
    used_points = Set{Int64}()

    for (section_idx, section) in enumerate(vsm_wing.unrefined_sections)
        # Map LE_point to closest unused structural point
        le_pos = section.LE_point
        min_dist = Inf
        closest_le_idx = wing_points[1].idx

        for point in wing_points
            if point.idx in used_points
                continue  # Already assigned
            end
            dist = norm(point.pos_cad - le_pos)
            if dist < min_dist
                min_dist = dist
                closest_le_idx = point.idx
            end
        end

        point_to_vsm_point[closest_le_idx] = (Int64(section_idx), :LE)
        push!(used_points, closest_le_idx)

        # Map TE_point to closest unused structural point
        te_pos = section.TE_point
        min_dist = Inf
        closest_te_idx = wing_points[1].idx

        for point in wing_points
            if point.idx in used_points
                continue  # Already assigned
            end
            dist = norm(point.pos_cad - te_pos)
            if dist < min_dist
                min_dist = dist
                closest_te_idx = point.idx
            end
        end

        point_to_vsm_point[closest_te_idx] = (Int64(section_idx), :TE)
        push!(used_points, closest_te_idx)
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
    d = r_le_rel - r_te_rel  # chord direction (LE→TE)
    d_norm_sq = dot(d, d)
    if d_norm_sq < 1e-12
        # Degenerate chord: just split equally and spanwise-preserve the torque
        F_le = 0.5 * F_panel
        F_te = 0.5 * F_panel
    else
        # Minimum-norm split that preserves moment about r_ref
        w_le = clamp(-dot(r_te_rel, d) / d_norm_sq, 0.0, 1.0)
        r_weighted = w_le * r_le_rel + (1 - w_le) * r_te_rel
        M_cp = M_panel - cross(r_ac - r_ref, F_panel)
        M_target = M_cp - cross(r_weighted, F_panel)
        ΔF = cross(M_target, d) / d_norm_sq
        F_le = w_le * F_panel + ΔF
        F_te = (1 - w_le) * F_panel - ΔF
    end

    # Spanwise split preserving moment about r_ref
    span_split = function (F::SVector{3}, r_mid::SVector{3}, r_left::SVector{3}, r_right::SVector{3})
        a = cross(r_left - r_ref, F)
        b = cross(r_right - r_ref, F)
        m_target = cross(r_mid - r_ref, F)
        ab = a - b
        denom = dot(ab, ab)
        w = denom < 1e-14 ? 0.5 : clamp(dot(ab, m_target - b) / denom, 0.0, 1.0)
        return (w * F, (1 - w) * F)
    end

    F_le_left, F_le_right = span_split(F_le, r_le_mid, nodes[1], nodes[2])
    # nodes[3] = TE right, nodes[4] = TE left → pass left/right accordingly
    F_te_left, F_te_right = span_split(F_te, r_te_mid, nodes[4], nodes[3])

    return (F_le_left, F_le_right, F_te_right, F_te_left, nodes)
end

"""
    distribute_panel_forces_to_points!(wing::VSMWing, points::AbstractVector{Point})

Distribute VSM forces to structural points using refined panel forces.

After VSM solve, each refined panel force/moment is split into corner-node
forces (moment-preserving about the chosen reference) and then aggregated to
the structural LE/TE points of the parent section (1:1 mapping).

# Algorithm
1. Initialize all WING point aero_forces to zero
2. Build inverse mapping from section → LE/TE structural point indices
3. For each refined panel of this wing:
   - Get panel force/moment from solver solution (body frame)
   - Map panel to its parent section using `refined_panel_mapping`
   - Split to LE/TE forces with `compute_aerostruc_loads`
   - Accumulate forces at the corresponding structural points

# Arguments
- `wing::VSMWing`: Wing with REFINE type and solved VSM state
- `points::AbstractVector{Point}`: All structural points (will filter for WING type)
"""
function distribute_panel_forces_to_points!(wing::VSMWing, points::AbstractVector{Point})
    @assert wing.wing_type == REFINE "Can only distribute forces for REFINE wings"

    sol = wing.vsm_solver.sol
    panels = wing.vsm_aero.panels
    panel_to_section = wing.vsm_wing.refined_panel_mapping
    f_body = sol.f_body_3D
    m_body = sol.m_body_3D

    # Initialize all WING point forces to zero
    for point in points
        if point.type == WING && point.wing_idx == wing.idx
            point.aero_force_b .= 0.0
        end
    end

    # Build inverse mapping: (section_idx, :LE/:TE) -> point_idx
    vsm_point_to_struct = Dict{Tuple{Int64, Symbol}, Int64}()
    for (point_idx, (section_idx, le_or_te)) in wing.point_to_vsm_point
        vsm_point_to_struct[(section_idx, le_or_te)] = point_idx
    end

    # Determine offset of this wing's panels in the solver arrays
    start_idx = 1
    if hasproperty(wing.vsm_solver, :body_aero)
        for w in wing.vsm_solver.body_aero.wings
            w === wing && break
            start_idx += length(w.vsm_aero.panels)
        end
    end

    # For each refined panel, split force/moment to LE/TE of its parent section
    n_panels_wing = length(panels)
    for local_panel_idx in 1:n_panels_wing
        panel_idx = start_idx + local_panel_idx - 1
        panel = panels[local_panel_idx]
        scale = 1.0 + (isfinite(wing.aero_scale_chord) ? wing.aero_scale_chord : AERO_SCALE_CHORD)
        Fp = scale .* SVector{3}(f_body[:, panel_idx])
        Mp = scale .* SVector{3}(m_body[:, panel_idx])

        section_idx = panel_to_section[local_panel_idx]
        section = wing.vsm_wing.unrefined_sections[section_idx]

        le_key = (Int64(section_idx), :LE)
        te_key = (Int64(section_idx), :TE)
        haskey(vsm_point_to_struct, le_key) || continue
        haskey(vsm_point_to_struct, te_key) || continue

        F_le_left, F_le_right, F_te_right, F_te_left, _ = compute_aerostruc_loads(panel, Fp, Mp)
        F_le = F_le_left + F_le_right
        F_te = F_te_left + F_te_right

        points[vsm_point_to_struct[le_key]].aero_force_b .+= F_le
        points[vsm_point_to_struct[te_key]].aero_force_b .+= F_te
    end

    return nothing
end

"""
    update_vsm_wing_from_structure!(wing::VSMWing, points::AbstractVector{Point})

Update VSM section points (LE/TE) directly from structural point positions using 1:1 mapping.

This creates two-way coupling: structural deformation → VSM sections → aero forces.

# Algorithm
Uses direct 1:1 correspondence between structural points and VSM section points:
1. For each structural WING point:
   - Calculate current position in body frame: pos_b = R_b_to_w' * (pos_w - origin)
   - Find corresponding VSM section point (LE or TE) via wing.point_to_vsm_point
   - Set VSM section point directly: section.LE_point = pos_b (or TE_point)

# Notes
- Section points are stored in body frame coordinates
- `wing.R_b_to_w` and `wing.pos_w` are updated each timestep from structural geometry (symbolic equations)
- To get world coordinates: `world_pos = wing.R_b_to_w * section.LE_point + wing.pos_w`

# Arguments
- `wing::VSMWing`: Wing with REFINE type
- `points::AbstractVector{Point}`: All structural points (will filter for WING type)
"""
function update_vsm_wing_from_structure!(wing::VSMWing, points::AbstractVector{Point})
    @assert wing.wing_type == REFINE "Can only update wing geometry for REFINE wings"

    # Get current R_b_to_w and origin from wing state
    # (These are updated during simulation from structural geometry)
    R_b_to_w = wing.R_b_to_w
    origin = wing.pos_w

    # Update each VSM section point directly from its corresponding structural point
    for (point_idx, (section_idx, le_or_te)) in wing.point_to_vsm_point
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
    # Do NOT call reinit! on wing - only modify sections!
    # body_aero reinit! will update panels from modified sections (called in update_vsm!)
    return nothing
end
