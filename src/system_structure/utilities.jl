# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

"""
Utility functions for SystemStructure.

This file contains:
- Validation functions
- reinit! for SystemStructure
- State copy and update functions
- Damping setters
- Segment statistics
"""

# ==================== SHARED HELPERS ==================== #

"""
    segment_cad_length(segment::Segment, points)

Compute segment length from endpoint `pos_cad` positions.
"""
function segment_cad_length(segment::Segment, points)
    point1 = points[segment.point_idxs[1]]
    point2 = points[segment.point_idxs[2]]
    return norm(point1.pos_cad - point2.pos_cad)
end

"""
    segment_world_length(segment::Segment, points)

Compute segment length from endpoint `pos_w` positions.
"""
function segment_world_length(segment::Segment, points)
    point1 = points[segment.point_idxs[1]]
    point2 = points[segment.point_idxs[2]]
    return norm(point1.pos_w - point2.pos_w)
end




# ==================== VALIDATION ==================== #

"""
    validate_sys_struct(sys_struct::SystemStructure)

Validate a `SystemStructure` for common configuration errors.

This function checks for issues that can cause initialization failures or
numerical problems during simulation. It emits warnings for suspicious
configurations and throws assertions for definite errors.

# Validations Performed

## Point Validations
- NaN extra_mass (error)
- Negative extra_mass (warning)
- Non-positive total_mass for DYNAMIC points (error) - checked before NaN position
- NaN position (error) - often caused by zero mass

## Wing Validations
- Non-positive mass (error) - checked before NaN position
- Zero or near-zero principal inertia components on RIGID_DYNAMICS wings (error/warning)
- NaN inertia values (error)
- Empty twist_surface list for RIGID_DYNAMICS wings (warning)
- NaN position (error) - often caused by zero mass/inertia

## Winch Validations
- Zero or negative inertia_total (error)
- Very small inertia_total (warning)
- NaN inertia_total (error)
- Non-positive drum_radius (error)
- Non-positive gear_ratio (error)

## Segment Validations
- Unusual diameter outside (0, 1) m range (warning)
- Non-positive rest length l0 (error)
- Zero or negative stiffness (warning)
- Negative damping (warning)

## Pulley Validations
- Zero total length constraint (error)

## TwistSurface Validations
- Inconsistent moment_frac across twist_surfaces (error)
"""
function validate_sys_struct(sys_struct::SystemStructure)
    (; points, twist_surfaces, segments, pulleys, wings, winches) = sys_struct

    # ==================== POINT VALIDATIONS ==================== #
    for point in points
        # Check for NaN extra_mass
        if isnan(point.extra_mass)
            error("Point $(point.name) has NaN extra_mass")
        end

        # Warn about negative extra_mass (physically nonsensical but still works)
        if point.extra_mass < 0
            @warn "Point $(point.name) has negative extra_mass $(point.extra_mass) kg. " *
                  "This is physically nonsensical."
        end

        # Check mass before NaN position: NaN pos is often caused by zero mass.
        if point.type == DYNAMIC && point.total_mass <= 0
            error("Point $(point.name) has non-positive total_mass ($(point.total_mass)). " *
                  "This will cause division by zero in acceleration calculations.")
        end

        # Check for NaN position (often a symptom of zero mass)
        if any(isnan.(point.pos_w))
            error("Point $(point.name) has NaN position: pos_w = $(point.pos_w)")
        end
    end

    # ==================== WING VALIDATIONS ==================== #
    for wing in wings
        # Check mass/inertia before NaN position: NaN pos is often caused by zero mass.
        if wing.mass <= 0
            error("Wing $(wing.name) has non-positive mass ($(wing.mass)). " *
                  "This will cause division by zero in acceleration calculations.")
        end

        if wing.dynamics_type == RIGID_DYNAMICS
            I_b = wing.inertia_principal

            # Check for zero or suspiciously small inertia (before NaN checks)
            for i in 1:3
                if I_b[i] ≈ 0.0
                    error("Wing $(wing.name) has zero inertia component " *
                          "I_b[$i] = $(I_b[i]). " *
                          "All principal inertia components must be non-zero.")
                elseif I_b[i] < 1e-6
                    @warn "Wing $(wing.name) has very small inertia component " *
                          "I_b[$i] = $(I_b[i]) kg⋅m², may cause numerical issues"
                end
            end

            # Check for NaN inertia
            if any(isnan.(I_b))
                error("Wing $(wing.name) has NaN inertia: I_b = $I_b")
            end

            # AeroNone does not couple to sections, so missing twist_surfaces is fine.
            if isempty(wing.twist_surface_idxs) &&
               couples_to_sections(wing.aero)
                @warn "Wing $(wing.name) (RIGID_DYNAMICS)" *
                    " has no twist_surfaces"
            end
        end

        # Check for NaN position (often a symptom of zero mass/inertia)
        if any(isnan.(wing.pos_w))
            error("Wing $(wing.name) has NaN position: pos_w = $(wing.pos_w)")
        end
    end

    # ==================== WINCH VALIDATIONS ==================== #
    for winch in winches
        # Check for NaN inertia
        if isnan(winch.inertia_total)
            error("Winch $(winch.name) has NaN inertia_total")
        end

        # Check for zero or negative inertia
        if winch.inertia_total ≈ 0.0
            error("Winch $(winch.name) has zero inertia_total. " *
                  "All winches must have non-zero inertia.")
        elseif winch.inertia_total < 0
            error("Winch $(winch.name) has negative inertia_total " *
                  "$(winch.inertia_total) kg⋅m². Inertia must be positive.")
        elseif winch.inertia_total < 1e-6
            @warn "Winch $(winch.name) has very small inertia_total " *
                  "$(winch.inertia_total) kg⋅m², may cause numerical issues"
        end

        # Check for non-positive drum_radius
        if winch.drum_radius <= 0
            error("Winch $(winch.name) has non-positive drum_radius ($(winch.drum_radius)). " *
                  "This will cause division by zero in torque conversions.")
        end

        # Check for non-positive gear_ratio
        if winch.gear_ratio <= 0
            error("Winch $(winch.name) has non-positive gear_ratio ($(winch.gear_ratio)). " *
                  "This will cause division by zero in speed/torque calculations.")
        end
    end

    # ==================== SEGMENT VALIDATIONS ==================== #
    for segment in segments
        # Wing structural segments don't use diameter (stiffness explicit, drag from VSM)
        wing_structural = all(points[i].type == WING for i in segment.point_idxs)
        if !wing_structural && !(0 < segment.diameter < 1)
            @warn "Segment $(segment.name) has unusual diameter " *
                  "$(segment.diameter) m (expected range: 0 to 1 m)"
        end

        # Rest length must be positive
        if segment.l0 <= 0
            error("Segment $(segment.name) has non-positive rest length " *
                  "l0 = $(segment.l0) m. This will cause division by zero.")
        end

        # Warn about zero or negative stiffness/damping (callable laws self-validate).
        if segment.unit_stiffness isa Real
            if segment.unit_stiffness ≈ 0.0
                @warn "Segment $(segment.name) has zero stiffness"
            elseif segment.unit_stiffness < 0
                @warn "Segment $(segment.name) has negative stiffness " *
                      "$(segment.unit_stiffness) N"
            end
        end

        if segment.unit_damping < 0
            @warn "Segment $(segment.name) has negative damping " *
                  "$(segment.unit_damping) N⋅s"
        end
    end

    # ==================== PULLEY VALIDATIONS ==================== #
    for pulley in pulleys
        if pulley.sum_len ≈ 0
            error("Pulley $(pulley.name) has zero total length constraint " *
                  "(sum_len = $(pulley.sum_len) m). " *
                  "Pulley constraints must have non-zero total length.")
        end
    end

    # ==================== TWIST_SURFACE VALIDATIONS ==================== #
    if length(twist_surfaces) > 0
        first_moment_frac = twist_surfaces[1].moment_frac
        for twist_surface in twist_surfaces
            if !(twist_surface.moment_frac ≈ first_moment_frac)
                error("TwistSurface $(twist_surface.name) has moment_frac = " *
                      "$(twist_surface.moment_frac), but all twist_surfaces must have the " *
                      "same moment_frac (first twist_surface has $(first_moment_frac))")
            end
        end
    end

    return nothing
end

# ==================== TETHER INIT LEN ==================== #

"""
    tether_ordered_point_idxs(tether, segments)

Point indices along the tether, ordered from `start_point_idx`
through each segment's far endpoint to the end point.
"""
function tether_ordered_point_idxs(tether, segments)
    idxs = Int64[tether.start_point_idx]
    for seg_idx in tether.segment_idxs
        push!(idxs, segments[seg_idx].point_idxs[2])
    end
    return idxs
end

"""
    tether_anchor_free(tether, boundary)

Return `(anchor_idx, free_idx)` for a root tether: the endpoint in
`boundary` (`STATIC`/winch points) is the anchor, the other is free.
Returns `(nothing, nothing)` if neither endpoint is on a boundary, or if
both are (a both-fixed tether cannot be placed; the caller warns and skips it).
"""
function tether_anchor_free(tether, boundary)
    start_on_boundary = tether.start_point_idx in boundary
    end_on_boundary = tether.end_point_idx in boundary
    if start_on_boundary && end_on_boundary
        return nothing, nothing
    elseif start_on_boundary
        return tether.start_point_idx, tether.end_point_idx
    elseif end_on_boundary
        return tether.end_point_idx, tether.start_point_idx
    end
    return nothing, nothing
end

"""
    rigid_point_siblings(points, wings)

Map each point index that rides a rigid body to the set of all points sharing
that body, so downstream traversal moves them as one unit. Covers `WING` points
of a `RIGID_DYNAMICS` wing (grouped by `wing_idx`) and `BODY_STATIC` points
(grouped by `body_idx`). These points carry no inter-point segments, so the set
captures their connectivity.
"""
function rigid_point_siblings(points, wings)
    siblings = Dict{Int64, Set{Int64}}()
    for wing in wings
        wing.dynamics_type == RIGID_DYNAMICS || continue
        members = Set{Int64}(point.idx for point in points
            if point.type == WING && point.wing_idx == wing.idx)
        for member in members
            siblings[member] = members
        end
    end
    body_idxs = unique(point.body_idx for point in points
        if point.type == BODY_STATIC)
    for body_idx in body_idxs
        members = Set{Int64}(point.idx for point in points
            if point.type == BODY_STATIC && point.body_idx == body_idx)
        for member in members
            siblings[member] = members
        end
    end
    return siblings
end

"""
    tether_downstream_idxs(tether, segments, boundary, from_idx,
                           anchor_idx, rigid_siblings)

Breadth-first set of point indices reachable from `from_idx` (the
tether's free end) through segments outside this tether and through
`rigid_siblings`, stopping at boundary points. These are the points
that must translate with the free end when the tether is repositioned.
Errors if traversal reaches `anchor_idx` (a loop back to the anchor).
"""
function tether_downstream_idxs(tether, segments, boundary,
                                from_idx, anchor_idx, rigid_siblings)
    own = Set{Int64}(tether_ordered_point_idxs(tether, segments))
    tether_segment_set = Set(tether.segment_idxs)
    visited = copy(own)
    downstream = Set{Int64}()
    queue = Int64[from_idx]
    while !isempty(queue)
        current_idx = pop!(queue)
        neighbors = Int64[]
        for seg in segments
            seg.idx in tether_segment_set && continue
            point_idx1, point_idx2 = seg.point_idxs
            if point_idx1 == current_idx
                push!(neighbors, point_idx2)
            elseif point_idx2 == current_idx
                push!(neighbors, point_idx1)
            end
        end
        if haskey(rigid_siblings, current_idx)
            for sibling in rigid_siblings[current_idx]
                sibling == current_idx || push!(neighbors, sibling)
            end
        end
        for neighbor_idx in neighbors
            if neighbor_idx == anchor_idx
                error("Tether $(tether.name): downstream structure " *
                      "connects back to the anchor point. Cannot " *
                      "apply tether length scaling.")
            end
            neighbor_idx in visited && continue
            push!(visited, neighbor_idx)
            neighbor_idx in boundary && continue
            push!(downstream, neighbor_idx)
            push!(queue, neighbor_idx)
        end
    end
    return downstream
end

"""
    twist_surface_tethers_by_overlap(specified, reach)

Cluster the `specified` tethers with a union-find over `reach`
(point indices each tether touches): tethers whose reaches intersect
share structure and land in the same cluster. Returns a vector of
tether vectors, one per cluster.
"""
function twist_surface_tethers_by_overlap(specified, reach)
    n = length(specified)
    parent = collect(1:n)
    function find_root(i)
        while parent[i] != i
            parent[i] = parent[parent[i]]
            i = parent[i]
        end
        return i
    end
    for i in 1:n
        for j in i+1:n
            isempty(intersect(reach[specified[i].idx],
                              reach[specified[j].idx])) && continue
            root_i = find_root(i)
            root_j = find_root(j)
            root_i == root_j && continue
            parent[root_i] = root_j
        end
    end
    twist_surfaces = Dict{Int64, Vector{Tether}}()
    for i in 1:n
        push!(get!(() -> Tether[], twist_surfaces, find_root(i)), specified[i])
    end
    return collect(values(twist_surfaces))
end

"""
    tether_unit_stiffness(tether, segments)

Return the common per-unit-length stiffness `[N]` of the tether's
segments. Errors if the segments are not uniform, since the spring
inversion in `apply_tether_init_forces!` assumes a single stiffnesys_state.
"""
function tether_unit_stiffness(tether, segments)
    any(!(segments[i].unit_stiffness isa Real) for i in tether.segment_idxs) &&
        error("Tether $(tether.name): init_tether_force needs constant " *
              "unit_stiffness; a segment has a nonlinear force law. " *
              "Use init_stretch_frac instead.")
    stiffness_values = SimFloat[segments[seg_idx].unit_stiffness
                                for seg_idx in tether.segment_idxs]
    stiffness = first(stiffness_values)
    all(≈(stiffness), stiffness_values) ||
        error("Tether $(tether.name): requires uniform unit_stiffness " *
              "across its segments, got $stiffness_values")
    return stiffness
end

"""
    apply_cluster_init_stretched_len!(cluster, points, segments,
                                      downstream, boundary; prn=true)

Reposition one cluster of root tethers so each sits at its
`init_stretched_len` standoff. Each tether contributes the
displacement that would move its free end onto the target length
along the anchor→free direction; the free end and everything
downstream of it are translated by the mean of those displacements,
then interior points are redistributed proportionally along each
tether. For a multi-tether cluster, logs an `@info` when `prn`.
"""
function apply_cluster_init_stretched_len!(
    cluster, points, segments, bodies, downstream, boundary; prn=true)
    snaps = map(cluster) do tether
        anchor_idx, free_idx = tether_anchor_free(tether, boundary)
        anchor_pos = copy(points[anchor_idx].pos_w)
        free_pos = copy(points[free_idx].pos_w)
        ordered = tether_ordered_point_idxs(tether, segments)
        seg_lens = SimFloat[segment_world_length(segments[seg_idx], points)
                            for seg_idx in tether.segment_idxs]
        if ordered[1] != anchor_idx
            reverse!(ordered)
            reverse!(seg_lens)
        end
        path_len = sum(seg_lens)
        path_len > 0 || error("Tether $(tether.name): current length is " *
            "zero, cannot scale to its stretched length")
        (; tether, free_idx, anchor_pos, free_pos, ordered, seg_lens, path_len)
    end

    deltas = [(snap.tether.init_stretched_len::SimFloat / snap.path_len - 1) .*
              (snap.free_pos .- snap.anchor_pos) for snap in snaps]
    delta = sum(deltas) ./ length(deltas)

    if length(cluster) > 1 && prn
        names = join((string(snap.tether.name) for snap in snaps), ", ")
        @info "Tethers ($names) feed one structure; placing it to the " *
              "mean stretched length and direction of all."
    end
    norm(delta) ≈ 0 && return

    moved = Set{Int64}()
    for snap in snaps
        for idx in downstream[snap.tether.idx]
            idx in moved && continue
            push!(moved, idx)
            points[idx].pos_w .+= delta
        end
        if !(snap.free_idx in moved)
            push!(moved, snap.free_idx)
            points[snap.free_idx].pos_w .+= delta
        end
    end

    # Move the body, not its points: the pos~anchor constraint would snap them back.
    # WING points carry their body association in wing_idx (body_idx is 0);
    # BODY_STATIC points use body_idx.
    moved_bodies = Set{Int64}()
    for idx in moved
        point = points[idx]
        if point.type == WING && point.wing_idx != 0
            push!(moved_bodies, point.wing_idx)
        elseif point.body_idx != 0
            push!(moved_bodies, point.body_idx)
        end
    end
    for body_idx in moved_bodies
        bodies[body_idx].pos_w .+= delta
        bodies[body_idx].com_w .+= delta
    end

    for snap in snaps
        length(snap.ordered) <= 2 && continue
        line = (snap.free_pos .+ delta) .- snap.anchor_pos
        cum = 0.0
        for k in 2:length(snap.ordered)-1
            cum += snap.seg_lens[k-1]
            points[snap.ordered[k]].pos_w .=
                snap.anchor_pos .+ (cum / snap.path_len) .* line
        end
    end
end

"""
    apply_tether_init_stretched_lens!(sys_struct::SystemStructure; prn=true)

Scale `pos_w` so each tether with an explicit `init_stretched_len` sits at
that standoff. Call after `copy_cad_to_world!`. Rest length is derived
separately by `apply_tether_init_forces!`.

Only tethers with one endpoint on a boundary (`STATIC` or winch point) are
placed; that endpoint is the fixed anchor (start or end). Scaling runs from
the anchor toward the free end, translating everything downstream of it.
A tether with neither endpoint anchored is an error. Roots feeding one
structure form a cluster, placed by their mean displacement (length and
direction).

Errors if a downstream segment connects back to the anchor.
"""
function apply_tether_init_stretched_lens!(sys_struct::SystemStructure;
                                           prn=true)
    (; points, segments, tethers, winches, wings) = sys_struct

    specified = [tether for tether in tethers
                 if !isnothing(tether.init_stretched_len)]
    isempty(specified) && return

    rigid_siblings = rigid_point_siblings(points, wings)

    # Boundary = externally world-fixed points: STATIC, winch, and BODY_STATIC on a STATIC body.
    bodies = sys_struct.bodies
    boundary = Set{Int64}(w.winch_point_idx for w in winches)
    for point in points
        point.type == STATIC && push!(boundary, point.idx)
        point.type == BODY_STATIC && bodies[point.body_idx].type == STATIC &&
            push!(boundary, point.idx)
    end

    # Both-fixed tethers (both endpoints on a boundary) are warned and skipped.
    both_fixed(tether) = tether.start_point_idx in boundary &&
                         tether.end_point_idx in boundary
    placeable = filter(!both_fixed, specified)
    for tether in specified
        both_fixed(tether) && @warn "Tether $(tether.name): both endpoints " *
            "are fixed; skipping its length placement."
    end
    isempty(placeable) && return
    specified = placeable

    anchor_free = Dict(tether.idx => tether_anchor_free(tether, boundary)
                       for tether in specified)
    non_root = [tether for tether in specified
                if isnothing(anchor_free[tether.idx][1])]
    if !isempty(non_root)
        names = join((string(tether.name) for tether in non_root), ", ")
        error("tether length is only supported on tethers anchored at " *
              "a STATIC or winch point. Tether(s) ($names) have neither " *
              "endpoint anchored; their position rides the root tether.")
    end

    downstream = Dict(tether.idx => tether_downstream_idxs(
                          tether, segments, boundary, anchor_free[tether.idx][2],
                          anchor_free[tether.idx][1], rigid_siblings)
                      for tether in specified)
    reach = Dict(tether.idx => union(
        setdiff(Set{Int64}(tether_ordered_point_idxs(tether, segments)),
                boundary),
        downstream[tether.idx]) for tether in specified)

    for cluster in twist_surface_tethers_by_overlap(specified, reach)
        apply_cluster_init_stretched_len!(cluster, points, segments,
                                          sys_struct.bodies,
                                          downstream, boundary; prn)
    end
end

"""
    init_unstretched_len(tether, segments) -> SimFloat

Derived initial unstretched (rest) length from the tether's current
(placed) stretched length `stretched = Σ segment lengths`:
- `init_stretch_frac` set: `len = stretch_frac · stretched`
  (`< 1` pre-stretch, `1` neutral, `> 1` slack).
- otherwise from `init_tether_force` (default 0):
  `len = stretched · (1 − force / unit_stiffness)` (zero-velocity,
  tension branch). Force 0 gives `len = stretched` (no tension).

Errors if both `init_stretch_frac` and `init_tether_force` are set,
if `stretch_frac ≤ 0`, if `force < 0`, if `force ≥
unit_stiffness`, or if the segments have non-uniform `unit_stiffness`.
"""
function init_unstretched_len(tether, segments)
    stretched = sum(segments[seg_idx].len
                    for seg_idx in tether.segment_idxs)
    frac = tether.init_stretch_frac
    force = tether.init_tether_force
    if !isnothing(frac) && !isnothing(force)
        error("Tether $(tether.name): set only one of " *
              "init_stretch_frac and init_tether_force")
    end
    if !isnothing(frac)
        frac > 0 || error("Tether $(tether.name): " *
            "init_stretch_frac $frac must be positive")
        return stretched * frac
    end
    force = something(force, 0.0)
    force >= 0 || error("Tether $(tether.name): " *
        "init_tether_force $force N is negative; " *
        "compression is not supported")
    force == 0 && return stretched
    stiffness = tether_unit_stiffness(tether, segments)
    force < stiffness || error("Tether $(tether.name): " *
        "init_tether_force $force N ≥ unit_stiffness $stiffness N; " *
        "no positive rest length achieves this force")
    return stretched * (1 - force / stiffness)
end

"""
    apply_tether_init_forces!(sys_struct::SystemStructure)

Set every tether's `len` to its [`init_unstretched_len`](@ref).
Must be called after segment world lengths are current.
"""
function apply_tether_init_forces!(sys_struct::SystemStructure)
    (; segments, tethers) = sys_struct
    for tether in tethers
        isempty(tether.segment_idxs) && continue
        tether.len = init_unstretched_len(tether, segments)
    end
end

# ==================== REINIT! FOR SYSTEM STRUCTURE ==================== #

"""
    reinit!(sys_struct::SystemStructure, set::Settings; kwargs...)

Re-initialize a `SystemStructure` from a `Settings` object.

This function resets various component states (e.g., winch lengths, twist_surface twists,
pulley positions) to their initial values as defined in the `Settings` object. It
is typically called before starting a new simulation run.

Pulley lengths are initialized proportionally based on current segment lengths:
`pulley.len = segment1.len / (segment1.len+segment2.len) * pulley.sum_len`

# Keyword Arguments
- `ignore_l0::Bool=false`: If true, recalculate segment rest lengths from current positions
- `remake_vsm::Bool=false`: If true, recreate VSM wing, aerodynamics, and solver from settings.
  This is useful after modifying `aero_geometry.yaml` or other VSM-related configuration files.
  For PARTICLE_DYNAMICS wings, also rebuilds the `point_to_vsm_point` mapping.
- `apply_transforms::Bool=true`: If false, skip applying spatial transforms
  (translate, rotate, heading) during reinitialization.
- `apply_tether_lengths::Bool=true`: If false, skip scaling point positions
  to match `tether.init_stretched_len`.
- `prn::Bool=true`: If true, print info messages (e.g. when several root
  tethers are placed to their mean stretched length).
"""
function reinit!(sys_struct::SystemStructure, set::Settings;
                 ignore_l0::Bool=false, remake_vsm::Bool=false,
                 reset_vel::Bool=true, apply_transforms::Bool=true,
                 apply_tether_lengths::Bool=true, prn::Bool=true)
    (; points, twist_surfaces, segments, pulleys, tethers, winches, wings, transforms) = sys_struct

    for winch in winches
        winch.vel = winch.init_vel
    end

    # Reset body pose to CAD (idempotent placement); ODE state derived below.
    for rigid_body in sys_struct.bodies
        rigid_body.pos_w .= rigid_body.pos_cad
        rigid_body.Q_b_to_w .= rotation_matrix_to_quaternion(rigid_body.R_b_to_c)
        if reset_vel
            rigid_body.vel_w .= 0.0
            rigid_body.ω_b .= 0.0
        end
        init_rigid_body!(rigid_body)
    end

    for twist_surface in twist_surfaces
        twist_surface.type == STATIC && continue
        twist_surface.twist = 0.0
        twist_surface.twist_ω = 0.0
    end

    # Transforms are not updated from Settings; YAML structure geometry has priority.

    # Step 1: copy CAD geometry to world frame
    copy_cad_to_world!(points, sys_struct.bodies; update_vel=reset_vel)

    # Step 2: apply stretched lengths (scales pos_w)
    if apply_tether_lengths
        apply_tether_init_stretched_lens!(sys_struct; prn)
    end

    # Step 3: compute segment lengths from pos_w
    for segment in segments
        len = segment_world_length(segment, points)
        (segment.l0 ≈ 0) && (segment.l0 = len)
        segment.len = len
    end

    apply_tether_init_forces!(sys_struct)

    for tether in tethers
        n = length(tether.segment_idxs)
        n == 0 && continue
        l0 = tether.len / n
        for seg_idx in tether.segment_idxs
            segments[seg_idx].l0 = l0
        end
    end

    for pulley in pulleys
        segment1, segment2 = segments[pulley.segment_idxs[1]],
                             segments[pulley.segment_idxs[2]]
        pulley.sum_len = segment1.l0 + segment2.l0

        # Proportional to current segment lengths (accurate for asymmetric bridles).
        pulley.len = segment1.len / (segment1.len+segment2.len) *
                     pulley.sum_len

        pulley.vel = 0.0
    end

    # Step 5: apply transforms (translate/rotate/heading); pos_w already initialized.
    if apply_transforms
        reinit!(transforms, sys_struct; update_vel=reset_vel)
    end

    # Recreate each wing's aero engine from settings (no-op for engine-less modes).
    if remake_vsm
        for wing in wings
            remake_aero!(wing.aero, wing, set, sys_struct.vsm_set,
                         points, twist_surfaces)
        end
    end

    # Compute per-wing wind from settings
    wind_vec_gnd = set.wind_vec

    wind_factor = WindFactor(sys_struct.am, sys_struct.set.profile_law)
    for wing in wings
        # Calculate wind at wing position using atmospheric model
        wing.v_wind .= wind_factor(wing.pos_w[3]) * wind_vec_gnd

        R_b_to_w = wing.R_b_to_w::Matrix{SimFloat}
        if wing.dynamics_type == PARTICLE_DYNAMICS
            va_wing_w = wing.v_wind - wing.vel_w + wing.wind_disturb
            wing.va_b .= R_b_to_w' * va_wing_w
        else
            # Initialize the aero operating point from the initial wind
            init_aero_state!(wing.aero, wing, R_b_to_w' * wind_vec_gnd)
        end
    end

    # validate_sys_struct() runs later: total_mass needs the live integrator.

    # Recalculate segment rest lengths from current positions if requested
    if ignore_l0
        for segment in segments
            point1 = points[segment.point_idxs[1]]
            point2 = points[segment.point_idxs[2]]
            segment.l0 = norm(point2.pos_w - point1.pos_w)
        end
    end

    # Joint rest geometry, from the final placed body poses (as-placed = unstrained).
    init_joint_rest!.(sys_struct.elastic_joints, Ref(sys_struct.bodies))
    init_joint_rest!.(sys_struct.timoshenko_joints, Ref(sys_struct.bodies))

    return nothing
end

# ==================== STATE COPY ==================== #

"""
    copy!(sys1::SystemStructure, sys2::SystemStructure)

Copy the dynamic state from one `SystemStructure` (`sys1`) to another (`sys2`).

This function is designed to transfer the state (positions, velocities, etc.) between
two system models, which can have different levels of fidelity. For example, it can
copy the state from a detailed multi-segment tether model (`sys1`) to a simplified
single-segment model (`sys2`).

The function handles several cases:
- If `sys1` and `sys2` have the same structure, it performs a direct copy of all point states.
- If `sys2` is a simplified (1-segment per tether) version of `sys1`, it copies the
  positions and velocities of the tether endpoints.
- It also copies the state of wings, twist_surfaces, winches, and pulleys where applicable.
"""
function copy!(sys1::SystemStructure, sys2::SystemStructure)

    # copy point pos and vel
    if length(sys1.points) > 0
        if length(sys1.points) == length(sys2.points)
            for (point1, point2) in zip(sys1.points, sys2.points)
                point2.pos_w .= point1.pos_w
                point2.vel_w .= point1.vel_w
                point2.disturb .= point1.disturb
            end
        # if different number of points, copy only the tether points
        elseif length(sys1.tethers) > 1 && length(sys1.tethers) == length(sys2.tethers)
            for (tether1, tether2) in zip(sys1.tethers, sys2.tethers)
                if length(tether1.segment_idxs) == length(tether2.segment_idxs)
                    # copy the points of the segments of the tethers
                    for (segment_idx1, segment_idx2) in zip(tether1.segment_idxs, tether2.segment_idxs)
                        point_idxs1 = sys1.segments[segment_idx1].point_idxs
                        point_idxs2 = sys2.segments[segment_idx2].point_idxs
                        for (point_idx1, point_idx2) in zip(point_idxs1, point_idxs2)
                            sys2.points[point_idx2].pos_w .= sys1.points[point_idx1].pos_w
                            sys2.points[point_idx2].vel_w .= sys1.points[point_idx1].vel_w
                            sys2.points[point_idx2].disturb .= sys1.points[point_idx1].disturb
                        end
                    end
                elseif length(tether2.segment_idxs) == 1
                    # copy the first and last point of the tether
                    point_idxs1 = [sys1.segments[tether1.segment_idxs[1]].point_idxs[1],
                                   sys1.segments[tether1.segment_idxs[end]].point_idxs[2]]
                    point_idxs2 = sys2.segments[tether2.segment_idxs[1]].point_idxs
                    sys2.points[point_idxs2[1]].pos_w .= sys1.points[point_idxs1[1]].pos_w
                    sys2.points[point_idxs2[2]].pos_w .= sys1.points[point_idxs1[2]].pos_w
                    sys2.points[point_idxs2[1]].vel_w .= sys1.points[point_idxs1[1]].vel_w
                    sys2.points[point_idxs2[2]].vel_w .= sys1.points[point_idxs1[2]].vel_w
                    sys2.points[point_idxs2[1]].disturb .= sys1.points[point_idxs1[1]].disturb
                    sys2.points[point_idxs2[2]].disturb .= sys1.points[point_idxs1[2]].disturb
                end
            end
        end
    end

    # copy twist and twist_ω of twist_surfaces
    if length(sys1.twist_surfaces) > 0 && length(sys1.twist_surfaces) == length(sys2.twist_surfaces)
        for (twist_surface1, twist_surface2) in zip(sys1.twist_surfaces, sys2.twist_surfaces)
            twist_surface2.twist = twist_surface1.twist
            twist_surface2.twist_ω = twist_surface1.twist_ω
        end
    end

    # copy tether lengths
    if length(sys1.tethers) > 0 &&
       length(sys1.tethers) == length(sys2.tethers)
        for (tether2, tether1) in
                zip(sys2.tethers, sys1.tethers)
            tether2.len = tether1.len
        end
    end

    # copy winch velocities
    if length(sys1.winches) > 0 &&
       length(sys1.winches) == length(sys2.winches)
        for (winch2, winch1) in
                zip(sys2.winches, sys1.winches)
            winch2.vel = winch1.vel
        end
    end

    # copy pulley lengths and velocities
    if length(sys1.pulleys) > 0 && length(sys1.pulleys) == length(sys2.pulleys)
        for (pulley1, pulley2) in zip(sys1.pulleys, sys2.pulleys)
            pulley2.len = pulley1.len
            pulley2.vel = pulley1.vel
        end
    end

    # copy wing positions and velocities
    if length(sys1.wings) > 0 && length(sys1.wings) == length(sys2.wings)
        for (wing1, wing2) in zip(sys1.wings, sys2.wings)
            wing2.pos_w .= wing1.pos_w
            wing2.vel_w .= wing1.vel_w
            wing2.ω_b .= wing1.ω_b
            wing2.Q_b_to_w .= wing1.Q_b_to_w
        end
    end
end

# ==================== SYSSTATE INTEROP ==================== #

"""
    update_from_sysstate!(sys::SystemStructure, sys_state::SysState)

Update the dynamic state of a `SystemStructure` from a `SysState` snapshot.

This function copies the state variables that are present in `SysState` (such as point
positions, wing orientations, winch lengths, and twist angles) into an existing `SystemStructure`.
Fields that cannot be populated from `SysState` (such as aerodynamic forces, moments, and
segment forces) are set to `NaN` to prevent them from being plotted.

This is useful for visualizing a `SysLog` by extracting individual `SysState` snapshots
and applying them to a `SystemStructure` for plotting with the Makie extension.

# Arguments
- `sys::SystemStructure`: The system structure to update (must already exist with correct topology).
- `sys_state::SysState`: The state snapshot to copy from.

# Example
```julia
# Load a system log
sim_log = load_log(...)

# Create a SystemStructure with the same topology
sys = SystemStructure(se(), "ram")

# Update from a specific time step
update_from_sysstate!(sys, sim_log.syslog[100])

# Plot the system at that time step
plot(sys)
```

# Notes
- The `SystemStructure` must have been created with the same model configuration as the
  simulation that generated the `SysLog`.
- Aerodynamic and force fields are set to `NaN` and will not be plotted.
- The number of points in `sys` must match the parametric type `P` of `SysState{P}`.
"""
function update_from_sysstate!(sys::SystemStructure, sys_state::SysState{P}) where P
    (; points, twist_surfaces, tethers, winches, wings, bodies) = sys

    # Position slot layout (points, panel corners, wing origins, body origins).
    slots = position_slots(sys)
    n_points = length(points)
    n_panel_corners = count_aero_log_points(wings)
    n_wings = length(wings)
    total_without_wings = n_points + n_panel_corners
    has_wing_slots = P == slots.total

    if !has_wing_slots && P != total_without_wings
        error("SystemStructure expects $(slots.total) points " *
              "($n_points regular + $n_panel_corners corners + " *
              "$n_wings wings + $(length(bodies)) bodies) or " *
              "$total_without_wings without wing slots, but SysState has $P points")
    end

    # Update point positions (X, Y, Z from SysState)
    for point in points
        point.pos_w[1] = sys_state.X[point.idx]
        point.pos_w[2] = sys_state.Y[point.idx]
        point.pos_w[3] = sys_state.Z[point.idx]
        # Set velocity to zero (not available in basic SysState)
        point.vel_w .= 0.0
        # Set forces to NaN (not available in SysState)
        point.force .= NaN
    end

    # Update wing state if wings exist
    if length(wings) > 0 && length(wings) == 1  # Currently only support single-wing systems
        wing = wings[1]

        # Copy orientation quaternion
        wing.Q_b_to_w .= sys_state.orient

        # Copy spherical coordinates
        wing.elevation = Float64(sys_state.elevation)
        wing.azimuth = Float64(sys_state.azimuth)
        wing.heading = Float64(sys_state.heading)

        if has_wing_slots
            wing_slot = slots.wings[wing.idx]
            wing.pos_w[1] = sys_state.X[wing_slot]
            wing.pos_w[2] = sys_state.Y[wing_slot]
            wing.pos_w[3] = sys_state.Z[wing_slot]
        else
            wing.pos_w .= [mean(sys_state.X), mean(sys_state.Y), mean(sys_state.Z)]
        end

        # Copy velocity if available in vel_kite
        wing.vel_w .= sys_state.vel_kite

        # Set angular velocity to NaN (turn_rates in SysState, but need conversion)
        wing.ω_b .= sys_state.turn_rates

        # Set aerodynamic quantities to NaN (to prevent plotting)
        wing.aero_force_b .= NaN
        wing.aero_moment_b .= NaN
        wing.tether_force .= NaN
        wing.tether_moment .= NaN
        wing.va_b .= NaN
        wing.v_wind .= sys_state.v_wind_kite
        wing.aoa = Float64(sys_state.AoA)
        wing.course = Float64(sys_state.course)
        wing.acc_w .= 0.0
        wing.turn_rate .= sys_state.turn_rates
        wing.turn_acc .= 0.0
    end

    # Standalone bodies: origins after the wing slots, orientations after the wings.
    if has_wing_slots
        for rigid_body in bodies
            slot = slots.bodies[rigid_body.idx]
            rigid_body.pos_w[1] = sys_state.X[slot]
            rigid_body.pos_w[2] = sys_state.Y[slot]
            rigid_body.pos_w[3] = sys_state.Z[slot]
            rigid_body.Q_b_to_w .= sys_state.orients[n_wings + rigid_body.idx]
        end
    end

    # Update twist_surface twist angles
    n_twist_surfaces = min(length(twist_surfaces), 4)  # SysState stores up to 4 twist angles
    for i in 1:n_twist_surfaces
        if i <= length(twist_surfaces)
            twist_surfaces[i].twist = Float64(sys_state.twist_angles[i])
            twist_surfaces[i].twist_ω = 0.0  # Not available in SysState
            # Set forces/moments to NaN
            twist_surfaces[i].tether_force = NaN
            twist_surfaces[i].tether_moment = NaN
            twist_surfaces[i].aero_moment = NaN
        end
    end

    for wing in wings
        restore_aero_twist!(wing.aero, wing, twist_surfaces)
    end

    # Update tether lengths from SysState (per-tether)
    for (tether_idx, tether) in enumerate(tethers)
        tether_idx > 4 && break
        tether.len = Float64(sys_state.l_tether[tether_idx])
    end

    # Update winch state from SysState (per-winch)
    n_winches = min(length(winches), 4)
    for i in 1:n_winches
        winches[i].force .= NaN
        winches[i].friction = NaN
        winches[i].acc = 0.0
        winches[i].vel = Float64(sys_state.v_reelout[i])
        winches[i].set_value = Float64(sys_state.set_torque[i])
    end

    corner_idx = n_points
    for wing in wings
        corner_idx = read_aero_log_points!(wing.aero, wing, sys,
                                           sys_state, corner_idx)
    end

    # Update global wind vector (only if wind_vec mode is active)
    if sys.set.use_wind_vec
        sys.set.wind_vec = MVec3(sys_state.v_wind_gnd)
    end

    # Segment lengths/forces come from symbolic getters, not from the SysState.

    return nothing
end

# ==================== DAMPING SETTERS ==================== #

"""
    set_world_frame_damping(sys::SystemStructure, damping, point_idxs)

Set the world frame damping coefficient for specified points in the system structure.

World frame damping applies a velocity-dependent drag force in the global
reference frame: ``\\mathbf{F}_{damp} = -c_{damp} \\odot \\mathbf{v}``, where
``c_{damp}`` is the damping vector and ``\\odot`` is element-wise multiplication.

# Arguments
- `sys::SystemStructure`: The system structure to modify.
- `damping::Union{Real, AbstractVector}`: Damping coefficient(s) [N·s/m].
  Scalar applies same value to all 3 axes. Vector must have 3 elements for [x,y,z] damping.
- `point_idxs`: Indices of points to apply damping to.

# Returns
- `nothing`
"""
function set_world_frame_damping(sys::SystemStructure, damping::Union{Real, AbstractVector},
                                 point_idxs)
    damp_vec = damping isa Real ? SVector{3,SimFloat}(damping, damping, damping) : SVector{3,SimFloat}(damping)
    @assert length(damp_vec) == 3 "Damping must be scalar or 3-element vector"
    for idx in point_idxs
        sys.points[idx].world_frame_damping = damp_vec
    end
    return nothing
end

"""
    set_world_frame_damping(sys::SystemStructure, damping)

Set the world frame damping coefficient for all points in the system structure.

World frame damping applies a velocity-dependent drag force in the global
reference frame: ``\\mathbf{F}_{damp} = -c_{damp} \\odot \\mathbf{v}``, where
``c_{damp}`` is the damping vector and ``\\odot`` is element-wise multiplication.

# Arguments
- `sys::SystemStructure`: The system structure to modify.
- `damping::Union{Real, AbstractVector}`: Damping coefficient(s) [N·s/m].
  Scalar applies same value to all 3 axes. Vector must have 3 elements for [x,y,z] damping.

# Returns
- `nothing`
"""
function set_world_frame_damping(sys::SystemStructure, damping::Union{Real, AbstractVector})
    set_world_frame_damping(sys, damping, eachindex(sys.points))
end

"""
    set_body_frame_damping(sys::SystemStructure, damping, point_idxs)

Set the body frame damping coefficient for specified points in the system structure.

# Arguments
- `sys::SystemStructure`: The system structure to modify.
- `damping::Union{Real, AbstractVector}`: Damping coefficient(s) [N·s/m].
  Scalar applies same value to all 3 axes. Vector must have 3 elements for [x,y,z] damping.
- `point_idxs`: Indices of points to apply damping to.

# Returns
- `nothing`
"""
function set_body_frame_damping(sys::SystemStructure, damping::Union{Real, AbstractVector},
                                point_idxs)
    damp_vec = damping isa Real ? SVector{3,SimFloat}(damping, damping, damping) : SVector{3,SimFloat}(damping)
    @assert length(damp_vec) == 3 "Damping must be scalar or 3-element vector"
    for idx in point_idxs
        sys.points[idx].body_frame_damping = damp_vec
    end
    return nothing
end

"""
    set_body_frame_damping(sys::SystemStructure, damping)

Set the body frame damping coefficient for all points in the system structure.

# Arguments
- `sys::SystemStructure`: The system structure to modify.
- `damping::Union{Real, AbstractVector}`: Damping coefficient(s) [N·s/m].
  Scalar applies same value to all 3 axes. Vector must have 3 elements for [x,y,z] damping.

# Returns
- `nothing`
"""
function set_body_frame_damping(sys::SystemStructure, damping::Union{Real, AbstractVector})
    set_body_frame_damping(sys, damping, eachindex(sys.points))
end

# ==================== SEGMENT STATISTICS ==================== #

"""
    segment_stretch_stats(sys::SystemStructure)

Calculate segment stretch statistics for segments in tension.

Returns the maximum and mean relative stretch of segments where len > l0,
along with the index of the segment with maximum stretch.
Relative stretch is defined as (current_length - l0) / l0.
Only segments in tension (stretched) are included in the statistics.

For pulley segments, the combined length of both segments is used against
the pulley's sum_l0, since the pulley constraint distributes length between them.

# Arguments
- `sys::SystemStructure`: System structure with current segment states

# Returns
- `(max_stretch, mean_stretch, max_idx)`: Tuple of maximum stretch, mean stretch,
  and index of the segment with maximum stretch (or first pulley segment index)
"""
function segment_stretch_stats(sys::SystemStructure)
    if isempty(sys.segments)
        return (0.0, 0.0, 0)
    end

    # Build set of segment indices that belong to pulleys
    pulley_seg_idxs = Set{Int64}()
    for pulley in sys.pulleys
        push!(pulley_seg_idxs, pulley.segment_idxs[1])
        push!(pulley_seg_idxs, pulley.segment_idxs[2])
    end

    stretch_data = Tuple{Int64, Float64}[]

    # Add pulley stretches (combined length of both segments)
    for pulley in sys.pulleys
        seg1 = sys.segments[pulley.segment_idxs[1]]
        seg2 = sys.segments[pulley.segment_idxs[2]]
        combined_len = seg1.len + seg2.len
        sum_l0 = pulley.sum_len  # sum_len stores seg1.l0 + seg2.l0
        if combined_len > sum_l0
            push!(stretch_data, (pulley.segment_idxs[1],
                                 (combined_len - sum_l0) / sum_l0))
        end
    end

    # Add non-pulley segment stretches
    for seg in sys.segments
        if seg.idx ∉ pulley_seg_idxs && seg.len > seg.l0
            push!(stretch_data, (seg.idx, (seg.len - seg.l0) / seg.l0))
        end
    end

    if isempty(stretch_data)
        return (0.0, 0.0, 0)
    end

    stretches = [entry[2] for entry in stretch_data]
    max_stretch = maximum(stretches)
    mean_stretch = sum(stretches) / length(stretches)
    max_idx = stretch_data[argmax(stretches)][1]

    return (max_stretch, mean_stretch, max_idx)
end
