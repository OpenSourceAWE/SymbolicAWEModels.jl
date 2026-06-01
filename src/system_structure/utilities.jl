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
    p1 = points[segment.point_idxs[1]]
    p2 = points[segment.point_idxs[2]]
    return norm(p1.pos_cad - p2.pos_cad)
end

"""
    segment_world_length(segment::Segment, points)

Compute segment length from endpoint `pos_w` positions.
"""
function segment_world_length(segment::Segment, points)
    p1 = points[segment.point_idxs[1]]
    p2 = points[segment.point_idxs[2]]
    return norm(p1.pos_w - p2.pos_w)
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
- Empty group list for RIGID_DYNAMICS wings (warning)
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

## Group Validations
- Inconsistent moment_frac across groups (error)
"""
function validate_sys_struct(sys_struct::SystemStructure)
    (; points, groups, segments, pulleys, wings, winches) = sys_struct

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

        # Check total_mass for division by zero (DYNAMIC points only)
        # NOTE: Check mass before NaN position - NaN pos is often caused by zero mass
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
        # Check for non-positive mass (all wing types)
        # NOTE: Check mass/inertia before NaN position - NaN pos is often caused by zero mass
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

            # Warn if RIGID_DYNAMICS wing has no groups
            # (expected for AERO_NONE which skips auto-group creation)
            if isempty(wing.group_idxs) &&
               wing.aero_mode != AERO_NONE
                @warn "Wing $(wing.name) (RIGID_DYNAMICS)" *
                    " has no groups"
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
        # Diameter should be in valid range (warn only, not critical)
        if !(0 < segment.diameter < 1)
            @warn "Segment $(segment.name) has unusual diameter " *
                  "$(segment.diameter) m (expected range: 0 to 1 m)"
        end

        # Rest length must be positive
        if segment.l0 <= 0
            error("Segment $(segment.name) has non-positive rest length " *
                  "l0 = $(segment.l0) m. This will cause division by zero.")
        end

        # Warn about zero or negative stiffness/damping
        if segment.unit_stiffness ≈ 0.0
            @warn "Segment $(segment.name) has zero stiffness"
        elseif segment.unit_stiffness < 0
            @warn "Segment $(segment.name) has negative stiffness " *
                  "$(segment.unit_stiffness) N"
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

    # ==================== GROUP VALIDATIONS ==================== #
    if length(groups) > 0
        first_moment_frac = groups[1].moment_frac
        for group in groups
            if !(group.moment_frac ≈ first_moment_frac)
                error("Group $(group.name) has moment_frac = " *
                      "$(group.moment_frac), but all groups must have the " *
                      "same moment_frac (first group has $(first_moment_frac))")
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
Returns `(nothing, nothing)` if neither endpoint is on a boundary;
errors if both are.
"""
function tether_anchor_free(tether, boundary)
    s_in = tether.start_point_idx in boundary
    e_in = tether.end_point_idx in boundary
    if s_in && e_in
        error("Tether $(tether.name): both endpoints are " *
              "ground-fixed; cannot place it to a length.")
    elseif s_in
        return tether.start_point_idx, tether.end_point_idx
    elseif e_in
        return tether.end_point_idx, tether.start_point_idx
    end
    return nothing, nothing
end

"""
    rigid_point_siblings(points, wings)

Map each `WING`-type point index of a `RIGID_DYNAMICS` wing to the
set of all such points sharing that wing. These points move as one
rigid body without inter-point segments, so the set captures their
connectivity for downstream traversal.
"""
function rigid_point_siblings(points, wings)
    siblings = Dict{Int64, Set{Int64}}()
    for wing in wings
        wing.dynamics_type == RIGID_DYNAMICS || continue
        members = Set{Int64}(p.idx for p in points
            if p.type == WING && p.wing_idx == wing.idx)
        for m in members
            siblings[m] = members
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
            p1, p2 = seg.point_idxs
            if p1 == current_idx
                push!(neighbors, p2)
            elseif p2 == current_idx
                push!(neighbors, p1)
            end
        end
        if haskey(rigid_siblings, current_idx)
            for sib in rigid_siblings[current_idx]
                sib == current_idx || push!(neighbors, sib)
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
    group_tethers_by_overlap(specified, reach)

Cluster the `specified` tethers with a union-find over `reach`
(point indices each tether touches): tethers whose reaches intersect
share structure and land in the same cluster. Returns a vector of
tether vectors, one per cluster.
"""
function group_tethers_by_overlap(specified, reach)
    n = length(specified)
    parent = collect(1:n)
    function find_root(i)
        parent[i] == i && return i
        parent[i] = find_root(parent[i])
        return parent[i]
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
    groups = Dict{Int64, Vector{Tether}}()
    for i in 1:n
        push!(get!(() -> Tether[], groups, find_root(i)), specified[i])
    end
    return collect(values(groups))
end

"""
    tether_unit_stiffness(tether, segments)

Return the common per-unit-length stiffness `[N]` of the tether's
segments. Errors if the segments are not uniform, since the spring
inversion in `apply_tether_init_forces!` assumes a single stiffness.
"""
function tether_unit_stiffness(tether, segments)
    ks = SimFloat[segments[si].unit_stiffness
                  for si in tether.segment_idxs]
    k = first(ks)
    all(≈(k), ks) || error("Tether $(tether.name): requires " *
        "uniform unit_stiffness across its segments, got $ks")
    return k
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
    cluster, points, segments, downstream, boundary; prn=true)
    snaps = map(cluster) do t
        anchor_idx, free_idx = tether_anchor_free(t, boundary)
        anchor_pos = copy(points[anchor_idx].pos_w)
        free_pos = copy(points[free_idx].pos_w)
        ordered = tether_ordered_point_idxs(t, segments)
        seg_lens = SimFloat[segment_world_length(segments[si], points)
                            for si in t.segment_idxs]
        if ordered[1] != anchor_idx
            reverse!(ordered)
            reverse!(seg_lens)
        end
        path_len = sum(seg_lens)
        path_len > 0 || error("Tether $(t.name): current length is " *
            "zero, cannot scale to its stretched length")
        (; t, free_idx, anchor_pos, free_pos, ordered, seg_lens, path_len)
    end

    deltas = [(s.t.init_stretched_len::SimFloat / s.path_len - 1) .*
              (s.free_pos .- s.anchor_pos) for s in snaps]
    delta = sum(deltas) ./ length(deltas)

    if length(cluster) > 1 && prn
        names = join((string(s.t.name) for s in snaps), ", ")
        @info "Tethers ($names) feed one structure; placing it to the " *
              "mean stretched length and direction of all."
    end
    norm(delta) ≈ 0 && return

    moved = Set{Int64}()
    for s in snaps
        for idx in downstream[s.t.idx]
            idx in moved && continue
            push!(moved, idx)
            points[idx].pos_w .+= delta
        end
        if !(s.free_idx in moved)
            push!(moved, s.free_idx)
            points[s.free_idx].pos_w .+= delta
        end
    end

    for s in snaps
        length(s.ordered) <= 2 && continue
        line = (s.free_pos .+ delta) .- s.anchor_pos
        cum = 0.0
        for k in 2:length(s.ordered)-1
            cum += s.seg_lens[k-1]
            points[s.ordered[k]].pos_w .=
                s.anchor_pos .+ (cum / s.path_len) .* line
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

    specified = [t for t in tethers if !isnothing(t.init_stretched_len)]
    isempty(specified) && return

    rigid_siblings = rigid_point_siblings(points, wings)

    boundary = Set{Int64}(w.winch_point_idx for w in winches)
    for point in points
        point.type == STATIC && push!(boundary, point.idx)
    end

    anchor_free = Dict(t.idx => tether_anchor_free(t, boundary)
                       for t in specified)
    non_root = [t for t in specified if isnothing(anchor_free[t.idx][1])]
    if !isempty(non_root)
        names = join((string(t.name) for t in non_root), ", ")
        error("tether length is only supported on tethers anchored at " *
              "a STATIC or winch point. Tether(s) ($names) have neither " *
              "endpoint anchored; their position rides the root tether.")
    end

    downstream = Dict(t.idx => tether_downstream_idxs(
                          t, segments, boundary, anchor_free[t.idx][2],
                          anchor_free[t.idx][1], rigid_siblings)
                      for t in specified)
    reach = Dict(t.idx => union(
        setdiff(Set{Int64}(tether_ordered_point_idxs(t, segments)), boundary),
        downstream[t.idx]) for t in specified)

    for cluster in group_tethers_by_overlap(specified, reach)
        apply_cluster_init_stretched_len!(cluster, points, segments,
                                          downstream, boundary; prn)
    end
end

"""
    apply_tether_init_forces!(sys_struct::SystemStructure)

Derive every tether's unstretched length `len` from its current
(placed) stretched length so the initial spring force equals
`init_tether_force` (default 0):
`len = stretched · (1 − force / unit_stiffness)` (zero-velocity,
tension branch of the segment spring law). Force 0 gives
`len = stretched` (zero tension).

Must be called after segment world lengths are current. Errors
if `force < 0` (compression unsupported), if `force ≥
unit_stiffness` (no positive rest length achieves it), or if a
tether's segments have non-uniform `unit_stiffness`.
"""
function apply_tether_init_forces!(sys_struct::SystemStructure)
    (; segments, tethers) = sys_struct
    for tether in tethers
        isempty(tether.segment_idxs) && continue
        stretched = sum(segments[si].len
                        for si in tether.segment_idxs)
        force = something(tether.init_tether_force, 0.0)
        force >= 0 || error("Tether $(tether.name): " *
            "init_tether_force $force N is negative; " *
            "compression is not supported")
        if force == 0
            tether.len = stretched
        else
            k = tether_unit_stiffness(tether, segments)
            force < k || error("Tether $(tether.name): " *
                "init_tether_force $force N ≥ unit_stiffness $k N; " *
                "no positive rest length achieves this force")
            tether.len = stretched * (1 - force / k)
        end
        tether.init_unstretched_len = tether.len
    end
end

# ==================== REINIT! FOR SYSTEM STRUCTURE ==================== #

"""
    reinit!(sys_struct::SystemStructure, set::Settings; kwargs...)

Re-initialize a `SystemStructure` from a `Settings` object.

This function resets various component states (e.g., winch lengths, group twists,
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
    (; points, groups, segments, pulleys, tethers, winches, wings, transforms) = sys_struct

    for winch in winches
        winch.vel = winch.init_vel
    end

    for group in groups
        group.twist = 0.0
        group.twist_ω = 0.0
    end

    # Transforms are not updated from Settings -
    # YAML structure geometry has priority

    # Step 1: copy CAD geometry to world frame
    copy_cad_to_world!(points, wings; update_vel=reset_vel)

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
        for si in tether.segment_idxs
            segments[si].l0 = l0
        end
    end

    for pulley in pulleys
        segment1, segment2 = segments[pulley.segment_idxs[1]],
                             segments[pulley.segment_idxs[2]]
        pulley.sum_len = segment1.l0 + segment2.l0

        # Initialize pulley.len proportional to current segment lengths
        # More accurate for asymmetric bridle configurations
        pulley.len = segment1.len / (segment1.len+segment2.len) *
                     pulley.sum_len

        pulley.vel = 0.0
    end

    # Step 5: apply transforms (translate/rotate/heading);
    # pos_w already initialized by copy_cad_to_world! +
    # apply_tether_init_stretched_lens!
    if apply_transforms
        reinit!(transforms, sys_struct; update_vel=reset_vel)
    end

    # Recreate VSM wing and aero if requested
    if remake_vsm
        for wing in wings
            wing isa VSMWing || continue
            # Recreate VSM wing from settings
            vsm_set = sys_struct.vsm_set::VortexStepMethod.VSMSettings
            wing.vsm_wing = create_vsm_wing(set, vsm_set;
                prn=false, sort_sections=false)
            wing.vsm_aero = VortexStepMethod.BodyAerodynamics([wing.vsm_wing])
            wing.vsm_solver = VortexStepMethod.Solver(wing.vsm_aero, vsm_set)

            # Transform sections: CAD → body frame
            # (must match SystemStructure constructor)
            vsm_wing = wing.vsm_wing
            vsm_wing.T_cad_body .= wing.pos_cad
            adjust_vsm_panels_to_origin!(
                vsm_wing, wing.pos_cad)
            rotate_vsm_sections!(
                vsm_wing, wing.R_b_to_c')
            vsm_wing.R_cad_body .= wing.R_b_to_c
            if wing.dynamics_type != PARTICLE_DYNAMICS
                apply_aero_z_offset!(
                    vsm_wing, wing.aero_z_offset)
            end
            VortexStepMethod.reinit!(wing.vsm_aero)

            # Match aero sections to structure (all types)
            match_aero_sections_to_structure!(
                wing, points; groups=groups)

            # Recompute group→section mapping
            if wing.dynamics_type == RIGID_DYNAMICS &&
               !isempty(wing.group_idxs)
                compute_spatial_group_mapping!(
                    wing, groups, points)
            end

            # PARTICLE_DYNAMICS-only: rebuild point mapping
            if wing.dynamics_type == PARTICLE_DYNAMICS &&
               !isnothing(wing.point_to_vsm_point)
                wing_point_idxs = collect(
                    keys(something(wing.point_to_vsm_point)))
                wing_pts = [points[idx]
                    for idx in wing_point_idxs]
                wing.point_to_vsm_point =
                    build_point_to_vsm_point_mapping(
                        wing_pts, wing)
            end
        end
    end

    # Compute per-wing wind from settings
    wind_vec_gnd = set.wind_vec

    for wing in wings
        # Calculate wind at wing position using atmospheric model
        wind_factor = calc_wind_factor(sys_struct.am,
                                       wing.pos_w[1], wing.pos_w[2],
                                       wing.pos_w[3], sys_struct)
        wing.v_wind .= wind_factor * wind_vec_gnd

        R_b_to_w = wing.R_b_to_w::Matrix{SimFloat}
        if wing.dynamics_type == PARTICLE_DYNAMICS
            va_wing_w = wing.v_wind - wing.vel_w + wing.wind_disturb
            wing.va_b .= R_b_to_w' * va_wing_w
        else
            # Initialize aero_y operating point
            if length(wing.aero_y) >= 2
                va_b_init = R_b_to_w' * wind_vec_gnd
                wing.aero_y .= 0.0
                wing.aero_y[1] = atan(
                    va_b_init[3], va_b_init[1])
                wing.aero_y[2] = atan(va_b_init[2],
                    hypot(va_b_init[1], va_b_init[3]))
            end
        end
    end

    # NOTE: validate_sys_struct() is called from model_management.jl after update_sys_struct!
    # because total_mass is only computed after the integrator exists.

    # Recalculate segment rest lengths from current positions if requested
    if ignore_l0
        for segment in segments
            p1 = points[segment.point_idxs[1]]
            p2 = points[segment.point_idxs[2]]
            segment.l0 = norm(p2.pos_w - p1.pos_w)
        end
    end

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
- It also copies the state of wings, groups, winches, and pulleys where applicable.
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

    # copy twist and twist_ω of groups
    if length(sys1.groups) > 0 && length(sys1.groups) == length(sys2.groups)
        for (group1, group2) in zip(sys1.groups, sys2.groups)
            group2.twist = group1.twist
            group2.twist_ω = group1.twist_ω
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
    update_from_sysstate!(sys::SystemStructure, ss::SysState)

Update the dynamic state of a `SystemStructure` from a `SysState` snapshot.

This function copies the state variables that are present in `SysState` (such as point
positions, wing orientations, winch lengths, and twist angles) into an existing `SystemStructure`.
Fields that cannot be populated from `SysState` (such as aerodynamic forces, moments, and
segment forces) are set to `NaN` to prevent them from being plotted.

This is useful for visualizing a `SysLog` by extracting individual `SysState` snapshots
and applying them to a `SystemStructure` for plotting with the Makie extension.

# Arguments
- `sys::SystemStructure`: The system structure to update (must already exist with correct topology).
- `ss::SysState`: The state snapshot to copy from.

# Example
```julia
# Load a system log
lg = load_log(...)

# Create a SystemStructure with the same topology
sys = SystemStructure(se(), "ram")

# Update from a specific time step
update_from_sysstate!(sys, lg.syslog[100])

# Plot the system at that time step
plot(sys)
```

# Notes
- The `SystemStructure` must have been created with the same model configuration as the
  simulation that generated the `SysLog`.
- Aerodynamic and force fields are set to `NaN` and will not be plotted.
- The number of points in `sys` must match the parametric type `P` of `SysState{P}`.
"""
function update_from_sysstate!(sys::SystemStructure, ss::SysState{P}) where P
    (; points, groups, tethers, winches, wings) = sys

    # Total slots: structural points + panel corners + wings.
    # Wing pos_w is appended after panel corners (see
    # update_sys_state!).
    n_points = length(points)
    n_panel_corners = isempty(wings) ? 0 : sum(
        length(wing.vsm_aero.panels) * 4 for wing in wings if wing isa VSMWing;
        init=0
    )
    n_wings = length(wings)
    total_with_wings = n_points + n_panel_corners + n_wings
    total_without_wings = n_points + n_panel_corners
    has_wing_slots = P == total_with_wings

    if !has_wing_slots && P != total_without_wings
        error("SystemStructure expects $total_with_wings points " *
              "($n_points regular + $n_panel_corners corners + " *
              "$n_wings wings) or $total_without_wings without " *
              "wing slots, but SysState has $P points")
    end

    # Update point positions (X, Y, Z from SysState)
    for point in points
        point.pos_w[1] = ss.X[point.idx]
        point.pos_w[2] = ss.Y[point.idx]
        point.pos_w[3] = ss.Z[point.idx]
        # Set velocity to zero (not available in basic SysState)
        point.vel_w .= 0.0
        # Set forces to NaN (not available in SysState)
        point.force .= NaN
    end

    # Update wing state if wings exist
    if length(wings) > 0 && length(wings) == 1  # Currently only support single-wing systems
        wing = wings[1]

        # Copy orientation quaternion
        wing.Q_b_to_w .= ss.orient

        # Copy spherical coordinates
        wing.elevation = Float64(ss.elevation)
        wing.azimuth = Float64(ss.azimuth)
        wing.heading = Float64(ss.heading)

        if has_wing_slots
            wing_slot = n_points + n_panel_corners + wing.idx
            wing.pos_w[1] = ss.X[wing_slot]
            wing.pos_w[2] = ss.Y[wing_slot]
            wing.pos_w[3] = ss.Z[wing_slot]
        else
            wing.pos_w .= [mean(ss.X), mean(ss.Y), mean(ss.Z)]
        end

        # Copy velocity if available in vel_kite
        wing.vel_w .= ss.vel_kite

        # Set angular velocity to NaN (turn_rates in SysState, but need conversion)
        wing.ω_b .= ss.turn_rates

        # Set aerodynamic quantities to NaN (to prevent plotting)
        wing.aero_force_b .= NaN
        wing.aero_moment_b .= NaN
        wing.tether_force .= NaN
        wing.tether_moment .= NaN
        wing.va_b .= NaN
        wing.v_wind .= ss.v_wind_kite
        wing.aoa = Float64(ss.AoA)
        wing.course = Float64(ss.course)
        wing.acc_w .= 0.0
        wing.turn_rate .= ss.turn_rates
        wing.turn_acc .= 0.0
    end

    # Update group twist angles
    n_groups = min(length(groups), 4)  # SysState stores up to 4 twist angles
    for i in 1:n_groups
        if i <= length(groups)
            groups[i].twist = Float64(ss.twist_angles[i])
            groups[i].twist_ω = 0.0  # Not available in SysState
            # Set forces/moments to NaN
            groups[i].tether_force = NaN
            groups[i].tether_moment = NaN
            groups[i].aero_moment = NaN
        end
    end

    # Update tether lengths from SysState (per-tether)
    for (ti, tether) in enumerate(tethers)
        ti > 4 && break
        tether.len = Float64(ss.l_tether[ti])
    end

    # Update winch state from SysState (per-winch)
    n_winches = min(length(winches), 4)
    for i in 1:n_winches
        winches[i].force .= NaN
        winches[i].friction = NaN
        winches[i].acc = 0.0
        winches[i].vel = Float64(ss.v_reelout[i])
        winches[i].set_value = Float64(ss.set_torque[i])
    end

    # Update VSM panel corner positions from world frame back to body frame
    corner_idx = n_points
    for wing in wings
        wing isa VSMWing || continue
        R_w_to_b = (wing.R_b_to_w::Matrix{SimFloat})'  # Transpose to get world-to-body rotation
        for panel in wing.vsm_aero.panels
            for j in 1:4
                corner_idx += 1
                # Get corner position from SysState (world frame)
                corner_w = [ss.X[corner_idx], ss.Y[corner_idx], ss.Z[corner_idx]]
                # Transform from world frame to body frame
                corner_b = R_w_to_b * (corner_w - wing.pos_w)
                # Update panel corner
                panel.corner_points[:, j] .= corner_b
            end
        end
    end

    # Update global wind vector (only if wind_vec mode is active)
    if sys.set.use_wind_vec
        sys.set.wind_vec = MVec3(ss.v_wind_gnd)
    end

    # Calculate segment lengths and forces from current positions and velocities
    # Note: velocities are set to zero, so damping term will be zero
    update_segment_forces!(sys)

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

    stretches = [s[2] for s in stretch_data]
    max_stretch = maximum(stretches)
    mean_stretch = sum(stretches) / length(stretches)
    max_idx = stretch_data[argmax(stretches)][1]

    return (max_stretch, mean_stretch, max_idx)
end
