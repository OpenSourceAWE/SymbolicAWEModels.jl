# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

"""
Utility functions for SystemStructure.

This file contains:
- Tether creation helpers
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

"""
    autocalc_tether_len(winch::Winch, tethers, segments)

Average unstretched tether length across all tethers connected
to this winch (sum of segment `l0` per tether, then average).
"""
function autocalc_tether_len(winch::Winch, tethers, segments)
    n = length(winch.tether_idxs)
    return sum(segments[seg_idx].l0
               for tether_idx in winch.tether_idxs
               for seg_idx in tethers[tether_idx].segment_idxs) / n
end

# ==================== TETHER CREATION ==================== #

"""
    create_tether(tether_idx, set, points, segments,
                  tethers, attach_point, dynamics_type;
                  z, unit_stiffness, unit_damping)

Procedurally create a multi-segment tether.

This function builds a tether from a specified number of
segments, connecting a given `attach_point` on the kite to
a new anchor point on the ground.
"""
function create_tether(tether_idx, set, points, segments,
                       tethers, attach_point,
                       dynamics_type; z=[0,0,1],
                       unit_stiffness=NaN,
                       unit_damping=NaN, d_pos=zeros(3))
    winch_pos = find_axis_point(
        attach_point.pos_cad, set.l_tether, z) .+ d_pos
    dir = winch_pos - attach_point.pos_cad
    segment_idxs = Int64[]
    ground_point_idx = 0
    for i in 1:set.segments
        frac = i / set.segments
        pos = attach_point.pos_cad + frac * dir
        point_idx = length(points) + 1
        segment_idx = length(segments) + 1
        if i == 1
            last_idx = attach_point.idx
        else
            last_idx = point_idx - 1
        end
        if i == set.segments
            points = [points;
                Point(point_idx, pos, STATIC)]
            ground_point_idx = points[end].idx
        else
            points = [points;
                Point(point_idx, pos, dynamics_type)]
        end
        segments = [segments;
            Segment(segment_idx, set, last_idx,
                    point_idx;
                    unit_stiffness, unit_damping)]
        push!(segment_idxs, segment_idx)
    end
    tethers = [tethers;
        Tether(tether_idx, segment_idxs)]
    return (points, segments, tethers,
            tethers[end].idx, ground_point_idx)
end

"""
    find_axis_point(P, l, v=[0,0,1])

Calculate the coordinates of a point `Q` that lies on a line defined by vector `v`
and is at a distance `l` from a given point `P`.
"""
function find_axis_point(P, l, v=[0,0,1])
    # Compute discriminant
    D = (v ⋅ P)^2 - norm(v)^2 * (norm(P)^2 - l^2)
    D < 0 && error("No real solution: l is too small or parameters invalid")
    # Solve quadratic for t, choose solution for negative direction
    t = (v ⋅ P - √D) / norm(v)^2
    # Compute point Q = t * v
    return [t * v[1], t * v[2], t * v[3]]
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
- Zero or near-zero principal inertia components on QUATERNION wings (error/warning)
- NaN inertia values (error)
- Empty group list for QUATERNION wings (warning)
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
    @unpack points, groups, segments, pulleys, wings, winches = sys_struct

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

        if wing.wing_type == QUATERNION
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

            # Warn if QUATERNION wing has no groups
            # (expected for AERO_NONE which skips auto-group creation)
            if isempty(wing.group_idxs) &&
               wing.aero_mode != AERO_NONE
                @warn "Wing $(wing.name) (QUATERNION)" *
                    " has no groups"
            end
        end

        # Check for NaN position (often a symptom of zero mass/inertia)
        if any(isnan.(wing.pos_w))
            error("Wing $(wing.name) has NaN position: pos_w = $(wing.pos_w)")
        end
        # REFINE wings don't use rigid body inertia, skip
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
    apply_tether_init_lens!(sys_struct::SystemStructure)

Scale tether point positions in `pos_w` to match each tether's `init_len`.
Must be called after `copy_cad_to_world!` (so `pos_w == pos_cad` at entry).

For each tether with a non-nothing `init_len`:
1. Scales all tether points (except start) radially from the start point.
2. Updates `segment.l0` for the tether's segments.
3. Sets `winch.tether_len = init_len` for any connected winch.
4. Propagates the end-point displacement via BFS through non-tether segments,
   translating downstream `pos_w` by the same delta.

Raises an error if a downstream non-tether segment connects back to the
tether's start point (would create an unsolvable constraint).
"""
function apply_tether_init_lens!(sys_struct::SystemStructure)
    @unpack points, segments, tethers, winches = sys_struct

    for tether in tethers
        isnothing(tether.init_len) && continue

        # Ordered point list: start → intermediates → end
        tether_point_idxs = Int64[tether.start_point_idx]
        for seg_idx in tether.segment_idxs
            push!(tether_point_idxs, segments[seg_idx].point_idxs[2])
        end

        # Always set winch tether_len when init_len is specified,
        # even if no geometric scaling is needed
        for winch in winches
            tether.idx in winch.tether_idxs || continue
            winch.tether_len = tether.init_len
        end

        current_len = sum(
            segment_world_length(segments[si], points)
            for si in tether.segment_idxs)
        current_len ≈ tether.init_len && continue
        current_len > 0 || error(
            "Tether $(tether.name): current length " *
            "is zero, cannot scale to init_len")

        scale = tether.init_len / current_len
        start_pos = copy(points[tether.start_point_idx].pos_w)
        old_end_pos = copy(points[tether.end_point_idx].pos_w)

        # Scale non-start tether points in pos_w
        for point_idx in tether_point_idxs[2:end]
            pt = points[point_idx]
            pt.pos_w .= start_pos .+ scale .* (pt.pos_w .- start_pos)
        end

        # Update segment l0 to match new positions
        for seg_idx in tether.segment_idxs
            seg = segments[seg_idx]
            seg.l0 = segment_world_length(seg, points)
        end

        delta = points[tether.end_point_idx].pos_w .- old_end_pos

        # BFS from end point through non-tether segments;
        # translate downstream pos_w by delta
        tether_segment_set = Set(tether.segment_idxs)
        visited = Set{Int64}(tether_point_idxs)
        queue = [tether.end_point_idx]

        while !isempty(queue)
            current_idx = popfirst!(queue)
            for seg in segments
                seg.idx in tether_segment_set && continue
                p1, p2 = seg.point_idxs
                neighbor_idx = if p1 == current_idx
                    p2
                elseif p2 == current_idx
                    p1
                else
                    continue
                end
                # Check start point before visited: start is in visited
                # but we still want to detect loops back to it
                if neighbor_idx == tether.start_point_idx
                    error("Tether $(tether.name): downstream structure " *
                          "connects back to tether start point. " *
                          "Cannot apply init_len scaling.")
                end
                neighbor_idx in visited && continue
                points[neighbor_idx].pos_w .+= delta
                push!(visited, neighbor_idx)
                push!(queue, neighbor_idx)
            end
        end
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
  For REFINE wings, also rebuilds the `point_to_vsm_point` mapping.
"""
function reinit!(sys_struct::SystemStructure, set::Settings;
                 ignore_l0::Bool=false, remake_vsm::Bool=false,
                 reset_vel::Bool=true)
    @unpack points, groups, segments, pulleys, tethers, winches, wings, transforms = sys_struct

    for winch in winches
        winch.tether_vel = winch.init_vel
    end

    for group in groups
        group.twist = 0.0
        group.twist_ω = 0.0
    end

    # Transforms are not updated from Settings - YAML structure geometry has priority

    # Step 1: copy CAD geometry to world frame (all points and wings)
    copy_cad_to_world!(points, wings; update_vel=reset_vel)

    # Step 2: apply tether initial lengths (scales pos_w; pos_cad unchanged)
    apply_tether_init_lens!(sys_struct)

    # Step 3: compute segment lengths from pos_w
    for segment in segments
        len = segment_world_length(segment, points)
        (segment.l0 ≈ 0) && (segment.l0 = len)
        segment.len = len
    end

    # Step 4: set winch tether_len from segment l0s
    # (skip if init_len already set it in apply_tether_init_lens!)
    for winch in winches
        any(!isnothing(tethers[ti].init_len)
            for ti in winch.tether_idxs) && continue
        winch.tether_len = autocalc_tether_len(
            winch, tethers, segments)
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

    # Calculate ground-level wind vector BEFORE transforms (needed for heading calculation)
    # Matches symbolic equations in generate_system.jl:1259-1264
    upwind_dir = deg2rad(set.upwind_dir)
    wind_elevation = sys_struct.wind_elevation
    wind_scale_gnd = set.v_wind

    wind_vec_base = [0.0, -1.0, 0.0]
    wind_vec_elevated = rotate_around_x(wind_vec_base, wind_elevation)
    wind_vec_rotated = rotate_around_z(wind_vec_elevated, -upwind_dir)
    sys_struct.wind_vec_gnd .= max(wind_scale_gnd, 1e-6) * wind_vec_rotated

    # Step 5: apply transforms (translate/rotate/heading);
    # pos_w already initialized by copy_cad_to_world! + apply_tether_init_lens!
    reinit!(transforms, sys_struct; update_vel=reset_vel)

    # Recreate VSM wing and aero if requested
    if remake_vsm
        for wing in wings
            wing isa VSMWing || continue
            # Recreate VSM wing from settings
            vsm_set = sys_struct.vsm_set::VortexStepMethod.VSMSettings
            wing.vsm_wing = create_vsm_wing(set, vsm_set;
                prn=false, sort_sections=false)
            wing.vsm_aero = VortexStepMethod.BodyAerodynamics([wing.vsm_wing])
            wing.vsm_solver = VortexStepMethod.Solver(wing.vsm_aero;
                solver_type=VortexStepMethod.NONLIN,
                atol=2e-8, rtol=2e-8)

            # Transform sections: CAD → body frame
            # (must match SystemStructure constructor)
            vsm_wing = wing.vsm_wing
            vsm_wing.T_cad_body .= wing.pos_cad
            adjust_vsm_panels_to_origin!(
                vsm_wing, wing.pos_cad)
            rotate_vsm_sections!(
                vsm_wing, wing.R_b_to_c')
            vsm_wing.R_cad_body .= wing.R_b_to_c
            if wing.wing_type != REFINE
                apply_aero_z_offset!(
                    vsm_wing, wing.aero_z_offset)
            end
            VortexStepMethod.reinit!(wing.vsm_aero)

            # Match aero sections to structure (all types)
            match_aero_sections_to_structure!(
                wing, points; groups=groups)

            # Recompute group→section mapping
            if wing.wing_type == QUATERNION &&
               !isempty(wing.group_idxs)
                compute_spatial_group_mapping!(
                    wing, groups, points)
            end

            # REFINE-only: rebuild point mapping
            if wing.wing_type == REFINE &&
               !isnothing(wing.point_to_vsm_point)
                wing_point_idxs = collect(
                    keys(something(wing.point_to_vsm_point)))
                wing_pts = [points[idx]
                    for idx in wing_point_idxs]
                wing.point_to_vsm_point =
                    build_point_to_vsm_point_mapping(
                        wing_pts, vsm_wing)
            end
        end
    end

    # Calculate ground-level wind vector with direction rotations
    # Matches symbolic equations in generate_system.jl:1259-1264
    upwind_dir = deg2rad(set.upwind_dir)
    wind_elevation = sys_struct.wind_elevation
    wind_scale_gnd = set.v_wind

    # Base wind vector: [0, -1, 0] points upwind
    wind_vec_base = [0.0, -1.0, 0.0]
    # Rotate by elevation around x-axis (vertical tilt)
    wind_vec_elevated = rotate_around_x(wind_vec_base, wind_elevation)
    # Rotate by upwind direction around z-axis (negative for convention)
    wind_vec_rotated = rotate_around_z(wind_vec_elevated, -upwind_dir)
    # Scale by ground wind speed
    wind_vec_gnd = max(wind_scale_gnd, 1e-6) * wind_vec_rotated

    for wing in wings
        # Calculate wind at wing position using atmospheric model
        # Matches symbolic equations in generate_system.jl
        wind_factor = calc_wind_factor(sys_struct.am,
                                       wing.pos_w[1], wing.pos_w[2], wing.pos_w[3], set)
        wing.v_wind .= wind_factor * wind_vec_gnd

        R_b_to_w = wing.R_b_to_w::Matrix{SimFloat}
        if wing.wing_type == REFINE
            # Initialize apparent wind in body frame for REFINE wings
            # va_wing = wind_vel - wing_vel + wind_disturb
            # va_b = R_b_to_w' * va_wing
            # At initialization: wing_vel typically 0, wind_disturb typically 0
            va_wing_w = wing.v_wind - wing.vel_w + wing.wind_disturb
            wing.va_b .= R_b_to_w' * va_wing_w
        else
            # Initialize vsm_y for QUATERNION wings (REFINE wings have ny=0)
            if length(wing.vsm_y) >= 3
                wing.vsm_y .= 0.0
                wing.vsm_y[1:3] .= R_b_to_w' * [set.v_wind, 0., 0.]
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
    simple = false

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
                    simple = true
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

    # copy winch tether lengths and velocities
    if length(sys1.winches) > 0 && length(sys1.winches) == length(sys2.winches)
        for (winch2, winch1) in zip(sys2.winches, sys1.winches)
            if !simple
                winch2.tether_len = winch1.tether_len
                winch2.tether_vel = winch1.tether_vel
            else
                tether_len_acc = 0.0
                for tether_idx in winch1.tether_idxs
                    tether2 = sys2.tethers[tether_idx]
                    segment2 = sys2.segments[tether2.segment_idxs[1]]
                    point_idxs2 = segment2.point_idxs
                    slen = norm(sys2.points[point_idxs2[1]].pos_w .-
                                        sys2.points[point_idxs2[2]].pos_w)
                    stiffness = segment2.unit_stiffness / slen
                    nt = length(winch1.tether_idxs)
                    tether_len_acc += (slen - norm(winch1.force)/stiffness/nt) / nt
                end
                winch2.tether_len = tether_len_acc
                winch2.tether_vel = winch1.tether_vel
            end
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
    @unpack points, groups, winches, wings = sys

    # Calculate expected total points (regular points + panel corners)
    n_points = length(points)
    n_panel_corners = isempty(wings) ? 0 : sum(
        length(wing.vsm_aero.panels) * 4 for wing in wings if wing isa VSMWing;
        init=0
    )
    expected_total = n_points + n_panel_corners

    # Verify compatibility
    if expected_total != P
        error("SystemStructure expects $expected_total points " *
              "($n_points regular + $n_panel_corners corners) " *
              "but SysState has $P points")
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

        # Compute wing position from average of points (if wings exist)
        # For a typical system, the wing COM is near the bridle attachment
        wing.pos_w .= [mean(ss.X), mean(ss.Y), mean(ss.Z)]

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

    # Update winch state
    n_winches = min(length(winches), 4)  # SysState stores up to 4 winches
    for i in 1:n_winches
        if i <= length(winches)
            winches[i].tether_len = Float64(ss.l_tether[i])
            winches[i].tether_vel = Float64(ss.v_reelout[i])
            winches[i].force .= NaN  # Force not directly available
            winches[i].friction = NaN
            winches[i].tether_acc = 0.0
            winches[i].set_value = Float64(ss.set_torque[i])
        end
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

    # Update global wind vector
    sys.wind_vec_gnd .= ss.v_wind_gnd

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
