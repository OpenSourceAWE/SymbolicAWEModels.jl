# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

"""
SystemStructure type and main constructor.

This file contains:
- SystemStructure struct definition
- Property accessors (getproperty/setproperty!)
- Main SystemStructure constructor with initialization logic
"""

# ==================== SYSTEM STRUCTURE ==================== #

"""
    struct SystemStructure

A discrete mass-spring-damper representation of a kite system.

This struct holds all components of the physical model, including points, segments,
winches, and wings, forming a complete description of the kite system's structure.

# Components
- [`Point`](@ref): Point masses.
- [`TwistSurface`](@ref): Collections of points for wing deformation.
- [`Segment`](@ref): Spring-damper elements.
- [`Pulley`](@ref): Elements that redistribute line lengths.
- [`Tether`](@ref): TwistSurfaces of segments controlled by a winch.
- [`Winch`](@ref): Ground-based winches.
- [`Wing`](@ref): Rigid wing bodies.
- [`Transform`](@ref): Spatial transformations for initial positioning.
"""
mutable struct SystemStructure{W<:AbstractWing}
    const name::String
    set::Settings
    const points::NamedCollection{Point}
    const twist_surfaces::NamedCollection{TwistSurface}
    const segments::NamedCollection{Segment}
    const pulleys::NamedCollection{Pulley}
    const tethers::NamedCollection{Tether}
    const winches::NamedCollection{Winch}
    const wings::NamedCollection{W}
    const transforms::NamedCollection{Transform}

    const am::AtmosphericModel
    stabilize::Bool
    fix_wing::Bool
    vsm_set::Union{Nothing, VortexStepMethod.VSMSettings}
end

function Base.getproperty(sys::SystemStructure, sym::Symbol)
    if sym == :total_mass
        # Sum of all point total_mass values (computed during simulation)
        # Falls back to extra_mass if total_mass is 0
        total = 0.0
        for point in getfield(sys, :points)
            if point.total_mass > 0
                total += point.total_mass
            else
                total += point.extra_mass
            end
        end
        return total
    elseif sym == :diff_vars
        vars = SimFloat[]
        # points
        points = getfield(sys, :points)
        for point in points
            if point.type == DYNAMIC
                append!(vars, point.pos_w)
                append!(vars, point.vel_w)
            end
        end
        # wings (principal frame ODE state, RIGID_DYNAMICS only)
        wings = getfield(sys, :wings)
        for wing in wings
            wing.dynamics_type != RIGID_DYNAMICS && continue
            append!(vars, wing.com_w)
            append!(vars, wing.com_vel)
            append!(vars, wing.Q_p_to_w)
            append!(vars, wing.ω_p)
        end
        # twist_surfaces
        twist_surfaces = getfield(sys, :twist_surfaces)
        for twist_surface in twist_surfaces
            if twist_surface.type == DYNAMIC
                push!(vars, twist_surface.twist)
                push!(vars, twist_surface.twist_ω)
            end
        end
        # pulleys
        pulleys = getfield(sys, :pulleys)
        for pulley in pulleys
            if pulley.type == DYNAMIC
                push!(vars, pulley.len)
                push!(vars, pulley.vel)
            end
        end
        # tethers
        tethers = getfield(sys, :tethers)
        for tether in tethers
            push!(vars, tether.len)
        end
        # winches
        winches = getfield(sys, :winches)
        for winch in winches
            push!(vars, winch.vel)
        end
        return reshape(vars, :, 1) # Return as a column vector (2D array)
    else
        return getfield(sys, sym)
    end
end

function Base.setproperty!(sys::SystemStructure, sym::Symbol, value)
    if sym == :diff_vars
        flat_value = vec(value) # Ensure value is a flat vector
        offset = 1
        # points
        points = getfield(sys, :points)
        for point in points
            if point.type == DYNAMIC
                point.pos_w .= @view flat_value[offset:offset+2]
                offset += 3
                point.vel_w .= @view flat_value[offset:offset+2]
                offset += 3
            end
        end
        # wings (principal frame ODE state, RIGID_DYNAMICS only)
        wings = getfield(sys, :wings)
        for wing in wings
            wing.dynamics_type != RIGID_DYNAMICS && continue
            wing.com_w .= @view flat_value[offset:offset+2]
            offset += 3
            wing.com_vel .= @view flat_value[offset:offset+2]
            offset += 3
            wing.Q_p_to_w .= @view flat_value[offset:offset+3]
            offset += 4
            wing.ω_p .= @view flat_value[offset:offset+2]
            offset += 3
        end
        # twist_surfaces
        twist_surfaces = getfield(sys, :twist_surfaces)
        for twist_surface in twist_surfaces
            if twist_surface.type == DYNAMIC
                twist_surface.twist = flat_value[offset]
                offset += 1
                twist_surface.twist_ω = flat_value[offset]
                offset += 1
            end
        end
        # pulleys
        pulleys = getfield(sys, :pulleys)
        for pulley in pulleys
            if pulley.type == DYNAMIC
                pulley.len = flat_value[offset]
                offset += 1
                pulley.vel = flat_value[offset]
                offset += 1
            end
        end
        # tethers
        tethers = getfield(sys, :tethers)
        for tether in tethers
            tether.len = flat_value[offset]
            offset += 1
        end
        # winches
        winches = getfield(sys, :winches)
        for winch in winches
            winch.vel = flat_value[offset]
            offset += 1
        end
        return value
    else
        return setfield!(sys, sym, value)
    end
end

"""
    calc_heading(sys::SystemStructure)

Calculate heading angles for all wings using the tangential
sphere frame method. Returns a vector of heading angles, one
per wing.
"""
function calc_heading(sys::SystemStructure)
    return [calc_heading(wing.R_b_to_w, wing.pos_w)
            for wing in getfield(sys, :wings)]
end

"""
    build_name_dict(items::Vector) -> Dict{Symbol, Int64}

Build a name→index dictionary from a vector of items with optional `name` fields.
Items with `name=nothing` are skipped. Integer names are converted to Symbols.
"""
function build_name_dict(items::Vector)
    name_to_idx = Dict{Symbol, Int64}()
    for (i, item) in enumerate(items)
        # Use try-catch to handle both direct fields and delegated properties (e.g., VSMWing.name -> base.name)
        item_name = try
            item.name
        catch
            nothing
        end
        if !isnothing(item_name)
            name = item_name isa Symbol ? item_name : Symbol(item_name)
            if haskey(name_to_idx, name)
                error("Duplicate name '$name' found at indices $(name_to_idx[name]) and $i")
            end
            name_to_idx[name] = i
        end
    end
    return name_to_idx
end

"""
    resolve_ref(ref::NameRef, name_dict::Dict{Symbol, Int64}, component_type::String) -> Int64

Resolve a reference (name or index) to an index using the name dictionary.
If ref is an integer, returns it directly. If ref is a symbol, looks up in dictionary.
"""
function resolve_ref(ref::Union{Int, Symbol}, name_dict::Dict{Symbol, Int64}, component_type::String)
    if ref isa Int
        return Int64(ref)
    else
        name = ref isa Symbol ? ref : Symbol(ref)
        if haskey(name_dict, name)
            return name_dict[name]
        else
            error("Unknown $component_type name: $name")
        end
    end
end

function resolve_ref(::Nothing, ::Dict{Symbol, Int64}, ::String)
    return Int64(0)
end

"""
    resolve_ref_spec(spec, name_dict, component_type) -> Union{Int64, Vector{Int64}, Nothing}

Resolve a reference point specification (single ref or vector of refs) to indices.
"""
function resolve_ref_spec(spec::Union{Int, Symbol}, name_dict::Dict{Symbol, Int64}, component_type::String)
    return resolve_ref(spec, name_dict, component_type)
end

function resolve_ref_spec(spec::AbstractVector, name_dict::Dict{Symbol, Int64}, component_type::String)
    return Int64[resolve_ref(ref, name_dict, component_type) for ref in spec]
end

function resolve_ref_spec(::Nothing, ::Dict{Symbol, Int64}, ::String)
    return nothing
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
    expand_auto_tethers!(points, segments, tethers, set)

For Route 2 tethers (auto-generation), create intermediate DYNAMIC
points and segments. Must be called before `assign_indices_and_resolve!`.

Detects Route 2 tethers by checking `start_point_ref !== nothing`
and `segment_refs` names not yet present in `segments`.
"""
function expand_auto_tethers!(
    points::Vector{Point},
    segments::Vector{Segment},
    tethers::Vector{Tether},
    set::Settings
)
    # Build point name lookup for finding start/end points
    point_names = Dict{Symbol, Int}()
    for (i, point) in enumerate(points)
        if !isnothing(point.name)
            name = point.name isa Symbol ? point.name : Symbol(point.name)
            point_names[name] = i
        end
    end

    # Build segment name set for checking existing segments
    seg_names = Set{Symbol}()
    for segment in segments
        if !isnothing(segment.name)
            name = segment.name isa Symbol ? segment.name : Symbol(segment.name)
            push!(seg_names, name)
        end
    end

    for tether in tethers
        # Skip Route 1 tethers (no start_point_ref)
        isnothing(tether.start_point_ref) && continue

        # Skip if segments already exist (user pre-created)
        first_seg_name = tether.segment_refs[1]
        first_seg_sym = first_seg_name isa Symbol ?
            first_seg_name : Symbol(first_seg_name)
        first_seg_sym in seg_names && continue

        # Resolve start and end points by name
        start_point_sym = tether.start_point_ref isa Symbol ?
            tether.start_point_ref :
            Symbol(tether.start_point_ref)
        end_point_sym = tether.end_point_ref isa Symbol ?
            tether.end_point_ref :
            Symbol(tether.end_point_ref)
        haskey(point_names, start_point_sym) || error(
            "Tether $(tether.name): start_point " *
            ":$start_point_sym not found in points")
        haskey(point_names, end_point_sym) || error(
            "Tether $(tether.name): end_point " *
            ":$end_point_sym not found in points")
        start_point_idx = point_names[start_point_sym]
        end_point_idx = point_names[end_point_sym]
        start_pos = points[start_point_idx].pos_cad
        end_pos = points[end_point_idx].pos_cad

        # Inherit transform from endpoints
        start_transform = points[start_point_idx].transform_ref
        end_transform = points[end_point_idx].transform_ref
        start_has_transform = start_transform != 0
        end_has_transform = end_transform != 0
        if start_has_transform && end_has_transform &&
           start_transform != end_transform
            error("Tether $(tether.name): " *
                "start_point :$start_point_sym (transform=" *
                "$start_transform) and end_point :$end_point_sym " *
                "(transform=$end_transform) have different " *
                "transforms")
        end
        transform_idx = start_has_transform ? start_transform :
            end_has_transform ? end_transform : 0

        n = tether.n_segments

        # Derive segment properties from settings if NaN
        unit_stiffness = tether.unit_stiffness
        unit_damping = tether.unit_damping
        diameter = tether.diameter
        if isnan(diameter)
            diameter = set.d_tether * 0.001  # mm → m
        end
        if isnan(unit_stiffness)
            unit_stiffness = set.e_tether * (diameter / 2)^2 * π
        end
        if isnan(unit_damping)
            if hasproperty(set, :rel_damping) &&
               set.rel_damping != 0.0
                unit_damping = set.rel_damping * unit_stiffness
            else
                unit_damping = 0.0
            end
        end

        rope_len = isnothing(tether.init_stretched_len) ?
            norm(end_pos - start_pos) : tether.init_stretched_len
        seg_l0 = rope_len / n
        tether.len = rope_len

        # Generate n-1 intermediate DYNAMIC points
        # (placed along the straight line at geometric spacing)
        direction = end_pos - start_pos
        for i in 1:(n - 1)
            frac = i / n
            pos = start_pos + frac * direction
            point_name = Symbol("$(tether.name)_point_$i")
            transform_kw = transform_idx == 0 ? () :
                (transform=transform_idx,)
            push!(points, Point(point_name, pos,
                DYNAMIC; transform_kw...))
            point_names[point_name] = length(points)
        end

        # Generate n segments
        for i in 1:n
            seg_name = tether.segment_refs[i]
            seg_sym = seg_name isa Symbol ? seg_name :
                Symbol(seg_name)
            if i == 1
                start_ref = tether.start_point_ref
            else
                start_ref = Symbol(
                    "$(tether.name)_point_$(i - 1)")
            end
            if i == n
                end_ref = tether.end_point_ref
            else
                end_ref = Symbol(
                    "$(tether.name)_point_$i")
            end
            push!(segments, Segment(
                seg_sym, start_ref, end_ref,
                unit_stiffness, unit_damping, diameter;
                l0=seg_l0))
            push!(seg_names, seg_sym)
        end
    end
end

"""
    assign_indices_and_resolve!(components, name_dicts)

Assign indices to all components based on their position in the vectors,
and resolve all references to indices.
"""
function assign_indices_and_resolve!(
    points::Vector{Point},
    twist_surfaces::Vector{TwistSurface},
    segments::Vector{Segment},
    pulleys::Vector{Pulley},
    tethers::Vector{Tether},
    winches::Vector{Winch},
    wings::Vector{<:AbstractWing},
    transforms::Vector{Transform}
)
    # Build name dictionaries FIRST (using idx values after assignment)
    # First pass: assign indices based on position
    for (i, point) in enumerate(points)
        point.idx = i
    end
    for (i, twist_surface) in enumerate(twist_surfaces)
        twist_surface.idx = i
    end
    for (i, segment) in enumerate(segments)
        segment.idx = i
    end
    for (i, pulley) in enumerate(pulleys)
        pulley.idx = i
    end
    for (i, tether) in enumerate(tethers)
        tether.idx = i
    end
    for (i, winch) in enumerate(winches)
        winch.idx = i
    end
    for (i, wing) in enumerate(wings)
        wing.idx = i
    end
    for (i, transform) in enumerate(transforms)
        transform.idx = i
    end

    # Build name lookup dictionaries
    point_names = build_name_dict(points)
    twist_surface_names = build_name_dict(twist_surfaces)
    segment_names = build_name_dict(segments)
    pulley_names = build_name_dict(pulleys)
    tether_names = build_name_dict(tethers)
    winch_names = build_name_dict(winches)
    wing_names = build_name_dict(wings)
    transform_names = build_name_dict(transforms)

    # Resolve references for all components
    # Points: resolve wing_ref and transform_ref
    for point in points
        point.wing_idx = resolve_ref(point.wing_ref, wing_names, "wing")
        point.transform_idx = resolve_ref(point.transform_ref, transform_names, "transform")
    end

    # TwistSurfaces: resolve point_refs
    for twist_surface in twist_surfaces
        twist_surface.point_idxs = Int64[resolve_ref(ref, point_names, "point") for ref in twist_surface.point_refs]
    end

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

    # Transforms: resolve wing_ref, rot_point_ref, base_point_ref, base_transform_ref
    # Use resolve_ref_spec (returns nothing for nothing inputs) since Transform
    # fields are Union{Int64, Nothing}.
    for transform in transforms
        transform.wing_idx = resolve_ref_spec(transform.wing_ref, wing_names, "wing")
        transform.rot_point_idx = resolve_ref_spec(transform.rot_point_ref, point_names, "point")
        transform.base_point_idx = resolve_ref_spec(transform.base_point_ref, point_names, "point")
        transform.base_transform_idx = resolve_ref_spec(transform.base_transform_ref, transform_names, "transform")
    end

    # Wings: resolve twist_surface_refs, transform_ref, and PARTICLE_DYNAMICS-specific refs
    for wing in wings
        # BaseWing fields
        wing.twist_surface_idxs = Int64[resolve_ref(ref, twist_surface_names, "twist_surface") for ref in wing.twist_surface_refs]
        wing.transform_idx = resolve_ref(wing.transform_ref, transform_names, "transform")

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

        # Resize per-mode aero state now that twist_surface_idxs are resolved
        # (initial sizing used n_unrefined as proxy which may differ)
        resize_aero_state!(wing.aero, wing)
    end

    return (point_names, twist_surface_names, segment_names, pulley_names,
            tether_names, winch_names, wing_names, transform_names)
end

"""
    init_body_frame_from_ref_points!(wing, points; prn=true)

Initialize wing body frame (R_b_to_c, pos_cad) from z/y
reference points. Shared by VSMWing PARTICLE_DYNAMICS and PlateWing.
"""
function init_body_frame_from_ref_points!(
    wing, points; prn=true
)
    isnothing(wing.origin) && return
    isnothing(wing.z_ref_points) && return
    isnothing(wing.y_ref_points) && return

    origin_pos = get_ref_position_from_points(
        points, wing.origin; field=:pos_cad)
    wing.pos_cad .= origin_pos

    # Temporarily set pos_w = pos_cad so
    # calc_particle_dynamics_wing_frame can read positions
    for point in points
        point.type == WING && point.wing_idx == wing.idx &&
            (point.pos_w .= point.pos_cad)
    end
    R_b_to_c, _ = calc_particle_dynamics_wing_frame(
        points, wing.z_ref_points,
        wing.y_ref_points, wing.origin)
    wing.R_b_to_c .= R_b_to_c
    wing.R_b_to_p .= Matrix{SimFloat}(I, 3, 3)

    if prn
        origin_rounded = round.(origin_pos; digits=3)
        @info "Wing $(wing.name) ($(typeof(wing).name.name))" *
              ": origin=[$origin_rounded]"
    end
end

"""
    calc_inertia_y_rotation(I_tensor)

Find the Y-axis rotation that diagonalizes the XZ block
of the inertia tensor (zeros out `I[1,3]` and `I[3,1]`).

Returns `(I_diag, Ry)` where `I_diag = Ry * I_tensor * Ry'`
and `Ry` is a rotation about the Y axis by angle
`θ = atan(2·I₁₃, I₁₁ − I₃₃) / 2`.
"""
function calc_inertia_y_rotation(I_tensor)
    θ = atan(2 * I_tensor[1, 3],
             I_tensor[1, 1] - I_tensor[3, 3]) / 2
    cθ, sθ = cos(θ), sin(θ)
    Ry = [cθ 0 sθ; 0 1 0; -sθ 0 cθ]
    I_diag = Ry * I_tensor * Ry'
    return I_diag, Ry
end

"""
    compute_spatial_twist_surface_mapping!(the_wing, twist_surfaces, points)

Partition the wing's unrefined VSM sections among its
twist_surfaces by spatial proximity: each unrefined section is
assigned to the single closest twist_surface (by distance between
section centre and twist_surface centre, both in body frame).

When `n_twist_surfaces == n_unrefined` this is the same 1:1
mapping as before. When `n_twist_surfaces < n_unrefined` a twist_surface
may own several adjacent sections; its single twist DOF
then drives all of them as a rigid unit. The case
`n_twist_surfaces > n_unrefined` is rejected — a twist DOF without
a section to drive would be undefined.
"""
function compute_spatial_twist_surface_mapping!(
    the_wing::Wing,
    twist_surfaces::AbstractVector{TwistSurface},
    points::AbstractVector{Point}
)
    the_vsm_wing = the_wing.vsm_wing
    n_unrefined = the_vsm_wing.n_unrefined_sections
    n_twist_surfaces = length(the_wing.twist_surface_idxs)

    n_twist_surfaces <= n_unrefined || error(
        "Wing $(the_wing.idx): n_twist_surfaces " *
        "($n_twist_surfaces) > n_unrefined sections " *
        "($n_unrefined). Reduce twist_surfaces or increase " *
        "aero resolution.")

    # Compute twist_surface centers in body frame
    twist_surface_centers = Vector{MVec3}(undef, n_twist_surfaces)
    for (local_idx, twist_surface_idx) in
            enumerate(the_wing.twist_surface_idxs)
        twist_surface = twist_surfaces[twist_surface_idx]
        center = zeros(3)
        for pt_idx in twist_surface.point_idxs
            center += the_wing.R_b_to_c' *
                (points[pt_idx].pos_cad -
                 the_wing.pos_cad)
        end
        twist_surface_centers[local_idx] =
            center / length(twist_surface.point_idxs)
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

    # Reset section lists (we rebuild the partition)
    for twist_surface_idx in the_wing.twist_surface_idxs
        empty!(twist_surfaces[twist_surface_idx].unrefined_section_idxs)
    end

    # Assign each unrefined section to nearest twist_surface
    for section_idx in 1:n_unrefined
        min_dist = Inf
        closest_local = 1
        for local_idx in 1:n_twist_surfaces
            dist = norm(unrefined_centers[section_idx] -
                     twist_surface_centers[local_idx])
            if dist < min_dist
                min_dist = dist
                closest_local = local_idx
            end
        end
        g_idx = the_wing.twist_surface_idxs[closest_local]
        push!(twist_surfaces[g_idx].unrefined_section_idxs,
              Int64(section_idx))
    end

    # Every twist_surface must claim at least one section
    for twist_surface_idx in the_wing.twist_surface_idxs
        twist_surface = twist_surfaces[twist_surface_idx]
        isempty(twist_surface.unrefined_section_idxs) && error(
            "Wing $(the_wing.idx): twist_surface " *
            "$(twist_surface.name) claims no unrefined " *
            "sections (likely coincident twist_surface centres).")
    end
end

# ==================== CONSTRUCTOR ==================== #

"""
    setup_wing_frame!(wing, points; prn=true)

Compute a wing's body frame (`R_b_to_c`, `pos_cad`) and, for `RIGID_DYNAMICS`, its
COM offset and principal inertia, from the WING points and ref points. This is
dynamics/geometry only — independent of the aero mode, which does its own
mode-specific setup afterwards in [`setup_aero!`](@ref).
"""
function setup_wing_frame!(wing, points; prn=true)
    if wing.dynamics_type == RIGID_DYNAMICS
        any(point.type == WING && point.wing_idx == wing.idx
            for point in points) || return nothing

        com_cad, inertia_normalized = normalized_inertia(wing.aero, wing, points)
        if !isnothing(inertia_normalized)
            # The hook returns per-unit-mass inertia [m²]; scale once here.
            I_cad = wing.mass .* inertia_normalized
            I_diag, Ry = calc_inertia_y_rotation(I_cad)
            wing.R_p_to_c .= Ry'  # principal→CAD
            wing.inertia_principal .= diag(I_diag)
        end

        # Body frame from ref points (else body = principal, origin = COM)
        if !isnothing(wing.origin) &&
           !isnothing(wing.z_ref_points) &&
           !isnothing(wing.y_ref_points)
            origin_cad = get_ref_position_from_points(
                points, wing.origin; field=:pos_cad)
            wing.pos_cad .= origin_cad
            for point in points
                point.type == WING &&
                    point.wing_idx == wing.idx &&
                    (point.pos_w .= point.pos_cad)
            end
            R_b_to_c, _ = calc_particle_dynamics_wing_frame(
                points, wing.z_ref_points, wing.y_ref_points, wing.origin)
            wing.R_b_to_c .= R_b_to_c
            wing.com_offset_b .= R_b_to_c' * (com_cad - origin_cad)
        else
            wing.pos_cad .= com_cad
            wing.R_b_to_c .= wing.R_p_to_c
            wing.com_offset_b .= 0.0
        end

        wing.R_b_to_p .= wing.R_p_to_c' * wing.R_b_to_c  # body → principal

        if prn
            I_rnd = round.(wing.inertia_principal; digits=4)
            offset_rounded = round.(wing.com_offset_b; digits=4)
            @info "RIGID_DYNAMICS wing $(wing.idx):" *
                " COM=[$(round.(com_cad; digits=3))]" *
                ", I=$I_rnd, com_offset_b=$offset_rounded"
        end
    else  # PARTICLE_DYNAMICS (VSM or flat-plate)
        init_body_frame_from_ref_points!(wing, points; prn)
    end
    return nothing
end

"""
    SystemStructure(name, set; points, twist_surfaces, segments, pulleys, tethers, winches, wings, transforms)

Constructs a `SystemStructure` object representing a complete kite system.

## Physical Models
- **"ram"**: A model with 4 deformable wing twist_surfaces and a complex pulley bridle system.
- **"simple_ram"**: A model with 4 deformable wing twist_surfaces and direct bridle connections.

# Arguments
- `name::String`: Model identifier ("ram", "simple_ram", or a custom name).
- `set::Settings`: Configuration parameters from `KiteUtils.jl`.

# Keyword Arguments
- `points`, `twist_surfaces`, `segments`, etc.: Vectors of the system components.
- `prn::Bool=true`: If true, print info messages about auto-generated components.

# Returns
- `SystemStructure`: A complete system ready for building a `SymbolicAWEModel`.
"""
function SystemStructure(name, set;
        points=Point[],
        twist_surfaces=TwistSurface[],
        segments=Segment[],
        pulleys=Pulley[],
        tethers=Tether[],
        winches=Winch[],
        wings=AbstractWing[],
        transforms=Transform[],
        ignore_l0::Bool=false,
        vsm_set=nothing,
        prn::Bool=true,
    )
    # Load VSMSettings if not provided and VSM wings exist
    has_vsm_wings = any(has_vsm_engine(wing.aero) for wing in wings)
    if isnothing(vsm_set) && has_vsm_wings
        model_dir = get_data_path()
        vsm_set_path = joinpath(model_dir, "vsm_settings.yaml")
        if isfile(vsm_set_path)
            vsm_set = VortexStepMethod.VSMSettings(
                vsm_set_path; data_prefix=false)
        end
    end

    # Validate all wings are the same concrete type
    # and narrow from AbstractWing[] to concrete type
    if length(wings) > 0
        W = typeof(wings[1])
        for i in 2:length(wings)
            @assert typeof(wings[i]) === W (
                "All wings must be the same concrete " *
                "type, got $(typeof(wings[i])) at " *
                "index $i, expected $W")
        end
        if eltype(wings) !== W
            wings = convert(Vector{W}, wings)
        end
    end

    # Expand Route 2 tethers (auto-generate points/segments)
    expand_auto_tethers!(points, segments, tethers, set)

    # Assign indices and resolve all references
    # This converts symbolic names to numeric indices
    (point_names_dict, twist_surface_names_dict,
     segment_names_dict, pulley_names_dict,
     tether_names_dict, winch_names_dict,
     wing_names_dict, transform_names_dict) =
        assign_indices_and_resolve!(
            points, twist_surfaces, segments, pulleys,
            tethers, winches, wings, transforms)

    # If no wings defined, convert WING points to STATIC
    if isempty(wings)
        wing_point_idxs = findall(point -> point.type == WING, points)
        if !isempty(wing_point_idxs)
            @warn "No wings provided but " *
                  "$(length(wing_point_idxs)) WING type " *
                  "points found. Converting to STATIC."
            for idx in wing_point_idxs
                points[idx] = Point(
                    points[idx].name,
                    points[idx].pos_cad,
                    STATIC;
                    extra_mass = points[idx].extra_mass,
                    body_frame_damping =
                        points[idx].body_frame_damping,
                    world_frame_damping =
                        points[idx].world_frame_damping,
                    transform = points[idx].transform_ref
                )
                points[idx].idx = idx  # Reassign idx after recreation
                points[idx].pos_w .= points[idx].pos_cad
                points[idx].vel_w .= 0.0
            end
        end
    end

    # Validate indices (now assigned by assign_indices_and_resolve!)
    for (i, point) in enumerate(points)
        @assert point.idx == i "Point $(point.idx) != $i"
        # Allow transform_idx=0 (no transform) or valid index
        @assert point.transform_idx == 0 ||
                point.transform_idx <= length(transforms)
    end
    for (i, twist_surface) in enumerate(twist_surfaces)
        @assert twist_surface.idx == i
    end
    for (i, segment) in enumerate(segments)
        @assert segment.idx == i
        (segment.l0 ≈ 0) && (segment.l0 = segment_cad_length(segment, points))
    end
    for (i, pulley) in enumerate(pulleys)
        @assert pulley.idx == i
    end
    for (i, tether) in enumerate(tethers)
        @assert tether.idx == i
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
            if point.type == WING && point.wing_idx == wing.idx]
        point_mass_sum = sum(
            point.extra_mass for point in points
            if point.type == WING && point.wing_idx == wing.idx; init=0.0)
        set_mass = hasproperty(set, :mass) ? set.mass : 0.0

        if set_mass > 0 && point_mass_sum > 0
            @warn "Both set.mass ($set_mass) and wing point masses " *
                  "($point_mass_sum) specified for wing $(wing.idx). " *
                  "Using wing point masses (sys_struct priority)."
            wing.mass = point_mass_sum
        elseif point_mass_sum > 0
            wing.mass = point_mass_sum
        elseif set_mass > 0
            nwing_points = length(wing_point_idxs)
            if nwing_points > 0
                mass_per_point = set_mass / nwing_points
                for point_idx in wing_point_idxs
                    points[point_idx].extra_mass = mass_per_point
                end
            end
            wing.mass = set_mass
        else
            wing.mass = 0.0
        end
    end

    # Body frame + COM/principal inertia (dynamics, independent of aero mode).
    for wing in wings
        setup_wing_frame!(wing, points; prn)
    end

    # Per-mode aero construction (dispatched; no-op for modes without one, e.g.
    # flat-plate). Each VSM wing's pipeline — panel transform to body frame,
    # auto-twist, section matching, twist-surface/point mappings — runs inside
    # setup_aero! (see aero_modes/common.jl + vsm_refine.jl).
    for (i, wing) in enumerate(wings)
        @assert wing.idx == i
        setup_aero!(wing.aero, wing, points, twist_surfaces; prn)
    end

    for (i, transform) in enumerate(transforms)
        @assert transform.idx == i
    end
    set.physical_model = name

    # Name dictionaries were already built by assign_indices_and_resolve!
    sys_struct = SystemStructure(name, set,
        NamedCollection{Point}(points, point_names_dict),
        NamedCollection{TwistSurface}(twist_surfaces, twist_surface_names_dict),
        NamedCollection{Segment}(segments, segment_names_dict),
        NamedCollection{Pulley}(pulleys, pulley_names_dict),
        NamedCollection{Tether}(tethers, tether_names_dict),
        NamedCollection{Winch}(winches, winch_names_dict),
        NamedCollection{eltype(wings)}(wings, wing_names_dict),
        NamedCollection{Transform}(transforms, transform_names_dict),
        AtmosphericModel(set), false, false, vsm_set)
    reinit!(sys_struct, set)

    # Recalculate segment rest lengths from current positions if requested
    if ignore_l0
        for segment in sys_struct.segments
            point1 = sys_struct.points[segment.point_idxs[1]]
            point2 = sys_struct.points[segment.point_idxs[2]]
            segment.l0 = norm(point2.pos_w - point1.pos_w)
        end
    end

    return sys_struct
end
