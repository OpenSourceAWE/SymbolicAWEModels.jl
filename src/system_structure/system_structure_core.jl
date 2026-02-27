# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

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
- [`Group`](@ref): Collections of points for wing deformation.
- [`Segment`](@ref): Spring-damper elements.
- [`Pulley`](@ref): Elements that redistribute line lengths.
- [`Tether`](@ref): Groups of segments controlled by a winch.
- [`Winch`](@ref): Ground-based winches.
- [`Wing`](@ref): Rigid wing bodies.
- [`Transform`](@ref): Spatial transformations for initial positioning.
"""
mutable struct SystemStructure{W<:AbstractWing}
    const name::String
    const set::Settings
    const points::NamedCollection{Point}
    const groups::NamedCollection{Group}
    const segments::NamedCollection{Segment}
    const pulleys::NamedCollection{Pulley}
    const tethers::NamedCollection{Tether}
    const winches::NamedCollection{Winch}
    const wings::NamedCollection{W}
    const transforms::NamedCollection{Transform}

    const y::Array{Float64, 2}
    const x::Array{Float64, 2}
    const jac::Array{Float64, 3}
    const wind_vec_gnd::KVec3
    const am::AtmosphericModel
    wind_elevation::SimFloat
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
        # wings (principal frame ODE state, QUATERNION only)
        wings = getfield(sys, :wings)
        for wing in wings
            wing.wing_type != QUATERNION && continue
            append!(vars, wing.com_w)
            append!(vars, wing.com_vel)
            append!(vars, wing.Q_p_to_w)
            append!(vars, wing.ω_p)
        end
        # groups
        groups = getfield(sys, :groups)
        for group in groups
            if group.type == DYNAMIC
                push!(vars, group.twist)
                push!(vars, group.twist_ω)
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
        # winches
        winches = getfield(sys, :winches)
        for winch in winches
            push!(vars, winch.tether_len)
            push!(vars, winch.tether_vel)
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
        # wings (principal frame ODE state, QUATERNION only)
        wings = getfield(sys, :wings)
        for wing in wings
            wing.wing_type != QUATERNION && continue
            wing.com_w .= @view flat_value[offset:offset+2]
            offset += 3
            wing.com_vel .= @view flat_value[offset:offset+2]
            offset += 3
            wing.Q_p_to_w .= @view flat_value[offset:offset+3]
            offset += 4
            wing.ω_p .= @view flat_value[offset:offset+2]
            offset += 3
        end
        # groups
        groups = getfield(sys, :groups)
        for group in groups
            if group.type == DYNAMIC
                group.twist = flat_value[offset]
                offset += 1
                group.twist_ω = flat_value[offset]
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
        # winches
        winches = getfield(sys, :winches)
        for winch in winches
            winch.tether_len = flat_value[offset]
            offset += 1
            winch.tether_vel = flat_value[offset]
            offset += 1
        end
        return value
    else
        return setfield!(sys, sym, value)
    end
end

"""
    calc_heading(sys::SystemStructure)

Calculate heading angles for all wings in the system structure.
Returns a vector of heading angles, one per wing.
"""
function calc_heading(sys::SystemStructure)
    wind_norm = normalize(getfield(sys, :wind_vec_gnd))
    return [calc_heading(wing.R_b_to_w, wind_norm) for wing in getfield(sys, :wings)]
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

function resolve_ref(ref::Nothing, name_dict::Dict{Symbol, Int64}, component_type::String)
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
    return Int64[resolve_ref(r, name_dict, component_type) for r in spec]
end

function resolve_ref_spec(spec::Nothing, name_dict::Dict{Symbol, Int64}, component_type::String)
    return nothing
end

"""
    assign_indices_and_resolve!(components, name_dicts)

Assign indices to all components based on their position in the vectors,
and resolve all references to indices.
"""
function assign_indices_and_resolve!(
    points::Vector{Point},
    groups::Vector{Group},
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
    for (i, group) in enumerate(groups)
        group.idx = i
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
    group_names = build_name_dict(groups)
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

    # Groups: resolve point_refs
    for group in groups
        group.point_idxs = Int64[resolve_ref(r, point_names, "point") for r in group.point_refs]
    end

    # Segments: resolve point_refs
    for segment in segments
        p1 = resolve_ref(segment.point_refs[1], point_names, "point")
        p2 = resolve_ref(segment.point_refs[2], point_names, "point")
        segment.point_idxs = (p1, p2)
    end

    # Pulleys: resolve segment_refs
    for pulley in pulleys
        s1 = resolve_ref(pulley.segment_refs[1], segment_names, "segment")
        s2 = resolve_ref(pulley.segment_refs[2], segment_names, "segment")
        pulley.segment_idxs = (s1, s2)
    end

    # Tethers: resolve segment_refs and winch_point_ref
    for tether in tethers
        tether.segment_idxs = Int64[resolve_ref(r, segment_names, "segment") for r in tether.segment_refs]
        tether.winch_point_idx = resolve_ref(tether.winch_point_ref, point_names, "point")
    end

    # Winches: resolve tether_refs
    for winch in winches
        winch.tether_idxs = Int64[resolve_ref(r, tether_names, "tether") for r in winch.tether_refs]
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

    # Wings: resolve group_refs, transform_ref, and REFINE-specific refs
    for wing in wings
        # BaseWing fields
        wing.group_idxs = Int64[resolve_ref(r, group_names, "group") for r in wing.group_refs]
        wing.transform_idx = resolve_ref(wing.transform_ref, transform_names, "transform")

        # VSMWing-specific REFINE fields
        if isa(wing, VSMWing)
            if !isnothing(wing.origin_ref)
                wing.origin_idx = resolve_ref(wing.origin_ref, point_names, "point")
            end
            if !isnothing(wing.z_ref_points_ref)
                z1 = resolve_ref_spec(wing.z_ref_points_ref[1], point_names, "point")
                z2 = resolve_ref_spec(wing.z_ref_points_ref[2], point_names, "point")
                wing.z_ref_points = (z1, z2)
            end
            if !isnothing(wing.y_ref_points_ref)
                y1 = resolve_ref_spec(wing.y_ref_points_ref[1], point_names, "point")
                y2 = resolve_ref_spec(wing.y_ref_points_ref[2], point_names, "point")
                wing.y_ref_points = (y1, y2)
            end
        end
    end

    return (point_names, group_names, segment_names, pulley_names,
            tether_names, winch_names, wing_names, transform_names)
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
    compute_spatial_group_mapping!(the_wing, groups, points)

Map groups to unrefined sections using spatial proximity.
Each group is assigned to the closest unrefined section
based on distance between centers.
"""
function compute_spatial_group_mapping!(
    the_wing::VSMWing,
    groups::AbstractVector{Group},
    points::AbstractVector{Point}
)
    the_vsm_wing = the_wing.vsm_wing
    n_unrefined = the_vsm_wing.n_unrefined_sections
    n_groups = length(the_wing.base.group_idxs)

    # Compute group centers in body frame
    group_centers = Vector{MVec3}(undef, n_groups)
    for (local_idx, group_idx) in
            enumerate(the_wing.base.group_idxs)
        group = groups[group_idx]
        center = zeros(3)
        for pt_idx in group.point_idxs
            center += the_wing.base.R_b_to_c' *
                (points[pt_idx].pos_cad -
                 the_wing.base.pos_cad)
        end
        group_centers[local_idx] =
            center / length(group.point_idxs)
    end

    # Compute unrefined section centers
    unrefined_centers = Vector{MVec3}(
        undef, n_unrefined)
    for i in 1:n_unrefined
        le_point =
            the_vsm_wing.unrefined_sections[i].LE_point
        te_point =
            the_vsm_wing.unrefined_sections[i].TE_point
        unrefined_centers[i] =
            (le_point + te_point) / 2
    end

    # Map each group to closest unrefined section
    for (local_idx, group_idx) in
            enumerate(the_wing.base.group_idxs)
        group = groups[group_idx]
        min_dist = Inf
        closest_idx = 1
        for unrefined_idx in 1:n_unrefined
            dist = norm(
                group_centers[local_idx] -
                unrefined_centers[unrefined_idx])
            if dist < min_dist
                min_dist = dist
                closest_idx = unrefined_idx
            end
        end
        group.unrefined_section_idxs =
            Int64[closest_idx]
    end

    # Validate: check all sections are covered
    assigned = Set{Int64}()
    for group_idx in the_wing.base.group_idxs
        union!(assigned,
            groups[group_idx].unrefined_section_idxs)
    end
    if length(assigned) != n_unrefined
        unassigned = setdiff(1:n_unrefined, assigned)
        @warn "Wing $(the_wing.base.idx): " *
            "$(length(unassigned)) unrefined sections " *
            "not assigned to any group: $unassigned"
    end
end

# ==================== CONSTRUCTOR ==================== #

"""
    SystemStructure(name, set; points, groups, segments, pulleys, tethers, winches, wings, transforms)

Constructs a `SystemStructure` object representing a complete kite system.

## Physical Models
- **"ram"**: A model with 4 deformable wing groups and a complex pulley bridle system.
- **"simple_ram"**: A model with 4 deformable wing groups and direct bridle connections.

# Arguments
- `name::String`: Model identifier ("ram", "simple_ram", or a custom name).
- `set::Settings`: Configuration parameters from `KiteUtils.jl`.

# Keyword Arguments
- `points`, `groups`, `segments`, etc.: Vectors of the system components.
- `prn::Bool=true`: If true, print info messages about auto-generated components.

# Returns
- `SystemStructure`: A complete system ready for building a `SymbolicAWEModel`.
"""
function SystemStructure(name, set;
        points=Point[],
        groups=Group[],
        segments=Segment[],
        pulleys=Pulley[],
        tethers=Tether[],
        winches=Winch[],
        wings=VSMWing[],
        transforms=Transform[],
        ignore_l0::Bool=false,
        vsm_set=nothing,
        prn::Bool=true,
    )
    # Load VSMSettings if not provided and wings exist
    if isnothing(vsm_set) && !isempty(wings)
        model_dir = get_data_path()
        vsm_set_path = joinpath(model_dir, "vsm_settings.yaml")
        if isfile(vsm_set_path)
            vsm_set = VortexStepMethod.VSMSettings(vsm_set_path; data_prefix=false)
        end
    end

    # Validate all wings are the same concrete type
    if length(wings) > 1
        W = typeof(wings[1])
        for i in 2:length(wings)
            @assert typeof(wings[i]) === W (
                "All wings must be the same concrete " *
                "type, got $(typeof(wings[i])) at " *
                "index $i, expected $W")
        end
    end

    # Assign indices and resolve all references FIRST
    # This converts symbolic names to numeric indices
    (point_names_dict, group_names_dict, segment_names_dict, pulley_names_dict,
     tether_names_dict, winch_names_dict, wing_names_dict, transform_names_dict) =
        assign_indices_and_resolve!(points, groups, segments, pulleys, tethers, winches, wings, transforms)

    # If no wings defined, convert WING points to STATIC
    if isempty(wings)
        wing_point_idxs = findall(p -> p.type == WING, points)
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
    for (i, group) in enumerate(groups)
        @assert group.idx == i
    end
    for (i, segment) in enumerate(segments)
        @assert segment.idx == i
        (segment.l0 ≈ 0) && (segment.l0 = norm(points[segment.point_idxs[1]].pos_cad - points[segment.point_idxs[2]].pos_cad))
    end
    for (i, pulley) in enumerate(pulleys)
        @assert pulley.idx == i
    end
    for (i, tether) in enumerate(tethers)
        @assert tether.idx == i
    end
    for (i, winch) in enumerate(winches)
        @assert winch.idx == i
        if iszero(winch.tether_len)
            for segment_idx in tethers[winch.tether_idxs[1]].segment_idxs
                winch.tether_len += segments[segment_idx].l0
            end
        end
    end
    # Compute body frame (COM + principal axes) and
    # transform VSM panels from CAD → body frame.
    # QUATERNION: COM from point masses, Y-axis rotation
    #   to diagonalize inertia tensor.
    # REFINE: origin from origin_idx, R_b_to_c from
    #   z/y_ref_points (no inertia needed).
    for wing in wings
        isa(wing, VSMWing) || continue
        vsm_wing = wing.vsm_wing

        if wing.wing_type == QUATERNION
            wing_pts = [p for p in points
                if p.type == WING &&
                   p.wing_idx == wing.idx]
            isempty(wing_pts) && continue

            masses = [p.extra_mass for p in wing_pts]
            total_m = sum(masses)
            com_cad = if total_m > 0
                sum(masses[j] .* wing_pts[j].pos_cad
                    for j in eachindex(wing_pts)) /
                    total_m
            else
                mean([p.pos_cad for p in wing_pts])
            end

            # Inertia tensor about COM in CAD frame
            if total_m > 0
                I_cad = zeros(3, 3)
                for (m, p) in zip(masses, wing_pts)
                    r = p.pos_cad - com_cad
                    I_cad += m * (dot(r, r) * I(3) -
                                  r * r')
                end
                I_diag, Ry =
                    calc_inertia_y_rotation(I_cad)
                wing.R_p_to_c .= Ry'  # principal→CAD
                wing.inertia_principal .= diag(I_diag)
            end

            # Compute body frame from ref points
            if !isnothing(wing.origin_idx) &&
               !isnothing(wing.z_ref_points) &&
               !isnothing(wing.y_ref_points)
                origin_cad =
                    points[wing.origin_idx].pos_cad
                wing.pos_cad .= origin_cad

                # Temporarily set pos_w = pos_cad
                for p in points
                    p.type == WING &&
                        p.wing_idx == wing.idx &&
                        (p.pos_w .= p.pos_cad)
                end
                R_b_to_c, _ = calc_refine_wing_frame(
                    points, wing.z_ref_points,
                    wing.y_ref_points,
                    wing.origin_idx)
                wing.R_b_to_c .= R_b_to_c

                # COM offset from body origin in body
                wing.com_offset_b .=
                    R_b_to_c' * (com_cad - origin_cad)
            else
                # No ref points: body = principal,
                # origin = COM
                wing.pos_cad .= com_cad
                wing.R_b_to_c .= wing.R_p_to_c
                wing.com_offset_b .= 0.0
            end

            # Transform VSM sections: CAD → body
            vsm_wing.T_cad_body .= wing.pos_cad
            adjust_vsm_panels_to_origin!(
                vsm_wing, wing.pos_cad)
            rotate_vsm_sections!(
                vsm_wing, wing.R_b_to_c')
            vsm_wing.R_cad_body .= wing.R_b_to_c
            apply_aero_z_offset!(
                vsm_wing, wing.aero_z_offset)
            VortexStepMethod.reinit!(wing.vsm_aero)

            # Body → principal (constant rotation)
            wing.R_b_to_p .= wing.R_p_to_c' * wing.R_b_to_c

            if prn
                I_rnd = round.(wing.inertia_principal;
                               digits=4)
                off = round.(wing.com_offset_b;
                             digits=4)
                @info "QUATERNION wing $(wing.idx):" *
                    " COM=[$(round.(com_cad; digits=3))]" *
                    ", I=$I_rnd" *
                    ", com_offset_b=$off"
            end

        elseif wing.wing_type == REFINE
            # Body frame from structural ref points.
            # Points are in CAD frame at construction
            # time (pos_w = pos_cad before transforms).
            if !isnothing(wing.origin_idx) &&
               !isnothing(wing.z_ref_points) &&
               !isnothing(wing.y_ref_points)
                origin_pos = points[
                    wing.origin_idx].pos_cad
                wing.pos_cad .= origin_pos

                # Temporarily set pos_w = pos_cad so
                # calc_refine_wing_frame can read them
                for p in points
                    p.type == WING &&
                        p.wing_idx == wing.idx &&
                        (p.pos_w .= p.pos_cad)
                end
                R_b_to_c, _ = calc_refine_wing_frame(
                    points, wing.z_ref_points,
                    wing.y_ref_points,
                    wing.origin_idx)
                wing.R_b_to_c .= R_b_to_c

                # Transform VSM sections: CAD → body
                vsm_wing.T_cad_body .= origin_pos
                adjust_vsm_panels_to_origin!(
                    vsm_wing, origin_pos)
                rotate_vsm_sections!(
                    vsm_wing, wing.R_b_to_c')
                vsm_wing.R_cad_body .= wing.R_b_to_c
                VortexStepMethod.reinit!(wing.vsm_aero)

                if prn
                    o = round.(origin_pos; digits=3)
                    @info "REFINE wing " *
                        "$(wing.idx): origin=[$o]"
                end
            end
        end
    end

    # Auto-create groups for QUATERNION wings if needed (before geometry initialization)
    # Skip for AERO_NONE — no aerodynamics means no twist DOFs needed.
    for (i, wing) in enumerate(wings)
        if wing.wing_type == QUATERNION &&
           isempty(wing.group_idxs) &&
           wing.aero_mode != AERO_NONE
            # Get WING-type points for this wing
            wing_point_idxs = findall(
                p -> p.type == WING && p.wing_idx == wing.idx, points)
            wing_points = [points[idx] for idx in wing_point_idxs]

            # Identify LE/TE pairs
            wing_segments = identify_wing_segments(wing_points)

            # Create a group for each section (LE/TE pair)
            # n_groups = n_unrefined_sections (one group per section)
            new_group_idxs = Int64[]

            for (le_idx, te_idx) in wing_segments
                group_idx = length(groups) + 1
                # Use integer as name for auto-created groups
                group_name = group_idx

                # Both LE and TE points (matches YAML convention)
                new_group = Group(group_name,
                    [le_idx, te_idx], DYNAMIC, 0.0)

                # Assign idx and resolve point_refs since
                # these are dynamically created
                new_group.idx = group_idx
                new_group.point_idxs = [le_idx, te_idx]

                push!(groups, new_group)
                push!(new_group_idxs, Int64(group_idx))
            end

            # Update wing with new groups and resize vsm arrays
            wing.group_idxs = new_group_idxs

            # Resize vsm arrays based on number of unrefined sections
            n_unrefined = wing.vsm_wing.n_unrefined_sections
            ny = 3 + n_unrefined + 3  # va(3) + twist(n_unrefined) + ω(3)
            nx = 3 + 3 + n_unrefined  # force(3) + moment(3) + unrefined_moments(n_unrefined)
            wing.vsm_y = zeros(SimFloat, ny)
            wing.vsm_x = zeros(SimFloat, nx)
            wing.vsm_jac = zeros(SimFloat, nx, ny)

            prn && @info "Auto-created $(length(new_group_idxs)) groups " *
                  "for QUATERNION wing $(wing.idx)"
        end
    end

    # Match aero sections to structural LE/TE for ALL
    # VSMWing types (runs after auto-group creation so
    # identify_wing_segments can use groups).
    for wing in wings
        isa(wing, VSMWing) || continue
        wing.aero_mode == AERO_NONE && continue
        match_aero_sections_to_structure!(
            wing, points; groups=groups)
    end

    # Clear REFINE wing.group_idxs — groups were used
    # for LE/TE identification but REFINE doesn't use
    # them for aerodynamics.  Groups stay in sys_struct
    # (useful for structural info / future linearization).
    for wing in wings
        if wing.wing_type == REFINE &&
           !isempty(wing.group_idxs)
            empty!(wing.group_idxs)
        end
    end

    # Initialize group-to-unrefined-section mapping for QUATERNION wings
    # Do this BEFORE y_airf calculation so the mapping is available
    for the_wing in wings
        if isa(the_wing, VSMWing) && the_wing.base.wing_type == QUATERNION && !isempty(the_wing.base.group_idxs)
            compute_spatial_group_mapping!(the_wing, groups, points)
        end
    end

    # Initialize group geometries from closest VSM panel
    for group in groups
        if iszero(group.chord)
            # Find which wing this group belongs to
            for wing in wings
                if group.idx in wing.group_idxs
                    # Compute group center in body frame (average of all attach points)
                    center = zeros(3)
                    for pt_idx in group.point_idxs
                        center += wing.R_b_to_c' * (points[pt_idx].pos_cad - wing.pos_cad)
                    end
                    center ./= length(group.point_idxs)

                    # Find closest panel (panels are in body
                    # frame after reinit!)
                    panels = wing.vsm_aero.panels
                    min_dist = Inf
                    closest_panel = panels[1]
                    for panel in panels
                        pc = (panel.LE_point_1 +
                              panel.LE_point_2 +
                              panel.TE_point_1 +
                              panel.TE_point_2) / 4
                        dist = norm(center - pc)
                        if dist < min_dist
                            min_dist = dist
                            closest_panel = panel
                        end
                    end

                    # Panel geometry already in body frame
                    group.le_pos .=
                        (closest_panel.LE_point_1 +
                         closest_panel.LE_point_2) / 2
                    group.chord .=
                        closest_panel.x_airf *
                        closest_panel.chord
                    group.y_airf .=
                        closest_panel.y_airf

                    break
                end
            end
        end
    end

    # Translate group le_pos from body origin to COM
    # (body frame). chord and y_airf are direction
    # vectors already in body frame from VSM panels.
    for wing in wings
        wing.wing_type != QUATERNION && continue
        for group_idx in wing.group_idxs
            group = groups[group_idx]
            group.le_pos .-= wing.com_offset_b
        end
    end

    for (i, wing) in enumerate(wings)
        @assert wing.idx == i
        # For REFINE wings, set defaults if not provided
        if wing.wing_type == REFINE
            # Build point_to_vsm_point mapping if not provided
            if isnothing(wing.point_to_vsm_point)
                # Get WING-type points for this wing
                wing_point_idxs = findall(
                    p -> p.type == WING && p.wing_idx == wing.idx, points)
                wing_points = [points[idx]
                    for idx in wing_point_idxs]
                wing.point_to_vsm_point =
                    build_point_to_vsm_point_mapping(
                        wing_points, wing.vsm_wing)
            end

            wing_point_idxs = collect(keys(
                wing.point_to_vsm_point))
            wing_points = [points[idx]
                for idx in wing_point_idxs]

            # For REFINE wings, pos_cad should be user-specified (KCU position)
            # or default to vsm_wing.T_cad_body (set in VSMWing constructor)
            # DO NOT calculate as centroid - that would misalign VSM panels

            # Identify wing segments (LE/TE pairs)
            if isnothing(wing.wing_segments)
                wing.wing_segments =
                    identify_wing_segments(
                        wing_points; groups=groups,
                        wing_group_idxs=wing.group_idxs)
            end

            # Set default reference points if not provided
            if isnothing(wing.z_ref_points) ||
               isnothing(wing.y_ref_points)
                segs = wing.wing_segments

                if isnothing(wing.z_ref_points)
                    # Use first segment (center LE-TE) for Z (normal)
                    wing.z_ref_points = segs[1]
                end

                if isnothing(wing.y_ref_points)
                    # Use center LE and mid-span LE for Y (spanwise)
                    mid = length(segs) ÷ 2
                    wing.y_ref_points = (segs[1][1],
                        segs[mid][1])
                end
            end

        end
    end

    # Calculate wing mass with either/or priority logic:
    # - If wing points have extra_mass specified, use point masses (priority)
    # - If only set.mass specified, distribute to wing points
    # - Warn if both sources have nonzero values
    for wing in wings
        wing_point_idxs = [p.idx for p in points if p.type == WING && p.wing_idx == wing.idx]

        # Sum of user-specified WING point masses
        point_mass_sum = sum(
            p.extra_mass for p in points if p.idx in wing_point_idxs;
            init=0.0
        )

        # Check for conflict between set.mass and point masses
        set_mass = hasproperty(set, :mass) ? set.mass : 0.0

        if set_mass > 0 && point_mass_sum > 0
            @warn "Both set.mass ($set_mass) and wing point masses ($point_mass_sum) " *
                  "specified for wing $(wing.idx). Using wing point masses (sys_struct " *
                  "priority)."
            wing.mass = point_mass_sum
        elseif point_mass_sum > 0
            # Wing point masses specified, no set.mass
            wing.mass = point_mass_sum
        elseif set_mass > 0
            # set.mass specified, distribute equally to wing points
            n_wing_points = length(wing_point_idxs)
            if n_wing_points > 0
                mass_per_point = set_mass / n_wing_points
                for point_idx in wing_point_idxs
                    points[point_idx].extra_mass = mass_per_point  # ASSIGN, not add
                end
            end
            wing.mass = set_mass
        else
            # Neither specified - wing has no mass
            wing.mass = 0.0
        end

        # Mass and inertia validation is done in validate_sys_struct()
    end

    for (i, transform) in enumerate(transforms)
        @assert transform.idx == i

        # Check for conflict with Settings (only warn if values differ)
        set_elev = hasproperty(set, :elevations) && i <= length(set.elevations) ?
                   set.elevations[i] : 0.0
        set_azim = hasproperty(set, :azimuths) && i <= length(set.azimuths) ?
                   set.azimuths[i] : 0.0
        set_head = hasproperty(set, :headings) && i <= length(set.headings) ?
                   set.headings[i] : 0.0

        sys_elev = rad2deg(transform.elevation)
        sys_azim = rad2deg(transform.azimuth)
        sys_head = rad2deg(transform.heading)

        # Only warn if both have nonzero values AND they differ
        elev_conflict = (set_elev != 0.0 && sys_elev != 0.0 && !isapprox(set_elev, sys_elev))
        azim_conflict = (set_azim != 0.0 && sys_azim != 0.0 && !isapprox(set_azim, sys_azim))
        head_conflict = (set_head != 0.0 && sys_head != 0.0 && !isapprox(set_head, sys_head))

        if elev_conflict || azim_conflict || head_conflict
            @warn "Transform $(transform.name): Settings and sys_struct have different " *
                  "angles. Settings: (elev=$(set_elev)°, azim=$(set_azim)°, " *
                  "head=$(set_head)°). Using sys_struct: (elev=$(sys_elev)°, " *
                  "azim=$(sys_azim)°, head=$(sys_head)°)."
        end

        # sys_struct values take priority - update Settings to match
        set.elevations[i] = sys_elev
        set.azimuths[i]   = sys_azim
        set.headings[i]   = sys_head
    end
    if length(wings) > 0
        # Use number of unrefined sections
        n_unrefined = wings[1].vsm_wing.n_unrefined_sections
        ny = 3 + n_unrefined + 3
        nx = 3 + 3 + n_unrefined
    else
        ny = 0
        nx = 0
    end
    y = zeros(length(wings), ny)
    x = zeros(length(wings), nx)
    jac = zeros(length(wings), nx, ny)
    set.physical_model = name

    # Name dictionaries were already built by assign_indices_and_resolve!
    sys_struct = SystemStructure(name, set,
        NamedCollection{Point}(points, point_names_dict),
        NamedCollection{Group}(groups, group_names_dict),
        NamedCollection{Segment}(segments, segment_names_dict),
        NamedCollection{Pulley}(pulleys, pulley_names_dict),
        NamedCollection{Tether}(tethers, tether_names_dict),
        NamedCollection{Winch}(winches, winch_names_dict),
        NamedCollection{eltype(wings)}(wings, wing_names_dict),
        NamedCollection{Transform}(transforms, transform_names_dict),
        y, x, jac, zeros(KVec3), AtmosphericModel(set), 0.0, false, false, vsm_set)
    reinit!(sys_struct, set)

    # Recalculate segment rest lengths from current positions if requested
    if ignore_l0
        for segment in sys_struct.segments
            p1 = sys_struct.points[segment.point_idxs[1]]
            p2 = sys_struct.points[segment.point_idxs[2]]
            segment.l0 = norm(p2.pos_w - p1.pos_w)
        end
    end

    return sys_struct
end
