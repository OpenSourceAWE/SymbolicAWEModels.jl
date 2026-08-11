# Copyright (c) 2026 Uwe Fechner, Bart van de Lint
# SPDX-License-Identifier: MIT

"""
    yaml_float_cell(value) -> Any

A float ready for a YAML cell: `nothing` when it is `NaN`, so the loader's
"absent means unset" rule reproduces it, and the value itself otherwise.
"""
yaml_float_cell(value) = isnan(value) ? nothing : SimFloat(value)

"""
    yaml_cell(value) -> String

Render one YAML cell. Floats print at full round-trip precision because
[`save_sys_struct_to_yaml`](@ref) is expected to reproduce an ODE state exactly,
so no value may be rounded on the way out.
"""
function yaml_cell(value)
    isnothing(value) && return "nothing"
    value isa Bool && return value ? "true" : "false"
    value isa Symbol && return String(value)
    value isa AbstractString && return String(value)
    value isa Integer && return string(value)
    if value isa Real
        isnan(value) && return ".nan"
        isinf(value) && return value > 0 ? ".inf" : "-.inf"
        return repr(SimFloat(value))
    end
    value isa Tuple && return yaml_cell(collect(value))
    if value isa AbstractVector
        return "[" * join(yaml_cell.(value), ", ") * "]"
    end
    error("save_sys_struct_to_yaml: cannot write a $(typeof(value)) to YAML")
end

"""
    yaml_ref(idx, collection) -> Any

The name of `collection[idx]`, falling back to the raw index when the component
is unnamed and to `nothing` when `idx` does not address one. Reference columns
are written by name so a saved structure stays readable and stays valid when
rows are reordered.
"""
function yaml_ref(idx, collection)
    (isnothing(idx) || idx <= 0 || idx > length(collection)) && return nothing
    name = collection[idx].name
    return isnothing(name) ? idx : name
end

"""
    dynamics_type_string(type::DynamicsType) -> String

The YAML spelling [`parse_dynamics_type`](@ref) reads back.
"""
function dynamics_type_string(type::DynamicsType)
    type == STATIC && return "STATIC"
    type == DYNAMIC && return "DYNAMIC"
    type == BODY_STATIC && return "BODY_STATIC"
    type == KINEMATIC && return "KINEMATIC"
    error("Unknown DynamicsType: $type")
end

"""
    wing_type_string(type::WingType) -> String

The YAML spelling [`parse_wing_type`](@ref) reads back.
"""
function wing_type_string(type::WingType)
    type == PARTICLE_DYNAMICS && return "PARTICLE_DYNAMICS"
    type == RIGID_DYNAMICS && return "RIGID_DYNAMICS"
    error("Unknown WingType: $type")
end

"""
    principal_frame_method_string(method::PrincipalFrameMethod) -> String

The YAML spelling [`parse_principal_frame_method`](@ref) reads back.
"""
function principal_frame_method_string(method::PrincipalFrameMethod)
    method == EIGEN_DECOMP && return "EIGEN_DECOMP"
    method == Y_ROTATION && return "Y_ROTATION"
    error("Unknown PrincipalFrameMethod: $method")
end

"""
    aero_mode_string(aero::AbstractAeroModel) -> String

The YAML spelling [`parse_aero_mode`](@ref) reads back.
"""
function aero_mode_string(aero::AbstractAeroModel)
    aero isa AeroNone && return "AeroNone"
    aero isa AeroDirect && return "AeroDirect"
    aero isa AeroLinearized && return "AeroLinearized"
    aero isa AeroPlate && return "AeroPlate"
    aero isa ContinuousAero && return "ContinuousAero"
    aero isa AeroPressure && return "AeroPressure"
    error("Unknown aero model: $(typeof(aero))")
end

"""
    ref_points_cell(ref_points) -> Any

A `(z|y)_ref_points` pair, or a single `origin`, as the loader's reference spec:
a bare name for one point, a plain list when the weights are the equal-weight
average the loader assumes for one, and `[[name, weight], ...]` otherwise.
"""
function ref_points_cell(ref_points)
    isnothing(ref_points) && return nothing
    ref_points isa Tuple && length(ref_points) == 2 &&
        return [ref_points_cell(ref_points[1]), ref_points_cell(ref_points[2])]
    ref_points isa WeightedRefPoints || return ref_points
    names = isempty(ref_points.refs) ? collect(ref_points.ids) :
        collect(ref_points.refs)
    length(names) == 1 && return names[1]
    count = length(names)
    all(weight -> abs(weight - 1 / count) < 1e-10, ref_points.weights) &&
        return names
    return [[names[k], ref_points.weights[k]] for k in eachindex(names)]
end

"""
    write_yaml_table(io, key, headers, rows; title=nothing)

Write one `key: {headers, data}` block, padding each column so the table stays
readable. Nothing is written for an empty `rows`.

A column every row leaves unset is dropped rather than written as a column of
`nothing`: the loader treats an absent column and an unset one identically, and
some reference columns reject an explicit `nothing` that they accept as absent.
"""
function write_yaml_table(io, key, headers, rows; title = nothing)
    isempty(rows) && return nothing
    cells = [[yaml_cell(value) for value in row] for row in rows]
    keep = [any(row[k] != "nothing" for row in cells) for k in eachindex(headers)]
    headers = headers[keep]
    cells = [row[keep] for row in cells]
    widths = [maximum(length(row[k]) for row in cells) for k in eachindex(headers)]
    isnothing(title) || write_yaml_title(io, title)
    println(io, key, ":")
    println(io, "  headers: [", join(String.(headers), ", "), "]")
    println(io, "  data:")
    for row in cells
        padded = [rpad(row[k], widths[k]) for k in eachindex(row)]
        println(io, "    - [", rstrip(join(padded, ", ")), "]")
    end
    println(io)
    return nothing
end

"""
    write_yaml_mappings(io, key, rows; title=nothing)

Write one `key: {data}` block whose rows are mappings rather than a fixed header
row. Used for wings, whose optional fields are too sparse for a table. Entries
whose value is `nothing` are left out.
"""
function write_yaml_mappings(io, key, rows; title = nothing)
    isempty(rows) && return nothing
    isnothing(title) || write_yaml_title(io, title)
    println(io, key, ":")
    println(io, "  data:")
    for row in rows
        first_field = true
        for (field, value) in row
            isnothing(value) && continue
            prefix = first_field ? "    - " : "      "
            println(io, prefix, field, ": ", yaml_cell(value))
            first_field = false
        end
    end
    println(io)
    return nothing
end

"""
    write_yaml_title(io, title)

Write the banner comment that separates blocks in a structural YAML.
"""
function write_yaml_title(io, title)
    bar = "#"^30
    println(io, bar)
    println(io, "## ", rpad(title * " ", 25, "#"))
    println(io, bar)
    return nothing
end

"""
    point_rows(sys) -> (headers, rows)

The `points` block, in the CAD frame the loader expects. `vel_w` travels with it
because the loader zeroes velocities when the column is absent.
"""
function point_rows(sys)
    headers = [:name, :pos_cad, :type, :wing_idx, :transform_idx, :body_idx,
               :joint, :vel_w, :extra_mass, :body_frame_damping,
               :world_frame_damping, :area, :drag_coeff, :fix_sphere,
               :fix_static]
    rows = Vector{Any}[]
    for point in sys.points
        push!(rows, Any[
            point.name, collect(point.pos_cad),
            dynamics_type_string(point.type),
            yaml_ref(point.wing_idx, sys.wings),
            yaml_ref(point.transform_idx, sys.transforms),
            yaml_ref(point.body_idx, sys.bodies),
            yaml_ref(point.joint_idx, sys.timoshenko_joints),
            collect(point.vel_w), point.extra_mass,
            collect(point.body_frame_damping),
            collect(point.world_frame_damping), point.area, point.drag_coeff,
            point.fix_sphere, point.fix_static])
    end
    return headers, rows
end

"""
    segment_rows(sys) -> (headers, rows)

The `segments` block, written with the resolved `unit_stiffness`/`unit_damping`
rather than a material so no modulus has to be re-derived on load. Errors on a
callable force law, which has no YAML spelling.
"""
function segment_rows(sys)
    headers = [:name, :point_i, :point_j, :l0, :diameter_mm, :unit_stiffness,
               :unit_damping, :compression_frac, :density]
    rows = Vector{Any}[]
    for segment in sys.segments
        segment.unit_stiffness isa Real || error("save_sys_struct_to_yaml: " *
            "segment $(segment.name) has a callable unit_stiffness, which " *
            "cannot be written to YAML; supply it programmatically instead.")
        push!(rows, Any[
            segment.name, yaml_ref(segment.point_idxs[1], sys.points),
            yaml_ref(segment.point_idxs[2], sys.points), segment.l0,
            yaml_float_cell(segment.diameter * 1000),
            segment.unit_stiffness, segment.unit_damping,
            segment.compression_frac, yaml_float_cell(segment.density)])
    end
    return headers, rows
end

"""
    pulley_rows(sys) -> (headers, rows)

The `pulleys` block. `sum_len`, `len` and `vel` are written because the pair of
rest lengths a pulley owns is not recoverable from the segments alone.
"""
function pulley_rows(sys)
    headers = [:name, :segment_i, :segment_j, :type, :efficiency, :damping,
               :brake, :friction_epsilon, :sum_len, :len, :vel]
    rows = Vector{Any}[]
    for pulley in sys.pulleys
        push!(rows, Any[
            pulley.name, yaml_ref(pulley.segment_idxs[1], sys.segments),
            yaml_ref(pulley.segment_idxs[2], sys.segments),
            dynamics_type_string(pulley.type), pulley.efficiency,
            pulley.damping, pulley.brake, pulley.friction_epsilon,
            pulley.sum_len, pulley.len, pulley.vel])
    end
    return headers, rows
end

"""
    twist_surface_rows(sys) -> (headers, rows)

The `twist_surfaces` block, carrying the `twist`/`twist_ω` ODE state.
"""
function twist_surface_rows(sys)
    headers = [:name, :points, :type, :moment_frac, :damping, :stiffness,
               :wing, :bodies, :flap_bodies, :flap_axis, :twist, :twist_vel]
    rows = Vector{Any}[]
    for surface in sys.twist_surfaces
        push!(rows, Any[
            surface.name,
            [yaml_ref(idx, sys.points) for idx in surface.point_idxs],
            dynamics_type_string(surface.type), surface.moment_frac,
            surface.damping, surface.stiffness,
            yaml_ref(surface.wing_idx, sys.wings),
            [yaml_ref(idx, sys.bodies) for idx in surface.body_idxs],
            [yaml_ref(idx, sys.bodies) for idx in surface.flap_body_idxs],
            collect(surface.flap_axis), surface.twist, surface.twist_ω])
    end
    return headers, rows
end

"""
    tether_rows(sys) -> (headers, rows)

The `tethers` block, written on the explicit-segments route so the segments a
tether owns are the ones already in the `segments` block. `len` carries the
unstretched-length ODE state.
"""
function tether_rows(sys)
    headers = [:name, :segment_idxs, :start_point, :end_point,
               :init_stretched_length, :len]
    rows = Vector{Any}[]
    for tether in sys.tethers
        push!(rows, Any[
            tether.name,
            [yaml_ref(idx, sys.segments) for idx in tether.segment_idxs],
            yaml_ref(tether.start_point_idx, sys.points),
            yaml_ref(tether.end_point_idx, sys.points),
            tether.stretched_len, tether.len])
    end
    return headers, rows
end

"""
    winch_rows(sys) -> (headers, rows)

The `winches` block, carrying the drum-velocity ODE state and the live
`set_value`.
"""
function winch_rows(sys)
    headers = [:name, :tether_idxs, :winch_point, :init_vel, :brake,
               :speed_controlled, :vel, :set_value]
    rows = Vector{Any}[]
    for winch in sys.winches
        push!(rows, Any[
            winch.name,
            [yaml_ref(idx, sys.tethers) for idx in winch.tether_idxs],
            yaml_ref(winch.winch_point_idx, sys.points), winch.init_vel,
            winch.brake, winch.speed_controlled, winch.vel, winch.set_value])
    end
    return headers, rows
end

"""
    body_mappings(sys; wings) -> Vector

The `bodies` or `wings` rows as mappings. `wings=true` selects the aero-carrying
bodies and adds the wing-only fields; `wings=false` selects the plain ones.
Positions and orientation come from the live `pos_w`/`Q_b_to_w` so the rigid
state round-trips.
"""
function body_mappings(sys; wings::Bool)
    rows = Vector{Pair{Symbol, Any}}[]
    for body in sys.bodies
        is_wing(body) == wings || continue
        state = Pair{Symbol, Any}[
            :vel => collect(body.vel_w),
            :Q_b_to_w => collect(body.Q_b_to_w),
            :omega_b => collect(body.ω_b),
            :principal_frame_method =>
                principal_frame_method_string(body.principal_frame_method)]
        transform = yaml_ref(body.transform_idx, sys.transforms)
        origin = ref_points_cell(body.origin)
        row = if wings
            Pair{Symbol, Any}[
                :name => body.name,
                :transform_idx => transform,
                :dynamics_type => wing_type_string(body.dynamics_type),
                :aero_mode => aero_mode_string(body.aero),
                :twist_surfaces => [yaml_ref(idx, sys.twist_surfaces)
                                    for idx in body.twist_surface_idxs],
                :z_ref_points => ref_points_cell(body.z_ref_points),
                :y_ref_points => ref_points_cell(body.y_ref_points),
                :origin_idx => origin,
                :pos_cad => isnothing(origin) ? collect(body.pos_cad) : nothing,
                :y_damping => body.damping[2] - body.damping[1],
                :angular_damping => body.damping[1],
                :principal_frame_method =>
                    principal_frame_method_string(body.principal_frame_method)]
        else
            Pair{Symbol, Any}[
                :name => body.name,
                :transform_idx => transform,
                :mass => body.mass,
                :pos => collect(body.pos_cad),
                :inertia_principal => collect(body.inertia_principal),
                :type => dynamics_type_string(body.type),
                :com_offset_b => collect(body.com_offset_b),
                :damping => collect(body.damping),
                :fix_sphere => body.fix_sphere,
                :ext_force_w => collect(body.ext_force_w),
                :ext_force_b => collect(body.ext_force_b),
                :ext_moment_b => collect(body.ext_moment_b),
                state...]
        end
        push!(rows, row)
    end
    return rows
end

"""
    joint_rows(sys, joints, scalar_fields) -> (headers, rows)

A `timoshenko_joints` or `elastic_joints` block: the two bodies, both anchors and
every `scalar_fields` stiffness/damping entry.
"""
function joint_rows(sys, joints, scalar_fields)
    headers = [:name, :body_a, :body_b, :anchor_a, :anchor_b, scalar_fields...]
    rows = Vector{Any}[]
    for joint in joints
        row = Any[joint.name, yaml_ref(joint.body_a_idx, sys.bodies),
                  yaml_ref(joint.body_b_idx, sys.bodies),
                  collect(joint.anchor_a_b), collect(joint.anchor_b_b)]
        for field in scalar_fields
            value = getfield(joint, field)
            push!(row, value isa Real ? yaml_float_cell(value) : value)
        end
        push!(rows, row)
    end
    return headers, rows
end

"""
    transform_rows(sys) -> (headers, rows)

The `transforms` block. Angles are written in degrees, matching the loader and
every other angle in a YAML file.
"""
function transform_rows(sys)
    headers = [:name, :elevation, :azimuth, :heading, :base_point_idx,
               :base_pos, :base_transform_idx, :wing_idx, :rot_point_idx,
               :elevation_vel, :azimuth_vel, :turn_rate]
    rows = Vector{Any}[]
    for transform in sys.transforms
        push!(rows, Any[
            transform.name, rad2deg(transform.elevation),
            rad2deg(transform.azimuth), rad2deg(transform.heading),
            yaml_ref(transform.base_point_idx, sys.points),
            isnothing(transform.base_pos) ? nothing :
                collect(transform.base_pos),
            yaml_ref(transform.base_transform_idx, sys.transforms),
            yaml_ref(transform.wing_idx, sys.wings),
            yaml_ref(transform.rot_point_idx, sys.points),
            rad2deg(transform.elevation_vel), rad2deg(transform.azimuth_vel),
            rad2deg(transform.turn_rate)])
    end
    return headers, rows
end

"""
    save_sys_struct_to_yaml(sys::SystemStructure, yaml_path; prn=true)

Write `sys` to `yaml_path` in the component-based format
[`load_sys_struct_from_yaml`](@ref) reads, and return the path.

Floats are written at full precision, so every parameter a model is built from —
rest lengths, stiffnesses, masses, pulley `sum_len`, damping, joint properties —
reloads bit for bit. The live non-positional state travels with them: pulley
`len`/`vel`, tether `len`, winch `vel`/`set_value`, twist `twist`/`twist_ω`,
point `vel_w` and body `vel`/`Q_b_to_w`/`ω_b`. A callable segment force law has
no YAML spelling and raises.

Positions are written as `pos_cad`, the CAD frame, not as the settled `pos_w`.
`pos_cad` is not merely an initial condition: `compute_spatial_twist_surface_mapping!`
matches structural geometry against the aero sections in that frame, and
`reinit!` re-derives `pos_w` from it by applying the transforms. Writing world
coordinates there would both place the structure twice and break the aero
mapping. A saved file therefore reproduces the geometry it was built from, and
settled positions have to be mapped back into the CAD frame by the caller.
"""
function save_sys_struct_to_yaml(sys::SystemStructure,
                                 yaml_path::AbstractString; prn::Bool = true)
    open(yaml_path, "w") do io
        println(io, "#"^30)
        println(io, "## System Structure ##########")
        println(io, "#"^30)
        println(io, "# Written by save_sys_struct_to_yaml from system ",
                repr(sys.name), ".")
        println(io, "# Positions and velocities are the live state, not the ",
                "design geometry.")
        println(io)

        write_yaml_table(io, "points", point_rows(sys)...; title = "Points")
        write_yaml_table(io, "segments", segment_rows(sys)...;
                         title = "Segments")
        write_yaml_table(io, "pulleys", pulley_rows(sys)...; title = "Pulleys")
        write_yaml_table(io, "twist_surfaces", twist_surface_rows(sys)...;
                         title = "Twist surfaces")
        write_yaml_table(io, "tethers", tether_rows(sys)...; title = "Tethers")
        write_yaml_table(io, "winches", winch_rows(sys)...; title = "Winches")
        write_yaml_mappings(io, "wings", body_mappings(sys; wings = true);
                            title = "Wings")
        write_yaml_mappings(io, "bodies", body_mappings(sys; wings = false);
                            title = "Bodies")
        write_yaml_table(io, "timoshenko_joints",
                         joint_rows(sys, sys.timoshenko_joints,
                                    (:EA, :GA, :GJ, :EIy, :EIz, :shear_coeff,
                                     :damping_trans, :damping_rot,
                                     :rest_length, :radius))...;
                         title = "Timoshenko joints")
        write_yaml_table(io, "elastic_joints",
                         joint_rows(sys, sys.elastic_joints,
                                    (:stiffness_axial, :stiffness_shear,
                                     :stiffness_torsion, :stiffness_bending,
                                     :damping_trans, :damping_rot,
                                     :radius))...;
                         title = "Elastic joints")
        write_yaml_table(io, "transforms", transform_rows(sys)...;
                         title = "Transforms")
    end
    prn && @info "System structure saved to $yaml_path"
    return yaml_path
end
