# Copyright (c) 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

"""
    system_definition(sys::SystemStructure) -> SystemDefinition

The KiteGeometry definition of `sys` in its initial pose: its components as rows, every
reference an index, a station only where it is a wing's, what only SymbolicAWEModels
reads in the models it registers ([`register_definition_models`](@ref)), and each body at
the centre of its own mass.
"""
function system_definition(sys::SystemStructure)
    materials = KiteGeometry.Material[]
    return SystemDefinition(;
        metadata=Metadata(; name=sys.name, n_points=length(sys.points)),
        points=[point_row(sys, point) for point in sys.points],
        segments=[segment_row!(materials, segment) for segment in sys.segments],
        spanwise_stations=[station_row(sys, station) for station in sys.stations
                           if !isnothing(station_wing_idx(sys, station))],
        pulleys=[pulley_row(pulley) for pulley in sys.pulleys],
        tethers=[tether_row(tether) for tether in sys.tethers],
        winches=[winch_row(winch) for winch in sys.winches],
        wings=[wing_row(sys, wing) for wing in sys.wings],
        bodies=[body_row(sys, body) for body in sys.bodies],
        tubes=[tube_row!(materials, sys, tube) for tube in sys.tubes],
        materials)
end

"""The name a definition gives component `component`: its own, else its index."""
definition_name(component) = Symbol(something(component.name, component.idx))

"""Index `idx` as a reference, `nothing` for the absent reference 0."""
optional_index(idx) = idx == 0 ? nothing : idx

"""The index of the material `law` with `parameters` among `materials`, added where none
matches."""
function material!(materials, law, parameters)
    index = findfirst(material -> material.law == law &&
                                  material.parameters == parameters, materials)
    isnothing(index) || return index
    push!(materials, KiteGeometry.Material(; name="$(law)_$(length(materials) + 1)", law,
                                           parameters))
    return length(materials)
end

function point_row(sys, point)
    body = optional_index(point.body_idx)
    return KiteGeometry.Point(; name=definition_name(point), fixed=point.type == STATIC,
        body, pos_ENU=point.pos_ENU, point.extra_mass,
        model=PointModel(string(point.type), Vector(point.body_frame_damping),
            Vector(point.world_frame_damping), point.area, point.drag_coeff,
            optional_index(point.tube_idx), point.fix_sphere, point.fix_static,
            optional_index(point.wing_idx),
            isnothing(body) ? nothing : Vector(point.anchor_b)))
end

function segment_row!(materials, segment)
    area = π * (segment.diameter / 2)^2
    stiffness = segment.unit_stiffness
    law, parameters = stiffness isa Real && (area > 0 || iszero(stiffness)) ?
        ("linear", OrderedDict{String, Any}("density" => segment.density,
            "youngs_modulus" => iszero(stiffness) ? 0.0 : stiffness / area)) :
        ("unit_stiffness", OrderedDict{String, Any}("density" => segment.density,
                                                     "unit_stiffness" => stiffness))
    return KiteGeometry.Segment(; name=definition_name(segment),
        points=segment.point_idxs, rest_length=segment.l0, segment.diameter,
        material=material!(materials, law, parameters),
        model=SegmentSpring(segment.unit_damping, segment.compression_frac,
                            segment.compression_damping_frac))
end

function station_row(sys, station)
    return KiteGeometry.SpanwiseStation(; name=definition_name(station),
        wing=station_wing_idx(sys, station), points=station.point_idxs,
        model=StationModel(string(station.type), station.moment_frac, station.damping,
            station.stiffness, optional_indices(station.body_idxs),
            optional_indices(station.flap_body_idxs),
            optional_indices(station.flap_point_idxs), Vector(station.flap_axis)))
end

"""Indices as references, `nothing` for none."""
optional_indices(idxs) = isempty(idxs) ? nothing : Vector{KiteGeometry.NameRef}(idxs)

"""The wing whose twist `station` carries: its own, else the one listing it, `nothing`
for a station of a body that is no wing."""
function station_wing_idx(sys, station)
    station.wing_idx in eachindex(sys.wings) && return station.wing_idx
    return findfirst(wing -> station.idx in wing.station_idxs, sys.wings)
end

pulley_row(pulley) = KiteGeometry.Pulley(; name=definition_name(pulley),
    segments=pulley.segment_idxs,
    model=PulleyModel(string(pulley.type), pulley.efficiency, pulley.damping,
                      pulley.brake, pulley.friction_epsilon))

tether_row(tether) = KiteGeometry.Tether(; name=definition_name(tether),
    start_point=tether.start_point_idx, end_point=tether.end_point_idx,
    segments=tether.segment_idxs,
    model=TetherModel(tether.init_stretched_len, tether.init_tether_force,
                      tether.init_stretch_frac))

winch_row(winch) = KiteGeometry.Winch(; name=definition_name(winch),
    tethers=winch.tether_idxs, winch_point=winch.winch_point_idx, winch.gear_ratio,
    winch.drum_radius,
    model=WinchModel(winch.init_vel, winch.brake, winch.speed_controlled,
                     winch_friction_epsilon(winch.model)))

"""The friction smoothing width of a winch model, `nothing` for one without."""
winch_friction_epsilon(model) =
    hasproperty(model, :friction_epsilon) ? model.friction_epsilon : nothing

function wing_row(sys, wing)
    return KiteGeometry.Wing(; name=definition_name(wing), pos_ENU=wing.mesh_frame_pos,
        Q_KA2ENU=rotation_matrix_to_quaternion(wing.R_mesh_to_w),
        model=WingModel(string(nameof(typeof(wing.aero))), Vector(wing.angular_damping),
            wing_engine_field(wing, :aero_scale_chord),
            wing_engine_field(wing, :aero_z_offset), string(wing.principal_frame_method),
            string(wing.dynamics_type), ref_points_cell(sys, wing.origin),
            ref_pair_cell(sys, wing.z_ref_points), ref_pair_cell(sys, wing.y_ref_points),
            nothing, nothing, nothing, wing_drag_corr(wing.aero),
            optional_indices(wing.station_idxs)))
end

"""A VSM engine field of `wing`, `nothing` for a wing without one."""
wing_engine_field(wing, field) =
    has_vsm_engine(wing.aero) ? Float64(getproperty(wing, field)) : nothing

"""The drag correction of a flat-plate aero model, `nothing` for another."""
wing_drag_corr(aero) = hasproperty(aero, :drag_corr) ? aero.drag_corr : nothing

"""Weighted reference points as the cell an authoring YAML gives them in."""
ref_points_cell(sys, ref_points::WeightedRefPoints) =
    [Any[definition_name(sys.points[id]), weight]
     for (id, weight) in zip(ref_points.ids, ref_points.weights)]
ref_points_cell(sys, ::Nothing) = nothing

ref_pair_cell(sys, pair) = isnothing(pair) ? nothing :
    [ref_points_cell(sys, pair[1]), ref_points_cell(sys, pair[2])]

function body_row(sys, body)
    R = initial_rotation(body)
    return KiteGeometry.Body(; name=definition_name(body),
        pos_ENU=body.pos_ENU + R * body.extra_com_offset_b, Q_KA2ENU=body.Q_KA_to_ENU,
        body.extra_mass, extra_inertia_KA=body.extra_inertia_b,
        model=BodyModel(string(body.type), optional_index(body.wing_idx),
            Vector(body.angular_damping), Vector(body.world_frame_damping),
            Vector(body.body_frame_damping), body.fix_sphere, body.fix_static,
            Vector(body.ext_force_w), Vector(body.ext_force_b), Vector(body.ext_moment_b),
            string(body.principal_frame_method)))
end

"""The parameters of the material of a tube's rigidity `law`: KiteGeometry's constants
for `breukels2011`, none for another."""
tube_law_parameters(law) = law == :breukels2011 ? copy(KiteGeometry.BREUKELS2011) :
    OrderedDict{String, Any}()

function tube_row!(materials, sys, tube)
    body_a, body_b = sys.bodies[tube.body_a_idx], sys.bodies[tube.body_b_idx]
    anchors = (Vector(tube.anchor_a_b - body_a.extra_com_offset_b),
               Vector(tube.anchor_b_b - body_b.extra_com_offset_b))
    model = tube.model
    linear(field) = (value = getfield(model, field); value isa Real ? value : nothing)
    row_model = model isa ElasticTube ?
        ElasticTubeModel(linear(:stiffness_axial), linear(:stiffness_shear),
                         linear(:stiffness_torsion), linear(:stiffness_bending),
                         model.damping, anchors...) :
        TimoshenkoTubeModel(linear(:EA), linear(:GA), linear(:GJ), linear(:EIy),
                            linear(:EIz), model.shear_coeff, model.damping, anchors...)
    return KiteGeometry.Tube(; name=definition_name(tube),
        bodies=(tube.body_a_idx, tube.body_b_idx), tube.diameter, tube.pressure,
        material=material!(materials, string(tube.law), tube_law_parameters(tube.law)),
        model=row_model)
end

"""
    store_definition_pose!(state::SystemState)

Write the pose `state` holds into its structure's definition, where it has one: each
point's position, each body's at the centre of its own mass, each wing's KA frame carried
with its body, and each segment's rest length.
"""
function store_definition_pose!(state::SystemState)
    sys = getfield(state, :structure)
    definition = sys.definition
    isnothing(definition) && return nothing
    for (row, point) in zip(definition.points, state.points)
        row.pos_ENU .= point.pos_w
    end
    for (row, wing) in zip(definition.wings, state.wings)
        R = quaternion_to_rotation_matrix(wing.Q_b_to_w) * initial_rotation(wing.row)'
        row.pos_ENU .= wing.pos_w + R * (row.pos_ENU - wing.row.pos_ENU)
        row.Q_KA2ENU .= rotation_matrix_to_quaternion(
            R * quaternion_to_rotation_matrix(row.Q_KA2ENU))
    end
    for (row, body) in zip(definition.bodies, state.bodies)
        R = quaternion_to_rotation_matrix(body.Q_b_to_w)
        row.pos_ENU .= body.pos_w + R * body.row.extra_com_offset_b
        row.Q_KA2ENU .= body.Q_b_to_w
    end
    for (row, segment) in zip(definition.segments, state.segments)
        row.rest_length = segment.l0
    end
    return nothing
end
