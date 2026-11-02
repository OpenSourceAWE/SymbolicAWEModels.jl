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
        points=collect(sys.points),
        segments=[segment_row!(materials, segment) for segment in sys.segments],
        spanwise_stations=[station_row(sys, station) for station in sys.stations
                           if !isnothing(station_wing_idx(sys, station))],
        pulleys=collect(sys.pulleys),
        tethers=[with_ends!(tether) for tether in sys.tethers],
        winches=[winch_row!(winch) for winch in sys.winches],
        wings=[wing_row(sys, wing) for wing in sys.wings],
        bodies=[body_row(sys, body) for body in sys.bodies],
        tubes=[tube_row!(materials, sys, tube) for tube in sys.tubes],
        materials)
end

"""
    attach_definition!(sys::SystemStructure)

Give `sys` the definition [`system_definition`](@ref) writes from its components, each
point, segment, station, pulley, tether, winch, body and wing of `sys` then being its
row there.
"""
function attach_definition!(sys::SystemStructure)
    definition = sys.definition = system_definition(sys)
    for (components, rows) in ((sys.points, definition.points),
                               (sys.segments, definition.segments),
                               (sys.stations, definition.spanwise_stations),
                               (sys.pulleys, definition.pulleys),
                               (sys.tethers, definition.tethers),
                               (sys.winches, definition.winches),
                               (sys.bodies, definition.bodies),
                               (sys.wings, definition.bodies))
        replace_with_rows!(components.items, rows)
    end
    for (wing, row) in zip(sys.wings, definition.wings)
        body_model(wing).wing_definition = row
    end
    return definition
end

"""Replace each of `components` by the row of `rows` that holds its model."""
function replace_with_rows!(components, rows)
    position = IdDict(component.model => i for (i, component) in enumerate(components))
    for row in rows
        i = get(position, row.model, nothing)
        isnothing(i) || (components[i] = row)
    end
    return components
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

"""`segment` with its stiffness and density as a material it appends to `materials`."""
function segment_row!(materials, segment)
    area = π * (segment.diameter / 2)^2
    stiffness = segment.unit_stiffness
    law, parameters = stiffness isa Real && (area > 0 || iszero(stiffness)) ?
        ("linear", OrderedDict{String, Any}("density" => segment.density,
            "youngs_modulus" => iszero(stiffness) ? 0.0 : stiffness / area)) :
        ("unit_stiffness", OrderedDict{String, Any}("density" => segment.density,
                                                     "unit_stiffness" => stiffness))
    segment.material = material!(materials, law, parameters)
    return segment
end

"""Indices as references, `nothing` for none."""
optional_indices(idxs) = isempty(idxs) ? nothing : Vector{KiteGeometry.NameRef}(idxs)

"""`station` with the wing whose twist it carries."""
function station_row(sys, station)
    station.wing_ref = station_wing_idx(sys, station)
    return station
end

"""The wing whose twist `station` carries: its own, else the one listing it, `nothing`
for a station of a body that is no wing."""
function station_wing_idx(sys, station)
    station.wing_idx in eachindex(sys.wings) && return station.wing_idx
    return findfirst(wing -> station.idx in wing.station_idxs, sys.wings)
end

"""`tether` naming the ends it resolves to, those its segments give where it names
none."""
function with_ends!(tether)
    tether.start_point = tether.start_point_idx
    tether.end_point = tether.end_point_idx
    return tether
end

"""`winch` naming its tethers and ground point by the indices they resolve to."""
function winch_row!(winch)
    winch.tethers = winch.tether_idxs
    winch.winch_point = winch.winch_point_idx
    return winch
end

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

"""`body` with its row at the centre of its own mass in its initial pose, its inertia in
that frame."""
function body_row(sys, body)
    R = initial_rotation(body)
    getfield(body, :pos_ENU) .= body.pos_ENU + R * body.extra_com_offset_b
    getfield(body, :Q_KA2ENU) .= body.Q_KA_to_ENU
    setfield!(body, :extra_inertia_KA, MMatrix{3, 3, Float64, 9}(body.extra_inertia_b))
    return body
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
body's at the centre of its own mass and each wing's KA frame carried with its body.
"""
function store_definition_pose!(state::SystemState)
    sys = getfield(state, :structure)
    definition = sys.definition
    isnothing(definition) && return nothing
    for (row, wing) in zip(definition.wings, state.wings)
        R = quaternion_to_rotation_matrix(wing.Q_b_to_w) * initial_rotation(wing.row)'
        row.pos_ENU .= wing.pos_w + R * (row.pos_ENU - wing.row.pos_ENU)
        row.Q_KA2ENU .= rotation_matrix_to_quaternion(
            R * quaternion_to_rotation_matrix(row.Q_KA2ENU))
    end
    for (row, body) in zip(definition.bodies, state.bodies)
        R = quaternion_to_rotation_matrix(body.Q_b_to_w)
        getfield(row, :pos_ENU) .= body.pos_w + R * body.row.extra_com_offset_b
        getfield(row, :Q_KA2ENU) .= body.Q_b_to_w
    end
    return nothing
end
