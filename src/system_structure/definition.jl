# Copyright (c) 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

"""
    system_definition(sys::SystemStructure) -> SystemDefinition

The KiteGeometry definition of `sys` in its initial pose: its components as rows, every
reference an index, a station only where it is a wing's, what only SymbolicAWEModels
reads in the models it registers ([`register_definition_models`](@ref)), and each body at
the centre of its own mass.
"""
function system_definition(sys)
    materials = sys.materials
    for station in sys.stations
        isnothing(station_wing_idx(sys, station)) || continue
        push!(body_model(sys.bodies[station.wing_idx]).stations_without_wing, station)
    end
    return SystemDefinition(;
        metadata=Metadata(; name=sys.name, n_points=length(sys.points)),
        points=[point_row!(point) for point in sys.points],
        segments=collect(sys.segments),
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

"""The name a definition gives component `component`: its own, else its index."""
definition_name(component) = Symbol(something(component.name, component.idx))

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

"""`point` naming the body and the tube it rides by the indices they resolve to."""
function point_row!(point)
    point.body = reference_column(point.body_idx)
    point.tube = reference_column(point.tube_idx)
    return point
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
    origin, R = mesh_frame(wing)
    return KiteGeometry.Wing(; name=definition_name(wing), pos_ENU=origin,
        Q_KA2ENU=rotation_matrix_to_quaternion(R), body=wing.idx,
        model=WingModel(aero_mode_name(wing.aero), Vector(wing.angular_damping),
            wing_engine_field(wing, :aero_scale_chord),
            wing_engine_field(wing, :aero_z_offset), string(wing.principal_frame_method),
            string(wing.dynamics_type), ref_points_cell(sys, wing.origin),
            ref_pair_cell(sys, wing.z_ref_points), ref_pair_cell(sys, wing.y_ref_points),
            nothing, nothing, nothing, wing_drag_corr(wing.aero),
            references_column(wing.station_idxs)))
end

"""The name of an aero mode a wing row gives."""
aero_mode_name(mode) = string(nameof(typeof(mode)))
aero_mode_name(::AeroSpec{M}) where {M} = string(nameof(M))

"""A VSM adjustment of `wing`, `nothing` for a wing flown without VSM."""
wing_engine_field(wing, field) =
    wing.aero isa AeroSpec ? Float64(getfield(body_model(wing), field)) : nothing

"""The drag correction of a flat-plate aero model, `nothing` for another."""
wing_drag_corr(aero) = hasproperty(aero, :drag_corr) ? aero.drag_corr : nothing

"""Weighted reference points as the cell an authoring YAML gives them in."""
ref_points_cell(sys, ref_points::WeightedRefPoints) =
    [Any[definition_name(sys.points[id]), weight]
     for (id, weight) in zip(ref_points.ids, ref_points.weights)]
ref_points_cell(sys, ::Nothing) = nothing

ref_pair_cell(sys, pair) = isnothing(pair) ? nothing :
    [ref_points_cell(sys, pair[1]), ref_points_cell(sys, pair[2])]

"""`body` with its inertia about its own centre of mass."""
function body_row(sys, body)
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
    row_model.simulated = tube
    return KiteGeometry.Tube(; name=definition_name(tube),
        bodies=(tube.body_a_idx, tube.body_b_idx), tube.diameter, tube.pressure,
        material=material!(materials, string(tube.law), tube_law_parameters(tube.law)),
        model=row_model)
end
