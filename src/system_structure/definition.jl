# Copyright (c) 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

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

"""The name of an aero mode a wing row gives."""
aero_mode_name(mode) = string(nameof(typeof(mode)))
aero_mode_name(::AeroSpec{M}) where {M} = string(nameof(M))

"""The `tubes` row of a tube built in Julia among `bodies`, its anchors measured from each
body's own centre of mass and its material, of its law, added to `materials`."""
function tube_row!(materials, bodies, tube)
    body(ref) = ref isa Integer ? bodies[ref] : bodies[findfirst(b -> b.name == ref, bodies)]
    anchors = (Vector(tube.anchor_a_b - own_com_offset(body(tube.bodies[1]))),
               Vector(tube.anchor_b_b - own_com_offset(body(tube.bodies[2]))))
    model = tube.model
    rigidities = map(field -> getfield(model, field), rigidity_fields(model))
    row_model = model isa ElasticTube ?
        ElasticTubeModel(rigidities..., model.damping, anchors...) :
        TimoshenkoTubeModel(rigidities..., model.shear_coeff, model.damping, anchors...)
    return KiteGeometry.Tube(; tube.name, tube.bodies, tube.diameter, tube.pressure,
        material=material!(materials, string(tube.material),
                           KiteGeometry.tube_law_parameters(tube.material)),
        model=row_model)
end
