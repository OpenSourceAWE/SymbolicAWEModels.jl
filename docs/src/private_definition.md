```@meta
CurrentModule = SymbolicAWEModels
```

# Private API — definition

The internals that read a KiteGeometry definition into a `SystemDefinition` and
write one back: the loader, the models SymbolicAWEModels registers for each block, the
components as definition rows, and the writer. Not part of the public API.

## Definition loader internals

```@docs
SymbolicAWEModels.optional_project_file
SymbolicAWEModels.set_wing_column!
SymbolicAWEModels.wings_first
SymbolicAWEModels.definition_tube
```

## KiteGeometry models

```@docs
SymbolicAWEModels.WingModel
SymbolicAWEModels.wing_aero_spec
SymbolicAWEModels.wing_dynamics_type
SymbolicAWEModels.wing_station_idxs
SymbolicAWEModels.unit_inertia_tensor
KiteGeometry.wing_spec(::SymbolicAWEModels.WingModel)
SymbolicAWEModels.BodyModel
SymbolicAWEModels.register_definition_models
```

## Component definition rows

```@docs
SymbolicAWEModels.per_axis
SymbolicAWEModels.material_unit_stiffness
SymbolicAWEModels.to_float64
SymbolicAWEModels.wing_column
SymbolicAWEModels.body_wing
SymbolicAWEModels.own_com_offset
SymbolicAWEModels.mesh_frame
SymbolicAWEModels.column_vector
SymbolicAWEModels.principal_axes
SymbolicAWEModels.axis_column
SymbolicAWEModels.body_principal_frame_method
SymbolicAWEModels.body_angular_damping
SymbolicAWEModels.origin_ENU
SymbolicAWEModels.set_origin_ENU!
```

## Definition writer internals

```@docs
SymbolicAWEModels.torque_winch
```

## Segments, tubes and a run's reads

```@docs
SymbolicAWEModels.BuiltSegment
SymbolicAWEModels.PointPositions
SymbolicAWEModels.body_apparent_wind
SymbolicAWEModels.built
SymbolicAWEModels.copy_run_columns!
SymbolicAWEModels.model_aero
SymbolicAWEModels.model_vsm_set
SymbolicAWEModels.copy_unsteady
SymbolicAWEModels.numeric_unit_stiffness
SymbolicAWEModels.point_apparent_wind
SymbolicAWEModels.point_segments
SymbolicAWEModels.point_wings
SymbolicAWEModels.segment_tethers
SymbolicAWEModels.segment_pulleys
SymbolicAWEModels.run_value
SymbolicAWEModels.set_run_value!
SymbolicAWEModels.segment_rest_length
SymbolicAWEModels.BodyGeometry
SymbolicAWEModels.body_geometry
SymbolicAWEModels.aero_mode
SymbolicAWEModels.BODY_GEOMETRY_FIELDS
SymbolicAWEModels.STATION_GEOMETRY_FIELDS
SymbolicAWEModels.StationGeometry
SymbolicAWEModels.station_geometry
SymbolicAWEModels.station_aero_moment
SymbolicAWEModels.run_station_aero_moment
SymbolicAWEModels.station_aero_moment_reader
SymbolicAWEModels.hold_station_moments!
SymbolicAWEModels.point_is_wing_node
SymbolicAWEModels.point_anchor_b
SymbolicAWEModels.pulley_sum_len
SymbolicAWEModels.winch_setting
SymbolicAWEModels.point_undeformed_b
SymbolicAWEModels.point_beam_anchor
SymbolicAWEModels.beam_anchor
SymbolicAWEModels.collection
SymbolicAWEModels.start_lengths!
SymbolicAWEModels.point_total_mass
SymbolicAWEModels.resolve_segment
SymbolicAWEModels.resolved_unit_damping
SymbolicAWEModels.run_tubes
SymbolicAWEModels.segment_spring
SymbolicAWEModels.station_tether_load
SymbolicAWEModels.typed_rows
SymbolicAWEModels.winch_force_vector
SymbolicAWEModels.wing_bodies
SymbolicAWEModels.wing_rates
```
