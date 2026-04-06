```@meta
CurrentModule = SymbolicAWEModels
```

# Private API


This page documents the internal functions and types of `SymbolicAWEModels.jl`. These are not
part of the public API and may change without notice. They are listed here for
developers and for those interested in the model's internal workings.

## Core types and constructors

```@docs
SymbolicAWEModels.SerializedModel
SymbolicAWEModels.SimFloat
SymbolicAWEModels.KVec3
SymbolicAWEModels.SVec3
VortexStepMethod.Wing
SymbolicAWEModels.create_tether
SymbolicAWEModels.create_vsm_wing
```

## State management and model simplification

```@docs
SymbolicAWEModels.copy!
SymbolicAWEModels.reinit!
SymbolicAWEModels.reposition!
SymbolicAWEModels.update_sys_struct!
SymbolicAWEModels.get_set_hash
SymbolicAWEModels.get_sys_struct_hash
```

## Physics and geometry helpers

```@docs
SymbolicAWEModels.calc_angle_of_attack
SymbolicAWEModels.calc_heading
SymbolicAWEModels.calc_R_t_to_w
SymbolicAWEModels.calc_R_v_to_w
SymbolicAWEModels.cad_to_body_frame
SymbolicAWEModels.calc_pos
SymbolicAWEModels.calc_winch_force
SymbolicAWEModels.find_axis_point
SymbolicAWEModels.quaternion_to_rotation_matrix
SymbolicAWEModels.rotation_matrix_to_quaternion
SymbolicAWEModels.rotate_v_around_k
SymbolicAWEModels.sym_normalize
SymbolicAWEModels.apply_heading
SymbolicAWEModels.get_rot_pos
SymbolicAWEModels.get_base_pos
SymbolicAWEModels.calc_aoa
```

## Equations and system management

```@docs
SymbolicAWEModels.create_sys!
SymbolicAWEModels.scalar_eqs!
SymbolicAWEModels.wing_eqs!
SymbolicAWEModels.vsm_eqs!
SymbolicAWEModels.point_eqs!
SymbolicAWEModels.segment_eqs!
SymbolicAWEModels.update_vsm!
SymbolicAWEModels.jacobian
SymbolicAWEModels.load_serialized_model!
SymbolicAWEModels.maybe_create_lin_prob!
SymbolicAWEModels.maybe_create_control_functions!
SymbolicAWEModels.maybe_create_prob!
SymbolicAWEModels.generate_control_funcs
SymbolicAWEModels.generate_lin_getters
SymbolicAWEModels.generate_prob_getters
SymbolicAWEModels.LinProbWithAttributes
SymbolicAWEModels.ProbWithAttributes
SymbolicAWEModels.ControlFuncWithAttributes
```

## Utility and internal functions

```@docs
SymbolicAWEModels.get_model_name
SymbolicAWEModels.calc_height
SymbolicAWEModels.set_depower_steering!
SymbolicAWEModels.min_chord_len
SymbolicAWEModels.pos
SymbolicAWEModels.spring_forces
SymbolicAWEModels.calc_spring_props
SymbolicAWEModels.set_v_wind_ground!
SymbolicAWEModels.in_percent_band
SymbolicAWEModels.step
SymbolicAWEModels.create_model_archive
SymbolicAWEModels.filecmp
SymbolicAWEModels.extract_model_archive
SymbolicAWEModels.copy_bin
SymbolicAWEModels.copy_examples
SymbolicAWEModels.copy_data
SymbolicAWEModels.copy_dir
SymbolicAWEModels.get_example_packages
SymbolicAWEModels.make_lin_sys_state
```

## Base overloads (internal use)

```@docs
SymbolicAWEModels._SAM_FIELDS
Base.getindex
Base.getproperty
Base.setproperty!
```

## YAML loader internals

```@docs
SymbolicAWEModels.get_field_or_nothing
SymbolicAWEModels.convert_to_type
SymbolicAWEModels.resolve_references
SymbolicAWEModels.calculate_derived_properties!
SymbolicAWEModels._extract_args
SymbolicAWEModels.call_yaml_constructor
SymbolicAWEModels.update_yaml_from_sys_struct!
SymbolicAWEModels.update_aero_yaml_from_struc_yaml!
```

## SystemStructure internals

```@docs
SymbolicAWEModels.segment_cad_length
SymbolicAWEModels.segment_world_length
SymbolicAWEModels.autocalc_tether_len
SymbolicAWEModels.apply_tether_init_lens!
SymbolicAWEModels.assign_indices_and_resolve!
SymbolicAWEModels.resolve_ref
SymbolicAWEModels.resolve_ref_spec
SymbolicAWEModels.validate_sys_struct
SymbolicAWEModels.build_name_dict
SymbolicAWEModels.identify_wing_segments
SymbolicAWEModels.match_aero_sections_to_structure!
SymbolicAWEModels.compute_spatial_group_mapping!
SymbolicAWEModels.copy_cad_to_world!
SymbolicAWEModels.adjust_vsm_panels_to_origin!
SymbolicAWEModels.apply_aero_z_offset!
SymbolicAWEModels.calc_refine_wing_frame
SymbolicAWEModels.calc_inertia_y_rotation
SymbolicAWEModels.rotate_vsm_sections!
SymbolicAWEModels.expand_auto_tethers!
SymbolicAWEModels.WeightedRefPoints
SymbolicAWEModels.resolve!
SymbolicAWEModels._validate_weights!
SymbolicAWEModels.SegmentType
```

## NamedCollection internals

```@docs
SymbolicAWEModels.names
SymbolicAWEModels.get_idx
SymbolicAWEModels.get_name
Base.keys(::NamedCollection)
Base.values(::NamedCollection)
Base.haskey(::NamedCollection, ::Symbol)
Base.setindex!(::NamedCollection, ::Any, ::Integer)
Base.setindex!(::NamedCollection, ::Any, ::Symbol)
```

## Equation builders

```@docs
SymbolicAWEModels.tether_eqs!
SymbolicAWEModels.pulley_eqs!
SymbolicAWEModels.winch_eqs!
SymbolicAWEModels.group_eqs!
```

## VSM and aerodynamics internals

```@docs
SymbolicAWEModels.build_point_to_vsm_point_mapping
SymbolicAWEModels.update_vsm_wing_from_structure!
SymbolicAWEModels.distribute_panel_forces_to_points!
SymbolicAWEModels.get_aero_force_override
SymbolicAWEModels.get_aero_moment_override
SymbolicAWEModels.get_group_moment_override
```

## Heading and geometry

```@docs
SymbolicAWEModels.solve_heading_rotation
SymbolicAWEModels.get_ref_position_from_points
SymbolicAWEModels.sym_calc_R_t_to_w
SymbolicAWEModels.wrap_to_pi
```

## Transform internals

```@docs
SymbolicAWEModels._apply_azimuth_elevation!
SymbolicAWEModels._apply_heading!
SymbolicAWEModels._finalize_transforms!
```

## Other internals

```@docs
SymbolicAWEModels.init_principal_frame!
SymbolicAWEModels.update_segment_forces!
KiteUtils.Logger(::SymbolicAWEModel, ::Int64)
```