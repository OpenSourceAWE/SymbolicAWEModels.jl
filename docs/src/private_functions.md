```@meta
CurrentModule = SymbolicAWEModels
```

# Private API

This page documents the internal functions and types of `SymbolicAWEModels.jl`. These are not
part of the public API and may change without notice. They are listed here for
developers and for those interested in the model's internal workings.

## Core Types and Constructors

```@docs
SymbolicAWEModels.SerializedModel
SymbolicAWEModels.SimFloat
SymbolicAWEModels.KVec3
SymbolicAWEModels.SVec3
VortexStepMethod.RamAirWing
SymbolicAWEModels.create_4_attach_ram_sys_struct
SymbolicAWEModels.create_tether
```

## State Management and Model Simplification

```@docs
SymbolicAWEModels.getstate
SymbolicAWEModels.setstate!
SymbolicAWEModels.set_measured!
SymbolicAWEModels.copy!
SymbolicAWEModels.reinit!
SymbolicAWEModels.update_sys_struct!
SymbolicAWEModels.get_set_hash
SymbolicAWEModels.get_sys_struct_hash
```

## Physics and Geometry Helpers

```@docs
SymbolicAWEModels.calc_speed_acc
SymbolicAWEModels.calc_moment_acc
SymbolicAWEModels.calc_angle_of_attack
SymbolicAWEModels.calc_heading
SymbolicAWEModels.calc_R_t_w
SymbolicAWEModels.calc_R_v_w
SymbolicAWEModels.cad_to_body_frame
SymbolicAWEModels.calc_pos
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

## Equations and System Management

```@docs
SymbolicAWEModels.create_sys!
SymbolicAWEModels.scalar_eqs!
SymbolicAWEModels.wing_eqs!
SymbolicAWEModels.linear_vsm_eqs!
SymbolicAWEModels.force_eqs!
SymbolicAWEModels.linearize_vsm!
SymbolicAWEModels.jacobian
SymbolicAWEModels.load_serialized_model!
SymbolicAWEModels.maybe_create_lin_prob!
SymbolicAWEModels.maybe_create_control_functions!
SymbolicAWEModels.maybe_create_prob!
SymbolicAWEModels.maybe_create_simple_lin_model!
SymbolicAWEModels.generate_control_funcs
SymbolicAWEModels.generate_simple_lin_model
SymbolicAWEModels.generate_lin_getters
SymbolicAWEModels.generate_prob_getters
SymbolicAWEModels.LinProbWithAttributes
SymbolicAWEModels.ProbWithAttributes
SymbolicAWEModels.SimpleLinModelWithAttributes
SymbolicAWEModels.ControlFuncWithAttributes
```

## Utility and Internal Functions

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
SymbolicAWEModels.create_default_models
SymbolicAWEModels.copy_bin
SymbolicAWEModels.copy_examples
SymbolicAWEModels.copy_data
SymbolicAWEModels.copy_dir
```

## Base Overloads (Internal Use)

```@docs
Base.getindex
Base.getproperty
Base.setproperty!
```

## Plotting Recipes (Internal Use)

```@docs
RecipesBase.apply_recipe
```
