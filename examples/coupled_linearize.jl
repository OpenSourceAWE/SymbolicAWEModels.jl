# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

"""
Linearize the 2-plate kite model at a steady-state operating point
and build a state-space representation.
"""

using SymbolicAWEModels, VortexStepMethod
using ModelingToolkit
using ModelingToolkit: t_nounits
using ControlSystemsBase

set_data_path("data/2plate_kite")
struc_yaml = joinpath(get_data_path(),
                      "quat_struc_geometry.yaml")
aero_yaml = joinpath(get_data_path(), "aero_geometry.yaml")
update_aero_yaml_from_struc_yaml!(struc_yaml, aero_yaml)

set = Settings("system.yaml")
vsm_set = VortexStepMethod.VSMSettings(
    joinpath(get_data_path(), "vsm_settings.yaml");
    data_prefix=false)

sys = load_sys_struct_from_yaml(struc_yaml;
    system_name="2plate_kite", set, vsm_set)
sam = SymbolicAWEModel(set, sys)

@variables begin
    heading(t_nounits)[1:1]
    angle_of_attack(t_nounits)[1:1]
    tether_len(t_nounits)[1:1]
    winch_force(t_nounits)[1:1]
end
outputs = [heading[1], angle_of_attack[1], tether_len[1],
           winch_force[1]]

init!(sam; outputs, create_lin_prob=true)
find_steady_state!(sam)

(; A, B, C, D) = SymbolicAWEModels.linearize!(sam)
@info "Linearized" A=size(A) B=size(B) C=size(C) D=size(D)

lin_sys = ss(A, B, C, D)
@info "State-space model" lin_sys
