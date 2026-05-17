# Copyright (c) 2025 Bart van de Lint, Uwe Fechner
# SPDX-License-Identifier: LGPL-3.0-only

"""
Linearize the 2-plate kite model at a steady-state operating point
and build a state-space representation.
"""

using Pkg
if Base.active_project() != joinpath(@__DIR__, "Project.toml")
    Pkg.activate(@__DIR__)
end

using KiteUtils: init!
using SymbolicAWEModels, VortexStepMethod
using ControlSystemsBase
using ModelingToolkit: @variables, t_nounits

set_data_path("data/2plate_kite")
struc_yaml = joinpath(get_data_path(),
                      "quat_struc_geometry.yaml")

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
@info "State-space model: $(size(A, 1)) states, " *
    "$(size(B, 2)) inputs, $(size(C, 1)) outputs"
