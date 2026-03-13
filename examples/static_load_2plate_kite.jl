# Copyright (c) 2025 Bart van de Lint, Jelle Poland
# SPDX-License-Identifier: MPL-2.0

"""
Static load test: apply external aerodynamic forces to the 2-plate
kite structure and simulate a few steps under those loads.
"""

using GLMakie
using SymbolicAWEModels, VortexStepMethod, KiteUtils
using LinearAlgebra
using SymbolicAWEModels: Point

MODEL_NAME = "2plate_kite"
n_steps = 3

# External aerodynamic loads per wing point [Fx, Fy, Fz]
F_AERO = [
     0     0     0;    # ground
   -50    50   225;    # right LE
    20    50   100;    # right TE
  -135     0   400;    # center LE
    50     0   250;    # center TE
   -50   -50   225;    # left LE
    20   -50   100;    # left TE
]

pkg_root = dirname(@__DIR__)
set_data_path(joinpath(pkg_root, "data", MODEL_NAME))

struc_yaml = joinpath(get_data_path(),
                      "quat_struc_geometry.yaml")
aero_yaml = joinpath(get_data_path(), "aero_geometry.yaml")
update_aero_yaml_from_struc_yaml!(struc_yaml, aero_yaml)

set = Settings("system.yaml")
vsm_set = VSMSettings(joinpath(get_data_path(),
    "vsm_settings.yaml"); data_prefix=false)

sys = load_sys_struct_from_yaml(struc_yaml;
    system_name=MODEL_NAME, set, vsm_set,
    aero_mode=AERO_NONE)

function apply_loads!(points, F)
    for (i, point) in enumerate(points)
        if i <= size(F, 1)
            point.disturb .= F[i, :]
        end
    end
end

apply_loads!(sys.points, F_AERO)

sam = SymbolicAWEModel(set, sys)
init!(sam; remake=false)
apply_loads!(sam.sys_struct.points, F_AERO)

logger = Logger(sam, n_steps)
sys_state = SysState(sam)

for step in 1:n_steps
    apply_loads!(sam.sys_struct.points, F_AERO)
    next_step!(sam)
    update_sys_state!(sys_state, sam)
    sys_state.time = step / set.sample_freq
    log!(logger, sys_state)
end

save_log(logger, "static_load")
syslog = load_log("static_load")
scene = replay(syslog, sam.sys_struct)
display(scene)
@info "Static load test complete" steps=n_steps
