# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# DefaultWinchModel: torque-controlled drum with Coulomb + viscous friction.
# `set_value` is interpreted as motor torque [N·m]. See common.jl for the
# interface.

"""
    DefaultWinchModel(; friction_epsilon=6.0)

Torque-controlled winch motor. `set_value` is the motor torque [N·m]. Coulomb
friction is smoothed by `friction_epsilon` (the `smooth_sign` transition width).
The drum parameters (`gear_ratio`, `drum_radius`, `f_coulomb`, `c_vf`,
`inertia_total`) live on the [`Winch`](@ref) struct; `friction_epsilon` is a
numerical property of this model and lives here (mutable, live-tunable).

# Equations
```
ω_motor   = vel / ratio,   ratio = drum_radius / gear_ratio
friction  = smooth_sign(ω_motor, friction_epsilon) * f_coulomb * ratio +
            c_vf * ω_motor * ratio^2
tau_total = set_value + ratio * force - friction
acc       = ifelse(brake > 0.5, 0, ratio * tau_total / inertia_total)
```
"""
mutable struct DefaultWinchModel <: AbstractWinchModel
    "Smoothing width for the Coulomb-friction sign function."
    friction_epsilon::SimFloat
end
DefaultWinchModel(; friction_epsilon=6.0) =
    DefaultWinchModel(SimFloat(friction_epsilon))

is_builtin_winch(::DefaultWinchModel) = true

function winch_component(::DefaultWinchModel, sys_struct, winch_idx; name, params)
    @variables begin
        vel(t)
        len(t)
        force(t)
        set_value(t)
        brake(t)
        acc(t)
        friction(t)
        ω_motor(t)
        tau_total(t)
        α_motor(t)
    end

    gear_ratio    = params.winches[winch_idx].gear_ratio
    drum_radius   = params.winches[winch_idx].drum_radius
    f_coulomb     = params.winches[winch_idx].f_coulomb
    c_vf          = params.winches[winch_idx].c_vf
    inertia_total = params.winches[winch_idx].inertia_total
    friction_eps  = params.winches[winch_idx].model.friction_epsilon
    smooth_sign(x, eps) = x / sqrt(x * x + eps * eps)
    ratio = drum_radius / gear_ratio

    eqs = [
        ω_motor   ~ vel / ratio
        friction  ~ smooth_sign(ω_motor, friction_eps) * f_coulomb * ratio +
                    c_vf * ω_motor * ratio^2
        tau_total ~ set_value + ratio * force - friction
        α_motor   ~ tau_total / inertia_total
        acc       ~ ifelse(brake > 0.5, 0.0, ratio * α_motor)
    ]
    return System(eqs, t,
                  [vel, len, force, set_value, brake, acc, friction,
                   ω_motor, tau_total, α_motor],
                  [gear_ratio, drum_radius, f_coulomb, c_vf, inertia_total,
                   friction_eps]; name)
end
