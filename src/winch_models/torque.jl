# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# TorqueWinch: torque-controlled drum with Coulomb + viscous friction.

"""
    TorqueWinch(; friction_epsilon=6.0)

Torque-controlled winch motor (the default winch model). `set_value` is the
motor torque [N·m]. Coulomb friction is smoothed by `friction_epsilon` (the
`smooth_sign` transition width). The drum parameters (`gear_ratio`,
`drum_radius`, `coulomb_friction`, `viscous_coefficient`, `inertia_total`) live
on the [`Winch`](@ref) struct; `friction_epsilon` is a numerical property of this
model and lives here (mutable, live-tunable).

# Equations
```
ω_motor   = vel / ratio,   ratio = drum_radius / gear_ratio
friction  = smooth_sign(ω_motor, friction_epsilon) * coulomb_friction * ratio +
            viscous_coefficient * ω_motor * ratio^2
tau_total = set_value + ratio * force - friction
acc       = ifelse(brake > 0.5, 0, ratio * tau_total / inertia_total)
```
"""
mutable struct TorqueWinch <: AbstractWinchModel
    "Smoothing width for the Coulomb-friction sign function."
    friction_epsilon::SimFloat
end
TorqueWinch(; friction_epsilon=6.0) =
    TorqueWinch(SimFloat(friction_epsilon))

is_builtin_winch(::TorqueWinch) = true

function winch_component(::TorqueWinch, sys_struct, winch_idx; name, params)
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

    gear_ratio          = params.winches[winch_idx].gear_ratio
    drum_radius         = params.winches[winch_idx].drum_radius
    coulomb_friction    = params.winches[winch_idx].coulomb_friction
    viscous_coefficient = params.winches[winch_idx].viscous_coefficient
    inertia_total       = params.winches[winch_idx].inertia_total
    friction_eps        = params.winches[winch_idx].model.friction_epsilon
    ratio = drum_radius / gear_ratio

    eqs = [
        ω_motor   ~ vel / ratio
        friction  ~ coulomb_viscous_friction(ω_motor, coulomb_friction,
                                             viscous_coefficient, friction_eps,
                                             ratio)
        tau_total ~ set_value + ratio * force - friction
        α_motor   ~ tau_total / inertia_total
        acc       ~ ifelse(brake > 0.5, 0.0, ratio * α_motor)
    ]
    return System(eqs, t,
                  [vel, len, force, set_value, brake, acc, friction,
                   ω_motor, tau_total, α_motor],
                  [gear_ratio, drum_radius, coulomb_friction, viscous_coefficient,
                   inertia_total, friction_eps]; name)
end
