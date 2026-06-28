# Copyright (c) 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# CascadedLengthWinch: length-controlled drum with a reel-speed cap.

"""
    CascadedLengthWinch(; v_max, position_gain, velocity_gain, friction_epsilon=6.0)

Length-controlled winch motor. `set_value` is the target tether length [m]. An
outer proportional law turns the length error into a velocity reference, hard-
clamped at `±v_max`; an inner proportional law with load feedforward tracks it.
The feedforward (`+friction − ratio·force`) cancels the load disturbance so
tether spring-mass bouncing does not perturb the reel speed. No integrator, so
there is no windup at the `v_max` clamp and no overshoot; the closed-loop speed
response is first-order.

Drum parameters (`gear_ratio`, `drum_radius`, `f_coulomb`, `c_vf`,
`inertia_total`) live on the [`Winch`](@ref); the gains and `friction_epsilon`
are this model's own fields (mutable, live-tunable). A low `v_max` gives a slow
quasi-static reel-in, useful for displacement-controlled load sweeps.

# Equations
```
ω_motor       = vel / ratio,   ratio = drum_radius / gear_ratio
friction      = smooth_sign(ω_motor, friction_epsilon) * f_coulomb * ratio +
                c_vf * ω_motor * ratio^2
vel_unclamped = position_gain * (set_value − len)
vel_ref       = clamp(vel_unclamped, −v_max, v_max)
tau_cmd       = velocity_gain * (vel_ref − vel) + friction − ratio * force
tau_net       = tau_cmd + ratio * force − friction
acc           = ifelse(brake > 0.5, 0, ratio * tau_net / inertia_total)
```
"""
mutable struct CascadedLengthWinch <: AbstractWinchModel
    "Reel-speed cap [m/s]; the velocity reference is clamped to `±v_max`."
    v_max::SimFloat
    "Outer-loop gain [1/s]: length error → velocity reference."
    position_gain::SimFloat
    "Inner-loop gain [N·m·s/m]: velocity error → motor torque."
    velocity_gain::SimFloat
    "Smoothing width for the Coulomb-friction sign function."
    friction_epsilon::SimFloat
end
CascadedLengthWinch(; v_max, position_gain, velocity_gain, friction_epsilon=6.0) =
    CascadedLengthWinch(SimFloat(v_max), SimFloat(position_gain),
                        SimFloat(velocity_gain), SimFloat(friction_epsilon))

is_builtin_winch(::CascadedLengthWinch) = true

function winch_component(::CascadedLengthWinch, sys_struct, winch_idx; name, params)
    @variables begin
        vel(t)
        len(t)
        force(t)
        set_value(t)
        brake(t)
        acc(t)
        friction(t)
        ω_motor(t)
        vel_unclamped(t)
        vel_ref(t)
        tau_cmd(t)
        tau_net(t)
    end

    gear_ratio    = params.winches[winch_idx].gear_ratio
    drum_radius   = params.winches[winch_idx].drum_radius
    f_coulomb     = params.winches[winch_idx].f_coulomb
    c_vf          = params.winches[winch_idx].c_vf
    inertia_total = params.winches[winch_idx].inertia_total
    friction_eps  = params.winches[winch_idx].model.friction_epsilon
    v_max         = params.winches[winch_idx].model.v_max
    position_gain = params.winches[winch_idx].model.position_gain
    velocity_gain = params.winches[winch_idx].model.velocity_gain
    smooth_sign(x, eps) = x / sqrt(x * x + eps * eps)
    ratio = drum_radius / gear_ratio

    eqs = [
        ω_motor       ~ vel / ratio
        friction      ~ smooth_sign(ω_motor, friction_eps) * f_coulomb * ratio +
                        c_vf * ω_motor * ratio^2
        vel_unclamped ~ position_gain * (set_value - len)
        vel_ref       ~ max(-v_max, min(v_max, vel_unclamped))
        tau_cmd       ~ velocity_gain * (vel_ref - vel) + friction - ratio * force
        tau_net       ~ tau_cmd + ratio * force - friction
        acc           ~ ifelse(brake > 0.5, 0.0,
                               ratio * tau_net / inertia_total)
    ]
    return System(eqs, t,
                  [vel, len, force, set_value, brake, acc, friction,
                   ω_motor, vel_unclamped, vel_ref, tau_cmd, tau_net],
                  [gear_ratio, drum_radius, f_coulomb, c_vf, inertia_total,
                   friction_eps, v_max, position_gain, velocity_gain]; name)
end
