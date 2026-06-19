# SPDX-FileCopyrightText: 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

"""
Custom winch model: cascaded length-velocity control with a v_max cap.

`set_value` is interpreted as a target tether length [m]. An outer
proportional law turns the length error into a velocity reference,
hard-clamped at ±v_max. An inner P-on-velocity controller with load
feedforward tracks that reference. The feedforward (`+friction -
ratio·force`) cancels the load disturbance so tether spring-mass
bouncing doesn't perturb `vel`. No integrator → no windup at
saturation → no overshoot. Closed-loop velocity response is 1st-order
with time constant `I / (ratio · K_p)`.
"""

using Pkg
if Base.active_project() != joinpath(@__DIR__, "Project.toml")
    Pkg.activate(@__DIR__)
end

using GLMakie
using KiteUtils: init!, next_step!
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using SymbolicAWEModels
using SymbolicAWEModels: AbstractWinchModel
import SymbolicAWEModels: Point  # resolve ambiguity with GLMakie

set_data_path(joinpath(dirname(@__DIR__), "data"))
set = Settings("base/system.yaml")
set.v_wind = 0.0

"""
Custom winch model: `set_value` is a target tether length. Drum params come from
the flat `params` view (mirrors `DefaultWinchModel`); the controller gains and
`friction_epsilon` are the model's own fields.
"""
mutable struct CascadedLengthWinch <: AbstractWinchModel
    v_max::Float64
    K_pos::Float64
    K_p::Float64
    friction_epsilon::Float64
end
CascadedLengthWinch(; v_max, K_pos, K_p, friction_epsilon=6.0) =
    CascadedLengthWinch(v_max, K_pos, K_p, friction_epsilon)

function SymbolicAWEModels.winch_component(model::CascadedLengthWinch,
        sys_struct, winch_idx; name, params)
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
    smooth_sign(x, eps) = x / sqrt(x * x + eps * eps)
    ratio = drum_radius / gear_ratio

    eqs = [
        ω_motor       ~ vel / ratio
        friction      ~ smooth_sign(ω_motor, friction_eps) *
                        f_coulomb * ratio +
                        c_vf * ω_motor * ratio^2
        vel_unclamped ~ model.K_pos * (set_value - len)
        vel_ref       ~ max(-model.v_max, min(model.v_max, vel_unclamped))
        tau_cmd       ~ model.K_p * (vel_ref - vel) + friction - ratio * force
        tau_net       ~ tau_cmd + ratio * force - friction
        acc           ~ ifelse(brake > 0.5, 0.0,
                               ratio * tau_net / inertia_total)
    ]
    return System(eqs, t,
                  [vel, len, force, set_value, brake, acc, friction,
                   ω_motor, vel_unclamped, vel_ref, tau_cmd, tau_net],
                  [gear_ratio, drum_radius, f_coulomb, c_vf, inertia_total,
                   friction_eps]; name)
end

points = [
    Point(:ground, [0.0, 0.0, 0.0], STATIC),
    Point(:mass, [0.0, 0.0, -50.0], DYNAMIC; extra_mass=10.0),
]
segments = [
    Segment(:line, :ground, :mass, 50_000.0, 500.0, 0.005; l0=50.0),
]
tethers = [Tether(:main, [:line], 50.0)]
winches = [Winch(:winch, set, [:main]; winch_point=:ground,
                 model=CascadedLengthWinch(v_max=1.0, K_pos=50.0, K_p=10.0))]
winches[1].inertia_total = 0.001  # tiny rotor → near-instant velocity tracking

sys_struct = SystemStructure("custom_winch_vel", set;
                             points, segments, tethers, winches)
sam = SymbolicAWEModel(set, sys_struct)
init!(sam; remake=true, prn=false)

winch  = sam.sys_struct.winches[:winch]
tether = sam.sys_struct.tethers[:main]

dt = 0.02
n_steps = 1500
times    = Float64[]
ref_len  = Float64[]
meas_len = Float64[]
meas_vel = Float64[]

for k in 1:n_steps
    tnow = k * dt
    l_target =
        tnow < 2.0  ? 50.0 :
        tnow < 15.0 ? 45.0 :
        55.0
    next_step!(sam; set_values=[l_target], dt=dt, vsm_interval=0)
    push!(times,    tnow)
    push!(ref_len,  l_target)
    push!(meas_len, tether.len)
    push!(meas_vel, winch.vel)
end

fig = Figure(size=(900, 600))
ax1 = Axis(fig[1, 1]; ylabel="tether length [m]",
           title="Cascaded length→velocity winch (v_max=1.0 m/s)")
lines!(ax1, times, ref_len;  label="target",   linestyle=:dash)
lines!(ax1, times, meas_len; label="measured")
axislegend(ax1; position=:rb)

ax2 = Axis(fig[2, 1]; xlabel="t [s]",
           ylabel="reel-out velocity [m/s]")
hlines!(ax2, [+1.0, -1.0]; color=:gray, linestyle=:dot)
lines!(ax2, times, meas_vel)

display(fig)
