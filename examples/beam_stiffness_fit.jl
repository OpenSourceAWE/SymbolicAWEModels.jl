# Copyright (c) 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Fit the bending stiffness of a 10-segment rigid-body beam to a measured
# root-torque-vs-tip-angle curve (e.g. a cantilever bend test of an inflatable
# tube). Clamp the root, apply a pure couple at the free tip, with no gravity:
# static equilibrium makes the internal bending moment constant — equal to the
# applied couple = the reaction torque at the root — at every joint. With N =
# n_segments - 1 identical joints each bending by Δθ, the tip angle is
# γ_tip = N·Δθ. So a measured law M = p(γ_tip) becomes the per-joint moment law
# f_joint(Δθ) = p(N·Δθ): the SAME polynomial, with the angle pre-scaled by N.
# Only the bending DOF is fitted; axial/shear/torsion are left stiff and linear.

using Pkg
Pkg.activate(@__DIR__)
using SymbolicAWEModels
using KiteUtils
using GLMakie

# ----- beam parameters -----
n_segments = 10
seg_length = 0.5                 # [m] per segment
seg_mass = 0.5                   # [kg] per segment
seg_radius = 0.02                # [m] equivalent rod radius (axial inertia)
n_joints = n_segments - 1        # angle scale between per-joint and tip angle

inertia = [0.5 * seg_mass * seg_radius^2,
           seg_mass * seg_length^2 / 12,
           seg_mass * seg_length^2 / 12]

# Non-bending DOF stay stiff/linear so a pure tip couple bends only the joints.
stiffness_axial = 1.0e5
stiffness_shear = 1.0e5
stiffness_torsion = 5.0e3
damping_trans = 50.0
damping_rot = 60.0             # transient-only; does not affect the static fit

# ----- measured data (virtual): root torque vs tip angle -----
# Replace with real cantilever-bend measurements. Gently stiffening, 0–40 N·m
# over 0–60°.
angle_deg = [0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0]
torque_data = [0.0, 5.3, 10.8, 16.8, 23.5, 31.2, 40.0]   # [N·m]
angle_rad = deg2rad.(angle_deg)

"""
    BendingPolynomial(coeffs, angle_scale)

Callable bending-stiffness law. Evaluates the odd polynomial
`p(γ) = c₁·γ + c₃·γ³ + c₅·γ⁵` (fitted to root-torque-vs-tip-angle data) at
`γ = angle_scale·joint_angle`, returning the restoring moment [N·m] for a single
joint. `angle_scale = 1` recovers the root-torque curve itself.
"""
struct BendingPolynomial
    coeffs::NTuple{3, Float64}
    angle_scale::Float64
end

function (law::BendingPolynomial)(joint_angle)
    γ = law.angle_scale * joint_angle
    c1, c3, c5 = law.coeffs
    return c1 * γ + c3 * γ^3 + c5 * γ^5
end

"""
    fit_bending_polynomial(angle, torque) -> NTuple{3,Float64}

Least-squares fit of the odd polynomial `c₁·γ + c₃·γ³ + c₅·γ⁵` (through the
origin: zero angle ⇒ zero moment, symmetric beam) to torque(angle) data, with
`angle` in radians.
"""
function fit_bending_polynomial(angle, torque)
    basis = hcat(angle, angle .^ 3, angle .^ 5)
    coeffs = basis \ torque
    return (coeffs[1], coeffs[2], coeffs[3])
end

coeffs = fit_bending_polynomial(angle_rad, torque_data)
root_law = BendingPolynomial(coeffs, 1.0)          # M = p(γ_tip)
joint_law = BendingPolynomial(coeffs, n_joints)    # per-joint moment(Δθ)
println("Fitted p(γ) = $(round(coeffs[1]; digits=3))·γ + " *
        "$(round(coeffs[2]; digits=3))·γ³ + $(round(coeffs[3]; digits=3))·γ⁵")

# ----- settings: reuse the committed beam data dir, no gravity -----
set_data_path(joinpath(dirname(@__DIR__), "data", "beam"))
set = Settings("system.yaml")
set.g_earth = 0.0                       # pure static bend test, no sag
set.physical_model = "beam_stiffness_fit"

# ----- build the beam with the fitted bending law -----
bodies = RigidBody[]
for i in 1:n_segments
    origin_x = (i - 0.5) * seg_length
    push!(bodies, RigidBody(Symbol("seg_$i");
        mass=seg_mass, inertia_principal=inertia,
        pos=[origin_x, 0.0, 0.0],
        fixed=(i == 1)))
end

joints = ElasticJoint[]
for i in 1:n_joints
    push!(joints, ElasticJoint(Symbol("joint_$i"),
        Symbol("seg_$i"), Symbol("seg_$(i+1)");
        anchor_a=[seg_length / 2, 0.0, 0.0],
        anchor_b=[-seg_length / 2, 0.0, 0.0],
        stiffness_axial, stiffness_shear, stiffness_torsion,
        stiffness_bending=joint_law,       # nonlinear, fitted
        damping_trans, damping_rot))
end

sys = SystemStructure("beam_stiffness_fit", set;
    rigid_bodies=bodies, elastic_joints=joints)
sam = SymbolicAWEModel(set, sys)
init!(sam)

# ----- validate: sweep the applied couple, settle, read tip angle -----
# Each torque level is held until the tip angle is stable for `stable_needed`
# consecutive steps — a static-equilibrium check robust to the beam's transient
# oscillation (an instantaneous zero-velocity is not equilibrium).
dt = 0.05
settle_max_time = 30.0
settle_steps = round(Int, settle_max_time / dt)
angle_tol = 1.0e-4               # [deg] per-step tip-angle change at equilibrium
stable_needed = 40

tip = sam.sys_struct.rigid_bodies[Symbol("seg_$n_segments")]

"""
    tip_angle_deg(body) -> Float64

Tip bending angle [deg]: the pitch of the body's x-axis (the beam axis) about the
world y-axis, read from its body→world rotation.
"""
function tip_angle_deg(body)
    R = SymbolicAWEModels.quaternion_to_rotation_matrix(body.Q_b_to_w)
    return rad2deg(atan(-R[3, 1], R[1, 1]))
end

logger = Logger(sam, length(torque_data) * settle_steps + 1)
sys_state = SysState(sam)
model_angle_deg = Float64[]
sim_time = 0.0

for target_torque in torque_data
    tip.ext_moment_b .= [0.0, target_torque, 0.0]
    prev_angle = NaN
    stable = 0
    for _ in 1:settle_steps
        next_step!(sam; dt, vsm_interval=0)
        global sim_time += dt
        update_sys_state!(sys_state, sam)
        sys_state.time = sim_time
        log!(logger, sys_state)
        angle = tip_angle_deg(tip)
        if abs(angle - prev_angle) < angle_tol
            stable += 1
            stable >= stable_needed && break
        else
            stable = 0
        end
        prev_angle = angle
    end
    push!(model_angle_deg, tip_angle_deg(tip))
end

# ----- report -----
println("\n torque[N·m]  data γ[deg]  model γ[deg]   error[deg]")
for i in eachindex(torque_data)
    err = model_angle_deg[i] - angle_deg[i]
    println(lpad(round(torque_data[i]; digits=1), 9), "  ",
            lpad(round(angle_deg[i]; digits=2), 10), "  ",
            lpad(round(model_angle_deg[i]; digits=2), 11), "  ",
            lpad(round(err; digits=3), 11))
end
rms = sqrt(sum((model_angle_deg .- angle_deg) .^ 2) / length(angle_deg))
println("RMS tip-angle error: $(round(rms; digits=3)) deg")

# ----- plot the fit and the model equilibria -----
fig = Figure()
ax = Axis(fig[1, 1]; xlabel="tip angle γ [deg]", ylabel="root torque [N·m]",
    title="Bending-stiffness fit: $n_segments-segment beam")
γ_curve = range(0, 60; length=100)
lines!(ax, γ_curve, root_law.(deg2rad.(γ_curve)); label="fitted p(γ)")
scatter!(ax, angle_deg, torque_data; markersize=12,
    label="virtual measurements")
scatter!(ax, model_angle_deg, torque_data; marker=:xcross, markersize=16,
    color=:red, label="beam model equilibria")
axislegend(ax; position=:lt)
display(GLMakie.Screen(), fig)

# ----- replay the bend sweep -----
save_log(logger, "beam_stiffness_fit")
syslog = load_log("beam_stiffness_fit")
scene = replay(syslog, sam.sys_struct; vector_scale=0.3)
display(GLMakie.Screen(), scene)

nothing
