# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Generic rigid body 6-DOF dynamics equation generation.

"""
    rigid_body_eqs!(eqs, defaults, idx; kwargs...)

Append the 6-DOF rigid body equations for body `idx` to `eqs` and its
initial-condition `defaults`. Given a total `force_w` at the center of mass and
`moment_w` about it (both world frame), integrate quaternion attitude and COM
translation, and emit the body-frame output.

This generator knows nothing about aerodynamics, pinning constraints, or
damping. Those are the caller's concern: a caller imposes them by passing the
`ω_kinematic`/`d_ω_p`/`d_com_w`/`d_com_vel` integration overrides. With the
overrides left at their defaults the body integrates freely.

# State (principal frame, integrated)
`com_w`, `com_vel`, `Q_p_to_w`, `ω_p`.

# Required keyword arguments
- `force_w`, `moment_w`: length-3 `Num` vectors, total force at / moment about
  the COM in world frame.
- `inertia_p`: length-3 principal inertia; `mass`: scalar.
- `apparent_mass`: entrained air resisting acceleration without weight (default 0).
- `R_b_to_p`: constant body→principal rotation; `com_offset_b`: COM offset in
  the body frame (origin→COM).
- State / output array variables (indexed `[.., idx]` internally): `com_w`,
  `com_vel`, `Q_p_to_w`, `ω_p`, `com_acc`, `α_p`, `R_p_to_w`, `moment_p`,
  `Q_p_vel`, `R_b_to_w`, `wing_pos`, `wing_vel`, `wing_acc`, `ω_b`, `α_b`,
  `Q_b_to_w`.
- Initial conditions: `initial_com_w`, `initial_com_vel`, `initial_Q_p_to_w`,
  `initial_ω_p` — `initial.*` view paths bound to the integrated state.

# Optional integration overrides (default to the unconstrained body)
- `ω_kinematic`: angular velocity used in quaternion kinematics (default `ω_p`).
- `d_ω_p`: RHS of `D(ω_p)` (default `α_p`).
- `d_com_w`: RHS of `D(com_w)` (default `com_vel`).
- `d_com_vel`: RHS of `D(com_vel)` (default `com_acc`).
"""
function rigid_body_eqs!(
    eqs, defaults, idx;
    force_w, moment_w, inertia_p, mass, R_b_to_p, com_offset_b,
    apparent_mass=0.0,
    com_w, com_vel, Q_p_to_w, ω_p,
    com_acc, α_p, R_p_to_w, moment_p, Q_p_vel,
    R_b_to_w, wing_pos, wing_vel, wing_acc, ω_b, α_b, Q_b_to_w,
    initial_com_w, initial_com_vel, initial_Q_p_to_w, initial_ω_p,
    ω_kinematic=nothing, d_ω_p=nothing, d_com_w=nothing, d_com_vel=nothing,
)
    ex = rigid_body_pose_expressions(
        force_w, moment_w, inertia_p, mass, R_b_to_p, apparent_mass, com_offset_b,
        com_w[:, idx], com_vel[:, idx], Q_p_to_w[:, idx], ω_p[:, idx];
        ω_kinematic, d_ω_p, d_com_w, d_com_vel)

    eqs = [
        eqs
        # === Quaternion kinematics ===
        [D(Q_p_to_w[i, idx]) ~ Q_p_vel[i, idx] for i = 1:4]
        [Q_p_vel[i, idx] ~ ex.Q_p_vel[i] for i = 1:4]

        # === Angular acceleration integration ===
        D(ω_p[:, idx]) ~ ex.d_ω

        [R_p_to_w[:, i, idx] ~ ex.R_p_to_w[:, i] for i = 1:3]

        # === Euler equations (principal frame, diagonal inertia) ===
        [α_p[k, idx] ~ ex.α_p[k] for k = 1:3]
        moment_p[:, idx] ~ ex.moment_p

        # === Translational dynamics ===
        D(com_w[:, idx]) ~ ex.d_com_w
        D(com_vel[:, idx]) ~ ex.d_com_vel
        com_acc[:, idx] ~ ex.com_acc

        # === Body frame output ===
        [R_b_to_w[:, i, idx] ~ ex.R_b_to_w[:, i] for i = 1:3]
        wing_pos[:, idx] ~ ex.pos_w
        wing_vel[:, idx] ~ ex.vel_w
        wing_acc[:, idx] ~ ex.acc_w
        ω_b[:, idx] ~ ex.ω_b
        α_b[:, idx] ~ ex.α_b
        [Q_b_to_w[k, idx] ~ ex.Q_b_to_w[k] for k = 1:4]
    ]

    defaults = [
        defaults
        bind_initial!(initial_com_w, collect(com_w[:, idx]))
        bind_initial!(initial_com_vel, collect(com_vel[:, idx]))
        bind_initial!(initial_Q_p_to_w, collect(Q_p_to_w[:, idx]))
        bind_initial!(initial_ω_p, collect(ω_p[:, idx]))
    ]

    return eqs, defaults
end
