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
    com_w, com_vel, Q_p_to_w, ω_p,
    com_acc, α_p, R_p_to_w, moment_p, Q_p_vel,
    R_b_to_w, wing_pos, wing_vel, wing_acc, ω_b, α_b, Q_b_to_w,
    initial_com_w, initial_com_vel, initial_Q_p_to_w, initial_ω_p,
    ω_kinematic=nothing, d_ω_p=nothing, d_com_w=nothing, d_com_vel=nothing,
)
    # Skew-symmetric matrix for quaternion kinematics
    Ω(ω) = [
        0 -ω[1] -ω[2] -ω[3]
        ω[1] 0 ω[3] -ω[2]
        ω[2] -ω[3] 0 ω[1]
        ω[3] ω[2] -ω[1] 0
    ]

    ω_kin = ω_kinematic === nothing ? ω_p[:, idx] : ω_kinematic
    d_ω = d_ω_p === nothing ? α_p[:, idx] : d_ω_p
    d_cw = d_com_w === nothing ? com_vel[:, idx] : d_com_w
    d_cv = d_com_vel === nothing ? com_acc[:, idx] : d_com_vel

    com_off_b = collect(com_offset_b)

    eqs = [
        eqs
        # === Quaternion kinematics ===
        [D(Q_p_to_w[i, idx]) ~ Q_p_vel[i, idx] for i = 1:4]
        [Q_p_vel[i, idx] ~ 0.5 * sum(
            Ω(ω_kin)[i, j] * Q_p_to_w[j, idx] for j = 1:4)
            for i = 1:4]

        # === Angular acceleration integration ===
        D(ω_p[:, idx]) ~ d_ω

        # R_p_to_w from Q_p_to_w
        [R_p_to_w[:, i, idx] ~
            quaternion_to_rotation_matrix(
                Q_p_to_w[:, idx])[:, i] for i = 1:3]

        # === Euler equations (principal frame, diagonal inertia) ===
        α_p[1, idx] ~ (moment_p[1, idx] +
            (inertia_p[2] - inertia_p[3]) *
            ω_p[2, idx] * ω_p[3, idx]) / inertia_p[1]
        α_p[2, idx] ~ (moment_p[2, idx] +
            (inertia_p[3] - inertia_p[1]) *
            ω_p[3, idx] * ω_p[1, idx]) / inertia_p[2]
        α_p[3, idx] ~ (moment_p[3, idx] +
            (inertia_p[1] - inertia_p[2]) *
            ω_p[1, idx] * ω_p[2, idx]) / inertia_p[3]

        # Total moment rotated to principal frame
        moment_p[:, idx] ~
            collect(R_p_to_w[:, :, idx])' * moment_w

        # === Translational dynamics ===
        D(com_w[:, idx]) ~ d_cw
        D(com_vel[:, idx]) ~ d_cv
        com_acc[:, idx] ~ force_w / mass
    ]

    # === Body frame output ===
    R_p_to_w_mat = collect(R_p_to_w[:, :, idx])
    R_b_to_w_mat = R_p_to_w_mat * R_b_to_p
    r_w = -(R_b_to_w_mat * com_off_b)              # COM→origin, world
    ω_w = R_p_to_w_mat * ω_p[:, idx]

    eqs = [
        eqs
        [R_b_to_w[:, i, idx] ~ R_b_to_w_mat[:, i] for i = 1:3]

        # Rigid body kinematics for origin pos/vel/acc
        wing_pos[:, idx] ~ com_w[:, idx] .+ r_w
        wing_vel[:, idx] ~ com_vel[:, idx] .+ (ω_w × r_w)
        wing_acc[:, idx] ~ com_acc[:, idx] .+
            ((R_p_to_w_mat * α_p[:, idx]) × r_w) .+
            (ω_w × (ω_w × r_w))

        # ω_b = R_b_to_p' * ω_p (constant rotation)
        ω_b[:, idx] ~ R_b_to_p' * ω_p[:, idx]
        α_b[:, idx] ~ R_b_to_p' * α_p[:, idx]
    ]

    # Q_b_to_w from R_b_to_w (one symbolic conversion, CSE-shared)
    R_body = R_b_to_w[:, :, idx]
    q_body = rotation_matrix_to_quaternion(R_body)
    eqs = [
        eqs
        Q_b_to_w[1, idx] ~ q_body[1]
        Q_b_to_w[2, idx] ~ q_body[2]
        Q_b_to_w[3, idx] ~ q_body[3]
        Q_b_to_w[4, idx] ~ q_body[4]
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
