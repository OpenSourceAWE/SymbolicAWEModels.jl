# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Wing rigid body dynamics equation generation

"""The columns of the symbolic `3×n` matrix `matrix`, each read element by element so no
symbolic slice is scalarized, indexed as KiteGeometry reads point positions."""
struct PointColumns{M}
    matrix::M
end
Base.getindex(columns::PointColumns, i) = [columns.matrix[j, i] for j in 1:3]

"""
    wing_eqs!(s, eqs, defaults, params, initial; kwargs...)

Generate the differential equations for the wing's
rigid body dynamics.

For RIGID_DYNAMICS wings:
- ODE state: `com_w`, `com_vel`, `Q_p_to_w`, `ω_p` (principal frame)
- Wing-specific loads (aero transport, tether, damping) and pinning
  constraints are assembled here, then the generic 6-DOF integration is
  delegated to `rigid_body_eqs!`.

For PARTICLE_DYNAMICS wings:
- No rigid body dynamics (handled by DYNAMIC points)
- `R_b_to_w` from structural ref points
- Principal frame variables set to zero/aliases
"""
function wing_eqs!(
    s, eqs, defaults, params, initial;
    ω_b, α_b, R_b_to_w, R_p_to_w,
    wing_pos, wing_vel, wing_acc,
    com_w, com_vel, com_acc, Q_p_to_w, ω_p, α_p,
    Q_b_to_w, moment_p, Q_p_vel,
    fix_wing, pos, vel, acc
)
    get_ref_position(pos, refs) = KiteGeometry.ref_position(PointColumns(pos), refs)

    for body_idx in eachindex(s.definition.bodies)
        # KINEMATIC wings: body frame fitted from ref points; DYNAMIC via body_eqs!.
        body_type(s, body_idx) == KINEMATIC || continue
        wing_idx = body_wing(s, body_idx)
        frame = wing_frame(s, wing_idx)
        (z_p1, z_p2), (y_p1, y_p2) = frame.z_ref_points, frame.y_ref_points
        pos_z1 = collect(get_ref_position(pos, z_p1))
        pos_z2 = collect(get_ref_position(pos, z_p2))
        pos_y1 = collect(get_ref_position(pos, y_p1))
        pos_y2 = collect(get_ref_position(pos, y_p2))

        R_wing = R_b_to_w[:, :, body_idx]
        q_wing = rotation_matrix_to_quaternion(R_wing)
        xaxis, yaxis, zaxis = wing_frame_columns(pos_z1, pos_z2, pos_y1, pos_y2;
            torn_frame=R_wing)

        frame_axes = (collect(R_wing[:, 1]), collect(R_wing[:, 2]),
                      collect(R_wing[:, 3]))
        omega_b = fitted_frame_omega(pos_z1, pos_z2, pos_y1, pos_y2,
            (collect(get_ref_position(vel, z_p1)),
             collect(get_ref_position(vel, z_p2)),
             collect(get_ref_position(vel, y_p1)),
             collect(get_ref_position(vel, y_p2))), frame_axes)
        eqs = [
            eqs
            R_b_to_w[:, 3, body_idx] ~ zaxis
            R_b_to_w[:, 1, body_idx] ~ xaxis
            R_b_to_w[:, 2, body_idx] ~ yaxis
            # Body frame output from origin ref points
            wing_pos[:, body_idx] ~ get_ref_position(pos, wing_origin(s, wing_idx))
            wing_vel[:, body_idx] ~ get_ref_position(vel, wing_origin(s, wing_idx))
            wing_acc[:, body_idx] ~ get_ref_position(acc, wing_origin(s, wing_idx))
            Q_b_to_w[1, body_idx] ~ q_wing[1]
            Q_b_to_w[2, body_idx] ~ q_wing[2]
            Q_b_to_w[3, body_idx] ~ q_wing[3]
            Q_b_to_w[4, body_idx] ~ q_wing[4]
            ω_b[:, body_idx] ~ omega_b
            α_b[:, body_idx] ~ zeros(3)
            # Principal frame aliases/zeros (no rigid-body rotation)
            [R_p_to_w[:, i, body_idx] ~ R_b_to_w[:, i, body_idx] for i = 1:3]
            com_w[:, body_idx] ~ wing_pos[:, body_idx]
            com_vel[:, body_idx] ~ wing_vel[:, body_idx]
            com_acc[:, body_idx] ~ wing_acc[:, body_idx]
            Q_p_to_w[:, body_idx] ~ Q_b_to_w[:, body_idx]
            ω_p[:, body_idx] ~ zeros(3)
            α_p[:, body_idx] ~ zeros(3)
            Q_p_vel[:, body_idx] ~ zeros(4)
            moment_p[:, body_idx] ~ zeros(3)
        ]
    end

    return eqs, defaults
end
