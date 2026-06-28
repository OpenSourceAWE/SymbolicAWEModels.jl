# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Wing rigid body dynamics equation generation

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
    # Weighted ref point position (symbolic)
    function get_ref_position(
        pos, ref_pt::WeightedRefPoints
    )
        id1 = ref_pt.ids[1]
        if length(ref_pt.ids) == 1
            return pos[:, id1]
        end
        # Element-wise access avoids symbolic slice scalarization issues.
        result = [
            sum(ref_pt.weights[i] *
                pos[j, ref_pt.ids[i]]
                for i in eachindex(ref_pt.ids))
            for j in 1:3
        ]
        return result
    end

    for wing in s.sys_struct.bodies
        # KINEMATIC wings: body frame fitted from ref points; DYNAMIC via body_eqs!.
        wing.type == KINEMATIC || continue
        z_p1, z_p2 = wing.z_ref_points
        y_p1, y_p2 = wing.y_ref_points
        pos_z1 = get_ref_position(pos, z_p1)
        pos_z2 = get_ref_position(pos, z_p2)
        pos_y1 = get_ref_position(pos, y_p1)
        pos_y2 = get_ref_position(pos, y_p2)

        R_wing = R_b_to_w[:, :, wing.idx]
        q_wing = rotation_matrix_to_quaternion(R_wing)
        eqs = [
            eqs
            # R_b_to_w from structural ref points
            R_b_to_w[:, 3, wing.idx] ~ smooth_normalize(pos_z2 - pos_z1)
            R_b_to_w[:, 1, wing.idx] ~ smooth_normalize(
                smooth_normalize(pos_y2 - pos_y1) × R_b_to_w[:, 3, wing.idx])
            R_b_to_w[:, 2, wing.idx] ~
                R_b_to_w[:, 3, wing.idx] × R_b_to_w[:, 1, wing.idx]
            # Body frame output from origin ref points
            wing_pos[:, wing.idx] ~ get_ref_position(pos, wing.origin)
            wing_vel[:, wing.idx] ~ get_ref_position(vel, wing.origin)
            wing_acc[:, wing.idx] ~ get_ref_position(acc, wing.origin)
            Q_b_to_w[1, wing.idx] ~ q_wing[1]
            Q_b_to_w[2, wing.idx] ~ q_wing[2]
            Q_b_to_w[3, wing.idx] ~ q_wing[3]
            Q_b_to_w[4, wing.idx] ~ q_wing[4]
            ω_b[:, wing.idx] ~ zeros(3)
            α_b[:, wing.idx] ~ zeros(3)
            # Principal frame aliases/zeros (no rigid-body rotation)
            [R_p_to_w[:, i, wing.idx] ~ R_b_to_w[:, i, wing.idx] for i = 1:3]
            com_w[:, wing.idx] ~ wing_pos[:, wing.idx]
            com_vel[:, wing.idx] ~ wing_vel[:, wing.idx]
            com_acc[:, wing.idx] ~ wing_acc[:, wing.idx]
            Q_p_to_w[:, wing.idx] ~ Q_b_to_w[:, wing.idx]
            ω_p[:, wing.idx] ~ zeros(3)
            α_p[:, wing.idx] ~ zeros(3)
            Q_p_vel[:, wing.idx] ~ zeros(4)
            moment_p[:, wing.idx] ~ zeros(3)
        ]
    end

    return eqs, defaults
end
