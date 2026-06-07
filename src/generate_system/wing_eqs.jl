# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Wing rigid body dynamics equation generation

"""
    wing_eqs!(s, eqs, psys, defaults; kwargs...)

Generate the differential equations for the wing's
rigid body dynamics.

For RIGID_DYNAMICS wings:
- ODE state: `com_w`, `com_vel`, `Q_p_to_w`, `ω_p` (principal frame)
- Euler rotation equations in principal frame (diagonal I)
- Newton's 2nd law for COM translation
- Body frame output (`R_b_to_w`, `wing_pos`, `ω_b`) computed
  algebraically via `R_b_to_w` = `R_p_to_w` * `R_b_to_p` (constant)

For PARTICLE_DYNAMICS wings:
- No rigid body dynamics (handled by DYNAMIC points)
- `R_b_to_w` from structural ref points
- Principal frame variables set to zero/aliases
"""
function wing_eqs!(
    s, eqs, psys, defaults;
    tether_wing_force, tether_wing_moment,
    aero_force_b, aero_moment_b,
    ω_b, α_b, R_b_to_w, R_p_to_w,
    wing_pos, wing_vel, wing_acc,
    com_w, com_vel, com_acc, Q_p_to_w, ω_p, α_p,
    fix_wing, pos, vel, acc
)
    wings = s.sys_struct.wings

    @variables begin
        # Principal frame intermediates
        α_p_damped(t)[1:3, eachindex(wings)]
        ω_p_stable(t)[1:3, eachindex(wings)]
        Q_p_vel(t)[1:4, eachindex(wings)]
        moment_p(t)[1:3, eachindex(wings)]
        # Body frame quaternion (algebraic for all)
        Q_b_to_w(t)[1:4, eachindex(wings)]
        # Forces, moments, mass
        moment_tether_wing(t)[1:3, eachindex(wings)]
        force_tether_wing(t)[1:3, eachindex(wings)]
        wing_mass(t)[eachindex(wings)]
        fix_wing_sphere(t)[eachindex(wings)]
    end
    moment_tether_wing = collect(moment_tether_wing)

    # Skew-symmetric matrix for quaternion kinematics
    Ω(ω) = [
        0 -ω[1] -ω[2] -ω[3]
        ω[1] 0 ω[3] -ω[2]
        ω[2] -ω[3] 0 ω[1]
        ω[3] ω[2] -ω[1] 0
    ]

    # Weighted ref point position (symbolic)
    function get_ref_position(
        pos, ref_pt::WeightedRefPoints
    )
        w1 = ref_pt.weights[1]
        id1 = ref_pt.ids[1]
        if length(ref_pt.ids) == 1
            return pos[:, id1]
        end
        # Use element-wise access to avoid
        # symbolic slice scalarization issues
        result = [
            sum(ref_pt.weights[i] *
                pos[j, ref_pt.ids[i]]
                for i in eachindex(ref_pt.ids))
            for j in 1:3
        ]
        return result
    end

    for wing in wings
        # ============= PARTICLE_DYNAMICS WINGS ============= #
        if wing.dynamics_type == PARTICLE_DYNAMICS
            z_p1, z_p2 = wing.z_ref_points
            y_p1, y_p2 = wing.y_ref_points
            pos_z1 = get_ref_position(pos, z_p1)
            pos_z2 = get_ref_position(pos, z_p2)
            pos_y1 = get_ref_position(pos, y_p1)
            pos_y2 = get_ref_position(pos, y_p2)

            R_wing = R_b_to_w[:, :, wing.idx]
            eqs = [
                eqs
                # R_b_to_w from structural ref points
                R_b_to_w[:, 3, wing.idx] ~
                    smooth_normalize(pos_z2 - pos_z1)
                R_b_to_w[:, 1, wing.idx] ~ smooth_normalize(
                    smooth_normalize(pos_y2 - pos_y1) ×
                    R_b_to_w[:, 3, wing.idx])
                R_b_to_w[:, 2, wing.idx] ~
                    R_b_to_w[:, 3, wing.idx] ×
                    R_b_to_w[:, 1, wing.idx]

                # Body frame output from origin ref points
                wing_pos[:, wing.idx] ~
                    get_ref_position(pos, wing.origin)
                wing_vel[:, wing.idx] ~
                    get_ref_position(vel, wing.origin)
                wing_acc[:, wing.idx] ~
                    get_ref_position(acc, wing.origin)

                # Q_b_to_w from R_b_to_w
                Q_b_to_w[1, wing.idx] ~
                    rotation_matrix_to_quaternion_w(
                        R_wing)
                Q_b_to_w[2, wing.idx] ~
                    rotation_matrix_to_quaternion_x(
                        R_wing)
                Q_b_to_w[3, wing.idx] ~
                    rotation_matrix_to_quaternion_y(
                        R_wing)
                Q_b_to_w[4, wing.idx] ~
                    rotation_matrix_to_quaternion_z(
                        R_wing)

                # Body frame angular state (zero for
                # PARTICLE_DYNAMICS — no rigid body rotation)
                ω_b[:, wing.idx] ~ zeros(3)
                α_b[:, wing.idx] ~ zeros(3)

                # Principal frame aliases/zeros
                [R_p_to_w[:, i, wing.idx] ~
                    R_b_to_w[:, i, wing.idx] for i = 1:3]
                com_w[:, wing.idx] ~
                    wing_pos[:, wing.idx]
                com_vel[:, wing.idx] ~
                    wing_vel[:, wing.idx]
                com_acc[:, wing.idx] ~
                    wing_acc[:, wing.idx]
                Q_p_to_w[:, wing.idx] ~
                    Q_b_to_w[:, wing.idx]
                ω_p[:, wing.idx] ~ zeros(3)
                α_p[:, wing.idx] ~ zeros(3)

                # Zero intermediates
                Q_p_vel[:, wing.idx] ~ zeros(4)
                moment_p[:, wing.idx] ~ zeros(3)
                α_p_damped[:, wing.idx] ~ zeros(3)
                ω_p_stable[:, wing.idx] ~ zeros(3)
                moment_tether_wing[:, wing.idx] ~
                    zeros(3)
                force_tether_wing[:, wing.idx] ~
                    zeros(3)
                wing_mass[wing.idx] ~ 0.0
                fix_wing_sphere[wing.idx] ~ false
            ]
            continue
        end

        # ============= RIGID_DYNAMICS WINGS ============= #

        I_p = get_inertia_principal(psys, wing.idx)
        com_axis = collect(
            smooth_normalize(com_w[:, wing.idx]))
        com_axis_p = collect(
            R_p_to_w[:, :, wing.idx]' * com_axis)

        eqs = [
            eqs
            fix_wing_sphere[wing.idx] ~
                get_fix_wing_sphere(psys, wing.idx)

            # === Principal frame quaternion
            #     kinematics ===
            # D(Q_p_to_w) = 0.5 * Ω(ω_p_stable) * Q_p_to_w
            [D(Q_p_to_w[i, wing.idx]) ~
                Q_p_vel[i, wing.idx] for i = 1:4]
            [Q_p_vel[i, wing.idx] ~ 0.5 * sum(
                Ω(ω_p_stable[:, wing.idx])[i, j] *
                Q_p_to_w[j, wing.idx]
                for j = 1:4) for i = 1:4]

            # Constrain ω for spherical joint
            ω_p_stable[:, wing.idx] ~ ifelse.(
                fix_wing == true,
                zeros(3),
                ifelse.(
                    fix_wing_sphere[wing.idx] == true,
                    ω_p[:, wing.idx] -
                    (ω_p[:, wing.idx] ⋅ com_axis_p) *
                    com_axis_p,
                    ω_p[:, wing.idx],
                ),
            )

            # Constrain angular acceleration
            D(ω_p[:, wing.idx]) ~ ifelse.(
                fix_wing == true,
                zeros(3),
                ifelse.(
                    fix_wing_sphere[wing.idx] == true,
                    α_p_damped[:, wing.idx] -
                    (α_p_damped[:, wing.idx] ⋅
                     com_axis_p) * com_axis_p,
                    α_p_damped[:, wing.idx],
                ),
            )

            # Damping in principal frame
            α_p_damped[:, wing.idx] ~ [
                α_p[1, wing.idx] -
                    get_angular_damping(
                        psys, wing.idx) *
                    ω_p[1, wing.idx],
                α_p[2, wing.idx] -
                    (get_y_damping(psys, wing.idx) +
                     get_angular_damping(
                         psys, wing.idx)) *
                    ω_p[2, wing.idx],
                α_p[3, wing.idx] +
                    get_z_disturb(psys, wing.idx) -
                    get_angular_damping(
                        psys, wing.idx) *
                    ω_p[3, wing.idx],
            ]

            # R_p_to_w from Q_p_to_w
            [R_p_to_w[:, i, wing.idx] ~
                quaternion_to_rotation_matrix(
                    Q_p_to_w[:, wing.idx])[:, i]
                for i = 1:3]

            # === Euler equations (principal frame,
            #     diagonal inertia) ===
            α_p[1, wing.idx] ~ (
                moment_p[1, wing.idx] +
                (I_p[2] - I_p[3]) *
                ω_p[2, wing.idx] *
                ω_p[3, wing.idx]) / I_p[1]
            α_p[2, wing.idx] ~ (
                moment_p[2, wing.idx] +
                (I_p[3] - I_p[1]) *
                ω_p[3, wing.idx] *
                ω_p[1, wing.idx]) / I_p[2]
            α_p[3, wing.idx] ~ (
                moment_p[3, wing.idx] +
                (I_p[1] - I_p[2]) *
                ω_p[1, wing.idx] *
                ω_p[2, wing.idx]) / I_p[3]

            # === Total moment about COM ===
            moment_tether_wing[:, wing.idx] ~
                tether_wing_moment[:, wing.idx]
        ]

        # Aero moment transport: body origin → COM
        com_off_b = collect(
            get_com_offset_b(psys, wing.idx))
        aero_moment_com_b =
            aero_moment_b[:, wing.idx] .+
            (aero_force_b[:, wing.idx] × com_off_b)
        # Total moment in world frame about COM
        moment_w =
            collect(R_b_to_w[:, :, wing.idx]) *
            aero_moment_com_b .+
            moment_tether_wing[:, wing.idx]

        eqs = [
            eqs
            # Rotate to principal frame
            moment_p[:, wing.idx] ~
                collect(R_p_to_w[:, :, wing.idx])' *
                moment_w

            # === Translational dynamics ===
            D(com_w[:, wing.idx]) ~ ifelse.(
                fix_wing == true,
                zeros(3),
                ifelse.(
                    fix_wing_sphere[wing.idx] == true,
                    (com_vel[:, wing.idx] ⋅ com_axis) *
                    com_axis,
                    com_vel[:, wing.idx],
                ),
            )
            D(com_vel[:, wing.idx]) ~ ifelse.(
                fix_wing == true,
                zeros(3),
                ifelse.(
                    fix_wing_sphere[wing.idx] == true,
                    (com_acc[:, wing.idx] ⋅ com_axis) *
                    com_axis,
                    com_acc[:, wing.idx],
                ),
            )
            wing_mass[wing.idx] ~
                get_wing_mass(psys, wing.idx)
            force_tether_wing[:, wing.idx] ~
                tether_wing_force[:, wing.idx]
            com_acc[:, wing.idx] ~ (
                force_tether_wing[:, wing.idx] .+
                collect(R_b_to_w[:, :, wing.idx]) *
                aero_force_b[:, wing.idx]
            ) / wing_mass[wing.idx]
        ]

        # === Body frame output (unified) ===
        # R_b_to_w = R_p_to_w * R_b_to_p (constant R_b_to_p)
        # When R_b_to_p = I (no ref points): R_b_to_w = R_p_to_w
        R_b_to_p = get_R_b_to_p(psys, wing.idx)
        R_p_to_w_mat = collect(R_p_to_w[:, :, wing.idx])
        R_b_to_w_mat = R_p_to_w_mat * R_b_to_p

        # COM→origin in world frame
        r_w = -(R_b_to_w_mat * com_off_b)
        ω_w = R_p_to_w_mat * ω_p[:, wing.idx]

        eqs = [
            eqs
            # R_b_to_w from constant body-principal rotation
            [R_b_to_w[:, i, wing.idx] ~
                R_b_to_w_mat[:, i] for i = 1:3]

            # wing_pos = com_w + R_b_to_w * (-com_offset_b)
            wing_pos[:, wing.idx] ~
                com_w[:, wing.idx] .+ r_w
            # Rigid body kinematics for origin vel/acc
            wing_vel[:, wing.idx] ~
                com_vel[:, wing.idx] .+
                (ω_w × r_w)
            wing_acc[:, wing.idx] ~
                com_acc[:, wing.idx] .+
                ((R_p_to_w_mat *
                  α_p[:, wing.idx]) × r_w) .+
                (ω_w × (ω_w × r_w))

            # ω_b = R_b_to_p' * ω_p (constant rotation)
            ω_b[:, wing.idx] ~
                R_b_to_p' * ω_p[:, wing.idx]
            α_b[:, wing.idx] ~
                R_b_to_p' * α_p[:, wing.idx]
        ]

        # Q_b_to_w from R_b_to_w (both cases)
        R_wing = R_b_to_w[:, :, wing.idx]
        eqs = [
            eqs
            Q_b_to_w[1, wing.idx] ~
                rotation_matrix_to_quaternion_w(
                    R_wing)
            Q_b_to_w[2, wing.idx] ~
                rotation_matrix_to_quaternion_x(
                    R_wing)
            Q_b_to_w[3, wing.idx] ~
                rotation_matrix_to_quaternion_y(
                    R_wing)
            Q_b_to_w[4, wing.idx] ~
                rotation_matrix_to_quaternion_z(
                    R_wing)
        ]

        # Defaults: principal frame ODE state
        defaults = [
            defaults
            [com_w[i, wing.idx] =>
                get_com_w(psys, wing.idx)[i]
                for i = 1:3]
            [com_vel[i, wing.idx] =>
                get_com_vel(psys, wing.idx)[i]
                for i = 1:3]
            [Q_p_to_w[i, wing.idx] =>
                get_Q_p_to_w(psys, wing.idx)[i]
                for i = 1:4]
            [ω_p[i, wing.idx] =>
                get_ω_p(psys, wing.idx)[i]
                for i = 1:3]
        ]
    end

    return eqs, defaults
end
