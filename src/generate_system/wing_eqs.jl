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
            q_wing = rotation_matrix_to_quaternion(R_wing)
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

                # Q_b_to_w from R_b_to_w (one symbolic conversion, CSE-shared)
                Q_b_to_w[1, wing.idx] ~ q_wing[1]
                Q_b_to_w[2, wing.idx] ~ q_wing[2]
                Q_b_to_w[3, wing.idx] ~ q_wing[3]
                Q_b_to_w[4, wing.idx] ~ q_wing[4]

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

        idx = wing.idx
        com_axis = collect(smooth_normalize(com_w[:, idx]))
        com_axis_p = collect(R_p_to_w[:, :, idx]' * com_axis)
        sphere = fix_wing_sphere[idx]

        # Wing-specific intermediates: pinning flag, damped angular
        # acceleration, tether loads, mass.
        eqs = [
            eqs
            fix_wing_sphere[idx] ~ params.wings[idx].fix_sphere

            ω_p_stable[:, idx] ~ ifelse.(
                fix_wing == true, zeros(3),
                ifelse.(sphere == true,
                    ω_p[:, idx] -
                    (ω_p[:, idx] ⋅ com_axis_p) * com_axis_p,
                    ω_p[:, idx]))

            α_p_damped[:, idx] ~ [
                α_p[1, idx] -
                    params.wings[idx].angular_damping * ω_p[1, idx],
                α_p[2, idx] -
                    (params.wings[idx].y_damping +
                     params.wings[idx].angular_damping) * ω_p[2, idx],
                α_p[3, idx] + params.wings[idx].z_disturb -
                    params.wings[idx].angular_damping * ω_p[3, idx],
            ]

            moment_tether_wing[:, idx] ~ tether_wing_moment[:, idx]
            force_tether_wing[:, idx] ~ tether_wing_force[:, idx]
            wing_mass[idx] ~ params.wings[idx].mass
        ]

        # Total force/moment at/about COM (world frame).
        # Aero moment transported body origin → COM.
        com_off_b = collect(params.wings[idx].com_offset_b)
        aero_moment_com_b = aero_moment_b[:, idx] .+
            (aero_force_b[:, idx] × com_off_b)
        moment_w = collect(R_b_to_w[:, :, idx]) * aero_moment_com_b .+
            moment_tether_wing[:, idx]
        force_w = force_tether_wing[:, idx] .+
            collect(R_b_to_w[:, :, idx]) * aero_force_b[:, idx]

        # Pinning constraints project the integrated derivatives:
        # fix_wing freezes the body, fix_wing_sphere keeps it on a sphere.
        d_ω_p = ifelse.(
            fix_wing == true, zeros(3),
            ifelse.(sphere == true,
                α_p_damped[:, idx] -
                (α_p_damped[:, idx] ⋅ com_axis_p) * com_axis_p,
                α_p_damped[:, idx]))
        d_com_w = ifelse.(
            fix_wing == true, zeros(3),
            ifelse.(sphere == true,
                (com_vel[:, idx] ⋅ com_axis) * com_axis,
                com_vel[:, idx]))
        d_com_vel = ifelse.(
            fix_wing == true, zeros(3),
            ifelse.(sphere == true,
                (com_acc[:, idx] ⋅ com_axis) * com_axis,
                com_acc[:, idx]))

        eqs, defaults = rigid_body_eqs!(
            eqs, defaults, idx;
            force_w, moment_w,
            inertia_p=params.wings[idx].inertia_principal,
            mass=wing_mass[idx],
            R_b_to_p=params.wings[idx].R_b_to_p,
            com_offset_b=com_off_b,
            com_w, com_vel, Q_p_to_w, ω_p,
            com_acc, α_p, R_p_to_w, moment_p, Q_p_vel,
            R_b_to_w, wing_pos, wing_vel, wing_acc, ω_b, α_b, Q_b_to_w,
            initial_com_w=initial.wings[idx].com_w,
            initial_com_vel=initial.wings[idx].com_vel,
            initial_Q_p_to_w=initial.wings[idx].Q_p_to_w,
            initial_ω_p=initial.wings[idx].ω_p,
            ω_kinematic=ω_p_stable[:, idx],
            d_ω_p, d_com_w, d_com_vel,
        )
    end

    return eqs, defaults
end
