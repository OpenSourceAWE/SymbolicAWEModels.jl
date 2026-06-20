# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Standalone rigid body dynamics equation generation. A thin wrapper that
# assembles the body's loads (gravity + settable external wrench) and delegates
# the 6-DOF integration to the shared `rigid_body_eqs!`.

"""
    body_eqs!(eqs, defaults, rigid_bodies, params, initial; kwargs...)

Generate the differential equations for each standalone `RigidBody`. Loads are
the accumulated joint wrench (`body_force`/`body_moment`, filled by `joint_eqs!`)
plus gravity (`-g·mass` at the COM, world frame) and the external wrench read
live from the struct (`ext_force_w` world, `ext_moment_b` body). Isotropic
angular damping is applied through the `d_ω_p` integration override.
"""
function body_eqs!(
    eqs, defaults, rigid_bodies, params, initial;
    body_force, body_moment,
    body_com_w, body_com_vel, body_com_acc, body_Q_p_to_w, body_ω_p, body_α_p,
    body_pos_w, body_vel_w, body_acc_w, body_ω_b, body_α_b, body_Q_b_to_w,
    body_R_b_to_w, body_R_p_to_w, body_moment_p, body_Q_p_vel,
)
    for rigid_body in rigid_bodies
        idx = rigid_body.idx
        mass = params.rigid_bodies[idx].mass

        gravity_w = Num[0, 0, -params.set.g_earth * mass]
        force_w = collect(body_force[:, idx]) .+
            collect(params.rigid_bodies[idx].ext_force_w) .+ gravity_w
        moment_w = collect(body_moment[:, idx]) .+
            collect(body_R_b_to_w[:, :, idx]) *
            collect(params.rigid_bodies[idx].ext_moment_b)

        # A fixed body freezes all DOF: zero every integrated derivative so the
        # state stays at its initial pose. Otherwise apply isotropic angular
        # damping in the principal frame.
        if rigid_body.fixed
            overrides = (ω_kinematic=zeros(3), d_ω_p=zeros(3),
                         d_com_w=zeros(3), d_com_vel=zeros(3))
        else
            d_ω_p = body_α_p[:, idx] .-
                params.rigid_bodies[idx].angular_damping * body_ω_p[:, idx]
            overrides = (; d_ω_p)
        end

        eqs, defaults = rigid_body_eqs!(
            eqs, defaults, idx;
            force_w, moment_w,
            inertia_p=params.rigid_bodies[idx].inertia_principal,
            mass,
            R_b_to_p=params.rigid_bodies[idx].R_b_to_p,
            com_offset_b=params.rigid_bodies[idx].com_offset_b,
            com_w=body_com_w, com_vel=body_com_vel,
            Q_p_to_w=body_Q_p_to_w, ω_p=body_ω_p,
            com_acc=body_com_acc, α_p=body_α_p, R_p_to_w=body_R_p_to_w,
            moment_p=body_moment_p, Q_p_vel=body_Q_p_vel,
            R_b_to_w=body_R_b_to_w,
            wing_pos=body_pos_w, wing_vel=body_vel_w, wing_acc=body_acc_w,
            ω_b=body_ω_b, α_b=body_α_b, Q_b_to_w=body_Q_b_to_w,
            initial_com_w=initial.rigid_bodies[idx].com_w,
            initial_com_vel=initial.rigid_bodies[idx].com_vel,
            initial_Q_p_to_w=initial.rigid_bodies[idx].Q_p_to_w,
            initial_ω_p=initial.rigid_bodies[idx].ω_p,
            overrides...,
        )
    end
    return eqs, defaults
end
