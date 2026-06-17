# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Standalone rigid body dynamics equation generation. A thin wrapper that
# assembles the body's loads (gravity + settable external wrench) and delegates
# the 6-DOF integration to the shared `rigid_body_eqs!`.

"""
    body_eqs!(eqs, defaults, psys, rigid_bodies; kwargs...)

Generate the differential equations for each standalone `RigidBody`. Loads are
the accumulated joint wrench (`body_force`/`body_moment`, filled by `joint_eqs!`)
plus gravity (`-g·mass` at the COM, world frame) and the external wrench read
live from the struct (`ext_force_w` world, `ext_moment_b` body). Isotropic
angular damping is applied through the `d_ω_p` integration override.
"""
function body_eqs!(
    eqs, defaults, psys, rigid_bodies;
    body_force, body_moment,
    body_com_w, body_com_vel, body_com_acc, body_Q_p_to_w, body_ω_p, body_α_p,
    body_pos_w, body_vel_w, body_acc_w, body_ω_b, body_α_b, body_Q_b_to_w,
    body_R_b_to_w, body_R_p_to_w, body_moment_p, body_Q_p_vel,
)
    for rigid_body in rigid_bodies
        idx = rigid_body.idx
        mass = get_body_mass(psys, idx)

        gravity_w = Num[0, 0, -get_g_earth(psys) * mass]
        force_w = collect(body_force[:, idx]) .+
            collect(get_body_ext_force_w(psys, idx)) .+ gravity_w
        moment_w = collect(body_moment[:, idx]) .+
            collect(body_R_b_to_w[:, :, idx]) *
            collect(get_body_ext_moment_b(psys, idx))

        # A fixed body freezes all DOF: zero every integrated derivative so the
        # state stays at its initial pose. Otherwise apply isotropic angular
        # damping in the principal frame.
        if rigid_body.fixed
            overrides = (ω_kinematic=zeros(3), d_ω_p=zeros(3),
                         d_com_w=zeros(3), d_com_vel=zeros(3))
        else
            d_ω_p = body_α_p[:, idx] .-
                get_body_angular_damping(psys, idx) * body_ω_p[:, idx]
            overrides = (; d_ω_p)
        end

        eqs, defaults = rigid_body_eqs!(
            eqs, defaults, idx;
            force_w, moment_w,
            inertia_p=get_body_inertia_principal(psys, idx),
            mass,
            R_b_to_p=get_body_R_b_to_p(psys, idx),
            com_offset_b=get_body_com_offset_b(psys, idx),
            com_w=body_com_w, com_vel=body_com_vel,
            Q_p_to_w=body_Q_p_to_w, ω_p=body_ω_p,
            com_acc=body_com_acc, α_p=body_α_p, R_p_to_w=body_R_p_to_w,
            moment_p=body_moment_p, Q_p_vel=body_Q_p_vel,
            R_b_to_w=body_R_b_to_w,
            wing_pos=body_pos_w, wing_vel=body_vel_w, wing_acc=body_acc_w,
            ω_b=body_ω_b, α_b=body_α_b, Q_b_to_w=body_Q_b_to_w,
            com_w_0=get_body_com_w(psys, idx),
            com_vel_0=get_body_com_vel(psys, idx),
            Q_p_to_w_0=get_body_Q_p_to_w(psys, idx),
            ω_p_0=get_body_ω_p(psys, idx),
            overrides...,
        )
    end
    return eqs, defaults
end
