# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Body loads assembly + integration overrides; delegates to `rigid_body_eqs!`.

"""
    body_eqs!(eqs, defaults, bodies, params, initial; kwargs...)

Generate the differential equations for each plain `Body` (no aero). Loads are
the accumulated joint wrench (`body_force`/`body_moment`, filled by `joint_eqs!`)
plus gravity (`-g·mass` at the COM) and the external wrench (`ext_force_w` world,
`ext_force_b`/`ext_moment_b` body). `STATIC` bodies are frozen; `fix_sphere`
confines the COM to a sphere about the world origin; `damping` is per-axis
angular damping.
"""
function body_eqs!(
    eqs, defaults, bodies, params, initial;
    body_force, body_moment,
    body_com_w, body_com_vel, body_com_acc, body_Q_p_to_w, body_ω_p, body_α_p,
    body_pos_w, body_vel_w, body_acc_w, body_ω_b, body_α_b, body_Q_b_to_w,
    body_R_b_to_w, body_R_p_to_w, body_moment_p, body_Q_p_vel,
)
    for rigid_body in bodies
        idx = rigid_body.idx
        # Bind body-frame outputs for every body (mtkcompile may integrate them, not principal).
        defaults = [
            defaults
            bind_initial!(initial.bodies[idx].pos_w, collect(body_pos_w[:, idx]))
            bind_initial!(initial.bodies[idx].vel_w, collect(body_vel_w[:, idx]))
            bind_initial!(initial.bodies[idx].R_b_to_w, body_R_b_to_w[:, :, idx])
        ]
        # KINEMATIC bodies (particle wings) are pose-fitted by wing_eqs!, not here.
        # Their slaved COM can still survive as a tearing iteration variable, so it
        # needs a u0 guess (build_initializeprob=false runs no DAE init to fill it).
        if rigid_body.type == KINEMATIC
            defaults = [
                defaults
                bind_initial!(initial.bodies[idx].com_w, collect(body_com_w[:, idx]))
                bind_initial!(initial.bodies[idx].com_vel, collect(body_com_vel[:, idx]))
            ]
            continue
        end
        mass = params.bodies[idx].mass
        R_b_to_w = collect(body_R_b_to_w[:, :, idx])

        # Loads at / about the COM (world frame).
        gravity_w = Num[0, 0, -params.set.g_earth * mass]
        force_w = collect(body_force[:, idx]) .+ gravity_w .+
            collect(params.bodies[idx].ext_force_w) .+
            R_b_to_w * collect(params.bodies[idx].ext_force_b)
        moment_w = collect(body_moment[:, idx]) .+
            R_b_to_w * collect(params.bodies[idx].ext_moment_b)

        (; ω_kinematic, d_ω_p, d_com_w, d_com_vel) = body_integration(
            params, idx, body_com_w[:, idx], body_com_vel[:, idx],
            body_ω_p[:, idx], body_α_p[:, idx], body_com_acc[:, idx],
            collect(body_R_p_to_w[:, :, idx]);
            frozen=(rigid_body.type == STATIC))

        eqs, defaults = rigid_body_eqs!(
            eqs, defaults, idx;
            force_w, moment_w,
            inertia_p=params.bodies[idx].inertia_principal, mass,
            R_b_to_p=params.bodies[idx].R_b_to_p,
            com_offset_b=params.bodies[idx].com_offset_b,
            com_w=body_com_w, com_vel=body_com_vel,
            Q_p_to_w=body_Q_p_to_w, ω_p=body_ω_p,
            com_acc=body_com_acc, α_p=body_α_p, R_p_to_w=body_R_p_to_w,
            moment_p=body_moment_p, Q_p_vel=body_Q_p_vel,
            R_b_to_w=body_R_b_to_w,
            wing_pos=body_pos_w, wing_vel=body_vel_w, wing_acc=body_acc_w,
            ω_b=body_ω_b, α_b=body_α_b, Q_b_to_w=body_Q_b_to_w,
            initial_com_w=initial.bodies[idx].com_w,
            initial_com_vel=initial.bodies[idx].com_vel,
            initial_Q_p_to_w=initial.bodies[idx].Q_p_to_w,
            initial_ω_p=initial.bodies[idx].ω_p,
            ω_kinematic, d_ω_p, d_com_w, d_com_vel,
        )
    end
    return eqs, defaults
end
