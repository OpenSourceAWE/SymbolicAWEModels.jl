# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# 6-DOF elastic joint equation generation.

"""
    joint_stiffness_term(joint, params, kind, Δ)

Restoring force/moment for one joint DOF, read as a flat parameter: a `Real`
stiffness is a numeric scalar param (`k·Δ`); an interpolation is a callable param
applied as `k(Δ)`. `kind`: 1=axial, 2=shear, 3=torsion, 4=bending.
"""
function joint_stiffness_term(joint, params, kind::Int, Δ)
    field = (:stiffness_axial, :stiffness_shear,
             :stiffness_torsion, :stiffness_bending)[kind]
    k = getproperty(params.elastic_joints[joint.idx], field)
    return getfield(joint, field) isa Real ? k * Δ : k(Δ)
end

"""
    joint_eqs!(eqs, elastic_joints, params; kwargs...)

For each `ElasticJoint`, compute the restoring wrench from the relative pose of
the two anchors (in body A's frame) and accumulate it — equal and opposite —
into `body_force`/`body_moment` (the same accumulators `body_eqs!` reads). The
relative rotation uses the small-angle vector extraction, exact for the small
per-joint rotations of a stiff chain.
"""
function joint_eqs!(
    eqs, elastic_joints, params;
    body_force, body_moment,
    body_com_w, body_pos_w, body_com_vel, body_ω_b, body_R_b_to_w,
)
    @variables begin
        joint_force_w(t)[1:3, eachindex(elastic_joints)]
        joint_torque_w(t)[1:3, eachindex(elastic_joints)]
    end

    for joint in elastic_joints
        j = joint.idx
        a = joint.body_a_idx
        b = joint.body_b_idx
        R_a = collect(body_R_b_to_w[:, :, a])
        R_b = collect(body_R_b_to_w[:, :, b])
        anchor_a = collect(params.elastic_joints[j].anchor_a_b)
        anchor_b = collect(params.elastic_joints[j].anchor_b_b)
        pos_a = collect(body_pos_w[:, a])
        pos_b = collect(body_pos_w[:, b])
        com_a = collect(body_com_w[:, a])
        com_b = collect(body_com_w[:, b])

        # Anchor world positions (body_pos_w is the body origin).
        pos_anchor_a = pos_a .+ R_a * anchor_a
        pos_anchor_b = pos_b .+ R_b * anchor_b

        # Rest references captured at init: the as-placed geometry is unstrained.
        rest_offset = collect(params.elastic_joints[j].rest_offset_a)
        R_rel0 = collect(params.elastic_joints[j].R_rel0)

        # Relative displacement, body A frame: [axial, shear_y, shear_z].
        Δr_a = R_a' * (pos_anchor_b .- pos_anchor_a) .- rest_offset

        # Relative rotation, body A frame: [torsion, bend_y, bend_z].
        R_rel = R_rel0' * (R_a' * R_b)
        Δθ_a = [
            0.5 * (R_rel[3, 2] - R_rel[2, 3]),
            0.5 * (R_rel[1, 3] - R_rel[3, 1]),
            0.5 * (R_rel[2, 1] - R_rel[1, 2]),
        ]

        # Anchor velocities: v = com_vel + ω_w × (anchor − com).
        ω_a_w = R_a * collect(body_ω_b[:, a])
        ω_b_w = R_b * collect(body_ω_b[:, b])
        vel_anchor_a = collect(body_com_vel[:, a]) .+
            (ω_a_w × (pos_anchor_a .- com_a))
        vel_anchor_b = collect(body_com_vel[:, b]) .+
            (ω_b_w × (pos_anchor_b .- com_b))
        Δv_a = R_a' * (vel_anchor_b .- vel_anchor_a)
        Δω_a = R_a' * (ω_b_w .- ω_a_w)

        damp_trans = params.elastic_joints[j].damping_trans
        damp_rot = params.elastic_joints[j].damping_rot

        # Built element-wise: symbolic-array broadcasting is fragile here.
        force_a = [
            -joint_stiffness_term(joint, params, 1, Δr_a[1]) - damp_trans * Δv_a[1],
            -joint_stiffness_term(joint, params, 2, Δr_a[2]) - damp_trans * Δv_a[2],
            -joint_stiffness_term(joint, params, 2, Δr_a[3]) - damp_trans * Δv_a[3],
        ]
        torque_a = [
            -joint_stiffness_term(joint, params, 3, Δθ_a[1]) - damp_rot * Δω_a[1],
            -joint_stiffness_term(joint, params, 4, Δθ_a[2]) - damp_rot * Δω_a[2],
            -joint_stiffness_term(joint, params, 4, Δθ_a[3]) - damp_rot * Δω_a[3],
        ]

        eqs = [
            eqs
            joint_force_w[:, j] ~ R_a * force_a
            joint_torque_w[:, j] ~ R_a * torque_a
        ]

        force_on_b = collect(joint_force_w[:, j])
        torque_on_b = collect(joint_torque_w[:, j])
        arm_b = pos_anchor_b .- com_b
        arm_a = pos_anchor_a .- com_a

        # Equal-and-opposite wrench, moments transported to each COM.
        moment_on_b = arm_b × force_on_b + torque_on_b
        moment_on_a = arm_a × (-force_on_b) - torque_on_b
        body_force[:, b] .+= force_on_b
        body_force[:, a] .-= force_on_b
        body_moment[:, b] .+= moment_on_b
        body_moment[:, a] .+= moment_on_a
    end
    return eqs
end
