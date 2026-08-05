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
        ex = elastic_joint_wrench(joint, params;
            force_w = joint_force_w[:, j], torque_w = joint_torque_w[:, j],
            pos_a = collect(body_pos_w[:, a]), R_a, com_a = collect(body_com_w[:, a]),
            com_vel_a = collect(body_com_vel[:, a]), omega_a_w = R_a * collect(body_ω_b[:, a]),
            pos_b = collect(body_pos_w[:, b]), R_b, com_b = collect(body_com_w[:, b]),
            com_vel_b = collect(body_com_vel[:, b]), omega_b_w = R_b * collect(body_ω_b[:, b]))
        eqs = [eqs; ex.tear_eqs]
        body_force[:, a] .+= ex.force_on_a
        body_force[:, b] .+= ex.force_on_b
        body_moment[:, a] .+= ex.moment_on_a
        body_moment[:, b] .+= ex.moment_on_b
    end
    return eqs
end
