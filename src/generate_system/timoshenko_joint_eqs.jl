# Copyright (c) 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Timoshenko joint equation generation. Each joint is the element of a 2-node
# corotational Timoshenko beam: it applies an equal-and-opposite restoring wrench
# to two bodies' load accumulators, from a consistent element stiffness that
# couples transverse displacement to rotation (transverse shear) — the
# distributed-compliance counterpart of joint_eqs.jl. A chain of bodies joined by
# these forms a beam.

"""
    timoshenko_rigidity(joint, params, field, arg)

Effective rigidity for one [`TimoshenkoJoint`](@ref) mode, read as a flat
parameter: a `Real` rigidity is a numeric scalar param used directly; a callable
is a callable param evaluated at the mode's strain/curvature `arg`. `field` is one
of `:EA`, `:GA`, `:GJ`, `:EIy`, `:EIz`.
"""
function timoshenko_rigidity(joint, params, field::Symbol, arg)
    rigidity = getproperty(params.timoshenko_joints[joint.idx], field)
    return getfield(joint, field) isa Real ? rigidity : rigidity(arg)
end

"""
    timoshenko_joint_eqs!(eqs, timoshenko_joints, params; kwargs...)

For each [`TimoshenkoJoint`](@ref), build a corotational element frame, extract the
small per-node deformations (axial stretch, chord-relative rotations) relative to
the rest geometry, evaluate the consistent Timoshenko stiffness (axial, torsion,
and two bending planes with the shear reduction `Φ = 12·EI/(k·GA·L²)`) — each
rigidity either constant or a callable of its strain/curvature ([`timoshenko_rigidity`](@ref)) — and
accumulate the restoring wrench — equal and opposite, transported to each COM —
into `body_force`/`body_moment` (the same accumulators `body_eqs!` reads).
Damping resists the axial stretch rate and the relative node spin.
"""
function timoshenko_joint_eqs!(
    eqs, timoshenko_joints, params;
    body_force, body_moment,
    body_com_w, body_pos_w, body_com_vel, body_ω_b, body_R_b_to_w,
)
    isempty(timoshenko_joints) && return eqs

    @variables begin
        timoshenko_force_a_w(t)[1:3, eachindex(timoshenko_joints)]
        timoshenko_force_b_w(t)[1:3, eachindex(timoshenko_joints)]
        timoshenko_moment_a_w(t)[1:3, eachindex(timoshenko_joints)]
        timoshenko_moment_b_w(t)[1:3, eachindex(timoshenko_joints)]
        # Torn intermediates so the reused frame subtree is not re-embedded/re-scalarized per nesting level.
        timoshenko_frame(t)[1:3, 1:3, eachindex(timoshenko_joints)]
        timoshenko_theta_a(t)[1:3, eachindex(timoshenko_joints)]
        timoshenko_theta_b(t)[1:3, eachindex(timoshenko_joints)]
    end

    for joint in timoshenko_joints
        j = joint.idx
        a = joint.body_a_idx
        b = joint.body_b_idx
        R_a = collect(body_R_b_to_w[:, :, a])
        R_b = collect(body_R_b_to_w[:, :, b])
        ex = timoshenko_element_wrench(joint, params;
            frame = timoshenko_frame[:, :, j],
            theta_a = timoshenko_theta_a[:, j], theta_b = timoshenko_theta_b[:, j],
            force_a = timoshenko_force_a_w[:, j], force_b = timoshenko_force_b_w[:, j],
            moment_a = timoshenko_moment_a_w[:, j], moment_b = timoshenko_moment_b_w[:, j],
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
