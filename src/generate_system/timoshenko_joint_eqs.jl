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
Damping resists the relative node velocity and spin.
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
    end

    for joint in timoshenko_joints
        j = joint.idx
        a = joint.body_a_idx
        b = joint.body_b_idx
        R_a = collect(body_R_b_to_w[:, :, a])
        R_b = collect(body_R_b_to_w[:, :, b])
        anchor_a = collect(params.timoshenko_joints[j].anchor_a_b)
        anchor_b = collect(params.timoshenko_joints[j].anchor_b_b)
        pos_a = collect(body_pos_w[:, a])
        pos_b = collect(body_pos_w[:, b])
        com_a = collect(body_com_w[:, a])
        com_b = collect(body_com_w[:, b])

        x_a = pos_a .+ R_a * anchor_a
        x_b = pos_b .+ R_b * anchor_b
        e1, e2, e3, len = timoshenko_element_frame(x_a, x_b, R_a)
        element_frame = [e1[1] e2[1] e3[1];
                         e1[2] e2[2] e3[2];
                         e1[3] e2[3] e3[3]]

        L0 = params.timoshenko_joints[j].rest_length
        R_a_rel0 = collect(params.timoshenko_joints[j].R_a_rel0)
        R_b_rel0 = collect(params.timoshenko_joints[j].R_b_rel0)

        # Deformational rotation of each node, in the element frame.
        Da = (element_frame' * R_a) * R_a_rel0'
        Db = (element_frame' * R_b) * R_b_rel0'
        θ_a = [0.5 * (Da[3, 2] - Da[2, 3]),
               0.5 * (Da[1, 3] - Da[3, 1]),
               0.5 * (Da[2, 1] - Da[1, 2])]
        θ_b = [0.5 * (Db[3, 2] - Db[2, 3]),
               0.5 * (Db[1, 3] - Db[3, 1]),
               0.5 * (Db[2, 1] - Db[1, 2])]
        δ = len - L0

        kshear = params.timoshenko_joints[j].shear_coeff

        # Per-mode strain/curvature: axial, twist rate, bending curvatures, and the
        # two transverse shear angles. Nonlinear rigidities are evaluated at these.
        ε = δ / L0
        κt = (θ_b[1] - θ_a[1]) / L0
        κy = (θ_b[2] - θ_a[2]) / L0
        κz = (θ_b[3] - θ_a[3]) / L0
        γy = 0.5 * (θ_a[2] + θ_b[2])
        γz = 0.5 * (θ_a[3] + θ_b[3])

        EA_eff = timoshenko_rigidity(joint, params, :EA, ε)
        GJ_eff = timoshenko_rigidity(joint, params, :GJ, κt)
        EIy_eff = timoshenko_rigidity(joint, params, :EIy, κy)
        EIz_eff = timoshenko_rigidity(joint, params, :EIz, κz)
        GAy_eff = timoshenko_rigidity(joint, params, :GA, γy)
        GAz_eff = timoshenko_rigidity(joint, params, :GA, γz)

        Φy = 12 * EIy_eff / (kshear * GAy_eff * L0^2)
        Φz = 12 * EIz_eff / (kshear * GAz_eff * L0^2)
        by = EIy_eff / (L0 * (1 + Φy))
        bz = EIz_eff / (L0 * (1 + Φz))
        shy = 6 * EIy_eff / (L0^2 * (1 + Φy))
        shz = 6 * EIz_eff / (L0^2 * (1 + Φz))
        Mt = GJ_eff / L0
        f_axial = EA_eff * δ / L0

        # Restoring nodal forces/moments in the element frame. The `shy/shz` terms
        # are the transverse shear coupling absent from a lumped joint.
        F_a_local = [f_axial,
                     -shz * (θ_a[3] + θ_b[3]),
                      shy * (θ_a[2] + θ_b[2])]
        M_a_local = [-Mt * (θ_a[1] - θ_b[1]),
                     -by * ((4 + Φy) * θ_a[2] + (2 - Φy) * θ_b[2]),
                     -bz * ((4 + Φz) * θ_a[3] + (2 - Φz) * θ_b[3])]
        F_b_local = [-f_axial,
                      shz * (θ_a[3] + θ_b[3]),
                     -shy * (θ_a[2] + θ_b[2])]
        M_b_local = [-Mt * (θ_b[1] - θ_a[1]),
                     -by * ((2 - Φy) * θ_a[2] + (4 + Φy) * θ_b[2]),
                     -bz * ((2 - Φz) * θ_a[3] + (4 + Φz) * θ_b[3])]

        # Damping resists relative node velocity/spin (rigid-motion-free).
        ω_a_w = R_a * collect(body_ω_b[:, a])
        ω_b_w = R_b * collect(body_ω_b[:, b])
        vel_a = collect(body_com_vel[:, a]) .+ (ω_a_w × (x_a .- com_a))
        vel_b = collect(body_com_vel[:, b]) .+ (ω_b_w × (x_b .- com_b))
        Δv = vel_b .- vel_a
        Δω = ω_b_w .- ω_a_w
        c_t = params.timoshenko_joints[j].damping_trans
        c_r = params.timoshenko_joints[j].damping_rot

        eqs = [
            eqs
            timoshenko_force_a_w[:, j] ~ element_frame * F_a_local .+ c_t .* Δv
            timoshenko_force_b_w[:, j] ~ element_frame * F_b_local .- c_t .* Δv
            timoshenko_moment_a_w[:, j] ~ element_frame * M_a_local .+ c_r .* Δω
            timoshenko_moment_b_w[:, j] ~ element_frame * M_b_local .- c_r .* Δω
        ]

        force_on_a = collect(timoshenko_force_a_w[:, j])
        force_on_b = collect(timoshenko_force_b_w[:, j])
        moment_on_a = (x_a .- com_a) × force_on_a .+ collect(timoshenko_moment_a_w[:, j])
        moment_on_b = (x_b .- com_b) × force_on_b .+ collect(timoshenko_moment_b_w[:, j])

        body_force[:, a] .+= force_on_a
        body_force[:, b] .+= force_on_b
        body_moment[:, a] .+= moment_on_a
        body_moment[:, b] .+= moment_on_b
    end
    return eqs
end
