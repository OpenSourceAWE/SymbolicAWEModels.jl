# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Lumped 6-DOF elastic tube equation generation.

"""
    tube_stiffness_term(tube_idx, params, kind, Δ)

Restoring force/moment for one DOF of elastic tube `tube_idx` of the run `params` reads,
as a flat parameter: a `Real` stiffness is a numeric scalar param (`k·Δ`); an
interpolation is a callable param applied as `k(Δ)`. `kind`: 1=axial, 2=shear,
3=torsion, 4=bending.
"""
function tube_stiffness_term(tube_idx, params, kind::Int, Δ)
    field = rigidity_fields(ElasticTubeModel())[kind]
    k = getproperty(params.tubes[tube_idx], field)
    return run_value(params.reg.sam, :tubes, tube_idx, field) isa Real ? k * Δ : k(Δ)
end

"""
    tube_rayleigh_term(tube_idx, params, kind, Δ, rate, beta)

Rayleigh stiffness-proportional damping for one tube DOF: the restoring map
evaluated at `Δ + beta*rate` minus at `Δ`. That is `beta·K_tangent·rate` to first
order and exact `beta·k·rate` for a `Real` stiffness, and it vanishes identically
when `rate` is zero, so rigid motion stays undamped whatever the stiffness law.
"""
function tube_rayleigh_term(tube_idx, params, kind::Int, Δ, rate, beta)
    return tube_stiffness_term(tube_idx, params, kind, Δ + beta * rate) -
           tube_stiffness_term(tube_idx, params, kind, Δ)
end

"""
    elastic_tube_eqs!(sam, eqs, params; kwargs...)

For each elastic tube ([`ElasticTubeModel`](@ref)) of the run `sam`, compute the
restoring wrench from the relative pose of the two anchors (in body A's frame) and
accumulate it — equal and opposite — into `body_force`/`body_moment` (the same
accumulators `body_eqs!` reads). The
relative rotation uses the small-angle vector extraction, exact for the small
per-tube rotations of a stiff chain.
"""
function elastic_tube_eqs!(
    sam, eqs, params;
    body_force, body_moment,
    body_com_w, body_pos_w, body_com_vel, body_ω_b, body_R_b_to_w,
)
    elastic_tubes = [tube_idx for (tube_idx, tube) in enumerate(sam.definition.tubes)
                     if tube_columns(tube) isa ElasticTubeModel]
    @variables begin
        joint_force_w(t)[1:3, eachindex(elastic_tubes)]
        joint_torque_w(t)[1:3, eachindex(elastic_tubes)]
    end

    for (j, tube_idx) in enumerate(elastic_tubes)
        a, b = sam.definition.tubes[tube_idx].bodies
        R_a = collect(body_R_b_to_w[:, :, a])
        R_b = collect(body_R_b_to_w[:, :, b])
        ex = elastic_tube_wrench(tube_idx, params;
            force_w = joint_force_w[:, j], torque_w = joint_torque_w[:, j],
            pos_a = collect(body_pos_w[:, a]), R_a, com_a = collect(body_com_w[:, a]),
            com_vel_a = collect(body_com_vel[:, a]),
            omega_a_w = R_a * collect(body_ω_b[:, a]),
            pos_b = collect(body_pos_w[:, b]), R_b, com_b = collect(body_com_w[:, b]),
            com_vel_b = collect(body_com_vel[:, b]),
            omega_b_w = R_b * collect(body_ω_b[:, b]))
        eqs = [eqs; ex.tear_eqs]
        body_force[:, a] .+= ex.force_on_a
        body_force[:, b] .+= ex.force_on_b
        body_moment[:, a] .+= ex.moment_on_a
        body_moment[:, b] .+= ex.moment_on_b
    end
    return eqs
end
