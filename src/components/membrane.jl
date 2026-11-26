# Copyright (c) 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# The membrane triangle: a constant-strain triangle of Saint-Venant–Kirchhoff fabric in
# plane stress whose energy is relaxed for wrinkling (Pipkin). The energy is the one
# definition; both backends take its corner forces as its symbolic gradient.

"""Floor under the squared radius of the principal strains' Mohr circle, keeping their
derivative finite under equal biaxial strain."""
const MEMBRANE_STRAIN_FLOOR = 1e-12

"""
    membrane_rest_metric(rest_lengths) -> (G11, G12, G22)

The metric of a triangle at rest in its edge basis, corner 1 to 2 and corner 1 to 3, from
the rest lengths of its edges 1→2, 2→3 and 3→1 [m²].
"""
function membrane_rest_metric(rest_lengths)
    l12, l23, l31 = rest_lengths[1], rest_lengths[2], rest_lengths[3]
    return l12^2, (l12^2 + l31^2 - l23^2) / 2, l31^2
end

"""
    membrane_rest_area(rest_lengths)

Area of a triangle at rest from the rest lengths of its edges [m²].
"""
function membrane_rest_area(rest_lengths)
    G11, G12, G22 = membrane_rest_metric(rest_lengths)
    return sqrt(max(G11 * G22 - G12^2, 0.0)) / 2
end

"""
    membrane_strain_energy(x1, x2, x3, rest_lengths, youngs_modulus, poisson_ratio,
                           thickness, slack_frac)

Strain energy [J] of a membrane triangle with corners at `x1`, `x2`, `x3` [m]. Taut,
where the principal Green strains satisfy `e2 + ν e1 ≥ 0`, it is Saint-Venant–Kirchhoff in
plane stress; wrinkled, where only `e1 > 0`, the uniaxial `Y e1² / 2`; slack, zero. A
`slack_frac` of the taut energy is kept everywhere.
"""
function membrane_strain_energy(x1, x2, x3, rest_lengths, youngs_modulus, poisson_ratio,
                                thickness, slack_frac)
    G11, G12, G22 = membrane_rest_metric(rest_lengths)
    det_G = G11 * G22 - G12^2
    edge_a = x2 .- x1
    edge_b = x3 .- x1
    g11, g12, g22 = edge_a ⋅ edge_a, edge_a ⋅ edge_b, edge_b ⋅ edge_b
    trace_C = (G22 * g11 - 2 * G12 * g12 + G11 * g22) / det_G
    det_C = (g11 * g22 - g12^2) / det_G
    trace_E = (trace_C - 2) / 2
    det_E = (det_C - trace_C + 1) / 4
    radius = sqrt((trace_E / 2)^2 - det_E + MEMBRANE_STRAIN_FLOOR)
    e1 = trace_E / 2 + radius
    e2 = trace_E / 2 - radius
    ν = poisson_ratio
    taut = youngs_modulus / (2 * (1 - ν^2)) * (trace_E^2 - 2 * (1 - ν) * det_E)
    relaxed = ifelse(e2 + ν * e1 >= 0, taut,
                     ifelse(e1 > 0, youngs_modulus * e1^2 / 2, zero(taut)))
    return membrane_rest_area(rest_lengths) * thickness *
           ((1 - slack_frac) * relaxed + slack_frac * taut)
end

"""
    membrane_dissipation(x1, x2, x3, v1, v2, v3, rest_lengths, thickness, damping)

Rayleigh dissipation function [W/2] of a membrane triangle with corners at `x1`, `x2`,
`x3` [m] moving at `v1`, `v2`, `v3` [m/s]: `damping` [Pa·s] times half the squared Green
strain rate, over the fabric's volume.
"""
function membrane_dissipation(x1, x2, x3, v1, v2, v3, rest_lengths, thickness, damping)
    G11, G12, G22 = membrane_rest_metric(rest_lengths)
    det_G = G11 * G22 - G12^2
    edge_a, edge_b = x2 .- x1, x3 .- x1
    rate_a, rate_b = v2 .- v1, v3 .- v1
    rate_11 = 2 * (edge_a ⋅ rate_a)
    rate_12 = edge_a ⋅ rate_b + rate_a ⋅ edge_b
    rate_22 = 2 * (edge_b ⋅ rate_b)
    M11 = (G22 * rate_11 - G12 * rate_12) / det_G
    M12 = (G22 * rate_12 - G12 * rate_22) / det_G
    M21 = (G11 * rate_12 - G12 * rate_11) / det_G
    M22 = (G11 * rate_22 - G12 * rate_12) / det_G
    strain_rate_sq = (M11^2 + 2 * M12 * M21 + M22^2) / 4
    return membrane_rest_area(rest_lengths) * thickness * damping * strain_rate_sq / 2
end

"""
    membrane_corner_forces(params, idx, positions, velocities) -> Matrix

Force on each corner of membrane `idx` [N], one column per corner: minus the gradient of
its [`membrane_strain_energy`](@ref) in the corner positions and of its
[`membrane_dissipation`](@ref) in the corner velocities, both differentiated symbolically.
`positions` and `velocities` hold one 3-vector per corner; the fabric's properties are
read from `params.membranes[idx]`.
"""
function membrane_corner_forces(params, idx, positions, velocities)
    fabric = params.membranes[idx]
    rest_lengths = [fabric.rest_lengths[i] for i in 1:3]
    @variables corner_pos[1:3, 1:3] corner_vel[1:3, 1:3]
    x, v = collect(corner_pos), collect(corner_vel)
    energy = membrane_strain_energy(x[:, 1], x[:, 2], x[:, 3], rest_lengths,
        fabric.youngs_modulus, fabric.poisson_ratio, fabric.thickness, fabric.slack_frac)
    dissipation = membrane_dissipation(x[:, 1], x[:, 2], x[:, 3], v[:, 1], v[:, 2],
        v[:, 3], rest_lengths, fabric.thickness, fabric.damping)
    force = .-(Symbolics.gradient(energy, vec(x)) .+
               Symbolics.gradient(dissipation, vec(v)))
    corners = Dict([vec(x) .=> reduce(vcat, collect.(positions));
                    vec(v) .=> reduce(vcat, collect.(velocities))])
    return reshape(substitute.(force, (corners,)), 3, 3)
end

"""
    membrane_variables()

The variables a membrane kernel declares, as a named tuple: each corner's `pos`/`vel`
inputs and the `force` output on it.
"""
function membrane_variables()
    vars = @variables begin
        p1_pos(t)[1:3], [input = true]
        p1_vel(t)[1:3], [input = true]
        p2_pos(t)[1:3], [input = true]
        p2_vel(t)[1:3], [input = true]
        p3_pos(t)[1:3], [input = true]
        p3_vel(t)[1:3], [input = true]
        p1_force(t)[1:3], [output = true]
        p2_force(t)[1:3], [output = true]
        p3_force(t)[1:3], [output = true]
    end
    return (; pos = vars[[1, 3, 5]], vel = vars[[2, 4, 6]], force = vars[7:9],
            all = vars)
end

"""
    MembraneTriangle(params, idx; name)

The kernel component of membrane `idx`: its corner forces
([`membrane_corner_forces`](@ref)) from its corners' positions and velocities.
"""
function MembraneTriangle(params, idx; name)
    io = membrane_variables()
    forces = membrane_corner_forces(params, idx, io.pos, io.vel)
    eqs = reduce(vcat, [collect(io.force[k]) .~ forces[:, k] for k in 1:3])
    return System(eqs, t, io.all, param_unknowns(params); name)
end
