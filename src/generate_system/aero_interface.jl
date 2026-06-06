# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

function angle_of_attack_sideslip(va)
    drag_dir = va ./ norm(va)
    alpha = atan(drag_dir[3], drag_dir[1])
    beta = atan(drag_dir[2], hypot(drag_dir[1], drag_dir[3]))
    return alpha, beta, drag_dir
end

function wind_axis_basis(drag_dir)
    span_body = SVector(0.0, 1.0, 0.0)
    lift_dir = normalize(cross(drag_dir, span_body))
    side_dir = cross(lift_dir, drag_dir)
    return lift_dir, side_dir
end

function dynamic_pressure_area(psys, wing, va)
    rho = calc_rho(psys.am, wing.pos_w[3])
    0.5 * rho * dot(va, va) * wing.vsm_aero.projected_area
end

function operating_point_component(j, alpha, beta, ω, twists)
    j == 1 && return alpha
    j == 2 && return beta
    j <= 5 && return ω[j - 2]
    return twists[j - 5]
end

function linearized_coefficient(wing, i, alpha, beta, ω, twists)
    acc = wing.aero_x[i]
    @inbounds for j in 1:length(wing.aero_y)
        delta = operating_point_component(j, alpha, beta, ω, twists) -
                wing.aero_y[j]
        acc += wing.aero_jac[i, j] * delta
    end
    return acc
end

function aero_force_moment(::DiscreteAero, psys, wing_idx, va, ω, twists)
    wing = psys.wings[wing_idx]
    f, m = wing.aero_force_b, wing.aero_moment_b
    SVector(f[1], f[2], f[3], m[1], m[2], m[3])
end

function aero_force_moment(::LinearizedAero, psys, wing_idx, va, ω, twists)
    wing = psys.wings[wing_idx]
    alpha, beta, drag_dir = angle_of_attack_sideslip(va)
    lift_dir, side_dir = wind_axis_basis(drag_dir)
    CL = linearized_coefficient(wing, 1, alpha, beta, ω, twists)
    CD = linearized_coefficient(wing, 2, alpha, beta, ω, twists)
    CS = linearized_coefficient(wing, 3, alpha, beta, ω, twists)
    qA = dynamic_pressure_area(psys, wing, va)
    drag_frac = get_drag_frac(psys, wing.idx)
    force = qA * (CL * lift_dir + CD * drag_frac * drag_dir + CS * side_dir)
    qA_cref = qA * wing.vsm_aero.c_ref
    CM1 = linearized_coefficient(wing, 4, alpha, beta, ω, twists)
    CM2 = linearized_coefficient(wing, 5, alpha, beta, ω, twists)
    CM3 = linearized_coefficient(wing, 6, alpha, beta, ω, twists)
    SVector(force[1], force[2], force[3],
            qA_cref * CM1, qA_cref * CM2, qA_cref * CM3)
end

function aero_force_moment(::ContinuousAero, psys, wing_idx, va, ω, twists)
    error("ContinuousAero not yet implemented (needs the VSM frozen-Γ kernel)")
end

function aero_twist_moment(::DiscreteAero, psys, wing_idx, gi, va, ω, twists)
    wing = psys.wings[wing_idx]
    psys.groups[wing.group_idxs[gi]].aero_moment
end

function aero_twist_moment(::LinearizedAero, psys, wing_idx, gi, va, ω, twists)
    wing = psys.wings[wing_idx]
    alpha, beta, _ = angle_of_attack_sideslip(va)
    coeff = linearized_coefficient(wing, 6 + gi, alpha, beta, ω, twists)
    dynamic_pressure_area(psys, wing, va) * wing.vsm_aero.c_ref * coeff
end

function aero_twist_moment(::ContinuousAero, psys, wing_idx, gi, va, ω, twists)
    error("ContinuousAero not yet implemented (needs the VSM frozen-Γ kernel)")
end

function wing_force_moment(psys, wing_idx, va1, va2, va3, ω1, ω2, ω3, twists)
    wing = psys.wings[wing_idx]
    aero_force_moment(wing.aero_type, psys, wing_idx,
                      SVector(va1, va2, va3), SVector(ω1, ω2, ω3), twists)
end

@register_array_symbolic wing_force_moment(
    psys::SystemStructure, wing_idx::Int,
    va1::Real, va2::Real, va3::Real, ω1::Real, ω2::Real, ω3::Real,
    twists::AbstractVector
) begin
    size = (6,)
    eltype = eltype(twists)
    ndims = 1
end

function group_twist_moment(psys, wing_idx, gi, va1, va2, va3, ω1, ω2, ω3, twists)
    wing = psys.wings[wing_idx]
    aero_twist_moment(wing.aero_type, psys, wing_idx, gi,
                      SVector(va1, va2, va3), SVector(ω1, ω2, ω3), twists)
end

@register_symbolic group_twist_moment(
    psys::SystemStructure, wing_idx::Int, gi::Int,
    va1::Real, va2::Real, va3::Real, ω1::Real, ω2::Real, ω3::Real,
    twists::AbstractVector)
