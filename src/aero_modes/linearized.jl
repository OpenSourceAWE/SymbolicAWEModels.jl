# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# AeroLinearized: first-order Taylor of VSM coefficients; numerics in common.jl.

"""
    AeroLinearized()

First-order Taylor expansion using the Jacobian from VSM linearization
(`RIGID_DYNAMICS` only). Carries a [`VSMEngine`](@ref); the no-arg form is the
engine-less marker filled in during wing construction.
"""
mutable struct AeroLinearized{E} <: AbstractVSMAero
    engine::Union{Nothing, E}
    AeroLinearized{E}(engine) where {E} = new{E}(engine)
end
AeroLinearized() = AeroLinearized{VSMEngine}(nothing)
AeroLinearized(engine::VSMEngine) = AeroLinearized{typeof(engine)}(engine)
attach_engine!(::AeroLinearized, engine::VSMEngine) = AeroLinearized(engine)

is_builtin_aero(::AeroLinearized) = true
aero_mode_tag(::AeroLinearized) = "lin"

function aero_component(::AeroLinearized, wing::RigidWing, sys_struct; name, params)
    wing_idx = wing.idx
    # Aero coefficients as flat params (frozen between VSM refreshes, synced per refresh).
    aero_y_p = params.wings[wing_idx].aero_y
    aero_x_p = params.wings[wing_idx].aero_x
    aero_jac_p = params.wings[wing_idx].aero_jac
    drag_frac_p = params.wings[wing_idx].drag_frac

    twist_surfaces = sys_struct.twist_surfaces
    num_twist_surfaces = length(wing.twist_surface_idxs)
    num_aero_inputs = length(wing.aero_y)
    area = wing.vsm_aero.projected_area
    c_ref = wing.vsm_aero.c_ref

    connectors = rigid_aero_connectors(num_twist_surfaces)
    @variables aero_input(t)[1:num_aero_inputs]

    apparent_wind = collect(connectors.va)
    omega = collect(connectors.omega)
    drag_dir = collect(apparent_wind ./ smooth_norm(apparent_wind))
    alpha = atan(drag_dir[3], drag_dir[1])
    beta = atan(drag_dir[2], smooth_norm((drag_dir[1], drag_dir[3])))

    twist_inputs = num_twist_surfaces > 0 ? collect(connectors.twist) : Num[]
    input_rhs = [alpha; beta; omega[1]; omega[2]; omega[3]; twist_inputs]

    delta(input_idx) = aero_input[input_idx] - aero_y_p[input_idx]
    coeff(output_idx) = aero_x_p[output_idx] +
        sum(aero_jac_p[output_idx, input_idx] * delta(input_idx)
            for input_idx in 1:num_aero_inputs)

    q_inf = 0.5 * connectors.rho * (apparent_wind ⋅ apparent_wind)
    qA = q_inf * area
    CL = coeff(1)
    CD = coeff(2)
    CS = coeff(3)

    crossed = collect(drag_dir × [0.0, 1.0, 0.0])
    lift_dir = collect(crossed ./ smooth_norm(crossed))
    side_dir = collect(lift_dir × drag_dir)
    drag_frac = drag_frac_p

    force_rhs = collect(qA * (CL * lift_dir +
        CD * drag_frac * drag_dir + CS * side_dir))
    moment_rhs = [qA * c_ref * coeff(3 + i) for i in 1:3]

    eqs = [collect(aero_input) .~ input_rhs
           collect(connectors.force) .~ force_rhs
           collect(connectors.moment) .~ moment_rhs]
    for twist_surface_pos in 1:num_twist_surfaces
        isempty(twist_surfaces[wing.twist_surface_idxs[twist_surface_pos]].unrefined_section_idxs) ?
            (eqs = [eqs; connectors.twist_moment[twist_surface_pos] ~ 0]) :
            (eqs = [eqs; connectors.twist_moment[twist_surface_pos] ~
                qA * c_ref * coeff(6 + twist_surface_pos)])
    end

    vars = rigid_unknowns(connectors)
    push!(vars, aero_input)
    return System(eqs, t, vars,
        [aero_y_p, aero_x_p, aero_jac_p, drag_frac_p]; name)
end

"""
    refresh_rigid_aero!(::AeroLinearized, wing, am, twist_surfaces; vsm_min_wind=0.5)

Linearized rigid-wing refresh. Computes the baseline wind-axis coefficients at the
operating point ([`rigid_aero_baseline!`](@ref)), then the `ForwardDiff` Jacobian
`d(coeffs)/d(inputs)` and stores it in `wing.aero_jac`. The compiled RHS uses that
Jacobian to reconstruct forces via a first-order Taylor expansion about the
operating point.
"""
function refresh_rigid_aero!(::AeroLinearized, wing, am, twist_surfaces;
                             vsm_min_wind=0.5)
    ctx = rigid_aero_baseline!(wing, twist_surfaces; vsm_min_wind)
    gamma0 = copy(wing.vsm_solver.sol.gamma_distribution)
    f_dual = y -> vsm_aero_coeffs(wing, y, ctx.va_mag, ctx.n_unrefined,
        ctx.n_twist_surfaces, ctx.twist_surface_idxs, twist_surfaces,
        ctx.moment_frac, ctx.shadow_ref; gamma_init=gamma0)
    ForwardDiff.jacobian!(wing.aero_jac, f_dual, ctx.y0)
    return nothing
end

"""
    refresh_particle_aero!(::AeroLinearized, wing, points, va_point_b_vals; vsm_min_wind=0.5)

Unsupported: `AeroLinearized` is not implemented for `PARTICLE_DYNAMICS` wings and
errors. Use `AeroDirect` (per-point nonlinear VSM) for particle wings.
"""
refresh_particle_aero!(::AeroLinearized, wing, points, va_point_b_vals;
                       vsm_min_wind=0.5) = error(
    "PARTICLE_DYNAMICS + AeroLinearized not yet implemented")
