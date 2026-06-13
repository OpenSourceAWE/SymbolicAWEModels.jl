# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# AeroDirect: frozen forces from the nonlinear VSM solve, held piecewise-constant
# between refreshes. Works for both RIGID_DYNAMICS and PARTICLE_DYNAMICS. Shared
# VSM numerics (rigid_aero_baseline!, the solve helpers) live in common.jl.

"""
    AeroDirect()

Stored forces from the nonlinear VSM solve, piecewise-constant between updates.
Carries a [`VSMEngine`](@ref); the no-arg form is the engine-less marker filled
in during wing construction.
"""
mutable struct AeroDirect <: AbstractVSMAero
    engine::Union{Nothing, VSMEngine}
end
AeroDirect() = AeroDirect(nothing)

is_builtin_aero(::AeroDirect) = true
aero_mode_tag(::AeroDirect) = "dir"
provides_aero_override(::AeroDirect) = true

function aero_component(::AeroDirect, sys_struct, wing_idx; name)
    psys = system_struct_param(sys_struct)
    wing = sys_struct.wings[wing_idx]

    if wing.dynamics_type == PARTICLE_DYNAMICS
        points = wing_points(sys_struct, wing)
        num_points = length(points)
        connectors = particle_aero_connectors(num_points)
        eqs = Equation[]
        for (point_num, point) in enumerate(points)
            eqs = [eqs
                   collect(connectors.point_force[:, point_num]) .~
                       [get_point_aero_force(psys, point.idx, i)
                        for i in 1:3]]
        end
        return System(eqs, t, particle_unknowns(connectors), [psys]; name)
    elseif wing.dynamics_type == RIGID_DYNAMICS
        twist_surfaces = sys_struct.twist_surfaces
        num_twist_surfaces = length(wing.twist_surface_idxs)
        connectors = rigid_aero_connectors(num_twist_surfaces)
        eqs = [collect(connectors.force) .~
                   [get_aero_force_override(psys, wing.idx, i) for i in 1:3]
               collect(connectors.moment) .~
                   [get_aero_moment_override(psys, wing.idx, i) for i in 1:3]]
        for (twist_surface_pos, twist_surface_idx) in enumerate(wing.twist_surface_idxs)
            rhs = isempty(twist_surfaces[twist_surface_idx].unrefined_section_idxs) ? 0 :
                get_twist_surface_moment_override(psys, wing.idx, Int64(twist_surface_idx))
            eqs = [eqs; connectors.twist_moment[twist_surface_pos] ~ rhs]
        end
        return System(eqs, t, rigid_unknowns(connectors), [psys]; name)
    else
        error("Unknown dynamics_type $(wing.dynamics_type) for wing $wing_idx.")
    end
end

"""
    refresh_rigid_aero!(::AeroDirect, wing, am, twist_surfaces; vsm_min_wind=0.5)

Direct rigid-wing refresh. Computes the baseline coefficients and applies the
resulting frozen body-frame force/moment ([`apply_direct_forces!`](@ref)), which
the RHS holds constant until the next refresh. Below `vsm_min_wind` the
coefficients, Jacobian, force, moment, and per-twist-surface moments are zeroed.
"""
function refresh_rigid_aero!(::AeroDirect, wing, am, twist_surfaces;
                             vsm_min_wind=0.5)
    if norm(wing.va_b) < vsm_min_wind
        fill!(wing.aero_x, 0.0)
        fill!(wing.aero_jac, 0.0)
        fill!(wing.aero_force_b, 0.0)
        fill!(wing.aero_moment_b, 0.0)
        for gidx in wing.twist_surface_idxs
            twist_surfaces[gidx].aero_moment = 0.0
        end
        return nothing
    end
    rigid_aero_baseline!(wing, twist_surfaces; vsm_min_wind)
    apply_direct_forces!(wing, am, wing.aero_x)
    return nothing
end

"""
    refresh_particle_aero!(::AeroDirect, wing, points, va_point_b_vals; vsm_min_wind=0.5)

Direct particle-wing refresh. Runs the full nonlinear VSM solve using each
section's apparent wind (averaged from its LE/TE point velocities in
`va_point_b_vals`), then distributes the resulting panel forces onto the wing's
structural points ([`distribute_panel_forces_to_points!`](@ref)). Below
`vsm_min_wind` the point forces are zeroed.
"""
function refresh_particle_aero!(::AeroDirect, wing, points, va_point_b_vals;
                                vsm_min_wind=0.5)
    if norm(wing.va_b) < vsm_min_wind
        for point in points
            if point.type == WING && point.wing_idx == wing.idx
                fill!(point.aero_force_b, 0.0)
            end
        end
        return nothing
    end

    update_vsm_wing_from_structure!(wing, points)
    set_particle_panel_va!(wing, va_point_b_vals)

    if !safe_vsm_solve!(wing.vsm_solver, wing.vsm_aero)
        throw(AssertionError("PARTICLE_DYNAMICS VSM solve failed (non-converged or non-finite) on wing $(wing.idx)"))
    end
    distribute_panel_forces_to_points!(wing, points)
    for point in points
        if point.type == WING && point.wing_idx == wing.idx &&
                any(!isfinite, point.aero_force_b)
            throw(AssertionError("PARTICLE_DYNAMICS: non-finite point force on wing $(wing.idx) point $(point.idx)"))
        end
    end
    return nothing
end

"""Apply direct forces from wind-axis coefficients."""
function apply_direct_forces!(wing, am, x0)
    va_b = wing.va_b
    if any(!isfinite, x0) || any(!isfinite, va_b)
        throw(AssertionError("AeroDirect: non-finite input on wing $(wing.idx)"))
    end
    va_sq = dot(va_b, va_b)
    rho = calc_rho(am, wing.pos_w[3])
    q_inf = 0.5 * rho * va_sq
    area = wing.vsm_aero.projected_area
    c_ref = wing.vsm_aero.c_ref

    CL, CD, CS = x0[1], x0[2], x0[3]
    span = SVector(0.0, 1.0, 0.0)
    drag_dir = va_b / norm(va_b)
    lift_dir = smooth_normalize(cross(drag_dir, span))
    side_dir = cross(lift_dir, drag_dir)

    wing.aero_force_b .= q_inf * area * (
        CL .* lift_dir .+
        CD * wing.drag_frac .* drag_dir .+
        CS .* side_dir)
    wing.aero_moment_b .= q_inf * area * c_ref .*
        x0[4:6]
end
