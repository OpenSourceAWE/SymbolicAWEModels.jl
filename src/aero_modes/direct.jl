# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# AeroDirect: frozen piecewise-constant forces from the nonlinear VSM solve.

"""
    AeroDirect()

Stored forces from the nonlinear VSM solve, piecewise-constant between updates.
Carries a [`VSMEngine`](@ref); `AeroDirect()` selects the mode for a wing
([`AeroSpec`](@ref)).
"""
mutable struct AeroDirect{E} <: AbstractVSMAero
    engine::E
end
AeroDirect() = AeroSpec{AeroDirect}()

is_builtin_aero(::AeroDirect) = true
aero_mode_tag(::Type{<:AeroDirect}) = "dir"

"""
    aero_component(::AeroDirect, ::ParticleDynamics, sam, wing_idx; name, params)

Frozen per-point VSM forces of wing `wing_idx` of `sam` (synced each refresh, so emitted as
flat params) on the particle connector contract.
"""
aero_component(::AeroDirect, ::ParticleDynamics, sam, wing_idx; name, params=nothing) =
    frozen_point_force_component(sam, wing_idx; name, params)

"""
    aero_component(mode::AeroDirect, ::RigidDynamics, sam, wing_idx; name, params)

Frozen lumped VSM force/moment of wing `wing_idx` of `sam` (synced each refresh, so emitted
as flat params) on the rigid connector contract; a station with sections takes its held
twist moment.
"""
function aero_component(mode::AeroDirect, ::RigidDynamics, sam, wing_idx;
                        name, params=nothing)
    station_idxs = wing_station_idxs(sam, wing_idx)
    connectors = rigid_aero_connectors(length(station_idxs))
    force_p = params.wings[wing_idx].aero_force_b
    moment_p = params.wings[wing_idx].aero_moment_b
    flat_ps = Any[force_p, moment_p]
    eqs = [collect(connectors.force) .~ collect(force_p)
           collect(connectors.moment) .~ collect(moment_p)]
    for (k, station_idx) in enumerate(station_idxs)
        if isempty(mode.engine.station_sections[k])
            push!(eqs, connectors.twist_moment[k] ~ 0)
        else
            moment_ts_p = params.spanwise_stations[station_idx].held_aero_moment
            push!(flat_ps, moment_ts_p)
            push!(eqs, connectors.twist_moment[k] ~ moment_ts_p)
        end
    end
    return System(eqs, t, rigid_unknowns(connectors), flat_ps; name)
end

"""
    refresh_rigid_aero!(mode::AeroDirect, sam, wing_idx; vsm_min_wind=0.5)

Direct rigid-wing refresh of wing `wing_idx` of `sam`. Computes the baseline coefficients
and applies the resulting frozen body-frame force/moment ([`apply_direct_forces!`](@ref)),
which the RHS holds constant until the next refresh. Below `vsm_min_wind` the
coefficients, Jacobian, force, moment, and per-station moments are zeroed.
"""
function refresh_rigid_aero!(mode::AeroDirect, sam, wing_idx; vsm_min_wind=0.5)
    if norm(wing_apparent_wind(sam, wing_idx)) < vsm_min_wind
        fill!(mode.aero_x, 0.0)
        fill!(mode.aero_jac, 0.0)
        set_run_value!(sam, :wings, wing_idx, :aero_force_b, zero(SVec3))
        set_run_value!(sam, :wings, wing_idx, :aero_moment_b, zero(SVec3))
    else
        rigid_aero_baseline!(mode, sam, wing_idx; vsm_min_wind)
        apply_direct_forces!(mode, sam, wing_idx, mode.aero_x)
    end
    hold_station_moments!(mode, sam, wing_idx)
    return nothing
end

"""Hold each station of wing `wing_idx` of `sam` that has sections at the twist moment [N·m]
the baseline coefficients `aero_x` of its aero `mode` give, until the next refresh."""
function hold_station_moments!(mode, sam, wing_idx)
    for (k, station_idx) in enumerate(wing_station_idxs(sam, wing_idx))
        isempty(mode.engine.station_sections[k]) && continue
        set_run_value!(sam, :spanwise_stations, station_idx, :held_aero_moment,
                       mode.aero_x[6 + k])
    end
    return nothing
end

"""Set the aerodynamic force of every node of wing `wing_idx` of `sam` to zero."""
function clear_point_aero_forces!(sam, wing_idx)
    for i in wing_points(sam, wing_idx)
        set_run_value!(sam, Val(:points), i, Val(:aero_force), zero(SVec3))
    end
    return nothing
end

"""
    refresh_particle_aero!(mode::AeroDirect, sam, wing_idx, va_point_b; vsm_min_wind=0.5)

Direct particle-wing refresh of wing `wing_idx` of `sam`. Runs the full nonlinear VSM solve
using each section's apparent wind (averaged from its LE/TE point columns of
`va_point_b`), then distributes the resulting panel forces onto the wing's structural
points ([`distribute_panel_forces_to_points!`](@ref)). Below `vsm_min_wind` the point
forces are zeroed. A failed solve throws ([`safe_vsm_solve!`](@ref)) before any point
force is written.
"""
function refresh_particle_aero!(mode::AeroDirect, sam, wing_idx, va_point_b;
vsm_min_wind=0.5)
    if norm(wing_apparent_wind(sam, wing_idx)) < vsm_min_wind
        clear_point_aero_forces!(sam, wing_idx)
        return nothing
    end
    update_vsm_wing_from_structure!(mode, sam, wing_idx)
    set_particle_panel_va!(mode, sam, wing_idx, va_point_b)
    solve_circulation!(sam, wing_idx, mode.vsm_solver, mode.vsm_aero)
    distribute_panel_forces_to_points!(mode, sam, wing_idx)
    for i in wing_points(sam, wing_idx)
        any(!isfinite, run_value(sam, Val(:points), i, Val(:aero_force_b))) &&
            throw(AssertionError("PARTICLE_DYNAMICS: non-finite point force on wing " *
                                 "$wing_idx point $i"))
    end
    return nothing
end

"""Set the frozen body-frame force and moment of wing `wing_idx` of `sam` from the wind-axis
coefficients `x0` of its aero `mode`."""
function apply_direct_forces!(mode, sam, wing_idx, x0)
    body_idx = wing_body(sam, wing_idx)
    va_b = run_value(sam, Val(:bodies), body_idx, Val(:va_b))
    if any(!isfinite, x0) || any(!isfinite, va_b)
        throw(AssertionError("AeroDirect: non-finite input on wing $wing_idx"))
    end
    height = run_value(sam, Val(:bodies), body_idx, Val(:pos_w))[3]
    q_inf = 0.5 * air_density(sam.am, height) * dot(va_b, va_b)
    area = mode.vsm_aero.projected_area
    c_ref = mode.vsm_aero.c_ref

    CL, CD, CS = x0[1], x0[2], x0[3]
    span = SVector(0.0, 1.0, 0.0)
    drag_dir = va_b / norm(va_b)
    lift_dir = smooth_normalize(cross(drag_dir, span))
    side_dir = cross(lift_dir, drag_dir)
    drag_frac = run_value(sam, Val(:bodies), body_idx, Val(:drag_frac))

    set_run_value!(sam, :wings, wing_idx, :aero_force_b,
                   q_inf * area * (CL .* lift_dir .+ CD * drag_frac .* drag_dir .+
                                   CS .* side_dir))
    set_run_value!(sam, :wings, wing_idx, :aero_moment_b, q_inf * area * c_ref .* x0[4:6])
    return nothing
end
