# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Aero coupling wiring (winch-style).
#
# Each wing's aerodynamics is a swappable subsystem built by
# `aero_component(wing.aero, …)` (see aero_modes/common.jl). This layer
# instantiates the component, validates its connector contract, drives the
# body-frame inputs, and reads the outputs back into the wing's aero variables.
# All built-in models (AeroNone / AeroDirect / AeroLinearized / AeroPlate) and
# any custom AbstractAeroModel go through the same wiring. Flat-plate aero uses
# the standard PARTICLE contract: it is the only mode that consumes the per-point
# `va`/`rho` inputs the wiring drives for every particle wing.

"""
    aero_eqs!(s, eqs; kwargs...)
        -> (eqs, aero_subsystems)

Instantiate and wire each wing's aero component. Returns the list of
component subsystems to attach to the parent `System`.
"""
function aero_eqs!(
    s, eqs, params;
    aero_force_b, aero_moment_b, twist_surface_aero_moment,
    twist_angle, twist_ω, va_wing_b, wing_pos, ω_b, R_b_to_w,
    pos, vel, va_point_b, height, aero_force_point_b=nothing
)
    (; twist_surfaces, wings, points) = s.sys_struct
    aero_subsystems = Any[]
    length(wings) == 0 && return eqs, aero_subsystems

    for wing in wings
        wing_idx = wing.idx
        subsys = aero_component(wing.aero, s.sys_struct, wing_idx;
                                name = Symbol("aero_$(wing_idx)"), params)
        push!(aero_subsystems, subsys)
        validate_aero_component(subsys, wing)

        if wing.dynamics_type == PARTICLE_DYNAMICS
            wing_points = [point for point in points
                           if point.type == WING && point.wing_idx == wing_idx]
            Rbw = R_b_to_w[:, :, wing_idx]
            aero_force_point = aero_force_point_b::AbstractArray
            for (k, point) in enumerate(wing_points)
                eqs = [eqs
                       collect(subsys.point_pos[:, k]) .~
                           collect(Rbw' * collect(pos[:, point.idx] -
                                                  wing_pos[:, wing_idx]))
                       collect(subsys.point_vel[:, k]) .~
                           collect(Rbw' * collect(vel[:, point.idx]))
                       collect(subsys.va[:, k]) .~
                           collect(va_point_b[:, point.idx])
                       subsys.rho[k] ~ calc_rho(s.am, height[point.idx])
                       collect(aero_force_point[:, point.idx]) .~
                           collect(subsys.point_force[:, k])]
            end
            eqs = [eqs
                   collect(aero_force_b[:, wing_idx]) .~
                       sum(collect(aero_force_point[:, point.idx])
                           for point in wing_points)
                   collect(aero_moment_b[:, wing_idx]) .~
                       sum(collect(cross(
                               collect(subsys.point_pos[:, k]),
                               collect(aero_force_point[:, wing_points[k].idx])))
                           for k in eachindex(wing_points))]
            continue
        end

        # RIGID_DYNAMICS
        num_twist_surfaces = length(wing.twist_surface_idxs)
        rho = calc_rho(s.am, wing_pos[3, wing_idx])
        eqs = [eqs
               collect(subsys.va) .~ collect(va_wing_b[:, wing_idx])
               subsys.rho ~ rho
               vec(collect(subsys.R_b_w)) .~ vec(R_b_to_w[:, :, wing_idx])
               collect(subsys.omega) .~ collect(ω_b[:, wing_idx])]
        if num_twist_surfaces > 0
            eqs = [eqs
                   collect(subsys.twist) .~
                       [twist_angle[twist_surfaces[gidx].idx]
                        for gidx in wing.twist_surface_idxs]
                   collect(subsys.twist_vel) .~
                       [twist_ω[twist_surfaces[gidx].idx]
                        for gidx in wing.twist_surface_idxs]]
        end

        eqs = [eqs
               collect(aero_force_b[:, wing_idx]) .~ collect(subsys.force)
               collect(aero_moment_b[:, wing_idx]) .~ collect(subsys.moment)]
        for (twist_surface_pos, gidx) in enumerate(wing.twist_surface_idxs)
            twist_surface = twist_surfaces[gidx]
            # FIXED surfaces without aero sections are zero-bound by
            # twist_surface_eqs!; everything else reads the component output.
            twist_surface.type == FIXED &&
                isempty(twist_surface.unrefined_section_idxs) && continue
            eqs = [eqs
                   twist_surface_aero_moment[twist_surface.idx] ~
                       subsys.twist_moment[twist_surface_pos]]
        end

    end
    return eqs, aero_subsystems
end
