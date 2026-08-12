# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Aero coupling wiring (winch-style): each wing's aero is a swappable subsystem.

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
    pos, vel, va_point_b, height, aero_force_point_b=nothing,
    twist_surface_delta=nothing
)
    (; twist_surfaces, wings) = s.sys_struct
    aero_subsystems = Any[]
    length(wings) == 0 && return eqs, aero_subsystems

    for wing in wings
        wing_idx = wing.idx
        subsys = aero_component(wing.aero, wing, s.sys_struct;
                                name = Symbol("aero_$(wing_idx)"), params)
        push!(aero_subsystems, subsys)
        validate_aero_component(subsys, wing)

        if wing.dynamics_type == PARTICLE_DYNAMICS
            nodes = wing_points(s.sys_struct, wing)
            aero_force_point = aero_force_point_b::AbstractArray
            wiring = particle_wing_aero_wiring(s, subsys;
                orientation = R_b_to_w[:, :, wing_idx],
                origin = collect(wing_pos[:, wing_idx]),
                positions = [pos[:, node.idx] for node in nodes],
                velocities = [vel[:, node.idx] for node in nodes],
                apparent_winds = [va_point_b[:, node.idx] for node in nodes],
                heights = [height[node.idx] for node in nodes])
            eqs = [eqs
                   wiring.eqs
                   flap_delta_eqs(wing, subsys,
                                  surface -> twist_surface_delta[surface])
                   collect(aero_force_b[:, wing_idx]) .~ wiring.force_b
                   collect(aero_moment_b[:, wing_idx]) .~ wiring.moment_b]
            for (k, node) in enumerate(nodes)
                eqs = [eqs
                       collect(aero_force_point[:, node.idx]) .~
                           collect(subsys.point_force[:, k])]
            end
            continue
        end

        # RIGID_DYNAMICS
        wiring = rigid_wing_aero_wiring(s, subsys, wing;
            apparent_wind_b = va_wing_b[:, wing_idx],
            height = wing_pos[3, wing_idx],
            frame = vec(R_b_to_w[:, :, wing_idx]),
            omega_b = ω_b[:, wing_idx],
            twist_angles = [twist_angle[twist_surfaces[gidx].idx]
                            for gidx in wing.twist_surface_idxs],
            twist_rates = [twist_ω[twist_surfaces[gidx].idx]
                           for gidx in wing.twist_surface_idxs])
        eqs = [eqs
               wiring.eqs
               collect(aero_force_b[:, wing_idx]) .~ wiring.force_b
               collect(aero_moment_b[:, wing_idx]) .~ wiring.moment_b]
        for (twist_surface_pos, gidx) in enumerate(wing.twist_surface_idxs)
            twist_surface = twist_surfaces[gidx]
            twist_surface_aero_driven(twist_surface) || continue
            eqs = [eqs
                   twist_surface_aero_moment[twist_surface.idx] ~
                       wiring.twist_moments[twist_surface_pos]]
        end

    end
    return eqs, aero_subsystems
end
