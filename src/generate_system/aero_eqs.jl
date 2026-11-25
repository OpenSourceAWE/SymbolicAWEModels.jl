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
    aero_force_b, aero_moment_b, station_aero_moment,
    twist_angle, twist_ω, va_wing_b, wing_pos, ω_b, R_b_to_w,
    pos, vel, va_point_b, height, aero_force_point_b=nothing,
    station_delta=nothing
)
    aero_subsystems = Any[]
    for (wing_idx, mode) in enumerate(s.aero)
        body_idx = wing_body(s, wing_idx)
        subsys = aero_component(mode, wing_dynamics(s, wing_idx), s, wing_idx;
                                name = Symbol("aero_$(wing_idx)"), params)
        push!(aero_subsystems, subsys)
        validate_aero_component(subsys, s, wing_idx)

        if body_dynamics_type(s, body_idx) == PARTICLE_DYNAMICS
            nodes = wing_points(s, wing_idx)
            aero_force_point = aero_force_point_b::AbstractArray
            wiring = particle_wing_aero_wiring(s, subsys;
                orientation = R_b_to_w[:, :, body_idx],
                origin = collect(wing_pos[:, body_idx]),
                positions = [pos[:, i] for i in nodes],
                velocities = [vel[:, i] for i in nodes],
                apparent_winds = [va_point_b[:, i] for i in nodes],
                heights = [pos[3, i] for i in nodes])
            eqs = [eqs
                   wiring.eqs
                   flap_delta_eqs(mode, subsys, surface -> station_delta[surface])
                   collect(aero_force_b[:, wing_idx]) .~ wiring.force_b
                   collect(aero_moment_b[:, wing_idx]) .~ wiring.moment_b]
            for (k, i) in enumerate(nodes)
                eqs = [eqs
                       collect(aero_force_point[:, i]) .~ collect(subsys.point_force[:, k])]
            end
            continue
        end

        # RIGID_DYNAMICS
        stations = wing_station_idxs(s, wing_idx)
        wiring = rigid_wing_aero_wiring(s, subsys, wing_idx;
            apparent_wind_b = va_wing_b[:, wing_idx],
            height = wing_pos[3, body_idx],
            frame = vec(R_b_to_w[:, :, body_idx]),
            omega_b = ω_b[:, body_idx],
            twist_angles = [twist_angle[station_idx] for station_idx in stations],
            twist_rates = [twist_ω[station_idx] for station_idx in stations])
        eqs = [eqs
               wiring.eqs
               collect(aero_force_b[:, wing_idx]) .~ wiring.force_b
               collect(aero_moment_b[:, wing_idx]) .~ wiring.moment_b]
        for (station_pos, station_idx) in enumerate(stations)
            station_aero_driven(s, station_idx) || continue
            eqs = [eqs
                   station_aero_moment[station_idx] ~ wiring.twist_moments[station_pos]]
        end
    end
    return eqs, aero_subsystems
end
