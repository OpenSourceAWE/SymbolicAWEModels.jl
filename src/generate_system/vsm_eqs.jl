# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Aero coupling wiring (winch-style).
#
# Each non-plate wing's aerodynamics is a swappable subsystem built by
# `wing.aero_model` (see aero_components.jl). This layer instantiates
# the component, validates its connector contract, drives the body-frame
# inputs, and reads the outputs back into the wing's aero variables.
# All built-in modes (AERO_NONE / AERO_DIRECT / AERO_LINEARIZED) and any
# AERO_CUSTOM model go through the same wiring.

"""
    vsm_eqs!(s, eqs, guesses, psys; kwargs...)
        -> (eqs, guesses, aero_subsystems)

Instantiate and wire each wing's aero component. Returns the list of
component subsystems to attach to the parent `System`.
"""
function vsm_eqs!(
    s, eqs, guesses, psys;
    aero_force_b, aero_moment_b, group_aero_moment,
    twist_angle, twist_ω, va_wing_b, wing_pos, ω_b, R_b_to_w,
    pos, vel, aero_force_point_b=nothing
)
    (; groups, wings, points) = s.sys_struct
    aero_subsystems = Any[]
    length(wings) == 0 && return eqs, guesses, aero_subsystems

    for wing in wings
        wing isa PlateWing && continue   # handled by plate_eqs!

        w = wing.idx
        subsys = wing.aero_model(s.sys_struct, w;
                                 name = Symbol("aero_$(w)"))
        validate_aero_component(subsys, wing)
        push!(aero_subsystems, subsys)

        if wing.dynamics_type == PARTICLE_DYNAMICS
            wing_points = [p for p in points
                           if p.type == WING && p.wing_idx == w]
            Rbw = R_b_to_w[:, :, w]
            afpb = aero_force_point_b::AbstractArray
            for (k, point) in enumerate(wing_points)
                eqs = [eqs
                       collect(subsys.point_pos[:, k]) .~
                           collect(Rbw' * collect(pos[:, point.idx] -
                                                  wing_pos[:, w]))
                       collect(subsys.point_vel[:, k]) .~
                           collect(Rbw' * collect(vel[:, point.idx]))
                       collect(afpb[:, point.idx]) .~
                           collect(subsys.point_force[:, k])]
            end
            eqs = [eqs
                   collect(aero_force_b[:, w]) .~
                       sum(collect(afpb[:, p.idx])
                           for p in wing_points)
                   collect(aero_moment_b[:, w]) .~ 0]
            continue
        end

        # RIGID_DYNAMICS
        ng = length(wing.group_idxs)
        rho = calc_rho(s.am, wing_pos[3, w])
        eqs = [eqs
               collect(subsys.va) .~ collect(va_wing_b[:, w])
               subsys.rho ~ rho
               vec(collect(subsys.R_b_w)) .~ vec(R_b_to_w[:, :, w])
               collect(subsys.omega) .~ collect(ω_b[:, w])]
        if ng > 0
            eqs = [eqs
                   collect(subsys.twist) .~
                       [twist_angle[groups[gidx].idx]
                        for gidx in wing.group_idxs]
                   collect(subsys.twist_vel) .~
                       [twist_ω[groups[gidx].idx]
                        for gidx in wing.group_idxs]]
        end

        eqs = [eqs
               collect(aero_force_b[:, w]) .~ collect(subsys.force)
               collect(aero_moment_b[:, w]) .~ collect(subsys.moment)]
        for (gi, gidx) in enumerate(wing.group_idxs)
            isempty(groups[gidx].unrefined_section_idxs) && continue
            eqs = [eqs
                   group_aero_moment[groups[gidx].idx] ~
                       subsys.twist_moment[gi]]
        end

        if s.set.quasi_static && hasproperty(subsys, :aero_input)
            ny = length(wing.aero_y)
            guesses = [guesses
                       [subsys.aero_input[iy] => get_aero_y(psys, w, iy)
                        for iy in 1:ny]]
        end
    end
    return eqs, guesses, aero_subsystems
end
