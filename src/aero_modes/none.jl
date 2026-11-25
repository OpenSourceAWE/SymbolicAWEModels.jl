# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# AeroNone: zero aerodynamic force; inherits AbstractAeroModel defaults (common.jl).

"""
    AeroNone()

No aerodynamic forces (returns zeros). A body with `aero = AeroNone()` is a plain
rigid/kinematic body. Needs no VSM geometry and carries no state.
"""
struct AeroNone <: AbstractAeroModel end

"""Whether a wing flown with aero `mode` has aerodynamics: one flown with `AeroNone` has
not."""
has_aero(::AbstractAeroModel) = true
has_aero(::AeroNone) = false

is_builtin_aero(::AeroNone) = true
aero_mode_tag(::Type{AeroNone}) = "none"

"""
    aero_component(::AeroNone, ::ParticleDynamics, sam, wing_idx; name)

Zero per-point forces on the nodes of wing `wing_idx` of `sam` (particle connector
contract). AeroNone supports both wing dynamics via one method each.
"""
function aero_component(::AeroNone, ::ParticleDynamics, sam, wing_idx; name, params=nothing)
    connectors = particle_aero_connectors(length(wing_points(sam, wing_idx)))
    eqs = vec(collect(connectors.point_force)) .~ 0
    return System(eqs, t, particle_unknowns(connectors), []; name)
end

"""
    aero_component(::AeroNone, ::RigidDynamics, sam, wing_idx; name)

Zero lumped force/moment of wing `wing_idx` of `sam` (rigid connector contract).
"""
function aero_component(::AeroNone, ::RigidDynamics, sam, wing_idx; name, params=nothing)
    num_stations = length(wing_station_idxs(sam, wing_idx))
    connectors = rigid_aero_connectors(num_stations)
    eqs = [collect(connectors.force) .~ 0
           collect(connectors.moment) .~ 0]
    num_stations > 0 && (eqs = [eqs; collect(connectors.twist_moment) .~ 0])
    return System(eqs, t, rigid_unknowns(connectors), []; name)
end
