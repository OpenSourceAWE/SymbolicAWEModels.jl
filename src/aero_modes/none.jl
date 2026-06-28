# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# AeroNone: zero aerodynamic force; inherits AbstractAeroModel defaults (common.jl).

"""
    AeroNone()

No aerodynamic forces (returns zeros). A body with `aero = AeroNone()` is a plain
rigid/kinematic body. Needs no VSM geometry and carries no state.
"""
struct AeroNone <: AbstractAeroModel end

"""A wing is a `Body` carrying aerodynamics; `Body{AeroNone}` is a plain body."""
is_wing(::Body{AeroNone}) = false
is_wing(::Body) = true

is_builtin_aero(::AeroNone) = true
aero_mode_tag(::AeroNone) = "none"
stores_point_force(::AeroNone) = false

"""
    aero_component(::AeroNone, wing::ParticleWing, sys_struct; name)

Zero per-point forces (particle connector contract). AeroNone supports both wing
dynamics via one method each.
"""
function aero_component(::AeroNone, wing::ParticleWing, sys_struct;
                        name, params=nothing)
    num_points = length(wing_points(sys_struct, wing))
    connectors = particle_aero_connectors(num_points)
    eqs = vec(collect(connectors.point_force)) .~ 0
    return System(eqs, t, particle_unknowns(connectors), []; name)
end

"""
    aero_component(::AeroNone, wing::RigidWing, sys_struct; name)

Zero lumped force/moment (rigid connector contract).
"""
function aero_component(::AeroNone, wing::RigidWing, sys_struct;
                        name, params=nothing)
    num_twist_surfaces = length(wing.twist_surface_idxs)
    connectors = rigid_aero_connectors(num_twist_surfaces)
    eqs = [collect(connectors.force) .~ 0
           collect(connectors.moment) .~ 0]
    num_twist_surfaces > 0 && (eqs = [eqs; collect(connectors.twist_moment) .~ 0])
    return System(eqs, t, rigid_unknowns(connectors), []; name)
end
