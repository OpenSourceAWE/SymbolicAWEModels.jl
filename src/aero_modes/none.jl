# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# AeroNone: zero aerodynamic force. VSM-backed (carries an engine so the geometry
# stays available) but produces no force. See common.jl for the interface.

"""
    AeroNone()

No aerodynamic forces (returns zeros). For debugging rigid body dynamics or a
wing with no aero coupling. A VSM-backed mode ([`AbstractVSMAero`](@ref)): it
carries a [`VSMEngine`](@ref) built from the wing's VSM geometry (so the geometry
stays available) but produces zero force.
"""
mutable struct AeroNone <: AbstractVSMAero
    engine::Union{Nothing, VSMEngine}
end
AeroNone() = AeroNone(nothing)

is_builtin_aero(::AeroNone) = true
aero_mode_tag(::AeroNone) = "none"
couples_to_sections(::AeroNone) = false
stores_point_force(::AeroNone) = false
calc_aoa(::AeroNone, wing) = SimFloat(NaN)

function aero_component(::AeroNone, sys_struct, wing_idx; name)
    psys = system_struct_param(sys_struct)
    wing = sys_struct.wings[wing_idx]

    if wing.dynamics_type == PARTICLE_DYNAMICS
        num_points = length(wing_points(sys_struct, wing))
        connectors = particle_aero_connectors(num_points)
        eqs = vec(collect(connectors.point_force)) .~ 0
        return System(eqs, t, particle_unknowns(connectors), [psys]; name)
    elseif wing.dynamics_type == RIGID_DYNAMICS
        num_twist_surfaces = length(wing.twist_surface_idxs)
        connectors = rigid_aero_connectors(num_twist_surfaces)
        eqs = [collect(connectors.force) .~ 0
               collect(connectors.moment) .~ 0]
        num_twist_surfaces > 0 && (eqs = [eqs; collect(connectors.twist_moment) .~ 0])
        return System(eqs, t, rigid_unknowns(connectors), [psys]; name)
    else
        error("Unknown dynamics_type $(wing.dynamics_type) for wing $wing_idx.")
    end
end
