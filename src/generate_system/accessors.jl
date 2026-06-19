# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Symbolic accessors: registered with SystemStructure (UnionAll)
# because @register_symbolic erases type parameters, causing
# duplicate methods if registered per concrete wing type.

# ==================== GENERIC ACCESSORS ==================== #
# These access SystemStructure fields shared by all wing types:
# points, segments, pulleys, winches, tethers, settings, BaseWing
# fields.

# ---- Points ----
get_pos_w(sys::SystemStructure, idx::Int64) =
    sys.points[idx].pos_w
@register_array_symbolic get_pos_w(
    sys::SystemStructure, idx::Int64) begin
    size = (3,)
    eltype = SimFloat
end

# ---- Rigid bodies ----
get_body_R_b_to_p(sys::SystemStructure, idx::Int64) =
    sys.rigid_bodies[idx].R_b_to_p
@register_array_symbolic get_body_R_b_to_p(
    sys::SystemStructure, idx::Int64) begin
    size = (3, 3)
    eltype = SimFloat
end

# ---- Elastic joints ----
# Restoring force/moment from one joint stiffness, given the relative DOF `Δ`.
# `kind`: 1=axial, 2=shear, 3=torsion, 4=bending. A `Real` stiffness is the
# linear law `k·Δ`; a callable interpolation gives the (possibly saturating)
# force directly as `f(Δ)`. The `stiffness_force` function barrier dispatches on
# the concrete stiffness type so the RHS stays allocation-free; `Δ` flows through
# as a ForwardDiff `Dual` for the Jacobian.
stiffness_force(k::Real, Δ) = k * Δ
stiffness_force(interpolation, Δ) = interpolation(Δ)

joint_stiffness(joint::ElasticJoint, kind::Int) =
    kind == 1 ? joint.stiffness_axial :
    kind == 2 ? joint.stiffness_shear :
    kind == 3 ? joint.stiffness_torsion : joint.stiffness_bending

get_joint_force(sys::SystemStructure, idx::Int64, kind::Int, Δ) =
    stiffness_force(joint_stiffness(sys.elastic_joints[idx], kind), Δ)
@register_symbolic get_joint_force(
    sys::SystemStructure, idx::Int64, kind::Int, Δ)

# ---- Winches ----
get_set_value(sys::SystemStructure, idx::Int64) =
    sys.winches[idx].set_value
@register_symbolic get_set_value(
    sys::SystemStructure, idx::Int64)
get_winch_gear_ratio(sys::SystemStructure, idx::Int64) =
    sys.winches[idx].gear_ratio
@register_symbolic get_winch_gear_ratio(
    sys::SystemStructure, idx::Int64)
get_winch_drum_radius(sys::SystemStructure, idx::Int64) =
    sys.winches[idx].drum_radius
@register_symbolic get_winch_drum_radius(
    sys::SystemStructure, idx::Int64)
get_winch_f_coulomb(sys::SystemStructure, idx::Int64) =
    sys.winches[idx].f_coulomb
@register_symbolic get_winch_f_coulomb(
    sys::SystemStructure, idx::Int64)
get_winch_c_vf(sys::SystemStructure, idx::Int64) =
    sys.winches[idx].c_vf
@register_symbolic get_winch_c_vf(
    sys::SystemStructure, idx::Int64)
get_winch_inertia_total(sys::SystemStructure, idx::Int64) =
    sys.winches[idx].inertia_total
@register_symbolic get_winch_inertia_total(
    sys::SystemStructure, idx::Int64)
get_winch_friction_epsilon(sys::SystemStructure, idx::Int64) =
    sys.winches[idx].friction_epsilon
@register_symbolic get_winch_friction_epsilon(
    sys::SystemStructure, idx::Int64)

const ZERO_WIND_FALLBACK = KVec3(1e-10, 0.0, 0.0)
"""
    get_wind_vec(sys)

Ground wind vector [m/s] from settings, with a tiny x-axis fallback for
exactly-zero wind (avoids normalize-by-zero). Read by the flat `wind_vec`
parameter sync (no longer a registered symbolic — it is the computed reader).
"""
function get_wind_vec(sys::SystemStructure)
    wv = sys.set.wind_vec
    if wv[1]^2 + wv[2]^2 + wv[3]^2 < 1e-20
        return ZERO_WIND_FALLBACK
    end
    return wv
end
