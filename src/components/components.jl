# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Point/Segment component Systems — the single source of truth both backends
# assemble from. Each is a small MTK `System` with array-valued I/O whose force law
# is the shared kernels (`kernels.jl`); air density and wind enter through the same
# env functions, so no physics is duplicated. Environment is captured from `s` at
# build time; parameters + defaults come from `params.<kind>[idx]` (flat_params.jl).
#
# Aggregation uses explicit force *inputs* (not an MTK `Flow`/`connect` connector),
# so array-valued I/O (`pos(t)[1:3]`) is preserved: a point receives the summed
# force of its incident segments through `force_in`, and each backend supplies that
# sum its own way — NetworkDynamics sums each edge output into the vertex input; the
# monolith writes the sum as an explicit equation over the `@named` subsystems. This
# is the one interface both an ND-graph and an MTK-`connect`-free assembly can share.

"""
    point_io()

The array-valued I/O variables of a point component: state `pos`/`vel`, the
aggregated `force_in`/`mass_in` a point receives from its incident segments (summed
by the assembly), plus `tension_in`/`pulley_len_out` used only by pulley/winch
assembly (0 for a plain point). Returns `(vars, pos, vel, force_in, mass_in,
tension_in, pulley_len_out)`.
"""
function point_io()
    vars = @variables begin
        pos(t)[1:3]
        vel(t)[1:3]
        force_in(t)[1:3], [input = true]
        mass_in(t), [input = true]
        tension_in(t), [input = true]
        pulley_len_out(t), [output = true]
    end
    return vars, pos, vel, force_in, mass_in, tension_in, pulley_len_out
end

"""
    segment_io()

The array-valued I/O variables of a segment component: the two endpoints' `pos`/`vel`
read as inputs and the endpoint loads written as outputs (`src_force`/`src_mass`,
`dst_force`/`dst_mass`, positive force-on-point sign). Returns `(vars, src_pos,
src_vel, dst_pos, dst_vel, src_force, src_mass, dst_force, dst_mass)`.
"""
function segment_io()
    vars = @variables begin
        src_pos(t)[1:3], [input = true]
        src_vel(t)[1:3], [input = true]
        dst_pos(t)[1:3], [input = true]
        dst_vel(t)[1:3], [input = true]
        src_force(t)[1:3], [output = true]
        src_mass(t), [output = true]
        dst_force(t)[1:3], [output = true]
        dst_mass(t), [output = true]
    end
    return vars, src_pos, src_vel, dst_pos, dst_vel,
           src_force, src_mass, dst_force, dst_mass
end

"""
    point_acceleration(s, pos, vel, structural_force, mass, drag_coeff, area,
                       world_damping, wind_gnd)

World-frame acceleration of a point mass: the structural force gathered from its
segments plus its own aerodynamic drag and gravity, per unit `mass`, minus
world-frame damping. `structural_force` is the physical net force on the point
(positive sign). Shared by the connector [`DynamicPoint`](@ref) (monolith) and the
network vertex (ext) so the point physics lives in one place (D2); each backend
supplies `structural_force` in its own aggregation convention.
"""
function point_acceleration(s, pos, vel, structural_force, mass, drag_coeff, area,
                            world_damping, wind_gnd)
    wind = WindFactor(s.am, s.set.profile_law)
    rho = calc_rho(s.am, max(0.0, pos[3]))
    va = wind(pos[3]) .* wind_gnd .- vel
    drag = point_drag_force(va, rho, drag_coeff, area)
    gravity = [0.0, 0.0, -s.set.g_earth * mass]
    return (structural_force .+ drag .+ gravity) ./ mass .- world_damping .* vel
end

"""
    wing_structural_segment(sys_struct, idx)

Whether segment `idx` is an internal wing-structural link — both endpoints are wing
nodes. Such a segment carries no tether drag (its aerodynamic load is owned by the
wing's VSM), so the assembly gives it the drag-free segment type rather than passing
a per-segment drag coefficient.
"""
function wing_structural_segment(sys_struct, idx)
    seg = sys_struct.segments[idx]
    return sys_struct.points[seg.point_idxs[1]].is_wing_node &&
           sys_struct.points[seg.point_idxs[2]].is_wing_node
end

"""
    segment_endpoint_loads(s, src_pos, src_vel, dst_pos, dst_vel, unit_stiffness,
                           unit_damping, compression_frac, l0, diameter, density,
                           cd_tether, wind_gnd; with_drag=true)

Physical loads a segment exerts: `(force_on_src, force_on_dst, half_mass, spring)`,
with forces in the *positive* (force-on-point) sign convention, `half_mass` the mass
each endpoint carries, and `spring` the signed scalar spring-damper tension along
the axis (positive in tension). Spring acts along the axis (opposite signs at the
ends); tether drag is split equally and omitted entirely when `with_drag=false` (the
drag-free wing-structural segment type, so `cd_tether`/`wind_gnd` are then unused).
Shared by [`SpringDamperSegment`](@ref) (monolith) and the network edge (ext), which
uses these signs directly and the scalar `spring` for pulley/winch tension aggregation.
"""
function segment_endpoint_loads(s, src_pos, src_vel, dst_pos, dst_vel,
                                unit_stiffness, unit_damping, compression_frac,
                                l0, diameter, density, cd_tether, wind_gnd;
                                with_drag = true)
    _, len, unit_vec, spring_vel =
        segment_geometry(src_pos, dst_pos, src_vel, dst_vel)
    spring = segment_spring_force(len, l0, spring_vel, unit_stiffness,
                                  unit_damping, compression_frac)
    spring_vec = spring .* unit_vec
    half_mass = segment_half_mass(l0, diameter, density)
    with_drag || return spring_vec, -spring_vec, half_mass, spring
    wind = WindFactor(s.am, s.set.profile_law)
    seg_pos_z = 0.5 * (src_pos[3] + dst_pos[3])
    rho = calc_rho(s.am, max(0.0, seg_pos_z))
    seg_vel = 0.5 .* (src_vel .+ dst_vel)
    va = wind(seg_pos_z) .* wind_gnd .- seg_vel
    drag = segment_perp_drag(va, unit_vec, rho, cd_tether, len * diameter)
    return spring_vec .+ 0.5 .* drag, -spring_vec .+ 0.5 .* drag, half_mass, spring
end

"""
    ground_wind_vec(params)

Ground-level wind vector as symbolic parameters — `params.set.wind_vec`, with a tiny
x-axis fallback when it is exactly zero (avoids normalize-by-zero, matching
`scalar_eqs!`). Returns a length-3 `Vector{Num}`.
"""
function ground_wind_vec(params)
    raw = collect(params.set.wind_vec)
    near_zero = sum(abs2, raw) < 1e-20
    fallback = (1e-10, 0.0, 0.0)
    return [ifelse(near_zero, fallback[k], raw[k]) for k in 1:3]
end

"""
    DynamicPoint(s, params, idx; name)

Particle vertex: a point mass integrating `pos`/`vel` under the net force gathered
at its `node`, its own aerodynamic drag, gravity, and world-frame damping. This is
the free-tether particle and the base for PARTICLE_DYNAMICS wing nodes; the wing
case adds an aero-force input on top of this same integrator (added with wing-node
wiring, Phase 4/5) rather than duplicating it. Translational mass is `extra_mass`
plus the incident segments' aggregated half-masses (`node.mass`). Its parameters are
read from `params.points[idx]` (mass/drag/area/damping) plus the computed
[`ground_wind_vec`](@ref) — the same param+defaults source both backends assemble
from.
"""
function DynamicPoint(s, params, idx; name)
    vars, pos, vel, force_in, mass_in, _, pulley_len_out = point_io()
    point = params.points[idx]
    mass = point.extra_mass + mass_in
    accel = point_acceleration(s, collect(pos), collect(vel), collect(force_in), mass,
        point.drag_coeff, point.area, collect(point.world_frame_damping),
        ground_wind_vec(params))
    eqs = [
        D.(collect(pos)) .~ collect(vel)
        D.(collect(vel)) .~ accel
        pulley_len_out ~ 0.0
    ]
    return System(eqs, t, vars, param_unknowns(params); name)
end

"""
    StaticPoint(s, params, idx; name)

Ground-anchored vertex with no dynamic state: `pos` is pinned to `params.points[idx]`'s
`pos_w` and `vel` is zero. Its `force_in`/`mass_in`/`tension_in` inputs are declared
(so segments may deliver to it) but ignored — the point does not move — matching the
STATIC branch of `point_eqs!`.
"""
function StaticPoint(s, params, idx; name)
    vars, pos, vel, _, _, _, pulley_len_out = point_io()
    eqs = [
        collect(pos) .~ collect(params.points[idx].pos_w)
        collect(vel) .~ zeros(3)
        pulley_len_out ~ 0.0
    ]
    return System(eqs, t, vars, param_unknowns(params); name)
end

"""
    SpringDamperSegment(s, params, idx; name)

Stateless edge reading its two endpoints' `pos`/`vel` and writing the spring-damper
force + tether drag on each (`src_force`/`dst_force`) and each endpoint's half-mass
(`src_mass`/`dst_mass`), in the **positive** force-on-point sign — the assembly sums
these into each point's `force_in`/`mass_in` (ND edge→vertex, or the monolith's
explicit sum). A [`wing_structural_segment`](@ref) has no tether drag: it drops the
drag term (a distinct compiled type) instead of reading `cd_tether`, so drag stays a
single global setting. Its spring/geometry parameters are read from
`params.segments[idx]`; `cd_tether` from `params.set` and the ground wind from
[`ground_wind_vec`](@ref).
"""
function SpringDamperSegment(s, params, idx; name)
    vars, src_pos, src_vel, dst_pos, dst_vel,
        src_force, src_mass, dst_force, dst_mass = segment_io()
    seg = params.segments[idx]
    with_drag = !wing_structural_segment(params.reg.sys_struct, idx)
    cd_tether = with_drag ? params.set.cd_tether : 0.0
    wind = with_drag ? ground_wind_vec(params) : zeros(3)
    force_on_src, force_on_dst, half_mass, _ = segment_endpoint_loads(
        s, collect(src_pos), collect(src_vel), collect(dst_pos), collect(dst_vel),
        seg.unit_stiffness, seg.unit_damping, seg.compression_frac, seg.l0,
        seg.diameter, seg.density, cd_tether, wind; with_drag)
    eqs = [
        collect(src_force) .~ force_on_src
        collect(dst_force) .~ force_on_dst
        src_mass ~ half_mass
        dst_mass ~ half_mass
    ]
    return System(eqs, t, vars, param_unknowns(params); name)
end
