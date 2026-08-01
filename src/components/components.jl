# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Point/Segment component Systems — the single source of truth (D2) both backends
# assemble from. Each is a small MTK `System` whose force law is the shared
# kernels (`kernels.jl`); air density and wind enter through the same env
# functions the monolith uses, so no physics is duplicated. Environment is
# captured from the model `s` at build time.
#
# The `Node` connector keeps array-valued variables (clean interface, few
# symbols). Component *equations* collect those to `Vector{Num}` and use broadcast
# `~`, because the reductions in the kernels (`smooth_norm`) and the `Vector{Num}`
# gravity term do not compose with raw symbolic-array terms (O3 policy: scalarize
# only where the array form errors).

"""
    Node()

Mechanical connector joining points and segments. Position and velocity are
`across` (potential) variables — equal for everything joined at the node — while
force and mass are `Flow` variables that sum over a connection, giving a point the
net force and total half-segment mass of its incident segments. The `Flow` sum is
the monolith's counterpart to the network backend's edge→vertex sum-aggregation
(D3). Components are scalarized (`pos1..3`, not `pos[1:3]`) because array potential
variables do not expand through `connect` in the current MTK (O3).
"""
@connector function Node(; name)
    vars = @variables begin
        pos1(t)
        pos2(t)
        pos3(t)
        vel1(t)
        vel2(t)
        vel3(t)
        force1(t), [connect = Flow]
        force2(t), [connect = Flow]
        force3(t), [connect = Flow]
        mass(t), [connect = Flow]
    end
    return System(Equation[], t, vars, []; name)
end

"""
    node_state(node)

Collect a [`Node`](@ref)'s scalarized components into `(pos, vel, force)` vectors
so the component force laws can work with 3-vectors while the connector stays
scalar (O3). `force` is the endpoint force to write (segments) or the aggregated
force gathered by a point.
"""
node_state(node) = ([node.pos1, node.pos2, node.pos3],
                    [node.vel1, node.vel2, node.vel3],
                    [node.force1, node.force2, node.force3])

"""
    DynamicPoint(s; name)

Particle vertex: a point mass integrating `pos`/`vel` under the net force gathered
at its `node`, its own aerodynamic drag, gravity, and world-frame damping. This is
the free-tether particle and the base for PARTICLE_DYNAMICS wing nodes; the wing
case adds an aero-force input on top of this same integrator (added with wing-node
wiring, Phase 4/5) rather than duplicating it. Translational mass is `extra_mass`
plus the incident segments' aggregated half-masses (`node.mass`).
"""
function DynamicPoint(s; name)
    @named node = Node()
    pars = @parameters begin
        extra_mass = 0.0
        drag_coeff = 0.0
        area = 0.0
        world_damping[1:3] = zeros(3)
        wind_gnd[1:3] = zeros(3)
    end
    g_earth = s.set.g_earth
    wind = WindFactor(s.am, s.set.profile_law)
    pos, vel, force = node_state(node)
    rho = calc_rho(s.am, max(0.0, pos[3]))
    va = wind(pos[3]) .* collect(wind_gnd) .- vel
    drag = point_drag_force(va, rho, drag_coeff, area)
    mass = extra_mass + node.mass
    gravity = [0.0, 0.0, -g_earth * mass]
    accel = (force .+ drag .+ gravity) ./ mass .- collect(world_damping) .* vel
    eqs = [
        D.(pos) .~ vel
        D.(vel) .~ accel
    ]
    return System(eqs, t, [], pars; systems = [node], name)
end

"""
    StaticPoint(s; name)

Ground-anchored vertex with no dynamic state: `pos` is pinned to `pos_w` and `vel`
is zero. Any force gathered at its `node` is absorbed as a reaction (the point does
not move), matching the STATIC branch of `point_eqs!`.
"""
function StaticPoint(s; name)
    @named node = Node()
    pars = @parameters pos_w[1:3] = zeros(3)
    pos, vel, _ = node_state(node)
    eqs = [
        pos .~ collect(pos_w)
        vel .~ zeros(3)
    ]
    return System(eqs, t, [], pars; systems = [node], name)
end

"""
    SpringDamperSegment(s; name)

Stateless edge between two `node`s carrying the spring-damper force and tether
drag. Writes the force on each endpoint (`+spring + ½ drag` on the source,
`−spring + ½ drag` on the destination) and each endpoint's half-mass into the
`Flow` variables, so the connected points aggregate force and mass. A rest length
`l0` and `unit_stiffness = 0` reproduce the wing-structural special cases without a
separate component (see §5.2 of the plan).
"""
function SpringDamperSegment(s; name)
    @named src = Node()
    @named dst = Node()
    pars = @parameters begin
        unit_stiffness = 0.0
        unit_damping = 0.0
        compression_frac = 0.1
        l0 = 1.0
        diameter = 0.0
        density = 0.0
        cd_tether = 1.0
        wind_gnd[1:3] = zeros(3)
    end
    wind = WindFactor(s.am, s.set.profile_law)
    src_pos, src_vel, src_force = node_state(src)
    dst_pos, dst_vel, dst_force = node_state(dst)
    _, len, unit_vec, spring_vel =
        segment_geometry(src_pos, dst_pos, src_vel, dst_vel)
    spring = segment_spring_force(len, l0, spring_vel, unit_stiffness,
                                  unit_damping, compression_frac)
    spring_vec = spring .* unit_vec
    seg_pos_z = 0.5 * (src_pos[3] + dst_pos[3])
    rho = calc_rho(s.am, max(0.0, seg_pos_z))
    seg_vel = 0.5 .* (src_vel .+ dst_vel)
    va = wind(seg_pos_z) .* collect(wind_gnd) .- seg_vel
    drag = segment_perp_drag(va, unit_vec, rho, cd_tether, len * diameter)
    half_mass = segment_half_mass(l0, diameter, density)
    # `Flow` variables sum to zero at a `connect`, so a point gathers −Σ(edge flow).
    # The edge therefore contributes the *negative* of the force on / half-mass at
    # each endpoint, leaving the point with +force and +mass. (The network backend
    # aggregates edge outputs into the vertex input directly, so its wrapper applies
    # the opposite sign — the shared physics lives in the kernels, D2.)
    eqs = [
        src_force .~ -spring_vec .- 0.5 .* drag
        dst_force .~ spring_vec .- 0.5 .* drag
        src.mass ~ -half_mass
        dst.mass ~ -half_mass
    ]
    return System(eqs, t, [], pars; systems = [src, dst], name)
end
