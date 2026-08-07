# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# The component physics — the single source of truth both backends assemble from.
# The force laws themselves are the shared kernels (`kernels.jl`); this file adds the
# symbolic layer over them: the parameter reads (`params.<kind>[idx]`,
# flat_params.jl), the environment (captured from `s` at build time, so air density
# and wind enter identically everywhere), and the per-component expression bundles.
#
# The first half is pure expression builders that the monolith's equation generators
# in `src/generate_system/` call directly. The second half wraps them into the small
# MTK `System`s the `ScheduledBackend` compiles one kernel from, each declaring
# exactly the inputs it reads and the outputs someone wires — no uniform width to pad
# to and nothing baked per endpoint orientation.
#
# Aggregation is by explicit force *inputs*, not an MTK `Flow`/`connect` connector,
# which is what keeps array-valued I/O (`pos(t)[1:3]`) intact: a point receives the
# summed force of its incident segments through `force_in`, and each backend supplies
# that sum its own way — the scheduled runtime gathers it from the wiring, the
# monolith writes it as an explicit equation.

"""`vector` with its component along the unit `axis` removed."""
remove_along(vector, axis) = vector .- (vector ⋅ axis) .* axis

"""Only `vector`'s component along the unit `axis`."""
keep_along(vector, axis) = (vector ⋅ axis) .* axis

"""
    point_acceleration(s, pos, vel, structural_force, mass, drag_coeff, area,
                       world_damping, wind_gnd)

World-frame acceleration of a point mass: the structural force gathered from its
segments plus its own aerodynamic drag and gravity, per unit `mass`, minus
world-frame damping. `structural_force` is the physical net force on the point
(positive sign). Shared by the monolith's `point_eqs!` and the [`Particle`](@ref)
component, so the point physics lives in one place; each backend supplies
`structural_force` in its own aggregation convention.
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
    point_particle_params(params, idx)

The shared DYNAMIC-particle parameters read from `params.points[idx]` — mass, drag,
area and world-frame damping as the point's own struct fields — plus the computed
ground wind `wind_gnd` ([`ground_wind_vec`](@ref)). Returns a named tuple consumed by
[`dynamic_point_dynamics`](@ref); each read registers the parameter on `params`.
"""
function point_particle_params(params, idx)
    point = params.points[idx]
    wind_gnd = ground_wind_vec(params)
    return (; extra_mass = point.extra_mass, drag_coeff = point.drag_coeff,
            area = point.area, world_damping = point.world_frame_damping,
            fix_sphere = point.fix_sphere, fix_static = point.fix_static, wind_gnd)
end

"""
    confined_derivatives(pos, vel, accel, pars)

`(D(pos), D(vel))` for a point that may be pinned: `fix_static` freezes it where it
is, and `fix_sphere` confines it to a sphere about the world origin by keeping only
the radial part of its velocity and acceleration. Matches the pair of `ifelse`s in
`point_eqs!`.
"""
function confined_derivatives(pos, vel, accel, pars)
    axis = collect(smooth_normalize(collect(pos)))
    velocity = ifelse.(pars.fix_sphere == true, keep_along(collect(vel), axis),
                       collect(vel))
    acceleration = ifelse.(pars.fix_sphere == true, keep_along(accel, axis), accel)
    frozen = pars.fix_static == true
    return ifelse.(frozen, zeros(3), velocity), ifelse.(frozen, zeros(3), acceleration)
end

"""
    dynamic_point_dynamics(s, pos, vel, force, mass, pars)

Shared body of the DYNAMIC point/pulley vertices: `D(pos)=vel` and
`D(vel)=point_acceleration(...)` from the shared kernel, reading the point's
drag/damping/wind parameters `pars` (a [`point_particle_params`](@ref) named tuple).
"""
function dynamic_point_dynamics(s, pos, vel, force, mass, pars)
    accel = point_acceleration(s, collect(pos), collect(vel), collect(force),
        mass, pars.drag_coeff, pars.area, collect(pars.world_damping),
        collect(pars.wind_gnd))
    velocity, acceleration = confined_derivatives(pos, vel, accel, pars)
    return [D.(collect(pos)) .~ velocity; D.(collect(vel)) .~ acceleration]
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
    segment_spring_params(params, idx; with_drag=true)

The spring-damper parameters read from `params.segments[idx]` (stiffness, damping,
compression fraction, diameter, density as the segment's own struct fields), plus
the global tether drag `cd_tether` (`params.set.cd_tether`) and the ground wind
`wind_gnd` ([`ground_wind_vec`](@ref)). With `with_drag=false` (the
[`wing_structural_segment`](@ref) edge) `cd_tether` is a literal `0` and unused.
`nonlinear` marks a callable `unit_stiffness` force law. Returns
`(spring_named_tuple, wind_gnd)`; each read registers the parameter on `params`.
"""
function segment_spring_params(params, idx; with_drag = true)
    seg = params.segments[idx]
    cd_tether = with_drag ? params.set.cd_tether : 0.0
    wind_gnd = ground_wind_vec(params)
    nonlinear = !(params.reg.sys_struct.segments[idx].unit_stiffness isa Real)
    spring = (; unit_stiffness = seg.unit_stiffness, unit_damping = seg.unit_damping,
              compression_frac = seg.compression_frac, diameter = seg.diameter,
              density = seg.density, cd_tether, nonlinear)
    return spring, wind_gnd
end

"""
    segment_load_terms(s, src_pos, src_vel, dst_pos, dst_vel, unit_stiffness,
                       unit_damping, compression_frac, l0, diameter, density,
                       cd_tether, wind_gnd; with_drag=true, nonlinear=false)

Every load term a segment produces, as a named tuple: the geometry (`len`,
`unit_vec`), the signed scalar spring-damper tension `spring` and its vector
`spring_vec`, the `half_mass` and `half_drag` each endpoint carries, and the total
`force_on_src`/`force_on_dst` in the positive force-on-point sign. With `nonlinear`
the `unit_stiffness` is a callable force law of strain
([`segment_nonlinear_force`](@ref)) rather than a linear rate; with
`with_drag = false` the tether drag is dropped entirely, so `cd_tether`/`wind_gnd`
are unused.
"""
function segment_load_terms(s, src_pos, src_vel, dst_pos, dst_vel,
                            unit_stiffness, unit_damping, compression_frac,
                            l0, diameter, density, cd_tether, wind_gnd;
                            with_drag = true, nonlinear = false)
    _, len, unit_vec, spring_vel =
        segment_geometry(src_pos, dst_pos, src_vel, dst_vel)
    spring = nonlinear ?
        segment_nonlinear_force(len, l0, spring_vel, unit_stiffness, unit_damping) :
        segment_spring_force(len, l0, spring_vel, unit_stiffness, unit_damping,
                             compression_frac)
    spring_vec = spring .* unit_vec
    half_mass = segment_half_mass(l0, diameter, density)
    half_drag = zeros(Num, 3)
    if with_drag
        wind = WindFactor(s.am, s.set.profile_law)
        seg_pos_z = 0.5 * (src_pos[3] + dst_pos[3])
        rho = calc_rho(s.am, max(0.0, seg_pos_z))
        seg_vel = 0.5 .* (src_vel .+ dst_vel)
        va = wind(seg_pos_z) .* wind_gnd .- seg_vel
        half_drag = 0.5 .*
            segment_perp_drag(va, unit_vec, rho, cd_tether, len * diameter)
    end
    return (; len, unit_vec, spring, spring_vec, half_mass, half_drag,
            force_on_src = spring_vec .+ half_drag,
            force_on_dst = .-spring_vec .+ half_drag)
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
    timoshenko_element_wrench(joint, params; frame, theta_a, theta_b, force_a, force_b,
        moment_a, moment_b, pos_a, R_a, com_a, com_vel_a, omega_a_w,
        pos_b, R_b, com_b, com_vel_b, omega_b_w)

Corotational Timoshenko element wrench shared by both backends (the monolith
`timoshenko_joint_eqs!` loop body and the scheduled joint component), so the
beam-element physics lives in one place. Given the two nodes' world poses (`pos`,
`R_b_to_w`, `com`, `com_vel`, world spin `omega_w`) and the joint's rest
geometry/rigidities, it builds the
element frame and per-node deformations, evaluates the consistent Timoshenko stiffness
(axial, torsion, two bending planes with shear reduction `Φ`) and damping, and returns
`(tear_eqs, force_on_a, moment_on_a, force_on_b, moment_on_b)` — the restoring wrench on
each node (world frame, transported to each COM). `frame`/`theta_a`/`theta_b`/`force_a`/
`force_b`/`moment_a`/`moment_b` are caller-supplied **torn** variables (array slices for
the monolith, standalone vars for the edge) so the reused frame/force subtrees are not
re-embedded; `tear_eqs` binds them.
"""
function timoshenko_element_wrench(joint, params;
        frame, theta_a, theta_b, force_a, force_b, moment_a, moment_b,
        pos_a, R_a, com_a, com_vel_a, omega_a_w,
        pos_b, R_b, com_b, com_vel_b, omega_b_w)
    j = joint.idx
    jp = params.timoshenko_joints[j]
    anchor_a = collect(jp.anchor_a_b)
    anchor_b = collect(jp.anchor_b_b)
    Ra = collect(R_a)
    Rb = collect(R_b)
    x_a = collect(pos_a) .+ Ra * anchor_a
    x_b = collect(pos_b) .+ Rb * anchor_b
    e1, e2, e3, len = timoshenko_element_frame(x_a, x_b, Ra)
    frame_expr = [e1[1] e2[1] e3[1];
                  e1[2] e2[2] e3[2];
                  e1[3] e2[3] e3[3]]
    element_frame = collect(frame)
    L0 = jp.rest_length
    R_a_rel0 = collect(jp.R_a_rel0)
    R_b_rel0 = collect(jp.R_b_rel0)
    Da = (element_frame' * Ra) * R_a_rel0'
    Db = (element_frame' * Rb) * R_b_rel0'
    θ_a_expr = [0.5 * (Da[3, 2] - Da[2, 3]),
                0.5 * (Da[1, 3] - Da[3, 1]),
                0.5 * (Da[2, 1] - Da[1, 2])]
    θ_b_expr = [0.5 * (Db[3, 2] - Db[2, 3]),
                0.5 * (Db[1, 3] - Db[3, 1]),
                0.5 * (Db[2, 1] - Db[1, 2])]
    θ_a = collect(theta_a)
    θ_b = collect(theta_b)
    δ = len - L0
    kshear = jp.shear_coeff
    ε = δ / L0
    κt = (θ_b[1] - θ_a[1]) / L0
    κy = (θ_b[2] - θ_a[2]) / L0
    κz = (θ_b[3] - θ_a[3]) / L0
    γy = 0.5 * (θ_a[2] + θ_b[2])
    γz = 0.5 * (θ_a[3] + θ_b[3])
    EA_eff = timoshenko_rigidity(joint, params, :EA, ε)
    GJ_eff = timoshenko_rigidity(joint, params, :GJ, κt)
    EIy_eff = timoshenko_rigidity(joint, params, :EIy, κy)
    EIz_eff = timoshenko_rigidity(joint, params, :EIz, κz)
    GAy_eff = timoshenko_rigidity(joint, params, :GA, γy)
    GAz_eff = timoshenko_rigidity(joint, params, :GA, γz)
    Φy = 12 * EIy_eff / (kshear * GAy_eff * L0^2)
    Φz = 12 * EIz_eff / (kshear * GAz_eff * L0^2)
    by = EIy_eff / (L0 * (1 + Φy))
    bz = EIz_eff / (L0 * (1 + Φz))
    shy = 6 * EIy_eff / (L0^2 * (1 + Φy))
    shz = 6 * EIz_eff / (L0^2 * (1 + Φz))
    Mt = GJ_eff / L0
    f_axial = EA_eff * δ / L0
    F_a_local = [f_axial, -shz * (θ_a[3] + θ_b[3]), shy * (θ_a[2] + θ_b[2])]
    M_a_local = [-Mt * (θ_a[1] - θ_b[1]),
                 -by * ((4 + Φy) * θ_a[2] + (2 - Φy) * θ_b[2]),
                 -bz * ((4 + Φz) * θ_a[3] + (2 - Φz) * θ_b[3])]
    F_b_local = [-f_axial, shz * (θ_a[3] + θ_b[3]), -shy * (θ_a[2] + θ_b[2])]
    M_b_local = [-Mt * (θ_b[1] - θ_a[1]),
                 -by * ((2 - Φy) * θ_a[2] + (4 + Φy) * θ_b[2]),
                 -bz * ((2 - Φz) * θ_a[3] + (4 + Φz) * θ_b[3])]
    ω_a_w = collect(omega_a_w)
    ω_b_w = collect(omega_b_w)
    vel_a = collect(com_vel_a) .+ (ω_a_w × (x_a .- collect(com_a)))
    vel_b = collect(com_vel_b) .+ (ω_b_w × (x_b .- collect(com_b)))
    Δv = vel_b .- vel_a
    Δω = ω_b_w .- ω_a_w
    c_t = jp.damping_trans
    c_r = jp.damping_rot
    tear_eqs = [
        vec(collect(frame)) ~ vec(frame_expr)
        collect(theta_a) ~ θ_a_expr
        collect(theta_b) ~ θ_b_expr
        collect(force_a) ~ element_frame * F_a_local .+ c_t .* Δv
        collect(force_b) ~ element_frame * F_b_local .- c_t .* Δv
        collect(moment_a) ~ element_frame * M_a_local .+ c_r .* Δω
        collect(moment_b) ~ element_frame * M_b_local .- c_r .* Δω
    ]
    force_on_a = collect(force_a)
    force_on_b = collect(force_b)
    moment_on_a = (x_a .- collect(com_a)) × force_on_a .+ collect(moment_a)
    moment_on_b = (x_b .- collect(com_b)) × force_on_b .+ collect(moment_b)
    return (; tear_eqs, force_on_a, moment_on_a, force_on_b, moment_on_b)
end

"""
    elastic_joint_wrench(joint, params; force_w, torque_w, pos_a, R_a, com_a, com_vel_a,
        omega_a_w, pos_b, R_b, com_b, com_vel_b, omega_b_w)

Lumped 6-DOF `ElasticJoint` restoring wrench shared by both backends (the monolith
`joint_eqs!` loop body and the scheduled elastic-joint component). From the relative
pose of the two anchors (in body A's frame) it builds the per-DOF restoring
force/torque (axial,
shear, torsion, bending stiffness + damping) and returns
`(tear_eqs, force_on_a, moment_on_a, force_on_b, moment_on_b)` — the equal-and-opposite
wrench transported to each COM. `force_w`/`torque_w` are the caller's **torn** world-frame
wrench variables (array slices for the monolith, standalone vars for the edge).
"""
function elastic_joint_wrench(joint, params; force_w, torque_w,
        pos_a, R_a, com_a, com_vel_a, omega_a_w,
        pos_b, R_b, com_b, com_vel_b, omega_b_w)
    j = joint.idx
    jp = params.elastic_joints[j]
    Ra = collect(R_a)
    Rb = collect(R_b)
    anchor_a = collect(jp.anchor_a_b)
    anchor_b = collect(jp.anchor_b_b)
    ca = collect(com_a)
    cb = collect(com_b)
    pos_anchor_a = collect(pos_a) .+ Ra * anchor_a
    pos_anchor_b = collect(pos_b) .+ Rb * anchor_b
    rest_offset = collect(jp.rest_offset_a)
    R_rel0 = collect(jp.R_rel0)
    Δr_a = Ra' * (pos_anchor_b .- pos_anchor_a) .- rest_offset
    R_rel = R_rel0' * (Ra' * Rb)
    Δθ_a = [0.5 * (R_rel[3, 2] - R_rel[2, 3]),
            0.5 * (R_rel[1, 3] - R_rel[3, 1]),
            0.5 * (R_rel[2, 1] - R_rel[1, 2])]
    ω_a_w = collect(omega_a_w)
    ω_b_w = collect(omega_b_w)
    vel_anchor_a = collect(com_vel_a) .+ (ω_a_w × (pos_anchor_a .- ca))
    vel_anchor_b = collect(com_vel_b) .+ (ω_b_w × (pos_anchor_b .- cb))
    Δv_a = Ra' * (vel_anchor_b .- vel_anchor_a)
    Δω_a = Ra' * (ω_b_w .- ω_a_w)
    damp_trans = jp.damping_trans
    damp_rot = jp.damping_rot
    force_a = [
        -joint_stiffness_term(joint, params, 1, Δr_a[1]) - damp_trans * Δv_a[1],
        -joint_stiffness_term(joint, params, 2, Δr_a[2]) - damp_trans * Δv_a[2],
        -joint_stiffness_term(joint, params, 2, Δr_a[3]) - damp_trans * Δv_a[3]]
    torque_a = [
        -joint_stiffness_term(joint, params, 3, Δθ_a[1]) - damp_rot * Δω_a[1],
        -joint_stiffness_term(joint, params, 4, Δθ_a[2]) - damp_rot * Δω_a[2],
        -joint_stiffness_term(joint, params, 4, Δθ_a[3]) - damp_rot * Δω_a[3]]
    tear_eqs = [collect(force_w) ~ Ra * force_a; collect(torque_w) ~ Ra * torque_a]
    force_on_b = collect(force_w)
    torque_on_b = collect(torque_w)
    arm_a = pos_anchor_a .- ca
    arm_b = pos_anchor_b .- cb
    force_on_a = .-force_on_b
    moment_on_a = arm_a × force_on_a .- torque_on_b
    moment_on_b = arm_b × force_on_b .+ torque_on_b
    return (; tear_eqs, force_on_a, moment_on_a, force_on_b, moment_on_b)
end

"""
    beam_hermite_ride_expressions(joint, params, point_idx; pos_a, R_a, com_a, com_vel_a,
        omega_a_w, pos_b, R_b, com_b, com_vel_b, omega_b_w)

Kinematics of a point riding `joint`'s corotational cubic-Hermite centerline at the
point's `beam_frac`, shared by both backends. From the two end bodies' world poses it
builds the element frame, the two nodes' chord-relative rotations, the transverse Hermite
deflection (+ a frame-carried `beam_offset_b`) and returns `(pos_point, vel_point, sfrac,
ride_velocity)` — the ride position, its rigid-blend velocity, the axial fraction that
splits any force at the point onto the two end bodies (`(1−sfrac)` to A, `sfrac` to B), and
`ride_velocity(p)`: the rigid-blend velocity as a function of a (possibly **torn**) ride
position `p`, so a backend can tear `pos_point` first and avoid re-embedding the heavy
element-frame subtree in the velocity.
"""
function beam_hermite_ride_expressions(joint, params, point_idx;
        pos_a, R_a, com_a, com_vel_a, omega_a_w,
        pos_b, R_b, com_b, com_vel_b, omega_b_w,
        frame = nothing, theta_a = nothing, theta_b = nothing)
    jp = params.timoshenko_joints[joint.idx]
    Ra = collect(R_a)
    Rb = collect(R_b)
    x_a = collect(pos_a) .+ Ra * collect(jp.anchor_a_b)
    x_b = collect(pos_b) .+ Rb * collect(jp.anchor_b_b)
    e1, e2, e3, beam_len = timoshenko_element_frame(x_a, x_b, Ra)
    frame_expr = [e1[1] e2[1] e3[1];
                  e1[2] e2[2] e3[2];
                  e1[3] e2[3] e3[3]]
    element_frame = frame === nothing ? frame_expr : collect(frame)
    Da = (element_frame' * Ra) * collect(jp.R_a_rel0)'
    Db = (element_frame' * Rb) * collect(jp.R_b_rel0)'
    θ_a_expr = [0.5 * (Da[3, 2] - Da[2, 3]), 0.5 * (Da[1, 3] - Da[3, 1]),
                0.5 * (Da[2, 1] - Da[1, 2])]
    θ_b_expr = [0.5 * (Db[3, 2] - Db[2, 3]), 0.5 * (Db[1, 3] - Db[3, 1]),
                0.5 * (Db[2, 1] - Db[1, 2])]
    θ_a = theta_a === nothing ? θ_a_expr : collect(theta_a)
    θ_b = theta_b === nothing ? θ_b_expr : collect(theta_b)
    sfrac = params.points[point_idx].beam_frac
    ec1 = element_frame[:, 1]
    ec2 = element_frame[:, 2]
    ec3 = element_frame[:, 3]
    N2 = beam_len * (sfrac - 2sfrac^2 + sfrac^3)
    N4 = beam_len * (-sfrac^2 + sfrac^3)
    v_defl = N2 * θ_a[3] + N4 * θ_b[3]
    w_defl = -(N2 * θ_a[2] + N4 * θ_b[2])
    x_center = x_a .+ (sfrac * beam_len) .* ec1 .+ v_defl .* ec2 .+ w_defl .* ec3
    pos_point = x_center .+ element_frame * collect(params.points[point_idx].beam_offset_b)
    ω_a_w = collect(omega_a_w)
    ω_b_w = collect(omega_b_w)
    cvel_a = collect(com_vel_a)
    cvel_b = collect(com_vel_b)
    ca = collect(com_a)
    cb = collect(com_b)
    ride_velocity(p) = (1 - sfrac) .* (cvel_a .+ (ω_a_w × (collect(p) .- ca))) .+
        sfrac .* (cvel_b .+ (ω_b_w × (collect(p) .- cb)))
    tear_eqs = Equation[]
    if frame !== nothing
        append!(tear_eqs, vec(collect(frame)) .~ vec(frame_expr))
        append!(tear_eqs, collect(theta_a) .~ θ_a_expr)
        append!(tear_eqs, collect(theta_b) .~ θ_b_expr)
    end
    return (; pos_point, vel_point = ride_velocity(pos_point), sfrac,
            ride_velocity, tear_eqs)
end

"""
    rigid_body_pose_expressions(force_w, moment_w, inertia_p, mass, R_b_to_p,
                                com_offset_b, com_w, com_vel, Q_p_to_w, ω_p;
                                ω_kinematic, d_ω_p, d_com_w, d_com_vel)

Pure 6-DOF rigid-body derivative and body-frame output expressions (principal frame)
— the single math source both the monolith ([`rigid_body_eqs!`](@ref)) and the
scheduled body component assemble from. Given the world-frame load at / about the COM
(`force_w`, `moment_w`) and the principal state (`com_w`, `com_vel`, `Q_p_to_w`,
`ω_p` as length-3/4 `Num` vectors), returns a named tuple of expressions: the state
derivatives `d_com_w`/`d_com_vel`/`d_Q`/`d_ω`, the Euler angular accel `α_p`, the
quaternion rate `Q_p_vel`, the COM accel `com_acc`, the principal moment `moment_p`,
and the body-frame outputs `R_p_to_w`, `R_b_to_w`, `pos_w`, `vel_w`, `acc_w`, `ω_b`,
`α_b`, `Q_b_to_w`. The optional integration overrides (`ω_kinematic`, `d_ω_p`,
`d_com_w`, `d_com_vel`) reproduce `fix_sphere`/`STATIC`; left `nothing` the body
integrates freely.
"""
function rigid_body_pose_expressions(force_w, moment_w, inertia_p, mass, R_b_to_p,
                                     com_offset_b, com_w, com_vel, Q_p_to_w, ω_p;
                                     ω_kinematic = nothing, d_ω_p = nothing,
                                     d_com_w = nothing, d_com_vel = nothing)
    skew4(ω) = [
        0 -ω[1] -ω[2] -ω[3]
        ω[1] 0 ω[3] -ω[2]
        ω[2] -ω[3] 0 ω[1]
        ω[3] ω[2] -ω[1] 0
    ]
    quaternion = collect(Q_p_to_w)
    ω = collect(ω_p)
    com = collect(com_w)
    com_velocity = collect(com_vel)
    R_body_to_principal = collect(R_b_to_p)
    com_off = collect(com_offset_b)
    inertia = collect(inertia_p)

    ω_kin = ω_kinematic === nothing ? ω : collect(ω_kinematic)
    quat_norm_gain = 10.0
    quat_norm_error = 1 - sum(quaternion[j]^2 for j in 1:4)
    Q_p_vel = [0.5 * sum(skew4(ω_kin)[i, j] * quaternion[j] for j in 1:4) +
               quat_norm_gain * quat_norm_error * quaternion[i] for i in 1:4]

    R_p_to_w = quaternion_to_rotation_matrix(quaternion)
    moment_p = R_p_to_w' * collect(moment_w)
    α_p = [
        (moment_p[1] + (inertia[2] - inertia[3]) * ω[2] * ω[3]) / inertia[1],
        (moment_p[2] + (inertia[3] - inertia[1]) * ω[3] * ω[1]) / inertia[2],
        (moment_p[3] + (inertia[1] - inertia[2]) * ω[1] * ω[2]) / inertia[3],
    ]
    com_acc = collect(force_w) ./ mass

    R_b_to_w = R_p_to_w * R_body_to_principal
    arm_w = -(R_b_to_w * com_off)
    ω_w = R_p_to_w * ω
    pos_w = com .+ arm_w
    vel_w = com_velocity .+ (ω_w × arm_w)
    acc_w = com_acc .+ ((R_p_to_w * α_p) × arm_w) .+ (ω_w × (ω_w × arm_w))
    ω_b = R_body_to_principal' * ω
    α_b = R_body_to_principal' * α_p
    Q_b_to_w = rotation_matrix_to_quaternion(R_b_to_w)

    return (;
        d_com_w = d_com_w === nothing ? com_velocity : collect(d_com_w),
        d_com_vel = d_com_vel === nothing ? com_acc : collect(d_com_vel),
        d_Q = Q_p_vel,
        d_ω = d_ω_p === nothing ? α_p : collect(d_ω_p),
        α_p, Q_p_vel, com_acc, moment_p,
        R_p_to_w, R_b_to_w, pos_w, vel_w, acc_w, ω_b, α_b, Q_b_to_w)
end

"""
    pulley_split_eqs(pulley_len, pulley_vel, tension_in, pulley_mass, pulley_damp,
                     pulley_len_out)

The rope-split dynamics a pulley point owns on top of its particle motion:
`D(pulley_len)=pulley_vel`, `D(pulley_vel)=tension_in/pulley_mass − pulley_damp·vel`
(the aggregated `tension_in` being `spring[seg1] − spring[seg2]`), and
`pulley_len_out=pulley_len` exposed so the incident segments read it as their `l0`.
Shared by [`PulleyPoint`](@ref) and [`WingNodePulleyPoint`](@ref).
"""
function pulley_split_eqs(pulley_len, pulley_vel, tension_in, pulley_mass, pulley_damp,
                          pulley_len_out)
    return [
        D(pulley_len) ~ pulley_vel;
        D(pulley_vel) ~ tension_in / pulley_mass - pulley_damp * pulley_vel;
        pulley_len_out ~ pulley_len;
    ]
end

"""
    wing_frame_columns(zp1, zp2, yp1, yp2)

The particle-wing body→world rotation columns fitted from the four structural ref
points, matching `wing_eqs.jl`: `z = normalize(zp2−zp1)`,
`x = normalize(normalize(yp2−yp1) × z)`, `y = z × x`. Returns `(xaxis, yaxis,
zaxis)`, each a 3-vector (the columns of `R_b_to_w`).
"""
function wing_frame_columns(zp1, zp2, yp1, yp2)
    zaxis = smooth_normalize(zp2 .- zp1)
    xaxis = smooth_normalize(smooth_normalize(yp2 .- yp1) × zaxis)
    yaxis = zaxis × xaxis
    return xaxis, yaxis, zaxis
end

# ==================== component `System`s ==================== #
# One per component type; the `ScheduledBackend` compiles a kernel from each and
# instances it per component. Everything above is the physics they share.

"""
    point_variables()

The variables a point component may declare, as a named tuple: its world `pos`/`vel`
outputs, the summed `force_in`/`mass_in` of its incident segments, their `drag_in`
share, and the observed `total_drag`. Each component lists only the ones it uses.
"""
function point_variables()
    vars = @variables begin
        pos(t)[1:3], [output = true]
        vel(t)[1:3], [output = true]
        force_in(t)[1:3], [input = true]
        mass_in(t), [input = true]
        drag_in(t)[1:3], [input = true]
        total_drag(t)[1:3]
    end
    return (; pos = vars[1], vel = vars[2], force_in = vars[3], mass_in = vars[4],
            drag_in = vars[5], total_drag = vars[6])
end

"""
    total_drag_eq(s, params, idx, io)

Bind `total_drag` to the point's own aerodynamic drag plus the share its segments
deliver — the monolith's `total_drag`, which the state getter scatters into
`point.drag_force`.
"""
function total_drag_eq(s, params, idx, io)
    point = params.points[idx]
    wind = WindFactor(s.am, s.set.profile_law)
    height = collect(io.pos)[3]
    apparent = wind(height) .* ground_wind_vec(params) .- collect(io.vel)
    own = point_drag_force(apparent, calc_rho(s.am, max(0.0, height)),
                           point.drag_coeff, point.area)
    return collect(io.total_drag) .~ own .+ collect(io.drag_in)
end

"""
    Particle(s, params, idx; name)

A free point mass: integrates `pos`/`vel` under the force gathered from its incident
segments, its own drag, gravity and world-frame damping, through the shared
[`point_acceleration`](@ref). Its translational mass is `extra_mass` plus the
incident segments' half-masses (`mass_in`).
"""
function Particle(s, params, idx; name)
    io = point_variables()
    pars = point_particle_params(params, idx)
    eqs = [dynamic_point_dynamics(s, io.pos, io.vel, io.force_in,
                                  pars.extra_mass + io.mass_in, pars);
           total_drag_eq(s, params, idx, io)]
    vars = [io.pos, io.vel, io.force_in, io.mass_in, io.drag_in, io.total_drag]
    return System(eqs, t, vars, param_unknowns(params); name)
end

"""
    Anchor(s, params, idx; name)

A ground-anchored point: `pos` is pinned to `params.points[idx].pos_w` and `vel` is
zero, matching the STATIC branch of `point_eqs!`. It carries no force or mass input
— it does not move and its mass is not integrated — but still observes `total_drag`.
"""
function Anchor(s, params, idx; name)
    io = point_variables()
    eqs = [
        collect(io.pos) .~ collect(params.points[idx].pos_w)
        collect(io.vel) .~ zeros(3)
        total_drag_eq(s, params, idx, io)
    ]
    vars = [io.pos, io.vel, io.drag_in, io.total_drag]
    return System(eqs, t, vars, param_unknowns(params); name)
end

"""
    pulley_rope_mass(params, pulley_idx, segment_idx)

The rope mass `sum_len · ρ · π (d/2)²` driving a pulley's rope split, built in
equation from the pulley's `sum_len` and one of its segments' material, so it
follows the struct live.
"""
function pulley_rope_mass(params, pulley_idx, segment_idx)
    segment = params.segments[segment_idx]
    return params.pulleys[pulley_idx].sum_len * segment.density *
           π * (segment.diameter / 2)^2
end

"""
    PulleyParticle(s, params, idx, pulley_idx, segment_idx; name)

A [`Particle`](@ref) that also owns a pulley's rope split: `D(pulley_len) =
pulley_vel` and `D(pulley_vel) = tension_in / rope_mass − damping · pulley_vel`,
with `tension_in` the imbalance `spring[seg1] − spring[seg2]` its two segments
deliver. `pulley_len_out` exposes the split so those segments read it as their rest
length.
"""
function PulleyParticle(s, params, idx, pulley_idx, segment_idx; name)
    io = point_variables()
    extra = @variables begin
        tension_in(t), [input = true]
        pulley_len_out(t), [output = true]
        pulley_len(t)
        pulley_vel(t)
    end
    pars = point_particle_params(params, idx)
    damping = make_param(:pulley_damp, 5.0)
    eqs = [
        dynamic_point_dynamics(s, io.pos, io.vel, io.force_in,
                               pars.extra_mass + io.mass_in, pars)
        total_drag_eq(s, params, idx, io)
        pulley_split_eqs(extra[3], extra[4], extra[1],
                         pulley_rope_mass(params, pulley_idx, segment_idx),
                         damping, extra[2])
    ]
    vars = [io.pos, io.vel, io.force_in, io.mass_in, io.drag_in, io.total_drag,
            extra...]
    return System(eqs, t, vars, [param_unknowns(params); damping]; name)
end

"""
    WinchAnchor(s, params, winch, idx; name)

A reeling winch at the STATIC point `idx`. It owns the motor speed `winch_vel` and
one `tether_len_k` state per connected tether, each an output its tether's segments
read as their rest length. The motor law is [`winch_component`](@ref) reused
verbatim; it reads the tension gathered at the winch point (the vector sum of the
spring forces there, as `winch_eqs!` forms it), the mean tether length, the control
`set_value` and the `brake`, and returns the drum acceleration.
"""
function WinchAnchor(s, params, winch, idx; name)
    count = length(winch.tether_idxs)
    io = point_variables()
    extra = @variables begin
        tension_in(t)[1:3], [input = true]
        winch_vel(t)
        winch_force(t)
        winch_acc(t)
        winch_friction(t)
    end
    lengths = map(1:count) do k
        len_name = Symbol(:tether_len_, k)
        only(@variables $len_name(t), [output = true])
    end
    set_value = make_param(:set_value, 0.0)
    brake = params.winches[winch.idx].brake
    motor = winch_component(winch.model, s.sys_struct, winch.idx;
                            name = :motor, params)
    validate_winch_component(motor, winch)
    eqs = [
        collect(io.pos) .~ collect(params.points[idx].pos_w)
        collect(io.vel) .~ zeros(3)
        total_drag_eq(s, params, idx, io)
        extra[3] ~ smooth_norm(collect(extra[1]))
        motor.vel ~ extra[2]
        motor.len ~ sum(lengths) / count
        motor.force ~ extra[3]
        motor.set_value ~ set_value
        motor.brake ~ brake
        extra[5] ~ motor.friction
        extra[4] ~ ifelse(params.winches[winch.idx].speed_controlled == true,
                          0.0, motor.acc)
        D(extra[2]) ~ ifelse(brake > 0.5, 0.0, extra[4])
        [D(len) ~ ifelse(brake > 0.5, 0.0, extra[2]) for len in lengths]
    ]
    vars = [io.pos, io.vel, io.drag_in, io.total_drag, extra..., lengths...]
    return System(eqs, t, vars, [param_unknowns(params); set_value];
                  name, systems = [motor])
end

"""
    segment_variables()

The variables every segment kernel declares, as a named tuple: the two endpoints'
`pos`/`vel` inputs, the force on each endpoint (positive force-on-point sign), the
`half_mass` and `half_drag` each endpoint carries, and the `spring_force`/`len`/`l0`
diagnostics. Mass and drag are one output each because both endpoints receive the
same value; the wiring delivers it twice.
"""
function segment_variables()
    vars = @variables begin
        src_pos(t)[1:3], [input = true]
        src_vel(t)[1:3], [input = true]
        dst_pos(t)[1:3], [input = true]
        dst_vel(t)[1:3], [input = true]
        src_force(t)[1:3], [output = true]
        dst_force(t)[1:3], [output = true]
        half_mass(t), [output = true]
        half_drag(t)[1:3], [output = true]
        spring_force(t)
        len(t)
        l0(t)
    end
    return (; src_pos = vars[1], src_vel = vars[2], dst_pos = vars[3],
            dst_vel = vars[4], src_force = vars[5], dst_force = vars[6],
            half_mass = vars[7], half_drag = vars[8], spring_force = vars[9],
            len = vars[10], l0 = vars[11], all = vars)
end

"""
    segment_eqs(s, params, idx, io, rest_len; with_drag=true)

The equations every segment kernel shares: the shared [`segment_load_terms`](@ref)
evaluated at `rest_len`, bound to the [`segment_variables`](@ref) outputs and diagnostics.
Returns `(eqs, loads)`; `loads` carries the scalar and vector spring tension a
pulley or tether segment emits on top of these.
"""
function segment_eqs(s, params, idx, io, rest_len; with_drag = true)
    spring, wind = segment_spring_params(params, idx; with_drag)
    loads = segment_load_terms(s, collect(io.src_pos), collect(io.src_vel),
        collect(io.dst_pos), collect(io.dst_vel), spring.unit_stiffness,
        spring.unit_damping, spring.compression_frac, rest_len, spring.diameter,
        spring.density, spring.cd_tether, collect(wind);
        with_drag, nonlinear = spring.nonlinear)
    eqs = [
        collect(io.src_force) .~ loads.force_on_src
        collect(io.dst_force) .~ loads.force_on_dst
        io.half_mass ~ loads.half_mass
        collect(io.half_drag) .~ loads.half_drag
        io.spring_force ~ loads.spring
        io.len ~ loads.len
        io.l0 ~ rest_len
    ]
    return eqs, loads
end

"""
    SpringSegment(s, params, idx; name, with_drag=true)

A spring-damper segment at a fixed rest length (`params.segments[idx].l0`). With
`with_drag = false` it is the drag-free wing-structural link — a distinct compiled
type rather than a zeroed drag coefficient.
"""
function SpringSegment(s, params, idx; name, with_drag = true)
    io = segment_variables()
    eqs, _ = segment_eqs(s, params, idx, io, params.segments[idx].l0; with_drag)
    return System(eqs, t, io.all, param_unknowns(params); name)
end

"""
    PulleySegment(s, params, idx, pulley_idx; name)

One of a pulley's two segments. Its rest length comes from the pulley point's
`pulley_len_out`, read as `rest_len`: the first segment takes it directly, the
second `sum_len − rest_len`, selected by the assembly-set `pulley_side` (`±1`). It
emits `pulley_side · spring` as `tension`, so the pulley point aggregates
`spring[seg1] − spring[seg2]`.
"""
function PulleySegment(s, params, idx, pulley_idx; name)
    io = segment_variables()
    extra = @variables begin
        rest_len(t), [input = true]
        tension(t), [output = true]
    end
    side = make_param(:pulley_side, 1.0)
    sum_len = params.pulleys[pulley_idx].sum_len
    eqs, loads = segment_eqs(s, params, idx, io,
                             ifelse(side > 0.0, extra[1], sum_len - extra[1]))
    push!(eqs, extra[2] ~ side * loads.spring)
    return System(eqs, t, [io.all; extra], [param_unknowns(params); side]; name)
end

"""
    TetherSegment(s, params, idx; name)

A segment of a winched tether. Its rest length is the winch's `tether_len_k` output
divided by the assembly-set `segment_count`. Both endpoints' spring force are
emitted as `src_tension`/`dst_tension`; the assembly wires only the one at the winch
point, which is how the winch sees the force `winch_eqs!` sums there.
"""
function TetherSegment(s, params, idx; name)
    io = segment_variables()
    extra = @variables begin
        rest_len(t), [input = true]
        src_tension(t)[1:3], [output = true]
        dst_tension(t)[1:3], [output = true]
    end
    count = make_param(:segment_count, 1.0)
    eqs, loads = segment_eqs(s, params, idx, io, extra[1] / count)
    eqs = [eqs
           collect(extra[2]) .~ loads.spring_vec
           collect(extra[3]) .~ .-loads.spring_vec]
    return System(eqs, t, [io.all; extra], [param_unknowns(params); count]; name)
end

"""
    body_variables()

The variables a rigid-body component declares, as a named tuple: the aggregated
`force_in`/`moment_in` at and about its COM, the pose outputs a point riding the
body reads (`pos`, `vel`, `frame` — `R_b_to_w` column-major — `com`, `com_velocity`,
`omega_w`), and the observed `acc`/`orientation`/`omega_b` the state getter scatters
into the struct.
"""
function body_variables()
    vars = @variables begin
        pos(t)[1:3], [output = true]
        vel(t)[1:3], [output = true]
        frame(t)[1:9], [output = true]
        com(t)[1:3], [output = true]
        com_velocity(t)[1:3], [output = true]
        omega_w(t)[1:3], [output = true]
        force_in(t)[1:3], [input = true]
        moment_in(t)[1:3], [input = true]
        acc(t)[1:3]
        orientation(t)[1:4]
        omega_b(t)[1:3]
    end
    return (; pos = vars[1], vel = vars[2], frame = vars[3], com = vars[4],
            com_velocity = vars[5], omega_w = vars[6], force_in = vars[7],
            moment_in = vars[8], acc = vars[9], orientation = vars[10],
            omega_b = vars[11], all = vars)
end

"""
    body_integration(params, idx, com_w, com_vel, omega_p, alpha_p, com_acc,
                     orientation_p)

The four integration overrides `rigid_body_pose_expressions` takes, matching
`body_eqs!`: `fix_sphere` confines the body to a sphere about the world origin by
keeping only the radial part of its COM velocity and acceleration and dropping the
radial part of its spin. `alpha_p`/`com_acc` are the caller's torn variables, so the
overrides can name the accelerations they correct without a cycle. Angular damping
is folded in here, as the monolith does.
"""
function body_integration(params, idx, com_w, com_vel, omega_p, alpha_p, com_acc,
                          orientation_p)
    body = params.bodies[idx]
    sphere = body.fix_sphere
    spin = collect(omega_p)
    damped = collect(alpha_p) .- collect(body.damping) .* spin
    axis = collect(smooth_normalize(collect(com_w)))
    axis_p = orientation_p' * axis
    return (; ω_kinematic = ifelse.(sphere == true, remove_along(spin, axis_p), spin),
            d_ω_p = ifelse.(sphere == true, remove_along(damped, axis_p), damped),
            d_com_w = ifelse.(sphere == true, keep_along(collect(com_vel), axis),
                              collect(com_vel)),
            d_com_vel = ifelse.(sphere == true, keep_along(collect(com_acc), axis),
                                collect(com_acc)))
end

"""
    RigidBody(s, params, idx; name)

Free 6-DOF rigid body: integrates the 13-state principal pose (`com_w`, `com_vel`,
`Q_p_to_w`, `ω_p`) under gravity, the wrench gathered at `force_in`/`moment_in`, the
external wrench (`ext_force_w`/`ext_force_b`/`ext_moment_b`) and per-axis angular
`damping`. Shares its whole 6-DOF math with the monolith through
[`rigid_body_pose_expressions`](@ref); `fix_sphere` stays a parameter, as in
`body_eqs!`. A clamped body is [`StaticBody`](@ref) instead.
"""
function RigidBody(s, params, idx; name)
    io = body_variables()
    state = @variables com_w(t)[1:3] com_vel(t)[1:3] Q(t)[1:4] omega_p(t)[1:3]
    torn = @variables alpha_p(t)[1:3] com_acc(t)[1:3]
    com_w, com_vel, Q, omega_p = state
    alpha_p, com_acc = torn
    body = params.bodies[idx]
    orientation_p = quaternion_to_rotation_matrix(collect(Q))
    orientation = orientation_p * collect(body.R_b_to_p)
    gravity = Num[0, 0, -params.set.g_earth * body.mass]
    force_w = collect(io.force_in) .+ gravity .+ collect(body.ext_force_w) .+
        orientation * collect(body.ext_force_b)
    moment_w = collect(io.moment_in) .+ orientation * collect(body.ext_moment_b)
    ex = rigid_body_pose_expressions(force_w, moment_w, body.inertia_principal,
        body.mass, body.R_b_to_p, body.com_offset_b, com_w, com_vel, Q, omega_p;
        body_integration(params, idx, com_w, com_vel, omega_p, alpha_p,
                         com_acc, orientation_p)...)
    eqs = [
        collect(alpha_p) .~ ex.α_p
        collect(com_acc) .~ ex.com_acc
        [D(com_w[i]) ~ ex.d_com_w[i] for i in 1:3]
        [D(com_vel[i]) ~ ex.d_com_vel[i] for i in 1:3]
        [D(Q[i]) ~ ex.d_Q[i] for i in 1:4]
        [D(omega_p[i]) ~ ex.d_ω[i] for i in 1:3]
        collect(io.pos) .~ ex.pos_w
        collect(io.vel) .~ ex.vel_w
        collect(io.frame) .~ vec(ex.R_b_to_w)
        collect(io.com) .~ collect(com_w)
        collect(io.com_velocity) .~ collect(com_vel)
        collect(io.omega_w) .~ ex.R_b_to_w * ex.ω_b
        collect(io.acc) .~ ex.acc_w
        collect(io.orientation) .~ ex.Q_b_to_w
        collect(io.omega_b) .~ ex.ω_b
    ]
    return System(eqs, t, [io.all; state; torn], param_unknowns(params); name)
end

"""
    StaticBody(s, params, idx; name)

Clamped (`STATIC`) rigid body: no state at all, just the fixed pose
`params.bodies[idx]` already holds, with zero velocity and spin. The monolith
freezes such a body by holding its thirteen states at their initial values; here
there is nothing to hold, so the pose is emitted directly. Its `force_in`/`moment_in`
are declared, so joints and ride points may deliver to it, and ignored.
"""
function StaticBody(s, params, idx; name)
    io = body_variables()
    body = params.bodies[idx]
    orientation = quaternion_to_rotation_matrix(collect(body.Q_b_to_w))
    eqs = [
        collect(io.pos) .~ collect(body.pos_w)
        collect(io.vel) .~ zeros(3)
        collect(io.frame) .~ vec(orientation)
        collect(io.com) .~ collect(body.com_w)
        collect(io.com_velocity) .~ zeros(3)
        collect(io.omega_w) .~ zeros(3)
        collect(io.acc) .~ zeros(3)
        collect(io.orientation) .~ collect(body.Q_b_to_w)
        collect(io.omega_b) .~ zeros(3)
    ]
    return System(eqs, t, io.all, param_unknowns(params); name)
end

"""
    body_pose_variables()

The pose inputs a component riding a rigid body reads off that body's outputs:
the body origin `pose_pos`, its orientation `pose_frame` (`R_b_to_w`, column-major),
its `pose_com`/`pose_com_velocity` and its world spin `pose_omega`.
"""
function body_pose_variables()
    vars = @variables begin
        pose_pos(t)[1:3], [input = true]
        pose_frame(t)[1:9], [input = true]
        pose_com(t)[1:3], [input = true]
        pose_com_velocity(t)[1:3], [input = true]
        pose_omega(t)[1:3], [input = true]
    end
    return (; pose_pos = vars[1], pose_frame = vars[2], pose_com = vars[3],
            pose_com_velocity = vars[4], pose_omega = vars[5], all = vars)
end

"""
    RidePoint(s, params, idx; name)

The *kinematics* of a `BODY_STATIC` point: its body's pose in, and its world
`pos`/`vel`, its moment `arm` (`pos − com`) and its `height` out, exactly as
`body_ride_eqs` places it. It is deliberately only half of a ride point — the load
it feeds back to the body is [`RideWrench`](@ref) — because the two halves sit on
opposite sides of
the same dependency chain: the position must be known before the incident segments
can compute their forces, and the force is only known after. One component holding
both would be a cycle; two are a four-layer schedule.
"""
function RidePoint(s, params, idx; name)
    pose = body_pose_variables()
    io = @variables begin
        pos(t)[1:3], [output = true]
        vel(t)[1:3], [output = true]
        arm(t)[1:3], [output = true]
        height(t), [output = true]
    end
    orientation = reshape(collect(pose.pose_frame), 3, 3)
    anchor = collect(pose.pose_pos) .+ orientation * collect(params.points[idx].anchor_b)
    lever = anchor .- collect(pose.pose_com)
    eqs = [
        collect(io[1]) .~ anchor
        collect(io[3]) .~ lever
        collect(io[2]) .~ rigid_body_point_velocity(pose, lever)
        io[4] ~ anchor[3]
    ]
    return System(eqs, t, [pose.all; io], param_unknowns(params); name)
end

"""
    ride_wrench_variables()

The variables shared by every anchored point's *statics* half, as a named tuple:
its own `height` and `vel` from the kinematics half, the `force_in`/`mass_in`/
`drag_in` its segments deliver, and the observed `total_drag`. The moment arms
differ per anchor kind and are declared by the component itself.
"""
function ride_wrench_variables()
    vars = @variables begin
        height(t), [input = true]
        vel(t)[1:3], [input = true]
        force_in(t)[1:3], [input = true]
        mass_in(t), [input = true]
        drag_in(t)[1:3], [input = true]
        total_drag(t)[1:3]
    end
    return (; height = vars[1], vel = vars[2], force_in = vars[3],
            mass_in = vars[4], drag_in = vars[5], total_drag = vars[6],
            all = vars)
end

"""
    ride_load(s, params, idx, io; with_gravity) -> (load, drag)

The world load an anchored point delivers to whatever carries it: the force its
segments deliver, its own aerodynamic drag at its height, its gravity and its
external force. `with_gravity = false` is the point that rides its own wing body,
whose mass is already counted at that body's COM (`rides_own_wing` in
`point_eqs!`). The drag is returned separately because it is also half of the
point's observed `total_drag`.
"""
function ride_load(s, params, idx, io; with_gravity)
    point = params.points[idx]
    wind = WindFactor(s.am, s.set.profile_law)
    apparent = wind(io.height) .* ground_wind_vec(params) .- collect(io.vel)
    drag = point_drag_force(apparent, calc_rho(s.am, max(0.0, io.height)),
                            point.drag_coeff, point.area)
    mass = point.extra_mass + io.mass_in
    gravity = with_gravity ? Num[0, 0, -params.set.g_earth * mass] : zeros(Num, 3)
    load = collect(io.force_in) .+ drag .+ gravity .+ collect(point.ext_force_w)
    return load, drag
end

"""
    RideWrench(s, params, idx; name, with_gravity=true)

The *statics* of a `BODY_STATIC` point: everything that flows back to its body. It
takes the point's own `height`/`vel`/`arm` from [`RidePoint`](@ref) and the load
[`ride_load`](@ref) builds, and emits `force_out` and `moment_out = arm ×
force_out` into the body's `force_in`/`moment_in`. `total_drag` is observed, as for
any other point.
"""
function RideWrench(s, params, idx; name, with_gravity = true)
    io = ride_wrench_variables()
    vars = @variables begin
        arm(t)[1:3], [input = true]
        force_out(t)[1:3], [output = true]
        moment_out(t)[1:3], [output = true]
    end
    load, drag = ride_load(s, params, idx, io; with_gravity)
    eqs = [
        collect(vars[2]) .~ load
        collect(vars[3]) .~ collect(vars[1]) × load
        collect(io.total_drag) .~ drag .+ collect(io.drag_in)
    ]
    return System(eqs, t, [io.all; vars], param_unknowns(params); name)
end

"""
    HermiteRidePoint(s, params, idx; name)

The *kinematics* of a point riding a Timoshenko beam's deflected centerline: both
end bodies' poses in, and the point's world `pos`/`vel`, its two moment arms
(`arm_a`/`arm_b`, the offsets from each end body's COM) and its `height` out. The
beam-anchored counterpart of [`RidePoint`](@ref), and split for the same reason;
the placement itself is the shared [`beam_hermite_ride_expressions`](@ref) that also
sources the monolith's `beam_hermite_ride_eqs`. The velocity is evaluated at the
already-bound `pos` output, so the heavy element-frame subtree is built once.
"""
function HermiteRidePoint(s, params, idx; name)
    io = joint_pose_variables()
    vars = @variables begin
        pos(t)[1:3], [output = true]
        vel(t)[1:3], [output = true]
        arm_a(t)[1:3], [output = true]
        arm_b(t)[1:3], [output = true]
        height(t), [output = true]
    end
    sys_struct = params.reg.sys_struct
    joint = sys_struct.timoshenko_joints[sys_struct.points[idx].joint_idx]
    ex = beam_hermite_ride_expressions(joint, params, idx; joint_poses(io)...)
    eqs = [
        collect(vars[1]) .~ ex.pos_point
        collect(vars[2]) .~ ex.ride_velocity(vars[1])
        collect(vars[3]) .~ collect(vars[1]) .- collect(io.a_com)
        collect(vars[4]) .~ collect(vars[1]) .- collect(io.b_com)
        vars[5] ~ vars[1][3]
    ]
    return System(eqs, t, [io.all; vars], param_unknowns(params); name)
end

"""
    HermiteRideWrench(s, params, idx; name)

The *statics* of a point riding a Timoshenko beam: the load [`ride_load`](@ref)
builds, split along the element by the point's axial fraction `beam_frac` — `(1 −
s)` onto the first end body and `s` onto the second — with each half's moment about
that body's COM. The beam-anchored counterpart of [`RideWrench`](@ref).
"""
function HermiteRideWrench(s, params, idx; name)
    io = ride_wrench_variables()
    out = joint_wrench_variables()
    vars = @variables begin
        arm_a(t)[1:3], [input = true]
        arm_b(t)[1:3], [input = true]
    end
    load, drag = ride_load(s, params, idx, io; with_gravity = true)
    frac = params.points[idx].beam_frac
    force_a = (1 - frac) .* load
    force_b = frac .* load
    eqs = [
        collect(out.force_a) .~ force_a
        collect(out.moment_a) .~ collect(vars[1]) × force_a
        collect(out.force_b) .~ force_b
        collect(out.moment_b) .~ collect(vars[2]) × force_b
        collect(io.total_drag) .~ drag .+ collect(io.drag_in)
    ]
    return System(eqs, t, [io.all; out.all; vars], param_unknowns(params); name)
end

"""
    joint_pose_variables()

The two end bodies' poses as inputs (`a_pos`/`a_frame`/`a_com`/`a_com_velocity`/
`a_omega` and the `b_…` set), as a named tuple. Declared by everything that spans
an element: the joints themselves and a point riding one.
"""
function joint_pose_variables()
    vars = @variables begin
        a_pos(t)[1:3], [input = true]
        a_frame(t)[1:9], [input = true]
        a_com(t)[1:3], [input = true]
        a_com_velocity(t)[1:3], [input = true]
        a_omega(t)[1:3], [input = true]
        b_pos(t)[1:3], [input = true]
        b_frame(t)[1:9], [input = true]
        b_com(t)[1:3], [input = true]
        b_com_velocity(t)[1:3], [input = true]
        b_omega(t)[1:3], [input = true]
    end
    return (; a_pos = vars[1], a_frame = vars[2], a_com = vars[3],
            a_com_velocity = vars[4], a_omega = vars[5], b_pos = vars[6],
            b_frame = vars[7], b_com = vars[8], b_com_velocity = vars[9],
            b_omega = vars[10], all = vars)
end

"""
    joint_wrench_variables()

The wrench an element delivers to each of its two end bodies, as outputs.
"""
function joint_wrench_variables()
    vars = @variables begin
        force_a(t)[1:3], [output = true]
        moment_a(t)[1:3], [output = true]
        force_b(t)[1:3], [output = true]
        moment_b(t)[1:3], [output = true]
    end
    return (; force_a = vars[1], moment_a = vars[2], force_b = vars[3],
            moment_b = vars[4], all = vars)
end

"""
    joint_variables()

Everything a body-to-body joint component declares: both end bodies' poses in
([`joint_pose_variables`](@ref)) and the restoring wrench on each out
([`joint_wrench_variables`](@ref)).
"""
function joint_variables()
    poses = joint_pose_variables()
    wrench = joint_wrench_variables()
    return merge(poses, wrench, (; all = [poses.all; wrench.all]))
end

"""
    joint_poses(io)

The two end bodies' poses in the argument form the shared joint wrench builders
take, reading them off a [`joint_variables`](@ref) tuple. The world spin arrives
directly as `a_omega`/`b_omega`, so unlike the monolith — which carries `ω_b` and
rotates it — there is nothing to rotate here.
"""
function joint_poses(io)
    return (; pos_a = collect(io.a_pos), R_a = reshape(collect(io.a_frame), 3, 3),
            com_a = collect(io.a_com), com_vel_a = collect(io.a_com_velocity),
            omega_a_w = collect(io.a_omega),
            pos_b = collect(io.b_pos), R_b = reshape(collect(io.b_frame), 3, 3),
            com_b = collect(io.b_com), com_vel_b = collect(io.b_com_velocity),
            omega_b_w = collect(io.b_omega))
end

"""
    joint_wrench_eqs(io, ex)

Bind a joint component's four wrench outputs from a shared wrench builder's result
`ex`, together with the torn equations it needs.
"""
joint_wrench_eqs(io, ex) = [
    ex.tear_eqs
    collect(io.force_a) .~ ex.force_on_a
    collect(io.moment_a) .~ ex.moment_on_a
    collect(io.force_b) .~ ex.force_on_b
    collect(io.moment_b) .~ ex.moment_on_b
]

"""
    ElasticJointComponent(s, params, idx; name)

Lumped 6-DOF elastic joint between two bodies: reads both poses and emits the
restoring wrench on each, through the shared [`elastic_joint_wrench`](@ref) that
also sources the monolith's `joint_eqs!`.
"""
function ElasticJointComponent(s, params, idx; name)
    io = joint_variables()
    torn = @variables joint_force(t)[1:3] joint_torque(t)[1:3]
    joint = params.reg.sys_struct.elastic_joints[idx]
    ex = elastic_joint_wrench(joint, params; force_w = torn[1],
                              torque_w = torn[2], joint_poses(io)...)
    return System(joint_wrench_eqs(io, ex), t, [io.all; torn],
                  param_unknowns(params); name)
end

"""
    TimoshenkoJointComponent(s, params, idx; name)

Corotational Timoshenko beam element between two bodies: reads both poses and emits
the restoring wrench on each, through the shared
[`timoshenko_element_wrench`](@ref) that also sources the monolith's
`timoshenko_joint_eqs!`. The element frame, the two nodes' chord-relative rotations
and the element forces are torn, so the shared frame subtree is built once.
"""
function TimoshenkoJointComponent(s, params, idx; name)
    io = joint_variables()
    torn = @variables begin
        element_frame(t)[1:3, 1:3]
        theta_a(t)[1:3]
        theta_b(t)[1:3]
        element_force_a(t)[1:3]
        element_force_b(t)[1:3]
        element_moment_a(t)[1:3]
        element_moment_b(t)[1:3]
    end
    joint = params.reg.sys_struct.timoshenko_joints[idx]
    ex = timoshenko_element_wrench(joint, params; frame = torn[1],
        theta_a = torn[2], theta_b = torn[3], force_a = torn[4],
        force_b = torn[5], moment_a = torn[6], moment_b = torn[7],
        joint_poses(io)...)
    return System(joint_wrench_eqs(io, ex), t, [io.all; torn],
                  param_unknowns(params); name)
end

"""
    flap_delta_expression(twist_surface, R_main, R_flap)

The signed live deflection δ of a flapped twist surface: the angle between its two
flap bodies' reference chords about the world hinge axis, referenced to rest. The
axis, the reference chords and the rest angle are frozen rest geometry, so δ is a
function of the two bodies' orientations alone. Sources both the monolith's
`twist_surface_delta_eqs!` and the [`FlapDelta`](@ref) component.
"""
function flap_delta_expression(twist_surface, R_main, R_flap)
    main_w = collect(R_main) * collect(twist_surface.flap_chord_refs[1])
    flap_w = collect(R_flap) * collect(twist_surface.flap_chord_refs[2])
    normal = collect(R_main) * collect(twist_surface.flap_axis)
    projected_main = main_w .- (main_w ⋅ normal) .* normal
    projected_flap = flap_w .- (flap_w ⋅ normal) .* normal
    return atan(normal ⋅ (projected_main × projected_flap),
                projected_main ⋅ projected_flap) - twist_surface.flap_rest_delta
end

"""
    FlapDelta(s, params, idx; name)

The live flap deflection of a flapped `KINEMATIC` twist surface: its two flap
bodies' orientations in, δ out, through the shared
[`flap_delta_expression`](@ref). Deflection is all such a surface contributes —
it carries no twist DOF of its own — and the wing's aero component reads it per
panel.
"""
function FlapDelta(s, params, idx; name)
    vars = @variables begin
        main_frame(t)[1:9], [input = true]
        flap_frame(t)[1:9], [input = true]
        delta(t), [output = true]
    end
    twist_surface = params.reg.sys_struct.twist_surfaces[idx]
    eqs = [vars[3] ~ flap_delta_expression(twist_surface,
        reshape(collect(vars[1]), 3, 3), reshape(collect(vars[2]), 3, 3))]
    return System(eqs, t, vars, param_unknowns(params); name)
end

"""
    indexed_vector_variables(base, count; input=false) -> Vector

`count` three-component variables named `base_1 … base_count`, declared as inputs or
as outputs. A component whose I/O width follows the model — a wing's aero reads one
position per structural point — needs them named individually, because the wiring
addresses `base_k` as one vector.
"""
function indexed_vector_variables(base::Symbol, count::Int; input = false)
    return map(1:count) do k
        name = Symbol(base, :_, k)
        input ? only(@variables $name(t)[1:3], [input = true]) :
                only(@variables $name(t)[1:3], [output = true])
    end
end

"""
    flap_delta_inputs(wing, subsys) -> (; vars, eqs)

One `flap_delta_g` input per twist surface some panel of `wing` maps to, and the
equations binding the aero component's per-panel `delta` connector to them, exactly
as the `PARTICLE_DYNAMICS` branch of `aero_eqs!` wires it. Empty when the wing's
aero exposes no `delta` connector.
"""
function flap_delta_inputs(wing, subsys)
    surfaces = wing_flap_surfaces(wing)
    isempty(surfaces) && return (; vars = Any[], eqs = Equation[])
    inputs = Dict(surface => scalar_input(Symbol(:flap_delta_, surface))
                  for surface in surfaces)
    panel_map = wing.aero.panel_twist_surface
    eqs = [subsys.delta[i] ~ inputs[panel_map[i]] for i in eachindex(panel_map)]
    return (; vars = Any[inputs[surface] for surface in surfaces], eqs)
end

"""
    wing_flap_surfaces(wing) -> Vector{Int}

The twist surfaces `wing`'s aero panels deflect with, sorted and unique, or empty
when its aero mode has no per-panel flap coupling. One `flap_delta` input and one
[`FlapDelta`](@ref) instance exist per entry.
"""
function wing_flap_surfaces(wing)
    hasproperty(wing.aero, :panel_twist_surface) || return Int[]
    panel_map = wing.aero.panel_twist_surface
    (length(panel_map) == length(wing.vsm_aero.panels) &&
     any(!=(0), panel_map)) || return Int[]
    return sort!(unique(panel_map))
end

"""A scalar input variable named `name`."""
scalar_input(name::Symbol) = only(@variables $name(t), [input = true])

"""
    ParticleWingAero(s, params, idx; name)

The live aerodynamic force on a fitted wing's structural points: every point's world
`pos`/`vel` and the wing's own pose in, the world force on each point out. It
rebuilds each point's body-frame position, velocity, apparent wind and air density
exactly as the `PARTICLE_DYNAMICS` branch of `aero_eqs!` does, hands them to the
wing's aero component and rotates the per-point force it returns back into the
world. A frozen mode compiles down to one parameter read per point; a live mode
solves in place. The wing's lumped `aero_force_b`/`aero_moment_b` are observed, as
the monolith emits them.
"""
function ParticleWingAero(s, params, idx; name)
    sys_struct = params.reg.sys_struct
    wing = sys_struct.wings[idx]
    points = wing_points(sys_struct, wing)
    count = length(points)
    pose = @variables begin
        wing_pos(t)[1:3], [input = true]
        wing_frame(t)[1:9], [input = true]
        aero_force_b(t)[1:3]
        aero_moment_b(t)[1:3]
    end
    positions = indexed_vector_variables(:point_pos, count; input = true)
    velocities = indexed_vector_variables(:point_vel, count; input = true)
    forces = indexed_vector_variables(:point_force, count)
    subsys = aero_component(wing.aero, wing, sys_struct; name = :aero, params)
    validate_aero_component(subsys, wing)
    orientation = reshape(collect(pose[2]), 3, 3)
    origin = collect(pose[1])
    wind_factor = param_computed!(params.reg, :wind_factor, WindFactorReader())
    wind_gnd = ground_wind_vec(params)
    eqs = Equation[]
    for k in 1:count
        position = collect(positions[k])
        velocity = collect(velocities[k])
        apparent = wind_factor(position[3]) .* wind_gnd .- velocity
        append!(eqs, collect(subsys.point_pos[:, k]) .~
                     orientation' * (position .- origin))
        append!(eqs, collect(subsys.point_vel[:, k]) .~ orientation' * velocity)
        append!(eqs, collect(subsys.va[:, k]) .~ orientation' * apparent)
        push!(eqs, subsys.rho[k] ~ calc_rho(s.am, max(0.0, position[3])))
        append!(eqs, collect(forces[k]) .~
                     orientation * collect(subsys.point_force[:, k]))
    end
    append!(eqs, collect(pose[3]) .~
                 sum(collect(subsys.point_force[:, k]) for k in 1:count))
    append!(eqs, collect(pose[4]) .~
                 sum(collect(subsys.point_pos[:, k]) ×
                     collect(subsys.point_force[:, k]) for k in 1:count))
    flaps = flap_delta_inputs(wing, subsys)
    append!(eqs, flaps.eqs)
    vars = [pose; positions; velocities; forces; flaps.vars]
    return System(eqs, t, vars, param_unknowns(params); name, systems = [subsys])
end

"""
    rigid_body_point_velocity(pose, lever)

World velocity of a point rigidly attached to a body: `com_velocity + ω × lever`,
where `lever` is the point's offset from that body's COM.
"""
rigid_body_point_velocity(pose, lever) =
    collect(pose.pose_com_velocity) .+ (collect(pose.pose_omega) × collect(lever))

"""A scalar output variable named `name`."""
scalar_output(name::Symbol) = only(@variables $name(t), [output = true])

"""
    WingAero(s, params, idx; name)

The aerodynamic wrench of a `RIGID_DYNAMICS` wing: its body's pose in, the world
force and the moment about its COM out, plus one twist moment per twist surface it
carries. It rebuilds the wing's apparent wind and air density from that pose exactly
as `scalar_eqs!` and `aero_eqs!` do, hands them to the wing's aero component, and
transports the returned body-frame wrench to the COM as `create_sys` does. Like
[`ParticleWingAero`](@ref) it is a component of its own: a body's pose does not
depend on the force it receives, so aero after pose and before the body's derivative
is just another schedule layer.
"""
function WingAero(s, params, idx; name)
    sys_struct = params.reg.sys_struct
    wing = sys_struct.wings[idx]
    pose = body_pose_variables()
    io = @variables begin
        force_out(t)[1:3], [output = true]
        moment_out(t)[1:3], [output = true]
        aero_force_b(t)[1:3]
        aero_moment_b(t)[1:3]
        va_b(t)[1:3]
        wind_vel(t)[1:3]
    end
    surfaces = wing.twist_surface_idxs
    twists = [scalar_input(Symbol(:twist_angle_, surface)) for surface in surfaces]
    rates = [scalar_input(Symbol(:twist_vel_, surface)) for surface in surfaces]
    moments = [scalar_output(Symbol(:twist_moment_, surface)) for surface in surfaces]
    subsys = aero_component(wing.aero, wing, sys_struct; name = :aero, params)
    validate_aero_component(subsys, wing)
    orientation = reshape(collect(pose.pose_frame), 3, 3)
    origin = collect(pose.pose_pos)
    velocity = rigid_body_point_velocity(pose, origin .- collect(pose.pose_com))
    wind_factor = param_computed!(params.reg, :wind_factor, WindFactorReader())
    eqs = [
        collect(io[6]) .~ wind_factor(origin[3]) .* ground_wind_vec(params)
        collect(io[5]) .~ orientation' * (collect(io[6]) .- velocity .+
                                          collect(params.wings[idx].wind_disturb))
        collect(subsys.va) .~ collect(io[5])
        subsys.rho ~ calc_rho(s.am, origin[3])
        vec(collect(subsys.R_b_w)) .~ collect(pose.pose_frame)
        collect(subsys.omega) .~ orientation' * collect(pose.pose_omega)
        collect(io[3]) .~ collect(subsys.force)
        collect(io[4]) .~ collect(subsys.moment)
        collect(io[1]) .~ orientation * collect(io[3])
        collect(io[2]) .~ orientation * (collect(io[4]) .+
            (collect(io[3]) × collect(params.bodies[idx].com_offset_b)))
    ]
    if !isempty(surfaces)
        append!(eqs, collect(subsys.twist) .~ twists)
        append!(eqs, collect(subsys.twist_vel) .~ rates)
        append!(eqs, [moments[j] ~ subsys.twist_moment[j]
                      for j in eachindex(surfaces)])
    end
    vars = [pose.all; io; twists; rates; moments]
    return System(eqs, t, vars, param_unknowns(params); name, systems = [subsys])
end

"""
    TwistSurfaceDOF(s, params, idx; name)

The added twist degree of freedom of a `DYNAMIC` twist surface: a thin plate hinged
at its leading edge, driven by the aerodynamic moment its wing's aero returns and
the bridle couple its points deliver, restrained by the surface's own stiffness and
damping. Its inertia `⅓·m·L²` takes the mass from those same points as an input, so
the component reads only its own surface's parameters. The monolith's `fix_wing`
freeze is not carried over: it is a parameter nothing ever sets.
"""
function TwistSurfaceDOF(s, params, idx; name)
    vars = @variables begin
        aero_moment_in(t), [input = true]
        node_moment_in(t), [input = true]
        node_mass_in(t), [input = true]
        twist_angle(t), [output = true]
        twist_vel(t), [output = true]
    end
    state = @variables free_twist_angle(t) twist_omega(t)
    surface = params.twist_surfaces[idx]
    inertia = 1 / 3 * vars[3] * smooth_norm(collect(surface.chord))^2
    angle = clamp(state[1], -deg2rad(90), deg2rad(90))
    eqs = [
        vars[4] ~ angle
        vars[5] ~ state[2]
        D(state[1]) ~ state[2]
        D(state[2]) ~ (vars[1] + vars[2]) / inertia - surface.damping * state[2] -
                      surface.stiffness * angle / inertia
    ]
    return System(eqs, t, [vars; state], param_unknowns(params); name)
end

"""
    PrescribedTwist(s, params, idx; name)

A `STATIC` twist surface's prescribed section twist: no state and no inputs, just
the `twist` its parameters hold. It exists so a node reading a twist angle reads one
whether its surface twists dynamically ([`TwistSurfaceDOF`](@ref)) or not.
"""
function PrescribedTwist(s, params, idx; name)
    vars = @variables begin
        twist_angle(t), [output = true]
        twist_vel(t), [output = true]
    end
    eqs = [vars[1] ~ params.twist_surfaces[idx].twist
           vars[2] ~ 0]
    return System(eqs, t, vars, param_unknowns(params); name)
end

"""
    twist_deformed_offset(params, idx, surface_idx, angle)

The body-frame offset of a rigid wing's structural node, rotated about its surface's
leading edge by the section twist `angle` — the placement `point_eqs!` gives a wing
node whose surface twists. Without a surface (`surface_idx == 0`) it is the node's
own `pos_b`.
"""
function twist_deformed_offset(params, idx, surface_idx, angle)
    surface_idx == 0 && return collect(params.points[idx].pos_b)
    surface = params.twist_surfaces[surface_idx]
    leading_edge = collect(surface.le_pos)
    chord = collect(params.points[idx].pos_b) .- leading_edge
    normal = chord × collect(surface.y_airf)
    return leading_edge .+ cos(angle) .* chord .- sin(angle) .* normal
end

"""
    TwistNodePoint(s, params, idx; name, surface_idx=0)

The *kinematics* of a structural node on a `RIGID_DYNAMICS` wing: its body's pose
and its surface's twist in, its world `pos`, its moment `arm` about the body COM and
its `height` out. Such a node is placed from the body's COM by a twist-deformed
body-frame offset ([`twist_deformed_offset`](@ref)) rather than by a fixed anchor,
and `point_eqs!` gives it no velocity of its own, which is why it is not a
[`RidePoint`](@ref).
"""
function TwistNodePoint(s, params, idx; name, surface_idx = 0)
    pose = body_pose_variables()
    io = @variables begin
        pos(t)[1:3], [output = true]
        vel(t)[1:3], [output = true]
        arm(t)[1:3], [output = true]
        height(t), [output = true]
    end
    twist = scalar_input(:twist_angle)
    orientation = reshape(collect(pose.pose_frame), 3, 3)
    lever = orientation * twist_deformed_offset(params, idx, surface_idx, twist)
    anchor = collect(pose.pose_com) .+ lever
    eqs = [
        collect(io[1]) .~ anchor
        collect(io[2]) .~ zeros(3)
        collect(io[3]) .~ lever
        io[4] ~ anchor[3]
    ]
    return System(eqs, t, [pose.all; io; twist], param_unknowns(params); name)
end

"""
    TwistNodeWrench(s, params, idx; name, surface_idx=0, gated=false)

The *statics* of a structural node on a `RIGID_DYNAMICS` wing: the load its
segments, its own drag and its external force deliver — no gravity, because the wing
body already carries the node's mass — the moment that load makes about the body COM
and, when the node belongs to a twist surface, the bridle couple it exerts on that
surface's hinge and the mass it lends to the surface's inertia. `gated` is the wing's
`group_points_moment = false`, which drops an in-surface node's moment on the body.
"""
function TwistNodeWrench(s, params, idx; name, surface_idx = 0, gated = false)
    io = ride_wrench_variables()
    vars = @variables begin
        arm(t)[1:3], [input = true]
        frame(t)[1:9], [input = true]
        force_out(t)[1:3], [output = true]
        moment_out(t)[1:3], [output = true]
    end
    twist = scalar_input(:twist_angle)
    point = params.points[idx]
    wind = WindFactor(s.am, s.set.profile_law)
    apparent = wind(io.height) .* ground_wind_vec(params) .- collect(io.vel)
    drag = point_drag_force(apparent, calc_rho(s.am, max(0.0, io.height)),
                            point.drag_coeff, point.area)
    load = collect(io.force_in) .+ drag .+ collect(point.ext_force_w)
    eqs = [
        collect(vars[3]) .~ load
        collect(vars[4]) .~ (gated ? zeros(Num, 3) : collect(vars[1]) × load)
        collect(io.total_drag) .~ drag .+ collect(io.drag_in)
    ]
    extra = Any[twist]
    if surface_idx > 0
        surface = params.twist_surfaces[surface_idx]
        node_moment = scalar_output(:node_moment)
        node_mass = scalar_output(:node_mass)
        chord = collect(surface.chord)
        axis = collect(smooth_normalize(chord))
        section_normal = sin(twist) .* axis .+
                         cos(twist) .* (axis × collect(surface.y_airf))
        direction = reshape(collect(vars[2]), 3, 3) * (-1 .* section_normal)
        offset = collect(point.pos_b) .- (collect(surface.le_pos) .+
                                          surface.moment_frac .* chord)
        eqs = [eqs
               node_moment ~ (offset ⋅ axis) * (load ⋅ direction)
               node_mass ~ point.extra_mass]
        append!(extra, Any[node_moment, node_mass])
    end
    return System(eqs, t, [io.all; vars; extra], param_unknowns(params); name)
end

"""
    wing_frame_variables()

The reference-point inputs a fitted (`KINEMATIC`) wing body reads: the two z-axis
and two y-axis structural reference positions its frame is built from, plus the
origin's position and velocity. Each is a weighted blend of real points, which the
wiring supplies — see [`Wiring`](@ref).
"""
function wing_frame_variables()
    vars = @variables begin
        z1_pos(t)[1:3], [input = true]
        z2_pos(t)[1:3], [input = true]
        y1_pos(t)[1:3], [input = true]
        y2_pos(t)[1:3], [input = true]
        origin_pos(t)[1:3], [input = true]
        origin_vel(t)[1:3], [input = true]
    end
    return (; z1_pos = vars[1], z2_pos = vars[2], y1_pos = vars[3],
            y2_pos = vars[4], origin_pos = vars[5], origin_vel = vars[6],
            all = vars)
end

"""
    KinematicBody(s, params, idx; name)

A `KINEMATIC` (particle-wing) body: no state at all. Its orientation is *fitted*
from four structural reference points through the shared
[`wing_frame_columns`](@ref), and its origin pose is read from its origin reference
point, exactly as the `KINEMATIC` branch of `wing_eqs!` does. The principal frame
is aliased to the body frame and the spin is zero, because such a wing has no rigid
rotation of its own — the particles carry the motion.
"""
function KinematicBody(s, params, idx; name)
    frame = wing_frame_variables()
    io = body_variables()
    columns = wing_frame_columns(collect(frame.z1_pos), collect(frame.z2_pos),
                                collect(frame.y1_pos), collect(frame.y2_pos))
    orientation = [columns[1][1] columns[2][1] columns[3][1];
                   columns[1][2] columns[2][2] columns[3][2];
                   columns[1][3] columns[2][3] columns[3][3]]
    eqs = [
        collect(io.pos) .~ collect(frame.origin_pos)
        collect(io.vel) .~ collect(frame.origin_vel)
        collect(io.frame) .~ vec(orientation)
        collect(io.com) .~ collect(frame.origin_pos)
        collect(io.com_velocity) .~ collect(frame.origin_vel)
        collect(io.omega_w) .~ zeros(3)
        collect(io.acc) .~ zeros(3)
        collect(io.orientation) .~ rotation_matrix_to_quaternion(orientation)
        collect(io.omega_b) .~ zeros(3)
    ]
    return System(eqs, t, [io.all; frame.all], param_unknowns(params); name)
end

"""
    body_frame_damp_accel(vel, body_damp, orientation, wing_vel)

Body-frame damping acceleration `R·(coeff ⊙ (Rᵀ·(vel − wing_vel)))`, the term
`point_damping_accel` adds for a point that damps against its wing's frame rather
than the world.
"""
body_frame_damp_accel(vel, body_damp, orientation, wing_vel) =
    orientation * (collect(body_damp) .* (orientation' * (collect(vel) .- wing_vel)))

"""
    WingNodePoint(s, params, idx; name, with_damping=true)

A particle belonging to a fitted wing: the shared [`Particle`](@ref) motion plus the
body-frame damping `point_eqs!` adds for a point that damps against its wing's frame
rather than the world, which reads the `wing_frame`/`wing_velocity` the fitted body
supplies. Its aerodynamic force arrives at `force_in` from
[`ParticleWingAero`](@ref), like any other load. `with_damping = false` drops the
damping term whose coefficient the struct leaves unset, so no zero-valued parameter
is generated, matching the monolith.
"""
function WingNodePoint(s, params, idx; name, with_damping = true)
    io = point_variables()
    extra = @variables begin
        wing_frame(t)[1:9], [input = true]
        wing_velocity(t)[1:3], [input = true]
    end
    orientation = reshape(collect(extra[1]), 3, 3)
    point = params.points[idx]
    pars = point_particle_params(params, idx)
    mass = pars.extra_mass + io.mass_in
    accel = point_acceleration(s, collect(io.pos), collect(io.vel),
        collect(io.force_in), mass, pars.drag_coeff, pars.area,
        collect(pars.world_damping), collect(pars.wind_gnd))
    with_damping && (accel = accel .- body_frame_damp_accel(io.vel,
        point.body_frame_damping, orientation, collect(extra[2])))
    velocity, acceleration = confined_derivatives(io.pos, io.vel, accel, pars)
    eqs = [
        D.(collect(io.pos)) .~ velocity
        D.(collect(io.vel)) .~ acceleration
        total_drag_eq(s, params, idx, io)
    ]
    vars = [io.pos, io.vel, io.force_in, io.mass_in, io.drag_in, io.total_drag,
            extra...]
    return System(eqs, t, vars, param_unknowns(params); name)
end
