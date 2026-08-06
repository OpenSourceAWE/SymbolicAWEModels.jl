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
    vertex_pose_io()

The extra wide-interface variables a non-body (point) vertex needs so it shares the
body-containing network's uniform I/O width (§8.5): the pose outputs
`pose_R`(9)/`pose_com`/`pose_com_vel`/`pose_omega` (zero for a point) and the unused
`moment_in` input. Returns `(vars, pose_R, pose_com, pose_com_vel, pose_omega,
moment_in)`.
"""
function vertex_pose_io()
    vars = @variables begin
        pose_R(t)[1:9]
        pose_com(t)[1:3]
        pose_com_vel(t)[1:3]
        pose_omega(t)[1:3]
        moment_in(t)[1:3], [input = true]
    end
    return vars, pose_R, pose_com, pose_com_vel, pose_omega, moment_in
end

"""
    finish_vertex(vars, eqs, params; name, wide)

Assemble a point vertex `System` from its narrow `vars`/`eqs`. When `wide`, append the
zero-valued pose outputs and the unused `moment_in` input ([`vertex_pose_io`](@ref)) so
the vertex matches the wide superset a body-containing network uses; otherwise build the
narrow point vertex unchanged (no regression on point-only models).
"""
function finish_vertex(vars, eqs, params; name, wide)
    if wide
        pvars, pose_R, pose_com, pose_com_vel, pose_omega, _ = vertex_pose_io()
        eqs = [eqs;
               collect(pose_R) .~ 0; collect(pose_com) .~ 0;
               collect(pose_com_vel) .~ 0; collect(pose_omega) .~ 0]
        vars = [vars; pvars]
    end
    return System(eqs, t, vars, param_unknowns(params); name)
end

"""
    segment_io()

The array-valued I/O variables of a segment component, the uniform edge interface
both backends share: the two endpoints' `pos`/`vel`/`pulley_len` read as inputs and
the endpoint loads written as outputs (`src_force`/`src_mass`/`src_tension`,
`dst_…`, positive force-on-point sign). `pulley_len` inputs and `tension` outputs
are used only by pulley/winch assembly (a plain segment ignores the inputs and emits
zero tension), but every edge declares them so the graph has one edge width. Returns
`(vars, src_pos, src_vel, src_pulley_len, dst_pos, dst_vel, dst_pulley_len,
src_force, src_mass, src_tension, dst_force, dst_mass, dst_tension)`.
"""
function segment_io()
    vars = @variables begin
        src_pos(t)[1:3], [input = true]
        src_vel(t)[1:3], [input = true]
        src_pulley_len(t), [input = true]
        dst_pos(t)[1:3], [input = true]
        dst_vel(t)[1:3], [input = true]
        dst_pulley_len(t), [input = true]
        src_force(t)[1:3], [output = true]
        src_mass(t), [output = true]
        src_tension(t), [output = true]
        dst_force(t)[1:3], [output = true]
        dst_mass(t), [output = true]
        dst_tension(t), [output = true]
    end
    return vars, src_pos, src_vel, src_pulley_len, dst_pos, dst_vel, dst_pulley_len,
           src_force, src_mass, src_tension, dst_force, dst_mass, dst_tension
end

"""
    segment_io_wide()

The **wide** edge I/O a body-containing network shares (§8.5): the narrow
[`segment_io`](@ref) tuple (positions 1–13 unchanged, so [`segment_loads`](@ref)/
[`endpoint_load_eqs`](@ref) read it directly), extended with each endpoint's pose
inputs (`src_pose_R`(9)/`src_pose_com`/`src_pose_com_vel`/`src_pose_omega`, `dst_…`) a
wrench edge reads off a body vertex, and the `src_moment`/`dst_moment` outputs (zero for
a plain segment). The first tuple element holds *all* variables (base + extras) for the
`System`. Returns `(io, extras)` where `io` is the 13-tuple and `extras` a NamedTuple of
the added handles.
"""
function segment_io_wide()
    base = segment_io()
    extravars = @variables begin
        src_pose_R(t)[1:9], [input = true]
        src_pose_com(t)[1:3], [input = true]
        src_pose_com_vel(t)[1:3], [input = true]
        src_pose_omega(t)[1:3], [input = true]
        dst_pose_R(t)[1:9], [input = true]
        dst_pose_com(t)[1:3], [input = true]
        dst_pose_com_vel(t)[1:3], [input = true]
        dst_pose_omega(t)[1:3], [input = true]
        src_moment(t)[1:3], [output = true]
        dst_moment(t)[1:3], [output = true]
    end
    allvars = [base[1]; extravars]
    io = (allvars, base[2:end]...)
    extras = (; src_pose_R = extravars[1], src_pose_com = extravars[2],
              src_pose_com_vel = extravars[3], src_pose_omega = extravars[4],
              dst_pose_R = extravars[5], dst_pose_com = extravars[6],
              dst_pose_com_vel = extravars[7], dst_pose_omega = extravars[8],
              src_moment = extravars[9], dst_moment = extravars[10])
    return io, extras
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
            area = point.area, world_damping = point.world_frame_damping, wind_gnd)
end

"""
    dynamic_point_dynamics(s, pos, vel, force, mass, pars)

Shared body of the DYNAMIC point/pulley vertices: `D(pos)=vel` and
`D(vel)=point_acceleration(...)` from the shared kernel, reading the vertex's
drag/damping/wind parameters `pars` (a [`point_particle_params`](@ref) named tuple).
Both backends build their particle vertices on this so the integrator lives once.
"""
function dynamic_point_dynamics(s, pos, vel, force, mass, pars)
    accel = point_acceleration(s, collect(pos), collect(vel), collect(force),
        mass, pars.drag_coeff, pars.area, collect(pars.world_damping),
        collect(pars.wind_gnd))
    return [D.(collect(pos)) .~ collect(vel);
            D.(collect(vel)) .~ accel]
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
    segment_spring_params(params, idx; with_drag=true)

The spring-damper parameters read from `params.segments[idx]` (stiffness, damping,
compression fraction, diameter, density as the segment's own struct fields), plus
the global tether drag `cd_tether` (`params.set.cd_tether`) and the ground wind
`wind_gnd` ([`ground_wind_vec`](@ref)). With `with_drag=false` (the
[`wing_structural_segment`](@ref) edge) `cd_tether` is a literal `0` and unused.
Returns `(spring_named_tuple, wind_gnd)`; each read registers the parameter on
`params`.
"""
function segment_spring_params(params, idx; with_drag = true)
    seg = params.segments[idx]
    cd_tether = with_drag ? params.set.cd_tether : 0.0
    wind_gnd = ground_wind_vec(params)
    spring = (; unit_stiffness = seg.unit_stiffness, unit_damping = seg.unit_damping,
              compression_frac = seg.compression_frac, diameter = seg.diameter,
              density = seg.density, cd_tether)
    return spring, wind_gnd
end

"""
    segment_loads(s, io, l0, spring, wind; with_drag=true)

Compute the positive endpoint forces, half-mass and scalar spring tension from the
shared [`segment_endpoint_loads`](@ref), reading the array-valued endpoint states out
of the [`segment_io`](@ref) tuple `io`. Returns `(force_on_src, force_on_dst,
half_mass, spring)`.
"""
function segment_loads(s, io, l0, spring, wind; with_drag = true)
    (_, src_pos, src_vel, _, dst_pos, dst_vel, _) = io
    return segment_endpoint_loads(
        s, collect(src_pos), collect(src_vel), collect(dst_pos), collect(dst_vel),
        spring.unit_stiffness, spring.unit_damping, spring.compression_frac, l0,
        spring.diameter, spring.density, spring.cd_tether, collect(wind); with_drag)
end

"""
    endpoint_load_eqs(io, force_on_src, force_on_dst, half_mass,
                      src_tension_val, dst_tension_val)

Bind the edge output variables (`src_force`/`src_mass`/`src_tension`, `dst_…`) of the
[`segment_io`](@ref) tuple `io` from the computed endpoint loads and the role-signed
tensions each endpoint reads. A plain segment passes `0` for both tensions.
"""
function endpoint_load_eqs(io, force_on_src, force_on_dst, half_mass,
                           src_tension_val, dst_tension_val)
    (_, _, _, _, _, _, _, src_force, src_mass, src_tension,
     dst_force, dst_mass, dst_tension) = io
    return [
        collect(src_force) .~ force_on_src;
        src_mass ~ half_mass; src_tension ~ src_tension_val;
        collect(dst_force) .~ force_on_dst;
        dst_mass ~ half_mass; dst_tension ~ dst_tension_val;
    ]
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
`timoshenko_joint_eqs!` loop body and the network joint edge), so the beam-element
physics lives in one place. Given the two nodes' world poses (`pos`, `R_b_to_w`, `com`,
`com_vel`, world spin `omega_w`) and the joint's rest geometry/rigidities, it builds the
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
`joint_eqs!` loop body and the network elastic-joint edge). From the relative pose of the
two anchors (in body A's frame) it builds the per-DOF restoring force/torque (axial,
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
deflection (+ a frame-carried `beam_offset_b`) and returns `(pos_point, vel_point, sfrac)`
— the ride position, its rigid-blend velocity, and the axial fraction that splits any
force applied at the point onto the two end bodies (`(1−sfrac)` to A, `sfrac` to B).
"""
function beam_hermite_ride_expressions(joint, params, point_idx;
        pos_a, R_a, com_a, com_vel_a, omega_a_w,
        pos_b, R_b, com_b, com_vel_b, omega_b_w)
    jp = params.timoshenko_joints[joint.idx]
    Ra = collect(R_a)
    Rb = collect(R_b)
    x_a = collect(pos_a) .+ Ra * collect(jp.anchor_a_b)
    x_b = collect(pos_b) .+ Rb * collect(jp.anchor_b_b)
    e1, e2, e3, beam_len = timoshenko_element_frame(x_a, x_b, Ra)
    element_frame = [e1[1] e2[1] e3[1];
                     e1[2] e2[2] e3[2];
                     e1[3] e2[3] e3[3]]
    Da = (element_frame' * Ra) * collect(jp.R_a_rel0)'
    Db = (element_frame' * Rb) * collect(jp.R_b_rel0)'
    θ_a = [0.5 * (Da[3, 2] - Da[2, 3]), 0.5 * (Da[1, 3] - Da[3, 1]),
           0.5 * (Da[2, 1] - Da[1, 2])]
    θ_b = [0.5 * (Db[3, 2] - Db[2, 3]), 0.5 * (Db[1, 3] - Db[3, 1]),
           0.5 * (Db[2, 1] - Db[1, 2])]
    sfrac = params.points[point_idx].beam_frac
    N2 = beam_len * (sfrac - 2sfrac^2 + sfrac^3)
    N4 = beam_len * (-sfrac^2 + sfrac^3)
    v_defl = N2 * θ_a[3] + N4 * θ_b[3]
    w_defl = -(N2 * θ_a[2] + N4 * θ_b[2])
    x_center = x_a .+ (sfrac * beam_len) .* e1 .+ v_defl .* e2 .+ w_defl .* e3
    pos_point = x_center .+ element_frame * collect(params.points[point_idx].beam_offset_b)
    ω_a_w = collect(omega_a_w)
    ω_b_w = collect(omega_b_w)
    vel_a = collect(com_vel_a) .+ (ω_a_w × (pos_point .- collect(com_a)))
    vel_b = collect(com_vel_b) .+ (ω_b_w × (pos_point .- collect(com_b)))
    vel_point = (1 - sfrac) .* vel_a .+ sfrac .* vel_b
    return (; pos_point, vel_point, sfrac)
end

"""
    rigid_body_pose_expressions(force_w, moment_w, inertia_p, mass, R_b_to_p,
                                com_offset_b, com_w, com_vel, Q_p_to_w, ω_p;
                                ω_kinematic, d_ω_p, d_com_w, d_com_vel)

Pure 6-DOF rigid-body derivative and body-frame output expressions (principal frame)
— the single math source both the monolith ([`rigid_body_eqs!`](@ref)) and the
network body vertex assemble from. Given the world-frame load at / about the COM
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
    body_io()

Declare the **wide** vertex I/O a body-containing network shares (§8.5): the uniform
output superset `pos`/`vel` (body origin), `pulley_len_out`, and the pose
`pose_R`(9)/`pose_com`/`pose_com_vel`/`pose_omega` an incident wrench edge reads to place
any ride anchor and transport its moment; the input superset is the aggregated wrench
`force_in`/`moment_in` plus the `mass_in`/`tension_in` a point vertex uses (declared but
unused here, so every vertex shares one input width). Returns `(vars, pos, vel,
pulley_len_out, pose_R, pose_com, pose_com_vel, pose_omega, force_in, moment_in)`.
"""
function body_io()
    vars = @variables begin
        pos(t)[1:3]
        vel(t)[1:3]
        pulley_len_out(t)
        pose_R(t)[1:9]
        pose_com(t)[1:3]
        pose_com_vel(t)[1:3]
        pose_omega(t)[1:3]
        force_in(t)[1:3], [input = true]
        mass_in(t), [input = true]
        tension_in(t), [input = true]
        moment_in(t)[1:3], [input = true]
    end
    return vars, pos, vel, pulley_len_out, pose_R, pose_com, pose_com_vel,
        pose_omega, force_in, moment_in
end

"""
    BodyVertex(s, params, idx; name)

Free rigid-body vertex: integrates the 13-state principal pose (`com_w`, `com_vel`,
`Q_p_to_w`, `ω_p`) under gravity, the aggregated wrench at `force_in`/`moment_in`,
the external wrench (`ext_force_w`/`ext_force_b`/`ext_moment_b`), and per-axis
angular `damping` — all read from `params.bodies[idx]`. Shares the 6-DOF math with
the monolith through [`rigid_body_pose_expressions`](@ref); emits the wide pose via
[`body_io`](@ref) (`pulley_len_out = 0`; `mass_in`/`tension_in` unused). `fix_sphere`
confinement is not built here. A clamped (`STATIC`) body uses [`StaticBody`](@ref) instead.
"""
function BodyVertex(s, params, idx; name)
    vars, pos, vel, pulley_len_out, pose_R, pose_com, pose_com_vel, pose_omega,
        force_in, moment_in = body_io()
    state = @variables com_w(t)[1:3] com_vel(t)[1:3] Q(t)[1:4] omega_p(t)[1:3]
    body = params.bodies[idx]
    R_b_to_w = quaternion_to_rotation_matrix(collect(Q)) * collect(body.R_b_to_p)
    gravity_w = Num[0, 0, -params.set.g_earth * body.mass]
    force_w = collect(force_in) .+ gravity_w .+ collect(body.ext_force_w) .+
        R_b_to_w * collect(body.ext_force_b)
    moment_w = collect(moment_in) .+ R_b_to_w * collect(body.ext_moment_b)
    ex = rigid_body_pose_expressions(force_w, moment_w, body.inertia_principal,
        body.mass, body.R_b_to_p, body.com_offset_b, com_w, com_vel, Q, omega_p)
    damping = collect(body.damping)
    eqs = [
        [D(com_w[i]) ~ ex.d_com_w[i] for i in 1:3]
        [D(com_vel[i]) ~ ex.d_com_vel[i] for i in 1:3]
        [D(Q[i]) ~ ex.d_Q[i] for i in 1:4]
        [D(omega_p[i]) ~ ex.α_p[i] - damping[i] * omega_p[i] for i in 1:3]
        pos ~ ex.pos_w
        vel ~ ex.vel_w
        pulley_len_out ~ 0.0
        [pose_R[k] ~ vec(ex.R_b_to_w)[k] for k in 1:9]
        pose_com ~ com_w
        pose_com_vel ~ com_vel
        pose_omega ~ ex.R_b_to_w * ex.ω_b
    ]
    return System(eqs, t, [vars; state], param_unknowns(params); name)
end

"""
    StaticBody(s, params, idx; name)

Clamped (`STATIC`) rigid-body vertex: no dynamic state, a fixed wide pose emitted from
`params.bodies[idx]` (`pos_w`, `com_w`, and `R_b_to_w` from the fixed `Q_b_to_w`), with
zero velocity/spin. Matches the monolith's frozen `STATIC` body (whose 13 states are held
at their initial pose) while staying stateless so the network needs no state for it. Its
`force_in`/`moment_in` inputs are declared (so joints/segments may deliver to it) but
ignored — the body does not move.
"""
function StaticBody(s, params, idx; name)
    vars, pos, vel, pulley_len_out, pose_R, pose_com, pose_com_vel, pose_omega,
        force_in, moment_in = body_io()
    body = params.bodies[idx]
    R_b_to_w = quaternion_to_rotation_matrix(collect(body.Q_b_to_w))
    eqs = [
        pos ~ collect(body.pos_w)
        vel ~ zeros(3)
        pulley_len_out ~ 0.0
        [pose_R[k] ~ vec(R_b_to_w)[k] for k in 1:9]
        pose_com ~ collect(body.com_w)
        pose_com_vel ~ zeros(3)
        pose_omega ~ zeros(3)
    ]
    return System(eqs, t, vars, param_unknowns(params); name)
end

"""
    DynamicPoint(s, params, idx; name)

Particle vertex: a point mass integrating `pos`/`vel` under the net force gathered
at its `force_in`, its own aerodynamic drag, gravity, and world-frame damping. This
is the free-tether particle and the base for PARTICLE_DYNAMICS wing nodes; the wing
case adds an aero-force input on top of this same integrator (added with wing-node
wiring) rather than duplicating it. Translational mass is `extra_mass` plus the
incident segments' aggregated half-masses (`mass_in`). Its parameters are read from
`params.points[idx]` ([`point_particle_params`](@ref)) — the same param+defaults
source both backends assemble from.
"""
function DynamicPoint(s, params, idx; name, wide = false)
    vars, pos, vel, force_in, mass_in, _, pulley_len_out = point_io()
    pars = point_particle_params(params, idx)
    eqs = [dynamic_point_dynamics(s, pos, vel, force_in, pars.extra_mass + mass_in,
                                  pars);
           pulley_len_out ~ 0.0]
    return finish_vertex(vars, eqs, params; name, wide)
end

"""
    StaticPoint(s, params, idx; name)

Ground-anchored vertex with no dynamic state: `pos` is pinned to `params.points[idx]`'s
`pos_w` and `vel` is zero. Its `force_in`/`mass_in`/`tension_in` inputs are declared
(so segments may deliver to it) but ignored — the point does not move — matching the
STATIC branch of `point_eqs!`.
"""
function StaticPoint(s, params, idx; name, wide = false)
    vars, pos, vel, _, _, _, pulley_len_out = point_io()
    eqs = [
        collect(pos) .~ collect(params.points[idx].pos_w)
        collect(vel) .~ zeros(3)
        pulley_len_out ~ 0.0
    ]
    return finish_vertex(vars, eqs, params; name, wide)
end

"""
    SpringDamperSegment(s, params, idx; name)

Stateless edge reading its two endpoints' `pos`/`vel` and writing the spring-damper
force + tether drag on each (`src_force`/`dst_force`) and each endpoint's half-mass
(`src_mass`/`dst_mass`), in the **positive** force-on-point sign — the assembly sums
these into each point's `force_in`/`mass_in` (ND edge→vertex, or the monolith's
explicit sum). Emits zero tension (a plain/structural segment feeds no pulley/winch).
A [`wing_structural_segment`](@ref) has no tether drag: it drops the drag term (a
distinct compiled type) instead of reading `cd_tether`, so drag stays a single global
setting. Its rest length is the frozen `params.segments[idx].l0`; spring/geometry
parameters come from `params.segments[idx]`, `cd_tether` from `params.set`.
"""
function SpringDamperSegment(s, params, idx; name, wide = false)
    io, extras = wide ? segment_io_wide() : (segment_io(), nothing)
    vars = io[1]
    with_drag = !wing_structural_segment(params.reg.sys_struct, idx)
    spring, wind = segment_spring_params(params, idx; with_drag)
    l0 = params.segments[idx].l0
    force_on_src, force_on_dst, half_mass, _ =
        segment_loads(s, io, l0, spring, wind; with_drag)
    eqs = endpoint_load_eqs(io, force_on_src, force_on_dst, half_mass, 0.0, 0.0)
    if wide
        eqs = [eqs; collect(extras.src_moment) .~ 0; collect(extras.dst_moment) .~ 0]
    end
    return System(eqs, t, vars, param_unknowns(params); name)
end

"""
    pulley_split_eqs(pulley_len, pulley_vel, tension_in, pulley_mass, pulley_damp,
                     pulley_len_out)

The rope-split dynamics a pulley vertex owns on top of its particle motion:
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
    PulleyPoint(s, params, idx, pulley_mass; name)

Dynamic pulley vertex: a particle ([`DynamicPoint`](@ref) motion) that additionally
owns the pulley rope split ([`pulley_split_eqs`](@ref)). `pulley_mass` is the rope
mass driving the split acceleration (supplied by the assembly, which knows the pulley
topology). Its `pulley_damp` is a fixed default parameter.
"""
function PulleyPoint(s, params, idx, pulley_mass; name)
    vars, pos, vel, force_in, mass_in, tension_in, pulley_len_out = point_io()
    extra = @variables pulley_len(t) pulley_vel(t)
    append!(vars, extra)
    pars = point_particle_params(params, idx)
    pulley_damp = make_param(:pulley_damp, 5.0)
    eqs = [
        dynamic_point_dynamics(s, pos, vel, force_in, pars.extra_mass + mass_in, pars);
        pulley_split_eqs(pulley_len, pulley_vel, tension_in, pulley_mass, pulley_damp,
                         pulley_len_out);
    ]
    return System(eqs, t, vars, [param_unknowns(params); pulley_damp]; name)
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

"""
    wing_frame_rotation(zp1, zp2, yp1, yp2)

The body→world rotation matrix `R_b_to_w` fitted from the four ref points, its columns
given by [`wing_frame_columns`](@ref). Shared by the wing-node body-frame damping and
the frozen aero-force rotation.
"""
function wing_frame_rotation(zp1, zp2, yp1, yp2)
    xaxis, yaxis, zaxis = wing_frame_columns(zp1, zp2, yp1, yp2)
    return [xaxis[1] yaxis[1] zaxis[1];
            xaxis[2] yaxis[2] zaxis[2];
            xaxis[3] yaxis[3] zaxis[3]]
end

"""
    body_frame_damp_accel(vel, body_damp, rot, ovel)

The body-frame damping acceleration `R·(coeff ⊙ (Rᵀ·(vel − wing_vel)))`, with the
wing frame `R = rot` ([`wing_frame_rotation`](@ref)) and the wing velocity taken as
the origin velocity `ovel`.
"""
function body_frame_damp_accel(vel, body_damp, rot, ovel)
    return rot * (collect(body_damp) .* (rot' * (collect(vel) .- ovel)))
end

"""
    wing_node_inputs()

The 15 external-input variables a wing node reads from its wing's ref points: the four
ref-point positions `zp1/zp2/yp1/yp2` (12) and the wing origin velocity `ovel` (3),
all as scalar `[input = true]` variables. Both backends supply them — the network
through NetworkDynamics `extin`, the monolith through equations wired to the wing
frame. Returns `(ext, zp1, zp2, yp1, yp2, ovel)` with each ref point a 3-vector.
"""
function wing_node_inputs()
    ext = @variables begin
        zp1x(t), [input = true]; zp1y(t), [input = true]; zp1z(t), [input = true]
        zp2x(t), [input = true]; zp2y(t), [input = true]; zp2z(t), [input = true]
        yp1x(t), [input = true]; yp1y(t), [input = true]; yp1z(t), [input = true]
        yp2x(t), [input = true]; yp2y(t), [input = true]; yp2z(t), [input = true]
        ovx(t), [input = true]; ovy(t), [input = true]; ovz(t), [input = true]
    end
    (zp1x, zp1y, zp1z, zp2x, zp2y, zp2z, yp1x, yp1y, yp1z,
     yp2x, yp2y, yp2z, ovx, ovy, ovz) = ext
    return ext, [zp1x, zp1y, zp1z], [zp2x, zp2y, zp2z], [yp1x, yp1y, yp1z],
           [yp2x, yp2y, yp2z], [ovx, ovy, ovz]
end

"""
    wing_node_extra_accel(point, rot, ovel, vel, mass)

The extra acceleration a KINEMATIC wing node adds on top of the shared
[`point_acceleration`](@ref): its frozen per-point aero force `aero_force_b` (body
frame, refreshed each VSM step) rotated to world by the fitted wing frame `rot`, minus
the body-frame damping ([`body_frame_damp_accel`](@ref)). Returns the world-frame
acceleration vector `aero − damp`.
"""
function wing_node_extra_accel(point, rot, ovel, vel, mass)
    damp = body_frame_damp_accel(vel, point.body_frame_damping, rot, ovel)
    aero = (rot * collect(point.aero_force_b)) ./ mass
    return aero .- damp
end

"""
    WingNodePoint(s, params, idx; name)

A DYNAMIC particle belonging to a `KINEMATIC` wing (`is_wing_node`): the shared
[`DynamicPoint`](@ref) motion plus [`wing_node_extra_accel`](@ref) (frozen aero +
body-frame damping). The wing frame and wing velocity are read through the shared
[`wing_node_inputs`](@ref) (the network supplies them via `extin`). The same kernel
serves aero-only nodes, whose `body_frame_damping` defaults to zero.
"""
function WingNodePoint(s, params, idx; name)
    vars, pos, vel, force_in, mass_in, _, pulley_len_out = point_io()
    ext, zp1, zp2, yp1, yp2, ovel = wing_node_inputs()
    append!(vars, ext)
    pars = point_particle_params(params, idx)
    point = params.points[idx]
    mass = pars.extra_mass + mass_in
    accel = point_acceleration(s, collect(pos), collect(vel), collect(force_in),
        mass, pars.drag_coeff, pars.area, collect(pars.world_damping),
        collect(pars.wind_gnd))
    rot = wing_frame_rotation(zp1, zp2, yp1, yp2)
    eqs = [
        D.(collect(pos)) .~ collect(vel);
        D.(collect(vel)) .~ accel .+ wing_node_extra_accel(point, rot, ovel, vel, mass);
        pulley_len_out ~ 0.0;
    ]
    return System(eqs, t, vars, param_unknowns(params); name)
end

"""
    WingNodePulleyPoint(s, params, idx, pulley_mass; name)

A dynamic pulley vertex ([`PulleyPoint`](@ref)) that also belongs to a `KINEMATIC`
wing, carrying the frozen aero force and body-frame damping of [`WingNodePoint`](@ref).
Used for pulley points that are also wing nodes.
"""
function WingNodePulleyPoint(s, params, idx, pulley_mass; name)
    vars, pos, vel, force_in, mass_in, tension_in, pulley_len_out = point_io()
    extra = @variables pulley_len(t) pulley_vel(t)
    ext, zp1, zp2, yp1, yp2, ovel = wing_node_inputs()
    append!(vars, extra); append!(vars, ext)
    pars = point_particle_params(params, idx)
    point = params.points[idx]
    pulley_damp = make_param(:pulley_damp, 5.0)
    mass = pars.extra_mass + mass_in
    accel = point_acceleration(s, collect(pos), collect(vel), collect(force_in),
        mass, pars.drag_coeff, pars.area, collect(pars.world_damping),
        collect(pars.wind_gnd))
    rot = wing_frame_rotation(zp1, zp2, yp1, yp2)
    eqs = [
        D.(collect(pos)) .~ collect(vel);
        D.(collect(vel)) .~ accel .+ wing_node_extra_accel(point, rot, ovel, vel, mass);
        pulley_split_eqs(pulley_len, pulley_vel, tension_in, pulley_mass, pulley_damp,
                         pulley_len_out);
    ]
    return System(eqs, t, vars, [param_unknowns(params); pulley_damp]; name)
end

"""
    live_aero_node_inputs(num_points)

The extra external inputs a live-aero wing node reads beyond the ref-point frame
([`wing_node_inputs`](@ref)): the wing origin position (`opx/opy/opz`, 3) and every
wing point's world position and velocity (`wpos_k_c`/`wvel_k_c`, `6·num_points`), all
`[input = true]`. The network supplies them via NetworkDynamics `extin`. Returns
`(extra_vars, wing_pos, pos_list, vel_list)` where `pos_list[k]`/`vel_list[k]` are the
3-vectors for wing point `k` (in the `wing_points` order the aero component expects).
"""
function live_aero_node_inputs(num_points)
    opos = @variables opx(t), [input = true]
    append!(opos, @variables opy(t), [input = true])
    append!(opos, @variables opz(t), [input = true])
    wing_pos = collect(opos)
    extra = Any[opos...]
    pos_list = Vector{Vector{Num}}(undef, num_points)
    vel_list = Vector{Vector{Num}}(undef, num_points)
    for k in 1:num_points
        pk = Num[]
        vk = Num[]
        for c in 1:3
            pnm = Symbol(:wpos_, k, :_, c)
            vnm = Symbol(:wvel_, k, :_, c)
            push!(pk, only(@variables $pnm(t), [input = true]))
            push!(vk, only(@variables $vnm(t), [input = true]))
        end
        append!(extra, pk)
        append!(extra, vk)
        pos_list[k] = pk
        vel_list[k] = vk
    end
    return extra, wing_pos, pos_list, vel_list
end

"""
    live_aero_connector_eqs(subsys, s, rot, wing_pos, wind_gnd, pos_list, vel_list)

Wire the shared aero component `subsys`'s particle connectors from the per-point world
positions/velocities (`pos_list`/`vel_list`, read via `extin`), transforming them into
the wing body frame by the fitted rotation `rot` and origin `wing_pos`. Mirrors the
particle wiring of `aero_eqs!` so the network forms the identical body-frame
`point_pos`/`point_vel`/`va`/`rho` the monolith does (apparent wind
`wind(z)·wind_gnd − vel`, density at the clamped height).
"""
function live_aero_connector_eqs(subsys, s, rot, wing_pos, wind_gnd, pos_list, vel_list)
    wind = WindFactor(s.am, s.set.profile_law)
    eqs = Equation[]
    for k in eachindex(pos_list)
        pos_k = pos_list[k]
        vel_k = vel_list[k]
        z_k = pos_k[3]
        va_w = wind(z_k) .* wind_gnd .- vel_k
        eqs = [eqs
               collect(subsys.point_pos[:, k]) .~ rot' * (pos_k .- wing_pos)
               collect(subsys.point_vel[:, k]) .~ rot' * vel_k
               collect(subsys.va[:, k]) .~ rot' * va_w
               subsys.rho[k] ~ calc_rho(s.am, max(0.0, z_k))]
    end
    return eqs
end

"""
    wing_aero_aggregate_vars()

Declare the six scalar observed variables (`wing_aero_force_b_1..3`,
`wing_aero_moment_b_1..3`) that a [`LiveAeroWingNodePoint`](@ref) exposes so the
network state getter can read the wing-level body-frame aero force and moment,
mirroring the monolith's `aero_force_b`/`aero_moment_b` observables. Returns the
force and moment component vectors.
"""
function wing_aero_aggregate_vars()
    force = Num[]
    moment = Num[]
    for c in 1:3
        fnm = Symbol(:wing_aero_force_b_, c)
        mnm = Symbol(:wing_aero_moment_b_, c)
        push!(force, only(@variables $fnm(t)))
        push!(moment, only(@variables $mnm(t)))
    end
    return force, moment
end

"""
    LiveAeroWingNodePoint(s, params, aero_params, idx, wing, slot, num_points; name)

A DYNAMIC particle wing node whose aero force is computed **live** (not frozen): the
shared [`DynamicPoint`](@ref) motion plus body-frame damping, plus the live per-point
aero force from the shared [`aero_component`](@ref) nested as the `aero` subsystem. The
wing frame is fitted from the ref points ([`wing_node_inputs`](@ref)) and the aero
component reads every wing point's world state ([`live_aero_node_inputs`](@ref)), all
via `extin`; the connectors are wired by [`live_aero_connector_eqs`](@ref). Because
NetworkDynamics computes a vertex force in the f-pass (the only place `extin` is
available) the whole wing's aero is evaluated in each wing-node kernel, which selects
its own point force by the per-instance `aero_slot` (so one kernel serves the whole
wing). Reuses `aero_component` unchanged (the same source both backends assemble), with
`aero_params` a separate registry so the namespaced aero parameters do not collide with
the point parameters. Covers every live particle mode ([`ContinuousAero`](@ref),
[`AeroPressure`](@ref) without flap, [`AeroPlate`](@ref)); the frozen
[`AeroDirect`](@ref) and force-free [`AeroNone`](@ref) keep [`WingNodePoint`](@ref).
"""
function LiveAeroWingNodePoint(s, params, aero_params, idx, wing, slot, num_points; name)
    vars, pos, vel, force_in, mass_in, _, pulley_len_out = point_io()
    ext, zp1, zp2, yp1, yp2, ovel = wing_node_inputs()
    append!(vars, ext)
    aero_ext, wing_pos, pos_list, vel_list = live_aero_node_inputs(num_points)
    append!(vars, aero_ext)
    pars = point_particle_params(params, idx)
    point = params.points[idx]
    mass = pars.extra_mass + mass_in
    base = point_acceleration(s, collect(pos), collect(vel), collect(force_in),
        mass, pars.drag_coeff, pars.area, collect(pars.world_damping),
        collect(pars.wind_gnd))
    rot = wing_frame_rotation(zp1, zp2, yp1, yp2)
    wind_gnd = ground_wind_vec(params)
    subsys = aero_component(wing.aero, wing, params.reg.sys_struct;
                            name = :aero, params = aero_params)
    eqs = live_aero_connector_eqs(subsys, s, rot, wing_pos, wind_gnd, pos_list, vel_list)
    aero_slot = make_param(:aero_slot, Float64(slot))
    my_force_b = zeros(Num, 3)
    for k in 1:num_points
        pick = ifelse(abs(aero_slot - k) < 0.5, 1.0, 0.0)
        my_force_b = my_force_b .+ pick .* collect(subsys.point_force[:, k])
    end
    damp = body_frame_damp_accel(vel, point.body_frame_damping, rot, ovel)
    aero = (rot * my_force_b) ./ mass
    wing_force, wing_moment = wing_aero_aggregate_vars()
    append!(vars, [wing_force; wing_moment])
    agg_force = sum(collect(subsys.point_force[:, k]) for k in 1:num_points)
    agg_moment = sum(cross(collect(subsys.point_pos[:, k]),
                           collect(subsys.point_force[:, k])) for k in 1:num_points)
    eqs = [eqs
           D.(collect(pos)) .~ collect(vel)
           D.(collect(vel)) .~ base .+ aero .- damp
           wing_force .~ agg_force
           wing_moment .~ agg_moment
           pulley_len_out ~ 0.0]
    return System(eqs, t, vars, [param_unknowns(params); aero_slot];
                  name, systems = [subsys])
end

"""
    WinchPoint(s, winch, winch_point; name)

Reeling winch vertex at a `STATIC` winch point. Owns the motor speed `winch_vel` and
one `tether_len` state per connected tether (`tether_len_1`, …). The winch motor law
is `winch_component(winch.model, …)` reused verbatim with a fresh `ParamView` (drum
parameters baked as defaults); it reads the summed tether tension
`smooth_norm(tension_in)`, the mean tether length, the control `set_value` and the
`brake`, and returns the drum acceleration `acc`. Integrates
`D(winch_vel) = brake·0 + acc` and `D(tether_len_k) = brake·0 + winch_vel`; each
`tether_len_k` is read by that tether's segments (the network wires it through an
`extin`).
"""
function WinchPoint(s, winch, winch_point; name)
    winch_point.type == STATIC || error(
        "NetworkBackend: winch $(winch.name) is at a non-STATIC point; only " *
        "STATIC winch points are supported so far.")
    n_tethers = length(winch.tether_idxs)
    vars, pos, vel, _, _, tension_in, pulley_len_out = point_io()
    winch_vel, winch_force = @variables winch_vel(t) winch_force(t)
    tether_lens = map(1:n_tethers) do k
        nm = Symbol(:tether_len_, k)
        only(@variables $nm(t))
    end
    append!(vars, [winch_vel, winch_force])
    append!(vars, tether_lens)
    pos_w = make_array_param(:pos_w, zeros(3))
    set_value = make_param(:set_value, 0.0)
    brake = make_param(:brake, 0.0)
    speed_controlled = make_param(:speed_controlled, 0.0)
    pars = [pos_w, set_value, brake, speed_controlled]

    view = ParamView(ParamRegistry(s.sys_struct))
    motor = winch_component(winch.model, s.sys_struct, winch.idx;
                            name = :motor, params = view)
    validate_winch_component(motor, winch)
    winch_acc = ifelse(speed_controlled > 0.5, 0.0, motor.acc)

    eqs = [
        collect(pos) .~ collect(pos_w);
        collect(vel) .~ zeros(3);
        pulley_len_out ~ 0.0;
        winch_force ~ smooth_norm(tension_in);
        motor.vel ~ winch_vel;
        motor.len ~ sum(tether_lens) / n_tethers;
        motor.force ~ winch_force;
        motor.set_value ~ set_value;
        motor.brake ~ brake;
        D(winch_vel) ~ ifelse(brake > 0.5, 0.0, winch_acc);
    ]
    for tl in tether_lens
        push!(eqs, D(tl) ~ ifelse(brake > 0.5, 0.0, winch_vel))
    end
    return System(eqs, t, vars, pars; name, systems = [motor])
end
