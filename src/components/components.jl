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
function DynamicPoint(s, params, idx; name)
    vars, pos, vel, force_in, mass_in, _, pulley_len_out = point_io()
    pars = point_particle_params(params, idx)
    eqs = [dynamic_point_dynamics(s, pos, vel, force_in, pars.extra_mass + mass_in,
                                  pars);
           pulley_len_out ~ 0.0]
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
explicit sum). Emits zero tension (a plain/structural segment feeds no pulley/winch).
A [`wing_structural_segment`](@ref) has no tether drag: it drops the drag term (a
distinct compiled type) instead of reading `cd_tether`, so drag stays a single global
setting. Its rest length is the frozen `params.segments[idx].l0`; spring/geometry
parameters come from `params.segments[idx]`, `cd_tether` from `params.set`.
"""
function SpringDamperSegment(s, params, idx; name)
    io = segment_io()
    vars = io[1]
    with_drag = !wing_structural_segment(params.reg.sys_struct, idx)
    spring, wind = segment_spring_params(params, idx; with_drag)
    l0 = params.segments[idx].l0
    force_on_src, force_on_dst, half_mass, _ =
        segment_loads(s, io, l0, spring, wind; with_drag)
    eqs = endpoint_load_eqs(io, force_on_src, force_on_dst, half_mass, 0.0, 0.0)
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
