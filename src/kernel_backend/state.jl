# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Scattering the runtime's buffers into a run's `SysState`, and
# writing the control setpoints into the parameter buffer. Every value has a fixed
# global slot, resolved once at assembly, so both are plain indexed loops.

"""
    PointReadout(point, pos, vel, drag, wind, force, aero)

Where one point's results live: its `pos`/`vel` in the output buffer and its
`total_drag`, `wind_vec` and `net_force` in the observable buffer — all three from
the instance that owns its loads, which for a point riding a body is its wrench
half — and its world aerodynamic force in the output buffer. A name that instance
does not observe gets no slots and is left alone.
"""
struct PointReadout
    point::Int
    pos::Vector{Int}
    vel::Vector{Int}
    drag::Vector{Int}
    wind::Vector{Int}
    force::Vector{Int}
    aero::Vector{Int}
end

"""
    SegmentReadout(segment, spring_force)

Where one segment's diagnostics live in the observable buffer.
"""
struct SegmentReadout
    segment::Int
    spring_force::Int
end

"""
    PulleyReadout(pulley, len, vel)

Where one pulley's rope split lives in the state vector.
"""
struct PulleyReadout
    pulley::Int
    len::Int
    vel::Int
end

"""
    WinchReadout(winch, vel, force, acc, friction, set_value, tether_lengths)

Where one winch's results live: `vel` and its per-tether lengths in the state
vector, the force/acceleration/friction in the observable buffer, and the control
setpoint in the parameter buffer.
"""
struct WinchReadout
    winch::Int
    vel::Int
    force::Int
    acc::Int
    friction::Int
    set_value::Int
    tether_lengths::Vector{Tuple{Int, Int}}
end

"""
    BodyReadout(body, pos, vel, acc, orientation, omega_b)

Where one body's results live: its pose in the output buffer, and its acceleration,
body-frame orientation and spin in the observable buffer. Its principal frame follows
from these.
"""
struct BodyReadout
    body::Int
    pos::Vector{Int}
    vel::Vector{Int}
    acc::Vector{Int}
    orientation::Vector{Int}
    omega_b::Vector{Int}
end

"""
    KinematicWingReadout(wing, body, frame, aero_points)

The reference points the kinematics of fitted (`KINEMATIC`) wing `wing`, whose body is
`body`, are rebuilt from after a step. Such a wing has no state of its own — its
frame, apparent wind and per-point apparent wind are functions of where its points
ended up, which is what [`wing_kinematics_from_points!`](@ref) computes. `frame` is
the wing's KiteGeometry `WingFrame`, so a blend of several points is rebuilt at the
same weights the [`Wiring`](@ref) feeds the kinematic body kernel.
"""
struct KinematicWingReadout
    wing::Int
    body::Int
    frame::KiteGeometry.WingFrame
    aero_points::Vector{Int}
end

"""
    StationReadout(surface, angle, rate, tether_force, tether_moment, aero_moment)

Where one station's results live: its twist and rate in the output buffer, and
the bridle couple and aerodynamic hinge moment in the observable buffer. A surface
whose twist is prescribed reports the same names, with the couple bound to zero.
"""
struct StationReadout
    surface::Int
    angle::Int
    rate::Int
    tether_force::Int
    tether_moment::Int
    aero_moment::Int
end

"""
    WingAeroReadout(wing, body, force, moment, wind)

Where the aero component of wing `wing`, whose body is `body`, observes what the run
holds of it: its lumped body-frame `aero_force_b`/`aero_moment_b`, and — for a rigid
wing, whose wind the component computes rather than the getter — its body's
`wind_vec`. A name the component does not observe gets no slots and is left alone.
"""
struct WingAeroReadout
    wing::Int
    body::Int
    force::Vector{Int}
    moment::Vector{Int}
    wind::Vector{Int}
end

"""
    RigidWingReadout(wing, alpha_b, base_point)

Where the reported scalars of `RIGID_DYNAMICS` wing `wing` come from: the observed
angular acceleration of its body, plus the transform base point the elevation and
azimuth are measured against. Its pose, spin and apparent wind are already scattered
by the body and aero readouts, so these are all that is missing to fill the same
scalars the monolith reads out of the integrator.
"""
struct RigidWingReadout
    wing::Int
    alpha_b::Vector{Int}
    base_point::Int
end

"""
    KernelStateGetter(model, rhs, sam)

Callable `(integrator, sam)` that scatters the runtime's results into the `SysState`
of the run `sam`, mirroring the monolith's `get_all_state`: point positions,
velocities, wind, loads and drag, segment tension, pulley splits, winch state, tether
lengths, body poses, station twist, and every wing's reported scalars. It re-runs the
output and observable maps for the integrator's current state first, since the last
RHS evaluation need not correspond to it.
"""
struct KernelStateGetter{R}
    rhs::R
    points::Vector{PointReadout}
    segments::Vector{SegmentReadout}
    pulleys::Vector{PulleyReadout}
    winches::Vector{WinchReadout}
    bodies::Vector{BodyReadout}
    kinematic::Vector{KinematicWingReadout}
    rigid::Vector{RigidWingReadout}
    aero::Vector{WingAeroReadout}
    twist::Vector{StationReadout}
    alpha_b::KVec3
end

function KernelStateGetter(model::KernelModel, rhs, sam)
    system = model.system
    points = [PointReadout(i,
                  buffer_slots(system, instance, :outputs, :pos),
                  buffer_slots(system, instance, :outputs, :vel),
                  buffer_slots(system, drag_source(model, i), :observables,
                               :total_drag),
                  observed_slots(system, drag_source(model, i), :wind_vec),
                  observed_slots(system, drag_source(model, i), :net_force),
                  aero_force_slots(model, i))
              for (i, instance) in enumerate(model.point_instances)]
    segments = [SegmentReadout(k,
                    only(buffer_slots(system, instance, :observables, :spring_force)))
                for (k, instance) in enumerate(model.segment_instances)]
    pulleys = PulleyReadout[]
    for (i, role) in enumerate(model.point_roles)
        role.kind === :pulley || continue
        instance = model.point_instances[i]
        push!(pulleys, PulleyReadout(role.pulley_idx,
            only(buffer_slots(system, instance, :states, :pulley_len)),
            only(buffer_slots(system, instance, :states, :pulley_vel))))
    end
    return KernelStateGetter(rhs, points, segments, pulleys, winch_readouts(model, sam),
                             body_readouts(model), kinematic_wing_readouts(sam),
                             rigid_wing_readouts(model, sam),
                             wing_aero_readouts(model, sam), station_readouts(model),
                             zero(KVec3))
end

"""
    rigid_wing_readouts(model, sam)

One [`RigidWingReadout`](@ref) per `RIGID_DYNAMICS` wing of `sam`. A fitted wing is
covered by [`kinematic_wing_readouts`](@ref) instead.
"""
function rigid_wing_readouts(model::KernelModel, sam)
    readouts = RigidWingReadout[]
    for (wing_idx, wing) in enumerate(sam.definition.wings)
        wing_dynamics_type(wing) == RIGID_DYNAMICS || continue
        body_idx = wing_body(sam, wing_idx)
        instance = model.body_instances[body_idx]
        instance == 0 && continue
        push!(readouts, RigidWingReadout(wing_idx,
            buffer_slots(model.system, instance, :observables, :alpha_b),
            transform_base_point(sam.transforms, body_idx)))
    end
    return readouts
end

"""
    kinematic_wing_readouts(sam)

One [`KinematicWingReadout`](@ref) per wing of `sam` whose body is `KINEMATIC`,
holding its weighted reference points and its aerodynamic surface points, resolved
once.
"""
function kinematic_wing_readouts(sam)
    readouts = KinematicWingReadout[]
    for wing_idx in eachindex(sam.definition.wings)
        body_idx = wing_body(sam, wing_idx)
        body_type(sam, body_idx) == KINEMATIC || continue
        push!(readouts, KinematicWingReadout(wing_idx, body_idx,
            wing_frame(sam, wing_idx), collect(Int, wing_points(sam, wing_idx))))
    end
    return readouts
end

"""
    station_readouts(model) -> Vector{StationReadout}

One [`StationReadout`](@ref) per station that has a twist instance. A
`KINEMATIC` surface has none — its deflection is a [`FlapDelta`](@ref), which the
aero reads directly and the run does not carry.
"""
function station_readouts(model::KernelModel)
    readouts = StationReadout[]
    for (station_idx, instance) in enumerate(model.twist_instances)
        instance == 0 && continue
        output(name) = only(buffer_slots(model.system, instance, :outputs, name))
        observed(name) = only(buffer_slots(model.system, instance, :observables, name))
        push!(readouts, StationReadout(station_idx, output(:twist_angle),
            output(:twist_vel), observed(:tether_force), observed(:tether_moment),
            observed(:aero_moment)))
    end
    return readouts
end

"""
    wing_aero_readouts(model, sam) -> Vector{WingAeroReadout}

One [`WingAeroReadout`](@ref) per wing of `sam` that has an aero instance, resolving
whichever of its names that instance observes.
"""
function wing_aero_readouts(model::KernelModel, sam)
    readouts = WingAeroReadout[]
    for (wing_idx, instance) in enumerate(model.aero_instances)
        instance == 0 && continue
        observed(name) = observed_slots(model.system, instance, name)
        push!(readouts, WingAeroReadout(wing_idx, wing_body(sam, wing_idx),
            observed(:aero_force_b), observed(:aero_moment_b), observed(:wind_vel)))
    end
    return readouts
end

"""The slots `name` occupies in `instance`'s observables, or none if it has no such
observable."""
function observed_slots(system, instance::Int, name::Symbol)
    kernel = system.kernels[system.instances[instance].kernel]
    has_slot(kernel.observables, name) || return Int[]
    return buffer_slots(system, instance, :observables, name)
end

"""The instance that observes point `i`'s `total_drag`: its own, or its wrench half
when the point rides a body."""
drag_source(model::KernelModel, i) =
    model.wrench_instances[i] == 0 ? model.point_instances[i] : model.wrench_instances[i]

"""The output slots of point `i`'s world aerodynamic force, or none when it has no
[`AeroPointForce`](@ref) — a point outside a panel-decomposed wing."""
function aero_force_slots(model::KernelModel, i)
    instance, output = model.aero_force_outputs[i]
    instance == 0 && return Int[]
    return buffer_slots(model.system, instance, :outputs, output)
end

"""One [`BodyReadout`](@ref) per body."""
function body_readouts(model::KernelModel)
    system = model.system
    return [BodyReadout(body_idx,
                buffer_slots(system, instance, :outputs, :pos),
                buffer_slots(system, instance, :outputs, :vel),
                buffer_slots(system, instance, :observables, :acc),
                buffer_slots(system, instance, :observables, :orientation),
                buffer_slots(system, instance, :observables, :omega_b))
            for (body_idx, instance) in enumerate(model.body_instances)]
end

"""One [`WinchReadout`](@ref) per winch of `sam`, resolving its per-tether length
slots."""
function winch_readouts(model::KernelModel, sam)
    system = model.system
    readouts = WinchReadout[]
    for (i, role) in enumerate(model.point_roles)
        role.kind === :winch || continue
        instance = model.point_instances[i]
        tethers = sam.definition.winches[role.winch_idx].tethers
        alias = length(tethers) == 1 ? Symbol("motor₊len") : nothing
        lengths = [(tether, only(winch_state_slots(system, instance,
                                                   Symbol(:tether_len_, k), alias)))
                   for (k, tether) in enumerate(tethers)]
        push!(readouts, WinchReadout(role.winch_idx,
            only(winch_state_slots(system, instance, :winch_vel,
                                   Symbol("motor₊vel"))),
            only(buffer_slots(system, instance, :observables, :winch_force)),
            only(buffer_slots(system, instance, :observables, :winch_acc)),
            only(buffer_slots(system, instance, :observables, :winch_friction)),
            only(buffer_slots(system, instance, :params, :set_value)), lengths))
    end
    return readouts
end

"""
    read_station_aero_moment(getter::KernelStateGetter, prob, integrator, s)

The aerodynamic twist moment on station `s` [N·m] at `integrator`'s state, from the
kernels `getter` reads its run's state with; `0.0` for a station without twist.
"""
function read_station_aero_moment(getter::KernelStateGetter, prob, integrator, s)
    readout = findfirst(readout -> readout.surface == s, getter.twist)
    isnothing(readout) && return 0.0
    scratch = refresh_outputs!(getter.rhs, integrator.u, integrator.p, integrator.t)
    return scratch.observable[getter.twist[readout].aero_moment]
end

function (getter::KernelStateGetter)(integrator, sam)
    scratch = refresh_outputs!(getter.rhs, integrator.u, integrator.p, integrator.t)
    output, observable, u = scratch.output, scratch.observable, integrator.u
    columns = run_columns(sam)
    for readout in getter.points
        i = readout.point
        copy_slots!(xyz(columns, :X, :Y, :Z), i, output, readout.pos)
        copy_slots!(xyz(columns, :VX, :VY, :VZ), i, output, readout.vel)
        copy_slots!(xyz(columns, :drag_force_x, :drag_force_y, :drag_force_z), i,
                    observable, readout.drag)
        copy_slots!(xyz(columns, :wind_x, :wind_y, :wind_z), i, observable,
                    readout.wind)
        copy_slots!(xyz(columns, :force_x, :force_y, :force_z), i, observable,
                    readout.force)
        copy_slots!(xyz(columns, :aero_force_x, :aero_force_y, :aero_force_z), i,
                    output, readout.aero)
    end
    for readout in getter.segments
        columns.spring_force[readout.segment] = observable[readout.spring_force]
    end
    for readout in getter.pulleys
        columns.pulley_len[readout.pulley] = u[readout.len]
        columns.pulley_vel[readout.pulley] = u[readout.vel]
    end
    for readout in getter.winches
        columns.v_reelout[readout.winch] = u[readout.vel]
        columns.winch_acc[readout.winch] = observable[readout.acc]
        columns.winch_friction[readout.winch] = observable[readout.friction]
        columns.set_torque[readout.winch] = integrator.p.numeric[readout.set_value]
        for (tether, slot) in readout.tether_lengths
            columns.l_tether[tether] = u[slot]
        end
    end
    for readout in getter.bodies
        slot = body_position_slot(sam, readout.body)
        frame = body_frame(sam, readout.body)
        copy_slots!(xyz(columns, :X, :Y, :Z), slot, output, readout.pos)
        copy_slots!(xyz(columns, :VX, :VY, :VZ), slot, output, readout.vel)
        copy_slots!(xyz(columns, :acc_x, :acc_y, :acc_z), slot, observable, readout.acc)
        copy_slots!(xyz(columns, :Qw, :Qx, :Qy, :Qz), frame, observable,
                    readout.orientation)
        copy_slots!(xyz(columns, :turn_rate_x, :turn_rate_y, :turn_rate_z), frame,
                    observable, readout.omega_b)
    end
    for readout in getter.twist
        columns.twist_angles[readout.surface] = output[readout.angle]
        columns.twist_vel[readout.surface] = output[readout.rate]
    end
    for readout in getter.aero
        copy_slots!(xyz(columns, :aero_force_KA_x, :aero_force_KA_y, :aero_force_KA_z),
                    readout.wing, observable, readout.force)
        copy_slots!(xyz(columns, :aero_moment_KA_x, :aero_moment_KA_y,
                        :aero_moment_KA_z), readout.wing, observable, readout.moment)
        copy_slots!(xyz(columns, :wind_x, :wind_y, :wind_z),
                    body_position_slot(sam, readout.body), observable, readout.wind)
    end
    for readout in getter.kinematic
        wing_kinematics_from_points!(sam, readout.wing, readout.frame,
            readout.aero_points;
            base_point = transform_base_point(sam.transforms, readout.body))
    end
    for readout in getter.rigid
        alpha_b = getter.alpha_b
        @inbounds for k in eachindex(readout.alpha_b)
            alpha_b[k] = observable[readout.alpha_b[k]]
        end
        write_wing_scalars!(sam, readout.wing; base_point = readout.base_point, alpha_b)
    end
    return nothing
end

"""
    copy_slots!(columns, entry, buffer, slots)

Copy `buffer[slots[k]]` into entry `entry` of `columns[k]`, one column of a run's
`SysState` per component; nothing for no slots.
"""
@inline function copy_slots!(columns::Tuple, entry, buffer, slots)
    @inbounds for k in eachindex(slots)
        columns[k][entry] = buffer[slots[k]]
    end
    return nothing
end

"""
    KernelControlSetter(model)

Callable `(target, values)` writing each winch's control setpoint into the
parameter buffer, mirroring the monolith's `set_set_values`. `values[i]` is the
setpoint for winch `i`.
"""
struct KernelControlSetter
    slots::Vector{Tuple{Int, Int}}
end

function KernelControlSetter(model::KernelModel, sam)
    slots = [(readout.winch, readout.set_value)
             for readout in winch_readouts(model, sam)]
    return isempty(slots) ? nothing : KernelControlSetter(slots)
end

function (setter::KernelControlSetter)(target, values)
    for (winch, slot) in setter.slots
        target.p.numeric[slot] = values[winch]
    end
    return nothing
end
