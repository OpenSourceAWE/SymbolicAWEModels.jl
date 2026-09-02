# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Scattering the runtime's buffers back into the `SystemStructure`, and
# writing the control setpoints into the parameter buffer. Every value has a fixed
# global slot, resolved once at assembly, so both are plain indexed loops.

"""
    PointReadout(point, pos, vel, drag, wind, force, mass, va, va_frame, aero)

Where one point's results live: its `pos`/`vel` in the output buffer and its
`total_drag`, `wind_vec` and `net_force` in the observable buffer — all three from
the instance that owns its loads, which for a point riding a body is its wrench
half. A name that instance does not observe gets no slots and is left alone.

`mass` is the slot its segments sum their half-masses into, which with `extra_mass`
is the monolith's `point_mass`.

`va` is the point's `va_b` from its [`AeroInflowPoint`](@ref), the apparent wind the
aero is solved on, and is empty for a point that has no aero instance. `va_frame` is
then the body whose frame to rebuild it in ([`va_frame_body`](@ref)).
"""
struct PointReadout
    point::Int
    pos::Vector{Int}
    vel::Vector{Int}
    drag::Vector{Int}
    wind::Vector{Int}
    force::Vector{Int}
    mass::Int
    va::Vector{Int}
    va_frame::Int
    aero::Vector{Int}
end

"""
    SegmentReadout(segment, spring_force, len, l0)

Where one segment's diagnostics live in the observable buffer.
"""
struct SegmentReadout
    segment::Int
    spring_force::Int
    len::Int
    l0::Int
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
    BodyReadout(body, pos, vel, acc, com, com_velocity, orientation, omega_b,
                principal, spin_p)

Where one body's results live: its pose in the output buffer, its acceleration,
body-frame orientation and spin in the observable buffer, and its principal
attitude/spin in the state vector.
"""
struct BodyReadout
    body::Int
    pos::Vector{Int}
    vel::Vector{Int}
    acc::Vector{Int}
    com::Vector{Int}
    com_velocity::Vector{Int}
    orientation::Vector{Int}
    omega_b::Vector{Int}
    principal::Vector{Int}
    spin_p::Vector{Int}
end

"""
    KinematicWingReadout(body, z1, z2, y1, y2, origin, aero_points)

The reference points a fitted (`KINEMATIC`) wing's kinematics are rebuilt from after
a step. Such a wing has no state of its own — its frame, apparent wind and per-point
apparent wind are functions of where its points ended up, which is what
[`wing_kinematics_from_points!`](@ref) computes. Each reference is the wing's own
[`WeightedRefPoints`](@ref), so a blend of several points is rebuilt at the same
weights the [`Wiring`](@ref) feeds the kinematic body kernel.
"""
struct KinematicWingReadout
    body::Int
    z1::WeightedRefPoints
    z2::WeightedRefPoints
    y1::WeightedRefPoints
    y2::WeightedRefPoints
    origin::WeightedRefPoints
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
    WingAeroReadout(wing, force, moment, apparent, wind)

`wing` is the body index, as `Body.idx` is. Where one wing's aero component observes
the quantities the struct carries: its
lumped body-frame `aero_force_b`/`aero_moment_b`, and — for a rigid wing, whose
apparent wind the component computes rather than the getter — its `va_b` and
`wind_vec`. A name the component does not observe gets no slots and is left alone.
"""
struct WingAeroReadout
    wing::Int
    force::Vector{Int}
    moment::Vector{Int}
    apparent::Vector{Int}
    wind::Vector{Int}
end

"""
    RigidWingReadout(wing, frame, alpha_b, base_point)

`wing` is the body index, as `Body.idx` is. Where a `RIGID_DYNAMICS` wing's
reported scalars come from: its body's `frame` output and its observed angular
acceleration, plus the transform base point the elevation and azimuth are measured
against. Its pose, spin and apparent wind are already scattered by the body and
aero readouts, so these are all that is missing to fill the same scalars the
monolith reads out of the integrator.
"""
struct RigidWingReadout
    wing::Int
    frame::Vector{Int}
    alpha_b::Vector{Int}
    base_point::Int
end

"""
    KernelStateGetter(model)

Callable `(integrator, sys_struct)` that scatters the runtime's results back into
the struct, mirroring the monolith's `get_all_state`: point positions, velocities,
apparent wind and drag, segment tension/length/rest length, pulley splits, winch
state, tether lengths, and every wing's reported scalars. It re-runs the output and
observable maps for the integrator's current state first, since the last RHS
evaluation need not correspond to it. A point's `total_mass` comes from the
half-masses its segments gather into it, so it follows their rest lengths as the
monolith's `point_mass` does.

A point's `va_b` is copied from its [`AeroInflowPoint`](@ref) after the fitted-wing
pass, so the value the aero is actually solved on is the compiled model's rather
than a refit of it. A point with no such instance gets the monolith's fallback, its
apparent wind in its own wing's frame, or the first body's if it is not a wing node.
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

function KernelStateGetter(model::KernelModel, rhs, sys_struct)
    system = model.system
    points = [PointReadout(idx,
                  buffer_slots(system, instance, :outputs, :pos),
                  buffer_slots(system, instance, :outputs, :vel),
                  buffer_slots(system, drag_source(model, idx), :observables,
                               :total_drag),
                  observed_slots(system, drag_source(model, idx), :wind_vec),
                  observed_slots(system, drag_source(model, idx), :net_force),
                  mass_slot(system, drag_source(model, idx)),
                  inflow_slots(model, idx), va_frame_body(sys_struct, idx),
                  aero_force_slots(model, idx))
              for (idx, instance) in enumerate(model.point_instances)]
    segments = [SegmentReadout(idx,
                    only(buffer_slots(system, instance, :observables, :spring_force)),
                    only(buffer_slots(system, instance, :observables, :len)),
                    only(buffer_slots(system, instance, :observables, :l0)))
                for (idx, instance) in enumerate(model.segment_instances)]
    pulleys = PulleyReadout[]
    for (idx, role) in enumerate(model.point_roles)
        role.kind === :pulley || continue
        instance = model.point_instances[idx]
        push!(pulleys, PulleyReadout(role.pulley_idx,
            only(buffer_slots(system, instance, :states, :pulley_len)),
            only(buffer_slots(system, instance, :states, :pulley_vel))))
    end
    return KernelStateGetter(rhs, points, segments, pulleys,
                                winch_readouts(model, sys_struct),
                                body_readouts(model, sys_struct),
                                kinematic_wing_readouts(sys_struct),
                                rigid_wing_readouts(model, sys_struct),
                                wing_aero_readouts(model, sys_struct),
                                station_readouts(model), zero(KVec3))
end

"""
    rigid_wing_readouts(model, sys_struct)

One [`RigidWingReadout`](@ref) per `RIGID_DYNAMICS` wing. A fitted wing is covered
by [`kinematic_wing_readouts`](@ref) instead.
"""
function rigid_wing_readouts(model::KernelModel, sys_struct)
    readouts = RigidWingReadout[]
    transforms = sys_struct.transforms
    for wing in sys_struct.wings
        wing.dynamics_type == RIGID_DYNAMICS || continue
        instance = model.body_instances[wing.idx]
        instance == 0 && continue
        base_point = (wing.transform_idx != 0 &&
                      wing.transform_idx <= length(transforms)) ?
            transforms[wing.transform_idx].base_point_idx : 0
        push!(readouts, RigidWingReadout(wing.idx,
            buffer_slots(model.system, instance, :outputs, :frame),
            buffer_slots(model.system, instance, :observables, :alpha_b),
            base_point))
    end
    return readouts
end

"""
    kinematic_wing_readouts(sys_struct)

One [`KinematicWingReadout`](@ref) per `KINEMATIC` body, holding its weighted
reference points and its aerodynamic surface points, resolved once.
"""
function kinematic_wing_readouts(sys_struct)
    readouts = KinematicWingReadout[]
    for (idx, body) in enumerate(sys_struct.bodies)
        body.type == KINEMATIC || continue
        aero = [i for (i, point) in enumerate(sys_struct.points)
                if point.is_wing_node && point.wing_idx == idx]
        push!(readouts, KinematicWingReadout(idx, body.z_ref_points[1],
            body.z_ref_points[2], body.y_ref_points[1],
            body.y_ref_points[2], body.origin, aero))
    end
    return readouts
end

"""
    station_readouts(model) -> Vector{StationReadout}

One [`StationReadout`](@ref) per station that has a twist instance. A
`KINEMATIC` surface has none — its deflection is a [`FlapDelta`](@ref), which the
aero reads directly and the struct does not carry.
"""
function station_readouts(model::KernelModel)
    readouts = StationReadout[]
    for (idx, instance) in enumerate(model.twist_instances)
        instance == 0 && continue
        output(name) = only(buffer_slots(model.system, instance, :outputs, name))
        observed(name) = only(buffer_slots(model.system, instance, :observables, name))
        push!(readouts, StationReadout(idx, output(:twist_angle),
            output(:twist_vel), observed(:tether_force), observed(:tether_moment),
            observed(:aero_moment)))
    end
    return readouts
end

"""
    wing_aero_readouts(model, sys_struct) -> Vector{WingAeroReadout}

One [`WingAeroReadout`](@ref) per wing that has an aero instance, resolving whichever
of the four names that instance observes.
"""
function wing_aero_readouts(model::KernelModel, sys_struct)
    readouts = WingAeroReadout[]
    for (idx, instance) in enumerate(model.aero_instances)
        instance == 0 && continue
        observed(name) = observed_slots(model.system, instance, name)
        push!(readouts, WingAeroReadout(idx, observed(:aero_force_b),
            observed(:aero_moment_b), observed(:va_b), observed(:wind_vel)))
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

"""The body whose frame point `idx`'s `va_b` is expressed in when the point has no
[`AeroInflowPoint`](@ref) to read it from: its own wing, or the first body if it is
not a wing node, which is the fallback `point_eqs!` uses. 0 for a model with no wing,
where the monolith leaves `va_b` at zero. Resolved once here because the choice is
fixed for the life of the model."""
function va_frame_body(sys_struct, idx)
    isempty(sys_struct.wings) && return 0
    point = sys_struct.points[idx]
    return point.is_wing_node ? point.wing_idx : 1
end

"""The slot `instance` gathers its incident segments' half-masses into, or 0 if it
has no `mass_in` input. With `extra_mass` this is the monolith's `point_mass`, which
moves as the segments' rest lengths do."""
function mass_slot(system, instance::Int)
    kernel = system.kernels[system.instances[instance].kernel]
    has_slot(kernel.inputs, :mass_in) || return 0
    return only(buffer_slots(system, instance, :inputs, :mass_in))
end

"""The instance that observes point `idx`'s `total_drag`: its own, or its wrench
half when the point rides a body."""
drag_source(model::KernelModel, idx) =
    model.wrench_instances[idx] == 0 ? model.point_instances[idx] :
    model.wrench_instances[idx]

"""The output slots of point `idx`'s body-frame aerodynamic force, or none when it
has no [`AeroPointForce`](@ref) — a point outside a panel-decomposed wing."""
function aero_force_slots(model::KernelModel, idx)
    instance = model.aero_force_instances[idx]
    instance == 0 && return Int[]
    return buffer_slots(model.system, instance, :outputs, :force_b)
end

"""The output slots of point `idx`'s `va_b`, or none when it has no
[`AeroInflowPoint`](@ref). This is the apparent wind the compiled model solves the
aero on; reading it keeps the struct on the same value rather than a refit."""
function inflow_slots(model::KernelModel, idx)
    instance = model.inflow_instances[idx]
    instance == 0 && return Int[]
    return buffer_slots(model.system, instance, :outputs, :va_b)
end

"""One [`BodyReadout`](@ref) per body. Only a `DYNAMIC` body integrates, so a
clamped or fitted one keeps whatever principal attitude and spin the struct holds —
which for a fitted wing is what its own pose output already implies."""
function body_readouts(model::KernelModel, sys_struct)
    system = model.system
    principal(instance, idx, name) =
        sys_struct.bodies[idx].type == DYNAMIC ?
        buffer_slots(system, instance, :states, name) : Int[]
    return [BodyReadout(idx,
                buffer_slots(system, instance, :outputs, :pos),
                buffer_slots(system, instance, :outputs, :vel),
                buffer_slots(system, instance, :observables, :acc),
                buffer_slots(system, instance, :outputs, :com),
                buffer_slots(system, instance, :outputs, :com_velocity),
                buffer_slots(system, instance, :observables, :orientation),
                buffer_slots(system, instance, :observables, :omega_b),
                principal(instance, idx, :Q),
                principal(instance, idx, :omega_p))
            for (idx, instance) in enumerate(model.body_instances)]
end

"""One [`WinchReadout`](@ref) per winch, resolving its per-tether length slots."""
function winch_readouts(model::KernelModel, sys_struct)
    system = model.system
    readouts = WinchReadout[]
    for (idx, role) in enumerate(model.point_roles)
        role.kind === :winch || continue
        instance = model.point_instances[idx]
        winch = sys_struct.winches[role.winch_idx]
        alias = length(winch.tether_idxs) == 1 ? Symbol("motor₊len") : nothing
        lengths = [(tether, only(winch_state_slots(system, instance,
                                                   Symbol(:tether_len_, k), alias)))
                   for (k, tether) in enumerate(winch.tether_idxs)]
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

function (getter::KernelStateGetter)(integrator, sys_struct::SystemStructure)
    scratch = refresh_outputs!(getter.rhs, integrator.u, integrator.p, integrator.t)
    for readout in getter.points
        point = sys_struct.points[readout.point]
        copy_slots!(point.pos_w, scratch.output, readout.pos)
        copy_slots!(point.vel_w, scratch.output, readout.vel)
        copy_slots!(point.drag_force, scratch.observable, readout.drag)
        copy_slots!(point.wind_vec, scratch.observable, readout.wind)
        copy_slots!(point.force, scratch.observable, readout.force)
        isempty(readout.aero) ||
            copy_slots!(point.aero_force_b, scratch.output, readout.aero)
        readout.mass == 0 || (point.total_mass =
            point.extra_mass + scratch.input[readout.mass])
    end
    for readout in getter.segments
        segment = sys_struct.segments[readout.segment]
        segment.force = scratch.observable[readout.spring_force]
        segment.len = scratch.observable[readout.len]
        segment.l0 = scratch.observable[readout.l0]
    end
    for readout in getter.pulleys
        pulley = sys_struct.pulleys[readout.pulley]
        pulley.len = integrator.u[readout.len]
        pulley.vel = integrator.u[readout.vel]
    end
    for readout in getter.winches
        winch = sys_struct.winches[readout.winch]
        winch.vel = integrator.u[readout.vel]
        winch.acc = scratch.observable[readout.acc]
        winch.friction = scratch.observable[readout.friction]
        winch.set_value = integrator.p.numeric[readout.set_value]
        fill!(winch.force, 0.0)
        winch.force[1] = scratch.observable[readout.force]
        for (tether, slot) in readout.tether_lengths
            sys_struct.tethers[tether].len = integrator.u[slot]
        end
    end
    for readout in getter.bodies
        body = sys_struct.bodies[readout.body]
        copy_slots!(body.pos_w, scratch.output, readout.pos)
        copy_slots!(body.vel_w, scratch.output, readout.vel)
        copy_slots!(body.com_w, scratch.output, readout.com)
        copy_slots!(body.com_vel, scratch.output, readout.com_velocity)
        copy_slots!(body.acc_w, scratch.observable, readout.acc)
        copy_slots!(body.Q_b_to_w, scratch.observable, readout.orientation)
        copy_slots!(body.ω_b, scratch.observable, readout.omega_b)
        copy_slots!(body.Q_p_to_w, integrator.u, readout.principal)
        copy_slots!(body.ω_p, integrator.u, readout.spin_p)
        if body.type == KINEMATIC
            body.Q_p_to_w .= body.Q_b_to_w
            fill!(body.ω_p, 0.0)
        end
    end
    for readout in getter.twist
        surface = sys_struct.stations[readout.surface]
        surface.twist = scratch.output[readout.angle]
        surface.twist_ω = scratch.output[readout.rate]
        surface.tether_force = scratch.observable[readout.tether_force]
        surface.tether_moment = scratch.observable[readout.tether_moment]
        surface.aero_moment = scratch.observable[readout.aero_moment]
    end
    for readout in getter.aero
        wing = sys_struct.bodies[readout.wing]
        copy_slots!(wing.aero_force_b, scratch.observable, readout.force)
        copy_slots!(wing.aero_moment_b, scratch.observable, readout.moment)
        copy_slots!(wing.va_b, scratch.observable, readout.apparent)
        copy_slots!(wing.wind_vec, scratch.observable, readout.wind)
    end
    for readout in getter.kinematic
        wing = sys_struct.bodies[readout.body]
        transforms = sys_struct.transforms
        base_point = (wing.transform_idx != 0 &&
                      wing.transform_idx <= length(transforms)) ?
            transforms[wing.transform_idx].base_point_idx : 0
        wing_kinematics_from_points!(wing, sys_struct.points, sys_struct.set,
            sys_struct.am, sys_struct.wind_mode;
            zp1 = readout.z1, zp2 = readout.z2, yp1 = readout.y1,
            yp2 = readout.y2, origin = readout.origin,
            aero_points = readout.aero_points, base_point,
            stations = sys_struct.stations)
    end
    for readout in getter.rigid
        wing = sys_struct.bodies[readout.wing]
        copy_slots!(wing.R_b_to_w, scratch.output, readout.frame)
        copy_slots!(getter.alpha_b, scratch.observable, readout.alpha_b)
        write_wing_scalars!(wing, sys_struct.points;
            base_point = readout.base_point, alpha_b = getter.alpha_b,
            stations = sys_struct.stations)
    end
    for readout in getter.points
        point = sys_struct.points[readout.point]
        if !isempty(readout.va)
            copy_slots!(point.va_b, scratch.output, readout.va)
        elseif readout.va_frame != 0
            frame = sys_struct.bodies[readout.va_frame].R_b_to_w
            mul!(point.va_b, frame', point.wind_vec .- point.vel_w)
        end
    end
    write_stretched_lengths!(sys_struct)
    return nothing
end

"""Copy `buffer[slots]` into `target`, component by component."""
function copy_slots!(target, buffer, slots)
    @inbounds for k in eachindex(slots)
        target[k] = buffer[slots[k]]
    end
    return nothing
end

"""
    write_stretched_lengths!(sys_struct)

Each tether's stretched length is the sum of its segments' current lengths, which
is what `tether_eqs!` computes symbolically.
"""
function write_stretched_lengths!(sys_struct)
    for tether in sys_struct.tethers
        tether.stretched_len =
            sum(sys_struct.segments[idx].len for idx in tether.segment_idxs)
    end
    return nothing
end

"""
    KernelControlSetter(model)

Callable `(target, values)` writing each winch's control setpoint into the
parameter buffer, mirroring the monolith's `set_set_values`. `values[winch.idx]` is
the setpoint for that winch.
"""
struct KernelControlSetter
    slots::Vector{Tuple{Int, Int}}
end

function KernelControlSetter(model::KernelModel, sys_struct)
    slots = [(readout.winch, readout.set_value)
             for readout in winch_readouts(model, sys_struct)]
    return isempty(slots) ? nothing : KernelControlSetter(slots)
end

function (setter::KernelControlSetter)(target, values)
    for (winch, slot) in setter.slots
        target.p.numeric[slot] = values[winch]
    end
    return nothing
end
