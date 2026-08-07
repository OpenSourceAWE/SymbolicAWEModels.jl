# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Scattering the scheduled runtime's buffers back into the `SystemStructure`, and
# writing the control setpoints into the parameter buffer. Every value has a fixed
# global slot, resolved once at assembly, so both are plain indexed loops.

"""
    PointReadout(point, pos, vel, drag)

Where one point's results live: its `pos`/`vel` in the output buffer and its
`total_drag` in the observable buffer.
"""
struct PointReadout
    point::Int
    pos::Vector{Int}
    vel::Vector{Int}
    drag::Vector{Int}
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
    ScheduledStateGetter(model)

Callable `(integrator, sys_struct)` that scatters the runtime's results back into
the struct, mirroring the monolith's `get_all_state`: point positions, velocities
and drag, segment tension/length/rest length, pulley splits, winch state and tether
lengths. It re-runs the output and observable maps for the integrator's current
state first, since the last RHS evaluation need not correspond to it.
"""
struct ScheduledStateGetter{R}
    rhs::R
    points::Vector{PointReadout}
    segments::Vector{SegmentReadout}
    pulleys::Vector{PulleyReadout}
    winches::Vector{WinchReadout}
    bodies::Vector{BodyReadout}
end

function ScheduledStateGetter(model::ScheduledModel, rhs, sys_struct)
    system = model.system
    points = [PointReadout(idx,
                  buffer_slots(system, instance, :outputs, :pos),
                  buffer_slots(system, instance, :outputs, :vel),
                  buffer_slots(system, drag_source(model, idx), :observables,
                               :total_drag))
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
    return ScheduledStateGetter(rhs, points, segments, pulleys,
                                winch_readouts(model, sys_struct),
                                body_readouts(model, sys_struct))
end

"""The instance that observes point `idx`'s `total_drag`: its own, or its wrench
half when the point rides a body."""
drag_source(model::ScheduledModel, idx) =
    model.wrench_instances[idx] == 0 ? model.point_instances[idx] :
    model.wrench_instances[idx]

"""One [`BodyReadout`](@ref) per body. A clamped body has no state, so its
principal attitude and spin keep whatever the struct already holds."""
function body_readouts(model::ScheduledModel, sys_struct)
    system = model.system
    principal(instance, idx, name) =
        sys_struct.bodies[idx].type == STATIC ? Int[] :
        buffer_slots(system, instance, :states, name)
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
function winch_readouts(model::ScheduledModel, sys_struct)
    system = model.system
    readouts = WinchReadout[]
    for (idx, role) in enumerate(model.point_roles)
        role.kind === :winch || continue
        instance = model.point_instances[idx]
        winch = sys_struct.winches[role.winch_idx]
        lengths = [(tether, only(buffer_slots(system, instance, :states,
                                              Symbol(:tether_len_, k))))
                   for (k, tether) in enumerate(winch.tether_idxs)]
        push!(readouts, WinchReadout(role.winch_idx,
            only(buffer_slots(system, instance, :states, :winch_vel)),
            only(buffer_slots(system, instance, :observables, :winch_force)),
            only(buffer_slots(system, instance, :observables, :winch_acc)),
            only(buffer_slots(system, instance, :observables, :winch_friction)),
            only(buffer_slots(system, instance, :params, :set_value)), lengths))
    end
    return readouts
end

function (getter::ScheduledStateGetter)(integrator, sys_struct::SystemStructure)
    scratch = refresh_outputs!(getter.rhs, integrator.u, integrator.p, integrator.t)
    for readout in getter.points
        point = sys_struct.points[readout.point]
        copy_slots!(point.pos_w, scratch.output, readout.pos)
        copy_slots!(point.vel_w, scratch.output, readout.vel)
        copy_slots!(point.drag_force, scratch.observable, readout.drag)
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
    ScheduledControlSetter(model)

Callable `(target, values)` writing each winch's control setpoint into the
parameter buffer, mirroring the monolith's `set_set_values`. `values[winch.idx]` is
the setpoint for that winch.
"""
struct ScheduledControlSetter
    slots::Vector{Tuple{Int, Int}}
end

function ScheduledControlSetter(model::ScheduledModel, sys_struct)
    slots = [(readout.winch, readout.set_value)
             for readout in winch_readouts(model, sys_struct)]
    return isempty(slots) ? nothing : ScheduledControlSetter(slots)
end

function (setter::ScheduledControlSetter)(target, values)
    for (winch, slot) in setter.slots
        target.p.numeric[slot] = values[winch]
    end
    return nothing
end
