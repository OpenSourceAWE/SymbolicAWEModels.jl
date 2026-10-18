# Copyright (c) 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

"""
    ComponentState

Supertype of what a run holds for one component of a [`SystemStructure`](@ref): the
component itself, its `row`, and the state the run gives it. A property the record does
not hold is the row's.
"""
abstract type ComponentState end

function Base.getproperty(state::ComponentState, name::Symbol)
    hasfield(typeof(state), name) && return getfield(state, name)
    return getproperty(getfield(state, :row), name)
end

function Base.setproperty!(state::ComponentState, name::Symbol, value)
    hasfield(typeof(state), name) || return setproperty!(getfield(state, :row), name, value)
    return setfield!(state, name, convert(fieldtype(typeof(state), name), value))
end

Base.propertynames(state::ComponentState) =
    (fieldnames(typeof(state))..., propertynames(getfield(state, :row))...)

"""
    PointState{P}

A point as a run holds it: where it is and what acts on it.

$(TYPEDFIELDS)
"""
mutable struct PointState{P} <: ComponentState
    "the structure's point"
    const row::P
    "Position in world frame [m]."
    const pos_w::KVec3
    "Velocity in world frame [m/s]."
    const vel_w::KVec3
    "Position relative to the wing's centre of mass in its principal frame [m]: the
    undeformed one carried by the live twist."
    const pos_b::KVec3
    "Net force on the point [N]."
    const force::KVec3
    "Aerodynamic force in body frame [N] (wing nodes of a particle wing)."
    const aero_force_b::KVec3
    "Drag force in world frame [N], the point's own and the segments' given to it."
    const drag_force::KVec3
    "Apparent velocity in body frame [m/s]."
    const va_b::KVec3
    "Wind velocity at the point, world frame [m/s]: a height-profile output, or the
    settable input under `PerPointWind`."
    const wind_vec::KVec3
    "Mass [kg]: `extra_mass` and the half-masses of its segments at their rest lengths."
    total_mass::SimFloat
    "Entrained-air mass [kg] resisting acceleration without adding weight."
    apparent_mass::SimFloat
end
PointState(point) = PointState(point, copy(point.pos_ENU), zeros(KVec3),
                               copy(point.pos_undeformed_b), zeros(KVec3), zeros(KVec3),
                               zeros(KVec3), zeros(KVec3), zeros(KVec3), point.total_mass,
                               point.apparent_mass)

"""
    StationState{S}

A station as a run holds it: its twist and the moments on it.

$(TYPEDFIELDS)
"""
mutable struct StationState{S} <: ComponentState
    "the structure's station"
    const row::S
    "Twist angle [rad]."
    twist::SimFloat
    "Twist angular velocity [rad/s]."
    twist_ω::SimFloat
    "Tether force on the station [N]."
    tether_force::SimFloat
    "Tether moment on the station [N·m]."
    tether_moment::SimFloat
    "Aerodynamic moment on the station [N·m]."
    aero_moment::SimFloat
end
StationState(station) = StationState(station, station.twist, 0.0, 0.0, 0.0, 0.0)

"""
    SegmentState{S}

A segment as a run holds it: its rest length, length and force.

$(TYPEDFIELDS)
"""
mutable struct SegmentState{S} <: ComponentState
    "the structure's segment"
    const row::S
    "Rest length [m]; a winch reels it."
    l0::SimFloat
    "Length [m]."
    len::SimFloat
    "Force [N]."
    force::SimFloat
end
SegmentState(segment) = SegmentState(segment, segment.l0, 0.0, 0.0)

"""
    PulleyState{P}

A pulley as a run holds it: how the rope over it is split.

$(TYPEDFIELDS)
"""
mutable struct PulleyState{P} <: ComponentState
    "the structure's pulley"
    const row::P
    "Length of rope on the first segment's side [m]."
    len::SimFloat
    "Velocity of the rope over the pulley [m/s]."
    vel::SimFloat
end
PulleyState(pulley) = PulleyState(pulley, 0.0, 0.0)

"""
    TetherState{T}

A tether as a run holds it: its reeled-out length.

$(TYPEDFIELDS)
"""
mutable struct TetherState{T} <: ComponentState
    "the structure's tether"
    const row::T
    "Unstretched length [m], the sum of its segments' rest lengths."
    len::SimFloat
    "Stretched length [m]."
    stretched_len::SimFloat
end
TetherState(tether) = TetherState(tether, tether.len, 0.0)

"""
    WinchState{W}

A winch as a run holds it: its reel-out speed, its setpoint and the forces on it.

$(TYPEDFIELDS)
"""
mutable struct WinchState{W} <: ComponentState
    "the structure's winch"
    const row::W
    "Reel-out velocity [m/s]; set it for a speed-controlled winch."
    vel::SimFloat
    "Reel-out acceleration [m/s²]."
    acc::SimFloat
    "Setpoint the winch's model is driven by, the motor torque [N·m] by default."
    set_value::SimFloat
    "Force at the winch point [N]."
    const force::KVec3
    "Friction force [N]."
    friction::SimFloat
end
WinchState(winch) = WinchState(winch, winch.init_vel, 0.0, 0.0, zeros(KVec3), 0.0)

"""
    BodyState{B}

A body as a run holds it, a wing among them: its pose and motion in its body and
principal frames, and what the air does to it.

$(TYPEDFIELDS)
"""
mutable struct BodyState{B} <: ComponentState
    "the structure's body"
    const row::B
    "Body→world orientation quaternion."
    const Q_b_to_w::Vector{SimFloat}
    "Angular velocity in body frame [rad/s]."
    const ω_b::KVec3
    "Body-origin position in world frame [m]."
    const pos_w::KVec3
    "Body-origin velocity in world frame [m/s]."
    const vel_w::KVec3
    "Body-origin acceleration in world frame [m/s²]."
    const acc_w::KVec3
    "Centre-of-mass position in world frame [m]."
    const com_w::KVec3
    "Centre-of-mass velocity in world frame [m/s]."
    const com_vel::KVec3
    "Principal→world orientation quaternion."
    const Q_p_to_w::Vector{SimFloat}
    "Angular velocity in principal frame [rad/s]."
    const ω_p::KVec3
    "Apparent wind in body frame [m/s]."
    const va_b::KVec3
    "Wind velocity at the body, world frame [m/s]: a height-profile output, or the
    settable input under `PerPointWind`."
    const wind_vec::KVec3
    "Aerodynamic force in body frame [N]."
    const aero_force_b::KVec3
    "Aerodynamic moment in body frame [N·m]."
    const aero_moment_b::KVec3
    "Elevation of the body seen from its transform's base [rad]."
    elevation::SimFloat
    "Elevation rate [rad/s]."
    elevation_vel::SimFloat
    "Elevation acceleration [rad/s²]."
    elevation_acc::SimFloat
    "Azimuth of the body seen from its transform's base [rad]."
    azimuth::SimFloat
    "Azimuth rate [rad/s]."
    azimuth_vel::SimFloat
    "Azimuth acceleration [rad/s²]."
    azimuth_acc::SimFloat
    "Heading [rad]."
    heading::SimFloat
    "Turn rate in body frame [rad/s]."
    const turn_rate::KVec3
    "Turn acceleration in body frame [rad/s²]."
    const turn_acc::KVec3
    "Course [rad]."
    course::SimFloat
    "Angle of attack [rad]."
    aoa::SimFloat
    "Entrained-air mass [kg] resisting acceleration without adding weight."
    apparent_mass::SimFloat
end
BodyState(body) = BodyState(body, copy(body.Q_KA_to_ENU), zeros(KVec3),
                            copy(body.pos_ENU), zeros(KVec3), zeros(KVec3), zeros(KVec3),
                            zeros(KVec3), [1.0, 0.0, 0.0, 0.0], zeros(KVec3),
                            zeros(KVec3), zeros(KVec3), zeros(KVec3), zeros(KVec3),
                            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, zeros(KVec3),
                            zeros(KVec3), 0.0, 0.0, body.apparent_mass)

function Base.getproperty(body::BodyState, name::Symbol)
    name === :R_b_to_w && return quaternion_to_rotation_matrix(getfield(body, :Q_b_to_w))
    name === :R_p_to_w && return quaternion_to_rotation_matrix(getfield(body, :Q_p_to_w))
    hasfield(BodyState, name) && return getfield(body, name)
    return getproperty(getfield(body, :row), name)
end

function Base.setproperty!(body::BodyState, name::Symbol, value)
    name === :R_b_to_w || return invoke(setproperty!, Tuple{ComponentState, Symbol, Any},
                                        body, name, value)
    value isa AbstractMatrix ||
        error("R_b_to_w takes a rotation matrix, not a $(typeof(value))")
    getfield(body, :Q_b_to_w) .= rotation_matrix_to_quaternion(value)
    return value
end

"""A record of a rigid wing."""
const RigidWingState = BodyState{<:RigidWing}

"""A record of a particle wing."""
const ParticleWingState = BodyState{<:ParticleWing}

"""
    SystemState{SS}

What a run holds of a [`SystemStructure`](@ref): one record per point, station,
segment, pulley, tether, winch and body, each found by index or name as in the
structure. `wings` holds the records of the bodies that are wings, the same objects as
the first of `bodies`, and `tubes` the structure's tubes. Any other property is the
structure's.
"""
mutable struct SystemState{SS}
    const structure::SS
    const points::NamedCollection{PointState{Point}}
    const stations::NamedCollection{StationState{Station}}
    const segments::NamedCollection{SegmentState{Segment}}
    const pulleys::NamedCollection{PulleyState{Pulley}}
    const tethers::NamedCollection{TetherState{Tether}}
    const winches::NamedCollection{WinchState{Winch}}
    const bodies::NamedCollection{BodyState}
    const wings::NamedCollection{BodyState}
    const tubes::NamedCollection{Tube}
end

function Base.getproperty(state::SystemState, name::Symbol)
    name === :state_vars && return state_vars(state)
    hasfield(SystemState, name) && return getfield(state, name)
    return getproperty(getfield(state, :structure), name)
end

function Base.setproperty!(state::SystemState, name::Symbol, value)
    name === :state_vars || error("SystemState has no settable property $name")
    set_state_vars!(state, vec(value))
    return value
end

"""The run's integrated state as one column: every point's position and velocity,
every body's principal-frame state, every station's twist and its rate, every pulley's
split and its rate, every tether's length and every winch's reel-out speed."""
function state_vars(state::SystemState)
    vars = SimFloat[]
    foreach(point -> append!(vars, point.pos_w, point.vel_w), state.points)
    foreach(body -> append!(vars, body.com_w, body.com_vel, body.Q_p_to_w, body.ω_p),
            state.bodies)
    foreach(station -> push!(vars, station.twist, station.twist_ω), state.stations)
    foreach(pulley -> push!(vars, pulley.len, pulley.vel), state.pulleys)
    foreach(tether -> push!(vars, tether.len), state.tethers)
    foreach(winch -> push!(vars, winch.vel), state.winches)
    return reshape(vars, :, 1)
end

"""Write `vars`, laid out as [`state_vars`](@ref) gives it, into the run's state."""
function set_state_vars!(state::SystemState, vars)
    values = Iterators.Stateful(vars)
    for point in state.points
        point.pos_w .= next_values!(values, 3)
        point.vel_w .= next_values!(values, 3)
    end
    for body in state.bodies
        body.com_w .= next_values!(values, 3)
        body.com_vel .= next_values!(values, 3)
        body.Q_p_to_w .= next_values!(values, 4)
        body.ω_p .= next_values!(values, 3)
    end
    for station in state.stations
        station.twist, station.twist_ω = next_values!(values, 2)
    end
    for pulley in state.pulleys
        pulley.len, pulley.vel = next_values!(values, 2)
    end
    foreach(tether -> tether.len = popfirst!(values), state.tethers)
    foreach(winch -> winch.vel = popfirst!(values), state.winches)
    return state
end

"""The next `n` of the stateful iterator `values`."""
next_values!(values, n) = collect(Iterators.take(values, n))

"""
    SystemState(sys::SystemStructure)

A record for every component of `sys`, holding the component, at the start of a run
from the pose `sys` is placed in ([`start_run!`](@ref)).
"""
function SystemState(sys)
    bodies = NamedCollection(BodyState[BodyState(body) for body in sys.bodies])
    wings = NamedCollection(BodyState[bodies[wing.idx] for wing in sys.wings])
    state = SystemState(sys, records(PointState, sys.points),
                        records(StationState, sys.stations),
                        records(SegmentState, sys.segments),
                        records(PulleyState, sys.pulleys),
                        records(TetherState, sys.tethers),
                        records(WinchState, sys.winches), bodies, wings, sys.tubes)
    return start_run!(state)
end

"""One `State` record per row of `rows`, by index and name."""
records(State, rows) = NamedCollection([State(row) for row in rows])
