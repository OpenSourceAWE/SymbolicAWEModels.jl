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

A point as a run holds it.
"""
mutable struct PointState{P} <: ComponentState
    "the structure's point"
    const row::P
end

"""
    StationState{S}

A station as a run holds it.
"""
mutable struct StationState{S} <: ComponentState
    "the structure's station"
    const row::S
end

"""
    SegmentState{S}

A segment as a run holds it.
"""
mutable struct SegmentState{S} <: ComponentState
    "the structure's segment"
    const row::S
end

"""
    PulleyState{P}

A pulley as a run holds it.
"""
mutable struct PulleyState{P} <: ComponentState
    "the structure's pulley"
    const row::P
end

"""
    TetherState{T}

A tether as a run holds it.
"""
mutable struct TetherState{T} <: ComponentState
    "the structure's tether"
    const row::T
end

"""
    WinchState{W}

A winch as a run holds it.
"""
mutable struct WinchState{W} <: ComponentState
    "the structure's winch"
    const row::W
end

"""
    BodyState{B}

A body as a run holds it, a wing among them.
"""
mutable struct BodyState{B} <: ComponentState
    "the structure's body"
    const row::B
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
the first of `bodies`. Any other property is the structure's.
"""
struct SystemState{SS}
    structure::SS
    points::NamedCollection{PointState{Point}}
    stations::NamedCollection{StationState{Station}}
    segments::NamedCollection{SegmentState{Segment}}
    pulleys::NamedCollection{PulleyState{Pulley}}
    tethers::NamedCollection{TetherState{Tether}}
    winches::NamedCollection{WinchState{Winch}}
    bodies::NamedCollection{BodyState}
    wings::NamedCollection{BodyState}
end

function Base.getproperty(state::SystemState, name::Symbol)
    hasfield(SystemState, name) && return getfield(state, name)
    return getproperty(getfield(state, :structure), name)
end

"""
    SystemState(sys::SystemStructure)

A record for every component of `sys`, holding the component.
"""
function SystemState(sys)
    bodies = NamedCollection(BodyState[BodyState(body) for body in sys.bodies])
    wings = NamedCollection(BodyState[bodies[wing.idx] for wing in sys.wings])
    return SystemState(sys, records(PointState, sys.points),
                       records(StationState, sys.stations),
                       records(SegmentState, sys.segments),
                       records(PulleyState, sys.pulleys), records(TetherState, sys.tethers),
                       records(WinchState, sys.winches), bodies, wings)
end

"""One `State` record per row of `rows`, by index and name."""
records(State, rows) = NamedCollection([State(row) for row in rows])
