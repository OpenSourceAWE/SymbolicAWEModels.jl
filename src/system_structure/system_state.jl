# Copyright (c) 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# ==================== COLUMN VECTORS ==================== #

"""
    ColumnVec3{V}

Entry `i` of three columns `x`, `y`, `z`, read and written as one 3-vector: a point's
position is `ColumnVec3(X, Y, Z, i)` over a `SysState`'s `X/Y/Z`.
"""
struct ColumnVec3{V <: AbstractVector{SimFloat}} <: StaticVector{3, SimFloat}
    x::V
    y::V
    z::V
    i::Int
end

Base.@propagate_inbounds Base.getindex(c::ColumnVec3, k::Int) =
    k == 1 ? c.x[c.i] : k == 2 ? c.y[c.i] : c.z[c.i]
Base.@propagate_inbounds function Base.setindex!(c::ColumnVec3, value, k::Int)
    k == 1 ? (c.x[c.i] = value) : k == 2 ? (c.y[c.i] = value) : (c.z[c.i] = value)
    return c
end
Base.Tuple(c::ColumnVec3) = (c[1], c[2], c[3])
Base.copy(c::ColumnVec3) = MVector{3, SimFloat}(Tuple(c))
StaticArrays.similar_type(::Type{<:ColumnVec3}, ::Type{T}, ::Size{(3,)}) where {T} =
    SVector{3, T}

"""`ColumnVec3` over entry `i` of the three columns `xyz`."""
column_vec3(xyz::NTuple{3}, i) = ColumnVec3(xyz..., i)

"""
    ColumnQuat{V}

Entry `i` of four columns `w`, `x`, `y`, `z`, read and written as one scalar-first
quaternion: a body's orientation is `ColumnQuat(Qw, Qx, Qy, Qz, frame)`.
"""
struct ColumnQuat{V <: AbstractVector{SimFloat}} <: StaticVector{4, SimFloat}
    w::V
    x::V
    y::V
    z::V
    i::Int
end

Base.@propagate_inbounds Base.getindex(q::ColumnQuat, k::Int) =
    k == 1 ? q.w[q.i] : k == 2 ? q.x[q.i] : k == 3 ? q.y[q.i] : q.z[q.i]
Base.@propagate_inbounds function Base.setindex!(q::ColumnQuat, value, k::Int)
    k == 1 ? (q.w[q.i] = value) : k == 2 ? (q.x[q.i] = value) :
        k == 3 ? (q.y[q.i] = value) : (q.z[q.i] = value)
    return q
end
Base.Tuple(q::ColumnQuat) = (q[1], q[2], q[3], q[4])
Base.copy(q::ColumnQuat) = SimFloat[q[1], q[2], q[3], q[4]]
StaticArrays.similar_type(::Type{<:ColumnQuat}, ::Type{T}, ::Size{(4,)}) where {T} =
    SVector{4, T}

# ==================== RUN COLUMNS ==================== #

"""
    RunColumns

The outputs of a run a `SysState` has no column for yet
(OpenSourceAWE/KiteUtils.jl#166), in its layout: one vector per quantity, or three for a
vector quantity, indexed like the structure's points, segments, stations, winches and
bodies. A wing's are its body's.

$(TYPEDFIELDS)
"""
struct RunColumns
    "net force on each point, world frame [N]"
    point_force::NTuple{3, Vector{SimFloat}}
    "aerodynamic force on each point, its wing's body frame [N]"
    point_aero_force_b::NTuple{3, Vector{SimFloat}}
    "apparent wind at each point, its wing's body frame [m/s]"
    point_va_b::NTuple{3, Vector{SimFloat}}
    "wind at each point, world frame [m/s]"
    point_wind::NTuple{3, Vector{SimFloat}}
    "mass of each point: its own and its segments' half-masses at their rest lengths [kg]"
    point_mass::Vector{SimFloat}
    "entrained-air mass of each point [kg]"
    point_apparent_mass::Vector{SimFloat}
    "rest length of each segment, which a winch reels [m]"
    segment_l0::Vector{SimFloat}
    "tether force on each station [N]"
    station_tether_force::Vector{SimFloat}
    "tether moment on each station [N·m]"
    station_tether_moment::Vector{SimFloat}
    "aerodynamic moment on each station [N·m]"
    station_aero_moment::Vector{SimFloat}
    "force at each winch point, world frame [N]"
    winch_force::NTuple{3, Vector{SimFloat}}
    "reel-out acceleration of each winch [m/s²]"
    winch_acc::Vector{SimFloat}
    "friction force of each winch [N]"
    winch_friction::Vector{SimFloat}
    "acceleration of each body's origin, world frame [m/s²]"
    body_acc::NTuple{3, Vector{SimFloat}}
    "apparent wind at each body, its body frame [m/s]"
    body_va_b::NTuple{3, Vector{SimFloat}}
    "wind at each body, world frame [m/s]"
    body_wind::NTuple{3, Vector{SimFloat}}
    "entrained-air mass of each body [kg]"
    body_apparent_mass::Vector{SimFloat}
    "aerodynamic force on each body, its body frame [N]"
    body_aero_force_b::NTuple{3, Vector{SimFloat}}
    "aerodynamic moment on each body, its body frame [N·m]"
    body_aero_moment_b::NTuple{3, Vector{SimFloat}}
    "turn rate of each body, its body frame [rad/s]"
    body_turn_rate::NTuple{3, Vector{SimFloat}}
    "turn acceleration of each body, its body frame [rad/s²]"
    body_turn_acc::NTuple{3, Vector{SimFloat}}
    "elevation of each body seen from its transform's base [rad]"
    body_elevation::Vector{SimFloat}
    "elevation rate of each body [rad/s]"
    body_elevation_vel::Vector{SimFloat}
    "elevation acceleration of each body [rad/s²]"
    body_elevation_acc::Vector{SimFloat}
    "azimuth of each body seen from its transform's base [rad]"
    body_azimuth::Vector{SimFloat}
    "azimuth rate of each body [rad/s]"
    body_azimuth_vel::Vector{SimFloat}
    "azimuth acceleration of each body [rad/s²]"
    body_azimuth_acc::Vector{SimFloat}
    "heading of each body [rad]"
    body_heading::Vector{SimFloat}
    "course of each body [rad]"
    body_course::Vector{SimFloat}
    "angle of attack of each body [rad]"
    body_aoa::Vector{SimFloat}
end

"""
    RunColumns(; points, segments, stations, winches, bodies)

Zero columns for a run of that many components.
"""
function RunColumns(; points, segments, stations, winches, bodies)
    column(n) = zeros(SimFloat, n)
    columns3(n) = (column(n), column(n), column(n))
    return RunColumns(columns3(points), columns3(points), columns3(points),
        columns3(points), column(points), column(points), column(segments),
        column(stations), column(stations), column(stations), columns3(winches),
        column(winches), column(winches), columns3(bodies), columns3(bodies),
        columns3(bodies), column(bodies), columns3(bodies), columns3(bodies),
        columns3(bodies), columns3(bodies), column(bodies), column(bodies),
        column(bodies), column(bodies), column(bodies), column(bodies), column(bodies),
        column(bodies), column(bodies))
end

# ==================== SYSTEM STATE ==================== #

"""
    ComponentState

Supertype of a record of a run's state: a component of the structure, its `row`, the
run's state `state` and the component's index `idx`. A property the record does not
take from the state is the row's.
"""
abstract type ComponentState end

"""
    SystemState{SS, ST}

What a run holds of a [`SystemStructure`](@ref). Its state is a KiteUtils `SysState` of
`Float64`, `sysstate`, beside `columns` for the outputs that has no column for. Its
`points`, `stations`, `segments`, `pulleys`, `tethers`, `winches`, `bodies` and `wings`
are records of the structure's components, found by index or name, each reading its
state from those columns and every other property from its component. Any other
property, `tubes` among them, is the structure's.
"""
mutable struct SystemState{SS, ST <: SysState}
    const structure::SS
    const sysstate::ST
    const columns::RunColumns
    "the position slot before the first body's in `sysstate.X/Y/Z`"
    const body_slot0::Int
    "the number of wings, whose orientation frames come before the bodies'"
    const n_wings::Int
    "the records of the bodies, built once: a body's type is its row's, which differs
    from body to body"
    bodies::NamedCollection{ComponentState}
    "the records of the wings, the first of `bodies`"
    wings::NamedCollection{ComponentState}
end

"""
    sysstate_counts(sys) -> (P, counts)

The position slot count `P` and the keywords `counts` of `KiteUtils.SysState(P; ...)`
for `sys`, shared by a run's state and the log frames written from it.
"""
function sysstate_counts(sys)
    counts = (wings = length(sys.wings), orients = n_orient_frames(sys),
              deflections = n_flap_deflections(sys), pulleys = length(sys.pulleys),
              winches = length(sys.winches), tethers = length(sys.tethers),
              segments = length(sys.segments))
    return position_slots(sys).total, counts
end

"""
    SystemState(sys::SystemStructure)

The state of a run of `sys`, at its start from the pose `sys` is placed in
([`start_run!`](@ref)).
"""
function SystemState(sys)
    P, counts = sysstate_counts(sys)
    columns = RunColumns(; points=length(sys.points), segments=length(sys.segments),
                         stations=length(sys.stations), winches=length(sys.winches),
                         bodies=length(sys.bodies))
    no_records = NamedCollection(ComponentState[])
    state = SystemState(sys, SysState(P; counts..., precision=SimFloat), columns,
                        first(position_slots(sys).bodies) - 1, length(sys.wings),
                        no_records, no_records)
    bodies = NamedCollection(ComponentState[BodyState(body, state, i)
                                            for (i, body) in enumerate(sys.bodies)])
    setfield!(state, :bodies, bodies)
    setfield!(state, :wings, NamedCollection(ComponentState[bodies[i]
                                                            for i in eachindex(sys.wings)]))
    return start_run!(state)
end

"""
    Records{R, S, C}

The records of one kind of component of a run's state `state`: `R(row, state, i)` for
the `i`-th of the structure's `rows`, by index or name.
"""
struct Records{R, S, C} <: AbstractVector{R}
    state::S
    rows::C
end
Records{R}(state::S, rows::C) where {R, S, C} = Records{R, S, C}(state, rows)

Base.size(records::Records) = size(records.rows)
Base.@propagate_inbounds Base.getindex(records::Records{R}, i::Int) where {R} =
    R(records.rows[i], records.state, i)
Base.getindex(records::Records, name::Symbol) = records[get_idx(records.rows, name)]
Base.haskey(records::Records, name::Symbol) = haskey(records.rows, name)

@inline function Base.getproperty(state::SystemState, name::Symbol)
    name === :points && return records(PointState, state, name)
    name === :stations && return records(StationState, state, name)
    name === :segments && return records(SegmentState, state, name)
    name === :pulleys && return records(PulleyState, state, name)
    name === :tethers && return records(TetherState, state, name)
    name === :winches && return records(WinchState, state, name)
    name === :state_vars && return state_vars(state)
    hasfield(SystemState, name) && return getfield(state, name)
    return getproperty(getfield(state, :structure), name)
end

function Base.setproperty!(state::SystemState, name::Symbol, value)
    name === :state_vars || error("SystemState has no settable property $name")
    set_state_vars!(state, vec(value))
    return value
end

"""The records `R` of the structure's component collection `name`."""
@inline records(R, state::SystemState, name) =
    Records{R}(state, getproperty(getfield(state, :structure), name))

# ==================== RECORDS ==================== #


Base.propertynames(record::ComponentState) =
    (state_properties(typeof(record))..., propertynames(getfield(record, :row))...)

@inline function Base.getproperty(record::ComponentState, name::Symbol)
    name === :row && return getfield(record, :row)
    return state_property(record, Val(name))
end

@inline function Base.setproperty!(record::ComponentState, name::Symbol, value)
    name in state_properties(typeof(record)) ||
        return setproperty!(getfield(record, :row), name, value)
    set_state_property!(record, Val(name), value)
    return value
end

@inline state_property(record::ComponentState, ::Val{name}) where {name} =
    getproperty(getfield(record, :row), name)
set_state_property!(record::ComponentState, ::Val{name}, value) where {name} =
    (state_property(record, Val(name)) .= value)

"""The `SysState` a record's state is held in."""
sysstate(record::ComponentState) = getfield(getfield(record, :state), :sysstate)
"""The [`RunColumns`](@ref) a record's state is held in."""
columns(record::ComponentState) = getfield(getfield(record, :state), :columns)
"""The index of a record's component in its collection."""
index(record::ComponentState) = getfield(record, :idx)

"""
    @state_record Name docstring begin
        property => getter
        property => (getter, setter)
    end

A record type `Name{P, S}` of a run's state and the properties it takes from the state:
`getter` reads one from record `r`, and `setter` writes `value` to a scalar one; a
vector one is written in place through its getter.
"""
macro state_record(name, doc, block)
    properties = Symbol[]
    methods = Expr[]
    for line in block.args
        (line isa Expr && line.head === :call && line.args[1] === :(=>)) || continue
        property, access = QuoteNode(line.args[2]), line.args[3]
        push!(properties, line.args[2])
        getter, setter = access isa Expr && access.head === :tuple ?
            access.args : (access, nothing)
        push!(methods, :(@inline state_property(r::$name, ::Val{$property}) = $getter))
        isnothing(setter) || push!(methods,
            :(@inline set_state_property!(r::$name, ::Val{$property}, value) = $setter))
    end
    return esc(quote
        @doc $doc struct $name{P, S} <: ComponentState
            row::P
            state::S
            idx::Int
        end
        state_properties(::Type{<:$name}) = $(Tuple(properties))
        $(methods...)
    end)
end

"""The error a write to a property derived from others raises."""
derived(name) = error("$name follows from the rest of the run's state; set that instead")

@state_record PointState """
    PointState

A point of a run's state: its position and velocity and the forces on it.
""" begin
    pos_w => ColumnVec3(sysstate(r).X, sysstate(r).Y, sysstate(r).Z, index(r))
    vel_w => ColumnVec3(sysstate(r).VX, sysstate(r).VY, sysstate(r).VZ, index(r))
    drag_force => ColumnVec3(sysstate(r).drag_force_x, sysstate(r).drag_force_y,
                             sysstate(r).drag_force_z, index(r))
    force => column_vec3(columns(r).point_force, index(r))
    aero_force_b => column_vec3(columns(r).point_aero_force_b, index(r))
    va_b => column_vec3(columns(r).point_va_b, index(r))
    wind_vec => column_vec3(columns(r).point_wind, index(r))
    total_mass => (columns(r).point_mass[index(r)],
                   columns(r).point_mass[index(r)] = value)
    apparent_mass => (columns(r).point_apparent_mass[index(r)],
                      columns(r).point_apparent_mass[index(r)] = value)
end

@state_record StationState """
    StationState

A station of a run's state: its twist and the moments on it.
""" begin
    twist => (sysstate(r).twist_angles[index(r)],
              sysstate(r).twist_angles[index(r)] = value)
    twist_ω => (sysstate(r).twist_vel[index(r)], sysstate(r).twist_vel[index(r)] = value)
    tether_force => (columns(r).station_tether_force[index(r)],
                     columns(r).station_tether_force[index(r)] = value)
    tether_moment => (columns(r).station_tether_moment[index(r)],
                      columns(r).station_tether_moment[index(r)] = value)
    aero_moment => (columns(r).station_aero_moment[index(r)],
                    columns(r).station_aero_moment[index(r)] = value)
end

@state_record SegmentState """
    SegmentState

A segment of a run's state: its rest length, which a winch reels, its length between
its points and its spring force.
""" begin
    l0 => (columns(r).segment_l0[index(r)], columns(r).segment_l0[index(r)] = value)
    len => (segment_world_length(getfield(r, :row), getfield(r, :state).points),
            derived(:len))
    force => (sysstate(r).spring_force[index(r)],
              sysstate(r).spring_force[index(r)] = value)
end

@state_record PulleyState """
    PulleyState

A pulley of a run's state: how the rope over it is split.
""" begin
    len => (sysstate(r).pulley_len[index(r)], sysstate(r).pulley_len[index(r)] = value)
    vel => (sysstate(r).pulley_vel[index(r)], sysstate(r).pulley_vel[index(r)] = value)
end

@state_record TetherState """
    TetherState

A tether of a run's state: its unstretched length, and its stretched length along its
segments.
""" begin
    len => (sysstate(r).l_tether[index(r)], sysstate(r).l_tether[index(r)] = value)
    stretched_len => (tether_stretched_length(r), derived(:stretched_len))
end

"""The length of `tether`'s segments end to end in its run [m]."""
function tether_stretched_length(tether)
    segments = getfield(tether, :state).segments
    return sum(segments[idx].len for idx in getfield(tether, :row).segment_idxs; init=0.0)
end

@state_record WinchState """
    WinchState

A winch of a run's state: its reel-out speed, its setpoint and the forces on it.
""" begin
    vel => (sysstate(r).v_reelout[index(r)], sysstate(r).v_reelout[index(r)] = value)
    set_value => (sysstate(r).set_torque[index(r)],
                  sysstate(r).set_torque[index(r)] = value)
    force => column_vec3(columns(r).winch_force, index(r))
    acc => (columns(r).winch_acc[index(r)], columns(r).winch_acc[index(r)] = value)
    friction => (columns(r).winch_friction[index(r)],
                 columns(r).winch_friction[index(r)] = value)
end

"""The position slot of `body` in its run's `sysstate.X/Y/Z`."""
body_slot(body) = getfield(getfield(body, :state), :body_slot0) + index(body)

"""The orientation frame of `body` in its run's `sysstate`."""
body_frame(body) = getfield(getfield(body, :state), :n_wings) + index(body)

"""The centre of mass of body record `r` in its body frame, from its origin [m]."""
com_offset(r) = getfield(r, :row).com_offset_b

@state_record BodyState """
    BodyState

A body of a run's state, a wing among them: its pose and motion in its body frame, the
same in its principal frame, derived from them, and what the air does to it.
""" begin
    pos_w => ColumnVec3(sysstate(r).X, sysstate(r).Y, sysstate(r).Z, body_slot(r))
    vel_w => ColumnVec3(sysstate(r).VX, sysstate(r).VY, sysstate(r).VZ, body_slot(r))
    Q_b_to_w => ColumnQuat(sysstate(r).Qw, sysstate(r).Qx, sysstate(r).Qy,
                           sysstate(r).Qz, body_frame(r))
    ω_b => ColumnVec3(sysstate(r).turn_rate_x, sysstate(r).turn_rate_y,
                      sysstate(r).turn_rate_z, body_frame(r))
    R_b_to_w => (quaternion_to_rotation_matrix(r.Q_b_to_w),
                 r.Q_b_to_w .= rotation_matrix_to_quaternion(value))
    R_p_to_w => (r.R_b_to_w * getfield(r, :row).R_b_to_p', derived(:R_p_to_w))
    com_w => (r.pos_w + r.R_b_to_w * com_offset(r), derived(:com_w))
    com_vel => (r.vel_w + cross(r.R_b_to_w * r.ω_b, r.R_b_to_w * com_offset(r)),
                derived(:com_vel))
    Q_p_to_w => (rotation_matrix_to_quaternion(r.R_p_to_w), derived(:Q_p_to_w))
    ω_p => (getfield(r, :row).R_b_to_p * r.ω_b, derived(:ω_p))
    acc_w => column_vec3(columns(r).body_acc, index(r))
    va_b => column_vec3(columns(r).body_va_b, index(r))
    wind_vec => column_vec3(columns(r).body_wind, index(r))
    apparent_mass => (columns(r).body_apparent_mass[index(r)],
                      columns(r).body_apparent_mass[index(r)] = value)
    aero_force_b => column_vec3(columns(r).body_aero_force_b, index(r))
    aero_moment_b => column_vec3(columns(r).body_aero_moment_b, index(r))
    turn_rate => column_vec3(columns(r).body_turn_rate, index(r))
    turn_acc => column_vec3(columns(r).body_turn_acc, index(r))
    elevation => (columns(r).body_elevation[index(r)],
                  columns(r).body_elevation[index(r)] = value)
    elevation_vel => (columns(r).body_elevation_vel[index(r)],
                      columns(r).body_elevation_vel[index(r)] = value)
    elevation_acc => (columns(r).body_elevation_acc[index(r)],
                      columns(r).body_elevation_acc[index(r)] = value)
    azimuth => (columns(r).body_azimuth[index(r)],
                columns(r).body_azimuth[index(r)] = value)
    azimuth_vel => (columns(r).body_azimuth_vel[index(r)],
                    columns(r).body_azimuth_vel[index(r)] = value)
    azimuth_acc => (columns(r).body_azimuth_acc[index(r)],
                    columns(r).body_azimuth_acc[index(r)] = value)
    heading => (columns(r).body_heading[index(r)],
                columns(r).body_heading[index(r)] = value)
    course => (columns(r).body_course[index(r)], columns(r).body_course[index(r)] = value)
    aoa => (columns(r).body_aoa[index(r)], columns(r).body_aoa[index(r)] = value)
end

"""A record of a rigid wing."""
const RigidWingState = BodyState{<:RigidWing}

"""A record of a particle wing."""
const ParticleWingState = BodyState{<:ParticleWing}

# ==================== STATE VECTOR ==================== #

"""The run's integrated state as one column: every point's position and velocity, every
body's pose and motion in its body frame, every station's twist and its rate, every
pulley's split and its rate, every tether's length and every winch's reel-out speed."""
function state_vars(state::SystemState)
    vars = SimFloat[]
    for point in state.points
        append!(vars, point.pos_w, point.vel_w)
    end
    for body in state.bodies
        append!(vars, body.pos_w, body.vel_w, body.Q_b_to_w, body.ω_b)
    end
    for station in state.stations
        push!(vars, station.twist, station.twist_ω)
    end
    for pulley in state.pulleys
        push!(vars, pulley.len, pulley.vel)
    end
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
        body.pos_w .= next_values!(values, 3)
        body.vel_w .= next_values!(values, 3)
        body.Q_b_to_w .= next_values!(values, 4)
        body.ω_b .= next_values!(values, 3)
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
