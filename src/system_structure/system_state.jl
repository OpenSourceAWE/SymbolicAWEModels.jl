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

# ==================== RUN CACHE ==================== #

"""
    RunCache

What a run keeps that a `SysState` has no column for and that cannot be recomputed from
the rest of its state: the segments' rest lengths and the apparent masses a caller or the
air sets, and what the aerodynamic refresh keeps between refreshes. One vector per
quantity, or three for a vector quantity, indexed like the structure's components. A
wing's are its body's.

$(TYPEDFIELDS)
"""
struct RunCache
    "Each segment's rest length [m], which a winch reels."
    segment_l0::Vector{SimFloat}
    "Each point's apparent mass [kg]."
    point_apparent_mass::Vector{SimFloat}
    "Each body's apparent mass [kg]."
    body_apparent_mass::Vector{SimFloat}
    "Each point's aerodynamic force in its wing's body frame [N]."
    point_aero_force_b::NTuple{3, Vector{SimFloat}}
    "The wind at each body [m/s]."
    body_wind::NTuple{3, Vector{SimFloat}}
    "Each station's aerodynamic twist moment [N·m]."
    station_aero_moment::Vector{SimFloat}
end

"""The zeroed [`RunCache`](@ref) of a structure of `points`, `segments`, `stations` and
`bodies` components."""
function RunCache(; points, segments, stations, bodies)
    column(n) = zeros(SimFloat, n)
    columns3(n) = (column(n), column(n), column(n))
    return RunCache(column(segments), column(points), column(bodies), columns3(points),
                    columns3(bodies), column(stations))
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
`Float64`, `sysstate`, beside `cache` for what it has no column for. Its
`points`, `stations`, `segments`, `pulleys`, `tethers`, `winches`, `bodies` and `wings`
are records of the structure's components, found by index or name, each reading its
state from those columns and every other property from its component. It holds the
run's settings `set`, its atmosphere `am`, its `wind_mode` and its `transforms`. Any
other property, `tubes` among them, is the structure's.
"""
mutable struct SystemState{SS, ST <: SysState}
    const structure::SS
    const sysstate::ST
    const cache::RunCache
    "the segments as this run reads them"
    const built_segments::Vector{BuiltSegment}
    "each point's segments"
    const point_segments::Vector{Vector{Int}}
    "the position slot before the first body's in `sysstate.X/Y/Z`"
    const body_slot0::Int
    "the number of wings, whose orientation frames come before the bodies'"
    const n_wings::Int
    "the records of the bodies, built once: a body's type is its row's, which differs
    from body to body"
    bodies::NamedCollection{ComponentState}
    "the records of the wings, the first of `bodies`"
    wings::NamedCollection{ComponentState}
    "the KiteGeometry transforms placing the run"
    const transforms::Rows{Transform}
    "the run's settings; `nothing` for a state that only draws or places the structure"
    const set::Union{Nothing, Settings}
    "the atmosphere of `set`"
    const am::Union{Nothing, AtmosphericModel}
    "where the wind comes from: [`ProfileWind`](@ref) or [`PerPointWind`](@ref)"
    const wind_mode::WindMode
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
    SystemState(sys::SystemStructure; set=nothing, transforms=Transform[],
                wind_mode=ProfileWind())

The state of a run of `sys` with settings `set` and `wind_mode`, placed by KiteGeometry
`transforms`, at its start from the pose `sys` is placed in ([`start_run!`](@ref)).
"""
function SystemState(sys; set=nothing, transforms=Transform[],
                     wind_mode::WindMode=ProfileWind())
    P, counts = sysstate_counts(sys)
    cache = RunCache(; points=length(sys.points), segments=length(sys.segments),
                     stations=length(sys.stations), bodies=length(sys.bodies))
    no_records = NamedCollection(ComponentState[])
    built_segments = [resolve_segment(row, sys.materials, set) for row in sys.segments]
    state = SystemState(sys, SysState(P; counts..., precision=SimFloat), cache,
                        built_segments, point_segments(sys),
                        first(position_slots(sys).bodies) - 1, length(sys.wings),
                        no_records, no_records, convert(Rows{Transform}, transforms),
                        set, isnothing(set) ? nothing : AtmosphericModel(set), wind_mode)
    bodies = NamedCollection(ComponentState[BodyState(body, state, i)
                                            for (i, body) in enumerate(sys.bodies)])
    setfield!(state, :bodies, bodies)
    setfield!(state, :wings, NamedCollection(ComponentState[bodies[i]
                                                            for i in eachindex(sys.wings)]))
    return start_run!(state)
end

"""The indices of the segments on each point of `sys`."""
function point_segments(sys)
    on_point = [Int[] for _ in sys.points]
    for (i, segment) in enumerate(sys.segments), point_idx in segment.point_idxs
        push!(on_point[point_idx], i)
    end
    return on_point
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
"""The [`RunCache`](@ref) of a record's state."""
cache(record::ComponentState) = getfield(getfield(record, :state), :cache)
"""What a run reads of a segment record's segment."""
built(record) = getfield(getfield(record, :state), :built_segments)[index(record)]
"""The stiffness per unit length [N] of a segment record whose stiffness is a number,
typed `Float64`."""
numeric_unit_stiffness(record) = material_unit_stiffness(built(record).material,
    getfield(record, :row).diameter, to_float64)
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
    force => ColumnVec3(sysstate(r).force_x, sysstate(r).force_y, sysstate(r).force_z,
                        index(r))
    aero_force_b => column_vec3(cache(r).point_aero_force_b, index(r))
    wind_vec => ColumnVec3(sysstate(r).wind_x, sysstate(r).wind_y, sysstate(r).wind_z,
                           index(r))
    va_b => (point_apparent_wind(r), derived(:va_b))
    total_mass => (point_total_mass(r), derived(:total_mass))
    apparent_mass => (cache(r).point_apparent_mass[index(r)],
                      cache(r).point_apparent_mass[index(r)] = value)
end

@state_record StationState """
    StationState

A station of a run's state: its twist and the moments on it.
""" begin
    twist => (sysstate(r).twist_angles[index(r)],
              sysstate(r).twist_angles[index(r)] = value)
    twist_ω => (sysstate(r).twist_vel[index(r)], sysstate(r).twist_vel[index(r)] = value)
    tether_force => (station_tether_load(r).force, derived(:tether_force))
    tether_moment => (station_tether_load(r).moment, derived(:tether_moment))
    aero_moment => (cache(r).station_aero_moment[index(r)],
                    cache(r).station_aero_moment[index(r)] = value)
end

@state_record SegmentState """
    SegmentState

A segment of a run's state: its index, its rest length, which a winch reels, its length
between its points, its spring force and its spring, its stiffness and density those of
its material.
""" begin
    idx => (index(r), derived(:idx))
    unit_stiffness => (material_unit_stiffness(built(r).material, getfield(r, :row).diameter),
                       set_material_stiffness!(built(r).material, value,
                                               getfield(r, :row).diameter))
    density => (built(r).material.parameters["density"],
                built(r).material.parameters["density"] = SimFloat(value))
    unit_damping => (built(r).unit_damping, built(r).unit_damping = value)
    compression_frac => (built(r).compression_frac, built(r).compression_frac = value)
    compression_damping_frac => (built(r).compression_damping_frac,
                                 built(r).compression_damping_frac = value)
    l0 => (cache(r).segment_l0[index(r)], cache(r).segment_l0[index(r)] = value)
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
    force => (winch_force_vector(r), derived(:force))
    acc => (sysstate(r).winch_acc[index(r)], sysstate(r).winch_acc[index(r)] = value)
    friction => (sysstate(r).winch_friction[index(r)],
                 sysstate(r).winch_friction[index(r)] = value)
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
    acc_w => ColumnVec3(sysstate(r).acc_x, sysstate(r).acc_y, sysstate(r).acc_z,
                        body_slot(r))
    wind_vec => column_vec3(cache(r).body_wind, index(r))
    va_b => (body_apparent_wind(r), derived(:va_b))
    apparent_mass => (cache(r).body_apparent_mass[index(r)],
                      cache(r).body_apparent_mass[index(r)] = value)
    aero_force_b => ColumnVec3(sysstate(r).aero_force_KA_x, sysstate(r).aero_force_KA_y,
                               sysstate(r).aero_force_KA_z, index(r))
    aero_moment_b => ColumnVec3(sysstate(r).aero_moment_KA_x, sysstate(r).aero_moment_KA_y,
                                sysstate(r).aero_moment_KA_z, index(r))
    turn_rate => (wing_rates(r).turn_rate, derived(:turn_rate))
    turn_acc => ColumnVec3(sysstate(r).turn_acc_x, sysstate(r).turn_acc_y,
                           sysstate(r).turn_acc_z, body_frame(r))
    elevation => (sysstate(r).elevations[index(r)],
                  sysstate(r).elevations[index(r)] = value)
    elevation_vel => (wing_rates(r).elevation_vel, derived(:elevation_vel))
    elevation_acc => (wing_rates(r).elevation_acc, derived(:elevation_acc))
    azimuth => (sysstate(r).azimuths[index(r)], sysstate(r).azimuths[index(r)] = value)
    azimuth_vel => (wing_rates(r).azimuth_vel, derived(:azimuth_vel))
    azimuth_acc => (wing_rates(r).azimuth_acc, derived(:azimuth_acc))
    heading => (sysstate(r).headings[index(r)], sysstate(r).headings[index(r)] = value)
    course => (sysstate(r).courses[index(r)], sysstate(r).courses[index(r)] = value)
    aoa => (sysstate(r).AoAs[index(r)], sysstate(r).AoAs[index(r)] = value)
end

"""The apparent wind on `point` in the body frame of its wing [m/s], zero where it has
none."""
function point_apparent_wind(point)
    wing_idx = point.wing_idx
    wing_idx > 0 || return zero(SVector{3, SimFloat})
    wing = getfield(point, :state).bodies[wing_idx]
    return SVector{3}(wing.R_b_to_w' * (point.wind_vec - point.vel_w))
end

"""The apparent wind on `body` in its body frame [m/s]."""
body_apparent_wind(body) = SVector{3}(body.R_b_to_w' *
    (body.wind_vec - body.vel_w + getfield(body, :row).wind_disturb))

"""The mass of `point`, its own and half of each segment on it [kg]."""
function point_total_mass(point)
    state = getfield(point, :state)
    mass = point.extra_mass
    for i in getfield(state, :point_segments)[index(point)]
        segment = state.segments[i]
        mass += segment_half_mass(segment.l0, segment.diameter, segment.density)
    end
    return mass
end

"""The force along and moment about the twist axis of `station` from the loads on its
points, zero for a station that does not twist."""
function station_tether_load(station)
    state = getfield(station, :state)
    station.type == KINEMATIC && return (; force = 0.0, moment = 0.0)
    orientation = state.bodies[station.wing_idx].R_b_to_w
    force, moment = 0.0, 0.0
    for point_idx in station.point_idxs
        point = state.points[point_idx]
        couple = twist_bridle_couple(station, point.pos_undeformed_b, station.twist,
                                     orientation)
        point_force = dot(point.force, couple.direction)
        force += point_force
        moment += couple.arm * point_force
    end
    return (; force, moment)
end

"""The spring force of the tethers of `winch` on its point [N]: the pull on the winch
point, negated."""
function winch_force_vector(winch)
    state = getfield(winch, :state)
    winch_point = winch.winch_point_idx
    force = zero(SVector{3, SimFloat})
    for tether_idx in winch.tether_idxs, i in state.tethers[tether_idx].segment_idxs
        segment = state.segments[i]
        ends = segment.point_idxs
        winch_point in ends || continue
        other = state.points[ends[1] == winch_point ? ends[2] : ends[1]]
        away = SVector{3}(other.pos_w - state.points[winch_point].pos_w)
        force -= segment.force * away / norm(away)
    end
    return force
end

"""The turn rate and the rates of elevation and azimuth of `wing`, from its pose, motion
and base point."""
function wing_rates(wing)
    state = getfield(wing, :state)
    base_point = transform_base_point(state.transforms, index(wing))
    rel_pos = base_point == 0 ? collect(wing.pos_w) :
        collect(wing.pos_w) .- collect(state.points[base_point].pos_w)
    R_b_to_w = collect(wing.R_b_to_w)
    e_x = R_b_to_w[:, 1]
    return wing_scalar_kinematics(; rel_pos, e_x, R_t_to_w = sym_calc_R_t_to_w(rel_pos),
        R_v_to_w = calc_R_v_to_w(rel_pos, e_x), R_b_to_w, vel = collect(wing.vel_w),
        acc = collect(wing.acc_w), omega_b = collect(wing.ω_b), alpha_b = zeros(3),
        va_b = collect(wing.va_b), twist_offset = 0.0)
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
