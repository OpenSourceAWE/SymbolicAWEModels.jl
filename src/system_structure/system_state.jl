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

# ==================== SYSTEM STATE ==================== #

"""
    ComponentState

Supertype of a record of a run's state: a component of the structure, its `row`, the
run's state `state` and the component's index `idx`. A property the record does not
take from the state is the row's, read only: the run never writes its definition.
"""
abstract type ComponentState end

"""
    StationGeometry

What a run works out of a station's geometry: where its section sits against its wing's
VSM panels, and its flap at rest.

$(TYPEDFIELDS)
"""
@kwdef mutable struct StationGeometry
    "Leading edge position, body frame [m]."
    le_pos::KVec3 = zeros(KVec3)
    "Chord vector, body frame [m]."
    chord::KVec3 = zeros(KVec3)
    "Spanwise vector, local panel frame."
    y_airf::KVec3 = zeros(KVec3)
    "Indices of the unrefined VSM sections of the station."
    unrefined_section_idxs::Vector{Int64} = Int64[]
    "Indices of the VSM panels whose twist moment the station takes."
    panel_idxs::Vector{Int64} = Int64[]
    "Share of its wing's mass on no point, by the area of its panels [kg]."
    body_mass::SimFloat = 0.0
    "Reference chord directions `[main, flap]` of a flap, each body's frame."
    flap_chord_refs::Vector{KVec3} = KVec3[]
    "Deflection of the flap in its initial pose [rad]."
    flap_rest_delta::SimFloat = 0.0
end

"""
    BodyGeometry

What a run works out of a body: its mass with that of the points it carries
([`update_mass_properties!`](@ref)) and, for the body of a wing, the frame and stations its
wing's row gives.

$(TYPEDFIELDS)
"""
@kwdef mutable struct BodyGeometry
    "Mass of the body with the points it carries [kg]."
    total_mass::SimFloat = 0.0
    "Centre of that mass from the body-frame origin, body frame [m]."
    com_offset_b::KVec3 = zeros(KVec3)
    "Principal moments of inertia about that centre [kg·m²]."
    inertia_principal::KVec3 = zeros(KVec3)
    "Rotation from the body frame into the principal frame."
    R_b_to_p::Matrix{SimFloat} = Matrix{SimFloat}(I, 3, 3)
    "Centre of the body's own mass from the body-frame origin, body frame [m]."
    extra_com_offset_b::KVec3 = zeros(KVec3)
    "The index of the `wings` row whose body this is, 0 for none."
    wing_row::Int = 0
    "The dynamics of that wing, `RIGID_DYNAMICS` for a body that is no wing."
    dynamics_type::WingType = RIGID_DYNAMICS
    "Indices of the wing's stations, from +y to -y."
    station_idxs::Vector{Int} = Int[]
    "Reference points of the body origin, pairs of point index and weight."
    origin::Union{Nothing, Vector{Pair{Int, Float64}}} = nothing
    "Reference points of the body's origin and axes, KiteGeometry's `WingFrame`."
    wing_frame::Union{Nothing, KiteGeometry.WingFrame} = nothing
    "The body's origin in the frame its wing's mesh is authored in [m]."
    mesh_origin_b::KVec3 = zeros(KVec3)
    "The body's axes in the frame its wing's mesh is authored in."
    R_b_to_mesh::Matrix{SimFloat} = Matrix{SimFloat}(I, 3, 3)
end

"""The fields of a [`BodyGeometry`](@ref) a body record reads and writes as its own."""
const BODY_GEOMETRY_FIELDS = fieldnames(BodyGeometry)

"""
    BodyGeometry(sys, body)

What a run starts from of body `body` of definition `sys`: its own centre of mass and, for
the body of a wing, the wing's dynamics, stations and reference points and its pose in the
frame of the wing's mesh, its mass still to be worked out.
"""
function BodyGeometry(sys, body)
    row = sys.bodies[body]
    wing = something(findfirst(wing -> wing.body == body, sys.wings), 0)
    wing == 0 && return BodyGeometry(; extra_com_offset_b=KVec3(own_com_offset(row)))
    spec = KiteGeometry.wing_spec(sys.wings[wing].model)
    origin = spec.origin
    R_b_to_w = quaternion_to_rotation_matrix(row.Q_KA2ENU)
    origin_pos = isnothing(origin) ? KVec3(row.pos_ENU) :
        get_ref_position_from_points(sys.points, origin; field=:pos_ENU)
    mesh_origin, R_mesh_to_w = KVec3(sys.wings[wing].pos_ENU),
        quaternion_to_rotation_matrix(sys.wings[wing].Q_KA2ENU)
    return BodyGeometry(; wing_row=wing, dynamics_type=wing_dynamics_type(sys.wings[wing]),
        extra_com_offset_b=R_b_to_w' * (row.pos_ENU - origin_pos),
        station_idxs=wing_station_idxs(sys, wing), origin, wing_frame=spec.frame,
        mesh_origin_b=R_mesh_to_w' * (origin_pos - mesh_origin),
        R_b_to_mesh=R_mesh_to_w' * R_b_to_w)
end

"""The stations wing `wing` of `sys` lists, from +y to -y."""
wing_station_idxs(sys, wing) = Vector{Int}(something(sys.wings[wing].model.stations, Int[]))

"""The fields of a [`StationGeometry`](@ref) a station record reads and writes as its
own."""
const STATION_GEOMETRY_FIELDS = fieldnames(StationGeometry)

"""The geometry of station `row` a run starts from: the chord and flap references it
gives, the rest zero until the run works them out."""
StationGeometry(row) = StationGeometry(; chord=KVec3(something(row.chord, zeros(3))),
    y_airf=KVec3(something(row.y_airf, zeros(3))),
    flap_chord_refs=KVec3.(something(row.flap_chord_refs, KVec3[])))

"""
    SystemState{SS, ST}

What a run holds of a `SystemDefinition`. Its state is a KiteUtils `SysState` of
`Float64`, `sysstate`, beside the integrator's parameters no column holds. Its
`points`, `spanwise_stations`, `segments`, `pulleys`, `tethers`, `winches`, `bodies` and
`wings`
are records of the structure's components, found by index or name, each reading its
state from those columns and every other property from its component. It holds the
run's settings `set`, its atmosphere `am`, its `wind_mode` and its `transforms`. Any
other property, `tubes` among them, is the structure's.
"""
mutable struct SystemState{SS, ST <: SysState, R <: NamedTuple}
    const structure::SS
    const sysstate::ST
    "the segments as this run reads them where it holds no value of its own"
    const built_segments::Vector{BuiltSegment}
    "each point's segments"
    const point_segments::Vector{Vector{Int}}
    "the wing each point belongs to, `nothing` for none"
    const point_wings::Vector{Union{Nothing, Int}}
    "the tether each segment belongs to, 0 for none"
    const segment_tethers::Vector{Int}
    "the pulley each segment runs over and its side of it, `(0, 0)` for none"
    const segment_pulleys::Vector{Tuple{Int, Int}}
    "the position slot before the first body's in `sysstate.X/Y/Z`"
    const body_slot0::Int
    "the number of wings, whose orientation frames come before the bodies'"
    const n_wings::Int
    "the rows of the points, stations, segments, pulleys, tethers and winches, typed"
    const rows::R
    "the tubes as this run simulates them"
    const tubes::Rows{KiteGeometry.Tube}
    "what the run works out of each body's mass"
    const body_geometry::Vector{BodyGeometry}
    "what the run works out of each station's geometry"
    const station_geometry::Vector{StationGeometry}
    "each body's aero mode, with the model's VSM engine where it has one"
    const aero::Vector{AbstractAeroModel}
    "the records of the bodies, built once: a body's type is its row's, which differs
    from body to body"
    bodies::Rows{ComponentState}
    "the records of the wings, the first of `bodies`"
    wings::Rows{ComponentState}
    "the KiteGeometry transforms placing the run"
    const transforms::Rows{Transform}
    "the run's settings; `nothing` for a state that only draws or places the structure"
    const set::Union{Nothing, Settings}
    "the VSM settings the run's engines are built from"
    const vsm_set::Union{Nothing, VortexStepMethod.VSMSettings}
    "the atmosphere of `set`"
    const am::Union{Nothing, AtmosphericModel}
    "where the wind comes from: [`ProfileWind`](@ref) or [`PerPointWind`](@ref)"
    const wind_mode::WindMode
    "the parameters of the run's integrator, `nothing` before it has one"
    run::Any
end

"""
    sysstate_counts(state) -> (P, counts)

The position slot count `P` and the keywords `counts` of `KiteUtils.SysState(P; ...)`
for run `state`, shared by its `sysstate` and the log frames written from it: its
`gamma_distribution` holds the circulation of every VSM panel of its wings, in order.
"""
function sysstate_counts(sys, n_corners=count_aero_log_points(wing_bodies(sys)),
                         panels=sum(n_vsm_panels, sys.aero; init=0))
    counts = (wings = length(sys.wings), orients = n_orient_frames(sys),
              deflections = n_flap_deflections(sys), pulleys = length(sys.pulleys),
              winches = length(sys.winches), tethers = length(sys.tethers),
              segments = length(sys.segments), panels)
    return position_slots(sys, n_corners).total, counts
end

"""
    SystemState(sys::SystemDefinition; set=nothing, transforms=Transform[],
                vsm_set=nothing, aero=Dict(), wind_mode=ProfileWind())

The state of a run of `sys` with settings `set` and `wind_mode`, placed by KiteGeometry
`transforms`, at its start from the pose `sys` is placed in ([`start_run!`](@ref)), its
wings flown with the aero modes `aero` gives by name and else their rows name.
"""
function SystemState(sys; set=nothing, transforms=Transform[], vsm_set=nothing,
                     aero=Dict{Symbol, AbstractAeroModel}(),
                     wind_mode::WindMode=ProfileWind())
    geometry = [BodyGeometry(sys, body) for body in eachindex(sys.bodies)]
    specs = [wing_aero_spec(wing, aero) for wing in sys.wings]
    vsm_set = model_vsm_set(specs, set, vsm_set)
    modes = AbstractAeroModel[AeroNone() for _ in sys.bodies]
    for (spec, wing) in zip(specs, sys.wings)
        modes[wing.body] = model_aero(spec, wing,
            length(geometry[wing.body].station_idxs), set, vsm_set)
    end
    n_corners = sum(n_aero_log_points(mode, geometry[i])
                    for (i, mode) in enumerate(modes); init=0)
    P, counts = sysstate_counts(sys, n_corners, sum(n_vsm_panels, modes; init=0))
    no_records = Rows{ComponentState}()
    built_segments = [resolve_segment(row, sys.materials, set) for row in sys.segments]
    rows = (points=typed_rows(KiteGeometry.Point{PointModel}, sys.points),
            spanwise_stations=typed_rows(KiteGeometry.SpanwiseStation{StationModel},
                                         sys.spanwise_stations),
            segments=typed_rows(KiteGeometry.Segment{SegmentSpring}, sys.segments),
            pulleys=typed_rows(KiteGeometry.Pulley{PulleyModel}, sys.pulleys),
            tethers=typed_rows(KiteGeometry.Tether{KiteGeometry.NoModel}, sys.tethers),
            winches=typed_rows(KiteGeometry.Winch{<:WinchModel}, sys.winches))
    state = SystemState(sys, SysState(P; counts..., precision=SimFloat),
                        built_segments, point_segments(sys), point_wings(sys, geometry),
                        segment_tethers(sys),
                        segment_pulleys(sys),
                        first(position_slots(sys, n_corners).bodies) - 1,
                        length(sys.wings), rows,
                        run_tubes(sys, geometry), geometry,
                        [StationGeometry(row) for row in sys.spanwise_stations], modes,
                        no_records, no_records, convert(Rows{Transform}, transforms),
                        set, vsm_set,
                        isnothing(set) ? nothing : AtmosphericModel(set), wind_mode,
                        nothing)
    bodies = Rows(ComponentState[BodyState(body, state, i)
                                 for (i, body) in enumerate(sys.bodies)])
    setfield!(state, :bodies, bodies)
    setfield!(state, :wings, Rows(ComponentState[bodies[i] for i in eachindex(sys.wings)]))
    start_lengths!(state)
    init_derived_properties!(state; prn=false)
    foreach(wing -> setup_model_aero!(wing, state), state.wings)
    return start_run!(state)
end

"""`rows` as a collection of `T`."""
typed_rows(::Type{T}, rows) where {T} = Rows{T}(rows)

"""The tubes of `sys` as a run simulates them, the bodies' own centres of mass lying as
their `geometry` gives."""
function run_tubes(sys, geometry)
    tubes = Rows(KiteGeometry.Tube[definition_tube(sys, tube,
                          map(body -> geometry[body].extra_com_offset_b, tube.bodies))
                      for tube in sys.tubes])
    foreach(((i, tube),) -> tube.idx = i, enumerate(tubes))
    return tubes
end

"""`vsm_set`, else the project's `vsm_settings` for a run of `set` with wings flown with
the aero mode `specs` of which one is VSM's."""
function model_vsm_set(specs, set, vsm_set)
    isnothing(vsm_set) || return vsm_set
    (isnothing(set) || !any(spec -> spec isa AeroSpec, specs)) && return nothing
    path = optional_project_file("vsm_settings")
    return isnothing(path) ? nothing : VortexStepMethod.VSMSettings(path; data_prefix=false)
end

"""The indices of the segments on each point of `sys`."""
function point_segments(sys)
    on_point = [Int[] for _ in sys.points]
    for (i, segment) in enumerate(sys.segments), point_idx in segment.points
        push!(on_point[point_idx], i)
    end
    return on_point
end

"""The wing each point of `sys` belongs to, `nothing` for none: the wing of a station it is
on, else the wing whose body it rides, `geometry` giving each body's wing."""
function point_wings(sys, geometry)
    wings = Union{Nothing, Int}[nothing for _ in sys.points]
    for (i, point) in enumerate(sys.points)
        body = point.body
        isnothing(body) || geometry[body].wing_row == 0 ||
            (wings[i] = geometry[body].wing_row)
    end
    for station in sys.spanwise_stations, point in station.points
        wings[point] = station.wing
    end
    return wings
end

"""The pulley each segment of `sys` runs over and its side of it, `(0, 0)` for none."""
function segment_pulleys(sys)
    over = fill((0, 0), length(sys.segments))
    for (i, pulley) in enumerate(sys.pulleys), (side, segment) in enumerate(pulley.segments)
        over[segment] = (i, side)
    end
    return over
end

"""The index of the tether each segment of `sys` belongs to, 0 for none."""
function segment_tethers(sys)
    tether_of = zeros(Int, length(sys.segments))
    for (i, tether) in enumerate(sys.tethers), segment in tether.segments
        tether_of[segment] = i
    end
    return tether_of
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
Base.getindex(records::Records, name::Symbol) =
    records[KiteGeometry.row_index(records.rows, name)]
Base.haskey(records::Records, name::Symbol) = haskey(records.rows, name)

@inline function Base.getproperty(state::SystemState, name::Symbol)
    rows = getfield(state, :rows)
    name === :points && return Records{PointState}(state, rows.points)
    name === :spanwise_stations &&
        return Records{StationState}(state, rows.spanwise_stations)
    name === :segments && return Records{SegmentState}(state, rows.segments)
    name === :pulleys && return Records{PulleyState}(state, rows.pulleys)
    name === :tethers && return Records{TetherState}(state, rows.tethers)
    name === :winches && return Records{WinchState}(state, rows.winches)
    name === :state_vars && return state_vars(state)
    hasfield(SystemState, name) && return getfield(state, name)
    return getproperty(getfield(state, :structure), name)
end

function Base.setproperty!(state::SystemState, name::Symbol, value)
    name === :state_vars || error("SystemState has no settable property $name")
    set_state_vars!(state, vec(value))
    return value
end

# ==================== RECORDS ==================== #


Base.propertynames(record::ComponentState) =
    (state_properties(typeof(record))..., propertynames(getfield(record, :row))...)

@inline function Base.getproperty(record::ComponentState, name::Symbol)
    name === :row && return getfield(record, :row)
    name === :idx && return index(record)
    return state_property(record, Val(name))
end

@inline function Base.setproperty!(record::ComponentState, name::Symbol, value)
    name in state_properties(typeof(record)) || error("$(getfield(record, :row).name)'s " *
        "$name is its definition's: set it on the definition, not on a run's state")
    set_state_property!(record, Val(name), value)
    return value
end

@inline state_property(record::ComponentState, ::Val{name}) where {name} =
    getproperty(getfield(record, :row), name)
set_state_property!(record::ComponentState, ::Val{name}, value) where {name} =
    (state_property(record, Val(name)) .= value)

"""The `SysState` a record's state is held in."""
sysstate(record::ComponentState) = getfield(getfield(record, :state), :sysstate)
"""The run's state a record reads."""
run_state(record::ComponentState) = getfield(record, :state)
"""What a run reads of a segment record's segment."""
built(record) = getfield(getfield(record, :state), :built_segments)[index(record)]
"""The stiffness per unit length [N] of a segment record whose stiffness is a number,
typed `Float64`."""
numeric_unit_stiffness(record) = to_float64(record.unit_stiffness)
"""The index of a record's component in its collection."""
index(record::ComponentState) = getfield(record, :idx)

"""
    run_value(record, name, default)

The value of `record`'s `name` the run's integrator holds, else `default`: before the run
has an integrator, or where no parameter of it is read from `name`. A vector comes back
immutable, as writing into it would not reach the integrator.
"""
function run_value(record, name, default)
    run = getfield(run_state(record), :run)
    isnothing(run) && return default
    value = run_param(run, (collection(record), index(record), name))
    isnothing(value) && return default
    value isa AbstractVector{<:Real} && length(value) == 3 && return SVector{3}(value)
    return value
end

"""
    set_run_value!(record, name, value)

Write `value` to `record`'s `name` in the run's integrator where it is one of its
parameters, else to the column of its row, refusing a value only a run holds before it
has started.
"""
function set_run_value!(record, name, value)
    run = getfield(run_state(record), :run)
    isnothing(run) ||
        set_run_param!(run, (collection(record), index(record), name), value) &&
        return value
    row = getfield(record, :row)
    hasproperty(row, name) ||
        error("$(row.name)'s $name is a run's: set it once the run has started")
    setproperty!(row, name, value)
    return value
end

"""The rest length of segment record `r` [m]: its tether's unstretched length shared
over its segments, its side of the rope over its pulley, else what its run holds, else
its row's."""
function segment_rest_length(r)
    state = run_state(r)
    tether = getfield(state, :segment_tethers)[index(r)]
    if tether > 0
        segments = getfield(state, :structure).tethers[tether].segments
        return sysstate(r).l_tether[tether] / length(segments)
    end
    pulley, side = getfield(state, :segment_pulleys)[index(r)]
    if pulley > 0
        len = sysstate(r).pulley_len[pulley]
        return side == 1 ? len : state.pulleys[pulley].sum_len - len
    end
    return run_value(r, :rest_length, something(getfield(r, :row).rest_length, 0.0))
end

"""Whether point record `r` is a structural node of a wing ([`wing_node`](@ref))."""
point_is_wing_node(r) = wing_node(getfield(run_state(r), :structure), index(r))

"""Where point record `r` sits in the frame of the body it rides [m]: its `anchor_b`, else
where its initial position puts it; zero for a point on no body."""
function point_anchor_b(r)
    anchor = getfield(r, :row).anchor_b
    isnothing(anchor) || return KVec3(anchor)
    body = r.body
    isnothing(body) && return zero(KVec3)
    return KVec3(initial_body_position(run_state(r).bodies[body], getfield(r, :row)))
end

"""
    point_undeformed_b(r)

Where wing-node record `r` sits relative to its wing in the initial pose, wing frame [m]:
from a rigid wing's centre of mass, or from a particle wing's frame fitted to the initial
positions of its reference points; zero for a point that is no wing node.
"""
function point_undeformed_b(r)
    wing_idx = r.wing_idx
    (isnothing(wing_idx) || !r.is_wing_node) && return zero(KVec3)
    state = run_state(r)
    wing = state.bodies[wing_idx]
    pos = KVec3(getfield(r, :row).pos_ENU)
    if body_dynamics_type(wing) == RIGID_DYNAMICS
        R = initial_rotation(wing)
        return KVec3(R' * (pos - (origin_ENU(wing) .+ R * wing.com_offset_b)))
    end
    isnothing(wing.wing_frame) &&
        return zero(KVec3)
    R, origin = calc_particle_dynamics_wing_frame(state.points, wing.wing_frame;
                                                  field=:pos_ENU)
    return KVec3(R' * (pos - origin))
end

"""The position along the beam element of the tube point record `r` rides and its offset
off the centreline ([`beam_anchor`](@ref)), `(0, 0)` for a point on no tube."""
function point_beam_anchor(r)
    tube = r.tube
    isnothing(tube) && return (0.0, zero(KVec3))
    state = run_state(r)
    return beam_anchor(getfield(r, :row).pos_ENU, state.tubes[tube], state.bodies)
end

"""The length of the rope over pulley record `r` [m]: its two segments' rest lengths in the
definition."""
function pulley_sum_len(r)
    segments = getfield(run_state(r), :structure).segments
    return sum(something(segments[idx].rest_length, 0.0) for idx in r.segments)
end

"""Winch record `r`'s column `column`, else the run's setting `setting`, else `NaN`."""
function winch_setting(r, column, setting)
    value = getproperty(getfield(r, :row), column)
    isnothing(value) || return value
    set = getfield(run_state(r), :set)
    return isnothing(set) ? NaN : getproperty(set, setting)
end

"""The aerodynamic twist moment on station record `r` [N·m] at its run's integrator
state, `0.0` before the run has started."""
function station_aero_moment(r)
    run = getfield(run_state(r), :run)
    return isnothing(run) ? 0.0 : run_station_aero_moment(run, index(r))
end

"""
    @state_record Name docstring begin
        property => getter
        property => (getter, setter)
    end holder => FIELDS ...

A record type `Name{P, S}` of a run's state and the properties it takes from the state:
`getter` reads one from record `r`, and `setter` writes `value` to a scalar one; a
vector one is written in place through its getter. Each `holder => FIELDS` adds the
fields named in the constant `FIELDS`, read and written on `holder(r)`.
"""
macro state_record(name, doc, block, forwards...)
    accessors = Tuple{Symbol, Any}[(line.args[2], line.args[3]) for line in block.args
        if line isa Expr && line.head === :call && line.args[1] === :(=>)]
    for forward in forwards
        holder, fields = forward.args[2], getfield(__module__, forward.args[3])
        for field in fields
            push!(accessors, (field, :((getproperty($holder(r), $(QuoteNode(field))),
                setproperty!($holder(r), $(QuoteNode(field)), value)))))
        end
    end
    properties = first.(accessors)
    methods = Expr[]
    for (property, access) in accessors
        getter, setter = access isa Expr && access.head === :tuple ?
            access.args : (access, nothing)
        key = QuoteNode(property)
        push!(methods, :(@inline state_property(r::$name, ::Val{$key}) = $getter))
        isnothing(setter) || push!(methods,
            :(@inline set_state_property!(r::$name, ::Val{$key}, value) = $setter))
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
    aero_force => ColumnVec3(sysstate(r).aero_force_x, sysstate(r).aero_force_y,
                             sysstate(r).aero_force_z, index(r))
    aero_force_b => (point_aero_force_b(r), derived(:aero_force_b))
    wind_vec => ColumnVec3(sysstate(r).wind_x, sysstate(r).wind_y, sysstate(r).wind_z,
                           index(r))
    va_b => (point_apparent_wind(r), derived(:va_b))
    total_mass => (point_total_mass(r), derived(:total_mass))
    apparent_mass => (run_value(r, :apparent_mass, something(r.row.apparent_mass, 0.0)),
                      set_run_value!(r, :apparent_mass, value))
    extra_mass => (run_value(r, :extra_mass, something(r.row.extra_mass, 0.0)),
                   set_run_value!(r, :extra_mass, value))
    body_frame_damping => (run_value(r, :body_frame_damping,
                                     per_axis(r.row.body_frame_damping)),
                           set_run_value!(r, :body_frame_damping, value))
    world_frame_damping => (run_value(r, :world_frame_damping,
                                      per_axis(r.row.world_frame_damping)),
                            set_run_value!(r, :world_frame_damping, value))
    area => (run_value(r, :area, something(r.row.area, 0.0)),
             set_run_value!(r, :area, value))
    drag_coeff => (run_value(r, :drag_coeff, something(r.row.drag_coeff, 0.0)),
                   set_run_value!(r, :drag_coeff, value))
    fix_sphere => (run_value(r, :fix_sphere, something(r.row.fix_sphere, false)),
                   set_run_value!(r, :fix_sphere, value))
    fix_static => (run_value(r, :fix_static, something(r.row.fix_static, false)),
                   set_run_value!(r, :fix_static, value))
    ext_force_w => (run_value(r, :ext_force_w, zero(KVec3)),
                    set_run_value!(r, :ext_force_w, value))
    is_wing_node => (point_is_wing_node(r), derived(:is_wing_node))
    wing_idx => (getfield(run_state(r), :point_wings)[index(r)], derived(:wing_idx))
    anchor_b => (point_anchor_b(r), derived(:anchor_b))
    pos_undeformed_b => (point_undeformed_b(r), derived(:pos_undeformed_b))
    beam_frac => (first(point_beam_anchor(r)), derived(:beam_frac))
    beam_offset_b => (last(point_beam_anchor(r)), derived(:beam_offset_b))
end

@state_record StationState """
    StationState

A station of a run's state: its twist, the moments on it, and its geometry as the run
works it out.
""" begin
    twist => (sysstate(r).twist_angles[index(r)],
              sysstate(r).twist_angles[index(r)] = value)
    twist_ω => (sysstate(r).twist_vel[index(r)], sysstate(r).twist_vel[index(r)] = value)
    tether_force => (station_tether_load(r).force, derived(:tether_force))
    tether_moment => (station_tether_load(r).moment, derived(:tether_moment))
    aero_moment => (station_aero_moment(r), set_run_value!(r, :aero_moment, value))
    moment_frac => (run_value(r, :moment_frac, something(r.row.moment_frac, 0.0)),
                    set_run_value!(r, :moment_frac, value))
    damping => (run_value(r, :damping, something(r.row.damping, 50.0)),
                set_run_value!(r, :damping, value))
    stiffness => (run_value(r, :stiffness, something(r.row.stiffness, 0.0)),
                  set_run_value!(r, :stiffness, value))
    flap_axis => (run_value(r, :flap_axis,
                            something(r.row.flap_axis, SVector(0.0, 1.0, 0.0))),
                  set_run_value!(r, :flap_axis, value))
    area => (something(r.row.area, NaN), derived(:area))
end station_geometry => STATION_GEOMETRY_FIELDS

"""The [`StationGeometry`](@ref) of station record `r`."""
station_geometry(r) = getfield(run_state(r), :station_geometry)[index(r)]

@state_record SegmentState """
    SegmentState

A segment of a run's state: its index, its rest length, which a winch reels, its length
between its points, its spring force and its spring, its stiffness and density those of
its material.
""" begin
    unit_stiffness => (run_value(r, :unit_stiffness,
                                 material_unit_stiffness(built(r).material, r.diameter)),
                       set_run_value!(r, :unit_stiffness, value))
    density => (run_value(r, :density, built(r).material.parameters["density"]),
                set_run_value!(r, :density, value))
    unit_damping => (run_value(r, :unit_damping, built(r).unit_damping),
                     set_run_value!(r, :unit_damping, value))
    compression_frac => (run_value(r, :compression_frac, built(r).compression_frac),
                         set_run_value!(r, :compression_frac, value))
    compression_damping_frac => (run_value(r, :compression_damping_frac,
                                           built(r).compression_damping_frac),
                                 set_run_value!(r, :compression_damping_frac, value))
    rest_length => (segment_rest_length(r), set_run_value!(r, :rest_length, value))
    len => (segment_world_length(getfield(r, :row), getfield(r, :state).points),
            derived(:len))
    force => (sysstate(r).spring_force[index(r)],
              sysstate(r).spring_force[index(r)] = value)
    diameter => (run_value(r, :diameter, something(r.row.diameter, 0.0)),
                 set_run_value!(r, :diameter, value))
end

@state_record PulleyState """
    PulleyState

A pulley of a run's state: how the rope over it is split.
""" begin
    len => (sysstate(r).pulley_len[index(r)], sysstate(r).pulley_len[index(r)] = value)
    vel => (sysstate(r).pulley_vel[index(r)], sysstate(r).pulley_vel[index(r)] = value)
    efficiency => (run_value(r, :efficiency, something(r.row.efficiency, 0.95)),
                   set_run_value!(r, :efficiency, value))
    damping => (run_value(r, :damping, something(r.row.damping, 0.0)),
                set_run_value!(r, :damping, value))
    brake => (run_value(r, :brake, something(r.row.brake, false)),
              set_run_value!(r, :brake, value))
    friction_epsilon => (run_value(r, :friction_epsilon,
                                   something(r.row.friction_epsilon, 0.1)),
                         set_run_value!(r, :friction_epsilon, value))
    sum_len => (pulley_sum_len(r), derived(:sum_len))
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
    return sum(segments[idx].len for idx in getfield(tether, :row).segments; init=0.0)
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
    init_vel => (something(r.row.init_vel, 0.0), r.row.init_vel = value)
    brake => (run_value(r, :brake, something(r.row.brake, 0.0)),
              set_run_value!(r, :brake, value))
    speed_controlled => (run_value(r, :speed_controlled,
                                   something(r.row.speed_controlled, false)),
                         set_run_value!(r, :speed_controlled, value))
    coulomb_friction => (run_value(r, :coulomb_friction,
                                   winch_setting(r, :coulomb_friction, :f_coulomb)),
                         set_run_value!(r, :coulomb_friction, value))
    viscous_coefficient => (run_value(r, :viscous_coefficient,
                                      winch_setting(r, :viscous_coefficient, :c_vf)),
                            set_run_value!(r, :viscous_coefficient, value))
    inertia_total => (run_value(r, :inertia_total,
                                winch_setting(r, :inertia_total, :inertia_total)),
                      set_run_value!(r, :inertia_total, value))
end

"""The collection of a run's state a record is one of, as a parameter path names it."""
collection(::PointState) = :points
collection(::StationState) = :spanwise_stations
collection(::SegmentState) = :segments
collection(::PulleyState) = :pulleys
collection(::TetherState) = :tethers
collection(::WinchState) = :winches

"""The position slot of `body` in its run's `sysstate.X/Y/Z`."""
body_slot(body) = getfield(getfield(body, :state), :body_slot0) + index(body)

"""The orientation frame of `body` in its run's `sysstate`."""
body_frame(body) = getfield(getfield(body, :state), :n_wings) + index(body)

"""The centre of mass of body record `r` in its body frame, from its origin [m]."""
com_offset(r) = body_geometry(r).com_offset_b

"""The aero mode of body record `r`, with the model's VSM engine where it has one."""
aero_mode(r) = getfield(run_state(r), :aero)[index(r)]

"""The [`BodyGeometry`](@ref) of body record `r`."""
body_geometry(r) = getfield(run_state(r), :body_geometry)[index(r)]

@state_record BodyState """
    BodyState

A body of a run's state, a wing among them: its pose and motion in its body frame, the
same in its principal frame, derived from them, what the air does to it, and the fields of
its aero mode's VSM engine, read through to the model's.
""" begin
    aero => (aero_mode(r), derived(:aero))
    vsm_aoa => (mean(r.vsm_solver.sol.alpha_dist), derived(:vsm_aoa))
    pos_w => ColumnVec3(sysstate(r).X, sysstate(r).Y, sysstate(r).Z, body_slot(r))
    vel_w => ColumnVec3(sysstate(r).VX, sysstate(r).VY, sysstate(r).VZ, body_slot(r))
    Q_b_to_w => ColumnQuat(sysstate(r).Qw, sysstate(r).Qx, sysstate(r).Qy,
                           sysstate(r).Qz, body_frame(r))
    ω_b => ColumnVec3(sysstate(r).turn_rate_x, sysstate(r).turn_rate_y,
                      sysstate(r).turn_rate_z, body_frame(r))
    R_b_to_w => (quaternion_to_rotation_matrix(r.Q_b_to_w),
                 r.Q_b_to_w .= rotation_matrix_to_quaternion(value))
    R_p_to_w => (r.R_b_to_w * body_geometry(r).R_b_to_p', derived(:R_p_to_w))
    com_w => (r.pos_w + r.R_b_to_w * com_offset(r), derived(:com_w))
    com_vel => (r.vel_w + cross(r.R_b_to_w * r.ω_b, r.R_b_to_w * com_offset(r)),
                derived(:com_vel))
    Q_p_to_w => (rotation_matrix_to_quaternion(r.R_p_to_w), derived(:Q_p_to_w))
    ω_p => (body_geometry(r).R_b_to_p * r.ω_b, derived(:ω_p))
    acc_w => ColumnVec3(sysstate(r).acc_x, sysstate(r).acc_y, sysstate(r).acc_z,
                        body_slot(r))
    wind_vec => ColumnVec3(sysstate(r).wind_x, sysstate(r).wind_y, sysstate(r).wind_z,
                           body_slot(r))
    va_b => (body_apparent_wind(r), derived(:va_b))
    principal_frame_method => (body_principal_frame_method(r),
                               derived(:principal_frame_method))
    extra_mass => (something(r.row.extra_mass, 0.0), derived(:extra_mass))
    origin_ENU => (origin_ENU(r), derived(:origin_ENU))
    angular_damping => (run_value(r, :angular_damping, body_angular_damping(r)),
                        set_run_value!(r, :angular_damping, value))
    world_frame_damping => (run_value(r, :world_frame_damping,
                                      per_axis(r.row.world_frame_damping)),
                            set_run_value!(r, :world_frame_damping, value))
    body_frame_damping => (run_value(r, :body_frame_damping,
                                     per_axis(r.row.body_frame_damping)),
                           set_run_value!(r, :body_frame_damping, value))
    ext_force_w => (run_value(r, :ext_force_w, column_vector(r.row, :ext_force_w)),
                    set_run_value!(r, :ext_force_w, value))
    ext_force_b => (run_value(r, :ext_force_b, column_vector(r.row, :ext_force_b)),
                    set_run_value!(r, :ext_force_b, value))
    ext_moment_b => (run_value(r, :ext_moment_b, column_vector(r.row, :ext_moment_b)),
                     set_run_value!(r, :ext_moment_b, value))
    fix_sphere => (run_value(r, :fix_sphere, something(r.row.fix_sphere, false)),
                   set_run_value!(r, :fix_sphere, value))
    fix_static => (run_value(r, :fix_static, something(r.row.fix_static, false)),
                   set_run_value!(r, :fix_static, value))
    apparent_mass => (run_value(r, :apparent_mass, something(r.row.apparent_mass, 0.0)),
                      set_run_value!(r, :apparent_mass, value))
    group_points_moment => (run_value(r, :group_points_moment,
                                      something(wing_column(r, :group_points_moment),
                                                true)),
                            set_run_value!(r, :group_points_moment, value))
    wind_disturb => (run_value(r, :wind_disturb, zero(SVector{3, Float64})),
                     set_run_value!(r, :wind_disturb, value))
    drag_frac => (run_value(r, :drag_frac, 1.0), set_run_value!(r, :drag_frac, value))
    aero_force_b => ColumnVec3(sysstate(r).aero_force_KA_x, sysstate(r).aero_force_KA_y,
                               sysstate(r).aero_force_KA_z, index(r))
    aero_moment_b => ColumnVec3(sysstate(r).aero_moment_KA_x, sysstate(r).aero_moment_KA_y,
                                sysstate(r).aero_moment_KA_z, index(r))
    turn_rate => (wing_rates(r).turn_rate, derived(:turn_rate))
    turn_acc => ColumnVec3(sysstate(r).turn_acc_x, sysstate(r).turn_acc_y,
                           sysstate(r).turn_acc_z, index(r))
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
end aero_mode => VSM_ENGINE_FIELDS body_geometry => BODY_GEOMETRY_FIELDS

"""The apparent wind on `point` in the body frame of its wing [m/s], zero where it has
none."""
function point_apparent_wind(point)
    wing_idx = point.wing_idx
    isnothing(wing_idx) && return zero(SVector{3, SimFloat})
    wing = getfield(point, :state).bodies[wing_idx]
    return SVector{3}(wing.R_b_to_w' * (point.wind_vec - point.vel_w))
end
collection(::BodyState) = :bodies

"""The aerodynamic force on `point` in the body frame of its wing [N], zero where it has
none."""
function point_aero_force_b(point)
    wing_idx = point.wing_idx
    isnothing(wing_idx) && return zero(SVector{3, SimFloat})
    wing = getfield(point, :state).bodies[wing_idx]
    return SVector{3}(wing.R_b_to_w' * point.aero_force)
end

"""The apparent wind on `body` in its body frame [m/s]."""
body_apparent_wind(body) = SVector{3}(body.R_b_to_w' *
    (body.wind_vec - body.vel_w + body.wind_disturb))

"""The mass of `point`, its own and half of each segment on it [kg]."""
function point_total_mass(point)
    state = getfield(point, :state)
    mass = point.extra_mass
    for i in getfield(state, :point_segments)[index(point)]
        segment = state.segments[i]
        mass += segment_half_mass(segment.rest_length, segment.diameter, segment.density)
    end
    return mass
end

"""The force along and moment about the twist axis of `station` from the loads on its
points, zero for a station that does not twist."""
function station_tether_load(station)
    state = getfield(station, :state)
    station_type(station) == KINEMATIC && return (; force = 0.0, moment = 0.0)
    orientation = state.bodies[station.wing].R_b_to_w
    force, moment = 0.0, 0.0
    for point_idx in station.points
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
    winch_point = winch.winch_point
    force = zero(SVector{3, SimFloat})
    for tether_idx in winch.tethers, i in state.tethers[tether_idx].segments
        segment = state.segments[i]
        ends = segment.points
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
    for station in state.spanwise_stations
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
    for station in state.spanwise_stations
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
