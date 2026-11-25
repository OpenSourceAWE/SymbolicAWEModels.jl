# Copyright (c) 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# ==================== QUANTITY KINDS ==================== #

"""
    Columns(names...; slot)

A quantity a run holds in the columns `names` of its `SysState`, at the entry
`slot(sam, i)` gives item `i` (the item's own index by default): one name for a number,
three for a vector, four for a scalar-first quaternion.
"""
struct Columns{Names, S}
    slot::S
end
Columns(names::Symbol...; slot=item_slot) = Columns{names, typeof(slot)}(slot)

"""The entry of item `i` in the columns of its collection: `i` itself."""
item_slot(sam, i) = i

"""The position slot of body `body_idx` in the `SysState` of `sam`."""
body_position_slot(sam, body_idx) = getfield(sam, :topology).body_slot0 + body_idx

"""The orientation frame of body `body_idx` in the `SysState` of `sam`: after the wings'."""
body_frame(sam, body_idx) = getfield(sam, :topology).n_wings + body_idx

"""
    Parameter(default)

A quantity the integrator of a run holds as a parameter once it is built, and
`default(sam, i)` holds before, or where no parameter survived compilation.
"""
struct Parameter{D}
    default::D
end

"""
    Authored(column, fallback)

The default of a [`Parameter`](@ref): the item's row column `column`, else `fallback`.
"""
struct Authored{T}
    column::Symbol
    fallback::T
end

"""
    Derived(value)

A quantity that follows from the rest of a run, `value(sam, i)`, read only.
"""
struct Derived{F}
    value::F
end

# ==================== QUANTITY TABLES ==================== #

"""The run's quantities of each point, by name."""
const POINT_QUANTITIES = (
    pos_w = Columns(:X, :Y, :Z),
    vel_w = Columns(:VX, :VY, :VZ),
    drag_force = Columns(:drag_force_x, :drag_force_y, :drag_force_z),
    force = Columns(:force_x, :force_y, :force_z),
    aero_force = Columns(:aero_force_x, :aero_force_y, :aero_force_z),
    wind_vec = Columns(:wind_x, :wind_y, :wind_z),
    aero_force_b = Derived((sam, i) -> point_aero_force_b(sam, i)),
    va_b = Derived((sam, i) -> point_apparent_wind(sam, i)),
    total_mass = Derived((sam, i) -> point_total_mass(sam, i)),
    apparent_mass = Parameter(Authored(:apparent_mass, 0.0)),
    extra_mass = Parameter(Authored(:extra_mass, 0.0)),
    body_frame_damping = Parameter((sam, i) ->
        per_axis(sam.definition.points[i].body_frame_damping)),
    world_frame_damping = Parameter((sam, i) ->
        per_axis(sam.definition.points[i].world_frame_damping)),
    area = Parameter(Authored(:area, 0.0)),
    drag_coeff = Parameter(Authored(:drag_coeff, 0.0)),
    fix_sphere = Parameter(Authored(:fix_sphere, false)),
    fix_static = Parameter(Authored(:fix_static, false)),
    ext_force_w = Parameter((sam, i) -> zero(SVec3)),
    is_wing_node = Derived((sam, i) -> wing_node(sam.definition, i)),
    wing_idx = Derived((sam, i) -> point_wing(sam, i)),
    anchor_b = Derived((sam, i) -> point_anchor_b(sam, i)),
    pos_undeformed_b = Derived((sam, i) -> point_undeformed_b(sam, i)),
    beam_frac = Derived((sam, i) -> first(point_beam_anchor(sam, i))),
    beam_offset_b = Derived((sam, i) -> last(point_beam_anchor(sam, i))),
)

"""The run's quantities of each spanwise station, by name."""
const STATION_QUANTITIES = (
    twist = Columns(:twist_angles),
    twist_ω = Columns(:twist_vel),
    tether_force = Derived((sam, i) -> station_tether_load(sam, i).force),
    tether_moment = Derived((sam, i) -> station_tether_load(sam, i).moment),
    aero_moment = Derived((sam, i) -> station_aero_moment(sam, i)),
    held_aero_moment = Parameter((sam, i) -> 0.0),
    moment_frac = Parameter(Authored(:moment_frac, 0.0)),
    damping = Parameter(Authored(:damping, 50.0)),
    stiffness = Parameter(Authored(:stiffness, 0.0)),
    flap_axis = Parameter((sam, i) -> SVec3(something(
        sam.definition.spanwise_stations[i].flap_axis, (0.0, 1.0, 0.0)))),
    area = Derived((sam, i) -> something(sam.definition.spanwise_stations[i].area, NaN)),
    le_pos = Derived((sam, i) -> station_le_pos(sam, i)),
    chord = Derived((sam, i) -> station_chord(sam, i)),
    y_airf = Derived((sam, i) -> station_y_airf(sam, i)),
    unrefined_section_idxs = Derived((sam, i) ->
        something(station_vsm(sam, i, :station_sections), Int[])),
    panel_idxs = Derived((sam, i) ->
        something(station_vsm(sam, i, :station_panels), Int[])),
    body_mass = Derived((sam, i) ->
        something(station_vsm(sam, i, :station_body_mass), 0.0)),
    flap_chord_refs = Derived((sam, i) -> flap_chord_refs(sam, i)),
    flap_rest_delta = Derived((sam, i) -> flap_rest_delta(sam, i)),
)

"""The run's quantities of each segment, by name."""
const SEGMENT_QUANTITIES = (
    unit_stiffness = Parameter((sam, i) -> material_unit_stiffness(
        segment_material(sam, i), run_value(sam, Val(:segments), i, Val(:diameter)))),
    density = Parameter((sam, i) -> segment_material(sam, i).parameters["density"]),
    unit_damping = Parameter((sam, i) -> segment_unit_damping(
        sam.definition.segments[i], segment_material(sam, i), sam.set)),
    compression_frac = Parameter(Authored(:compression_frac, 0.1)),
    compression_damping_frac = Parameter(Authored(:compression_damping_frac, 1.0)),
    rest_length = Parameter((sam, i) -> segment_rest_length(sam, i)),
    len = Derived((sam, i) -> segment_world_length(sam, i)),
    force = Columns(:spring_force),
    diameter = Parameter((sam, i) -> segment_diameter(sam, i)),
)

"""The run's quantities of each pulley, by name."""
const PULLEY_QUANTITIES = (
    len = Columns(:pulley_len),
    vel = Columns(:pulley_vel),
    efficiency = Parameter(Authored(:efficiency, 0.95)),
    damping = Parameter(Authored(:damping, 0.0)),
    brake = Parameter(Authored(:brake, false)),
    friction_epsilon = Parameter(Authored(:friction_epsilon, 0.1)),
    sum_len = Derived((sam, i) -> pulley_sum_len(sam, i)),
)

"""The run's quantities of each tether, by name."""
const TETHER_QUANTITIES = (
    len = Columns(:l_tether),
    stretched_len = Derived((sam, i) -> tether_stretched_length(sam, i)),
)

"""The run's quantities of each winch, by name."""
const WINCH_QUANTITIES = (
    vel = Columns(:v_reelout),
    set_value = Columns(:set_torque),
    acc = Columns(:winch_acc),
    friction = Columns(:winch_friction),
    force = Derived((sam, i) -> winch_force_vector(sam, i)),
    init_vel = Derived((sam, i) -> something(sam.definition.winches[i].init_vel, 0.0)),
    brake = Parameter(Authored(:brake, 0.0)),
    speed_controlled = Parameter(Authored(:speed_controlled, false)),
    coulomb_friction = Parameter((sam, i) -> winch_setting(sam, i, :coulomb_friction,
                                                            :f_coulomb)),
    viscous_coefficient = Parameter((sam, i) -> winch_setting(sam, i, :viscous_coefficient,
                                                               :c_vf)),
    inertia_total = Parameter((sam, i) -> winch_setting(sam, i, :inertia_total,
                                                         :inertia_total)),
    drum_radius = Parameter(Authored(:drum_radius, NaN)),
    gear_ratio = Parameter(Authored(:gear_ratio, NaN)),
    model = Derived((sam, i) -> sam.definition.winches[i].model),
)

"""The run's quantities of each body, by name."""
const BODY_QUANTITIES = (
    pos_w = Columns(:X, :Y, :Z; slot=body_position_slot),
    vel_w = Columns(:VX, :VY, :VZ; slot=body_position_slot),
    acc_w = Columns(:acc_x, :acc_y, :acc_z; slot=body_position_slot),
    wind_vec = Columns(:wind_x, :wind_y, :wind_z; slot=body_position_slot),
    Q_b_to_w = Columns(:Qw, :Qx, :Qy, :Qz; slot=body_frame),
    ω_b = Columns(:turn_rate_x, :turn_rate_y, :turn_rate_z; slot=body_frame),
    R_b_to_w = Derived((sam, body_idx) -> body_rotation(sam, body_idx)),
    R_p_to_w = Derived((sam, body_idx) ->
        body_rotation(sam, body_idx) * body_R_b_to_p(sam, body_idx)'),
    Q_p_to_w = Derived((sam, body_idx) -> rotation_matrix_to_quaternion(
        body_rotation(sam, body_idx) * body_R_b_to_p(sam, body_idx)')),
    ω_p = Derived((sam, body_idx) -> body_R_b_to_p(sam, body_idx) *
        run_value(sam, Val(:bodies), body_idx, Val(:ω_b))),
    com_w = Derived((sam, body_idx) -> body_com_w(sam, body_idx)),
    com_vel = Derived((sam, body_idx) -> body_com_vel(sam, body_idx)),
    va_b = Derived((sam, body_idx) -> body_apparent_wind(sam, body_idx)),
    total_mass = Derived((sam, body_idx) -> body_mass_properties(sam, body_idx).total_mass),
    total_com_offset_b = Derived((sam, body_idx) ->
        body_mass_properties(sam, body_idx).com_offset_b),
    inertia_principal = Derived((sam, body_idx) ->
        body_mass_properties(sam, body_idx).inertia_principal),
    R_b_to_p = Derived((sam, body_idx) -> body_R_b_to_p(sam, body_idx)),
    extra_com_offset_b = Derived((sam, body_idx) -> extra_com_offset_b(sam, body_idx)),
    extra_mass = Derived((sam, body_idx) ->
        something(sam.definition.bodies[body_idx].extra_mass, 0.0)),
    origin_ENU = Derived((sam, body_idx) -> origin_ENU(sam, body_idx)),
    dynamics_type = Derived((sam, body_idx) -> body_dynamics_type(sam, body_idx)),
    principal_frame_method = Derived((sam, body_idx) ->
        body_principal_frame_method(sam, body_idx)),
    angular_damping = Parameter((sam, body_idx) -> body_angular_damping(sam, body_idx)),
    world_frame_damping = Parameter((sam, body_idx) ->
        per_axis(sam.definition.bodies[body_idx].world_frame_damping)),
    body_frame_damping = Parameter((sam, body_idx) ->
        per_axis(sam.definition.bodies[body_idx].body_frame_damping)),
    ext_force_w = Parameter((sam, body_idx) ->
        column_vector(sam.definition.bodies[body_idx], :ext_force_w)),
    ext_force_b = Parameter((sam, body_idx) ->
        column_vector(sam.definition.bodies[body_idx], :ext_force_b)),
    ext_moment_b = Parameter((sam, body_idx) ->
        column_vector(sam.definition.bodies[body_idx], :ext_moment_b)),
    fix_sphere = Parameter(Authored(:fix_sphere, false)),
    fix_static = Parameter(Authored(:fix_static, false)),
    apparent_mass = Parameter(Authored(:apparent_mass, 0.0)),
    group_points_moment = Parameter((sam, body_idx) ->
        something(wing_column(sam, body_idx, :group_points_moment), true)),
    wind_disturb = Parameter((sam, body_idx) -> zero(SVec3)),
    drag_frac = Parameter((sam, body_idx) -> 1.0),
)

"""The run's quantities of each wing that are not its body's, by name."""
const WING_QUANTITIES = (
    aero = Derived((sam, wing_idx) -> sam.aero[wing_idx]),
    aero_x = Derived((sam, wing_idx) -> sam.aero[wing_idx].aero_x),
    aero_y = Derived((sam, wing_idx) -> sam.aero[wing_idx].aero_y),
    aero_jac = Derived((sam, wing_idx) -> sam.aero[wing_idx].aero_jac),
    unsteady = Derived((sam, wing_idx) -> sam.aero[wing_idx].unsteady),
    vsm_aoa = Derived((sam, wing_idx) -> mean(circulation_alpha(sam, wing_idx))),
    aero_force_b = Columns(:aero_force_KA_x, :aero_force_KA_y, :aero_force_KA_z),
    aero_moment_b = Columns(:aero_moment_KA_x, :aero_moment_KA_y, :aero_moment_KA_z),
    turn_acc = Columns(:turn_acc_x, :turn_acc_y, :turn_acc_z),
    elevation = Columns(:elevations),
    azimuth = Columns(:azimuths),
    heading = Columns(:headings),
    course = Columns(:courses),
    aoa = Columns(:AoAs),
    turn_rate = Derived((sam, wing_idx) -> wing_rates(sam, wing_idx).turn_rate),
    elevation_vel = Derived((sam, wing_idx) -> wing_rates(sam, wing_idx).elevation_vel),
    elevation_acc = Derived((sam, wing_idx) -> wing_rates(sam, wing_idx).elevation_acc),
    azimuth_vel = Derived((sam, wing_idx) -> wing_rates(sam, wing_idx).azimuth_vel),
    azimuth_acc = Derived((sam, wing_idx) -> wing_rates(sam, wing_idx).azimuth_acc),
    station_idxs = Derived((sam, wing_idx) -> wing_station_idxs(sam, wing_idx)),
    origin = Derived((sam, wing_idx) -> wing_origin(sam, wing_idx)),
    wing_frame = Derived((sam, wing_idx) -> wing_frame(sam, wing_idx)),
    mesh_origin_b = Derived((sam, wing_idx) -> wing_mesh_pose(sam, wing_idx).origin_b),
    R_b_to_mesh = Derived((sam, wing_idx) -> wing_mesh_pose(sam, wing_idx).R_b_to_mesh),
)

"""The run's quantities of each tube, by name."""
const TUBE_QUANTITIES = (
    EA = Parameter((sam, tube_idx) -> tube_rigidity(sam, tube_idx, :EA)),
    GA = Parameter((sam, tube_idx) -> tube_rigidity(sam, tube_idx, :GA)),
    GJ = Parameter((sam, tube_idx) -> tube_rigidity(sam, tube_idx, :GJ)),
    EIy = Parameter((sam, tube_idx) -> tube_rigidity(sam, tube_idx, :EIy)),
    EIz = Parameter((sam, tube_idx) -> tube_rigidity(sam, tube_idx, :EIz)),
    stiffness_axial = Parameter((sam, tube_idx) ->
        tube_rigidity(sam, tube_idx, :stiffness_axial)),
    stiffness_shear = Parameter((sam, tube_idx) ->
        tube_rigidity(sam, tube_idx, :stiffness_shear)),
    stiffness_torsion = Parameter((sam, tube_idx) ->
        tube_rigidity(sam, tube_idx, :stiffness_torsion)),
    stiffness_bending = Parameter((sam, tube_idx) ->
        tube_rigidity(sam, tube_idx, :stiffness_bending)),
    shear_coeff = Parameter((sam, tube_idx) -> tube_column(sam, tube_idx, :shear_coeff,
                                                           5 / 6)),
    damping = Parameter((sam, tube_idx) -> tube_column(sam, tube_idx, :damping, 0.0)),
    anchor_a_b = Derived((sam, tube_idx) -> tube_anchors_b(sam, tube_idx)[1]),
    anchor_b_b = Derived((sam, tube_idx) -> tube_anchors_b(sam, tube_idx)[2]),
    rest_length = Derived((sam, tube_idx) -> tube_rest_value(sam, tube_idx, :rest_length)),
    R_a_rel0 = Derived((sam, tube_idx) -> tube_rest_value(sam, tube_idx, :R_a_rel0)),
    R_b_rel0 = Derived((sam, tube_idx) -> tube_rest_value(sam, tube_idx, :R_b_rel0)),
    rest_offset_a = Derived((sam, tube_idx) ->
        tube_rest_value(sam, tube_idx, :rest_offset_a)),
    R_rel0 = Derived((sam, tube_idx) -> tube_rest_value(sam, tube_idx, :R_rel0)),
)

"""The quantity tables of a run's collections, by collection."""
const QUANTITIES = (points = POINT_QUANTITIES, spanwise_stations = STATION_QUANTITIES,
                    segments = SEGMENT_QUANTITIES, pulleys = PULLEY_QUANTITIES,
                    tethers = TETHER_QUANTITIES, winches = WINCH_QUANTITIES,
                    bodies = BODY_QUANTITIES, wings = WING_QUANTITIES,
                    tubes = TUBE_QUANTITIES)

"""Whether `collection` names a collection of a run's items."""
is_collection(collection) =
    collection isa Symbol && hasfield(typeof(QUANTITIES), collection)

"""Whether the run's items of `collection` have a quantity `quantity`."""
has_quantity(collection::Symbol, quantity::Symbol) =
    hasfield(typeof(getfield(QUANTITIES, collection)), quantity)

# ==================== RUN VALUES ==================== #

"""The index of item `id` of `collection` of `sam`: `id` itself, or the index of the row
named `id`."""
item_index(sam, collection, id::Integer) = Int(id)
item_index(sam, collection, id::Symbol) =
    KiteGeometry.row_index(getproperty(sam.definition, collection), id)

"""
    run_value(sam, collection, id, quantity)

The value of `quantity` of item `id`, an index or a row name, of `collection` of the run
`sam` (`:points`, `:spanwise_stations`, `:segments`, `:pulleys`, `:tethers`, `:winches`,
`:bodies`, `:wings` or `:tubes`), as the tables [`POINT_QUANTITIES`](@ref) and its
siblings define it. A wing's quantity that is not its own is its body's.
"""
run_value(sam, collection::Symbol, id, quantity::Symbol) =
    run_value(sam, Val(collection), item_index(sam, collection, id), Val(quantity))

@inline function run_value(sam, ::Val{C}, i::Int, ::Val{Q}) where {C, Q}
    if C === :wings && !hasfield(typeof(WING_QUANTITIES), Q)
        return run_value(sam, Val(:bodies), wing_body(sam, i), Val(Q))
    end
    hasfield(typeof(getfield(QUANTITIES, C)), Q) ||
        error("$C have no quantity $Q in a run")
    return spec_value(getfield(getfield(QUANTITIES, C), Q), sam, Val(C), i, Val(Q))
end

"""
    set_run_value!(sam, collection, id, quantity, value)

Write `value` to `quantity` of item `id` of `collection` of the run `sam`: into the
columns of its `SysState`, or into its integrator's parameter. Refuses a quantity that
follows from others, and a parameter before `init!` or one compilation pruned.
"""
function set_run_value!(sam, collection::Symbol, id, quantity::Symbol, value)
    set_run_value!(sam, Val(collection), item_index(sam, collection, id), Val(quantity),
                   value)
    return value
end

function set_run_value!(sam, ::Val{C}, i::Int, ::Val{Q}, value) where {C, Q}
    if C === :wings && !hasfield(typeof(WING_QUANTITIES), Q)
        return set_run_value!(sam, Val(:bodies), wing_body(sam, i), Val(Q), value)
    end
    hasfield(typeof(getfield(QUANTITIES, C)), Q) ||
        error("$C have no quantity $Q in a run")
    return set_spec_value!(getfield(getfield(QUANTITIES, C), Q), sam, Val(C), i, Val(Q),
                           value)
end

@inline spec_value(spec::Columns{Names}, sam, ::Val, i, ::Val) where {Names} =
    column_value(getfield(sam, :sysstate), Names, spec.slot(sam, i))
spec_value(spec::Derived, sam, ::Val, i, ::Val) = spec.value(sam, i)
spec_value(spec::Parameter, sam, ::Val{C}, i, ::Val{Q}) where {C, Q} =
    parameter_value(sam, (C, i, Q), () -> default_value(spec.default, sam, Val(C), i))

"""The value `default` gives item `i` of collection `C` before a run holds its own."""
default_value(default, sam, ::Val, i) = default(sam, i)
default_value(default::Authored, sam, ::Val{C}, i) where {C} =
    something(getproperty(getproperty(getfield(sam, :definition), C)[i], default.column),
              default.fallback)

@inline column_value(sysstate, names::NTuple{1, Symbol}, slot) =
    getfield(sysstate, names[1])[slot]
@inline column_value(sysstate, names::NTuple{3, Symbol}, slot) =
    SVector(getfield(sysstate, names[1])[slot], getfield(sysstate, names[2])[slot],
            getfield(sysstate, names[3])[slot])
@inline column_value(sysstate, names::NTuple{4, Symbol}, slot) =
    SVector(getfield(sysstate, names[1])[slot], getfield(sysstate, names[2])[slot],
            getfield(sysstate, names[3])[slot], getfield(sysstate, names[4])[slot])

function set_spec_value!(spec::Columns{Names}, sam, ::Val, i, ::Val, value) where {Names}
    sysstate, slot = getfield(sam, :sysstate), spec.slot(sam, i)
    if length(Names) == 1
        getfield(sysstate, Names[1])[slot] = value
    else
        length(value) == length(Names) ||
            error("$(length(Names)) values expected, got $(length(value))")
        for (k, name) in enumerate(Names)
            getfield(sysstate, name)[slot] = value[k]
        end
    end
    return value
end
set_spec_value!(::Derived, sam, ::Val{C}, i, ::Val{Q}, value) where {C, Q} =
    error("$C $i: $Q follows from the rest of the run; set that instead")
function set_spec_value!(::Parameter, sam, ::Val{C}, i, ::Val{Q}, value) where {C, Q}
    integrator = getfield(sam, :integrator)
    isnothing(integrator) && error("$C $i: $Q is a parameter of a run's integrator, " *
        "which init! builds; set it on the definition before init!")
    access = parameter_access(sam, (C, i, Q))
    isnothing(access) && error("$C $i: $Q is not a parameter of this model")
    write_parameter!(access, integrator, value)
    return value
end

"""
    parameter_value(sam, path, default)

The value the parameter read from `path` has in the integrator of `sam`, else
`default()`: before `init!`, or where no parameter of it survived compilation. A
3-vector comes back immutable.
"""
function parameter_value(sam, path, default)
    integrator = getfield(sam, :integrator)
    isnothing(integrator) && return default()
    access = parameter_access(sam, path)
    isnothing(access) && return default()
    value = read_parameter(access, integrator)
    value isa AbstractVector{<:Real} && length(value) == 3 && return SVector{3}(value)
    return value
end

"""The access to the parameter read from `path` in the integrator of `sam`, built once
per integrator, `nothing` for none."""
parameter_access(sam, path) = get!(getfield(sam, :parameter_access), path) do
    build_parameter_access(sam.prob.param_sync, sam.prob, path)
end

# ==================== TOPOLOGY ==================== #

"""
    run_topology(sys, n_corners)

The lookups of the structure of definition `sys` a run reads its items by: each point's
segments and wing, each segment's tether and pulley with its side of it, each wing's
stations and body, each body's wing, and the slots of each group of `SysState`
positions, `n_corners` of them for the wings' aero log points.
"""
function run_topology(sys, n_corners)
    body_wings = zeros(Int, length(sys.bodies))
    for (wing_idx, wing) in enumerate(sys.wings)
        body_wings[wing.body] = wing_idx
    end
    slots = position_slots(sys, n_corners)
    return (; point_segments = point_segments(sys),
            point_wings = point_wings(sys, body_wings),
            segment_tethers = segment_tethers(sys),
            segment_pulleys = segment_pulleys(sys),
            wing_stations = [station_indices(sys, wing_idx)
                             for wing_idx in eachindex(sys.wings)],
            wing_bodies = [Int(wing.body) for wing in sys.wings], body_wings, slots,
            body_slot0 = first(slots.bodies) - 1,
            n_wings = length(sys.wings))
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
on, else the wing whose body it rides, `body_wings` giving each body's wing."""
function point_wings(sys, body_wings)
    wings = Union{Nothing, Int}[nothing for _ in sys.points]
    for (i, point) in enumerate(sys.points)
        body = point.body
        isnothing(body) || body_wings[body] == 0 || (wings[i] = body_wings[body])
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

"""The wing point `i` of `sam` belongs to, `nothing` for none."""
point_wing(sam, i) = sam.topology.point_wings[i]

"""The wing whose body is body `body_idx` of `sam`, 0 for none."""
body_wing(sam, body_idx) = sam.topology.body_wings[body_idx]

"""The body of wing `wing_idx` of `sam`."""
wing_body(sam, wing_idx) = getfield(sam, :topology).wing_bodies[wing_idx]

# ==================== POSITIONS ==================== #

"""
    ColumnPositions(x, y, z)

Vectors of a run's points read from three of its `SysState` columns, by point index, as
KiteGeometry reads positions.
"""
struct ColumnPositions{V}
    x::V
    y::V
    z::V
end
Base.getindex(positions::ColumnPositions, i) =
    SVector(positions.x[i], positions.y[i], positions.z[i])

"""The world positions of the points of `sam` [m]."""
run_positions(sam) = ColumnPositions(sam.sysstate.X, sam.sysstate.Y, sam.sysstate.Z)

"""The world velocities of the points of `sam` [m/s]."""
run_velocities(sam) = ColumnPositions(sam.sysstate.VX, sam.sysstate.VY, sam.sysstate.VZ)

"""The positions of the points of `sam` in the pose its definition is placed in [m]."""
initial_positions(sam) = KiteGeometry.Positions(sam.definition.points)

# ==================== BODIES ==================== #

"""The [`WingType`](@ref) of body `body_idx` of `sam`: its wing's, `RIGID_DYNAMICS` for a
body that is no wing."""
function body_dynamics_type(sam, body_idx)
    wing_idx = body_wing(sam, body_idx)
    wing_idx == 0 && return RIGID_DYNAMICS
    return wing_dynamics_type(sam.definition.wings[wing_idx])
end

"""The dynamics type of body `body_idx` of `sam`: `STATIC` where it is held, `KINEMATIC` for
the body of a particle wing and `DYNAMIC` for another."""
function body_type(sam, body_idx)
    sam.definition.bodies[body_idx].fixed && return STATIC
    return body_dynamics_type(sam, body_idx) == PARTICLE_DYNAMICS ? KINEMATIC : DYNAMIC
end

"""Column `column` of the model of body `body_idx` of `sam`, else the same column of its
wing's `wings` row, `nothing` where neither gives it."""
function wing_column(sam, body_idx, column)
    model = sam.definition.bodies[body_idx].model
    value = hasproperty(model, column) ? getproperty(model, column) : nothing
    wing_idx = body_wing(sam, body_idx)
    (isnothing(value) && wing_idx != 0) || return value
    return getproperty(sam.definition.wings[wing_idx].model, column)
end

"""How the principal frame of body `body_idx` of `sam` is found: its own column, else its
wing's, else `EIGEN_DECOMP`."""
body_principal_frame_method(sam, body_idx) = parse_principal_frame_method(
    something(wing_column(sam, body_idx, :principal_frame_method), "EIGEN_DECOMP"))

"""The spin damping of body `body_idx` of `sam` per axis [1/s]: its own column, else its
wing's, else none on a body that is no wing and 150 about y on a wing."""
body_angular_damping(sam, body_idx) =
    per_axis(something(wing_column(sam, body_idx, :angular_damping),
                       body_wing(sam, body_idx) == 0 ? 0.0 : [0.0, 150.0, 0.0]))

"""The body→world rotation of body `body_idx` of `sam` in its run."""
body_rotation(sam, body_idx) =
    quaternion_to_rotation_matrix(run_value(sam, Val(:bodies), body_idx, Val(:Q_b_to_w)))

"""The rotation from the body frame of body `body_idx` of `sam` into its principal frame."""
body_R_b_to_p(sam, body_idx) = body_mass_properties(sam, body_idx).R_b_to_p

"""The centre of mass of body `body_idx` of `sam` in the world [m]."""
body_com_w(sam, body_idx) = run_value(sam, Val(:bodies), body_idx, Val(:pos_w)) +
    body_rotation(sam, body_idx) * body_mass_properties(sam, body_idx).com_offset_b

"""The velocity of the centre of mass of body `body_idx` of `sam` in the world [m/s]."""
function body_com_vel(sam, body_idx)
    R = body_rotation(sam, body_idx)
    return run_value(sam, Val(:bodies), body_idx, Val(:vel_w)) +
        cross(R * run_value(sam, Val(:bodies), body_idx, Val(:ω_b)),
              R * body_mass_properties(sam, body_idx).com_offset_b)
end

"""The apparent wind on body `body_idx` of `sam` in its body frame [m/s]."""
body_apparent_wind(sam, body_idx) = SVector{3}(body_rotation(sam, body_idx)' *
    (run_value(sam, Val(:bodies), body_idx, Val(:wind_vec)) -
     run_value(sam, Val(:bodies), body_idx, Val(:vel_w)) +
     run_value(sam, Val(:bodies), body_idx, Val(:wind_disturb))))

"""The KiteGeometry `WingSpec` of wing `wing_idx` of `sam`: its origin and frame reference
points."""
wing_spec(sam, wing_idx) = KiteGeometry.wing_spec(sam.definition.wings[wing_idx].model)

"""Reference points of the origin of wing `wing_idx` of `sam`, pairs of point index and
weight, `nothing` for its body's position."""
wing_origin(sam, wing_idx) = wing_spec(sam, wing_idx).origin

"""Reference points of the origin and axes of wing `wing_idx` of `sam`, KiteGeometry's
`WingFrame`, `nothing` for none."""
wing_frame(sam, wing_idx) = wing_spec(sam, wing_idx).frame

"""The stations wing `wing_idx` of definition `sys` lists, from +y to -y."""
station_indices(sys, wing_idx) =
    Vector{Int}(something(sys.wings[wing_idx].model.stations, Int[]))

"""The stations wing `wing_idx` of `sam` lists, from +y to -y."""
wing_station_idxs(sam, wing_idx) = sam.topology.wing_stations[wing_idx]

"""The origin of the body frame of wing `wing_idx` of `sam` in its initial pose [m]: its
reference points, else its body's position."""
function wing_origin_ENU(sam, wing_idx)
    origin = wing_origin(sam, wing_idx)
    body = sam.definition.bodies[wing_body(sam, wing_idx)]
    isnothing(origin) && return KVec3(body.pos_ENU)
    return KVec3(KiteGeometry.ref_position(initial_positions(sam), origin))
end

"""The centre of the own mass of body `body_idx` of `sam` from the origin of its frame, body
frame [m]: a wing body's position from its wing's origin, else its `com_offset_b`."""
function extra_com_offset_b(sam, body_idx)
    row = sam.definition.bodies[body_idx]
    wing_idx = body_wing(sam, body_idx)
    wing_idx == 0 && return KVec3(own_com_offset(row))
    return KVec3(initial_rotation(row)' * (row.pos_ENU - wing_origin_ENU(sam, wing_idx)))
end

"""
    origin_ENU(sam, body_idx)

The origin of the frame of body `body_idx` of `sam` in the world [m] in its initial pose:
its `pos_ENU`, the centre of its own mass, less that centre's offset from the origin.
"""
origin_ENU(sam, body_idx) = KVec3(sam.definition.bodies[body_idx].pos_ENU) -
    initial_rotation(sam.definition.bodies[body_idx]) * extra_com_offset_b(sam, body_idx)

"""Where the position `pos_ENU` [m] sits in the frame of body `body_idx` of `sam`, from its
origin, in its initial pose."""
initial_body_position(sam, body_idx, pos_ENU) =
    initial_rotation(sam.definition.bodies[body_idx])' *
    (pos_ENU - origin_ENU(sam, body_idx))

"""The origin and axes of the body of wing `wing_idx` of `sam` in the frame its mesh is
authored in: `(; origin_b, R_b_to_mesh)`."""
function wing_mesh_pose(sam, wing_idx)
    wing = sam.definition.wings[wing_idx]
    R_mesh_to_w = quaternion_to_rotation_matrix(wing.Q_KA2ENU)
    R_b_to_w = initial_rotation(sam.definition.bodies[wing_body(sam, wing_idx)])
    origin = wing_origin_ENU(sam, wing_idx)
    return (; origin_b = KVec3(R_mesh_to_w' * (origin - wing.pos_ENU)),
            R_b_to_mesh = R_mesh_to_w' * R_b_to_w)
end

# ==================== MASS ==================== #

"""
    carrier_body_idx(sam, i) -> Int

The body whose `total_mass` holds point `i`'s mass: the DYNAMIC or STATIC body it rides,
or the body of the `RIGID_DYNAMICS` wing whose node it is. 0 when no integrated body
carries it, and the point then carries its own weight.
"""
function carrier_body_idx(sam, i)
    body = sam.definition.points[i].body
    !isnothing(body) && body_type(sam, body) != KINEMATIC && return Int(body)
    wing_idx = point_wing(sam, i)
    (isnothing(wing_idx) || !wing_node(sam.definition, i)) && return 0
    body_idx = wing_body(sam, wing_idx)
    return body_type(sam, body_idx) == KINEMATIC ? 0 : body_idx
end

"""Where point `i` of `sam`, carried by body `body_idx`, sits in its frame [m]: a rider's
`anchor_b`, else where its initial position puts it."""
function carried_position_b(sam, i, body_idx)
    sam.definition.points[i].body == body_idx && return point_anchor_b(sam, i)
    return initial_body_position(sam, body_idx, sam.definition.points[i].pos_ENU)
end

"""
    body_mass_properties(sam, body_idx)

The mass of body `body_idx` of `sam` with the points it carries ([`carrier_body_idx`](@ref))
as point masses of their `total_mass` at their anchors: `(; total_mass, com_offset_b,
inertia_principal, R_b_to_p)`, its centre of mass from its origin and its principal inertia
about it. The body of a particle wing has the mass of the wing's free points and of the
bodies sharing a tube with its structure ([`particle_wing_parts`](@ref)).
"""
function body_mass_properties(sam, body_idx)
    if body_type(sam, body_idx) == KINEMATIC
        wing_idx = body_wing(sam, body_idx)
        total_mass = wing_idx == 0 ? 0.0 : particle_wing_mass(sam, wing_idx)
        return (; total_mass, com_offset_b = zero(SVec3), inertia_principal = zero(SVec3),
                R_b_to_p = Matrix{SimFloat}(I, 3, 3))
    end
    row = sam.definition.bodies[body_idx]
    extra_mass = something(row.extra_mass, 0.0)
    extra_com = SVec3(extra_com_offset_b(sam, body_idx))
    carried = [(point_total_mass(sam, i), SVec3(carried_position_b(sam, i, body_idx)))
               for i in eachindex(sam.definition.points)
               if carrier_body_idx(sam, i) == body_idx]
    total_mass = extra_mass + sum(first, carried; init=0.0)
    first_moment = extra_mass * extra_com + sum(m * r for (m, r) in carried;
                                                init=zero(SVec3))
    com = total_mass > 0 ? first_moment / total_mass : extra_com
    inertia = something(row.extra_inertia_KA, zeros(3, 3)) .+
        point_mass_inertia(extra_mass, extra_com - com)
    for (mass, position) in carried
        inertia = inertia .+ point_mass_inertia(mass, position - com)
    end
    inertia_principal, R_b_to_p = principal_axes(inertia,
                                                 body_principal_frame_method(sam, body_idx))
    return (; total_mass, com_offset_b = com, inertia_principal = SVec3(inertia_principal),
            R_b_to_p = Matrix{SimFloat}(R_b_to_p))
end

"""The mass of particle wing `wing_idx` of `sam` [kg]: its member points no body carries and
the bodies sharing a tube with its structure."""
function particle_wing_mass(sam, wing_idx)
    point_idxs, body_idxs = particle_wing_parts(sam, wing_idx)
    return sum(point_total_mass(sam, i) for i in point_idxs
               if carrier_body_idx(sam, i) == 0; init=0.0) +
        sum(body_mass_properties(sam, body_idx).total_mass for body_idx in body_idxs;
            init=0.0)
end

"""
    particle_wing_parts(sam, wing_idx) -> (point_idxs, body_idxs)

The member points of particle wing `wing_idx` of `sam` and the bodies of no wing sharing a
tube (KiteGeometry's `connected_body_groups`) with its structure, which carry its mass. The
groups are seeded from its stations' member and flap bodies and the bodies its member points
ride. Members are its stations' points, else the points of its frame
([`wing_frame_member`](@ref)).
"""
function particle_wing_parts(sam, wing_idx)
    sys = sam.definition
    root = KiteGeometry.connected_body_groups(length(sys.bodies), sys.tubes)
    own_body = wing_body(sam, wing_idx)
    seed_roots, point_idxs = Set{Int}(), Set{Int}()
    for station_idx in wing_station_idxs(sam, wing_idx)
        station = sys.spanwise_stations[station_idx]
        union!(point_idxs, station.points)
        for body_idx in Iterators.flatten((something(station.bodies, ()),
                                    something(station.flap_bodies, ())))
            body_idx != 0 && body_idx != own_body && push!(seed_roots, root[body_idx])
        end
    end
    if isempty(point_idxs)
        for i in eachindex(sys.points)
            wing_frame_member(sam, i, wing_idx) && push!(point_idxs, i)
        end
    end
    for i in point_idxs
        body_idx = sys.points[i].body
        !isnothing(body_idx) && body_idx != own_body && push!(seed_roots, root[body_idx])
    end
    body_idxs = [body_idx for body_idx in eachindex(sys.bodies)
                 if body_wing(sam, body_idx) == 0 && root[body_idx] in seed_roots]
    return point_idxs, body_idxs
end

"""Whether point `i` of `sam` belongs to the frame of wing `wing_idx`: a node of its surface
or a `BODY_STATIC` point of the wing."""
wing_frame_member(sam, i, wing_idx) = (wing_node(sam.definition, i) ||
    point_type(sam.definition.points[i]) == BODY_STATIC) && point_wing(sam, i) == wing_idx

# ==================== POINTS ==================== #

"""The mass of point `i` of `sam`, its own and half of each segment on it [kg]."""
function point_total_mass(sam, i)
    mass = run_value(sam, Val(:points), i, Val(:extra_mass))
    for k in sam.topology.point_segments[i]
        mass += segment_half_mass(run_value(sam, Val(:segments), k, Val(:rest_length)),
                                  run_value(sam, Val(:segments), k, Val(:diameter)),
                                  run_value(sam, Val(:segments), k, Val(:density)))
    end
    return mass
end

"""The rotation of the body of the wing point `i` of `sam` belongs to, `nothing` for a
point of no wing."""
function point_wing_rotation(sam, i)
    wing_idx = point_wing(sam, i)
    return isnothing(wing_idx) ? nothing : body_rotation(sam, wing_body(sam, wing_idx))
end

"""The apparent wind on point `i` of `sam` in the body frame of its wing [m/s], zero
where it has none."""
function point_apparent_wind(sam, i)
    R = point_wing_rotation(sam, i)
    isnothing(R) && return zero(SVec3)
    return SVec3(R' * (run_value(sam, Val(:points), i, Val(:wind_vec)) -
                       run_value(sam, Val(:points), i, Val(:vel_w))))
end

"""The aerodynamic force on point `i` of `sam` in the body frame of its wing [N], zero
where it has none."""
function point_aero_force_b(sam, i)
    R = point_wing_rotation(sam, i)
    isnothing(R) && return zero(SVec3)
    return SVec3(R' * run_value(sam, Val(:points), i, Val(:aero_force)))
end

"""Where point `i` of `sam` sits in the frame of the body it rides [m]: its `anchor_b`,
else where its initial position puts it; zero for a point on no body."""
function point_anchor_b(sam, i)
    row = sam.definition.points[i]
    isnothing(row.anchor_b) || return KVec3(row.anchor_b)
    isnothing(row.body) && return zero(KVec3)
    return KVec3(initial_body_position(sam, row.body, row.pos_ENU))
end

"""
    point_undeformed_b(sam, i)

Where wing node `i` of `sam` sits relative to its wing in the initial pose, wing frame
[m]: from a rigid wing's centre of mass, or from a particle wing's frame fitted to the
initial positions of its reference points; zero for a point that is no wing node.
"""
function point_undeformed_b(sam, i)
    wing_idx = point_wing(sam, i)
    (isnothing(wing_idx) || !wing_node(sam.definition, i)) && return zero(KVec3)
    body_idx = wing_body(sam, wing_idx)
    pos = SVec3(sam.definition.points[i].pos_ENU)
    if body_dynamics_type(sam, body_idx) == RIGID_DYNAMICS
        R = initial_rotation(sam.definition.bodies[body_idx])
        com = body_mass_properties(sam, body_idx).com_offset_b
        return KVec3(R' * (pos - (origin_ENU(sam, body_idx) + R * com)))
    end
    frame = wing_frame(sam, wing_idx)
    isnothing(frame) && return zero(KVec3)
    R, origin = KiteGeometry.wing_frame(initial_positions(sam), frame)
    return KVec3(R' * (pos - origin))
end

"""The position along the beam element of the tube point `i` of `sam` rides and its
offset off the centreline ([`beam_anchor`](@ref)), `(0, 0)` for a point on no tube."""
function point_beam_anchor(sam, i)
    row = sam.definition.points[i]
    isnothing(row.tube) && return (0.0, zero(KVec3))
    return beam_anchor(sam, row.pos_ENU, row.tube)
end

# ==================== SEGMENTS, PULLEYS, TETHERS, WINCHES ==================== #

"""The material segment `i` of `sam` names."""
segment_material(sam, i) = sam.definition.materials[sam.definition.segments[i].material]

"""The diameter of segment `i` of `sam` [m], which its row must give."""
function segment_diameter(sam, i)
    segment = sam.definition.segments[i]
    isnothing(segment.diameter) && error("segment $(segment.name) gives no diameter")
    return segment.diameter
end

"""The stiffness per unit length [N] of segment `i` of `sam` whose stiffness is a number,
typed `Float64`."""
numeric_unit_stiffness(sam, i) =
    to_float64(run_value(sam, Val(:segments), i, Val(:unit_stiffness)))

"""The rest length of segment `i` of `sam` [m]: its tether's unstretched length shared
over its segments, its side of the rope over its pulley, else its row's."""
function segment_rest_length(sam, i)
    tether = sam.topology.segment_tethers[i]
    if tether > 0
        segments = sam.definition.tethers[tether].segments
        return sam.sysstate.l_tether[tether] / length(segments)
    end
    pulley, side = sam.topology.segment_pulleys[i]
    if pulley > 0
        len = sam.sysstate.pulley_len[pulley]
        return side == 1 ? len : pulley_sum_len(sam, pulley) - len
    end
    return something(sam.definition.segments[i].rest_length, 0.0)
end

"""The length of segment `i` of `sam` between its points in its run [m]."""
function segment_world_length(sam, i)
    first_point, second_point = sam.definition.segments[i].points
    return norm(run_value(sam, Val(:points), first_point, Val(:pos_w)) -
                run_value(sam, Val(:points), second_point, Val(:pos_w)))
end

"""The length of the rope over pulley `i` of `sam` [m]: its two segments' rest lengths in
the definition."""
pulley_sum_len(sam, i) = sum(something(sam.definition.segments[k].rest_length, 0.0)
                             for k in sam.definition.pulleys[i].segments)

"""The length of the segments of tether `i` of `sam` end to end in its run [m]."""
tether_stretched_length(sam, i) =
    sum(segment_world_length(sam, k) for k in sam.definition.tethers[i].segments;
        init=0.0)

"""Winch `i` of `sam`'s column `column`, else the run's setting `setting`, else `NaN`."""
function winch_setting(sam, i, column, setting)
    value = getproperty(sam.definition.winches[i], column)
    isnothing(value) || return value
    return isnothing(sam.set) ? NaN : getproperty(sam.set, setting)
end

"""The spring force of the tethers of winch `i` of `sam` on its point [N]: the pull on
the winch point, negated."""
function winch_force_vector(sam, i)
    winch = sam.definition.winches[i]
    winch_point = winch.winch_point
    force = zero(SVec3)
    for tether in winch.tethers, k in sam.definition.tethers[tether].segments
        ends = sam.definition.segments[k].points
        winch_point in ends || continue
        other = ends[1] == winch_point ? ends[2] : ends[1]
        away = run_value(sam, Val(:points), other, Val(:pos_w)) -
            run_value(sam, Val(:points), winch_point, Val(:pos_w))
        force -= run_value(sam, Val(:segments), k, Val(:force)) * away / norm(away)
    end
    return force
end

# ==================== STATIONS ==================== #

"""The wing of station `station_idx` of `sam` and the station's position among that wing's
stations."""
function station_place(sam, station_idx)
    wing_idx = sam.definition.spanwise_stations[station_idx].wing
    stations = wing_station_idxs(sam, wing_idx)
    return wing_idx, something(findfirst(==(station_idx), stations), 0)
end

"""What the aero mode of the wing of station `station_idx` of `sam` worked out of it,
`field` of its VSM engine ([`station_vsm_value`](@ref)), `nothing` where it worked out
none."""
function station_vsm(sam, station_idx, field)
    wing_idx, k = station_place(sam, station_idx)
    k == 0 && return nothing
    return station_vsm_value(sam.aero[wing_idx], field, k)
end

"""The leading edge of station `station_idx` of `sam` from the centre of mass of its wing's
body, body frame [m]; zero on a particle wing."""
function station_le_pos(sam, station_idx)
    wing_idx = sam.definition.spanwise_stations[station_idx].wing
    body_idx = wing_body(sam, wing_idx)
    body_type(sam, body_idx) == KINEMATIC && return zero(KVec3)
    le_b = something(station_vsm(sam, station_idx, :station_le_b), zero(KVec3))
    return KVec3(le_b - body_mass_properties(sam, body_idx).com_offset_b)
end

"""The chord vector of station `station_idx` of `sam`, body frame [m]: its row's, else its
VSM section's, else zero."""
function station_chord(sam, station_idx)
    chord = sam.definition.spanwise_stations[station_idx].chord
    isnothing(chord) || return KVec3(chord)
    return KVec3(something(station_vsm(sam, station_idx, :station_chord), zero(KVec3)))
end

"""The spanwise vector of station `station_idx` of `sam`: its row's, else its VSM section's,
else zero."""
function station_y_airf(sam, station_idx)
    y_airf = sam.definition.spanwise_stations[station_idx].y_airf
    isnothing(y_airf) || return KVec3(y_airf)
    return KVec3(something(station_vsm(sam, station_idx, :station_y_airf), zero(KVec3)))
end

"""The force along and moment about the twist axis of station `station_idx` of `sam` from
the loads on its points, zero for a station that does not twist."""
function station_tether_load(sam, station_idx)
    station = sam.definition.spanwise_stations[station_idx]
    station_type(station) == KINEMATIC && return (; force = 0.0, moment = 0.0)
    orientation = body_rotation(sam, wing_body(sam, station.wing))
    surface = (; chord = station_chord(sam, station_idx),
               y_airf = station_y_airf(sam, station_idx),
               le_pos = station_le_pos(sam, station_idx),
               moment_frac = run_value(sam, Val(:spanwise_stations), station_idx,
                                       Val(:moment_frac)))
    twist = run_value(sam, Val(:spanwise_stations), station_idx, Val(:twist))
    force, moment = 0.0, 0.0
    for i in station.points
        couple = twist_bridle_couple(surface, point_undeformed_b(sam, i), twist,
                                     orientation)
        point_force = dot(run_value(sam, Val(:points), i, Val(:force)), couple.direction)
        force += point_force
        moment += couple.arm * point_force
    end
    return (; force, moment)
end

"""The aerodynamic twist moment on station `station_idx` of `sam` [N·m] at its integrator's
state, `0.0` before `init!`."""
function station_aero_moment(sam, station_idx)
    integrator = getfield(sam, :integrator)
    isnothing(integrator) && return 0.0
    return read_station_aero_moment(sam.prob.get_all_state, sam.prob, integrator,
                                    station_idx)
end

"""The aerodynamic twist moment on station `station_idx` [N·m] at `integrator`'s state, from
the monolithic system of `prob`."""
read_station_aero_moment(_, prob, integrator, station_idx) =
    getu(prob.sys, prob.sys.station_aero_moment[station_idx])(integrator)

"""The reference chord directions `[main, flap]` of the flap of station `station_idx` of
`sam`, each body's frame: its row's, else each body's x axis."""
function flap_chord_refs(sam, station_idx)
    refs = sam.definition.spanwise_stations[station_idx].flap_chord_refs
    isnothing(refs) || return KVec3.(refs)
    return [KVec3(1.0, 0.0, 0.0), KVec3(1.0, 0.0, 0.0)]
end

"""
    flap_rest_delta(sam, station_idx)

The deflection of the flap of station `station_idx` of `sam` [rad] in the initial pose: what
the live deflection is measured from. Zero for a station without a flap.
"""
function flap_rest_delta(sam, station_idx)
    station = sam.definition.spanwise_stations[station_idx]
    if has_point_flap(station)
        positions = initial_positions(sam)
        fore, hinge, aft = station.flap_points
        return point_flap_angle(sam, station_idx, positions[fore], positions[hinge],
                                positions[aft],
                                initial_rotation(sam.definition.wings[station.wing]))
    end
    has_body_flap(station) || return 0.0
    bodies = sam.definition.bodies
    return body_flap_angle(sam, station_idx,
                           initial_rotation(bodies[station.flap_bodies[1]]),
                           initial_rotation(bodies[station.flap_bodies[2]]))
end

# ==================== TUBES ==================== #

"""
    tube_columns(tube)

The columns of the element the KiteGeometry `tubes` row `tube` is simulated as: its
model's, else those of a [`TimoshenkoTubeModel`](@ref) with every column unset.
"""
tube_columns(tube) =
    tube.model isa KiteGeometry.NoModel ? TimoshenkoTubeModel() : tube.model

"""
    tube_column(sam, tube_idx, column, fallback)

Column `column` of tube `tube_idx` of `sam`, else `fallback`; `nothing` for a column its
kind of element does not have.
"""
function tube_column(sam, tube_idx, column, fallback)
    columns = tube_columns(sam.definition.tubes[tube_idx])
    hasfield(typeof(columns), column) || return nothing
    return something(getfield(columns, column), fallback)
end

"""
    tube_rigidity(sam, tube_idx, column)

Rigidity `column` of tube `tube_idx` of `sam`, a number or a callable: its row's, else the
one of a Timoshenko tube's material at its diameter and pressure
([`tube_law_rigidities`](@ref)); `nothing` for a column its kind of element does not have.
"""
function tube_rigidity(sam, tube_idx, column)
    tube = sam.definition.tubes[tube_idx]
    columns = tube_columns(tube)
    hasfield(typeof(columns), column) || return nothing
    value = getfield(columns, column)
    isnothing(value) && return law_rigidity(columns, sam, tube_idx, column)
    return value isa Real ? SimFloat(value) : value
end

"""
    law_rigidity(columns, sam, tube_idx, column)

Rigidity `column` of tube `tube_idx` of `sam`, whose element has the row model `columns`,
from its material: an elastic tube has none, so it errors.
"""
function law_rigidity(::TimoshenkoTubeModel, sam, tube_idx, column)
    tube = sam.definition.tubes[tube_idx]
    EA, GA, EI, GJ = tube_law_rigidities(sam.definition.materials[tube.material],
                                         tube.diameter, tube.pressure)
    return getfield((; EA, GA, GJ, EIy = EI, EIz = EI), column)
end
law_rigidity(::ElasticTubeModel, sam, tube_idx, column) =
    error("tube $(sam.definition.tubes[tube_idx].name): $column is unset; an elastic " *
          "tube takes no stiffness from its material")

"""The ends of tube `tube_idx` of `sam` on its two bodies, each from its body's frame
origin, body frame [m]."""
function tube_anchors_b(sam, tube_idx)
    columns = tube_columns(sam.definition.tubes[tube_idx])
    return (KVec3(something(columns.anchor_a, zeros(3))),
            KVec3(something(columns.anchor_b, zeros(3))))
end

"""
    tube_initial_frames(sam, tube_idx) -> (R_a, R_b, anchor_a_w, anchor_b_w)

World rotations of the two bodies of tube `tube_idx` of `sam` and the world positions of its
two anchors, in the bodies' initial poses.
"""
function tube_initial_frames(sam, tube_idx)
    body_a, body_b = sam.definition.tubes[tube_idx].bodies
    R_a = initial_rotation(sam.definition.bodies[body_a])
    R_b = initial_rotation(sam.definition.bodies[body_b])
    anchor_a, anchor_b = tube_anchors_b(sam, tube_idx)
    return R_a, R_b, origin_ENU(sam, body_a) .+ R_a * anchor_a,
        origin_ENU(sam, body_b) .+ R_b * anchor_b
end

"""
    tube_rest(sam, tube_idx) -> NamedTuple

The rest state of tube `tube_idx` of `sam`, in which it is unstrained, from its bodies'
initial poses: for an elastic tube its anchor offset in body A's frame `rest_offset_a`
[m] and the relative rotation `R_rel0 = R_a' R_b`; for a Timoshenko tube its length
`rest_length` [m] and each node's orientation in the element frame, `R_a_rel0` and
`R_b_rel0`.
"""
tube_rest(sam, tube_idx) = tube_rest(tube_columns(sam.definition.tubes[tube_idx]),
                                     tube_initial_frames(sam, tube_idx)...)

function tube_rest(::ElasticTubeModel, R_a, R_b, anchor_a_w, anchor_b_w)
    return (; rest_offset_a = KVec3(R_a' * (anchor_b_w .- anchor_a_w)),
            R_rel0 = Matrix{SimFloat}(R_a' * R_b))
end

function tube_rest(::TimoshenkoTubeModel, R_a, R_b, x_a, x_b)
    e1, e2, e3, len = timoshenko_element_frame(x_a, x_b, R_a)
    element_frame = [e1 e2 e3]
    return (; rest_length = len, R_a_rel0 = Matrix{SimFloat}(element_frame' * R_a),
            R_b_rel0 = Matrix{SimFloat}(element_frame' * R_b))
end

"""Rest quantity `name` of tube `tube_idx` of `sam` ([`tube_rest`](@ref)); `nothing` for one
its kind of element does not have."""
tube_rest_value(sam, tube_idx, name) = get(tube_rest(sam, tube_idx), name, nothing)

"""
    beam_anchor(sam, pos_ENU, tube_idx) -> (beam_frac, beam_offset_b)

Where a point at `pos_ENU` [m] rides the beam element of tube `tube_idx` of `sam` in the
bodies' initial poses: its position `frac ∈ [0, 1]` along the element and its offset off the
centreline in the rest element frame [m].
"""
function beam_anchor(sam, pos_ENU, tube_idx)
    R_a, _, node_a, node_b = tube_initial_frames(sam, tube_idx)
    e1, e2, e3, len = timoshenko_element_frame(node_a, node_b, R_a)
    rel = pos_ENU .- node_a
    frac = clamp(dot(rel, e1) / len, 0.0, 1.0)
    perp = rel .- (frac * len) .* e1
    return frac, KVec3(dot(perp, e1), dot(perp, e2), dot(perp, e3))
end

# ==================== WINGS ==================== #

"""The turn rate and the rates of elevation and azimuth of wing `wing_idx` of `sam`, from
its body's pose and motion and its base point."""
function wing_rates(sam, wing_idx)
    body_idx = wing_body(sam, wing_idx)
    body(quantity) = collect(run_value(sam, Val(:bodies), body_idx, Val(quantity)))
    base_point = transform_base_point(sam.transforms, body_idx)
    rel_pos = base_point == 0 ? body(:pos_w) :
        body(:pos_w) .- collect(run_value(sam, Val(:points), base_point, Val(:pos_w)))
    R_b_to_w = collect(body_rotation(sam, body_idx))
    e_x = R_b_to_w[:, 1]
    return wing_scalar_kinematics(; rel_pos, e_x, R_t_to_w = sym_calc_R_t_to_w(rel_pos),
        R_v_to_w = calc_R_v_to_w(rel_pos, e_x), R_b_to_w, vel = body(:vel_w),
        acc = body(:acc_w), omega_b = body(:ω_b), alpha_b = zeros(3),
        va_b = collect(body_apparent_wind(sam, body_idx)), twist_offset = 0.0)
end

# ==================== STATE VECTOR ==================== #

"""The run's integrated state as one column: every point's position and velocity, every
body's pose and motion in its body frame, every station's twist and its rate, every
pulley's split and its rate, every tether's length and every winch's reel-out speed."""
function state_vars(sam)
    sys = sam.definition
    value(collection, i, quantity) = run_value(sam, Val(collection), i, Val(quantity))
    vars = SimFloat[]
    for i in eachindex(sys.points)
        append!(vars, value(:points, i, :pos_w), value(:points, i, :vel_w))
    end
    for body_idx in eachindex(sys.bodies)
        append!(vars, value(:bodies, body_idx, :pos_w), value(:bodies, body_idx, :vel_w),
                value(:bodies, body_idx, :Q_b_to_w), value(:bodies, body_idx, :ω_b))
    end
    for station_idx in eachindex(sys.spanwise_stations)
        push!(vars, value(:spanwise_stations, station_idx, :twist),
              value(:spanwise_stations, station_idx, :twist_ω))
    end
    for i in eachindex(sys.pulleys)
        push!(vars, value(:pulleys, i, :len), value(:pulleys, i, :vel))
    end
    append!(vars, (value(:tethers, i, :len) for i in eachindex(sys.tethers)))
    append!(vars, (value(:winches, i, :vel) for i in eachindex(sys.winches)))
    return reshape(vars, :, 1)
end

"""Write `vars`, laid out as [`state_vars`](@ref) gives it, into the run `sam`."""
function set_state_vars!(sam, vars)
    sys = sam.definition
    values = Iterators.Stateful(vec(vars))
    put!(collection, i, quantity, n) =
        set_run_value!(sam, Val(collection), i, Val(quantity),
                       n == 1 ? popfirst!(values) : collect(Iterators.take(values, n)))
    for i in eachindex(sys.points)
        put!(:points, i, :pos_w, 3)
        put!(:points, i, :vel_w, 3)
    end
    for body_idx in eachindex(sys.bodies)
        put!(:bodies, body_idx, :pos_w, 3)
        put!(:bodies, body_idx, :vel_w, 3)
        put!(:bodies, body_idx, :Q_b_to_w, 4)
        put!(:bodies, body_idx, :ω_b, 3)
    end
    for station_idx in eachindex(sys.spanwise_stations)
        put!(:spanwise_stations, station_idx, :twist, 1)
        put!(:spanwise_stations, station_idx, :twist_ω, 1)
    end
    for i in eachindex(sys.pulleys)
        put!(:pulleys, i, :len, 1)
        put!(:pulleys, i, :vel, 1)
    end
    foreach(i -> put!(:tethers, i, :len, 1), eachindex(sys.tethers))
    foreach(i -> put!(:winches, i, :vel, 1), eachindex(sys.winches))
    return sam
end
