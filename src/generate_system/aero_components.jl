# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Aero coupling components (winch-style swappable subsystems).
#
# A wing carries an `aero_model` builder selected by `aero_mode`
# (see `resolve_aero_model`). Each builder returns a `System` whose
# connectors are fixed by the wing's `dynamics_type`:
#
#   RIGID_DYNAMICS (num_groups = length(wing.group_idxs)):
#     inputs:  va[1:3], rho, R_b_w[1:3,1:3], omega[1:3],
#              twist[1:num_groups], twist_vel[1:num_groups]
#     outputs: force[1:3], moment[1:3], twist_moment[1:num_groups]
#
#   PARTICLE_DYNAMICS (num_points = number of WING points):
#     inputs:  point_pos[1:3,1:num_points], point_vel[1:3,1:num_points]
#     outputs: point_force[1:3,1:num_points]
#
# Everything is in the wing body frame. The wiring layer (vsm_eqs!)
# drives the inputs and reads the outputs.

# ==================== connector declarations ==================== #
#
# Connectors are declared as array variables and passed to `System`
# unflattened — MTK accepts array (and scalar) variables as unknowns.
# Listing them explicitly is required so that input connectors with no
# internal equation (a built-in ignores most of them) still exist for
# the wiring layer to bind, the same way the winch component lists
# `len`.

function _rigid_aero_connectors(num_groups::Int)
    @variables begin
        va(t)[1:3]
        rho(t)
        R_b_w(t)[1:3, 1:3]
        omega(t)[1:3]
        force(t)[1:3]
        moment(t)[1:3]
    end
    if num_groups > 0
        @variables twist(t)[1:num_groups] twist_vel(t)[1:num_groups] twist_moment(t)[1:num_groups]
    else
        twist = nothing
        twist_vel = nothing
        twist_moment = nothing
    end
    return (; va, rho, R_b_w, omega, force, moment,
            twist, twist_vel, twist_moment)
end

function _rigid_unknowns(connectors)
    vars = Any[connectors.va, connectors.rho, connectors.R_b_w,
               connectors.omega, connectors.force, connectors.moment]
    connectors.twist === nothing ||
        append!(vars, Any[connectors.twist, connectors.twist_vel,
                          connectors.twist_moment])
    return vars
end

function _particle_aero_connectors(num_points::Int)
    @variables begin
        point_pos(t)[1:3, 1:num_points]
        point_vel(t)[1:3, 1:num_points]
        point_force(t)[1:3, 1:num_points]
    end
    return (; point_pos, point_vel, point_force)
end

_particle_unknowns(connectors) =
    Any[connectors.point_pos, connectors.point_vel, connectors.point_force]

function _wing_points(sys_struct, wing)
    return [point for point in sys_struct.points
            if point.type == WING && point.wing_idx == wing.idx]
end

# ==================== NoAero ==================== #

function default_aero_none(sys_struct, wing_idx; name)
    SST = typeof(sys_struct)
    @parameters (psys::SST = sys_struct), [tunable = false]
    wing = sys_struct.wings[wing_idx]

    if wing.dynamics_type == PARTICLE_DYNAMICS
        num_points = length(_wing_points(sys_struct, wing))
        connectors = _particle_aero_connectors(num_points)
        eqs = vec(collect(connectors.point_force)) .~ 0
        return System(eqs, t, _particle_unknowns(connectors), [psys]; name)
    elseif wing.dynamics_type == RIGID_DYNAMICS
        num_groups = length(wing.group_idxs)
        connectors = _rigid_aero_connectors(num_groups)
        eqs = [collect(connectors.force) .~ 0
               collect(connectors.moment) .~ 0]
        num_groups > 0 && (eqs = [eqs; collect(connectors.twist_moment) .~ 0])
        return System(eqs, t, _rigid_unknowns(connectors), [psys]; name)
    else
        error("Unknown dynamics_type $(wing.dynamics_type) for wing $wing_idx.")
    end
end

# ==================== DiscreteAero (AERO_DIRECT) ==================== #

function default_aero_direct(sys_struct, wing_idx; name)
    SST = typeof(sys_struct)
    @parameters (psys::SST = sys_struct), [tunable = false]
    wing = sys_struct.wings[wing_idx]

    if wing.dynamics_type == PARTICLE_DYNAMICS
        points = _wing_points(sys_struct, wing)
        num_points = length(points)
        connectors = _particle_aero_connectors(num_points)
        eqs = Equation[]
        for (point_num, point) in enumerate(points)
            eqs = [eqs
                   collect(connectors.point_force[:, point_num]) .~
                       [get_point_aero_force(psys, point.idx, i)
                        for i in 1:3]]
        end
        return System(eqs, t, _particle_unknowns(connectors), [psys]; name)
    elseif wing.dynamics_type == RIGID_DYNAMICS
        groups = sys_struct.groups
        num_groups = length(wing.group_idxs)
        connectors = _rigid_aero_connectors(num_groups)
        eqs = [collect(connectors.force) .~
                   [get_aero_force_override(psys, wing.idx, i) for i in 1:3]
               collect(connectors.moment) .~
                   [get_aero_moment_override(psys, wing.idx, i) for i in 1:3]]
        for (group_pos, group_idx) in enumerate(wing.group_idxs)
            rhs = isempty(groups[group_idx].unrefined_section_idxs) ? 0 :
                get_group_moment_override(psys, wing.idx, Int64(group_idx))
            eqs = [eqs; connectors.twist_moment[group_pos] ~ rhs]
        end
        return System(eqs, t, _rigid_unknowns(connectors), [psys]; name)
    else
        error("Unknown dynamics_type $(wing.dynamics_type) for wing $wing_idx.")
    end
end

# ==================== LinearizedAero ==================== #

function default_aero_linearized(sys_struct, wing_idx; name)
    SST = typeof(sys_struct)
    @parameters (psys::SST = sys_struct), [tunable = false]
    wing = sys_struct.wings[wing_idx]

    wing.dynamics_type == PARTICLE_DYNAMICS && error(
        "AERO_LINEARIZED is not supported for PARTICLE_DYNAMICS " *
        "wings (wing $wing_idx); use AERO_DIRECT or a custom model.")
    wing isa VSMWing || error(
        "AERO_LINEARIZED wing $wing_idx is not a VSMWing.")

    groups = sys_struct.groups
    num_groups = length(wing.group_idxs)
    num_aero_inputs = length(wing.aero_y)
    area = wing.vsm_aero.projected_area
    c_ref = wing.vsm_aero.c_ref

    connectors = _rigid_aero_connectors(num_groups)
    @variables aero_input(t)[1:num_aero_inputs]

    apparent_wind = collect(connectors.va)
    omega = collect(connectors.omega)
    drag_dir = collect(apparent_wind ./ smooth_norm(apparent_wind))
    alpha = atan(drag_dir[3], drag_dir[1])
    beta = atan(drag_dir[2], smooth_norm((drag_dir[1], drag_dir[3])))

    twist_inputs = num_groups > 0 ? collect(connectors.twist) : Num[]
    input_rhs = [alpha; beta; omega[1]; omega[2]; omega[3]; twist_inputs]

    delta(input_idx) = aero_input[input_idx] - get_aero_y(psys, wing.idx, input_idx)
    coeff(output_idx) = get_aero_x(psys, wing.idx, output_idx) +
        sum(get_aero_jac(psys, wing.idx, output_idx, input_idx) * delta(input_idx)
            for input_idx in 1:num_aero_inputs)

    q_inf = 0.5 * connectors.rho * (apparent_wind ⋅ apparent_wind)
    qA = q_inf * area
    CL = coeff(1)
    CD = coeff(2)
    CS = coeff(3)

    crossed = collect(drag_dir × [0.0, 1.0, 0.0])
    lift_dir = collect(crossed ./ smooth_norm(crossed))
    side_dir = collect(lift_dir × drag_dir)
    drag_frac = get_drag_frac(psys, wing.idx)

    force_rhs = collect(qA * (CL * lift_dir +
        CD * drag_frac * drag_dir + CS * side_dir))
    moment_rhs = [qA * c_ref * coeff(3 + i) for i in 1:3]

    eqs = [collect(aero_input) .~ input_rhs
           collect(connectors.force) .~ force_rhs
           collect(connectors.moment) .~ moment_rhs]
    for group_pos in 1:num_groups
        isempty(groups[wing.group_idxs[group_pos]].unrefined_section_idxs) ?
            (eqs = [eqs; connectors.twist_moment[group_pos] ~ 0]) :
            (eqs = [eqs; connectors.twist_moment[group_pos] ~
                qA * c_ref * coeff(6 + group_pos)])
    end

    vars = _rigid_unknowns(connectors)
    push!(vars, aero_input)
    return System(eqs, t, vars, [psys]; name)
end

# ==================== PlateAero (not via component path) ==================== #

default_aero_plate(sys_struct, wing_idx; name) = error(
    "PlateWing aerodynamics use plate_eqs!, not the aero component path.")

# ==================== validation ==================== #

function validate_aero_component(subsys, wing)
    if wing.dynamics_type == RIGID_DYNAMICS
        required = Symbol[:va, :rho, :R_b_w, :omega, :force, :moment]
        length(wing.group_idxs) > 0 &&
            append!(required, [:twist, :twist_vel, :twist_moment])
    else
        required = Symbol[:point_pos, :point_vel, :point_force]
    end
    required_str = join(required, ", ")
    for con in required
        hasproperty(subsys, con) || error(
            "Wing $(wing.name): aero component is missing required " *
            "connector `$con`. Required: $required_str.")
    end
    return nothing
end
