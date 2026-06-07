# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Aero coupling components (winch-style swappable subsystems).
#
# A wing carries an `aero_model` builder selected by `aero_mode`
# (see `resolve_aero_model`). Each builder returns a `System` whose
# connectors are fixed by the wing's `dynamics_type`:
#
#   RIGID_DYNAMICS (ng = length(wing.group_idxs)):
#     inputs:  va[1:3], rho, R_b_w[1:3,1:3], omega[1:3],
#              twist[1:ng], twist_vel[1:ng]            (ng > 0 only)
#     outputs: force[1:3], moment[1:3], twist_moment[1:ng]
#
#   PARTICLE_DYNAMICS (np = number of WING points):
#     inputs:  point_pos[1:3,1:np], point_vel[1:3,1:np]
#     outputs: point_force[1:3,1:np]
#
# Everything is in the wing body frame. The wiring layer (vsm_eqs!)
# drives the inputs and reads the outputs.

_flatten_vars(x) = x === nothing ? Num[] : vec(collect(x))

# ==================== connector declarations ==================== #

function _rigid_aero_connectors(ng::Int)
    @variables begin
        va(t)[1:3]
        rho(t)
        R_b_w(t)[1:3, 1:3]
        omega(t)[1:3]
        force(t)[1:3]
        moment(t)[1:3]
    end
    if ng > 0
        @variables twist(t)[1:ng] twist_vel(t)[1:ng] twist_moment(t)[1:ng]
    else
        twist = nothing
        twist_vel = nothing
        twist_moment = nothing
    end
    return (; va, rho, R_b_w, omega, force, moment,
            twist, twist_vel, twist_moment)
end

function _rigid_unknowns(c)
    return [_flatten_vars(c.va); c.rho; _flatten_vars(c.R_b_w);
            _flatten_vars(c.omega); _flatten_vars(c.force);
            _flatten_vars(c.moment); _flatten_vars(c.twist);
            _flatten_vars(c.twist_vel); _flatten_vars(c.twist_moment)]
end

function _particle_aero_connectors(np::Int)
    @variables begin
        point_pos(t)[1:3, 1:np]
        point_vel(t)[1:3, 1:np]
        point_force(t)[1:3, 1:np]
    end
    return (; point_pos, point_vel, point_force)
end

_particle_unknowns(c) =
    [_flatten_vars(c.point_pos); _flatten_vars(c.point_vel);
     _flatten_vars(c.point_force)]

function _wing_points(sys_struct, wing)
    return [p for p in sys_struct.points
            if p.type == WING && p.wing_idx == wing.idx]
end

# ==================== NoAero ==================== #

function default_aero_none(sys_struct, wing_idx; name)
    SST = typeof(sys_struct)
    @parameters (psys::SST = sys_struct), [tunable = false]
    wing = sys_struct.wings[wing_idx]

    if wing.dynamics_type == PARTICLE_DYNAMICS
        np = length(_wing_points(sys_struct, wing))
        c = _particle_aero_connectors(np)
        eqs = [collect(c.point_force) .~ 0]
        return System(eqs, t, _particle_unknowns(c), [psys]; name)
    end

    ng = length(wing.group_idxs)
    c = _rigid_aero_connectors(ng)
    eqs = [collect(c.force) .~ 0
           collect(c.moment) .~ 0]
    ng > 0 && (eqs = [eqs; collect(c.twist_moment) .~ 0])
    return System(eqs, t, _rigid_unknowns(c), [psys]; name)
end

# ==================== DiscreteAero (AERO_DIRECT) ==================== #

function default_aero_direct(sys_struct, wing_idx; name)
    SST = typeof(sys_struct)
    @parameters (psys::SST = sys_struct), [tunable = false]
    wing = sys_struct.wings[wing_idx]

    if wing.dynamics_type == PARTICLE_DYNAMICS
        points = _wing_points(sys_struct, wing)
        np = length(points)
        c = _particle_aero_connectors(np)
        eqs = Equation[]
        for (k, point) in enumerate(points)
            eqs = [eqs
                   collect(c.point_force[:, k]) .~
                       [get_point_aero_force(psys, point.idx, i)
                        for i in 1:3]]
        end
        return System(eqs, t, _particle_unknowns(c), [psys]; name)
    end

    groups = sys_struct.groups
    ng = length(wing.group_idxs)
    c = _rigid_aero_connectors(ng)
    w = wing.idx
    eqs = [collect(c.force) .~
               [get_aero_force_override(psys, w, i) for i in 1:3]
           collect(c.moment) .~
               [get_aero_moment_override(psys, w, i) for i in 1:3]]
    for (gi, gidx) in enumerate(wing.group_idxs)
        rhs = isempty(groups[gidx].unrefined_section_idxs) ? 0 :
            get_group_moment_override(psys, w, Int64(gidx))
        eqs = [eqs; c.twist_moment[gi] ~ rhs]
    end
    return System(eqs, t, _rigid_unknowns(c), [psys]; name)
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
    ng = length(wing.group_idxs)
    ny = length(wing.aero_y)
    nx = length(wing.aero_x)
    w = wing.idx
    area = wing.vsm_aero.projected_area
    c_ref = wing.vsm_aero.c_ref

    c = _rigid_aero_connectors(ng)
    @variables aero_input(t)[1:ny]

    va = collect(c.va)
    omega = collect(c.omega)
    drag_dir = collect(va ./ smooth_norm(va))
    alpha = atan(drag_dir[3], drag_dir[1])
    beta = atan(drag_dir[2], smooth_norm((drag_dir[1], drag_dir[3])))

    twist_inputs = ng > 0 ? collect(c.twist) : Num[]
    input_rhs = [alpha; beta; omega[1]; omega[2]; omega[3]; twist_inputs]

    delta(iy) = aero_input[iy] - get_aero_y(psys, w, iy)
    coeff(ix) = get_aero_x(psys, w, ix) +
        sum(get_aero_jac(psys, w, ix, iy) * delta(iy) for iy in 1:ny)

    q_inf = 0.5 * c.rho * (va ⋅ va)
    qA = q_inf * area
    CL = coeff(1)
    CD = coeff(2)
    CS = coeff(3)

    crossed = collect(drag_dir × [0.0, 1.0, 0.0])
    lift_dir = collect(crossed ./ smooth_norm(crossed))
    side_dir = collect(lift_dir × drag_dir)
    drag_frac = get_drag_frac(psys, w)

    force_rhs = collect(qA * (CL * lift_dir +
        CD * drag_frac * drag_dir + CS * side_dir))
    moment_rhs = [qA * c_ref * coeff(3 + i) for i in 1:3]

    eqs = [collect(aero_input) .~ input_rhs
           collect(c.force) .~ force_rhs
           collect(c.moment) .~ moment_rhs]
    for gi in 1:ng
        isempty(groups[wing.group_idxs[gi]].unrefined_section_idxs) ?
            (eqs = [eqs; c.twist_moment[gi] ~ 0]) :
            (eqs = [eqs; c.twist_moment[gi] ~ qA * c_ref * coeff(6 + gi)])
    end

    unknowns = [_rigid_unknowns(c); collect(aero_input)]
    return System(eqs, t, unknowns, [psys]; name)
end

# ==================== PlateAero (not via component path) ==================== #

default_aero_plate(sys_struct, wing_idx; name) = error(
    "PlateWing aerodynamics use plate_eqs!, not the aero component path.")

# ==================== validation ==================== #

"""
    validate_aero_component(subsys, wing)

Check that `subsys` (built by `wing.aero_model`) exposes the
connector contract for the wing's `dynamics_type`.
"""
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
