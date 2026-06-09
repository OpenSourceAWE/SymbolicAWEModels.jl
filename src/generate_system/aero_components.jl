# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Aero coupling components (winch-style swappable subsystems).
#
# A wing carries an `aero::AbstractAeroModel`; `aero_component(mode, …)` is
# dispatched on its type and returns a `System` whose connectors are fixed by
# the wing's `dynamics_type`:
#
#   RIGID_DYNAMICS (num_twist_surfaces = length(wing.twist_surface_idxs)):
#     inputs:  va[1:3], rho, R_b_w[1:3,1:3], omega[1:3],
#              twist[1:num_twist_surfaces], twist_vel[1:num_twist_surfaces]
#     outputs: force[1:3], moment[1:3], twist_moment[1:num_twist_surfaces]
#
#   PARTICLE_DYNAMICS (num_points = number of WING points):
#     inputs:  point_pos[1:3,1:np], point_vel[1:3,1:np],
#              va[1:3,1:np], rho[1:np]
#     outputs: point_force[1:3,1:np]
#   (VSM particle modes read frozen forces and ignore va/rho; flat-plate uses
#    them to compute per-point forces symbolically.)
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

function _rigid_aero_connectors(num_twist_surfaces::Int)
    @variables begin
        va(t)[1:3]
        rho(t)
        R_b_w(t)[1:3, 1:3]
        omega(t)[1:3]
        force(t)[1:3]
        moment(t)[1:3]
    end
    if num_twist_surfaces > 0
        @variables twist(t)[1:num_twist_surfaces] twist_vel(t)[1:num_twist_surfaces] twist_moment(t)[1:num_twist_surfaces]
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
        va(t)[1:3, 1:num_points]
        rho(t)[1:num_points]
        point_force(t)[1:3, 1:num_points]
    end
    return (; point_pos, point_vel, va, rho, point_force)
end

_particle_unknowns(connectors) =
    Any[connectors.point_pos, connectors.point_vel, connectors.va,
        connectors.rho, connectors.point_force]

function _wing_points(sys_struct, wing)
    return [point for point in sys_struct.points
            if point.type == WING && point.wing_idx == wing.idx]
end

# ==================== dispatch ==================== #

"""
    aero_component(mode::AbstractAeroModel, sys_struct, wing_idx; name) -> System

Build the aero subsystem for `sys_struct.wings[wing_idx]`, selected by dispatch
on the wing's `aero` model. Returns a `System` exposing the connectors fixed by
the wing's `dynamics_type` (see this file's header). Add a method on a custom
`AbstractAeroModel` subtype to plug in your own aerodynamics.
"""
function aero_component end

# ==================== NoAero ==================== #

function aero_component(::AeroNone, sys_struct, wing_idx; name)
    psys = system_struct_param(sys_struct)
    wing = sys_struct.wings[wing_idx]

    if wing.dynamics_type == PARTICLE_DYNAMICS
        num_points = length(_wing_points(sys_struct, wing))
        connectors = _particle_aero_connectors(num_points)
        eqs = vec(collect(connectors.point_force)) .~ 0
        return System(eqs, t, _particle_unknowns(connectors), [psys]; name)
    elseif wing.dynamics_type == RIGID_DYNAMICS
        num_twist_surfaces = length(wing.twist_surface_idxs)
        connectors = _rigid_aero_connectors(num_twist_surfaces)
        eqs = [collect(connectors.force) .~ 0
               collect(connectors.moment) .~ 0]
        num_twist_surfaces > 0 && (eqs = [eqs; collect(connectors.twist_moment) .~ 0])
        return System(eqs, t, _rigid_unknowns(connectors), [psys]; name)
    else
        error("Unknown dynamics_type $(wing.dynamics_type) for wing $wing_idx.")
    end
end

# ==================== DiscreteAero (AeroDirect) ==================== #

function aero_component(::AeroDirect, sys_struct, wing_idx; name)
    psys = system_struct_param(sys_struct)
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
        twist_surfaces = sys_struct.twist_surfaces
        num_twist_surfaces = length(wing.twist_surface_idxs)
        connectors = _rigid_aero_connectors(num_twist_surfaces)
        eqs = [collect(connectors.force) .~
                   [get_aero_force_override(psys, wing.idx, i) for i in 1:3]
               collect(connectors.moment) .~
                   [get_aero_moment_override(psys, wing.idx, i) for i in 1:3]]
        for (twist_surface_pos, twist_surface_idx) in enumerate(wing.twist_surface_idxs)
            rhs = isempty(twist_surfaces[twist_surface_idx].unrefined_section_idxs) ? 0 :
                get_twist_surface_moment_override(psys, wing.idx, Int64(twist_surface_idx))
            eqs = [eqs; connectors.twist_moment[twist_surface_pos] ~ rhs]
        end
        return System(eqs, t, _rigid_unknowns(connectors), [psys]; name)
    else
        error("Unknown dynamics_type $(wing.dynamics_type) for wing $wing_idx.")
    end
end

# ==================== LinearizedAero ==================== #

function aero_component(::AeroLinearized, sys_struct, wing_idx; name)
    psys = system_struct_param(sys_struct)
    wing = sys_struct.wings[wing_idx]

    wing.dynamics_type == PARTICLE_DYNAMICS && error(
        "AeroLinearized is not supported for PARTICLE_DYNAMICS " *
        "wings (wing $wing_idx); use AeroDirect or a custom model.")
    is_vsm(wing) || error(
        "AeroLinearized wing $wing_idx has no VSM engine.")

    twist_surfaces = sys_struct.twist_surfaces
    num_twist_surfaces = length(wing.twist_surface_idxs)
    num_aero_inputs = length(wing.aero_y)
    area = wing.vsm_aero.projected_area
    c_ref = wing.vsm_aero.c_ref

    connectors = _rigid_aero_connectors(num_twist_surfaces)
    @variables aero_input(t)[1:num_aero_inputs]

    apparent_wind = collect(connectors.va)
    omega = collect(connectors.omega)
    drag_dir = collect(apparent_wind ./ smooth_norm(apparent_wind))
    alpha = atan(drag_dir[3], drag_dir[1])
    beta = atan(drag_dir[2], smooth_norm((drag_dir[1], drag_dir[3])))

    twist_inputs = num_twist_surfaces > 0 ? collect(connectors.twist) : Num[]
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
    for twist_surface_pos in 1:num_twist_surfaces
        isempty(twist_surfaces[wing.twist_surface_idxs[twist_surface_pos]].unrefined_section_idxs) ?
            (eqs = [eqs; connectors.twist_moment[twist_surface_pos] ~ 0]) :
            (eqs = [eqs; connectors.twist_moment[twist_surface_pos] ~
                qA * c_ref * coeff(6 + twist_surface_pos)])
    end

    vars = _rigid_unknowns(connectors)
    push!(vars, aero_input)
    return System(eqs, t, vars, [psys]; name)
end

# ==================== PlateAero ==================== #
#
# Flat-plate aero uses the same PARTICLE connector contract as the other
# per-point modes; it is the only one that consumes the `va`/`rho` inputs (the
# VSM particle modes read frozen forces and ignore them). Each WING point is a
# 1-point `FIXED` TwistSurface section; the per-point force is computed from the
# section's twisted body-frame axes, the point's apparent wind, and its density.

function aero_component(::AeroPlate, sys_struct, wing_idx; name)
    psys = system_struct_param(sys_struct)
    wing = sys_struct.wings[wing_idx]
    wing.dynamics_type == RIGID_DYNAMICS && error(
        "AeroPlate is only supported for PARTICLE_DYNAMICS wings (wing $wing_idx).")

    twist_surfaces = sys_struct.twist_surfaces
    points = _wing_points(sys_struct, wing)
    num_points = length(points)
    connectors = _particle_aero_connectors(num_points)

    eqs = Equation[]
    for (point_num, point) in enumerate(points)
        ts_idx = 0
        for gidx in wing.twist_surface_idxs
            if twist_surfaces[gidx].point_idxs[1] == point.idx
                ts_idx = gidx
                break
            end
        end
        ts_idx == 0 && error(
            "Wing $wing_idx: WING point $(point.idx) is not a flat-plate " *
            "section point.")

        x_airf = smooth_normalize(collect(get_twist_surface_chord(psys, ts_idx)))
        y_airf = collect(get_twist_surface_y_airf(psys, ts_idx))
        twist = get_twist(psys, ts_idx)
        x_twisted = cos(twist) * x_airf + sin(twist) * (y_airf × x_airf)
        z_twisted = x_twisted × y_airf

        apparent_wind = collect(connectors.va[:, point_num])
        v_tan = apparent_wind ⋅ x_twisted
        v_norm = apparent_wind ⋅ z_twisted
        alpha_deg = rad2deg(atan(v_norm, v_tan))

        cl = get_plate_cl(psys, wing_idx, alpha_deg)
        cd = get_plate_drag_corr(psys, wing_idx) *
             get_plate_cd(psys, wing_idx, alpha_deg)

        q = 0.5 * connectors.rho[point_num] * (v_tan^2 + v_norm^2)
        q_drag = 0.5 * connectors.rho[point_num] * (apparent_wind ⋅ apparent_wind)

        alpha_rad = atan(v_norm, v_tan)
        va_airf_dir = cos(alpha_rad) * x_twisted + sin(alpha_rad) * z_twisted
        lift_dir = smooth_normalize(va_airf_dir × y_airf)
        drag_dir = smooth_normalize(y_airf × lift_dir)

        area = get_twist_surface_area(psys, ts_idx)
        eqs = [eqs
               connectors.point_force[:, point_num] ~
                   q * area * cl * lift_dir + q_drag * area * cd * drag_dir]
    end

    return System(eqs, t, _particle_unknowns(connectors), [psys]; name)
end

# ==================== validation ==================== #

function validate_aero_component(subsys, wing)
    if wing.dynamics_type == RIGID_DYNAMICS
        required = Symbol[:va, :rho, :R_b_w, :omega, :force, :moment]
        length(wing.twist_surface_idxs) > 0 &&
            append!(required, [:twist, :twist_vel, :twist_moment])
    else
        required = Symbol[:point_pos, :point_vel, :va, :rho, :point_force]
    end
    required_str = join(required, ", ")
    for con in required
        hasproperty(subsys, con) || error(
            "Wing $(wing.name): aero component is missing required " *
            "connector `$con`. Required: $required_str.")
    end
    return nothing
end
