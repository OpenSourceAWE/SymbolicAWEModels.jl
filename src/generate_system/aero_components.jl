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

# ==================== connector declarations ==================== #
#
# Connectors are declared as array variables and passed to `System`
# unflattened — MTK accepts array (and scalar) variables as unknowns.
# Listing them explicitly is required so that input connectors with no
# internal equation (a built-in ignores most of them) still exist for
# the wiring layer to bind, the same way the winch component lists
# `len`.

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
    vars = Any[c.va, c.rho, c.R_b_w, c.omega, c.force, c.moment]
    c.twist === nothing ||
        append!(vars, Any[c.twist, c.twist_vel, c.twist_moment])
    return vars
end

function _particle_aero_connectors(np::Int)
    @variables begin
        point_pos(t)[1:3, 1:np]
        point_vel(t)[1:3, 1:np]
        point_force(t)[1:3, 1:np]
    end
    return (; point_pos, point_vel, point_force)
end

_particle_unknowns(c) = Any[c.point_pos, c.point_vel, c.point_force]

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
        eqs = vec(collect(c.point_force)) .~ 0
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

    vars = _rigid_unknowns(c)
    push!(vars, aero_input)
    return System(eqs, t, vars, [psys]; name)
end

# ==================== PlateAero (flat-plate polar groups) ==================== #
#
# Flat-plate aero for wings whose groups carry a polar (calc_cl/calc_cd).
# Each group is a 1-point FIXED panel; its twist is prescribed. The math is
# done entirely in the wing body frame.
#
# Extra (opportunistic) PARTICLE connectors beyond the base contract:
#   point_va  — body-frame apparent wind per point (bound to `va_point_b`)
#   point_rho — air density per point (bound to `calc_rho(am, height)`)

function _group_of_point(groups, group_idxs, point_idx)
    for gidx in group_idxs
        point_idx in groups[gidx].point_idxs && return gidx
    end
    error("Point $point_idx is not a member of any group on its wing.")
end

"""
    _plate_group_force_b(psys, gidx, x_airf, y_airf, twist, va_b, rho)

Body-frame flat-plate force on a polar group. `x_airf`/`y_airf` are the
build-time chord/span directions (body frame); `twist`, `va_b`, `rho` are
symbolic. Lift perpendicular to the in-plane flow, drag along it.
"""
function _plate_group_force_b(psys, gidx, x_airf, y_airf, twist, va_b, rho)
    ct = cos(twist)
    st = sin(twist)
    x_tw = ct .* x_airf .+ st .* cross(y_airf, x_airf)
    z_tw = cross(x_tw, y_airf)
    v_tan = dot(va_b, x_tw)
    v_norm = dot(va_b, z_tw)
    alpha = rad2deg(atan(v_norm, v_tan))
    cl = get_group_cl(psys, gidx, alpha)
    cd = get_group_drag_corr(psys, gidx) * get_group_cd(psys, gidx, alpha)
    q = 0.5 * rho * (v_tan^2 + v_norm^2)
    q_drag = 0.5 * rho * dot(va_b, va_b)
    a_rad = atan(v_norm, v_tan)
    va_airf_dir = cos(a_rad) .* x_tw .+ sin(a_rad) .* z_tw
    lift_dir = smooth_normalize(cross(va_airf_dir, y_airf))
    drag_dir = smooth_normalize(cross(y_airf, lift_dir))
    area = get_group_area(psys, gidx)
    lift = (q * area * cl) .* lift_dir
    drag = (q_drag * area * cd) .* drag_dir
    return collect(lift .+ drag)
end

function default_aero_plate(sys_struct, wing_idx; name)
    SST = typeof(sys_struct)
    @parameters (psys::SST = sys_struct), [tunable = false]
    wing = sys_struct.wings[wing_idx]
    groups = sys_struct.groups
    group_idxs = wing.group_idxs
    isempty(group_idxs) && error(
        "AERO_PLATE wing $wing_idx has no polar groups.")

    if wing.dynamics_type == PARTICLE_DYNAMICS
        points = _wing_points(sys_struct, wing)
        np = length(points)
        @variables begin
            point_pos(t)[1:3, 1:np]
            point_vel(t)[1:3, 1:np]
            point_force(t)[1:3, 1:np]
            point_va(t)[1:3, 1:np]
            point_rho(t)[1:np]
        end
        eqs = Equation[]
        for (k, point) in enumerate(points)
            gidx = _group_of_point(groups, group_idxs, point.idx)
            g = groups[gidx]
            x_airf = normalize(Vector(g.chord))
            y_airf = Vector(g.y_airf)
            tw = get_twist(psys, gidx)
            va_b = collect(point_va[:, k])
            f = _plate_group_force_b(psys, gidx, x_airf, y_airf,
                                     tw, va_b, point_rho[k])
            eqs = [eqs; collect(point_force[:, k]) .~ f]
        end
        vars = Any[point_pos, point_vel, point_force, point_va, point_rho]
        return System(eqs, t, vars, [psys]; name)
    end

    # RIGID_DYNAMICS: sum group forces and moments about the COM.
    ng = length(group_idxs)
    c = _rigid_aero_connectors(ng)
    va = collect(c.va)
    omega = collect(c.omega)
    force = zeros(Num, 3)
    moment = zeros(Num, 3)
    for (gi, gidx) in enumerate(group_idxs)
        g = groups[gidx]
        x_airf = normalize(Vector(g.chord))
        y_airf = Vector(g.y_airf)
        r_g = Vector(sys_struct.points[g.point_idxs[1]].pos_b)
        va_g = va .- cross(omega, r_g)
        f = _plate_group_force_b(psys, gidx, x_airf, y_airf,
                                 c.twist[gi], va_g, c.rho)
        force = force .+ f
        moment = moment .+ cross(r_g, f)
    end
    eqs = [collect(c.force) .~ force
           collect(c.moment) .~ moment]
    for gi in 1:ng
        eqs = [eqs; c.twist_moment[gi] ~ 0]
    end
    return System(eqs, t, _rigid_unknowns(c), [psys]; name)
end

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
