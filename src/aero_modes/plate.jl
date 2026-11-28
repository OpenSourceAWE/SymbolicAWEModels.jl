# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# AeroPlate: flat-plate CL/CD-lookup aerodynamics (PARTICLE_DYNAMICS only).

"""
    AeroPlate(calc_cl, calc_cd; drag_corr=1.0)

Flat-plate CL/CD lookup aerodynamics. Carries the shared polar lookups
(`calc_cl`/`calc_cd`: `α_deg → coefficient`) and drag correction used by all of
a wing's flat-plate (1-point `STATIC`) spanwise stations. One polar set per
wing.
"""
mutable struct AeroPlate{CL, CD} <: AbstractAeroModel
    calc_cl::CL
    calc_cd::CD
    drag_corr::SimFloat
end

AeroPlate(calc_cl, calc_cd; drag_corr=1.0) =
    AeroPlate{typeof(calc_cl), typeof(calc_cd)}(
        calc_cl, calc_cd, SimFloat(drag_corr))

"""
    AeroPlate()

A flat plate whose polars a model takes from its settings ([`model_aero`](@ref)).
"""
AeroPlate() = AeroPlate(nothing, nothing)

"""A flat plate with the polars of `set` where `mode` gives none, and the wing row's
`drag_corr`."""
function model_aero(mode::AeroPlate{Nothing}, wing, n_stations, set, vsm_set)
    isnothing(set) && return mode
    calc_cl, calc_cd = create_plate_interpolations(set.alpha_cl, set.cl_list, set.cd_list;
                                                   alpha_cd=set.alpha_cd)
    return AeroPlate(calc_cl, calc_cd; drag_corr=something(wing.model.drag_corr, 0.93))
end

is_builtin_aero(::AeroPlate) = true
aero_mode_tag(::Type{<:AeroPlate}) = "plate"
function calc_aoa(::AeroPlate, sam, wing_idx)
    va_b = wing_apparent_wind(sam, wing_idx)
    return atan(va_b[3], va_b[1])
end

# ==================== polar interpolation + accessors ==================== #

"""
    create_plate_interpolations(alpha_deg, cl_data, cd_data;
        alpha_cd=nothing, spline=:cubic)

Create CL and CD interpolation objects from polar data vectors.

# Arguments
- `alpha_deg`: angle of attack values [deg]
- `cl_data`: lift coefficient values
- `cd_data`: drag coefficient values
- `alpha_cd`: separate alpha values for CD (default: same
  as CL)
- `spline`: `:cubic` for cubic spline, `:linear` for
  piecewise linear

# Returns
- `(cl_interp, cd_interp)` tuple of interpolation objects
"""
function create_plate_interpolations(
    alpha_cl, cl_data, cd_data;
    alpha_cd=nothing, spline=:cubic
)
    alpha_cd_vec = isnothing(alpha_cd) ? alpha_cl : alpha_cd
    if spline == :cubic
        cl_interp = CubicSpline(
            Vector{Float64}(cl_data),
            Vector{Float64}(alpha_cl))
        cd_interp = CubicSpline(
            Vector{Float64}(cd_data),
            Vector{Float64}(alpha_cd_vec))
    elseif spline == :linear
        cl_interp = LinearInterpolation(
            Vector{Float64}(cl_data),
            Vector{Float64}(alpha_cl))
        cd_interp = LinearInterpolation(
            Vector{Float64}(cd_data),
            Vector{Float64}(alpha_cd_vec))
    else
        error("Unknown spline type: $spline. " *
              "Use :cubic or :linear.")
    end
    return (cl_interp, cd_interp)
end

# ==================== equation builder ==================== #

"""
    aero_component(::AeroPlate, ::ParticleDynamics, sam, wing_idx; name, params)

Flat-plate aero component of wing `wing_idx` of `sam`. Uses the same `PARTICLE_DYNAMICS`
connector contract as the other per-point modes; it is the only one that consumes the
`va`/`rho` inputs (the VSM particle modes read frozen forces and ignore them). Each wing
node is a 1-point `STATIC` spanwise station; the per-point force is computed from the
section's twisted body-frame axes, the point's apparent wind, and its air density.
"""
function aero_component(::AeroPlate, ::ParticleDynamics, sam, wing_idx; name,
params=nothing)
    stations = sam.definition.spanwise_stations
    station_idxs = wing_station_idxs(sam, wing_idx)
    point_idxs = wing_points(sam, wing_idx)
    connectors = particle_aero_connectors(length(point_idxs))

    # Airfoil polars as callable params (live splines, applied as `cl(α)`).
    calc_cl_p = params.wings[wing_idx].aero.calc_cl
    calc_cd_p = params.wings[wing_idx].aero.calc_cd
    flat_ps = Any[calc_cl_p, calc_cd_p]
    eqs = Equation[]
    for (point_num, i) in enumerate(point_idxs)
        position = findfirst(idx -> stations[idx].points[1] == i, station_idxs)
        isnothing(position) && error("Wing $wing_idx: wing node $i is not a flat-plate " *
                                     "section point.")
        station_idx = station_idxs[position]
        chord_p = params.spanwise_stations[station_idx].chord
        y_airf_p = params.spanwise_stations[station_idx].y_airf
        twist_p = params.spanwise_stations[station_idx].twist
        drag_corr_p = params.wings[wing_idx].aero.drag_corr
        area_p = params.spanwise_stations[station_idx].area
        append!(flat_ps, (chord_p, y_airf_p, twist_p, drag_corr_p, area_p))
        x_airf = smooth_normalize(collect(chord_p))
        y_airf = collect(y_airf_p)
        twist = twist_p
        x_twisted = cos(twist) * x_airf + sin(twist) * (y_airf × x_airf)
        z_twisted = x_twisted × y_airf

        apparent_wind = collect(connectors.va[:, point_num])
        v_tan = apparent_wind ⋅ x_twisted
        v_norm = apparent_wind ⋅ z_twisted
        alpha_deg = rad2deg(atan(v_norm, v_tan))

        cl = calc_cl_p(alpha_deg)
        cd = drag_corr_p * calc_cd_p(alpha_deg)

        q = 0.5 * connectors.rho[point_num] * (v_tan^2 + v_norm^2)
        q_drag = 0.5 * connectors.rho[point_num] * (apparent_wind ⋅ apparent_wind)

        alpha_rad = atan(v_norm, v_tan)
        va_airf_dir = cos(alpha_rad) * x_twisted + sin(alpha_rad) * z_twisted
        lift_dir = smooth_normalize(va_airf_dir × y_airf)
        drag_dir = smooth_normalize(y_airf × lift_dir)

        area = area_p
        push!(eqs, connectors.point_force[:, point_num] ~
                   q * area * cl * lift_dir + q_drag * area * cd * drag_dir)
    end

    return System(eqs, t, particle_unknowns(connectors), unique(flat_ps); name)
end

# ==================== log-point hooks ==================== #

"""
    plate_corners(sam, station_idx, R_b_to_w) -> NTuple{4, Vector}

World-frame corners of the display quad of flat-plate station `station_idx` of `sam`, its
wing's body turned by `R_b_to_w`. The station's structural point sits at quarter chord; the
quad is a square of side `sqrt(area)` spanned by the twisted chord direction and
`y_airf`, so the quad area matches the station's `area`.
"""
function plate_corners(sam, station_idx, R_b_to_w)
    station(quantity) = run_value(sam, :spanwise_stations, station_idx, quantity)
    x_airf = normalize(station(:chord))
    y_airf = station(:y_airf)
    twist = station(:twist)
    x_twisted = cos(twist) * x_airf + sin(twist) * (y_airf × x_airf)
    side = sqrt(station(:area))
    chord_w = R_b_to_w * (side * x_twisted)
    span_half_w = R_b_to_w * (0.5 * side * y_airf)
    point_pos_w = run_value(sam, Val(:points),
                            sam.definition.spanwise_stations[station_idx].points[1],
                            Val(:pos_w))
    le_mid = point_pos_w - 0.25 * chord_w
    te_mid = le_mid + chord_w
    return (le_mid + span_half_w, te_mid + span_half_w,
            te_mid - span_half_w, le_mid - span_half_w)
end

"""
    n_aero_log_points(::AeroPlate, sys, wing_idx) -> Int

4 quad corners per flat-plate station of wing `wing_idx` of definition `sys`
([`plate_corners`](@ref)).
"""
n_aero_log_points(::AeroPlate, sys, wing_idx) = 4 * length(station_indices(sys, wing_idx))

"""
    write_aero_log_points!(::AeroPlate, sam, wing_idx, sys_state, point_idx, zoom) -> Int

Log each flat-plate station's display quad of wing `wing_idx` of `sam`
([`plate_corners`](@ref)).
"""
function write_aero_log_points!(::AeroPlate, sam, wing_idx, sys_state, point_idx, zoom)
    R_b_to_w = body_rotation(sam, wing_body(sam, wing_idx))
    for station_idx in wing_station_idxs(sam, wing_idx),
            corner_w in plate_corners(sam, station_idx, R_b_to_w)
        point_idx += 1
        sys_state.X[point_idx] = corner_w[1] * zoom
        sys_state.Y[point_idx] = corner_w[2] * zoom
        sys_state.Z[point_idx] = corner_w[3] * zoom
    end
    return point_idx
end

"""
    read_aero_log_points!(::AeroPlate, sam, wing_idx, sys_state, point_idx) -> Int

Skip the quad-corner slots: plate corners are derived from the restored
structural point positions, so nothing is read back.
"""
read_aero_log_points!(mode::AeroPlate, sam, wing_idx, sys_state, point_idx) =
    point_idx + n_aero_log_points(mode, sam.definition, wing_idx)
