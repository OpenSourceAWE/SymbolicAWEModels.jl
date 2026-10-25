# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# AeroPlate: flat-plate CL/CD-lookup aerodynamics (PARTICLE_DYNAMICS only).

"""
    AeroPlate(calc_cl, calc_cd; drag_corr=1.0)

Flat-plate CL/CD lookup aerodynamics. Carries the shared polar lookups
(`calc_cl`/`calc_cd`: `α_deg → coefficient`) and drag correction used by all of
a wing's flat-plate (1-point `STATIC`) [`Station`](@ref)s. One polar set per
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

Polar-less marker, used only to select the flat-plate path when parsing a YAML
`aero_mode`; the real lookups are attached when the flat-plate wing is built.
"""
AeroPlate() = AeroPlate(nothing, nothing)

is_builtin_aero(::AeroPlate) = true
aero_mode_tag(::AeroPlate) = "plate"
calc_aoa(::AeroPlate, wing) = atan(wing.va_b[3], wing.va_b[1])

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
    aero_component(::AeroPlate, wing::ParticleWingState, sys_struct; name)

Flat-plate aero component. Uses the same `PARTICLE_DYNAMICS` connector contract as
the other per-point modes; it is the only one that consumes the `va`/`rho` inputs
(the VSM particle modes read frozen forces and ignore them). Each wing node is a
1-point `STATIC` [`Station`](@ref) section; the per-point force is computed
from the section's twisted body-frame axes, the point's apparent wind, and its
air density.
"""
function aero_component(::AeroPlate, wing::ParticleWingState, sys_struct;
                        name, params=nothing)
    wing_idx = wing.idx
    stations = sys_struct.stations
    points = wing_points(sys_struct, wing)
    num_points = length(points)
    connectors = particle_aero_connectors(num_points)

    # Airfoil polars as callable params (live splines, applied as `cl(α)`).
    calc_cl_p = params.wings[wing_idx].aero.calc_cl
    calc_cd_p = params.wings[wing_idx].aero.calc_cd
    flat_ps = Any[calc_cl_p, calc_cd_p]
    eqs = Equation[]
    for (point_num, point) in enumerate(points)
        ts_idx = 0
        for gidx in wing.station_idxs
            if stations[gidx].point_idxs[1] == point.idx
                ts_idx = gidx
                break
            end
        end
        ts_idx == 0 && error(
            "Wing $wing_idx: wing node $(point.idx) is not a flat-plate " *
            "section point.")

        chord_p = params.stations[ts_idx].chord
        y_airf_p = params.stations[ts_idx].y_airf
        twist_p = params.stations[ts_idx].twist
        drag_corr_p = params.wings[wing_idx].aero.drag_corr
        area_p = params.stations[ts_idx].area
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
        eqs = [eqs
               connectors.point_force[:, point_num] ~
                   q * area * cl * lift_dir + q_drag * area * cd * drag_dir]
    end

    return System(eqs, t, particle_unknowns(connectors), unique(flat_ps); name)
end

# ==================== log-point hooks ==================== #

"""
    plate_corners(station, point_pos_w, R_b_to_w) -> NTuple{4, Vector}

World-frame corners of a flat-plate section's display quad. The section's
structural point sits at quarter chord; the quad is a square of side
`sqrt(area)` spanned by the twisted chord direction and `y_airf`, so the quad
area matches the section's `area`.
"""
function plate_corners(station, point_pos_w, R_b_to_w)
    x_airf = normalize(station.chord)
    y_airf = station.y_airf
    twist = station.twist
    x_twisted = cos(twist) * x_airf + sin(twist) * (y_airf × x_airf)
    side = sqrt(station.area)
    chord_w = R_b_to_w * (side * x_twisted)
    span_half_w = R_b_to_w * (0.5 * side * y_airf)
    le_mid = point_pos_w - 0.25 * chord_w
    te_mid = le_mid + chord_w
    return (le_mid + span_half_w, te_mid + span_half_w,
            te_mid - span_half_w, le_mid - span_half_w)
end

"""
    n_aero_log_points(::AeroPlate, wing) -> Int

4 quad corners per flat-plate section ([`plate_corners`](@ref)).
"""
n_aero_log_points(::AeroPlate, wing) = 4 * length(wing.station_idxs)

"""
    write_aero_log_points!(::AeroPlate, wing, sys_struct, sys_state,
                           point_idx, zoom) -> Int

Log each flat-plate section's display quad ([`plate_corners`](@ref)).
"""
function write_aero_log_points!(::AeroPlate, wing, sys_struct, sys_state,
                                point_idx, zoom)
    R_b_to_w = wing.R_b_to_w::Matrix{SimFloat}
    for station_idx in wing.station_idxs
        station = sys_struct.stations[station_idx]
        point = sys_struct.points[station.point_idxs[1]]
        for corner_w in plate_corners(station, point.pos_w, R_b_to_w)
            point_idx += 1
            sys_state.X[point_idx] = corner_w[1] * zoom
            sys_state.Y[point_idx] = corner_w[2] * zoom
            sys_state.Z[point_idx] = corner_w[3] * zoom
        end
    end
    return point_idx
end

"""
    read_aero_log_points!(::AeroPlate, wing, sys_struct, sys_state,
                          point_idx) -> Int

Skip the quad-corner slots: plate corners are derived from the restored
structural point positions, so nothing is read back.
"""
read_aero_log_points!(mode::AeroPlate, wing, sys_struct, sys_state,
                      point_idx) = point_idx + n_aero_log_points(mode, wing)

# ==================== YAML construction ==================== #

function load_wing(mode::AeroPlate, definition, w, set, dynamics_type, vsm_set,
                   station_refs, stations; transform)
    return load_plate_wing(definition, w, set, dynamics_type, mode, stations;
                           transform)
end

"""
    load_plate_wing(definition, w, set, dynamics_type, aero_mode, stations; transform)

Load flat-plate wing `w` of `definition` from its wing row and the `surfaces` block of
its extras. Each surface becomes a 1-point `STATIC` [`Station`](@ref) appended to
`stations`; the wing references them by name. CL/CD interpolations come from `Settings`
polar data.
"""
function load_plate_wing(definition, w, set, dynamics_type, aero_mode, stations;
                         transform)
    wing = definition.wings[w]
    name = row_name(wing)
    cl_interp, cd_interp = create_plate_interpolations(
        set.alpha_cl, set.cl_list, set.cd_list;
        alpha_cd=set.alpha_cd)

    section_refs = NameRef[]
    surfaces = get(definition.extras, "surfaces", nothing)
    if !isnothing(surfaces) && !isnothing(get(surfaces, "data", nothing))
        for (si, surf_row) in enumerate(KiteGeometry.parse_table(surfaces))
            surf_name = haskey(surf_row, :name) && !isnothing(surf_row.name) ?
                Symbol(surf_row.name) : Symbol("$(name)_plate_$si")
            twist = haskey(surf_row, :twist) && !isnothing(surf_row.twist) ?
                float(surf_row.twist) : 0.0
            push!(stations, Station(
                surf_name, [cell_ref(surf_row.point_idx)], STATIC, 0.0;
                x_airf=collect(Float64, surf_row.x_airf),
                y_airf=collect(Float64, surf_row.y_airf),
                area=float(surf_row.area), twist))
            push!(section_refs, surf_name)
        end
    end

    PlateWing(name, section_refs, cl_interp, cd_interp;
              dynamics_type, transform,
              angular_damping=something(wing.angular_damping, KVec3(0.0, 150.0, 0.0)),
              drag_corr=something(wing.drag_corr, 0.93),
              z_ref_points=maybe(ref_pair_spec, wing.z_ref_points),
              y_ref_points=maybe(ref_pair_spec, wing.y_ref_points),
              origin=maybe(ref_points_spec, wing.origin_idx))
end
