# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# AeroPlate: flat-plate CL/CD-lookup aerodynamics (PARTICLE_DYNAMICS only). Each
# WING point is a 1-point FIXED TwistSurface; per-point forces are built
# symbolically from the section's twisted axes, apparent wind, and polar lookups.

"""
    AeroPlate(calc_cl, calc_cd; drag_corr=1.0)

Flat-plate CL/CD lookup aerodynamics. Carries the shared polar lookups
(`calc_cl`/`calc_cd`: `α_deg → coefficient`) and drag correction used by all of
a wing's flat-plate (1-point `FIXED`) [`TwistSurface`](@ref)s. One polar set per
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

get_plate_cl(sys::SystemStructure, wing_idx::Int64, alpha_deg) =
    sys.wings[wing_idx].aero.calc_cl(alpha_deg)
@register_symbolic get_plate_cl(
    sys::SystemStructure, wing_idx::Int64, alpha_deg)

get_plate_cd(sys::SystemStructure, wing_idx::Int64, alpha_deg) =
    sys.wings[wing_idx].aero.calc_cd(alpha_deg)
@register_symbolic get_plate_cd(
    sys::SystemStructure, wing_idx::Int64, alpha_deg)

get_plate_drag_corr(sys::SystemStructure, idx::Int64) =
    sys.wings[idx].aero.drag_corr
@register_symbolic get_plate_drag_corr(
    sys::SystemStructure, idx::Int64)

get_twist_surface_area(sys::SystemStructure, idx::Int64) =
    sys.twist_surfaces[idx].area
@register_symbolic get_twist_surface_area(
    sys::SystemStructure, idx::Int64)

# ==================== equation builder ==================== #

"""
    aero_component(::AeroPlate, sys_struct, wing_idx; name)

Flat-plate aero component. Uses the same `PARTICLE_DYNAMICS` connector contract as
the other per-point modes; it is the only one that consumes the `va`/`rho` inputs
(the VSM particle modes read frozen forces and ignore them). Each WING point is a
1-point `FIXED` [`TwistSurface`](@ref) section; the per-point force is computed
from the section's twisted body-frame axes, the point's apparent wind, and its
air density.
"""
function aero_component(::AeroPlate, sys_struct, wing_idx; name)
    psys = system_struct_param(sys_struct)
    wing = sys_struct.wings[wing_idx]

    twist_surfaces = sys_struct.twist_surfaces
    points = wing_points(sys_struct, wing)
    num_points = length(points)
    connectors = particle_aero_connectors(num_points)

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

    return System(eqs, t, particle_unknowns(connectors), [psys]; name)
end

# ==================== log-point hooks ==================== #

"""
    plate_corners(twist_surface, point_pos_w, R_b_to_w) -> NTuple{4, Vector}

World-frame corners of a flat-plate section's display quad. The section's
structural point sits at quarter chord; the quad is a square of side
`sqrt(area)` spanned by the twisted chord direction and `y_airf`, so the quad
area matches the section's `area`.
"""
function plate_corners(twist_surface, point_pos_w, R_b_to_w)
    x_airf = normalize(twist_surface.chord)
    y_airf = twist_surface.y_airf
    twist = twist_surface.twist
    x_twisted = cos(twist) * x_airf + sin(twist) * (y_airf × x_airf)
    side = sqrt(twist_surface.area)
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
n_aero_log_points(::AeroPlate, wing) = 4 * length(wing.twist_surface_idxs)

"""
    write_aero_log_points!(::AeroPlate, wing, sys_struct, sys_state,
                           point_idx, zoom) -> Int

Log each flat-plate section's display quad ([`plate_corners`](@ref)).
"""
function write_aero_log_points!(::AeroPlate, wing, sys_struct, sys_state,
                                point_idx, zoom)
    R_b_to_w = wing.R_b_to_w::Matrix{SimFloat}
    for twist_surface_idx in wing.twist_surface_idxs
        twist_surface = sys_struct.twist_surfaces[twist_surface_idx]
        point = sys_struct.points[twist_surface.point_idxs[1]]
        for corner_w in plate_corners(twist_surface, point.pos_w, R_b_to_w)
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

function load_wing(mode::AeroPlate, row, idx, data, set, wing_type, vsm_set,
                   yaml_to_ref, yaml_parse_ref_points, yaml_parse_origin,
                   twist_surfaces)
    return load_plate_wing(row, idx, data, set, wing_type, mode,
        yaml_to_ref, yaml_parse_ref_points, yaml_parse_origin, twist_surfaces)
end

"""
    load_plate_wing(row, idx, data, set, wing_type, aero_mode,
                    yaml_to_ref, yaml_parse_ref_points,
                    yaml_parse_origin, twist_surfaces)

Load a flat-plate wing from a YAML wing row + `surfaces` block. Each surface
becomes a 1-point `FIXED` [`TwistSurface`](@ref) appended to `twist_surfaces`; the
wing references them by name. CL/CD interpolations come from `Settings` polar data.
"""
function load_plate_wing(row, idx, data, set, wing_type, aero_mode,
                         yaml_to_ref, yaml_parse_ref_points,
                         yaml_parse_origin, twist_surfaces)
    name = if haskey(row, :name) && !isnothing(row.name)
        Symbol(row.name)
    else
        idx
    end

    cl_interp, cd_interp = create_plate_interpolations(
        set.alpha_cl, set.cl_list, set.cd_list;
        alpha_cd=set.alpha_cd)

    drag_corr = hasfield(typeof(row), :drag_corr) &&
        !isnothing(row.drag_corr) ? float(row.drag_corr) : 0.93
    y_damping = hasfield(typeof(row), :y_damping) &&
        !isnothing(row.y_damping) ?
        float(row.y_damping) : 150.0

    z_ref = yaml_parse_ref_points(row, :z_ref_points)
    y_ref = yaml_parse_ref_points(row, :y_ref_points)
    origin = yaml_parse_origin(row, :origin_idx)
    transform = if hasfield(typeof(row), :transform_idx) &&
                   !isnothing(row.transform_idx)
        yaml_to_ref(row.transform_idx)
    else
        nothing
    end

    section_refs = NameRef[]
    if haskey(data, "surfaces") &&
       haskey(data["surfaces"], "data") &&
       data["surfaces"]["data"] !== nothing
        surf_rows = parse_table(data["surfaces"])
        for (si, surf_row) in enumerate(surf_rows)
            surf_name = haskey(surf_row, :name) && !isnothing(surf_row.name) ?
                Symbol(surf_row.name) : Symbol("$(name)_plate_$si")
            x_airf = collect(Float64, surf_row.x_airf)
            y_airf = collect(Float64, surf_row.y_airf)
            area = float(surf_row.area)
            point = yaml_to_ref(surf_row.point_idx)
            twist = hasfield(typeof(surf_row), :twist) &&
                !isnothing(surf_row.twist) ?
                float(surf_row.twist) : 0.0
            push!(twist_surfaces, TwistSurface(
                surf_name, [point], FIXED, 0.0;
                x_airf, y_airf, area, twist))
            push!(section_refs, surf_name)
        end
    end

    PlateWing(name, section_refs, cl_interp, cd_interp;
              dynamics_type=wing_type, transform, y_damping,
              drag_corr, z_ref_points=z_ref, y_ref_points=y_ref,
              origin)
end
