# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Registered symbolic functions for flat-plate aerodynamics.
# The actual equation generation is in generate_system/plate_eqs.jl.

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

# Registered symbolic accessors for PlateWing fields.
# Single concrete type — no duplicate method risk.
get_plate_cl(sys::SystemStructure{PlateWing},
             wing_idx::Int64, alpha_deg) =
    sys.wings[wing_idx].calc_cl(alpha_deg)
@register_symbolic get_plate_cl(
    sys::SystemStructure{PlateWing},
    wing_idx::Int64, alpha_deg)

get_plate_cd(sys::SystemStructure{PlateWing},
             wing_idx::Int64, alpha_deg) =
    sys.wings[wing_idx].calc_cd(alpha_deg)
@register_symbolic get_plate_cd(
    sys::SystemStructure{PlateWing},
    wing_idx::Int64, alpha_deg)

get_plate_drag_corr(sys::SystemStructure{PlateWing},
                    idx::Int64) =
    sys.wings[idx].drag_corr
@register_symbolic get_plate_drag_corr(
    sys::SystemStructure{PlateWing}, idx::Int64)

get_surface_x_airf(sys::SystemStructure{PlateWing},
                   wing_idx::Int64, surf_idx::Int64) =
    sys.wings[wing_idx].surfaces[surf_idx].x_airf
@register_array_symbolic get_surface_x_airf(
    sys::SystemStructure{PlateWing}, wing_idx::Int64,
    surf_idx::Int64) begin
    size = (3,)
    eltype = SimFloat
end

get_surface_y_airf(sys::SystemStructure{PlateWing},
                   wing_idx::Int64, surf_idx::Int64) =
    sys.wings[wing_idx].surfaces[surf_idx].y_airf
@register_array_symbolic get_surface_y_airf(
    sys::SystemStructure{PlateWing}, wing_idx::Int64,
    surf_idx::Int64) begin
    size = (3,)
    eltype = SimFloat
end

get_surface_area(sys::SystemStructure{PlateWing},
                 wing_idx::Int64, surf_idx::Int64) =
    sys.wings[wing_idx].surfaces[surf_idx].area
@register_symbolic get_surface_area(
    sys::SystemStructure{PlateWing},
    wing_idx::Int64, surf_idx::Int64)


get_surface_twist(sys::SystemStructure{PlateWing},
                  wing_idx::Int64, surf_idx::Int64) =
    sys.wings[wing_idx].surfaces[surf_idx].twist
@register_symbolic get_surface_twist(
    sys::SystemStructure{PlateWing},
    wing_idx::Int64, surf_idx::Int64)
