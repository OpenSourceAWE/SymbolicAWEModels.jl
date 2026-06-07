# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Polar (CL/CD) interpolation construction for flat-plate aerodynamics.
# The equation generation is in generate_system/aero_components.jl
# (`default_aero_plate`).

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

# Flat-plate aerodynamics is now expressed as polar `Group`s (see
# `plate_group`) and computed by the `default_aero_plate` component builder
# (src/generate_system/aero_components.jl). Group-level polar accessors
# (`get_group_cl` / `get_group_cd` / `get_group_area` / `get_group_drag_corr`)
# live in src/generate_system/accessors.jl.
