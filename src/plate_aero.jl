# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

# Registered symbolic functions for flat-plate aerodynamics.
# The actual equation generation is in generate_system/plate_eqs.jl.

"""
    create_plate_interpolations(alpha_deg, cl_data, cd_data)

Create CL and CD interpolation objects from polar data vectors.
Uses linear interpolation with flat extrapolation for CL and
line extrapolation for CD, matching VortexStepMethod conventions.

# Arguments
- `alpha_deg`: angle of attack values [deg]
- `cl_data`: lift coefficient values
- `cd_data`: drag coefficient values

# Returns
- `(cl_interp, cd_interp)` tuple of interpolation objects
"""
function create_plate_interpolations(
    alpha_cl, cl_data, cd_data;
    alpha_cd=nothing
)
    alpha_cd_vec = isnothing(alpha_cd) ? alpha_cl : alpha_cd
    cl_interp = linear_interpolation(
        Vector{Float64}(alpha_cl),
        Vector{Float64}(cl_data);
        extrapolation_bc=Flat())
    cd_interp = linear_interpolation(
        Vector{Float64}(alpha_cd_vec),
        Vector{Float64}(cd_data);
        extrapolation_bc=Line())
    return (cl_interp, cd_interp)
end

# Registered symbolic accessors for PlateWing fields
for T in [PlateWing]
    ST = SystemStructure{T}
    @eval begin
        get_plate_cl(sys::$ST, wing_idx::Int64,
                     alpha_deg) =
            sys.wings[wing_idx].calc_cl(alpha_deg)
        @register_symbolic get_plate_cl(
            sys::$ST, wing_idx::Int64, alpha_deg)

        get_plate_cd(sys::$ST, wing_idx::Int64,
                     alpha_deg) =
            sys.wings[wing_idx].calc_cd(alpha_deg)
        @register_symbolic get_plate_cd(
            sys::$ST, wing_idx::Int64, alpha_deg)

        get_plate_drag_corr(sys::$ST, idx::Int64) =
            sys.wings[idx].drag_corr
        @register_symbolic get_plate_drag_corr(
            sys::$ST, idx::Int64)

        get_plate_steering(sys::$ST, idx::Int64) =
            sys.wings[idx].rel_steering
        @register_symbolic get_plate_steering(
            sys::$ST, idx::Int64)

        get_plate_depower(sys::$ST, idx::Int64) =
            sys.wings[idx].alpha_depower
        @register_symbolic get_plate_depower(
            sys::$ST, idx::Int64)

        get_surface_x_airf(sys::$ST, wing_idx::Int64,
                           surf_idx::Int64) =
            sys.wings[wing_idx].surfaces[surf_idx].x_airf
        @register_array_symbolic get_surface_x_airf(
            sys::$ST, wing_idx::Int64,
            surf_idx::Int64) begin
            size = (3,)
            eltype = SimFloat
        end

        get_surface_y_airf(sys::$ST, wing_idx::Int64,
                           surf_idx::Int64) =
            sys.wings[wing_idx].surfaces[surf_idx].y_airf
        @register_array_symbolic get_surface_y_airf(
            sys::$ST, wing_idx::Int64,
            surf_idx::Int64) begin
            size = (3,)
            eltype = SimFloat
        end

        get_surface_area(sys::$ST, wing_idx::Int64,
                         surf_idx::Int64) =
            sys.wings[wing_idx].surfaces[surf_idx].area
        @register_symbolic get_surface_area(
            sys::$ST, wing_idx::Int64, surf_idx::Int64)

        get_surface_alpha_offset(sys::$ST,
                                 wing_idx::Int64,
                                 surf_idx::Int64) =
            sys.wings[wing_idx].surfaces[surf_idx].alpha_offset
        @register_symbolic get_surface_alpha_offset(
            sys::$ST, wing_idx::Int64, surf_idx::Int64)

        get_surface_twist_a(sys::$ST, wing_idx::Int64,
                            surf_idx::Int64) =
            sys.wings[wing_idx].surfaces[surf_idx].twist_a
        @register_symbolic get_surface_twist_a(
            sys::$ST, wing_idx::Int64, surf_idx::Int64)

        get_surface_twist_b(sys::$ST, wing_idx::Int64,
                            surf_idx::Int64) =
            sys.wings[wing_idx].surfaces[surf_idx].twist_b
        @register_symbolic get_surface_twist_b(
            sys::$ST, wing_idx::Int64, surf_idx::Int64)

        get_surface_twist_c(sys::$ST, wing_idx::Int64,
                            surf_idx::Int64) =
            sys.wings[wing_idx].surfaces[surf_idx].twist_c
        @register_symbolic get_surface_twist_c(
            sys::$ST, wing_idx::Int64, surf_idx::Int64)
    end
end
