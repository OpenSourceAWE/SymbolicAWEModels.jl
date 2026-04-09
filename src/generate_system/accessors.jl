# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

# Symbolic accessor functions for ModelingToolkit integration
#
# The following `get_*` functions are registered for use within ModelingToolkit.jl's
# symbolic context. They act as symbolic "shims" or placeholders that allow the
# equation generation logic to access fields from the `SystemStructure` (`psys`) and
# `Settings` (`pset`) parameter objects during the symbolic construction of the model.
# This avoids hard-coding parameters into the equations, allowing them to be changed
# at simulation time.

get_pos_w(sys::SystemStructure{VSMWing}, idx::Int64) = sys.points[idx].pos_w
@register_array_symbolic get_pos_w(sys::SystemStructure{VSMWing}, idx::Int64) begin
    size = (3,)
    eltype = SimFloat
end
get_vel_w(sys::SystemStructure{VSMWing}, idx::Int64) = sys.points[idx].vel_w
@register_array_symbolic get_vel_w(sys::SystemStructure{VSMWing}, idx::Int64) begin
    size = (3,)
    eltype = SimFloat
end
get_pos_b(sys::SystemStructure{VSMWing}, idx::Int64) = sys.points[idx].pos_b
@register_array_symbolic get_pos_b(sys::SystemStructure{VSMWing}, idx::Int64) begin
    size = (3,)
    eltype = SimFloat
end
get_va_b(sys::SystemStructure{VSMWing}, idx::Int64) = sys.points[idx].va_b
@register_array_symbolic get_va_b(sys::SystemStructure{VSMWing}, idx::Int64) begin
    size = (3,)
    eltype = SimFloat
end
get_wing_pos_w(sys::SystemStructure{VSMWing}, idx::Int64) = sys.wings[idx].pos_w
@register_array_symbolic get_wing_pos_w(sys::SystemStructure{VSMWing}, idx::Int64) begin
    size = (3,)
    eltype = SimFloat
end
get_wing_vel_w(sys::SystemStructure{VSMWing}, idx::Int64) = sys.wings[idx].vel_w
@register_array_symbolic get_wing_vel_w(sys::SystemStructure{VSMWing}, idx::Int64) begin
    size = (3,)
    eltype = SimFloat
end
get_Q_b_to_w(sys::SystemStructure{VSMWing}, idx::Int64) = sys.wings[idx].Q_b_to_w
@register_array_symbolic get_Q_b_to_w(sys::SystemStructure{VSMWing}, idx::Int64) begin
    size = (4,)
    eltype = SimFloat
end
# Principal frame accessors
get_com_w(sys::SystemStructure{VSMWing}, idx::Int64) =
    sys.wings[idx].com_w
@register_array_symbolic get_com_w(
    sys::SystemStructure{VSMWing}, idx::Int64) begin
    size = (3,)
    eltype = SimFloat
end
get_com_vel(sys::SystemStructure{VSMWing}, idx::Int64) =
    sys.wings[idx].com_vel
@register_array_symbolic get_com_vel(
    sys::SystemStructure{VSMWing}, idx::Int64) begin
    size = (3,)
    eltype = SimFloat
end
get_Q_p_to_w(sys::SystemStructure{VSMWing}, idx::Int64) =
    sys.wings[idx].Q_p_to_w
@register_array_symbolic get_Q_p_to_w(
    sys::SystemStructure{VSMWing}, idx::Int64) begin
    size = (4,)
    eltype = SimFloat
end
get_ω_p(sys::SystemStructure{VSMWing}, idx::Int64) =
    sys.wings[idx].ω_p
@register_array_symbolic get_ω_p(
    sys::SystemStructure{VSMWing}, idx::Int64) begin
    size = (3,)
    eltype = SimFloat
end
get_com_offset_b(sys::SystemStructure{VSMWing}, idx::Int64) =
    sys.wings[idx].com_offset_b
@register_array_symbolic get_com_offset_b(
    sys::SystemStructure{VSMWing}, idx::Int64) begin
    size = (3,)
    eltype = SimFloat
end
get_R_b_to_p(sys::SystemStructure{VSMWing}, idx::Int64) =
    sys.wings[idx].R_b_to_p
@register_array_symbolic get_R_b_to_p(
    sys::SystemStructure{VSMWing}, idx::Int64) begin
    size = (3, 3)
    eltype = SimFloat
end
get_wing_mass(sys::SystemStructure{VSMWing}, idx::Int64) = sys.wings[idx].mass
@register_symbolic get_wing_mass(sys::SystemStructure{VSMWing}, idx::Int64)
get_inertia_principal(sys::SystemStructure{VSMWing}, idx::Int64) =
    sys.wings[idx].inertia_principal
@register_array_symbolic get_inertia_principal(
    sys::SystemStructure{VSMWing}, idx::Int64) begin
    size = (3,)
    eltype = SimFloat
end
get_ω_b(sys::SystemStructure{VSMWing}, idx::Int64) = sys.wings[idx].ω_b
@register_array_symbolic get_ω_b(sys::SystemStructure{VSMWing}, idx::Int64) begin
    size = (3,)
    eltype = SimFloat
end
get_wind_disturb(sys::SystemStructure{VSMWing}, idx::Int64) = sys.wings[idx].wind_disturb
@register_array_symbolic get_wind_disturb(sys::SystemStructure{VSMWing}, idx::Int64) begin
    size = (3,)
    eltype = SimFloat
end
get_disturb(sys::SystemStructure{VSMWing}, idx::Int64) = sys.points[idx].disturb
@register_array_symbolic get_disturb(sys::SystemStructure{VSMWing}, idx::Int64) begin
    size = (3,)
    eltype = SimFloat
end
get_le_pos(sys::SystemStructure{VSMWing}, idx::Int64) = sys.groups[idx].le_pos
@register_array_symbolic get_le_pos(sys::SystemStructure{VSMWing}, idx::Int64) begin
    size = (3,)
    eltype = SimFloat
end
get_vsm_y(sys::SystemStructure{VSMWing}, idx::Int64, iy::Int) =
    sys.wings[idx].vsm_y[iy]
@register_symbolic get_vsm_y(
    sys::SystemStructure{VSMWing}, idx::Int64, iy::Int)
get_vsm_x(sys::SystemStructure{VSMWing}, idx::Int64, ix::Int) =
    sys.wings[idx].vsm_x[ix]
@register_symbolic get_vsm_x(
    sys::SystemStructure{VSMWing}, idx::Int64, ix::Int)
get_vsm_jac(sys::SystemStructure{VSMWing}, idx::Int64,
             ix::Int, iy::Int) =
    sys.wings[idx].vsm_jac[ix, iy]
@register_symbolic get_vsm_jac(
    sys::SystemStructure{VSMWing}, idx::Int64,
    ix::Int, iy::Int)
get_pulley_len(sys::SystemStructure{VSMWing}, idx::Int64) = sys.pulleys[idx].len
@register_symbolic get_pulley_len(sys::SystemStructure{VSMWing}, idx::Int64)
get_pulley_vel(sys::SystemStructure{VSMWing}, idx::Int64) = sys.pulleys[idx].vel
@register_symbolic get_pulley_vel(sys::SystemStructure{VSMWing}, idx::Int64)
get_set_value(sys::SystemStructure{VSMWing}, idx::Int64) = sys.winches[idx].set_value
@register_symbolic get_set_value(sys::SystemStructure{VSMWing}, idx::Int64)
get_twist(sys::SystemStructure{VSMWing}, idx::Int64) = sys.groups[idx].twist
@register_symbolic get_twist(sys::SystemStructure{VSMWing}, idx::Int64)
get_group_damping(sys::SystemStructure{VSMWing}, idx::Int64) = sys.groups[idx].damping
@register_symbolic get_group_damping(sys::SystemStructure{VSMWing}, idx::Int64)
get_twist_ω(sys::SystemStructure{VSMWing}, idx::Int64) = sys.groups[idx].twist_ω
@register_symbolic get_twist_ω(sys::SystemStructure{VSMWing}, idx::Int64)
get_group_y_airf(sys::SystemStructure{VSMWing}, idx::Int64) = sys.groups[idx].y_airf
@register_array_symbolic get_group_y_airf(sys::SystemStructure{VSMWing}, idx::Int64) begin
    size = (3,)
    eltype = SimFloat
end
get_group_chord(sys::SystemStructure{VSMWing}, idx::Int64) = sys.groups[idx].chord
@register_array_symbolic get_group_chord(sys::SystemStructure{VSMWing}, idx::Int64) begin
    size = (3,)
    eltype = SimFloat
end
get_group_le_pos(sys::SystemStructure{VSMWing}, idx::Int64) = sys.groups[idx].le_pos
@register_array_symbolic get_group_le_pos(sys::SystemStructure{VSMWing}, idx::Int64) begin
    size = (3,)
    eltype = SimFloat
end
get_extra_mass(sys::SystemStructure{VSMWing}, idx::Int64) = sys.points[idx].extra_mass
@register_symbolic get_extra_mass(sys::SystemStructure{VSMWing}, idx::Int64)
get_l0(sys::SystemStructure{VSMWing}, idx::Int64) = sys.segments[idx].l0
@register_symbolic get_l0(sys::SystemStructure{VSMWing}, idx::Int64)
get_diameter(sys::SystemStructure{VSMWing}, idx::Int64) = sys.segments[idx].diameter
@register_symbolic get_diameter(sys::SystemStructure{VSMWing}, idx::Int64)
get_compression_frac(sys::SystemStructure{VSMWing}, idx::Int64) =
    sys.segments[idx].compression_frac
@register_symbolic get_compression_frac(sys::SystemStructure{VSMWing}, idx::Int64)
get_moment_frac(sys::SystemStructure{VSMWing}, idx::Int64) = sys.groups[idx].moment_frac
@register_symbolic get_moment_frac(sys::SystemStructure{VSMWing}, idx::Int64)
get_sum_len(sys::SystemStructure{VSMWing}, idx::Int64) = sys.pulleys[idx].sum_len
@register_symbolic get_sum_len(sys::SystemStructure{VSMWing}, idx::Int64)
get_tether_len(sys::SystemStructure{VSMWing}, idx::Int64) = sys.tethers[idx].len
@register_symbolic get_tether_len(sys::SystemStructure{VSMWing}, idx::Int64)
get_winch_vel(sys::SystemStructure{VSMWing}, idx::Int64) = sys.winches[idx].vel
@register_symbolic get_winch_vel(sys::SystemStructure{VSMWing}, idx::Int64)
get_unit_stiffness(sys::SystemStructure{VSMWing}, idx::Int64) =
    sys.segments[idx].unit_stiffness
@register_symbolic get_unit_stiffness(sys::SystemStructure{VSMWing}, idx::Int64)
get_unit_damping(sys::SystemStructure{VSMWing}, idx::Int64) =
    sys.segments[idx].unit_damping
@register_symbolic get_unit_damping(sys::SystemStructure{VSMWing}, idx::Int64)
get_body_frame_damping(sys::SystemStructure{VSMWing}, idx::Int64) = sys.points[idx].body_frame_damping
@register_array_symbolic get_body_frame_damping(sys::SystemStructure{VSMWing}, idx::Int64) begin
    size = (3,)
    eltype = SimFloat
end
get_world_frame_damping(sys::SystemStructure{VSMWing}, idx::Int64) = sys.points[idx].world_frame_damping
@register_array_symbolic get_world_frame_damping(sys::SystemStructure{VSMWing}, idx::Int64) begin
    size = (3,)
    eltype = SimFloat
end
get_point_area(sys::SystemStructure{VSMWing}, idx::Int64) = sys.points[idx].area
@register_symbolic get_point_area(sys::SystemStructure{VSMWing}, idx::Int64)
get_point_drag_coeff(sys::SystemStructure{VSMWing}, idx::Int64) = sys.points[idx].drag_coeff
@register_symbolic get_point_drag_coeff(sys::SystemStructure{VSMWing}, idx::Int64)
function get_point_aero_force(
    sys::SystemStructure{VSMWing}, idx::Int64, component::Int
)
    point = sys.points[idx]
    if point.wing_idx > 0
        wing = sys.wings[point.wing_idx]
        wing.aero_mode == AERO_NONE && return 0.0
    end
    return point.aero_force_b[component]
end
@register_symbolic get_point_aero_force(
    sys::SystemStructure{VSMWing}, idx::Int64, component::Int)
get_brake(sys::SystemStructure{VSMWing}, idx::Int64) = sys.winches[idx].brake
@register_symbolic get_brake(sys::SystemStructure{VSMWing}, idx::Int64)
get_speed_controlled(sys::SystemStructure{VSMWing}, idx::Int64) =
    sys.winches[idx].speed_controlled
@register_symbolic get_speed_controlled(
    sys::SystemStructure{VSMWing}, idx::Int64)
get_fix_point_sphere(sys::SystemStructure{VSMWing}, idx::Int64) =
    sys.points[idx].fix_sphere
@register_symbolic get_fix_point_sphere(sys::SystemStructure{VSMWing}, idx::Int64)
get_fix_static(sys::SystemStructure{VSMWing}, idx::Int64) =
    sys.points[idx].fix_static
@register_symbolic get_fix_static(sys::SystemStructure{VSMWing}, idx::Int64)
get_fix_wing_sphere(sys::SystemStructure{VSMWing}, idx::Int64) = sys.wings[idx].fix_sphere
@register_symbolic get_fix_wing_sphere(sys::SystemStructure{VSMWing}, idx::Int64)
get_drag_frac(sys::SystemStructure{VSMWing}, idx::Int64) = sys.wings[idx].drag_frac
@register_symbolic get_drag_frac(sys::SystemStructure{VSMWing}, idx::Int64)
get_y_damping(sys::SystemStructure{VSMWing}, idx::Int64) = sys.wings[idx].y_damping
@register_symbolic get_y_damping(sys::SystemStructure{VSMWing}, idx::Int64)
get_angular_damping(sys::SystemStructure{VSMWing}, idx::Int64) =
    sys.wings[idx].angular_damping
@register_symbolic get_angular_damping(sys::SystemStructure{VSMWing}, idx::Int64)
get_z_disturb(sys::SystemStructure{VSMWing}, idx::Int64) = sys.wings[idx].z_disturb
@register_symbolic get_z_disturb(sys::SystemStructure{VSMWing}, idx::Int64)
get_winch_gear_ratio(sys::SystemStructure{VSMWing}, idx::Int64) =
    sys.winches[idx].gear_ratio
@register_symbolic get_winch_gear_ratio(sys::SystemStructure{VSMWing}, idx::Int64)
get_winch_drum_radius(sys::SystemStructure{VSMWing}, idx::Int64) =
    sys.winches[idx].drum_radius
@register_symbolic get_winch_drum_radius(sys::SystemStructure{VSMWing}, idx::Int64)
get_winch_f_coulomb(sys::SystemStructure{VSMWing}, idx::Int64) =
    sys.winches[idx].f_coulomb
@register_symbolic get_winch_f_coulomb(sys::SystemStructure{VSMWing}, idx::Int64)
get_winch_c_vf(sys::SystemStructure{VSMWing}, idx::Int64) = sys.winches[idx].c_vf
@register_symbolic get_winch_c_vf(sys::SystemStructure{VSMWing}, idx::Int64)
get_winch_inertia_total(sys::SystemStructure{VSMWing}, idx::Int64) =
    sys.winches[idx].inertia_total
@register_symbolic get_winch_inertia_total(sys::SystemStructure{VSMWing}, idx::Int64)
get_winch_friction_epsilon(sys::SystemStructure{VSMWing}, idx::Int64) =
    sys.winches[idx].friction_epsilon
@register_symbolic get_winch_friction_epsilon(
    sys::SystemStructure{VSMWing}, idx::Int64)
get_wind_elevation(sys::SystemStructure{VSMWing}) = sys.wind_elevation
@register_symbolic get_wind_elevation(sys::SystemStructure{VSMWing})
get_rho_tether(set::Settings) = set.rho_tether
@register_symbolic get_rho_tether(set::Settings)
get_cd_tether(set::Settings) = set.cd_tether
@register_symbolic get_cd_tether(set::Settings)
get_v_wind(set::Settings) = set.v_wind
@register_symbolic get_v_wind(set::Settings)
get_upwind_dir(set::Settings) = set.upwind_dir
@register_symbolic get_upwind_dir(set::Settings)
get_g_earth(set::Settings) = set.g_earth
@register_symbolic get_g_earth(set::Settings)

# ==================== AERO MODE FUNCTIONS ==================== #
# Aero mode is a build-time decision (part of model SHA hash).
# For AERO_LINEARIZED: symbolic linearization equations in vsm_eqs.jl.
# For AERO_DIRECT: override functions below read stored forces.
# For AERO_NONE: equations are zeros (no registered functions needed).

"""
    get_aero_force_override(sys, wing_idx, component)

Return override aero force for non-LINEARIZED modes:
- `AERO_DIRECT`: stored `wing.aero_force_b[c]`
- `AERO_NONE`:   `0.0`
"""
function get_aero_force_override(
    sys::SystemStructure{VSMWing}, idx::Int64, component::Int
)
    wing = sys.wings[idx]
    wing.aero_mode == AERO_DIRECT &&
        return wing.aero_force_b[component]
    return 0.0
end
@register_symbolic get_aero_force_override(
    sys::SystemStructure{VSMWing}, idx::Int64, component::Int)

"""
    get_aero_moment_override(sys, wing_idx, component)

Return override aero moment for non-LINEARIZED modes.
"""
function get_aero_moment_override(
    sys::SystemStructure{VSMWing}, idx::Int64, component::Int
)
    wing = sys.wings[idx]
    wing.aero_mode == AERO_DIRECT &&
        return wing.aero_moment_b[component]
    return 0.0
end
@register_symbolic get_aero_moment_override(
    sys::SystemStructure{VSMWing}, idx::Int64, component::Int)

"""
    get_group_moment_override(sys, wing_idx, group_idx)

Return override group aero moment for non-LINEARIZED modes.
"""
function get_group_moment_override(
    sys::SystemStructure{VSMWing}, wing_idx::Int64,
    group_idx::Int64
)
    wing = sys.wings[wing_idx]
    wing.aero_mode == AERO_DIRECT &&
        return sys.groups[group_idx].aero_moment
    return 0.0
end
@register_symbolic get_group_moment_override(
    sys::SystemStructure{VSMWing}, wing_idx::Int64,
    group_idx::Int64)
