# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Symbolic accessors: registered with SystemStructure (UnionAll)
# because @register_symbolic erases type parameters, causing
# duplicate methods if registered per concrete wing type.

# ==================== GENERIC ACCESSORS ==================== #
# These access SystemStructure fields shared by all wing types:
# points, segments, pulleys, winches, tethers, settings, BaseWing
# fields.

# ---- Points ----
get_pos_w(sys::SystemStructure, idx::Int64) =
    sys.points[idx].pos_w
@register_array_symbolic get_pos_w(
    sys::SystemStructure, idx::Int64) begin
    size = (3,)
    eltype = SimFloat
end
get_vel_w(sys::SystemStructure, idx::Int64) =
    sys.points[idx].vel_w
@register_array_symbolic get_vel_w(
    sys::SystemStructure, idx::Int64) begin
    size = (3,)
    eltype = SimFloat
end
get_pos_b(sys::SystemStructure, idx::Int64) =
    sys.points[idx].pos_b
@register_array_symbolic get_pos_b(
    sys::SystemStructure, idx::Int64) begin
    size = (3,)
    eltype = SimFloat
end
get_va_b(sys::SystemStructure, idx::Int64) =
    sys.points[idx].va_b
@register_array_symbolic get_va_b(
    sys::SystemStructure, idx::Int64) begin
    size = (3,)
    eltype = SimFloat
end
get_disturb(sys::SystemStructure, idx::Int64) =
    sys.points[idx].disturb
@register_array_symbolic get_disturb(
    sys::SystemStructure, idx::Int64) begin
    size = (3,)
    eltype = SimFloat
end
get_extra_mass(sys::SystemStructure, idx::Int64) =
    sys.points[idx].extra_mass
@register_symbolic get_extra_mass(
    sys::SystemStructure, idx::Int64)
get_body_frame_damping(sys::SystemStructure, idx::Int64) =
    sys.points[idx].body_frame_damping
@register_array_symbolic get_body_frame_damping(
    sys::SystemStructure, idx::Int64) begin
    size = (3,)
    eltype = SimFloat
end
get_world_frame_damping(sys::SystemStructure, idx::Int64) =
    sys.points[idx].world_frame_damping
@register_array_symbolic get_world_frame_damping(
    sys::SystemStructure, idx::Int64) begin
    size = (3,)
    eltype = SimFloat
end
get_point_area(sys::SystemStructure, idx::Int64) =
    sys.points[idx].area
@register_symbolic get_point_area(
    sys::SystemStructure, idx::Int64)
get_point_drag_coeff(sys::SystemStructure, idx::Int64) =
    sys.points[idx].drag_coeff
@register_symbolic get_point_drag_coeff(
    sys::SystemStructure, idx::Int64)
function get_point_aero_force(
    sys::SystemStructure, idx::Int64, component::Int
)
    point = sys.points[idx]
    if point.wing_idx > 0
        wing = sys.wings[point.wing_idx]
        wing.aero_mode == AERO_NONE && return 0.0
    end
    return point.aero_force_b[component]
end
@register_symbolic get_point_aero_force(
    sys::SystemStructure, idx::Int64, component::Int)
get_fix_point_sphere(sys::SystemStructure, idx::Int64) =
    sys.points[idx].fix_sphere
@register_symbolic get_fix_point_sphere(
    sys::SystemStructure, idx::Int64)
get_fix_static(sys::SystemStructure, idx::Int64) =
    sys.points[idx].fix_static
@register_symbolic get_fix_static(
    sys::SystemStructure, idx::Int64)

# ---- Wings (BaseWing fields via delegation) ----
get_wing_pos_w(sys::SystemStructure, idx::Int64) =
    sys.wings[idx].pos_w
@register_array_symbolic get_wing_pos_w(
    sys::SystemStructure, idx::Int64) begin
    size = (3,)
    eltype = SimFloat
end
get_wing_vel_w(sys::SystemStructure, idx::Int64) =
    sys.wings[idx].vel_w
@register_array_symbolic get_wing_vel_w(
    sys::SystemStructure, idx::Int64) begin
    size = (3,)
    eltype = SimFloat
end
get_Q_b_to_w(sys::SystemStructure, idx::Int64) =
    sys.wings[idx].Q_b_to_w
@register_array_symbolic get_Q_b_to_w(
    sys::SystemStructure, idx::Int64) begin
    size = (4,)
    eltype = SimFloat
end
get_com_w(sys::SystemStructure, idx::Int64) =
    sys.wings[idx].com_w
@register_array_symbolic get_com_w(
    sys::SystemStructure, idx::Int64) begin
    size = (3,)
    eltype = SimFloat
end
get_com_vel(sys::SystemStructure, idx::Int64) =
    sys.wings[idx].com_vel
@register_array_symbolic get_com_vel(
    sys::SystemStructure, idx::Int64) begin
    size = (3,)
    eltype = SimFloat
end
get_Q_p_to_w(sys::SystemStructure, idx::Int64) =
    sys.wings[idx].Q_p_to_w
@register_array_symbolic get_Q_p_to_w(
    sys::SystemStructure, idx::Int64) begin
    size = (4,)
    eltype = SimFloat
end
get_ω_p(sys::SystemStructure, idx::Int64) =
    sys.wings[idx].ω_p
@register_array_symbolic get_ω_p(
    sys::SystemStructure, idx::Int64) begin
    size = (3,)
    eltype = SimFloat
end
get_com_offset_b(sys::SystemStructure, idx::Int64) =
    sys.wings[idx].com_offset_b
@register_array_symbolic get_com_offset_b(
    sys::SystemStructure, idx::Int64) begin
    size = (3,)
    eltype = SimFloat
end
get_R_b_to_p(sys::SystemStructure, idx::Int64) =
    sys.wings[idx].R_b_to_p
@register_array_symbolic get_R_b_to_p(
    sys::SystemStructure, idx::Int64) begin
    size = (3, 3)
    eltype = SimFloat
end
get_wing_mass(sys::SystemStructure, idx::Int64) =
    sys.wings[idx].mass
@register_symbolic get_wing_mass(
    sys::SystemStructure, idx::Int64)
get_inertia_principal(sys::SystemStructure, idx::Int64) =
    sys.wings[idx].inertia_principal
@register_array_symbolic get_inertia_principal(
    sys::SystemStructure, idx::Int64) begin
    size = (3,)
    eltype = SimFloat
end
get_ω_b(sys::SystemStructure, idx::Int64) =
    sys.wings[idx].ω_b
@register_array_symbolic get_ω_b(
    sys::SystemStructure, idx::Int64) begin
    size = (3,)
    eltype = SimFloat
end
get_wind_disturb(sys::SystemStructure, idx::Int64) =
    sys.wings[idx].wind_disturb
@register_array_symbolic get_wind_disturb(
    sys::SystemStructure, idx::Int64) begin
    size = (3,)
    eltype = SimFloat
end
get_fix_wing_sphere(sys::SystemStructure, idx::Int64) =
    sys.wings[idx].fix_sphere
@register_symbolic get_fix_wing_sphere(
    sys::SystemStructure, idx::Int64)
get_drag_frac(sys::SystemStructure, idx::Int64) =
    sys.wings[idx].drag_frac
@register_symbolic get_drag_frac(
    sys::SystemStructure, idx::Int64)
get_y_damping(sys::SystemStructure, idx::Int64) =
    sys.wings[idx].y_damping
@register_symbolic get_y_damping(
    sys::SystemStructure, idx::Int64)
get_angular_damping(sys::SystemStructure, idx::Int64) =
    sys.wings[idx].angular_damping
@register_symbolic get_angular_damping(
    sys::SystemStructure, idx::Int64)
get_z_disturb(sys::SystemStructure, idx::Int64) =
    sys.wings[idx].z_disturb
@register_symbolic get_z_disturb(
    sys::SystemStructure, idx::Int64)

# ---- Segments ----
get_l0(sys::SystemStructure, idx::Int64) =
    sys.segments[idx].l0
@register_symbolic get_l0(
    sys::SystemStructure, idx::Int64)
get_diameter(sys::SystemStructure, idx::Int64) =
    sys.segments[idx].diameter
@register_symbolic get_diameter(
    sys::SystemStructure, idx::Int64)
get_compression_frac(sys::SystemStructure, idx::Int64) =
    sys.segments[idx].compression_frac
@register_symbolic get_compression_frac(
    sys::SystemStructure, idx::Int64)
get_unit_stiffness(sys::SystemStructure, idx::Int64) =
    sys.segments[idx].unit_stiffness
@register_symbolic get_unit_stiffness(
    sys::SystemStructure, idx::Int64)
get_unit_damping(sys::SystemStructure, idx::Int64) =
    sys.segments[idx].unit_damping
@register_symbolic get_unit_damping(
    sys::SystemStructure, idx::Int64)

# ---- Pulleys ----
get_pulley_len(sys::SystemStructure, idx::Int64) =
    sys.pulleys[idx].len
@register_symbolic get_pulley_len(
    sys::SystemStructure, idx::Int64)
get_pulley_vel(sys::SystemStructure, idx::Int64) =
    sys.pulleys[idx].vel
@register_symbolic get_pulley_vel(
    sys::SystemStructure, idx::Int64)
get_sum_len(sys::SystemStructure, idx::Int64) =
    sys.pulleys[idx].sum_len
@register_symbolic get_sum_len(
    sys::SystemStructure, idx::Int64)

# ---- Winches ----
get_set_value(sys::SystemStructure, idx::Int64) =
    sys.winches[idx].set_value
@register_symbolic get_set_value(
    sys::SystemStructure, idx::Int64)
get_winch_vel(sys::SystemStructure, idx::Int64) =
    sys.winches[idx].vel
@register_symbolic get_winch_vel(
    sys::SystemStructure, idx::Int64)
get_brake(sys::SystemStructure, idx::Int64) =
    sys.winches[idx].brake
@register_symbolic get_brake(
    sys::SystemStructure, idx::Int64)
get_speed_controlled(sys::SystemStructure, idx::Int64) =
    sys.winches[idx].speed_controlled
@register_symbolic get_speed_controlled(
    sys::SystemStructure, idx::Int64)
get_winch_gear_ratio(sys::SystemStructure, idx::Int64) =
    sys.winches[idx].gear_ratio
@register_symbolic get_winch_gear_ratio(
    sys::SystemStructure, idx::Int64)
get_winch_drum_radius(sys::SystemStructure, idx::Int64) =
    sys.winches[idx].drum_radius
@register_symbolic get_winch_drum_radius(
    sys::SystemStructure, idx::Int64)
get_winch_f_coulomb(sys::SystemStructure, idx::Int64) =
    sys.winches[idx].f_coulomb
@register_symbolic get_winch_f_coulomb(
    sys::SystemStructure, idx::Int64)
get_winch_c_vf(sys::SystemStructure, idx::Int64) =
    sys.winches[idx].c_vf
@register_symbolic get_winch_c_vf(
    sys::SystemStructure, idx::Int64)
get_winch_inertia_total(sys::SystemStructure, idx::Int64) =
    sys.winches[idx].inertia_total
@register_symbolic get_winch_inertia_total(
    sys::SystemStructure, idx::Int64)
get_winch_friction_epsilon(sys::SystemStructure, idx::Int64) =
    sys.winches[idx].friction_epsilon
@register_symbolic get_winch_friction_epsilon(
    sys::SystemStructure, idx::Int64)

# ---- Tethers ----
get_tether_len(sys::SystemStructure, idx::Int64) =
    sys.tethers[idx].len
@register_symbolic get_tether_len(
    sys::SystemStructure, idx::Int64)

# ---- Settings ----
const _ZERO_WIND_FALLBACK = KVec3(1e-10, 0.0, 0.0)

function get_wind_vec(sys::SystemStructure)
    wv = sys.set.wind_vec
    if wv[1]^2 + wv[2]^2 + wv[3]^2 < 1e-20
        return _ZERO_WIND_FALLBACK
    end
    return wv
end
@register_array_symbolic get_wind_vec(
    sys::SystemStructure) begin
    size = (3,)
    eltype = SimFloat
end
get_rho_tether(sys::SystemStructure) = sys.set.rho_tether
@register_symbolic get_rho_tether(sys::SystemStructure)
get_cd_tether(sys::SystemStructure) = sys.set.cd_tether
@register_symbolic get_cd_tether(sys::SystemStructure)
get_g_earth(sys::SystemStructure) = sys.set.g_earth
@register_symbolic get_g_earth(sys::SystemStructure)

# ---- Aero overrides ----
function get_aero_force_override(
    sys::SystemStructure, idx::Int64, component::Int
)
    wing = sys.wings[idx]
    wing.aero_mode == AERO_DIRECT &&
        return wing.aero_force_b[component]
    return 0.0
end
@register_symbolic get_aero_force_override(
    sys::SystemStructure, idx::Int64, component::Int)
function get_aero_moment_override(
    sys::SystemStructure, idx::Int64, component::Int
)
    wing = sys.wings[idx]
    wing.aero_mode == AERO_DIRECT &&
        return wing.aero_moment_b[component]
    return 0.0
end
@register_symbolic get_aero_moment_override(
    sys::SystemStructure, idx::Int64, component::Int)

# ==================== VSM-SPECIFIC ACCESSORS ==================== #
# These access VSMWing-specific fields or Group fields.
# Registered with SystemStructure (UnionAll) because
# @register_symbolic erases type parameters.

get_le_pos(sys::SystemStructure, idx::Int64) =
    sys.groups[idx].le_pos
@register_array_symbolic get_le_pos(
    sys::SystemStructure, idx::Int64) begin
    size = (3,)
    eltype = SimFloat
end
get_aero_y(sys::SystemStructure, idx::Int64, iy::Int) =
    sys.wings[idx].aero_y[iy]
@register_symbolic get_aero_y(
    sys::SystemStructure, idx::Int64, iy::Int)
get_aero_x(sys::SystemStructure, idx::Int64, ix::Int) =
    sys.wings[idx].aero_x[ix]
@register_symbolic get_aero_x(
    sys::SystemStructure, idx::Int64, ix::Int)
get_aero_jac(sys::SystemStructure, idx::Int64,
             ix::Int, iy::Int) =
    sys.wings[idx].aero_jac[ix, iy]
@register_symbolic get_aero_jac(
    sys::SystemStructure, idx::Int64,
    ix::Int, iy::Int)
get_twist(sys::SystemStructure, idx::Int64) =
    sys.groups[idx].twist
@register_symbolic get_twist(
    sys::SystemStructure, idx::Int64)
get_group_damping(sys::SystemStructure, idx::Int64) =
    sys.groups[idx].damping
@register_symbolic get_group_damping(
    sys::SystemStructure, idx::Int64)
get_twist_ω(sys::SystemStructure, idx::Int64) =
    sys.groups[idx].twist_ω
@register_symbolic get_twist_ω(
    sys::SystemStructure, idx::Int64)
get_group_y_airf(sys::SystemStructure, idx::Int64) =
    sys.groups[idx].y_airf
@register_array_symbolic get_group_y_airf(
    sys::SystemStructure, idx::Int64) begin
    size = (3,)
    eltype = SimFloat
end
get_group_chord(sys::SystemStructure, idx::Int64) =
    sys.groups[idx].chord
@register_array_symbolic get_group_chord(
    sys::SystemStructure, idx::Int64) begin
    size = (3,)
    eltype = SimFloat
end
get_group_le_pos(sys::SystemStructure, idx::Int64) =
    sys.groups[idx].le_pos
@register_array_symbolic get_group_le_pos(
    sys::SystemStructure, idx::Int64) begin
    size = (3,)
    eltype = SimFloat
end
get_moment_frac(sys::SystemStructure, idx::Int64) =
    sys.groups[idx].moment_frac
@register_symbolic get_moment_frac(
    sys::SystemStructure, idx::Int64)
function get_group_moment_override(
    sys::SystemStructure, wing_idx::Int64,
    group_idx::Int64
)
    wing = sys.wings[wing_idx]
    wing.aero_mode == AERO_DIRECT &&
        return sys.groups[group_idx].aero_moment
    return 0.0
end
@register_symbolic get_group_moment_override(
    sys::SystemStructure, wing_idx::Int64,
    group_idx::Int64)
