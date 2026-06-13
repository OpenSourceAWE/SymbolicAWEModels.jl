# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Point dynamics equation generation

"""
    point_eqs!(s, eqs, defaults, guesses, points, segments, twist_surfaces, wings, psys;
               R_b_to_w, wing_vel, wind_vec_gnd, twist_angle,
               pos, vel, acc, point_force, point_mass, spring_force_vec, drag_force, l0,
               spring_sum_force, point_drag_force, total_drag,
               disturb_force, tether_r, chord_b, fixed_pos, normal, pos_b,
               fix_point_sphere, fix_static, body_frame_damping, world_frame_damping,
               va_point_b, va_point_w, wind_at_point, height,
               aero_force_point_b,
               twist_surface_y_airf, tether_wing_force, tether_wing_moment)

Generate equations for all point types (STATIC, DYNAMIC, QUASI_STATIC, WING).

# Arguments
- `s::SymbolicAWEModel`: The main model object (for atmospheric model).
- `eqs`, `defaults`, `guesses`: Accumulating vectors for the MTK system.
- `points`, `segments`, `twist_surfaces`, `wings`: System components.
- `psys`: Symbolic parameter representing the system structure.
- `R_b_to_w`: Symbolic rotation matrix (body to world).
- `wing_vel`: Symbolic wing center of mass velocity.
- `wind_vec_gnd`: Symbolic ground-level wind vector.
- `twist_angle`: Symbolic twist_surface twist angle.
- `pos`, `vel`, `acc`: Pre-declared point state variables.
- `point_force`, `point_mass`: Pre-declared point force and mass variables.
- `spring_force_vec`, `drag_force`, `l0`: Pre-declared segment force variables.
- `spring_sum_force`: Pre-declared accumulated spring/drag forces variable.
- Other variables: Various point-specific symbolic variables.
- `tether_wing_force`, `tether_wing_moment`: Mutable arrays to accumulate forces/moments.

# Returns
- Tuple `(eqs, defaults, guesses)` with updated equation vectors.
  Note: `tether_wing_force` and `tether_wing_moment` are modified in-place.
"""
function point_eqs!(s, eqs, defaults, guesses, points, segments, twist_surfaces, wings, psys;
                    R_b_to_w, com_w,
                    wing_vel, wind_vec_gnd, twist_angle,
                    pos, vel, acc, point_force, point_mass, spring_force_vec, drag_force, l0,
                    spring_sum_force, point_drag_force, total_drag,
                    disturb_force, tether_r, chord_b, fixed_pos, normal, pos_b,
                    fix_point_sphere, fix_static, body_frame_damping, world_frame_damping,
                    va_point_b, va_point_w, wind_at_point, height,
                    aero_force_point_b,
                    twist_surface_y_airf, tether_wing_force, tether_wing_moment)

    for point in points
        F::Vector{Num} = zeros(Num, 3)
        seg_drag::Vector{Num} = zeros(Num, 3)
        mass = get_extra_mass(psys, point.idx)
        for segment in segments
            if point.idx in segment.point_idxs
                mass_per_meter =
                    get_density(psys, segment.idx) * π *
                    (get_diameter(psys, segment.idx) / 2)^2
                inverted = segment.point_idxs[2] == point.idx
                if inverted
                    F .-= spring_force_vec[:, segment.idx]
                else
                    F .+= spring_force_vec[:, segment.idx]
                end
                mass += mass_per_meter * l0[segment.idx] / 2
                half_seg_drag = 0.5 * drag_force[:, segment.idx]
                F .+= half_seg_drag
                seg_drag .+= half_seg_drag
            end
        end

        # The net force on the point. This variable is used by other components.
        eqs = [
            eqs
            spring_sum_force[:, point.idx] ~ F  # Store accumulated spring/drag forces
            point_mass[point.idx] ~ mass
            disturb_force[:, point.idx] ~ get_disturb(psys, point.idx)
            body_frame_damping[:, point.idx] ~ get_body_frame_damping(psys, point.idx)
            world_frame_damping[:, point.idx] ~ get_world_frame_damping(psys, point.idx)
        ]

        # Calculate apparent velocity for ALL points (needed for PARTICLE_DYNAMICS wings and generally useful)
        # Get the wing's R_b_to_w for transforming to body frame
        wing_idx_for_transform = if point.type == WING
            point.wing_idx
        elseif length(wings) > 0
            # Use first wing for non-wing points
            Int64(1)
        else
            nothing
        end

        if !isnothing(wing_idx_for_transform)
            eqs = [
                eqs
                height[point.idx] ~ max(0.0, pos[3, point.idx])
                wind_at_point[:, point.idx] ~
                    calc_wind_factor(s.am, pos[1, point.idx], pos[2, point.idx],
                                     pos[3, point.idx], psys) * wind_vec_gnd
                va_point_w[:, point.idx] ~
                    wind_at_point[:, point.idx] - vel[:, point.idx]
                va_point_b[:, point.idx] ~
                    R_b_to_w[:, :, wing_idx_for_transform]' * va_point_w[:, point.idx]
                point_drag_force[:, point.idx] ~
                    0.5 * calc_rho(s.am, height[point.idx]) *
                    get_point_drag_coeff(psys, point.idx) *
                    smooth_norm(va_point_w[:, point.idx]) *
                    get_point_area(psys, point.idx) *
                    va_point_w[:, point.idx]
            ]
        else
            # No wings - still compute wind and drag in world frame
            eqs = [
                eqs
                height[point.idx] ~ max(0.0, pos[3, point.idx])
                wind_at_point[:, point.idx] ~
                    calc_wind_factor(s.am, pos[1, point.idx], pos[2, point.idx],
                                     pos[3, point.idx], psys) * wind_vec_gnd
                va_point_w[:, point.idx] ~
                    wind_at_point[:, point.idx] - vel[:, point.idx]
                va_point_b[:, point.idx] ~ zeros(3)  # No body frame without wing
                point_drag_force[:, point.idx] ~
                    0.5 * calc_rho(s.am, height[point.idx]) *
                    get_point_drag_coeff(psys, point.idx) *
                    smooth_norm(va_point_w[:, point.idx]) *
                    get_point_area(psys, point.idx) *
                    va_point_w[:, point.idx]
            ]
        end

        # Total drag: point aero drag + share of segment drag
        eqs = [
            eqs
            total_drag[:, point.idx] ~
                point_drag_force[:, point.idx] + seg_drag
        ]

        if point.type == WING
            # Find the wing for this point
            wing = wings[point.wing_idx]

            if wing.dynamics_type == PARTICLE_DYNAMICS
                # PARTICLE_DYNAMICS wing: Points are DYNAMIC and receive lumped
                # panel/plate forces. Similar to DYNAMIC points but
                # with aero forces included.

                # Add aerodynamic forces (calculated in aero_eqs!)
                aero_force_w = R_b_to_w[:, :, wing.idx] * aero_force_point_b[:, point.idx]

                eqs = [
                    eqs
                    point_force[:, point.idx] ~
                        spring_sum_force[:, point.idx] + aero_force_w + Num[0, 0, -get_g_earth(psys) * mass] + disturb_force[:, point.idx] + point_drag_force[:, point.idx]
                ]

                # Damping terms (applied in body frame, then transformed to world frame)
                vel_diff_w = vel[:, point.idx] - wing_vel[:, point.wing_idx]
                vel_diff_b = R_b_to_w[:, :, wing.idx]' * vel_diff_w
                body_frame_damp_b = body_frame_damping[:, point.idx] .* vel_diff_b
                body_frame_damp_vec = R_b_to_w[:, :, wing.idx] * body_frame_damp_b
                world_frame_damp_vec = world_frame_damping[:, point.idx] .* vel[:, point.idx]

                # DYNAMIC point equations
                axis = smooth_normalize(pos[:, point.idx])
                eqs = [
                    eqs
                    fix_point_sphere[point.idx] ~ get_fix_point_sphere(psys, point.idx)
                    fix_static[point.idx] ~ get_fix_static(psys, point.idx)
                    D(pos[:, point.idx]) ~ ifelse.(
                        fix_static[point.idx] == true,
                        zeros(3),
                        ifelse.(fix_point_sphere[point.idx]==true,
                                vel[:, point.idx] ⋅ axis * axis,
                                vel[:, point.idx]
                        )
                    )
                    D(vel[:, point.idx]) ~ ifelse.(
                        fix_static[point.idx] == true,
                        zeros(3),
                        ifelse.(fix_point_sphere[point.idx]==true,
                                acc[:, point.idx] ⋅ axis * axis,
                                acc[:, point.idx]
                        )
                    )
                    acc[:, point.idx] ~ point_force[:, point.idx] ./ mass - body_frame_damp_vec - world_frame_damp_vec
                ]
                defaults = [
                    defaults
                    [pos[j, point.idx] => get_pos_w(psys, point.idx)[j] for j = 1:3]
                    [vel[j, point.idx] => get_vel_w(psys, point.idx)[j] for j = 1:3]
                ]

            elseif wing.dynamics_type == RIGID_DYNAMICS
                # RIGID_DYNAMICS wing: rigid body constraint
                eqs = [
                    eqs
                    point_force[:, point.idx] ~
                        spring_sum_force[:, point.idx] + Num[0, 0, -get_g_earth(psys) * mass] + disturb_force[:, point.idx] + point_drag_force[:, point.idx]
                ]

                found = 0
                twist_surface = nothing
                for twist_surface_ in twist_surfaces
                    if point.idx in twist_surface_.point_idxs
                        twist_surface = twist_surface_
                        found += 1
                    end
                end
                in_group = found == 1
                !(found in [0, 1]) && error(
                    "Kite point number $(point.idx) is part of $found twist_surfaces, " *
                    "and should be part of exactly 0 or 1 twist_surfaces.",
                )

                if found == 1
                    found = 0
                    for wing_ in wings
                        if twist_surface.idx in wing_.twist_surface_idxs
                            found += 1
                        end
                    end
                    !(found == 1) && error(
                        "Kite twist_surface number $(twist_surface.idx) is part of $found wings, " *
                        "and should be part of exactly 1 wing.",
                    )

                    eqs = [
                        eqs
                        fixed_pos[:, point.idx] ~ get_le_pos(psys, twist_surface.idx)
                        chord_b[:, point.idx] ~
                            get_pos_b(psys, point.idx) .- fixed_pos[:, point.idx]
                        normal[:, point.idx] ~ chord_b[:, point.idx] × twist_surface_y_airf[:, twist_surface.idx]
                        pos_b[:, point.idx] ~
                            fixed_pos[:, point.idx] .+
                            cos(twist_angle[twist_surface.idx]) * chord_b[:, point.idx] -
                            sin(twist_angle[twist_surface.idx]) * normal[:, point.idx]
                    ]
                elseif found == 0
                    eqs = [eqs; pos_b[:, point.idx] ~ get_pos_b(psys, point.idx)]
                end
                # Moment arm about COM (world frame)
                eqs = [
                    eqs
                    tether_r[:, point.idx] ~
                        pos[:, point.idx] -
                        com_w[:, point.wing_idx]
                ]
                # In-group (twist_surface) points can be excluded from the
                # wing moment via the wing's group_points_moment flag, while
                # their force always contributes.
                point_moment = tether_r[:, point.idx] ×
                    point_force[:, point.idx]
                if in_group
                    point_moment = ifelse.(
                        get_group_points_moment(
                            psys, point.wing_idx) == true,
                        point_moment, zeros(3))
                end
                tether_wing_moment[:, point.wing_idx] .+= point_moment
                tether_wing_force[:, point.wing_idx] .+=
                    point_force[:, point.idx]

                # Rigid body constraint: COM + R_b_to_w * pos_b
                # (pos_b is offset from COM in body frame)
                eqs = [
                    eqs
                    pos[:, point.idx] ~
                        com_w[:, point.wing_idx] +
                        R_b_to_w[:, :, point.wing_idx] *
                        pos_b[:, point.idx]
                    vel[:, point.idx] ~ zeros(3)
                    acc[:, point.idx] ~ zeros(3)
                ]
            else
                error("Unsupported dynamics_type $(wing.dynamics_type) " *
                      "for WING point $(point.idx)")
            end
        elseif point.type == STATIC
            # Define point_force for STATIC points
            eqs = [
                eqs
                point_force[:, point.idx] ~
                    spring_sum_force[:, point.idx] + Num[0, 0, -get_g_earth(psys) * mass] + disturb_force[:, point.idx] + point_drag_force[:, point.idx]
                pos[:, point.idx] ~ get_pos_w(psys, point.idx)
                vel[:, point.idx] ~ zeros(3)
                acc[:, point.idx] ~ zeros(3)
            ]
        elseif point.type == DYNAMIC
            # Define point_force for DYNAMIC points
            eqs = [
                eqs
                point_force[:, point.idx] ~
                    spring_sum_force[:, point.idx] + Num[0, 0, -get_g_earth(psys) * mass] + disturb_force[:, point.idx] + point_drag_force[:, point.idx]
            ]

            if length(wings) > 0
                # Damping applied in body frame, then transformed to world frame
                vel_diff_w = vel[:, point.idx] - wing_vel[:, point.wing_idx]
                vel_diff_b = R_b_to_w[:, :, point.wing_idx]' * vel_diff_w
                body_frame_damp_b = body_frame_damping[:, point.idx] .* vel_diff_b
                body_frame_damp_vec = R_b_to_w[:, :, point.wing_idx] * body_frame_damp_b
            else
                body_frame_damp_vec = zeros(3)
            end
            world_frame_damp_vec = world_frame_damping[:, point.idx] .* vel[:, point.idx]

            axis = smooth_normalize(pos[:, point.idx])
            eqs = [
                eqs
                fix_point_sphere[point.idx] ~ get_fix_point_sphere(psys, point.idx)
                fix_static[point.idx] ~ get_fix_static(psys, point.idx)
                D(pos[:, point.idx]) ~ ifelse.(
                    fix_static[point.idx] == true,
                    zeros(3),
                    ifelse.(fix_point_sphere[point.idx]==true,
                            vel[:, point.idx] ⋅ axis * axis,
                            vel[:, point.idx]
                    )
                )
                D(vel[:, point.idx]) ~ ifelse.(
                    fix_static[point.idx] == true,
                    zeros(3),
                    ifelse.(fix_point_sphere[point.idx]==true,
                            acc[:, point.idx] ⋅ axis * axis,
                            acc[:, point.idx]
                    )
                )
                acc[:, point.idx]    ~ point_force[:, point.idx] ./ mass - body_frame_damp_vec - world_frame_damp_vec
            ]
            defaults = [
                defaults
                [pos[j, point.idx] => get_pos_w(psys, point.idx)[j] for j = 1:3]
                [vel[j, point.idx] => get_vel_w(psys, point.idx)[j] for j = 1:3]
            ]
        elseif point.type == QUASI_STATIC
            # Define point_force for QUASI_STATIC points
            eqs = [
                eqs
                point_force[:, point.idx] ~
                    spring_sum_force[:, point.idx] + Num[0, 0, -get_g_earth(psys) * mass] + disturb_force[:, point.idx] + point_drag_force[:, point.idx]
                fix_static[point.idx] ~ get_fix_static(psys, point.idx)
                vel[:, point.idx] ~ zeros(3)
                acc[:, point.idx] ~ zeros(3)
                # When fix_static=true: position is fixed at get_pos_w
                # When fix_static=false: force must balance to zero
                ifelse.(
                    fix_static[point.idx] == true,
                    pos[:, point.idx],
                    point_force[:, point.idx]
                ) ~ ifelse.(
                    fix_static[point.idx] == true,
                    get_pos_w(psys, point.idx),
                    zeros(3)
                )
            ]
            guesses = [
                guesses
                [acc[j, point.idx] => 0 for j = 1:3]
                [pos[j, point.idx] => get_pos_w(psys, point.idx)[j] for j = 1:3]
                [point_force[j, point.idx] => 0 for j = 1:3]
            ]
        else
            error("Unknown point type: $(typeof(point))")
        end
    end

    return eqs, defaults, guesses
end
