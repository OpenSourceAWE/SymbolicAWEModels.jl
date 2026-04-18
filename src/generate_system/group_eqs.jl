# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

# Group twist dynamics equation generation

"""
    group_eqs!(eqs, defaults, guesses, groups, wings, psys;
               R_b_to_w, fix_wing, twist_angle, twist_ω, group_aero_moment,
               point_force, tether_wing_moment, group_y_airf, group_chord, group_le_pos)

Generate equations for deformable wing group twist dynamics.

# Arguments
- `eqs`, `defaults`, `guesses`: Accumulating vectors for the MTK system.
- `groups`: Collection of Group objects (deformable wing sections).
- `wings`: Collection of Wing objects.
- `psys`: Symbolic parameter representing the system structure.
- `R_b_to_w`: Symbolic rotation matrix (body to world).
- `fix_wing`: Symbolic boolean for fixing wing dynamics.
- `twist_angle`, `twist_ω`: Symbolic twist state variables.
- `group_aero_moment`: Symbolic aerodynamic moment on groups.
- `point_force`: Symbolic point force variable.
- `tether_wing_moment`: Accumulated tether moments on wings (for validation).
- `group_y_airf`, `group_chord`, `group_le_pos`: Symbolic group geometry variables.

# Returns
- Tuple `(eqs, defaults, guesses)` with updated equation vectors.
"""
function group_eqs!(eqs, defaults, guesses, groups, wings, psys;
                    R_b_to_w, fix_wing, twist_angle, twist_ω, group_aero_moment,
                    point_force, tether_wing_moment, group_y_airf, group_chord, group_le_pos)

    length(groups) == 0 && return eqs, defaults, guesses

    @variables begin
        trailing_edge_angle(t)[eachindex(groups)]
        trailing_edge_ω(t)[eachindex(groups)]
        trailing_edge_α(t)[eachindex(groups)]
        free_twist_angle(t)[eachindex(groups)]
        twist_α(t)[eachindex(groups)]
        group_tether_force(t)[eachindex(groups)]
        group_tether_moment(t)[eachindex(groups)]
        tether_force(t)[eachindex(groups[1].point_idxs), eachindex(groups)]
        tether_moment(t)[eachindex(groups[1].point_idxs), eachindex(groups)]
        r_group(t)[eachindex(groups[1].point_idxs), eachindex(groups)]
        r_vec(t)[1:3, eachindex(groups[1].point_idxs), eachindex(groups)]
    end

    for group in groups
        found = 0
        wing = nothing
        for wing_ in wings
            if group.idx in wing_.group_idxs
                wing = wing_
                found += 1
            end
        end
        !(found == 1) && error(
            "Kite group $(group.idx) is in $found wings; must be in exactly 1.",
        )

        all(iszero.(tether_wing_moment[:, wing.idx])) && error(
            "Tether wing moment is zero. At least one wing connection point " *
            "should not be part of a deforming group.",
        )

        # Set group geometry from getters (allows runtime updates)
        eqs = [
            eqs
            group_y_airf[:, group.idx] ~ get_group_y_airf(psys, group.idx)
            group_chord[:, group.idx] ~ get_group_chord(psys, group.idx)
            group_le_pos[:, group.idx] ~ get_group_le_pos(psys, group.idx)
        ]

        gc = collect(group_chord[:, group.idx])
        x_airf = smooth_normalize(gc)
        gy = collect(group_y_airf[:, group.idx])
        init_z_airf = x_airf × gy
        z_airf = sin(twist_angle[group.idx]) * x_airf + cos(twist_angle[group.idx]) * init_z_airf
        Rbw = collect(R_b_to_w[:, :, wing.idx])
        Rz = Rbw * (-1 * z_airf)  # Note: -z_airf has a bug, use -1 * z_airf instead
        gl = collect(group_le_pos[:, group.idx])

        for (i, point_idx) in enumerate(group.point_idxs)
            pf = collect(point_force[:, point_idx])
            rv = collect(r_vec[:, i, group.idx])
            pos_offset = collect(
                get_pos_b(psys, point_idx) .-
                (gl + get_moment_frac(psys, group.idx) * gc)
            )
            eqs = [
                eqs
                [r_vec[j, i, group.idx] ~ pos_offset[j]
                 for j in 1:3]
                r_group[i, group.idx] ~ rv ⋅ smooth_normalize(gc)
                tether_force[i, group.idx] ~ pf ⋅ Rz
                tether_moment[i, group.idx] ~ r_group[i, group.idx] * tether_force[i, group.idx]
            ]
        end

        # Inertia of a thin rectangular plate rotating around one edge
        # I = 1/3 × m × L² where m is total mass of group points
        group_chord = collect(group_chord)
        group_mass = sum(get_extra_mass(psys, point_idx) for point_idx in group.point_idxs)
        inertia = 1 / 3 * group_mass * smooth_norm(group_chord[:, group.idx])^2
        max_twist = deg2rad(90)

        eqs = [
            eqs
            group_tether_force[group.idx] ~ sum(tether_force[:, group.idx])
            group_tether_moment[group.idx] ~ sum(tether_moment[:, group.idx])
            twist_α[group.idx] ~
                (group_aero_moment[group.idx] + group_tether_moment[group.idx]) /
                inertia
            twist_angle[group.idx] ~
                clamp(free_twist_angle[group.idx], -max_twist, max_twist)
        ]
        if group.type == DYNAMIC
            eqs = [
                eqs
                D(free_twist_angle[group.idx]) ~
                    ifelse(fix_wing == true, 0, twist_ω[group.idx])
                D(twist_ω[group.idx]) ~ ifelse(
                    fix_wing == true,
                    0,
                    twist_α[group.idx] -
                    get_group_damping(psys, group.idx) * twist_ω[group.idx],
                )
            ]
            defaults = [
                defaults
                free_twist_angle[group.idx] => get_twist(psys, group.idx)
                twist_ω[group.idx] => get_twist_ω(psys, group.idx)
            ]
        elseif group.type == QUASI_STATIC
            eqs = [eqs; twist_ω[group.idx] ~ 0; twist_α[group.idx] ~ 0]
            guesses = [
                guesses
                free_twist_angle[group.idx] => 0
                twist_angle[group.idx] => 0
            ]
        else
            error("Wrong group type.")
        end
    end

    return eqs, defaults, guesses
end
