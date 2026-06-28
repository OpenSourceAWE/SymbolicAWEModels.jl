# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Segment spring-damper equation generation

"""
    segment_eqs!(s, eqs, points, segments, pulleys, tethers, bodies, params;
                 pos, vel, wind_vec_gnd, spring_force_vec, drag_force, l0,
                 pulley_len, tether_len)

Generate equations for segment spring-damper forces and aerodynamic drag.

# Arguments
- `s::SymbolicAWEModel`: The main model object (for atmospheric model).
- `eqs`: Accumulating equation vector for the MTK system.
- `points`, `segments`, `pulleys`, `tethers`, `bodies`: System components.
- `pos`, `vel`: Symbolic point state variables.
- `wind_vec_gnd`: Symbolic ground-level wind vector.
- `spring_force_vec`, `drag_force`, `l0`: Pre-declared segment force variables.
- `pulley_len`, `tether_len`: Symbolic state variables for pulley and tether lengths.

# Returns
- Tuple `(eqs, len, spring_force)` with updated equation vector
  and the segment length and spring force variables for use by other components.
"""
function segment_eqs!(s, eqs, points, segments,
                      pulleys, tethers, bodies,
                      params; pos, vel, wind_vec_gnd,
                      spring_force_vec, drag_force, l0,
                      pulley_len, tether_len)
    wind_factor = param_computed!(params.reg, :wind_factor, WindFactorReader())
    @variables begin
        # Spring-damper model
        segment_vec(t)[1:3, eachindex(segments)]
        unit_vec(t)[1:3, eachindex(segments)]
        len(t)[eachindex(segments)]
        rel_vel(t)[1:3, eachindex(segments)]
        spring_vel(t)[eachindex(segments)]
        spring_force(t)[eachindex(segments)]
        stiffness(t)[eachindex(segments)]
        damping(t)[eachindex(segments)]
        # Aerodynamic drag model
        segment_height(t)[eachindex(segments)]
        segment_vel(t)[1:3, eachindex(segments)]
        segment_rho(t)[eachindex(segments)]
        wind_vel(t)[1:3, eachindex(segments)]
        va(t)[1:3, eachindex(segments)]
        area(t)[eachindex(segments)]
        app_perp_vel(t)[1:3, eachindex(segments)]
    end

    for segment in segments
        p1, p2 = segment.point_idxs[1], segment.point_idxs[2]

        # WING-WING segments: rigid wings skip spring+drag, particle wings skip drag.
        p1_obj = points[p1]
        p2_obj = points[p2]
        is_wing_structural_segment = (p1_obj.type == WING && p2_obj.type == WING)

        # Check if this is a RIGID_DYNAMICS wing structural segment
        is_rigid_dynamics_wing_segment = false
        if is_wing_structural_segment
            # Both points should belong to the same wing
            wing = bodies[p1_obj.wing_idx]
            is_rigid_dynamics_wing_segment = (wing.dynamics_type == RIGID_DYNAMICS)
        end

        in_pulley = 0
        for pulley in pulleys
            if segment.idx == pulley.segment_idxs[1]
                eqs = [eqs; l0[segment.idx] ~ pulley_len[pulley.idx]]
                in_pulley += 1
            end
            if segment.idx == pulley.segment_idxs[2]
                eqs = [
                    eqs
                    l0[segment.idx] ~
                        params.pulleys[pulley.idx].sum_len - pulley_len[pulley.idx]
                ]
                in_pulley += 1
            end
        end
        (in_pulley > 1) && error(
            "Bridle segment $(segment.idx) is in $in_pulley pulleys; " *
            "should be in 0 or 1.",
        )

        if in_pulley == 0
            in_tether = 0
            tether_idx = 0
            for tether in tethers
                if segment.idx in tether.segment_idxs
                    tether_idx = tether.idx
                    in_tether += 1
                end
            end
            !(in_tether in [0, 1]) && error(
                "Segment $(segment.idx) is in " *
                "$in_tether tethers; should be 0 or 1.",
            )

            if in_tether == 1
                # l0 = tether_len / n_segments (winched and winchless alike).
                n_segs = length(
                    tethers[tether_idx].segment_idxs)
                eqs = [
                    eqs
                    l0[segment.idx] ~
                        tether_len[tether_idx] / n_segs
                ]
            else
                eqs = [eqs;
                    l0[segment.idx] ~
                        params.segments[segment.idx].l0]
            end
        end

        # Geometric quantities (always needed)
        eqs = [
            eqs
            segment_vec[:, segment.idx] ~ pos[:, p2] - pos[:, p1]
            len[segment.idx] ~ smooth_norm(segment_vec[:, segment.idx])
            unit_vec[:, segment.idx] ~ segment_vec[:, segment.idx] / len[segment.idx]
            rel_vel[:, segment.idx] ~ vel[:, p1] - vel[:, p2]
            spring_vel[segment.idx] ~ rel_vel[:, segment.idx] ⋅ unit_vec[:, segment.idx]
        ]

        # Spring force: zero for RIGID_DYNAMICS wing segments (rigid body), computed otherwise
        if is_rigid_dynamics_wing_segment
            eqs = [
                eqs
                damping[segment.idx] ~ 0.0
                stiffness[segment.idx] ~ 0.0
                spring_force[segment.idx] ~ 0.0
                spring_force_vec[:, segment.idx] ~ zeros(3)
            ]
        else
            eqs = [
                eqs
                damping[segment.idx] ~
                    params.segments[segment.idx].unit_damping / len[segment.idx]
                stiffness[segment.idx] ~ ifelse(
                    len[segment.idx] > l0[segment.idx],
                    params.segments[segment.idx].unit_stiffness / len[segment.idx],
                    params.segments[segment.idx].compression_frac *
                    params.segments[segment.idx].unit_stiffness / len[segment.idx],
                )
                spring_force[segment.idx] ~ (
                    stiffness[segment.idx] * (len[segment.idx] - l0[segment.idx]) -
                    damping[segment.idx] * spring_vel[segment.idx]
                )
                spring_force_vec[:, segment.idx] ~
                    spring_force[segment.idx] * unit_vec[:, segment.idx]
            ]
        end

        # Aerodynamic properties for all segments
        segment_pos_z = 0.5 * (pos[3, p1] + pos[3, p2])
        eqs = [
            eqs
            segment_height[segment.idx] ~ max(0.0, segment_pos_z)
            segment_vel[:, segment.idx] ~ 0.5 * (vel[:, p1] + vel[:, p2])
            segment_rho[segment.idx] ~ calc_rho(s.am, segment_height[segment.idx])
            wind_vel[:, segment.idx] ~
                wind_factor(segment_pos_z) * wind_vec_gnd
            va[:, segment.idx] ~
                wind_vel[:, segment.idx] - segment_vel[:, segment.idx]
            area[segment.idx] ~
                len[segment.idx] * params.segments[segment.idx].diameter
            app_perp_vel[:, segment.idx] ~
                va[:, segment.idx] -
                (va[:, segment.idx] ⋅ unit_vec[:, segment.idx]) *
                unit_vec[:, segment.idx]
        ]

        # Drag force: zero for wing structural segments (forces from VSM), otherwise computed
        if is_wing_structural_segment
            eqs = [eqs; drag_force[:, segment.idx] ~ zeros(3)]
        else
            eqs = [
                eqs
                drag_force[:, segment.idx] ~
                    (
                        0.5 * segment_rho[segment.idx] * params.set.cd_tether *
                        smooth_norm(va[:, segment.idx]) * area[segment.idx]
                    ) * app_perp_vel[:, segment.idx]
            ]
        end
    end

    return eqs, len, spring_force
end
