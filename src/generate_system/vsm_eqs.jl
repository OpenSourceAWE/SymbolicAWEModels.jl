# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

# VSM aerodynamics equation generation
#
# Aero mode is a build-time decision (part of model SHA hash).
# Each mode generates only the equations it needs:
#   AERO_NONE:       zeros (no VSM calls)
#   AERO_DIRECT:     reads stored forces via registered functions
#   AERO_LINEARIZED: full symbolic linearization (q∞·A·(x₀+J·Δ))
#
# REFINE wings use per-point forces via get_point_aero_force.

"""
    vsm_eqs!(s, eqs, guesses, psys; kwargs...)

Generate aerodynamic equations for all wings.

Aero mode is resolved at build time — each wing's `aero_mode`
determines which equations are generated:
- `AERO_LINEARIZED`: symbolic linearization equations
  (q∞·A·(x₀ + J·Δ)) that the solver can differentiate through
- `AERO_DIRECT`: registered functions returning stored forces
- `AERO_NONE`: zeros

For REFINE wings, per-point forces come from
`get_point_aero_force` (which also respects `AERO_NONE`).
"""
function vsm_eqs!(
    s, eqs, guesses, psys;
    aero_force_b, aero_moment_b, group_aero_moment,
    twist_angle, va_wing_b, wing_pos, ω_b,
    aero_force_point_b=nothing
)
    (; groups, wings, points) = s.sys_struct
    length(wings) == 0 && return eqs, guesses

    # Predeclare symbolic arrays so static analysis sees these bindings
    # even though they are only populated when a linearized wing exists.
    vsm_input_state = nothing
    vsm_input_state_delta = nothing
    vsm_input_state_prev = nothing
    force_jacobian = nothing
    vsm_output_force_prev = nothing
    q_inf = nothing
    no_scale_aero_force_b = nothing

    has_linearized = any(
        w isa VSMWing &&
        w.wing_type == QUATERNION &&
        w.aero_mode == AERO_LINEARIZED
        for w in wings)

    # Declare symbolic variables only when needed
    # for AERO_LINEARIZED QUATERNION wings
    if has_linearized
        first_lin_wing = first(w for w in wings if w isa VSMWing &&
            w.wing_type == QUATERNION && w.aero_mode == AERO_LINEARIZED)
        n_unrefined = first_lin_wing.vsm_wing.n_unrefined_sections
        ny_quaternion = 3 + n_unrefined + 3
        nx_values = [
            3 + 3 + w.vsm_wing.n_unrefined_sections
            for w in wings
            if w isa VSMWing && w.wing_type == QUATERNION &&
               w.aero_mode == AERO_LINEARIZED]
        nx_max = maximum(nx_values)

        @variables begin
            vsm_input_state(t)[
                1:ny_quaternion, eachindex(wings)]
            vsm_input_state_delta(t)[
                1:ny_quaternion, eachindex(wings)]
            vsm_input_state_prev(t)[
                1:ny_quaternion, eachindex(wings)]
            force_jacobian(t)[
                1:nx_max, 1:ny_quaternion,
                eachindex(wings)]
            vsm_output_force_prev(t)[
                1:nx_max, eachindex(wings)]
            q_inf(t)[eachindex(wings)]
            no_scale_aero_force_b(t)[
                1:3, eachindex(wings)]
        end
    end

    for wing in wings
        if wing isa VSMWing && wing.wing_type == REFINE
            # ========== REFINE WING ==========
            afpb = aero_force_point_b::AbstractArray
            wing_points = [
                p for p in points
                if p.type == WING &&
                    p.wing_idx == wing.idx
            ]

            if wing.aero_mode == AERO_NONE
                for point in wing_points
                    eqs = [
                        eqs
                        afpb[:, point.idx] ~ zeros(3)
                    ]
                end
            else
                # AERO_DIRECT: per-point forces
                for point in wing_points
                    eqs = [
                        eqs
                        afpb[:, point.idx] ~ [
                            get_point_aero_force(
                                psys, point.idx, i)
                            for i in 1:3
                        ]
                    ]
                end
            end

            eqs = [
                eqs
                aero_force_b[:, wing.idx] ~
                    sum([afpb[:, p.idx]
                         for p in wing_points])
                aero_moment_b[:, wing.idx] ~ zeros(3)
            ]

        elseif wing.aero_mode == AERO_NONE
            # ========== QUATERNION + AERO_NONE ==========
            eqs = [
                eqs
                aero_force_b[:, wing.idx] ~ zeros(3)
                aero_moment_b[:, wing.idx] ~ zeros(3)
            ]
            for gidx in wing.group_idxs
                group = groups[gidx]
                isempty(
                    group.unrefined_section_idxs
                ) && continue
                eqs = [
                    eqs
                    group_aero_moment[group.idx] ~ 0
                ]
            end

        elseif wing.aero_mode == AERO_DIRECT
            # ========== QUATERNION + AERO_DIRECT ==========
            eqs = [
                eqs
                aero_force_b[:, wing.idx] ~ [
                    get_aero_force_override(
                        psys, wing.idx, c)
                    for c in 1:3]
                aero_moment_b[:, wing.idx] ~ [
                    get_aero_moment_override(
                        psys, wing.idx, c)
                    for c in 1:3]
            ]
            for gidx in wing.group_idxs
                group = groups[gidx]
                isempty(
                    group.unrefined_section_idxs
                ) && continue
                eqs = [
                    eqs
                    group_aero_moment[group.idx] ~
                        get_group_moment_override(
                            psys, wing.idx,
                            Int64(gidx))
                ]
            end

        else
            # ========== QUATERNION + AERO_LINEARIZED =====
            # Full symbolic linearization equations
            wing isa VSMWing || error("AERO_LINEARIZED wing $(wing.idx) is not a VSMWing")

            area = wing.vsm_aero.projected_area
            n_un = wing.vsm_wing.n_unrefined_sections
            ny_quat = 3 + n_un + 3
            nx_quat = 3 + 3 + n_un

            force_b = no_scale_aero_force_b[:, wing.idx]
            wind_dir_b = smooth_normalize(
                va_wing_b[:, wing.idx])
            drag_force_b =
                (force_b ⋅ wind_dir_b) * wind_dir_b

            # Build twist mapping
            unrefined_to_group_twist =
                Vector{Any}(undef, n_un)
            for gidx in wing.group_idxs
                group = groups[gidx]
                for ui in group.unrefined_section_idxs
                    unrefined_to_group_twist[ui] =
                        twist_angle[group.idx]
                end
            end

            local prev_state_eqs = Equation[]
            for iy in 1:ny_quat
                push!(prev_state_eqs,
                    vsm_input_state_prev[iy, wing.idx] ~
                    get_vsm_y(psys, wing.idx, iy))
            end

            local prev_force_eqs = Equation[]
            for ix in 1:nx_quat
                push!(prev_force_eqs,
                    vsm_output_force_prev[ix, wing.idx] ~
                    get_vsm_x(psys, wing.idx, ix))
            end

            local force_jacobian_eqs = Equation[]
            for ix in 1:nx_quat
                for iy in 1:ny_quat
                    push!(force_jacobian_eqs,
                        force_jacobian[ix, iy, wing.idx] ~
                        get_vsm_jac(
                            psys, wing.idx, ix, iy))
                end
            end

            eqs = [
                eqs
                # Dynamic pressure
                q_inf[wing.idx] ~
                    0.5 *
                    calc_rho(s.am,
                        wing_pos[3, wing.idx]) *
                    smooth_norm(va_wing_b[:, wing.idx])^2

                # Load linearization data from struct
                prev_state_eqs
                prev_force_eqs
                force_jacobian_eqs

                # Current input state (symbolic)
                vsm_input_state[:, wing.idx] ~ [
                    va_wing_b[:, wing.idx]
                    unrefined_to_group_twist
                    ω_b[:, wing.idx]
                ]

                # Δstate = state - state₀
                vsm_input_state_delta[:, wing.idx] ~
                    vsm_input_state[:, wing.idx] -
                    vsm_input_state_prev[:, wing.idx]
            ]

            # Symbolic linearized expressions
            delta = vsm_input_state_delta[:, wing.idx]
            J = force_jacobian[:, :, wing.idx]
            x0 = vsm_output_force_prev[:, wing.idx]
            qA = q_inf[wing.idx] * area

            # Linearized force (1:3), moment (4:6)
            lin_force = qA * (x0[1:3] +
                J[1:3, :] * delta)
            lin_moment = qA * (x0[4:6] +
                J[4:6, :] * delta)

            # Linearized group moments
            group_moment_eqs = Equation[]
            for gidx in wing.group_idxs
                group = groups[gidx]
                isempty(
                    group.unrefined_section_idxs
                ) && continue
                moment_terms = []
                for ui in group.unrefined_section_idxs
                    vix = 6 + ui
                    push!(moment_terms,
                        x0[vix] +
                        sum([J[vix, iy] * delta[iy]
                            for iy in 1:ny_quat
                        ]))
                end
                push!(group_moment_eqs,
                    group_aero_moment[group.idx] ~
                        sum(moment_terms))
            end

            # Drag correction on linearized force
            lin_force_corrected =
                force_b +
                drag_force_b *
                    (get_drag_frac(
                        psys, wing.idx) - 1)

            eqs = [
                eqs
                # Intermediate linearized force
                force_b ~ lin_force

                # Group moments
                group_moment_eqs

                # Final force (with drag correction)
                aero_force_b[:, wing.idx] ~
                    lin_force_corrected

                # Final moment
                aero_moment_b[:, wing.idx] ~
                    lin_moment
            ]

            if s.set.quasi_static
                local wing_guesses = []
                for iy in 1:ny_quat
                    push!(wing_guesses,
                        vsm_input_state[iy, wing.idx] =>
                        get_vsm_y(psys, wing.idx, iy))
                end
                guesses = [
                    guesses
                    wing_guesses
                ]
            end
        end
    end
    return eqs, guesses
end
