# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

# VSM aerodynamics equation generation
#
# Aero mode is a build-time decision (part of model SHA hash).
# Each mode generates only the equations it needs:
#   AERO_NONE:       zeros (no VSM calls)
#   AERO_DIRECT:     reads stored forces via registered functions
#   AERO_LINEARIZED: wind-axis coefficient linearization
#                    with q∞·A·(CL·lift - CD·drag + CS·side)
#
# REFINE wings use per-point forces via get_point_aero_force.

"""
    vsm_eqs!(s, eqs, guesses, psys; kwargs...)

Generate aerodynamic equations for all wings.

Aero mode is resolved at build time — each wing's `aero_mode`
determines which equations are generated:
- `AERO_LINEARIZED`: wind-axis coefficient equations with
  Jacobian-based linearization around the operating point
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

    # Predeclare symbolic arrays
    aero_input = nothing
    aero_input_delta = nothing
    aero_input_prev = nothing
    aero_jac_sym = nothing
    aero_coeffs_prev = nothing
    q_inf = nothing

    has_linearized = any(
        w isa VSMWing &&
        w.wing_type == QUATERNION &&
        w.aero_mode == AERO_LINEARIZED
        for w in wings)

    # Declare symbolic variables for AERO_LINEARIZED
    if has_linearized
        first_lin_wing = first(
            w for w in wings if w isa VSMWing &&
            w.wing_type == QUATERNION &&
            w.aero_mode == AERO_LINEARIZED)
        ny = length(first_lin_wing.aero_y)
        nx_values = [
            length(w.aero_x)
            for w in wings
            if w isa VSMWing &&
               w.wing_type == QUATERNION &&
               w.aero_mode == AERO_LINEARIZED]
        nx = maximum(nx_values)

        @variables begin
            aero_input(t)[1:ny, eachindex(wings)]
            aero_input_delta(t)[
                1:ny, eachindex(wings)]
            aero_input_prev(t)[
                1:ny, eachindex(wings)]
            aero_jac_sym(t)[
                1:nx, 1:ny, eachindex(wings)]
            aero_coeffs_prev(t)[
                1:nx, eachindex(wings)]
            q_inf(t)[eachindex(wings)]
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
            # ========== QUATERNION + AERO_NONE =====
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
            # ========== QUATERNION + AERO_DIRECT ===
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
            # ========== QUATERNION + AERO_LINEARIZED
            # Wind-axis coefficient linearization
            wing isa VSMWing || error(
                "AERO_LINEARIZED wing $(wing.idx)" *
                " is not a VSMWing")

            area = wing.vsm_aero.projected_area
            c_ref = wing.vsm_aero.c_ref
            n_groups = length(wing.group_idxs)
            ny_w = length(wing.aero_y)
            nx_w = length(wing.aero_x)

            # ── Load stored operating point ──────
            local prev_input_eqs = Equation[]
            for iy in 1:ny_w
                push!(prev_input_eqs,
                    aero_input_prev[iy, wing.idx] ~
                    get_aero_y(
                        psys, wing.idx, iy))
            end

            local prev_coeff_eqs = Equation[]
            for ix in 1:nx_w
                push!(prev_coeff_eqs,
                    aero_coeffs_prev[ix, wing.idx] ~
                    get_aero_x(
                        psys, wing.idx, ix))
            end

            local jac_eqs = Equation[]
            for ix in 1:nx_w
                for iy in 1:ny_w
                    push!(jac_eqs,
                        aero_jac_sym[
                            ix, iy, wing.idx] ~
                        get_aero_jac(
                            psys, wing.idx, ix, iy))
                end
            end

            # ── Current input state (symbolic) ───
            w = wing.idx
            vx = va_wing_b[1, w]
            vy = va_wing_b[2, w]
            vz = va_wing_b[3, w]
            va_sq = vx^2 + vy^2 + vz^2
            va_mag = sqrt(va_sq + 1e-24)
            alpha_sym = atan(vz, vx)
            beta_sym = asin(vy / va_mag)

            # Per-group twist inputs
            twist_inputs = [
                twist_angle[groups[gidx].idx]
                for gidx in wing.group_idxs]

            eqs = [
                eqs
                # Dynamic pressure
                q_inf[w] ~
                    0.5 *
                    calc_rho(s.am,
                        wing_pos[3, w]) * va_sq

                # Load stored state from struct
                prev_input_eqs
                prev_coeff_eqs
                jac_eqs

                # Assemble current input state
                aero_input[:, wing.idx] ~ [
                    alpha_sym
                    beta_sym
                    ω_b[1, wing.idx]
                    ω_b[2, wing.idx]
                    ω_b[3, wing.idx]
                    twist_inputs
                ]

                # Delta = current - operating point
                aero_input_delta[:, wing.idx] ~
                    aero_input[:, wing.idx] -
                    aero_input_prev[:, wing.idx]
            ]

            # ── Coefficient reconstruction ───────
            # coeff(ix) = x0[ix] + Σ J[ix,iy]*Δ[iy]
            delta = aero_input_delta[:, w]
            J = aero_jac_sym[:, :, w]
            x0 = aero_coeffs_prev[:, w]

            coeff(ix) = x0[ix] + sum(
                J[ix, iy] * delta[iy]
                for iy in 1:ny_w)

            CL = coeff(1)
            CD = coeff(2)
            CS = coeff(3)
            qA = q_inf[w] * area

            # ── Force reconstruction ─────────────
            # drag_dir = va_b / |va_b|  (scalar)
            # lift_dir = normalize(cross(drag_dir,
            #            [0,1,0]))
            #          = [-vz, 0, vx] / sqrt(vx²+vz²)
            inv_va = 1 / va_mag
            xz_mag = sqrt(vx^2 + vz^2 + 1e-24)
            inv_xz = 1 / xz_mag

            drag = [vx, vy, vz] .* inv_va
            lift = [-vz, 0, vx] .* inv_xz
            side = [0.0, 1.0, 0.0]

            drag_frac = get_drag_frac(psys, w)
            force_eq = [
                qA * (CL * lift[i] -
                      CD * drag_frac * drag[i] +
                      CS * side[i])
                for i in 1:3]

            moment_eq = [
                qA * c_ref * coeff(3 + i)
                for i in 1:3]

            # ── Group moments ────────────────────
            group_moment_eqs = Equation[]
            for (gi, gidx) in
                    enumerate(wing.group_idxs)
                group = groups[gidx]
                isempty(
                    group.unrefined_section_idxs
                ) && continue
                push!(group_moment_eqs,
                    group_aero_moment[group.idx] ~
                        qA * c_ref * coeff(6 + gi))
            end

            eqs = [
                eqs
                group_moment_eqs

                aero_force_b[:, wing.idx] ~ force_eq
                aero_moment_b[:, wing.idx] ~
                    moment_eq
            ]

            if s.set.quasi_static
                local wing_guesses = []
                for iy in 1:ny_w
                    push!(wing_guesses,
                        aero_input[iy, wing.idx] =>
                        get_aero_y(
                            psys, wing.idx, iy))
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
