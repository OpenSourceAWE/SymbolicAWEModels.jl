# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

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
            ny_w = length(wing.aero_y)
            nx_w = length(wing.aero_x)
            w = wing.idx

            # ── Load stored operating point ──────
            prev_input_eqs = [
                aero_input_prev[iy, w] ~
                    get_aero_y(psys, w, iy)
                for iy in 1:ny_w]

            prev_coeff_eqs = [
                aero_coeffs_prev[ix, w] ~
                    get_aero_x(psys, w, ix)
                for ix in 1:nx_w]

            jac_eqs = [
                aero_jac_sym[ix, iy, w] ~
                    get_aero_jac(psys, w, ix, iy)
                for ix in 1:nx_w for iy in 1:ny_w]

            # ── Current input state (symbolic) ───
            # collect() so smooth_norm's mapreduce scalarises
            va = collect(va_wing_b[:, w])
            drag_dir = collect(va ./ smooth_norm(va))
            alpha_sym = atan(drag_dir[3], drag_dir[1])
            beta_sym = asin(drag_dir[2])

            twist_inputs = [
                twist_angle[groups[gidx].idx]
                for gidx in wing.group_idxs]

            eqs = [
                eqs
                q_inf[w] ~
                    0.5 *
                    calc_rho(s.am, wing_pos[3, w]) *
                    (va ⋅ va)

                prev_input_eqs
                prev_coeff_eqs
                jac_eqs

                aero_input[:, w] ~ [
                    alpha_sym
                    beta_sym
                    ω_b[1, w]
                    ω_b[2, w]
                    ω_b[3, w]
                    twist_inputs
                ]

                aero_input_delta[:, w] ~
                    aero_input[:, w] -
                    aero_input_prev[:, w]
            ]

            # ── Coefficient reconstruction ───────
            # coeff(ix) = x0[ix] + Σ J[ix,iy]*Δ[iy]
            delta = aero_input_delta[:, w]
            J = aero_jac_sym[:, :, w]
            x0 = aero_coeffs_prev[:, w]

            coeff(ix) = x0[ix] + sum(
                J[ix, iy] * delta[iy] for iy in 1:ny_w)

            CL = coeff(1)
            CD = coeff(2)
            CS = coeff(3)
            qA = q_inf[w] * area

            # ── Wind-axis basis (matches VSM) ────
            # drag = va / |va|
            # lift = normalize(drag × span)
            # side = lift × drag   (orthonormal triad)
            crossed = collect(drag_dir × [0.0, 1.0, 0.0])
            lift_dir = collect(
                crossed ./ smooth_norm(crossed))
            side_dir = collect(lift_dir × drag_dir)

            drag_frac = get_drag_frac(psys, w)
            force_eq = collect(qA * (
                CL * lift_dir +
                CD * drag_frac * drag_dir +
                CS * side_dir))

            moment_eq = [
                qA * c_ref * coeff(3 + i) for i in 1:3]

            group_moment_eqs = [
                group_aero_moment[groups[gidx].idx] ~
                    qA * c_ref * coeff(6 + gi)
                for (gi, gidx) in
                    enumerate(wing.group_idxs)
                if !isempty(
                    groups[gidx].unrefined_section_idxs)]

            eqs = [
                eqs
                group_moment_eqs

                aero_force_b[:, w] ~ force_eq
                aero_moment_b[:, w] ~ moment_eq
            ]

            if s.set.quasi_static
                wing_guesses = [
                    aero_input[iy, w] =>
                        get_aero_y(psys, w, iy)
                    for iy in 1:ny_w]
                guesses = [
                    guesses
                    wing_guesses
                ]
            end
        end
    end
    return eqs, guesses
end
