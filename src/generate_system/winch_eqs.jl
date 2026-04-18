# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

# Winch motor dynamics and per-tether length equation generation

"""
    winch_eqs!(eqs, defaults, winches, tethers, points,
               psys;
               point_force, set_values,
               tether_len, winch_vel)

Generate equations for winch motor dynamics and per-tether
length state.

Each tether gets a differential equation for `tether_len`:
- **With winch:** `D(tether_len) = winch_vel` (shared
  velocity from the winch motor).
- **Without winch:** `D(tether_len) = 0` (constant).

Each winch gets a differential equation for `winch_vel`:
- `D(winch_vel) = winch_acc` (from motor dynamics),
  gated by `brake` and `speed_controlled`.

# Returns
- Tuple `(eqs, defaults)` with updated equation vectors.
"""
function winch_eqs!(eqs, defaults, winches, tethers,
                    points, psys;
                    point_force, set_values,
                    tether_len, winch_vel)
    @variables begin
        winch_acc(t)[eachindex(winches)]
        winch_force(t)[eachindex(winches)]
        winch_force_vec(t)[1:3, eachindex(winches)]
        brake(t)[eachindex(winches)]
        speed_controlled(t)[eachindex(winches)]
        # Winch motor and friction dynamics
        ω_motor(t)[eachindex(winches)]
        tau_friction(t)[eachindex(winches)]
        tau_motor(t)[eachindex(winches)]
        tau_total(t)[eachindex(winches)]
        α_motor(t)[eachindex(winches)]
    end

    # Build tether → winch lookup
    tether_winch = Dict{Int, Int}()
    for winch in winches
        for ti in winch.tether_idxs
            if haskey(tether_winch, ti)
                error("Tether $ti is connected " *
                    "to winch $(tether_winch[ti]) " *
                    "and winch $(winch.idx). Each " *
                    "tether can have at most one " *
                    "winch.")
            end
            tether_winch[ti] = winch.idx
        end
    end

    # --- Per-tether length equations ---
    for tether in tethers
        if haskey(tether_winch, tether.idx)
            wi = tether_winch[tether.idx]
            eqs = [
                eqs
                D(tether_len[tether.idx]) ~
                    ifelse(brake[wi] == true, 0,
                           winch_vel[wi])
            ]
        else
            # Winchless tether: constant length
            eqs = [
                eqs
                D(tether_len[tether.idx]) ~ 0
            ]
        end
        defaults = [
            defaults
            tether_len[tether.idx] =>
                get_tether_len(psys, tether.idx)
        ]
    end

    # --- Per-winch velocity and motor dynamics ---
    for winch in winches
        isempty(winch.tether_idxs) &&
            error("Winch $(winch.name): no connected " *
                  "tethers; at least one is required.")
        winch_point_idx = winch.winch_point_idx
        (winch_point_idx > length(points)) &&
            error("Winch $(winch.name): point " *
                  "$winch_point_idx does not exist.")
        F = point_force[:, winch_point_idx]

        gear_ratio = get_winch_gear_ratio(
            psys, winch.idx)
        drum_radius = get_winch_drum_radius(
            psys, winch.idx)
        f_coulomb = get_winch_f_coulomb(
            psys, winch.idx)
        c_vf = get_winch_c_vf(psys, winch.idx)
        inertia_total = get_winch_inertia_total(
            psys, winch.idx)
        friction_eps = get_winch_friction_epsilon(
            psys, winch.idx)

        # Smooth sign function to avoid discontinuities
        # at zero velocity. eps controls transition width.
        smooth_sign(x, eps) =
            x / sqrt(x * x + eps * eps)

        eqs = [
            eqs
            D(winch_vel[winch.idx]) ~
                ifelse(brake[winch.idx] == true, 0,
                    ifelse(
                        speed_controlled[winch.idx]
                            == true,
                        0, winch_acc[winch.idx]))
            brake[winch.idx] ~
                get_brake(psys, winch.idx)
            speed_controlled[winch.idx] ~
                get_speed_controlled(psys, winch.idx)

            # Winch motor, gear, and friction dynamics
            ω_motor[winch.idx] ~
                gear_ratio / drum_radius *
                winch_vel[winch.idx]
            tau_friction[winch.idx] ~
                smooth_sign(
                    ω_motor[winch.idx], friction_eps) *
                f_coulomb * drum_radius / gear_ratio +
                c_vf * ω_motor[winch.idx] *
                drum_radius^2 / gear_ratio^2
            tau_motor[winch.idx] ~ set_values[winch.idx]
            tau_total[winch.idx] ~
                tau_motor[winch.idx] +
                drum_radius / gear_ratio *
                winch_force[winch.idx] -
                tau_friction[winch.idx]
            α_motor[winch.idx] ~
                tau_total[winch.idx] / inertia_total
            winch_acc[winch.idx] ~
                drum_radius / gear_ratio *
                α_motor[winch.idx]

            winch_force_vec[:, winch.idx] ~ F
            winch_force[winch.idx] ~
                smooth_norm(winch_force_vec[:, winch.idx])
        ]
        defaults = [
            defaults
            winch_vel[winch.idx] =>
                get_winch_vel(psys, winch.idx)
        ]
    end
    return eqs, defaults
end
