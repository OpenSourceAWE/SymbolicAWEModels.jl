# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

# Winch motor dynamics and per-tether length equation generation

"""
    winch_eqs!(eqs, defaults, winches, tethers, points,
               psys, _pset;
               point_force, set_values, tether_len, tether_vel)

Generate equations for winch motor dynamics and per-tether
length/velocity state.

Each tether gets differential equations for `tether_len` and
`tether_vel`:
- **With winch:** driven by `winch_acc` from motor dynamics,
  gated by `brake` and `speed_controlled`.
- **Without winch:** `tether_len` is constant (zero velocity).

# Returns
- Tuple `(eqs, defaults)` with updated equation vectors.
"""
function winch_eqs!(eqs, defaults, winches, tethers,
                    points, psys, _pset;
                    point_force, set_values,
                    tether_len, tether_vel)
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
            tether_winch[ti] = winch.idx
        end
    end

    # --- Per-tether differential equations ---
    for tether in tethers
        if haskey(tether_winch, tether.idx)
            wi = tether_winch[tether.idx]
            # Tether with winch: driven by winch dynamics
            eqs = [
                eqs
                D(tether_len[tether.idx]) ~
                    ifelse(brake[wi] == true, 0,
                           tether_vel[tether.idx])
                D(tether_vel[tether.idx]) ~
                    ifelse(brake[wi] == true, 0,
                        ifelse(
                            speed_controlled[wi] == true,
                            0, winch_acc[wi]))
            ]
        else
            # Winchless tether: constant length
            eqs = [
                eqs
                D(tether_len[tether.idx]) ~ 0
                D(tether_vel[tether.idx]) ~ 0
            ]
        end
        defaults = [
            defaults
            tether_len[tether.idx] =>
                get_tether_len(psys, tether.idx)
            tether_vel[tether.idx] =>
                get_tether_vel(psys, tether.idx)
        ]
    end

    # --- Per-winch motor dynamics ---
    for winch in winches
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

        # Use first tether velocity for motor speed
        first_tether_idx = winch.tether_idxs[1]
        tv = tether_vel[first_tether_idx]

        # Smooth sign function to avoid discontinuities
        # at zero velocity. eps controls transition width.
        smooth_sign(x, eps) =
            x / sqrt(x * x + eps * eps)

        eqs = [
            eqs
            brake[winch.idx] ~
                get_brake(psys, winch.idx)
            speed_controlled[winch.idx] ~
                get_speed_controlled(psys, winch.idx)

            # Winch motor, gear, and friction dynamics
            ω_motor[winch.idx] ~
                gear_ratio / drum_radius * tv
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
                norm(winch_force_vec[:, winch.idx])
        ]
    end
    return eqs, defaults
end
