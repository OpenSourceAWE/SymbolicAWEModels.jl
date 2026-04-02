# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

# Winch motor dynamics equation generation

"""
    winch_eqs!(eqs, defaults, winches, tethers, points, psys, _pset;
               point_force, set_values, tether_len, tether_vel)

Generate equations for winch motor dynamics and tether reeling.

# Arguments
- `eqs`, `defaults`: Accumulating vectors for the MTK system.
- `winches`: Collection of Winch objects.
- `tethers`: Collection of Tether objects.
- `points`: Collection of Point objects.
- `psys`, `pset`: Symbolic parameters representing system and settings.
- `point_force`: Symbolic point force variable.
- `set_values`: Symbolic winch torque setpoint variable.
- `tether_len`, `tether_vel`: Symbolic tether state variables.

# Returns
- Tuple `(eqs, defaults)` with updated equation vectors.
"""
function winch_eqs!(eqs, defaults, winches, tethers, points, psys, _pset;
                    point_force, set_values, tether_len, tether_vel)
    @variables begin
        tether_acc(t)[eachindex(winches)]
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

    for winch in winches
        F = zeros(Num, 3)
        for tether_idx in winch.tether_idxs
            winch_point_idx = tethers[tether_idx].winch_point_idx
            (winch_point_idx > length(points)) &&
                error("Point number $winch_point_idx does not exist.")
            F .+= point_force[:, winch_point_idx]
        end

        gear_ratio = get_winch_gear_ratio(psys, winch.idx)
        drum_radius = get_winch_drum_radius(psys, winch.idx)
        f_coulomb = get_winch_f_coulomb(psys, winch.idx)
        c_vf = get_winch_c_vf(psys, winch.idx)
        inertia_total = get_winch_inertia_total(psys, winch.idx)
        friction_eps = get_winch_friction_epsilon(psys, winch.idx)

        # Smooth sign function to avoid discontinuities
        # at zero velocity. eps controls transition width.
        smooth_sign(x, eps) = x / sqrt(x * x + eps * eps)

        eqs = [
            eqs
            brake[winch.idx] ~ get_brake(psys, winch.idx)
            speed_controlled[winch.idx] ~
                get_speed_controlled(psys, winch.idx)
            D(tether_len[winch.idx]) ~
                ifelse(brake[winch.idx] == true, 0,
                       tether_vel[winch.idx])
            D(tether_vel[winch.idx]) ~
                ifelse(brake[winch.idx] == true, 0,
                    ifelse(speed_controlled[winch.idx] == true,
                           0, tether_acc[winch.idx]))

            # Winch motor, gear, and friction dynamics
            ω_motor[winch.idx] ~
                gear_ratio / drum_radius * tether_vel[winch.idx]
            tau_friction[winch.idx] ~
                smooth_sign(ω_motor[winch.idx], friction_eps) *
                f_coulomb * drum_radius / gear_ratio +
                c_vf * ω_motor[winch.idx] *
                drum_radius^2 / gear_ratio^2
            tau_motor[winch.idx] ~ set_values[winch.idx]
            tau_total[winch.idx] ~
                tau_motor[winch.idx] +
                drum_radius / gear_ratio * winch_force[winch.idx] -
                tau_friction[winch.idx]
            α_motor[winch.idx] ~ tau_total[winch.idx] / inertia_total
            tether_acc[winch.idx] ~
                drum_radius / gear_ratio * α_motor[winch.idx]

            winch_force_vec[:, winch.idx] ~ F
            winch_force[winch.idx] ~ norm(winch_force_vec[:, winch.idx])
        ]
        defaults = [
            defaults
            tether_len[winch.idx] => get_tether_len(psys, winch.idx)
            tether_vel[winch.idx] => get_tether_vel(psys, winch.idx)
        ]
    end
    return eqs, defaults
end
