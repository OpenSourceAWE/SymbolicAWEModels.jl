# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Pulley dynamics equation generation

"""
    pulley_eqs!(eqs, defaults, guesses, pulleys, segments, psys;
                spring_force, pulley_len, pulley_vel)

Generate equations for pulley dynamics (rope distribution over pulleys).

# Arguments
- `eqs`, `defaults`, `guesses`: Accumulating vectors for the MTK system.
- `pulleys`: Collection of Pulley objects.
- `segments`: Collection of Segment objects (for mass calculation).
- `psys`: Symbolic parameter representing the system structure.
- `spring_force`: Symbolic segment spring force variable.
- `pulley_len`, `pulley_vel`: Symbolic pulley state variables.

# Returns
- Tuple `(eqs, defaults, guesses)` with updated equation vectors.
"""
function pulley_eqs!(eqs, defaults, guesses, pulleys, segments, psys;
                     spring_force, pulley_len, pulley_vel)
    @variables begin
        pulley_force(t)[eachindex(pulleys)]
        pulley_acc(t)[eachindex(pulleys)]
    end
    @parameters pulley_damp = 5.0

    for pulley in pulleys
        segment = segments[pulley.segment_idxs[1]]
        mass_per_meter =
            get_density(psys, segment.idx) * π * (get_diameter(psys, segment.idx) / 2)^2
        mass = get_sum_len(psys, pulley.idx) * mass_per_meter
        eqs = [
            eqs
            pulley_force[pulley.idx] ~
                spring_force[pulley.segment_idxs[1]] -
                spring_force[pulley.segment_idxs[2]]
            pulley_acc[pulley.idx] ~ pulley_force[pulley.idx] / mass
        ]
        if pulley.type == DYNAMIC
            eqs = [
                eqs
                D(pulley_len[pulley.idx]) ~ pulley_vel[pulley.idx]
                D(pulley_vel[pulley.idx]) ~
                    pulley_acc[pulley.idx] - pulley_damp * pulley_vel[pulley.idx]
            ]
            defaults = [
                defaults
                pulley_len[pulley.idx] => get_pulley_len(psys, pulley.idx)
                pulley_vel[pulley.idx] => get_pulley_vel(psys, pulley.idx)
            ]
        elseif pulley.type == QUASI_STATIC
            eqs = [eqs; pulley_vel[pulley.idx] ~ 0; pulley_acc[pulley.idx] ~ 0]
            guesses = [
                guesses
                pulley_len[pulley.idx] => get_l0(psys, pulley.segment_idxs[1])
            ]
        else
            error("Wrong pulley type")
        end
    end
    return eqs, defaults, guesses
end
