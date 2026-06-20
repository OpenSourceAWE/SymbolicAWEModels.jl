# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Pulley dynamics equation generation

"""
    pulley_eqs!(eqs, defaults, pulleys, segments, params;
                spring_force, pulley_len, pulley_vel)

Generate equations for pulley dynamics (rope distribution over pulleys).

# Arguments
- `eqs`, `defaults`: Accumulating vectors for the MTK system.
- `pulleys`: Collection of Pulley objects.
- `segments`: Collection of Segment objects (for mass calculation).
- `spring_force`: Symbolic segment spring force variable.
- `pulley_len`, `pulley_vel`: Symbolic pulley state variables.

# Returns
- Tuple `(eqs, defaults)` with updated equation vectors.
"""
function pulley_eqs!(eqs, defaults, pulleys, segments, params, initial;
                     spring_force, pulley_len, pulley_vel)
    @variables begin
        pulley_force(t)[eachindex(pulleys)]
        pulley_acc(t)[eachindex(pulleys)]
    end
    @parameters pulley_damp = 5.0

    for pulley in pulleys
        segment = segments[pulley.segment_idxs[1]]
        mass_per_meter =
            params.segments[segment.idx].density * π * (params.segments[segment.idx].diameter / 2)^2
        mass = params.pulleys[pulley.idx].sum_len * mass_per_meter
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
                bind_initial!(initial.pulleys[pulley.idx].len, pulley_len[pulley.idx])
                bind_initial!(initial.pulleys[pulley.idx].vel, pulley_vel[pulley.idx])
            ]
        else
            error("Wrong pulley type")
        end
    end
    return eqs, defaults
end
