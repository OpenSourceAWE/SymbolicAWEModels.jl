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
- `l0`: Symbolic segment rest lengths (segment 1's `l0` is bound to the pulley length).

# Returns
- Tuple `(eqs, defaults)` with updated equation vectors.
"""
function pulley_eqs!(eqs, defaults, pulleys, segments, params, initial;
                     spring_force, pulley_len, pulley_vel, l0)
    @variables begin
        pulley_force(t)[eachindex(pulleys)]
        pulley_acc(t)[eachindex(pulleys)]
    end
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
            braked = params.pulleys[pulley.idx].brake > 0.5
            accel = pulley_acc[pulley.idx] -
                    pulley_friction_force(params.pulleys[pulley.idx],
                                          pulley_vel[pulley.idx],
                                          0.5 * (spring_force[pulley.segment_idxs[1]] +
                                                 spring_force[pulley.segment_idxs[2]])) /
                    mass
            eqs = [
                eqs
                D(pulley_len[pulley.idx]) ~
                    ifelse(braked, zero(accel), pulley_vel[pulley.idx])
                D(pulley_vel[pulley.idx]) ~ ifelse(braked, zero(accel), accel)
            ]
            defaults = [
                defaults
                bind_initial!(initial.pulleys[pulley.idx].len, pulley_len[pulley.idx])
                bind_initial!(initial.pulleys[pulley.idx].vel, pulley_vel[pulley.idx])
                # l0[seg1] == pulley_len; bind for the tearing that keeps l0.
                bind_initial!(initial.pulleys[pulley.idx].len, l0[pulley.segment_idxs[1]])
            ]
        else
            error("Wrong pulley type")
        end
    end
    return eqs, defaults
end
