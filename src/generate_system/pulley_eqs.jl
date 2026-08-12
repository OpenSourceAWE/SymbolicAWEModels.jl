# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Pulley dynamics equation generation

"""
    pulley_eqs!(eqs, defaults, pulleys, segments, params;
                spring_force, pulley_len, pulley_vel)

Generate equations for pulley dynamics (rope distribution over pulleys).

The rope mass and the split dynamics come from the shared
[`pulley_rope_mass`](@ref) and [`pulley_split_eqs`](@ref); `pulley_force` (the
tension imbalance driving the split) and `pulley_acc` (that imbalance over the
rope mass) are emitted as diagnostics on top.

# Arguments
- `eqs`, `defaults`: Accumulating vectors for the MTK system.
- `pulleys`: Collection of Pulley objects.
- `segments`: Collection of Segment objects.
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
        pulley.type == DYNAMIC || error("Wrong pulley type")
        mass = pulley_rope_mass(params, pulley.idx, pulley.segment_idxs[1])
        line_tension = 0.5 * (spring_force[pulley.segment_idxs[1]] +
                              spring_force[pulley.segment_idxs[2]])
        eqs = [
            eqs
            pulley_force[pulley.idx] ~
                spring_force[pulley.segment_idxs[1]] -
                spring_force[pulley.segment_idxs[2]]
            pulley_acc[pulley.idx] ~ pulley_force[pulley.idx] / mass
            pulley_split_eqs(pulley_len[pulley.idx], pulley_vel[pulley.idx],
                             pulley_force[pulley.idx], line_tension, mass,
                             params.pulleys[pulley.idx])
        ]
        defaults = [
            defaults
            bind_initial!(initial.pulleys[pulley.idx].len, pulley_len[pulley.idx])
            bind_initial!(initial.pulleys[pulley.idx].vel, pulley_vel[pulley.idx])
            # l0[seg1] == pulley_len; bind for the tearing that keeps l0.
            bind_initial!(initial.pulleys[pulley.idx].len, l0[pulley.segment_idxs[1]])
        ]
    end
    return eqs, defaults
end
