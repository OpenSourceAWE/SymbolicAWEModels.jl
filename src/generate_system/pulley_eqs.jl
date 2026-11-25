# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Pulley dynamics equation generation

"""
    pulley_eqs!(eqs, defaults, pulleys, params, initial;
                spring_force, pulley_len, pulley_vel, l0)

Generate equations for pulley dynamics (rope distribution over pulleys).

The rope mass and the split dynamics come from the shared
[`pulley_rope_mass`](@ref) and [`pulley_split_eqs`](@ref); `pulley_force` (the
tension imbalance driving the split) and `pulley_acc` (that imbalance over the
rope mass) are emitted as diagnostics on top.

# Arguments
- `eqs`, `defaults`: Accumulating vectors for the MTK system.
- `pulleys`: The pulley rows of the definition.
- `spring_force`: Symbolic segment spring force variable.
- `pulley_len`, `pulley_vel`: Symbolic pulley state variables.
- `l0`: Symbolic segment rest lengths (segment 1's `l0` is bound to the pulley length).

# Returns
- Tuple `(eqs, defaults)` with updated equation vectors.
"""
function pulley_eqs!(eqs, defaults, pulleys, params, initial;
                     spring_force, pulley_len, pulley_vel, l0)
    @variables begin
        pulley_force(t)[eachindex(pulleys)]
        pulley_acc(t)[eachindex(pulleys)]
    end
    for (i, pulley) in enumerate(pulleys)
        pulley_type(pulley) == DYNAMIC || error("Wrong pulley type")
        mass = pulley_rope_mass(params, i, pulley.segments[1])
        line_tension = 0.5 * (spring_force[pulley.segments[1]] +
                              spring_force[pulley.segments[2]])
        eqs = [
            eqs
            pulley_force[i] ~
                spring_force[pulley.segments[1]] -
                spring_force[pulley.segments[2]]
            pulley_acc[i] ~ pulley_force[i] / mass
            pulley_split_eqs(pulley_len[i], pulley_vel[i],
                             pulley_force[i], line_tension, mass,
                             params.pulleys[i])
        ]
        defaults = [
            defaults
            bind_initial!(initial.pulleys[i].len, pulley_len[i])
            bind_initial!(initial.pulleys[i].vel, pulley_vel[i])
            # l0[seg1] == pulley_len; bind for the tearing that keeps l0.
            bind_initial!(initial.pulleys[i].len, l0[pulley.segments[1]])
        ]
    end
    return eqs, defaults
end
