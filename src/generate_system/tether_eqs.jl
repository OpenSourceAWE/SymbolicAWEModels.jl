# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Tether aggregation equation generation

"""
    tether_eqs!(eqs, tethers; len, spring_force)

Generate equations for tether stretched length and average spring force.

# Arguments
- `eqs`: Accumulating equation vector.
- `tethers`: Collection of Tether objects.
- `len`: Symbolic segment length variable.
- `spring_force`: Symbolic segment spring force variable.

# Returns
- Updated `eqs` vector with tether equations.
"""
function tether_eqs!(eqs, tethers; len, spring_force)
    @variables begin
        stretched_len(t)[eachindex(tethers)]
        tether_spring_force(t)[eachindex(tethers)]
    end
    for tether in tethers
        slen = zero(Num)
        tforce = zero(Num)
        for segment_idx in tether.segment_idxs
            slen += len[segment_idx]
            tforce += spring_force[segment_idx]
        end
        tforce /= length(tether.segment_idxs)
        eqs = [
            eqs
            stretched_len[tether.idx] ~ slen
            tether_spring_force[tether.idx] ~ tforce
        ]
    end
    return eqs
end
