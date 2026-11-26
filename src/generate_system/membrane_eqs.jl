# Copyright (c) 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

"""
    membrane_eqs!(eqs, membranes, params; pos, vel, membrane_force)

Append the force on each corner of each membrane, [`membrane_corner_forces`](@ref) of
its corners' `pos` and `vel`, to `eqs` as `membrane_force[:, corner, idx]`.
"""
function membrane_eqs!(eqs, membranes, params; pos, vel, membrane_force)
    for membrane in membranes
        corners = membrane.point_idxs
        forces = membrane_corner_forces(params, membrane.idx,
            [collect(pos[:, idx]) for idx in corners],
            [collect(vel[:, idx]) for idx in corners])
        for k in 1:3
            eqs = [eqs; collect(membrane_force[:, k, membrane.idx]) .~ forces[:, k]]
        end
    end
    return eqs
end
