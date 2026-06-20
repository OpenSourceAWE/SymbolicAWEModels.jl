# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Thin caller for winch motor dynamics: dispatches to `winch_component` on the
# winch's `model` (see winch_models/). Mirrors aero_eqs.jl.

"""
    winch_eqs!(eqs, defaults, winches, tethers, segments, points,
               sys_struct, params;
               spring_force_vec, set_values, tether_len,
               winch_vel, winch_acc, winch_force_vec, winch_friction)

Generate equations for winch motor dynamics and per-tether length
state, and return the list of `ODESystem` subsystems to attach to the
parent system.

For each winch:
1. Sum spring force vectors of the segments meeting the winch point
   (sign-aware via segment `point_idxs`) → `winch_force_vec`.
2. Instantiate the winch component via `winch_component(winch.model, …)`.
3. Validate the connector contract with
   [`validate_winch_component`](@ref).
4. Bind connectors to the parent variables (including
   `subsys.len ~ mean(tether_len[tether_idx] for tether_idx in winch.tether_idxs)`)
   and integrate `D(winch_vel) = ifelse(brake > 0.5, 0, winch_acc)`.
   When `winch.speed_controlled` is true, `winch_acc` is forced to 0
   (ignoring `subsys.acc`) so velocity is prescribed via `winch.vel`.

For each tether:
- With winch:    `D(tether_len) = ifelse(brake > 0.5, 0, winch_vel)`.
- Without winch: `D(tether_len) = 0`.
"""
function winch_eqs!(eqs, defaults, winches, tethers, segments, points,
                    sys_struct, params, initial;
                    spring_force_vec, set_values, tether_len,
                    winch_vel, winch_acc, winch_force_vec, winch_force,
                    winch_friction)
    tether_winch = Dict{Int, Int}()
    for winch in winches
        for tether_idx in winch.tether_idxs
            haskey(tether_winch, tether_idx) && error(
                "Tether $tether_idx is connected to winch " *
                "$(tether_winch[tether_idx]) and winch $(winch.idx). " *
                "Each tether can have at most one winch.")
            tether_winch[tether_idx] = winch.idx
        end
    end

    for tether in tethers
        if haskey(tether_winch, tether.idx)
            winch_idx = tether_winch[tether.idx]
            eqs = [eqs
                   D(tether_len[tether.idx]) ~
                       ifelse(params.winches[winch_idx].brake > 0.5,
                              0, winch_vel[winch_idx])]
        else
            eqs = [eqs; D(tether_len[tether.idx]) ~ 0]
        end
        defaults = [defaults
                    bind_initial!(initial.tethers[tether.idx].len,
                                  tether_len[tether.idx])]
    end

    winch_subsystems = Any[]
    for winch in winches
        isempty(winch.tether_idxs) &&
            error("Winch $(winch.name): no connected tethers; " *
                  "at least one is required.")
        winch_point_idx = winch.winch_point_idx
        (winch_point_idx > length(points)) &&
            error("Winch $(winch.name): point $winch_point_idx does not exist.")

        winch_seg_idxs = Set{Int}()
        for tether_idx in winch.tether_idxs
            union!(winch_seg_idxs, tethers[tether_idx].segment_idxs)
        end
        force_vec = zeros(Num, 3)
        for segment in segments
            segment.idx in winch_seg_idxs || continue
            if segment.point_idxs[1] == winch_point_idx
                force_vec .+= spring_force_vec[:, segment.idx]
            elseif segment.point_idxs[2] == winch_point_idx
                force_vec .-= spring_force_vec[:, segment.idx]
            end
        end

        subsys = winch_component(winch.model, sys_struct, winch.idx;
                                 name=Symbol("winch_$(winch.idx)"), params)
        validate_winch_component(subsys, winch)
        push!(winch_subsystems, subsys)

        brake_p = params.winches[winch.idx].brake
        eqs = [eqs
               winch_force_vec[:, winch.idx] ~ force_vec
               winch_force[winch.idx] ~
                   smooth_norm(winch_force_vec[:, winch.idx])
               subsys.vel       ~ winch_vel[winch.idx]
               subsys.len       ~
                   sum(tether_len[tether_idx] for tether_idx in winch.tether_idxs) /
                   length(winch.tether_idxs)
               subsys.force     ~ winch_force[winch.idx]
               subsys.set_value ~ set_values[winch.idx]
               subsys.brake     ~ brake_p
               winch_acc[winch.idx]      ~
                   ifelse(params.winches[winch.idx].speed_controlled == true,
                          0.0, subsys.acc)
               winch_friction[winch.idx] ~ subsys.friction
               D(winch_vel[winch.idx]) ~
                   ifelse(brake_p > 0.5, 0, winch_acc[winch.idx])]
        defaults = [defaults
                    bind_initial!(initial.winches[winch.idx].vel,
                                  winch_vel[winch.idx])]
    end
    return eqs, defaults, winch_subsystems
end
