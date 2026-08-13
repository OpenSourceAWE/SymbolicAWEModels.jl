# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Segment spring-damper equation generation

"""
    segment_rest_length_eqs(segment, pulleys, tethers, params;
                            l0, pulley_len, pulley_vel, tether_len, tether_vel)

Rest-length equations for one segment, and the rate that rest length moves at.
A pulley member takes either the pulley state `pulley_len` or the remainder
`sum_len - pulley_len`, a tether member takes `tether_len / n_segments`, and any
other segment keeps its fixed `l0` parameter. Returns `(eqs, rest_len_rate)`, the
rate carrying the sign of the state it follows and being zero for a fixed rest
length; [`segment_load_terms`](@ref) needs it to damp the extension rather than the
endpoint separation. Errors when a segment belongs to more than one pulley or more
than one tether.
"""
function segment_rest_length_eqs(segment, pulleys, tethers, params;
                                 l0, pulley_len, pulley_vel, tether_len, tether_vel)
    idx = segment.idx
    eqs = Equation[]
    rate = Num(0.0)
    in_pulley = 0
    for pulley in pulleys
        split_rate = pulley_len_rate(params.pulleys[pulley.idx],
                                     pulley_vel[pulley.idx])
        if idx == pulley.segment_idxs[1]
            push!(eqs, l0[idx] ~ pulley_len[pulley.idx])
            rate = split_rate
            in_pulley += 1
        end
        if idx == pulley.segment_idxs[2]
            push!(eqs, l0[idx] ~
                params.pulleys[pulley.idx].sum_len - pulley_len[pulley.idx])
            rate = -split_rate
            in_pulley += 1
        end
    end
    (in_pulley > 1) && error(
        "Bridle segment $idx is in $in_pulley pulleys; should be in 0 or 1.",
    )
    in_pulley == 1 && return eqs, rate

    in_tether = 0
    tether_idx = 0
    for tether in tethers
        if idx in tether.segment_idxs
            tether_idx = tether.idx
            in_tether += 1
        end
    end
    !(in_tether in [0, 1]) && error(
        "Segment $idx is in $in_tether tethers; should be 0 or 1.",
    )
    if in_tether == 1
        n_segments = length(tethers[tether_idx].segment_idxs)
        push!(eqs, l0[idx] ~ tether_len[tether_idx] / n_segments)
        rate = tether_vel[tether_idx] / n_segments
    else
        push!(eqs, l0[idx] ~ params.segments[idx].l0)
    end
    return eqs, rate
end

"""
    segment_eqs!(s, eqs, points, segments, pulleys, tethers, bodies, params;
                 pos, vel, wind_vec_gnd, spring_force_vec, drag_force, l0,
                 pulley_len, pulley_vel, tether_len, tether_vel)

Generate equations for segment spring-damper forces and aerodynamic drag.

Every load term comes from the shared [`segment_load_terms`](@ref) with the
parameters read by [`segment_spring_params`](@ref), so the monolith and the
`KernelBackend` evaluate the same force law. A
[`wing_structural_segment`](@ref) gets `with_drag = false` (its aerodynamic load is
owned by the wing's VSM), and one whose wing has `RIGID_DYNAMICS` skips the spring
as well, keeping only the geometry.

# Arguments
- `s::SymbolicAWEModel`: The main model object (for atmospheric model).
- `eqs`: Accumulating equation vector for the MTK system.
- `points`, `segments`, `pulleys`, `tethers`, `bodies`: System components.
- `pos`, `vel`: Symbolic point state variables.
- `wind_vec_gnd`: Symbolic ground-level wind vector.
- `spring_force_vec`, `drag_force`, `l0`: Pre-declared segment force variables.
- `pulley_len`, `tether_len`: Symbolic state variables for pulley and tether lengths.
- `pulley_vel`, `tether_vel`: Their rates, which a moving rest length's damper reads
  through [`segment_rest_length_eqs`](@ref).

# Returns
- Tuple `(eqs, len, spring_force)` with updated equation vector
  and the segment length and spring force variables for use by other components.
"""
function segment_eqs!(s, eqs, points, segments,
                      pulleys, tethers, bodies,
                      params; pos, vel, wind_vec_gnd,
                      spring_force_vec, drag_force, l0,
                      pulley_len, pulley_vel, tether_len, tether_vel)
    @variables begin
        segment_vec(t)[1:3, eachindex(segments)]
        unit_vec(t)[1:3, eachindex(segments)]
        len(t)[eachindex(segments)]
        rel_vel(t)[1:3, eachindex(segments)]
        spring_vel(t)[eachindex(segments)]
        spring_force(t)[eachindex(segments)]
    end
    wind_gnd = collect(wind_vec_gnd)

    for segment in segments
        idx = segment.idx
        src, dst = segment.point_idxs[1], segment.point_idxs[2]
        src_pos, src_vel = collect(pos[:, src]), collect(vel[:, src])
        dst_pos, dst_vel = collect(pos[:, dst]), collect(vel[:, dst])

        rest_eqs, rest_len_rate = segment_rest_length_eqs(
            segment, pulleys, tethers, params;
            l0, pulley_len, pulley_vel, tether_len, tether_vel)
        eqs = [
            eqs
            rest_eqs
            rel_vel[:, idx] ~ src_vel - dst_vel
        ]

        with_drag = !wing_structural_segment(params.reg.sys_struct, idx)
        rigid_link = !with_drag &&
            bodies[points[src].wing_idx].dynamics_type == RIGID_DYNAMICS
        if rigid_link
            seg_vec, seg_len, axis, closing_vel =
                segment_geometry(src_pos, dst_pos, src_vel, dst_vel)
            eqs = [
                eqs
                segment_vec[:, idx] ~ seg_vec
                len[idx] ~ seg_len
                unit_vec[:, idx] ~ axis
                spring_vel[idx] ~ closing_vel
                spring_force[idx] ~ 0.0
                spring_force_vec[:, idx] ~ zeros(3)
                drag_force[:, idx] ~ zeros(3)
            ]
            continue
        end

        spring, _, wind_factor = segment_spring_params(params, idx; with_drag)
        loads = segment_load_terms(s, src_pos, src_vel, dst_pos, dst_vel,
            spring.unit_stiffness, spring.unit_damping, spring.compression_frac,
            spring.compression_damping_frac, l0[idx], spring.diameter,
            spring.density, spring.cd_tether, wind_gnd, wind_factor;
            with_drag, nonlinear = spring.nonlinear, rest_len_rate)
        eqs = [
            eqs
            segment_vec[:, idx] ~ loads.segment_vec
            len[idx] ~ loads.len
            unit_vec[:, idx] ~ loads.unit_vec
            spring_vel[idx] ~ loads.spring_vel
            spring_force[idx] ~ loads.spring
            spring_force_vec[:, idx] ~ loads.spring_vec
            drag_force[:, idx] ~ 2 .* loads.half_drag
        ]
    end

    return eqs, len, spring_force
end
