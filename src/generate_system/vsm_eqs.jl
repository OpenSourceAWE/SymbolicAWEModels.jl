# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

function vsm_eqs!(
    s, eqs, guesses, psys;
    aero_force_b, aero_moment_b, group_aero_moment,
    twist_angle, va_wing_b, wing_pos, ω_b,
    aero_force_point_b=nothing
)
    (; groups, wings, points) = s.sys_struct
    length(wings) == 0 && return eqs, guesses

    for wing in wings
        if wing isa PlateWing
            continue
        elseif wing isa VSMWing && wing.dynamics_type == PARTICLE_DYNAMICS
            eqs = particle_wing_aero_eqs!(
                eqs, wing, points, psys;
                aero_force_b, aero_moment_b, aero_force_point_b)
        elseif wing.aero_type isa NoAero
            eqs = no_aero_wing_eqs!(
                eqs, wing, groups;
                aero_force_b, aero_moment_b, group_aero_moment)
        else
            eqs = rigid_wing_aero_eqs!(
                eqs, wing, groups, psys;
                aero_force_b, aero_moment_b, group_aero_moment,
                va_wing_b, ω_b, twist_angle)
        end
    end
    return eqs, guesses
end

function no_aero_wing_eqs!(
    eqs, wing, groups;
    aero_force_b, aero_moment_b, group_aero_moment
)
    eqs = [
        eqs
        aero_force_b[:, wing.idx] ~ zeros(3)
        aero_moment_b[:, wing.idx] ~ zeros(3)
    ]
    for gidx in wing.group_idxs
        isempty(groups[gidx].unrefined_section_idxs) && continue
        eqs = [eqs; group_aero_moment[groups[gidx].idx] ~ 0]
    end
    return eqs
end

function rigid_wing_aero_eqs!(
    eqs, wing, groups, psys;
    aero_force_b, aero_moment_b, group_aero_moment,
    va_wing_b, ω_b, twist_angle
)
    w = wing.idx
    twists = [twist_angle[groups[gidx].idx] for gidx in wing.group_idxs]
    va1, va2, va3 = va_wing_b[1, w], va_wing_b[2, w], va_wing_b[3, w]
    ω1, ω2, ω3 = ω_b[1, w], ω_b[2, w], ω_b[3, w]

    force_moment = wing_force_moment(
        psys, w, va1, va2, va3, ω1, ω2, ω3, twists)
    eqs = [
        eqs
        aero_force_b[:, w] ~
            [force_moment[1], force_moment[2], force_moment[3]]
        aero_moment_b[:, w] ~
            [force_moment[4], force_moment[5], force_moment[6]]
    ]

    for (gi, gidx) in enumerate(wing.group_idxs)
        isempty(groups[gidx].unrefined_section_idxs) && continue
        eqs = [
            eqs
            group_aero_moment[groups[gidx].idx] ~
                group_twist_moment(
                    psys, w, gi, va1, va2, va3, ω1, ω2, ω3, twists)
        ]
    end
    return eqs
end

function particle_wing_aero_eqs!(
    eqs, wing, points, psys;
    aero_force_b, aero_moment_b, aero_force_point_b
)
    point_force_b = aero_force_point_b::AbstractArray
    wing_points = [
        p for p in points
        if p.type == WING && p.wing_idx == wing.idx]

    if wing.aero_type isa NoAero
        for point in wing_points
            eqs = [eqs; point_force_b[:, point.idx] ~ zeros(3)]
        end
    else
        for point in wing_points
            eqs = [
                eqs
                point_force_b[:, point.idx] ~ [
                    get_point_aero_force(psys, point.idx, i)
                    for i in 1:3]
            ]
        end
    end

    eqs = [
        eqs
        aero_force_b[:, wing.idx] ~
            sum([point_force_b[:, p.idx] for p in wing_points])
        aero_moment_b[:, wing.idx] ~ zeros(3)
    ]
    return eqs
end
