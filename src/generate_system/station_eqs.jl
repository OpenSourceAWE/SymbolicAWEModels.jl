# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Station twist dynamics equation generation

"""
    validate_station_modes(stations, bodies)

Check that each station's twist mode is coherent with its owning body's dynamics
and its point count. Errors loudly on an inconsistent combination:

- `DYNAMIC` twist is an added rigid-body deformation DOF and needs
  a bridle couple, so it requires a `RIGID_DYNAMICS` body and ≥2 points.
- A 1-point station has no bridle couple to oppose a twist moment, so its only
  coherent twist is prescribed → must be `STATIC`.
- `STATIC` twist is prescribed (default 0). On a rigid wing it imposes the
  section twist; on a `PARTICLE_DYNAMICS` wing the free points carry the
  deformation, so a multi-point `STATIC` surface is an inert aero-section
  membership marker.
"""
function validate_station_modes(stations, bodies)
    for station in stations
        if station.type == KINEMATIC
            # KINEMATIC surfaces are derived (no twist DOF, no owning-body couple);
            # a flap variant needs its owning wing and an ordered [main, flap] pair.
            has_flap(station) || isempty(station.flap_body_idxs) || error(
                "Station $(station.name): flap needs exactly 2 flap " *
                "bodies [main, flap], got $(length(station.flap_body_idxs)).")
            has_flap(station) && station.wing_idx == 0 && error(
                "Station $(station.name): KINEMATIC flap needs a `wing`.")
            continue
        end
        owners = [body for body in bodies if station.idx in body.station_idxs]
        length(owners) == 1 || error(
            "Station $(station.name) is in $(length(owners)) bodies; must be in exactly 1.")
        wing = owners[1]
        rigid = wing.dynamics_type == RIGID_DYNAMICS
        npoints = length(station.point_idxs)
        if station.type == DYNAMIC
            rigid || error(
                "Station $(station.name): $(station.type) twist requires a " *
                "RIGID_DYNAMICS wing (differential/algebraic twist is a rigid " *
                "body DOF), but wing $(wing.name) is $(wing.dynamics_type).")
            npoints >= 2 || error(
                "Station $(station.name): $(station.type) twist needs a bridle couple " *
                "(≥2 points), got $npoints. A 1-point station must use STATIC twist.")
        elseif station.type == STATIC
            # STATIC = prescribed twist (default 0). On a rigid wing it imposes
            # the section twist; on a PARTICLE wing the free points carry the
            # deformation, so a multi-point STATIC surface is an inert
            # aero-section membership marker (its prescribed twist is unused).
        else
            error("Station $(station.name): unsupported twist mode $(station.type).")
        end
    end
    return nothing
end

"""
    station_eqs!(eqs, defaults, stations, bodies, params;
               R_b_to_w, fix_wing, twist_angle, twist_ω, station_aero_moment,
               point_force, station_y_airf, station_chord, station_le_pos)

Generate equations for deformable wing station twist dynamics. The couple each
of a surface's points delivers to its hinge is built from the shared
[`twist_bridle_couple`](@ref).

# Arguments
- `eqs`, `defaults`: Accumulating vectors for the MTK system.
- `stations`: Collection of Station objects (deformable wing sections).
- `bodies`: Collection of `Body` objects (the station owners).
- `R_b_to_w`: Symbolic rotation matrix (body to world).
- `fix_wing`: Symbolic boolean for fixing wing dynamics.
- `twist_angle`, `twist_ω`: Symbolic twist state variables.
- `station_aero_moment`: Symbolic aerodynamic moment on stations.
- `point_force`: Symbolic point force variable.
- `station_y_airf`, `station_chord`, `station_le_pos`: Symbolic station geometry variables.

# Returns
- Tuple `(eqs, defaults)` with updated equation vectors.
"""
function station_eqs!(eqs, defaults, stations, bodies, params, initial;
                    R_b_to_w, fix_wing, twist_angle, twist_ω, station_aero_moment,
                    point_force, station_y_airf, station_chord, station_le_pos)

    length(stations) == 0 && return eqs, defaults

    # Stations may have differing point counts (e.g. left/right halves of a
    # bridle with an asymmetric number of attachment points). Size the per-point
    # arrays by the largest point count so every station fits, and fill the
    # unused tail slots with zero equations below so the system stays balanced.
    max_npoints = maximum(length(ts.point_idxs) for ts in stations)

    @variables begin
        trailing_edge_angle(t)[eachindex(stations)]
        trailing_edge_ω(t)[eachindex(stations)]
        trailing_edge_α(t)[eachindex(stations)]
        free_twist_angle(t)[eachindex(stations)]
        twist_α(t)[eachindex(stations)]
        station_tether_force(t)[eachindex(stations)]
        station_tether_moment(t)[eachindex(stations)]
        tether_force(t)[1:max_npoints, eachindex(stations)]
        tether_moment(t)[1:max_npoints, eachindex(stations)]
        r_station(t)[1:max_npoints, eachindex(stations)]
        r_vec(t)[1:3, 1:max_npoints, eachindex(stations)]
    end

    for station in stations
        # KINEMATIC surfaces have no twist DOF; zero-bind their vars (δ is emitted by station_delta_eqs!).
        if station.type == KINEMATIC
            eqs = [
                eqs
                twist_angle[station.idx] ~ 0
                twist_ω[station.idx] ~ 0
                station_aero_moment[station.idx] ~ 0
                station_y_airf[:, station.idx] ~ zeros(3)
                station_chord[:, station.idx] ~ zeros(3)
                station_le_pos[:, station.idx] ~ zeros(3)
                station_tether_force[station.idx] ~ 0
                station_tether_moment[station.idx] ~ 0
                [tether_force[i, station.idx] ~ 0 for i in 1:max_npoints]
                [tether_moment[i, station.idx] ~ 0 for i in 1:max_npoints]
                [r_station[i, station.idx] ~ 0 for i in 1:max_npoints]
                [r_vec[j, i, station.idx] ~ 0 for i in 1:max_npoints for j in 1:3]
            ]
            continue
        end

        found = 0
        wing = nothing
        for body in bodies
            if station.idx in body.station_idxs
                wing = body
                found += 1
            end
        end
        !(found == 1) && error(
            "Kite station $(station.idx) is in $found bodies; must be in exactly 1.",
        )

        # An AeroNone owner drives no aero moment; aero_eqs! never sets it.
        no_aero = !is_wing(wing)
        no_aero &&
            (eqs = [eqs; station_aero_moment[station.idx] ~ 0])

        # Set station geometry from getters (allows runtime updates)
        eqs = [
            eqs
            station_y_airf[:, station.idx] ~ params.stations[station.idx].y_airf
            station_chord[:, station.idx] ~ params.stations[station.idx].chord
            station_le_pos[:, station.idx] ~ params.stations[station.idx].le_pos
        ]

        if station.type == STATIC
            eqs = [
                eqs
                twist_angle[station.idx] ~ params.stations[station.idx].twist
                twist_ω[station.idx] ~ 0
                station_tether_force[station.idx] ~ 0
                station_tether_moment[station.idx] ~ 0
                [tether_force[i, station.idx] ~ 0 for i in 1:max_npoints]
                [tether_moment[i, station.idx] ~ 0 for i in 1:max_npoints]
                [r_station[i, station.idx] ~ 0 for i in 1:max_npoints]
                [r_vec[j, i, station.idx] ~ 0 for i in 1:max_npoints for j in 1:3]
            ]
            (!no_aero && isempty(station.unrefined_section_idxs)) &&
                (eqs = [eqs; station_aero_moment[station.idx] ~ 0])
            continue
        end

        orientation = collect(R_b_to_w[:, :, wing.idx])
        surface_params = params.stations[station.idx]

        for (i, point_idx) in enumerate(station.point_idxs)
            couple = twist_bridle_couple(surface_params,
                                         params.points[point_idx].pos_undeformed_b,
                                         twist_angle[station.idx], orientation)
            point_load = collect(point_force[:, point_idx])
            eqs = [
                eqs
                [r_vec[j, i, station.idx] ~ couple.offset[j]
                 for j in 1:3]
                r_station[i, station.idx] ~ couple.arm
                tether_force[i, station.idx] ~ point_load ⋅ couple.direction
                tether_moment[i, station.idx] ~
                    r_station[i, station.idx] *
                    tether_force[i, station.idx]
            ]
        end

        # Zero out the unused tail rows for stations with fewer than
        # max_npoints points, so every declared array element gets exactly one
        # equation regardless of how many points this particular surface has.
        npoints = length(station.point_idxs)
        if npoints < max_npoints
            eqs = [
                eqs
                [tether_force[i, station.idx] ~ 0 for i in (npoints+1):max_npoints]
                [tether_moment[i, station.idx] ~ 0 for i in (npoints+1):max_npoints]
                [r_station[i, station.idx] ~ 0 for i in (npoints+1):max_npoints]
                [r_vec[j, i, station.idx] ~ 0 for i in (npoints+1):max_npoints for j in 1:3]
            ]
        end

        station_chord = collect(station_chord)
        station_mass = sum(params.points[point_idx].extra_mass for point_idx in station.point_idxs)
        twist = station_dynamics(;
            free_angle = free_twist_angle[station.idx],
            twist_vel = twist_ω[station.idx],
            aero_moment = station_aero_moment[station.idx],
            node_moment = station_tether_moment[station.idx],
            mass = station_mass,
            chord = station_chord[:, station.idx],
            damping = params.stations[station.idx].damping,
            stiffness = params.stations[station.idx].stiffness)

        eqs = [
            eqs
            station_tether_force[station.idx] ~ sum(tether_force[:, station.idx])
            station_tether_moment[station.idx] ~ sum(tether_moment[:, station.idx])
            twist_α[station.idx] ~ twist.twist_acc
            twist_angle[station.idx] ~ twist.angle
        ]
        if station.type == DYNAMIC
            eqs = [
                eqs
                D(free_twist_angle[station.idx]) ~
                    ifelse(fix_wing == true, 0, twist_ω[station.idx])
                D(twist_ω[station.idx]) ~
                    ifelse(fix_wing == true, 0, twist.twist_vel_rate)
            ]
            defaults = [
                defaults
                bind_initial!(initial.stations[station.idx].twist,
                              free_twist_angle[station.idx])
                bind_initial!(initial.stations[station.idx].twist_ω,
                              twist_ω[station.idx])
            ]
        else
            error("Wrong station type.")
        end
    end

    return eqs, defaults
end

"""
    station_delta_eqs!(eqs, stations; station_delta, body_R_b_to_w)

Emit the live flap deflection δ for each station into `station_delta`.
A flapped surface ([`has_flap`](@ref)) gets the signed angle between its two flap
bodies' reference chords about the world hinge axis, referenced to rest (same
formula as [`flap_delta`](@ref), evaluated symbolically from the bodies'
orientations `body_R_b_to_w`); every other surface gets 0. The flap axis,
reference chords and rest angle are frozen rest geometry baked in as constants.
"""
function station_delta_eqs!(eqs, stations;
                                  station_delta, body_R_b_to_w)
    length(stations) == 0 && return eqs
    for station in stations
        j = station.idx
        if !has_flap(station)
            eqs = [eqs; station_delta[j] ~ 0]
            continue
        end
        R_main = collect(body_R_b_to_w[:, :, station.flap_body_idxs[1]])
        R_flap = collect(body_R_b_to_w[:, :, station.flap_body_idxs[2]])
        eqs = [eqs
               station_delta[j] ~
                   flap_delta_expression(station, R_main, R_flap)]
    end
    return eqs
end
