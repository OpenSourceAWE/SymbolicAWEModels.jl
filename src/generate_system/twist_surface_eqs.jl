# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# TwistSurface twist dynamics equation generation

"""
    validate_twist_surface_modes(twist_surfaces, bodies)

Check that each twist_surface's twist mode is coherent with its owning body's dynamics
and its point count. Errors loudly on an inconsistent combination:

- `DYNAMIC` twist is an added rigid-body deformation DOF and needs
  a bridle couple, so it requires a `RIGID_DYNAMICS` body and ≥2 points.
- A 1-point twist_surface has no bridle couple to oppose a twist moment, so its only
  coherent twist is prescribed → must be `STATIC`.
- `STATIC` twist on a `PARTICLE_DYNAMICS` body cannot move free particles, so it is
  only coherent for a single point.
"""
function validate_twist_surface_modes(twist_surfaces, bodies)
    for twist_surface in twist_surfaces
        owners = [body for body in bodies if twist_surface.idx in body.twist_surface_idxs]
        length(owners) == 1 || error(
            "TwistSurface $(twist_surface.name) is in $(length(owners)) bodies; must be in exactly 1.")
        wing = owners[1]
        rigid = wing.dynamics_type == RIGID_DYNAMICS
        npoints = length(twist_surface.point_idxs)
        if twist_surface.type == DYNAMIC
            rigid || error(
                "TwistSurface $(twist_surface.name): $(twist_surface.type) twist requires a " *
                "RIGID_DYNAMICS wing (differential/algebraic twist is a rigid " *
                "body DOF), but wing $(wing.name) is $(wing.dynamics_type).")
            npoints >= 2 || error(
                "TwistSurface $(twist_surface.name): $(twist_surface.type) twist needs a bridle couple " *
                "(≥2 points), got $npoints. A 1-point twist_surface must use STATIC twist.")
        elseif twist_surface.type == STATIC
            (rigid || npoints == 1) || error(
                "TwistSurface $(twist_surface.name): STATIC twist on a PARTICLE_DYNAMICS wing is " *
                "only coherent for a single point (imposed twist cannot move free " *
                "particles), got $npoints points.")
        else
            error("TwistSurface $(twist_surface.name): unsupported twist mode $(twist_surface.type).")
        end
    end
    return nothing
end

"""
    twist_surface_eqs!(eqs, defaults, twist_surfaces, bodies, params;
               R_b_to_w, fix_wing, twist_angle, twist_ω, twist_surface_aero_moment,
               point_force, twist_surface_y_airf, twist_surface_chord, twist_surface_le_pos)

Generate equations for deformable wing twist_surface twist dynamics.

# Arguments
- `eqs`, `defaults`: Accumulating vectors for the MTK system.
- `twist_surfaces`: Collection of TwistSurface objects (deformable wing sections).
- `bodies`: Collection of `Body` objects (the twist_surface owners).
- `R_b_to_w`: Symbolic rotation matrix (body to world).
- `fix_wing`: Symbolic boolean for fixing wing dynamics.
- `twist_angle`, `twist_ω`: Symbolic twist state variables.
- `twist_surface_aero_moment`: Symbolic aerodynamic moment on twist_surfaces.
- `point_force`: Symbolic point force variable.
- `twist_surface_y_airf`, `twist_surface_chord`, `twist_surface_le_pos`: Symbolic twist_surface geometry variables.

# Returns
- Tuple `(eqs, defaults)` with updated equation vectors.
"""
function twist_surface_eqs!(eqs, defaults, twist_surfaces, bodies, params, initial;
                    R_b_to_w, fix_wing, twist_angle, twist_ω, twist_surface_aero_moment,
                    point_force, twist_surface_y_airf, twist_surface_chord, twist_surface_le_pos)

    length(twist_surfaces) == 0 && return eqs, defaults

    # Twist surfaces may have differing point counts (e.g. left/right halves of a
    # bridle with an asymmetric number of attachment points). Size the per-point
    # arrays by the largest point count so every twist_surface fits, and fill the
    # unused tail slots with zero equations below so the system stays balanced.
    max_npoints = maximum(length(ts.point_idxs) for ts in twist_surfaces)

    @variables begin
        trailing_edge_angle(t)[eachindex(twist_surfaces)]
        trailing_edge_ω(t)[eachindex(twist_surfaces)]
        trailing_edge_α(t)[eachindex(twist_surfaces)]
        free_twist_angle(t)[eachindex(twist_surfaces)]
        twist_α(t)[eachindex(twist_surfaces)]
        twist_surface_tether_force(t)[eachindex(twist_surfaces)]
        twist_surface_tether_moment(t)[eachindex(twist_surfaces)]
        tether_force(t)[1:max_npoints, eachindex(twist_surfaces)]
        tether_moment(t)[1:max_npoints, eachindex(twist_surfaces)]
        r_twist_surface(t)[1:max_npoints, eachindex(twist_surfaces)]
        r_vec(t)[1:3, 1:max_npoints, eachindex(twist_surfaces)]
    end

    for twist_surface in twist_surfaces
        found = 0
        wing = nothing
        for body in bodies
            if twist_surface.idx in body.twist_surface_idxs
                wing = body
                found += 1
            end
        end
        !(found == 1) && error(
            "Kite twist_surface $(twist_surface.idx) is in $found bodies; must be in exactly 1.",
        )

        # An AeroNone owner drives no aero moment; aero_eqs! never sets it.
        no_aero = !is_wing(wing)
        no_aero &&
            (eqs = [eqs; twist_surface_aero_moment[twist_surface.idx] ~ 0])

        # Set twist_surface geometry from getters (allows runtime updates)
        eqs = [
            eqs
            twist_surface_y_airf[:, twist_surface.idx] ~ params.twist_surfaces[twist_surface.idx].y_airf
            twist_surface_chord[:, twist_surface.idx] ~ params.twist_surfaces[twist_surface.idx].chord
            twist_surface_le_pos[:, twist_surface.idx] ~ params.twist_surfaces[twist_surface.idx].le_pos
        ]

        if twist_surface.type == STATIC
            eqs = [
                eqs
                twist_angle[twist_surface.idx] ~ params.twist_surfaces[twist_surface.idx].twist
                twist_ω[twist_surface.idx] ~ 0
                twist_surface_tether_force[twist_surface.idx] ~ 0
                twist_surface_tether_moment[twist_surface.idx] ~ 0
                [tether_force[i, twist_surface.idx] ~ 0 for i in 1:max_npoints]
                [tether_moment[i, twist_surface.idx] ~ 0 for i in 1:max_npoints]
                [r_twist_surface[i, twist_surface.idx] ~ 0 for i in 1:max_npoints]
                [r_vec[j, i, twist_surface.idx] ~ 0 for i in 1:max_npoints for j in 1:3]
            ]
            (!no_aero && isempty(twist_surface.unrefined_section_idxs)) &&
                (eqs = [eqs; twist_surface_aero_moment[twist_surface.idx] ~ 0])
            continue
        end

        gc = collect(twist_surface_chord[:, twist_surface.idx])
        x_airf = smooth_normalize(gc)
        gy = collect(twist_surface_y_airf[:, twist_surface.idx])
        init_z_airf = x_airf × gy
        z_airf = sin(twist_angle[twist_surface.idx]) * x_airf + cos(twist_angle[twist_surface.idx]) * init_z_airf
        Rbw = collect(R_b_to_w[:, :, wing.idx])
        Rz = Rbw * (-1 * z_airf)  # Note: -z_airf has a bug, use -1 * z_airf instead
        gl = collect(twist_surface_le_pos[:, twist_surface.idx])

        for (i, point_idx) in enumerate(twist_surface.point_idxs)
            pf = collect(point_force[:, point_idx])
            rv = collect(r_vec[:, i, twist_surface.idx])
            pos_offset = collect(
                params.points[point_idx].pos_b .-
                (gl + params.twist_surfaces[twist_surface.idx].moment_frac * gc)
            )
            eqs = [
                eqs
                [r_vec[j, i, twist_surface.idx] ~ pos_offset[j]
                 for j in 1:3]
                r_twist_surface[i, twist_surface.idx] ~ rv ⋅ smooth_normalize(gc)
                tether_force[i, twist_surface.idx] ~ pf ⋅ Rz
                tether_moment[i, twist_surface.idx] ~ r_twist_surface[i, twist_surface.idx] * tether_force[i, twist_surface.idx]
            ]
        end

        # Zero out the unused tail rows for twist_surfaces with fewer than
        # max_npoints points, so every declared array element gets exactly one
        # equation regardless of how many points this particular surface has.
        npoints = length(twist_surface.point_idxs)
        if npoints < max_npoints
            eqs = [
                eqs
                [tether_force[i, twist_surface.idx] ~ 0 for i in (npoints+1):max_npoints]
                [tether_moment[i, twist_surface.idx] ~ 0 for i in (npoints+1):max_npoints]
                [r_twist_surface[i, twist_surface.idx] ~ 0 for i in (npoints+1):max_npoints]
                [r_vec[j, i, twist_surface.idx] ~ 0 for i in (npoints+1):max_npoints for j in 1:3]
            ]
        end

        # Thin-plate inertia about one edge: I = 1/3·m·L² (m = twist_surface mass).
        twist_surface_chord = collect(twist_surface_chord)
        twist_surface_mass = sum(params.points[point_idx].extra_mass for point_idx in twist_surface.point_idxs)
        inertia = 1 / 3 * twist_surface_mass * smooth_norm(twist_surface_chord[:, twist_surface.idx])^2
        max_twist = deg2rad(90)

        eqs = [
            eqs
            twist_surface_tether_force[twist_surface.idx] ~ sum(tether_force[:, twist_surface.idx])
            twist_surface_tether_moment[twist_surface.idx] ~ sum(tether_moment[:, twist_surface.idx])
            twist_α[twist_surface.idx] ~
                (twist_surface_aero_moment[twist_surface.idx] + twist_surface_tether_moment[twist_surface.idx]) /
                inertia
            twist_angle[twist_surface.idx] ~
                clamp(free_twist_angle[twist_surface.idx], -max_twist, max_twist)
        ]
        if twist_surface.type == DYNAMIC
            eqs = [
                eqs
                D(free_twist_angle[twist_surface.idx]) ~
                    ifelse(fix_wing == true, 0, twist_ω[twist_surface.idx])
                D(twist_ω[twist_surface.idx]) ~ ifelse(
                    fix_wing == true,
                    0,
                    twist_α[twist_surface.idx] -
                    params.twist_surfaces[twist_surface.idx].damping * twist_ω[twist_surface.idx] -
                    params.twist_surfaces[twist_surface.idx].stiffness * twist_angle[twist_surface.idx] /
                    inertia,
                )
            ]
            defaults = [
                defaults
                bind_initial!(initial.twist_surfaces[twist_surface.idx].twist,
                              free_twist_angle[twist_surface.idx])
                bind_initial!(initial.twist_surfaces[twist_surface.idx].twist_ω,
                              twist_ω[twist_surface.idx])
            ]
        else
            error("Wrong twist_surface type.")
        end
    end

    return eqs, defaults
end
