# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Shared force-law kernels — the single source of truth every backend's point and
# segment physics consumes. Pure functions: geometry and material properties in,
# forces out. The environment (air density, apparent wind) is injected by the
# caller, so these never read the atmosphere or wind profile themselves and stay
# identical between the monolith and the `KernelBackend`.

"""
    segment_geometry(pos_src, pos_dst, vel_src, vel_dst)

Kinematics of a segment spanning `pos_src → pos_dst`. Returns
`(segment_vec, len, unit_vec, spring_vel)`, where `len` is the smoothed length,
`unit_vec` the axis, and `spring_vel = (vel_src - vel_dst) · unit_vec` the
along-axis closing speed (positive when the endpoints approach).
"""
function segment_geometry(pos_src, pos_dst, vel_src, vel_dst)
    segment_vec = pos_dst .- pos_src
    len = smooth_norm(segment_vec)
    unit_vec = segment_vec ./ len
    spring_vel = (vel_src .- vel_dst) ⋅ unit_vec
    return segment_vec, len, unit_vec, spring_vel
end

"""
    segment_spring_force(len, l0, spring_vel, unit_stiffness, unit_damping,
                         compression_frac)

Scalar spring-damper force along a segment (the Real-`unit_stiffness` branch of
`segment_eqs!`). Stiffness is `unit_stiffness / len` in tension and softens to
`compression_frac` of that under compression; the damping term
`(unit_damping / len) · spring_vel` opposes the closing speed and is deliberately
left out of that softening. Softening it would make the force jump at `len == l0`,
where the stiffness term crosses continuously only because `len - l0` vanishes
there; a taut tether sits on that crossing, so the jump diverges the solver.
Multiply by the segment `unit_vec` for the force vector on the source endpoint.
"""
function segment_spring_force(len, l0, spring_vel, unit_stiffness, unit_damping,
                              compression_frac)
    damping = unit_damping / len
    stiffness = ifelse(len > l0, unit_stiffness / len,
                       compression_frac * unit_stiffness / len)
    return stiffness * (len - l0) - damping * spring_vel
end

"""
    segment_nonlinear_force(len, l0, spring_vel, force_law, unit_damping)

Scalar force along a segment whose `unit_stiffness` is a callable force law of the
strain `(len − l0) / l0` (the callable branch of `segment_eqs!`). The law owns the
whole curve, slack and compression included, so there is no `compression_frac`; the
damping term `(unit_damping / len) · spring_vel` is the same as in the linear case.
"""
function segment_nonlinear_force(len, l0, spring_vel, force_law, unit_damping)
    return force_law((len - l0) / l0) - (unit_damping / len) * spring_vel
end

"""
    segment_perp_drag(va, unit_vec, rho, cd_tether, area)

Aerodynamic drag vector on a tether segment. Only the component of the apparent
wind `va` perpendicular to the segment axis `unit_vec` contributes:
`0.5 ρ c_d |va| A` scaling the perpendicular apparent wind. `area = len · diameter`
is the segment's projected area.
"""
function segment_perp_drag(va, unit_vec, rho, cd_tether, area)
    app_perp_vel = va .- (va ⋅ unit_vec) .* unit_vec
    return (0.5 * rho * cd_tether * smooth_norm(va) * area) .* app_perp_vel
end

"""
    point_drag_force(va, rho, drag_coeff, area)

Aerodynamic drag on a point mass from apparent wind `va`:
`0.5 ρ c_d |va| A · va`. Unlike a segment, the full apparent wind acts (no
axis projection).
"""
function point_drag_force(va, rho, drag_coeff, area)
    return (0.5 * rho * drag_coeff * smooth_norm(va) * area) .* va
end

"""
    segment_half_mass(l0, diameter, density)

Mass [kg] of half a tether segment, `density · π (diameter/2)² · l0 / 2`. Each
endpoint of a segment carries this share, so the two halves sum to the full
segment mass; a point's translational mass is `extra_mass` plus the halves of all
incident segments.
"""
function segment_half_mass(l0, diameter, density)
    mass_per_meter = density * π * (diameter / 2)^2
    return mass_per_meter * l0 / 2
end
