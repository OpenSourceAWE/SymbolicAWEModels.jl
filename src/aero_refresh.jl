# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Aero refresh: the low-frequency VSM-update path (every `vsm_interval` steps).
# `refresh_aero\!` orchestrates; per-mode work is dispatched on the wing's aero
# mode via `refresh_rigid_aero\!` / `refresh_particle_aero\!`. This is NOT the
# compiled RHS, so dynamic dispatch on the abstract `wing.aero` field is free.

"""
    refresh_aero!(sam::SymbolicAWEModel, prob::ProbWithAttributes,
                  integ=sam.integrator; vsm_min_wind=0.5)

Refresh each wing's aerodynamic state, dispatching on the wing's aero mode
([`refresh_rigid_aero!`](@ref) / [`refresh_particle_aero!`](@ref)). Runs on the
low-frequency VSM-update schedule (`vsm_interval`), not the compiled RHS.

**RIGID_DYNAMICS VSM modes:** compute wind-axis coefficients (CL, CD, CS, CM, cm)
at the operating point, plus the `ForwardDiff` Jacobian over `[α, β, ω₁, ω₂, ω₃,
θ_twist…]` (`AeroLinearized`) or the frozen forces (`AeroDirect`).

**PARTICLE_DYNAMICS VSM modes:** full nonlinear VSM solve with per-point force
distribution. Non-VSM modes (`AeroNone`/`AeroPlate`) are no-ops.
"""
function refresh_aero!(sam::SymbolicAWEModel,
                       prob::ProbWithAttributes,
                       integ=sam.integrator;
                       vsm_min_wind=0.5)
    wings = sam.sys_struct.wings
    twist_surfaces = sam.sys_struct.twist_surfaces
    points = sam.sys_struct.points

    length(wings) == 0 && return nothing

    for wing in wings
        wing.dynamics_type == RIGID_DYNAMICS || continue
        refresh_rigid_aero!(wing.aero, wing, sam.am, twist_surfaces;
                            vsm_min_wind)
    end

    any(w.dynamics_type === PARTICLE_DYNAMICS for w in wings) ||
        return nothing
    point_state = prob.get_point_state(integ)
    va_point_b_vals = point_state[4]
    for wing in wings
        wing.dynamics_type == PARTICLE_DYNAMICS || continue
        refresh_particle_aero!(wing.aero, wing, points, va_point_b_vals;
                               vsm_min_wind)
    end

    nothing
end

"""
    refresh_rigid_aero!(mode, wing, am, twist_surfaces; vsm_min_wind=0.5)

Refresh a `RIGID_DYNAMICS` wing's aero state, dispatched on its aero `mode`:
- `AeroNone` / any non-VSM mode → no-op (fallback).
- `AeroLinearized` → compute the baseline coefficients ([`rigid_aero_baseline!`](@ref))
  and the `ForwardDiff` Jacobian `d(coeffs)/d(inputs)` into `wing.aero_jac`.
- `AeroDirect` → compute the baseline coefficients and apply the frozen body-frame
  force/moment; below `vsm_min_wind` everything is zeroed.
"""
refresh_rigid_aero!(::AbstractAeroModel, wing, am, twist_surfaces;
                    vsm_min_wind=0.5) = nothing

"""
    refresh_rigid_aero!(::AeroLinearized, wing, am, twist_surfaces; vsm_min_wind=0.5)

Linearized rigid-wing refresh. Computes the baseline wind-axis coefficients at the
operating point ([`rigid_aero_baseline!`](@ref)), then the `ForwardDiff` Jacobian
`d(coeffs)/d(inputs)` and stores it in `wing.aero_jac`. The compiled RHS uses that
Jacobian to reconstruct forces via a first-order Taylor expansion about the
operating point.
"""
function refresh_rigid_aero!(::AeroLinearized, wing, am, twist_surfaces;
                             vsm_min_wind=0.5)
    ctx = rigid_aero_baseline!(wing, twist_surfaces; vsm_min_wind)
    gamma0 = copy(wing.vsm_solver.sol.gamma_distribution)
    f_dual = y -> vsm_aero_coeffs(wing, y, ctx.va_mag, ctx.n_unrefined,
        ctx.n_twist_surfaces, ctx.twist_surface_idxs, twist_surfaces,
        ctx.moment_frac, ctx.shadow_ref; gamma_init=gamma0)
    ForwardDiff.jacobian!(wing.aero_jac, f_dual, ctx.y0)
    return nothing
end

"""
    refresh_rigid_aero!(::AeroDirect, wing, am, twist_surfaces; vsm_min_wind=0.5)

Direct rigid-wing refresh. Computes the baseline coefficients and applies the
resulting frozen body-frame force/moment ([`apply_direct_forces!`](@ref)), which
the RHS holds constant until the next refresh. Below `vsm_min_wind` the
coefficients, Jacobian, force, moment, and per-twist-surface moments are zeroed.
"""
function refresh_rigid_aero!(::AeroDirect, wing, am, twist_surfaces;
                             vsm_min_wind=0.5)
    if norm(wing.va_b) < vsm_min_wind
        fill!(wing.aero_x, 0.0)
        fill!(wing.aero_jac, 0.0)
        fill!(wing.aero_force_b, 0.0)
        fill!(wing.aero_moment_b, 0.0)
        for gidx in wing.twist_surface_idxs
            twist_surfaces[gidx].aero_moment = 0.0
        end
        return nothing
    end
    rigid_aero_baseline!(wing, twist_surfaces; vsm_min_wind)
    apply_direct_forces!(wing, am, wing.aero_x)
    return nothing
end

"""
    refresh_particle_aero!(mode, wing, points, va_point_b_vals; vsm_min_wind=0.5)

Refresh a `PARTICLE_DYNAMICS` wing's aero state, dispatched on its aero `mode`:
- `AeroNone` / any non-VSM mode → no-op (fallback).
- `AeroDirect` → full nonlinear VSM solve with per-section apparent wind, then
  distribute panel forces onto the wing's structural points
  ([`distribute_panel_forces_to_points!`](@ref)); below `vsm_min_wind` the point
  forces are zeroed.
- `AeroLinearized` → unsupported (errors).
"""
refresh_particle_aero!(::AbstractAeroModel, wing, points, va_point_b_vals;
                       vsm_min_wind=0.5) = nothing

"""
    refresh_particle_aero!(::AeroLinearized, wing, points, va_point_b_vals; vsm_min_wind=0.5)

Unsupported: `AeroLinearized` is not implemented for `PARTICLE_DYNAMICS` wings and
errors. Use `AeroDirect` (per-point nonlinear VSM) for particle wings.
"""
refresh_particle_aero!(::AeroLinearized, wing, points, va_point_b_vals;
                       vsm_min_wind=0.5) = error(
    "PARTICLE_DYNAMICS + AeroLinearized not yet implemented")

"""
    refresh_particle_aero!(::AeroDirect, wing, points, va_point_b_vals; vsm_min_wind=0.5)

Direct particle-wing refresh. Runs the full nonlinear VSM solve using each
section's apparent wind (averaged from its LE/TE point velocities in
`va_point_b_vals`), then distributes the resulting panel forces onto the wing's
structural points ([`distribute_panel_forces_to_points!`](@ref)). Below
`vsm_min_wind` the point forces are zeroed.
"""
function refresh_particle_aero!(::AeroDirect, wing, points, va_point_b_vals;
                                vsm_min_wind=0.5)
    if norm(wing.va_b) < vsm_min_wind
        for point in points
            if point.type == WING && point.wing_idx == wing.idx
                fill!(point.aero_force_b, 0.0)
            end
        end
        return nothing
    end

    update_vsm_wing_from_structure!(wing, points)

    if !isnothing(wing.point_to_vsm_point)
        n_sections = length(wing.vsm_wing.unrefined_sections)
        section_va = Vector{Vector{Float64}}(undef, n_sections)

        vsm_point_to_struct = Dict{Tuple{Int64, Symbol}, Int64}()
        for (point_idx, (section_idx, le_or_te)) in wing.point_to_vsm_point
            vsm_point_to_struct[(section_idx, le_or_te)] = point_idx
        end

        for section_idx in 1:n_sections
            le_pi = get(vsm_point_to_struct, (Int64(section_idx), :LE), nothing)
            te_pi = get(vsm_point_to_struct, (Int64(section_idx), :TE), nothing)
            if !isnothing(le_pi) && !isnothing(te_pi)
                va_le = va_point_b_vals[:, le_pi]
                va_te = va_point_b_vals[:, te_pi]
                section_va[section_idx] = 0.5 * (va_le + va_te)
            else
                section_va[section_idx] = wing.va_b
            end
        end

        n_panels = length(wing.vsm_aero.panels)
        va_dist = zeros(n_panels, 3)
        mapping = wing.vsm_wing.refined_panel_mapping
        for rpi in 1:n_panels
            va_dist[rpi, :] .= section_va[mapping[rpi]]
        end
        set_va!(wing.vsm_aero, va_dist)
    else
        set_va!(wing.vsm_aero, wing.va_b)
    end

    if !safe_vsm_solve!(wing.vsm_solver, wing.vsm_aero)
        throw(AssertionError("PARTICLE_DYNAMICS VSM solve failed (non-converged or non-finite) on wing $(wing.idx)"))
    end
    distribute_panel_forces_to_points!(wing, points)
    for point in points
        if point.type == WING && point.wing_idx == wing.idx &&
                any(!isfinite, point.aero_force_b)
            throw(AssertionError("PARTICLE_DYNAMICS: non-finite point force on wing $(wing.idx) point $(point.idx)"))
        end
    end
    return nothing
end

# ── RIGID_DYNAMICS aero helpers ──────────────────────────

"""
    finite_full(x) -> Bool

`true` if `x` is finite. For a `ForwardDiff.Dual` it also checks every partial, so
a NaN/Inf *derivative* is caught — used by [`safe_vsm_solve!`](@ref) to reject a
bad VSM solve during the Jacobian pass, not just the value pass.
"""
finite_full(x::Real) = isfinite(x)
finite_full(x::ForwardDiff.Dual) =
    isfinite(ForwardDiff.value(x)) &&
    all(isfinite, ForwardDiff.partials(x))

"""
NaN/Inf-guarded `solve!`. Checks both Dual value and partials. On
non-finite or non-converged result, zero gamma and return `false`.
"""
function safe_vsm_solve!(solver, body_aero,
                          gamma_init=nothing; moment_frac=0.1)
    if isnothing(gamma_init)
        VortexStepMethod.solve!(solver, body_aero;
            moment_frac, log=false)
    else
        VortexStepMethod.solve!(solver, body_aero, gamma_init;
            moment_frac, log=false)
    end
    force_coeffs = solver.sol.force_coeffs
    moment_coeffs = solver.sol.moment_coeffs
    if !solver.lr.converged ||
            any(!finite_full, force_coeffs) ||
            any(!finite_full, moment_coeffs)
        if !isnothing(solver.sol.gamma_distribution)
            fill!(solver.sol.gamma_distribution, 0)
        end
        return false
    end
    return true
end

"""
    vsm_solve_objects(wing, ::Type{T}, shadow_ref) -> (body_aero, solver, wing)

The VSM solve objects [`vsm_aero_coeffs`](@ref) runs on, selected by the input
eltype `T`. A value pass (`Float64`) uses the wing's real objects. A `ForwardDiff`
pass feeds `Dual` numbers through the solve to get the Jacobian, but the real
objects' buffers are `Float64` and can't hold Duals — so for a `Dual` eltype we
solve on a `Dual`-typed "shadow" of the solver/aero. The shadow is expensive, so
it is built lazily and cached in `shadow_ref`, keyed by the `Dual` eltype (rebuilt
if it changes); `use_gamma_prev` warm-starts each perturbed solve from the
previous circulation.
"""
vsm_solve_objects(wing, ::Type{Float64}, shadow_ref) =
    (wing.vsm_aero, wing.vsm_solver, wing.vsm_wing)

function vsm_solve_objects(wing, ::Type{T}, shadow_ref) where {T}
    shadow = shadow_ref[]
    if shadow === nothing || eltype(shadow[1]._va) !== T
        shadow = VortexStepMethod.make_dual_shadow(
            wing.vsm_solver, wing.vsm_aero, T)
        shadow[2].use_gamma_prev = true
        shadow_ref[] = shadow
    end
    body_aero, solver = shadow
    return body_aero, solver, body_aero.wings[1]
end

"""
    vsm_aero_coeffs(wing, y, va_mag, n_unrefined, n_twist_surfaces,
                     twist_surface_idxs, twist_surfaces, moment_frac, shadow_ref;
                     gamma_init=nothing) -> Vector

Run one VSM solve at operating-point input `y = [α, β, ω₁, ω₂, ω₃, θ_twist…]` and
return the wind-axis coefficient vector `[CL, CD, CS, CM₁, CM₂, CM₃, cm_twist…]`.
`ForwardDiff.Dual`-aware via `vsm_solve_objects`: for a Dual eltype it solves on a
cached dual shadow of the VSM solver, so the same routine yields the Jacobian
under AD.
"""
function vsm_aero_coeffs(wing, y::AbstractVector{T},
        va_mag, n_unrefined, n_twist_surfaces,
        twist_surface_idxs, twist_surfaces, moment_frac,
        shadow_ref::Ref;
        gamma_init=nothing) where {T}

    body_aero_c, solver_c, wing_c = vsm_solve_objects(wing, T, shadow_ref)

    α = y[1]
    β = y[2]
    ω = MVector{3, T}(y[3], y[4], y[5])

    # Body-frame apparent wind from (α, β, va_mag)
    cα, sα = cos(α), sin(α)
    cβ, sβ = cos(β), sin(β)
    va_b_local = MVector{3, T}(va_mag * cα * cβ,
                               va_mag * sβ,
                               va_mag * sα * cβ)

    # Per-twist_surface → per-section twist
    theta = zeros(T, n_unrefined)
    for (twist_surface_index, gidx) in enumerate(twist_surface_idxs)
        for unrefined_index in twist_surfaces[gidx].unrefined_section_idxs
            theta[unrefined_index] = y[5 + twist_surface_index]
        end
    end

    if n_unrefined > 0
        VortexStepMethod.unrefined_deform!(
            wing_c, theta; smooth=false)
        VortexStepMethod.reinit!(
            body_aero_c; init_aero=false)
    end
    set_va!(body_aero_c, va_b_local, ω)
    if !safe_vsm_solve!(solver_c, body_aero_c, gamma_init;
                         moment_frac)
        throw(AssertionError("VSM solve failed (non-converged or non-finite) on wing $(wing.idx) [eltype=$T]"))
    end

    sol = solver_c.sol
    force_coeffs = sol.force_coeffs
    cm_body = sol.moment_coeffs
    moment_coeff_unrefined = sol.moment_coeff_unrefined_dist

    # Wind-axis basis (matches VSM): drag along va,
    # lift = normalize(drag × span), side = lift × drag.
    span = SVector(zero(T), one(T), zero(T))
    drag_dir = va_b_local ./ va_mag
    lift_dir = smooth_normalize(cross(drag_dir, span))
    side_dir = cross(lift_dir, drag_dir)

    x = zeros(T, 6 + n_twist_surfaces)
    x[1] = dot(force_coeffs, lift_dir)
    x[2] = dot(force_coeffs, drag_dir)
    x[3] = dot(force_coeffs, side_dir)
    x[4] = cm_body[1]
    x[5] = cm_body[2]
    x[6] = cm_body[3]
    for (twist_surface_index, gidx) in enumerate(twist_surface_idxs)
        x[6 + twist_surface_index] = sum(
            moment_coeff_unrefined[unrefined_index]
            for unrefined_index in
                twist_surfaces[gidx].unrefined_section_idxs;
            init = zero(T))
    end
    return x
end

"""
    rigid_aero_baseline!(wing, twist_surfaces; vsm_min_wind=0.5)

Compute the operating point and baseline wind-axis coefficients for one wing:
writes `wing.aero_y` / `wing.aero_x` and updates `twist_surfaces[gidx].aero_moment`.
Returns the context (`va_mag`, section counts, `moment_frac`, `shadow_ref`, `y0`)
the mode-specific reduction (`refresh_rigid_aero!`) needs for the Jacobian.
"""
function rigid_aero_baseline!(wing, twist_surfaces;
                              vsm_min_wind=0.5)
    va_b = wing.va_b
    va_mag_actual = norm(va_b)
    omega_b = wing.ω_b

    twist_surface_idxs = wing.twist_surface_idxs
    n_twist_surfaces = length(twist_surface_idxs)
    n_unrefined = wing.vsm_wing.n_unrefined_sections

    moment_frac = isempty(twist_surface_idxs) ? 0.25 :
        twist_surfaces[first(twist_surface_idxs)].moment_frac

    va_mag = max(va_mag_actual, vsm_min_wind)
    alpha_0 = atan(va_b[3], va_b[1])
    beta_0 = atan(va_b[2], hypot(va_b[1], va_b[3]))
    if !isfinite(alpha_0)
        alpha_0 = 0.0
    end
    if !isfinite(beta_0)
        beta_0 = 0.0
    end

    # Operating-point input vector y₀ = [α, β, ω, θ_twist_surface]
    y0 = wing.aero_y
    y0[1] = alpha_0
    y0[2] = beta_0
    y0[3] = omega_b[1]
    y0[4] = omega_b[2]
    y0[5] = omega_b[3]
    for (twist_surface_index, gidx) in enumerate(twist_surface_idxs)
        y0[5 + twist_surface_index] = twist_surfaces[gidx].twist
    end

    shadow_ref = Ref{Any}(nothing)
    f_baseline = y -> vsm_aero_coeffs(wing, y, va_mag,
        n_unrefined, n_twist_surfaces, twist_surface_idxs, twist_surfaces,
        moment_frac, shadow_ref)

    wing.aero_x .= f_baseline(y0)
    for (twist_surface_index, gidx) in enumerate(twist_surface_idxs)
        twist_surfaces[gidx].aero_moment = wing.aero_x[6 + twist_surface_index]
    end

    return (; va_mag, n_unrefined, n_twist_surfaces,
            twist_surface_idxs, moment_frac, shadow_ref, y0)
end

"""Apply direct forces from wind-axis coefficients."""
function apply_direct_forces!(wing, am, x0)
    va_b = wing.va_b
    if any(!isfinite, x0) || any(!isfinite, va_b)
        throw(AssertionError("AeroDirect: non-finite input on wing $(wing.idx)"))
    end
    va_sq = dot(va_b, va_b)
    rho = calc_rho(am, wing.pos_w[3])
    q_inf = 0.5 * rho * va_sq
    area = wing.vsm_aero.projected_area
    c_ref = wing.vsm_aero.c_ref

    CL, CD, CS = x0[1], x0[2], x0[3]
    span = SVector(0.0, 1.0, 0.0)
    drag_dir = va_b / norm(va_b)
    lift_dir = smooth_normalize(cross(drag_dir, span))
    side_dir = cross(lift_dir, drag_dir)

    wing.aero_force_b .= q_inf * area * (
        CL .* lift_dir .+
        CD * wing.drag_frac .* drag_dir .+
        CS .* side_dir)
    wing.aero_moment_b .= q_inf * area * c_ref .*
        x0[4:6]
end
