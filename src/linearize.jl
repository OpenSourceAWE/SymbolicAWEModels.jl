# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

"""
    find_steady_state!(s::SymbolicAWEModel, integ=s.integrator; t=1.0, dt=1/s.set.sample_freq)

Run the simulation for a short period to allow the system to settle.

During this period, the winches are braked and the wing's elevation and azimuth
angles are fixed, but it is free to move radially (in distance). This allows the
dynamic components of the bridle and tethers to settle into a stable, steady-state
equilibrium before starting a maneuver or analysis.

# Arguments
- `s::SymbolicAWEModel`: The model to be stabilized.
- `integ`: The integrator to use. Defaults to `s.integrator`.

# Keywords
- `t::Float64=1.0`: The duration [s] for which to run the settling simulation.
- `dt::Float64`: The time step [s] for the settling simulation.
"""
function find_steady_state!(sam::SymbolicAWEModel; 
                            t=2.0, dt=t/10, vsm_interval=1)
    (; winches, wings) = sam.sys_struct
    old_brakes = [winch.brake for winch in winches]
    old_fixes = [wing.fix_sphere for wing in wings]
    [winch.brake=true for winch in winches]
    [wing.fix_sphere=true for wing in wings]
    for _ in 1:Int(round(t ÷ dt))
        next_step!(sam; dt, vsm_interval)
    end
    [winch.brake=old_brakes[winch.idx] for winch in winches]
    [wing.fix_sphere=old_fixes[wing.idx] for wing in wings]
    update_sys_struct!(sam.prob, sam.integrator, sam.sys_struct)
    return nothing
end

"""
    update_vsm!(s::SymbolicAWEModel, integ=s.integrator)

Update the aerodynamic model from the Vortex Step Method.

**For QUATERNION wings:**
Computes wind-axis coefficients (CL, CD, CS, CM, cm) at the
current operating point, plus a `ForwardDiff` Jacobian over the
input vector `[α, β, ω₁, ω₂, ω₃, θ_group₁…]`. Stores the dense
Jacobian `d(coeffs)/d(inputs)` in `wing.aero_jac`.

**For REFINE wings:**
Full nonlinear VSM solve with per-point force distribution.
"""
function update_vsm!(sam::SymbolicAWEModel,
                     prob::ProbWithAttributes,
                     integ=sam.integrator;
                     vsm_min_wind=0.5)
    wings = sam.sys_struct.wings
    groups = sam.sys_struct.groups
    points = sam.sys_struct.points

    length(wings) == 0 && return nothing

    for wing in wings
        wing.wing_type != QUATERNION && continue
        wing.aero_mode == AERO_NONE && continue
        if norm(wing.va_b) < vsm_min_wind
            fill!(wing.aero_x, 0.0)
            fill!(wing.aero_jac, 0.0)
            if wing.aero_mode == AERO_DIRECT
                fill!(wing.aero_force_b, 0.0)
                fill!(wing.aero_moment_b, 0.0)
            end
            for gidx in wing.group_idxs
                groups[gidx].aero_moment = 0.0
            end
            continue
        end
        _update_quaternion_wing!(wing, sam.am, groups)
    end

    has_refine_wings = any(
        w.wing_type === REFINE for w in wings)
    if has_refine_wings
        point_state = prob.get_point_state(integ)
        va_point_b_vals = point_state[4]

        for wing in wings
            wing isa VSMWing || continue
            wing.wing_type != REFINE && continue
            wing.aero_mode == AERO_NONE && continue

            if wing.aero_mode == AERO_LINEARIZED
                error(
                    "REFINE + AERO_LINEARIZED " *
                    "not yet implemented")
            end

            if norm(wing.va_b) < vsm_min_wind
                for point in points
                    if point.type == WING &&
                            point.wing_idx == wing.idx
                        fill!(point.aero_force_b, 0.0)
                    end
                end
                continue
            end

            update_vsm_wing_from_structure!(
                wing, points)

            if !isnothing(wing.point_to_vsm_point)
                n_sections = length(
                    wing.vsm_wing.unrefined_sections)
                section_va =
                    Vector{Vector{Float64}}(
                        undef, n_sections)

                vsm_point_to_struct =
                    Dict{Tuple{Int64, Symbol}, Int64}()
                for (point_idx, (section_idx, le_or_te)) in
                        wing.point_to_vsm_point
                    vsm_point_to_struct[
                        (section_idx, le_or_te)] =
                        point_idx
                end

                for section_idx in 1:n_sections
                    le_pi = get(vsm_point_to_struct,
                        (Int64(section_idx), :LE),
                        nothing)
                    te_pi = get(vsm_point_to_struct,
                        (Int64(section_idx), :TE),
                        nothing)

                    if !isnothing(le_pi) &&
                            !isnothing(te_pi)
                        va_le =
                            va_point_b_vals[:, le_pi]
                        va_te =
                            va_point_b_vals[:, te_pi]
                        section_va[section_idx] =
                            0.5 * (va_le + va_te)
                    else
                        section_va[section_idx] =
                            wing.va_b
                    end
                end

                n_panels =
                    length(wing.vsm_aero.panels)
                va_dist = zeros(n_panels, 3)

                mapping = wing.vsm_wing.refined_panel_mapping
                for rpi in 1:n_panels
                    va_dist[rpi, :] .=
                        section_va[mapping[rpi]]
                end

                set_va!(wing.vsm_aero, va_dist)
            else
                set_va!(wing.vsm_aero, wing.va_b)
            end

            if !_safe_vsm_solve!(wing.vsm_solver, wing.vsm_aero)
                error("REFINE VSM solve failed (non-converged or non-finite) on wing $(wing.idx)")
            end
            distribute_panel_forces_to_points!(
                wing, points)
            for point in points
                if point.type == WING &&
                        point.wing_idx == wing.idx &&
                        any(!isfinite, point.aero_force_b)
                    error("REFINE: non-finite point force on wing $(wing.idx) point $(point.idx)")
                end
            end
        end
    end

    nothing
end

# ── QUATERNION aero helpers ──────────────────────────

_finite_full(x::Real) = isfinite(x)
_finite_full(x::ForwardDiff.Dual) =
    isfinite(ForwardDiff.value(x)) &&
    all(isfinite, ForwardDiff.partials(x))

"""
NaN/Inf-guarded `solve!`. Checks both Dual value and partials. On
non-finite or non-converged result, zero gamma and return `false`.
"""
function _safe_vsm_solve!(solver, body_aero,
                          gamma_init=nothing; moment_frac=0.1)
    if isnothing(gamma_init)
        VortexStepMethod.solve!(solver, body_aero;
            moment_frac, log=false)
    else
        VortexStepMethod.solve!(solver, body_aero, gamma_init;
            moment_frac, log=false)
    end
    cf = solver.sol.force_coeffs
    cm = solver.sol.moment_coeffs
    if !solver.lr.converged ||
            any(!_finite_full, cf) || any(!_finite_full, cm)
        if !isnothing(solver.sol.gamma_distribution)
            fill!(solver.sol.gamma_distribution, 0)
        end
        return false
    end
    return true
end

function _vsm_aero_coeffs(wing, y::AbstractVector{T},
        va_mag, n_unrefined, n_groups,
        group_idxs, groups, moment_frac,
        shadow_ref::Ref;
        gamma_init=nothing) where {T}

    if T === Float64
        body_aero_c = wing.vsm_aero
        solver_c = wing.vsm_solver
        wing_c = wing.vsm_wing
    else
        sh = shadow_ref[]
        if sh === nothing || eltype(sh[1]._va) !== T
            shadow_ref[] = VortexStepMethod.make_dual_shadow(
                wing.vsm_solver, wing.vsm_aero, T)
            sh = shadow_ref[]
            sh[2].use_gamma_prev = true
        end
        body_aero_c, solver_c = sh
        wing_c = body_aero_c.wings[1]
    end

    α = y[1]
    β = y[2]
    ω = MVector{3, T}(y[3], y[4], y[5])

    # Body-frame apparent wind from (α, β, va_mag)
    cα, sα = cos(α), sin(α)
    cβ, sβ = cos(β), sin(β)
    va_b_local = MVector{3, T}(va_mag * cα * cβ,
                               va_mag * sβ,
                               va_mag * sα * cβ)

    # Per-group → per-section twist
    theta = zeros(T, n_unrefined)
    for (gi, gidx) in enumerate(group_idxs)
        for ui in groups[gidx].unrefined_section_idxs
            theta[ui] = y[5 + gi]
        end
    end

    if n_unrefined > 0
        VortexStepMethod.unrefined_deform!(
            wing_c, theta; smooth=false)
        VortexStepMethod.reinit!(
            body_aero_c; init_aero=false)
    end
    set_va!(body_aero_c, va_b_local, ω)
    if !_safe_vsm_solve!(solver_c, body_aero_c, gamma_init;
                         moment_frac)
        error("VSM solve failed (non-converged or non-finite) on wing $(wing.idx) [eltype=$T]")
    end

    sol = solver_c.sol
    cf = sol.force_coeffs
    cm_body = sol.moment_coeffs
    cm_unr = sol.cm_unrefined_dist

    # Wind-axis basis (matches VSM): drag along va,
    # lift = normalize(drag × span), side = lift × drag.
    span = SVector(zero(T), one(T), zero(T))
    drag_dir = va_b_local ./ va_mag
    lift_dir = smooth_normalize(cross(drag_dir, span))
    side_dir = cross(lift_dir, drag_dir)

    x = zeros(T, 6 + n_groups)
    x[1] = dot(cf, lift_dir)
    x[2] = dot(cf, drag_dir)
    x[3] = dot(cf, side_dir)
    x[4] = cm_body[1]
    x[5] = cm_body[2]
    x[6] = cm_body[3]
    for (gi, gidx) in enumerate(group_idxs)
        x[6 + gi] = sum(cm_unr[ui]
            for ui in groups[gidx].unrefined_section_idxs;
            init = zero(T))
    end
    return x
end

"""
    _update_quaternion_wing!(wing, am, groups)

Compute baseline wind-axis coefficients and the
ForwardDiff Jacobian `d(coeffs)/d(inputs)` for one wing.

Writes `wing.aero_y / aero_x / aero_jac`, updates
`groups[gidx].aero_moment`, and (in AERO_DIRECT mode) writes
`wing.aero_force_b` / `wing.aero_moment_b`.
"""
function _update_quaternion_wing!(wing, am, groups)
    va_b = wing.va_b
    va_mag = norm(va_b)
    omega_b = wing.ω_b

    group_idxs = wing.group_idxs
    n_groups = length(group_idxs)
    n_unrefined = wing.vsm_wing.n_unrefined_sections

    moment_frac = isempty(group_idxs) ? 0.25 :
        groups[first(group_idxs)].moment_frac

    alpha_0 = atan(va_b[3], va_b[1])
    beta_0 = atan(va_b[2], hypot(va_b[1], va_b[3]))

    # Operating-point input vector y₀ = [α, β, ω, θ_group]
    y0 = wing.aero_y
    y0[1] = alpha_0
    y0[2] = beta_0
    y0[3] = omega_b[1]
    y0[4] = omega_b[2]
    y0[5] = omega_b[3]
    for (gi, gidx) in enumerate(group_idxs)
        y0[5 + gi] = groups[gidx].twist
    end

    shadow_ref = Ref{Any}(nothing)
    f_baseline = y -> _vsm_aero_coeffs(wing, y, va_mag,
        n_unrefined, n_groups, group_idxs, groups,
        moment_frac, shadow_ref)

    wing.aero_x .= f_baseline(y0)
    for (gi, gidx) in enumerate(group_idxs)
        groups[gidx].aero_moment = wing.aero_x[6 + gi]
    end

    if wing.aero_mode == AERO_LINEARIZED
        gamma0 = copy(wing.vsm_solver.sol.gamma_distribution)
        f_dual = y -> _vsm_aero_coeffs(wing, y, va_mag,
            n_unrefined, n_groups, group_idxs, groups,
            moment_frac, shadow_ref; gamma_init=gamma0)
        ForwardDiff.jacobian!(wing.aero_jac, f_dual, y0)
    elseif wing.aero_mode == AERO_DIRECT
        _apply_direct_forces!(wing, am, wing.aero_x)
    end
    return nothing
end

"""Apply direct forces from wind-axis coefficients."""
function _apply_direct_forces!(wing, am, x0)
    va_b = wing.va_b
    if any(!isfinite, x0) || any(!isfinite, va_b)
        error("AERO_DIRECT: non-finite input on wing $(wing.idx)")
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

"""
    linearize!(s::SymbolicAWEModel; set_values=s.get_set_values(s.integrator)) -> LinType

Compute the full state-space linearization of the model around the current operating point.

This function uses the `LinearizationProblem` generated by `ModelingToolkit.jl` to
calculate the A, B, C, and D matrices for the complete, high-order system.

# Arguments
- `s::SymbolicAWEModel`: The model to linearize.

# Keywords
- `set_values`: The control input vector `u` around which to linearize.

# Returns
- `LinType`: A NamedTuple `(A, B, C, D)` containing the state-space matrices.
"""
function linearize!(sam::SymbolicAWEModel; set_values=nothing)
    isnothing(sam.lin_prob) && error("Run init! with create_lin_prob=true")
    lin_prob = sam.lin_prob
    prob = sam.prob

    # copy set values from prob to lin prob
    if !isnothing(prob) && !isnothing(prob.get_set_values)
        if isnothing(set_values)
            set_values = prob.get_set_values(sam.integrator)
        end
        lin_prob.set_set_values(lin_prob.prob, set_values)
    end

    lin_model = solve(lin_prob.prob)[1]
    return lin_model
end

"""
    jacobian(f::Function, x::AbstractVector, ϵ::AbstractVector) -> Matrix

Numerically compute the Jacobian of a vector-valued function `f` at point `x`.

This function uses a simple forward finite difference method to approximate the partial
derivatives of `f` with respect to each component of `x`.

# Arguments
- `f::Function`: The function to differentiate (`y = f(x)`).
- `x::AbstractVector`: The point at which to evaluate the Jacobian.
- `ϵ::AbstractVector`: A vector of perturbation sizes for each component of `x`.

# Returns
- `Matrix`: The Jacobian matrix `J`, where `J[i, j] = ∂f[i] / ∂x[j]`.
"""
function jacobian(f::Function, x::AbstractVector, ϵ::AbstractVector)
    n = length(x)
    fx = f(x)
    m = length(fx)
    J = zeros(m, n)
    for i in 1:n
        x_perturbed = copy(x)
        x_perturbed[i] += ϵ[i]
        J[:, i] = (f(x_perturbed) - fx) / ϵ[i]
    end
    return J
end

