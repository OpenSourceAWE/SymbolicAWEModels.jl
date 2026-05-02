# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

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
current operating point via a baseline VSM solve, plus
finite-difference perturbation solves for each input direction
(alpha, beta, omega, twist). Stores the dense Jacobian
`d(coeffs)/d(inputs)` in `wing.aero_jac`.

**For REFINE wings:**
Full nonlinear VSM solve with per-point force distribution.
"""
function update_vsm!(sam::SymbolicAWEModel,
                     prob::ProbWithAttributes,
                     integ=sam.integrator)
    wings = sam.sys_struct.wings
    groups = sam.sys_struct.groups
    points = sam.sys_struct.points

    length(wings) == 0 && return nothing

    # Handle QUATERNION wings
    has_quaternion_wings = any(
        w.wing_type === QUATERNION for w in wings)
    if has_quaternion_wings
        for wing in wings
            wing.wing_type != QUATERNION && continue
            wing.aero_mode == AERO_NONE && continue

            va_b = wing.va_b
            va_mag = norm(va_b)
            if va_mag < 1e-3
                @warn "Apparent wind too small " *
                    "(va_b=$va_b), skipping " *
                    "wing $(wing.idx)"
                continue
            end
            if any(isnan.(wing.vsm_solver.sol.force))
                wing.vsm_solver.prob = nothing
                @warn "Resetting vsm solver."
            end

            group_idxs = wing.group_idxs
            n_groups = length(group_idxs)
            n_unrefined =
                wing.vsm_wing.n_unrefined_sections

            moment_frac = if isempty(group_idxs)
                0.25
            elseif length(groups) >=
                    maximum(group_idxs)
                groups[first(group_idxs)].moment_frac
            else
                0.25
            end

            # Current twist per unrefined section
            theta_0 = zeros(n_unrefined)
            group_twist_0 = zeros(n_groups)
            for (gi, gidx) in enumerate(group_idxs)
                g = groups[gidx]
                group_twist_0[gi] = g.twist
                for ui in g.unrefined_section_idxs
                    theta_0[ui] = g.twist
                end
            end

            # Operating point
            omega_b = copy(wing.ω_b)
            alpha_0 = atan(va_b[3], va_b[1])
            beta_0 = asin(clamp(
                va_b[2] / va_mag, -1, 1))

            # Baseline solve
            x0 = _vsm_solve_coeffs!(
                wing, theta_0, va_b, omega_b,
                moment_frac, groups)

            # Store operating point
            _store_operating_point!(
                wing, alpha_0, beta_0, omega_b,
                group_twist_0, groups)

            # Store baseline coefficients
            wing.aero_x .= x0

            if wing.aero_mode == AERO_LINEARIZED
                # Perturbation solves → Jacobian
                _compute_aero_jacobian!(
                    wing, x0, alpha_0, beta_0,
                    va_mag, omega_b, theta_0,
                    group_twist_0, moment_frac,
                    groups)
            elseif wing.aero_mode == AERO_DIRECT
                # Compute physical forces from
                # baseline wind-axis coefficients
                _apply_direct_forces!(
                    wing, sam.am, x0)
            end

            # Update group aero moments
            for (gi, gidx) in enumerate(group_idxs)
                groups[gidx].aero_moment = x0[6 + gi]
            end
        end
    end

    # Handle REFINE wings (full nonlinear solve)
    has_refine_wings = any(
        w.wing_type === REFINE for w in wings)
    if has_refine_wings
        point_state = prob.get_point_state(integ)
        va_point_b_vals = point_state[4]

        for wing in wings
            wing.wing_type != REFINE && continue
            wing.aero_mode == AERO_NONE && continue

            if wing.aero_mode == AERO_LINEARIZED
                error(
                    "REFINE + AERO_LINEARIZED " *
                    "not yet implemented")
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

            VortexStepMethod.solve!(
                wing.vsm_solver, wing.vsm_aero;
                log=false)
            distribute_panel_forces_to_points!(
                wing, points)
        end
    end

    nothing
end

# ── QUATERNION aero helpers ──────────────────────────

"""Convert (alpha, beta, va_mag) to body-frame va_b."""
function _va_from_angles(alpha, beta, va_mag)
    return [va_mag * cos(alpha) * cos(beta),
            va_mag * sin(beta),
            va_mag * sin(alpha) * cos(beta)]
end

"""
    _vsm_solve_coeffs!(wing, theta, va_b, omega_b,
                       moment_frac, groups)

Single VSM solve → wind-axis coefficient vector:
`[CL, CD, CS, CM1, CM2, CM3, cm_1..cm_n]`
"""
function _vsm_solve_coeffs!(
    wing, theta, va_b, omega_b,
    moment_frac, groups
)
    VortexStepMethod.unrefined_deform!(
        wing.vsm_wing, theta; smooth=false)
    VortexStepMethod.reinit!(
        wing.vsm_aero; init_aero=false)
    set_va!(wing.vsm_aero, va_b, omega_b)
    VortexStepMethod.solve!(
        wing.vsm_solver, wing.vsm_aero;
        moment_frac, log=false)

    sol = wing.vsm_solver.sol
    cf = sol.force_coeffs
    cm_body = sol.moment_coeffs
    cm_unr = sol.cm_unrefined_dist

    # Wind-axis decomposition
    va_norm = va_b / norm(va_b)
    side_dir = SVector(0.0, 1.0, 0.0)
    lift_dir = normalize(cross(va_norm, side_dir))

    n_groups = length(wing.group_idxs)
    x = zeros(6 + n_groups)
    x[1] = dot(cf, lift_dir)   # CL
    x[2] = dot(cf, va_norm)    # CD
    x[3] = dot(cf, side_dir)   # CS
    x[4] = cm_body[1]          # CM1
    x[5] = cm_body[2]          # CM2
    x[6] = cm_body[3]          # CM3

    for (gi, gidx) in enumerate(wing.group_idxs)
        g = groups[gidx]
        x[6 + gi] = sum(
            cm_unr[ui]
            for ui in g.unrefined_section_idxs)
    end
    return x
end

"""Store operating point in `wing.aero_y`.

Layout: [α, β, ω₁, ω₂, ω₃, θ_group₁…θ_groupₙ]
"""
function _store_operating_point!(
    wing, alpha, beta, omega_b,
    group_twist, groups
)
    wing.aero_y[1] = alpha
    wing.aero_y[2] = beta
    wing.aero_y[3] = omega_b[1]
    wing.aero_y[4] = omega_b[2]
    wing.aero_y[5] = omega_b[3]

    for gi in eachindex(group_twist)
        wing.aero_y[5 + gi] = group_twist[gi]
    end
end

"""
Compute dense Jacobian via finite-difference
perturbation solves. Columns are input perturbation
directions, rows are all output coefficients.

Input layout: [α, β, ω₁, ω₂, ω₃, θ_group₁…θ_groupₙ]
"""
function _compute_aero_jacobian!(
    wing, x0, alpha_0, beta_0,
    va_mag, omega_b, theta_0,
    group_twist_0, moment_frac, groups
)
    delta_angle = 1e-2   # rad
    delta_omega = 1e-2   # rad/s
    n_groups = length(wing.group_idxs)
    ny = length(wing.aero_y)
    jac = wing.aero_jac

    for col in 1:ny
        va_b = _va_from_angles(
            alpha_0, beta_0, va_mag)
        omega = copy(omega_b)
        theta = copy(theta_0)

        if col == 1       # alpha
            va_b = _va_from_angles(
                alpha_0 + delta_angle,
                beta_0, va_mag)
            d = delta_angle
        elseif col == 2   # beta
            va_b = _va_from_angles(
                alpha_0,
                beta_0 + delta_angle, va_mag)
            d = delta_angle
        elseif col <= 5   # omega_b[col-2]
            omega[col - 2] += delta_omega
            d = delta_omega
        else              # per-group twist
            gi = col - 5
            gidx = wing.group_idxs[gi]
            for ui in groups[gidx].unrefined_section_idxs
                theta[ui] += delta_angle
            end
            d = delta_angle
        end

        x_pert = _vsm_solve_coeffs!(
            wing, theta, va_b, omega,
            moment_frac, groups)
        jac[:, col] .= (x_pert .- x0) ./ d
    end

    # Restore baseline state on VSM
    VortexStepMethod.unrefined_deform!(
        wing.vsm_wing, theta_0; smooth=false)
    VortexStepMethod.reinit!(
        wing.vsm_aero; init_aero=false)
    set_va!(wing.vsm_aero,
        _va_from_angles(alpha_0, beta_0, va_mag),
        omega_b)
end

"""Apply direct forces from wind-axis coefficients."""
function _apply_direct_forces!(wing, am, x0)
    va_b = wing.va_b
    va_sq = dot(va_b, va_b)
    rho = calc_rho(am, wing.pos_w[3])
    q_inf = 0.5 * rho * va_sq
    area = wing.vsm_aero.projected_area
    c_ref = wing.vsm_aero.c_ref

    CL, CD, CS = x0[1], x0[2], x0[3]
    va_norm = va_b / norm(va_b)
    side_dir = SVector(0.0, 1.0, 0.0)
    lift_dir = normalize(cross(va_norm, side_dir))

    wing.aero_force_b .= q_inf * area * (
        CL .* lift_dir .-
        CD * wing.drag_frac .* va_norm .+
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

