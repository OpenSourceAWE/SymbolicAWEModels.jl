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

Update the aerodynamic model from the Vortex Step Method (VSM).

This function updates the VSM aerodynamics for all wings, with wing-type-specific behavior:

**For QUATERNION wings:**
- Takes the current kinematic state (apparent wind, angular velocity, twist angles)
- Linearizes the VSM aerodynamics around this operating point
- Updates the Jacobian (`vsm_jac`) and steady-state forces (`vsm_x`)

**For REFINE wings:**
- Updates VSM panel positions from current structural deformation
- Solves the full nonlinear VSM system
- Distributes panel forces to structural points via `point.aero_force_b`

This is typically called periodically during simulation based on the `vsm_interval` parameter.
"""
function update_vsm!(sam::SymbolicAWEModel, prob::ProbWithAttributes,
                     integ=sam.integrator; vsm_min_wind=0.5)
    wings = sam.sys_struct.wings
    groups = sam.sys_struct.groups
    points = sam.sys_struct.points

    length(wings) == 0 && return nothing

    # Handle QUATERNION wings
    has_quaternion_wings = any(
        w.wing_type === QUATERNION for w in wings)
    if has_quaternion_wings && !isnothing(prob.get_vsm_y)
        vsm_y = prob.get_vsm_y(integ)

        for wing in wings
            wing.wing_type != QUATERNION && continue
            wing.aero_mode == AERO_NONE && continue

            wing.vsm_y .= vsm_y[:, wing.idx]
            if norm(wing.vsm_y[1:3]) < vsm_min_wind
                fill!(wing.vsm_x, 0.0)
                fill!(wing.vsm_jac, 0.0)
                if wing.aero_mode == AERO_DIRECT
                    fill!(wing.aero_force_b, 0.0)
                    fill!(wing.aero_moment_b, 0.0)
                end
                for gidx in wing.group_idxs
                    groups[gidx].aero_moment = 0.0
                end
                continue
            end
            if any(isnan.(wing.vsm_solver.sol.force))
                wing.vsm_solver.prob = nothing
                @warn "Resetting vsm solver."
            end

            n_unrefined =
                wing.vsm_wing.n_unrefined_sections
            group_idxs = wing.group_idxs

            moment_frac = if isempty(group_idxs)
                0.25
            elseif length(groups) >= maximum(group_idxs)
                groups[first(group_idxs)].moment_frac
            else
                0.25
            end

            theta_idxs = isempty(group_idxs) ?
                nothing : (4:(3 + n_unrefined))

            # Both modes call linearize to get the VSM
            # solution at the current operating point.
            # AERO_LINEARIZED also stores the Jacobian.
            res = VortexStepMethod.linearize(
                wing.vsm_solver,
                wing.vsm_aero,
                wing.vsm_y;
                va_idxs=1:3,
                theta_idxs=theta_idxs,
                omega_idxs=(4 + n_unrefined):(6 + n_unrefined),
                moment_frac=moment_frac,
                aero_coeffs=true
            )

            if wing.aero_mode == AERO_LINEARIZED
                # Store Jacobian and coefficients for
                # the symbolic linearization equations
                wing.vsm_jac .= res[1]
                wing.vsm_x .= res[2]

            elseif wing.aero_mode == AERO_DIRECT
                # Compute physical forces from baseline
                # coefficients: F = q∞·A·C₀
                va_sq = wing.va_b[1]^2 +
                    wing.va_b[2]^2 + wing.va_b[3]^2
                rho = calc_rho(
                    sam.am, wing.pos_w[3])
                q_inf = 0.5 * rho * va_sq
                area = wing.vsm_aero.projected_area
                coeffs = res[2]
                wing.aero_force_b .=
                    q_inf * area .* coeffs[1:3]
                wing.aero_moment_b .=
                    q_inf * area .* coeffs[4:6]
            end

            # Map unrefined moments back to groups
            # (same for both LINEARIZED and DIRECT)
            if !isempty(group_idxs)
                unrefined_moments = res[2][7:end]
                for gidx in group_idxs
                    g = groups[gidx]
                    g.aero_moment = sum(
                        unrefined_moments[
                            g.unrefined_section_idxs])
                end
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

            VortexStepMethod.solve!(
                wing.vsm_solver, wing.vsm_aero;
                log=false)
            distribute_panel_forces_to_points!(
                wing, points)
        end
    end

    nothing
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

