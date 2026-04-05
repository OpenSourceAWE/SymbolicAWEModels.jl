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
    @unpack winches, wings = sam.sys_struct
    old_brakes = [winch.brake for winch in winches]
    old_fixes = [wing.fix_sphere for wing in wings]
    [winch.brake=true for winch in winches]
    [wing.fix_sphere=true for wing in wings]
    for _ in 1:Int(round(t÷dt))
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
                     integ=sam.integrator)
    wings = sam.sys_struct.wings
    groups = sam.sys_struct.groups
    points = sam.sys_struct.points

    length(wings) == 0 && return nothing

    # Handle QUATERNION wings
    has_quaternion_wings = any(
        w.wing_type == QUATERNION for w in wings)
    if has_quaternion_wings && !isnothing(prob.get_vsm_y)
        vsm_y = prob.get_vsm_y(integ)

        for wing in wings
            wing.wing_type != QUATERNION && continue
            wing.aero_mode == AERO_NONE && continue

            wing.vsm_y .= vsm_y[:, wing.idx]
            if norm(wing.vsm_y[1:3]) < 1e-3
                @warn "Apparent wind too small for VSM " *
                    "linearization (va_b=$(wing.vsm_y[1:3]))" *
                    ", skipping wing $(wing.idx)"
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
        w.wing_type == REFINE for w in wings)
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

                set_va!(wing.vsm_aero, va_dist,
                    zeros(MVec3))
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

    # copy state and settings to lin prob
    lin_prob.set_sys(lin_prob.prob, sam.sys_struct)
    lin_prob.set_set(lin_prob.prob, sam.set)

    lin_model = solve(lin_prob.prob)[1]
    return lin_model
end

"""
    getstate(sys_struct::SystemStructure) -> Tuple

Capture and return a snapshot of the key dynamic states of the system.
"""
function getstate(sys_struct::SystemStructure)
    @unpack wings, winches = sys_struct
    c = copy
    wing = wings[1]
    tether_len = [winch.tether_len for winch in winches]
    tether_vel = [winch.tether_vel for winch in winches]
    return (c(wing.pos_w), c(wing.vel_w),
            c(wing.wind_disturb), c(wing.R_b_to_w),
            c(wing.ω_b), tether_len, tether_vel,
            c(wing.com_w), c(wing.com_vel),
            c(wing.Q_p_to_w), c(wing.ω_p))
end

"""
    setstate!(sys_struct::SystemStructure, state)

Set the key dynamic states of the system from a snapshot tuple.
"""
function setstate!(sys_struct::SystemStructure, state)
    @unpack wings, winches = sys_struct
    wing = wings[1]
    (pos_w, vel_w, wind_disturb, R_b_to_w, ω_b,
     tether_len, tether_vel,
     com_w, com_vel, Q_p_to_w, ω_p) = state
    wing.pos_w .= pos_w
    wing.vel_w .= vel_w
    wing.wind_disturb .= wind_disturb
    wing.R_b_to_w .= R_b_to_w
    wing.ω_b .= ω_b
    wing.com_w .= com_w
    wing.com_vel .= com_vel
    wing.Q_p_to_w .= Q_p_to_w
    wing.ω_p .= ω_p
    for winch in winches
        winch.tether_len = tether_len[winch.idx]
        winch.tether_vel = tether_vel[winch.idx]
    end
end

"""
    set_measured!(sys_struct, heading, turn_rate, tether_len, tether_vel)

Adjust the model's state to match a set of "measured" or target values.

This function is typically used in state estimation or stabilization loops. It takes
a set of target values (e.g., from a sensor or a reference trajectory) and forces
the corresponding states in the `SystemStructure` to match, calculating kinematically
consistent values for other related states.
"""
function set_measured!(sys_struct::SystemStructure, 
    heading, turn_rate,
    tether_len, tether_vel
)
    @unpack wings, winches = sys_struct
    wing = wings[1]

    # get variables from integrator
    distance = norm(wing.pos_w)
    R_t_to_w = calc_R_t_to_w(wing.pos_w)

    # Update wing position for new tether length
    wing.pos_w .= R_t_to_w * [0, 0,
        distance + tether_len[1] - winches[1].tether_len]
    wing.vel_w .= R_t_to_w *
        [-wing.elevation_vel, wing.azimuth_vel, 0.0]
    wing.wind_disturb .= R_t_to_w *
        [0.0, 0.0, -tether_vel[1]]
    # Rotate to target heading (tangential sphere frame)
    cur_heading = calc_heading(
        wing.R_b_to_w, wing.pos_w)
    d_heading = heading - cur_heading
    R_b_to_w = zeros(3, 3)
    for i in 1:3
        R_b_to_w[:, i] .= R_t_to_w * rotate_around_z(
            R_t_to_w' * wing.R_b_to_w[:, i],
            d_heading)
    end
    wing.R_b_to_w = R_b_to_w
    # adjust the turn rates for observed turn rate
    wing.ω_b .= wing.R_b_to_w' * R_t_to_w * [wing.turn_rate[1], wing.turn_rate[2], turn_rate]
    # directly set tether length
    for winch in winches
        winch.tether_len = tether_len[winch.idx]
        winch.tether_vel = tether_vel[winch.idx]
    end
    # Derive principal frame from body frame
    R_b_to_w_mat = wing.R_b_to_w
    wing.com_w .= wing.pos_w .+
        R_b_to_w_mat * wing.com_offset_b
    R_p_to_w = R_b_to_w_mat * wing.R_b_to_c' * wing.R_p_to_c
    wing.Q_p_to_w .= rotation_matrix_to_quaternion(R_p_to_w)
    ω_w = R_b_to_w_mat * wing.ω_b
    r_com_w = R_b_to_w_mat * wing.com_offset_b
    wing.com_vel .= wing.vel_w .+ cross(ω_w, r_com_w)
    wing.ω_p .= R_p_to_w' * ω_w
    return nothing
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

function jacobian(f::Function, x::AbstractVector,
                  ϵ::AbstractVector, arg2)
    n = length(x)
    fx = f(x, arg2)
    m = length(fx)
    J = zeros(m, n)
    for i in 1:n
        x_perturbed = copy(x)
        x_perturbed[i] += ϵ[i]
        J[:, i] = (f(x_perturbed, arg2) - fx) / ϵ[i]
    end
    return J
end

function jacobian(f::Function, ::Val{:second_arg},
                  x::AbstractVector, ϵ::AbstractVector, arg1)
    n = length(x)
    fx = f(arg1, x)
    m = length(fx)
    J = zeros(m, n)
    for i in 1:n
        x_perturbed = copy(x)
        x_perturbed[i] += ϵ[i]
        J[:, i] = (f(arg1, x_perturbed) - fx) / ϵ[i]
    end
    return J
end

"""
    simple_linearize!(s::SymbolicAWEModel; tstab=10.0) -> LinType

Compute a simplified, low-order state-space model by numerically linearizing the full simulation.

This function performs system identification by perturbing the states and inputs of the
full nonlinear model and observing the effect on the state derivatives and outputs. It
runs the simulation for a short duration (`tstab`) after each perturbation to find the
steady-state response.

# Arguments
- `s::SymbolicAWEModel`: The model to linearize.

# Keywords
- `tstab::Float64=10.0`: The simulation time [s] to run after each perturbation to reach a steady state.

# Returns
- `LinType`: A NamedTuple `(A, B, C, D)` containing the identified low-order state-space matrices.
"""
function simple_linearize!(sam::SymbolicAWEModel; tstab=10.0)
    @unpack segments, winches, tethers, wings = sam.sys_struct
    integ = sam.integrator::OrdinaryDiffEqCore.ODEIntegrator
    prob = sam.prob
    update_sys_struct!(sam.prob, sam.integrator, sam.sys_struct)
    state0 = getstate(sam.sys_struct)
    old_brakes = [winch.brake for winch in winches]
    old_fixes = [wing.fix_sphere for wing in wings]
    [winch.brake=true for winch in winches]
    [wing.fix_sphere=true for wing in wings]
    lin_x0 = sam.simple_lin_model.get_x(integ)
    u0 = [winch.set_value for winch in sam.sys_struct.winches]
    @unpack A, B, C, D = sam.simple_lin_model.model
    A .= 0.0
    B .= 0.0
    C .= 0.0
    D .= 0.0

    # TODO: add sparsity pattern for the known zeros
    function f(x, u)
        heading = x[1]
        turn_rate = x[2]
        tether_len = x[3:5]
        tether_vel = x[6:8]
        set_measured!(sam.sys_struct, heading, turn_rate,
                        tether_len, tether_vel)
        prob.set_set_values(integ, u)
        OrdinaryDiffEqCore.reinit!(integ)
        ModelingToolkit.SciMLBase.step!(integ, tstab, true)
        return sam.simple_lin_model.get_dx(integ)
    end

    # yes it looks weird to step in an output function, but this is a steady state finder rather than output
    function h(x, u)
        heading = x[1]
        turn_rate = x[2]
        tether_len = x[3:5]
        tether_vel = x[6:8]
        set_measured!(sam.sys_struct, heading, turn_rate,
                        tether_len, tether_vel)
        prob.set_set_values(integ, u)
        OrdinaryDiffEqCore.reinit!(integ)
        ModelingToolkit.SciMLBase.step!(integ, tstab, true)
        return sam.simple_lin_model.get_y(integ)
    end

    segment = segments[tethers[1].segment_idxs[1]]
    mass_per_meter = sam.set.rho_tether * π * (segment.diameter/2)^2
    mass = winches[1].tether_len * mass_per_meter + sam.set.mass

    # calculate jacobian
    ϵ_x = [0.001, 0.1, 0.001, 0.001, 0.001, 0.1, 0.1, 0.1]
    ϵ_u = [1.0, 0.1, 0.1]
    A .= jacobian(f, lin_x0, ϵ_x, u0)
    B .= jacobian(f, Val(:second_arg), u0, ϵ_u, lin_x0)
    C .= jacobian(h, lin_x0, ϵ_x, u0)
    D .= 0.0
    D[4,1] = -mass * B[6,1]
    A[:,1] .= 0.0 # Aero moment due to change in heading cannot be found in steady state
    C[4,1] = 0.0
    prob.set_set_values(integ, u0)
    [winch.brake=old_brakes[winch.idx] for winch in winches]
    [wing.fix_sphere=old_fixes[wing.idx] for wing in wings]
    setstate!(sam.sys_struct, state0)
    OrdinaryDiffEqCore.reinit!(integ)
    return sam.simple_lin_model.model
end
