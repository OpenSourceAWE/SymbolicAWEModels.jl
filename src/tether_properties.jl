# SPDX-FileCopyrightText: 2025 Bart van de Lint <bart@vandelint.net>
#
# SPDX-License-Identifier: LGPL-3.0-only

"""
    in_percent_band(x, steady, delta_x, i, p) -> Bool

Helper function to check if a time series has settled within a percentage band.

It checks if all values of the time series `x` from index `i` to the end are
within a tolerance band defined by `p` percent of the total change `delta_x`.
"""
function in_percent_band(x, steady, delta_x, i, p)
    tol = p/100 * abs(delta_x)
    # All subsequent points must be within steady ± p%
    all(abs.(x[i:end] .- steady) .<= tol)
end

"""
    calc_spring_props(sam, tether_sam; prn=false) -> (Vector, Vector, Matrix, Float64)

Calculate the equivalent stiffness and damping, and return the step response data.

This function orchestrates the process by performing a step response test on the
`tether_sam` model and then analyzing the resulting tether length data.

# Arguments
- `sam::SymbolicAWEModel`: The reference model, used for its physical properties.
- `tether_sam::SymbolicAWEModel`: A copy of the model to perform the step test on.

# Keywords
- `prn::Bool=false`: If true, enables printing of intermediate results.

# Returns
- `Tuple{Vector{Float64}, Vector{Float64}, Matrix{Float64}, Float64}`: A tuple containing:
    1.  `unit_stiffness` [N]
    2.  `unit_damping` [Ns]
    3.  `tether_lens` (the step response data)
    4.  `dt` (the simulation time step)
"""
function calc_spring_props(sam::SymbolicAWEModel, tether_sam::SymbolicAWEModel; 
                           F_step=-0.1, prn=false)
    find_steady_state!(sam; t=10.0, dt=10.0, vsm_interval=0)
    copy!(sam.sys_struct, tether_sam.sys_struct)
    integrator = tether_sam.integrator
    prob = tether_sam.prob
    isnothing(integrator) && error("tether_sam.integrator is not initialized")
    isnothing(prob) && error("tether_sam.prob is not initialized")
    OrdinaryDiffEqCore.reinit!(integrator; reinit_dae=true)
    update_sys_struct!(prob, integrator, tether_sam.sys_struct)

    F_0 = [-tether_sam.sys_struct.points[i].force for i in 1:4]
    steps = 200
    tether_lens = step(tether_sam, steps, F_step, F_0)
    k_values, c_values = calc_spring_props(sam, tether_lens, F_step; prn)
    
    dt = 1/sam.set.sample_freq
    return k_values .* tether_lens[:,1], c_values .* tether_lens[:,1], tether_lens, dt
end

"""
    calc_spring_props(sam, tether_lens, F_step; p=5, prn=false) -> (Vector, Vector)

Calculate spring constant `k` and damping coefficient `c` from a step response.

This function analyzes the time series of tether lengths (`tether_lens`) resulting
from a step force (`F_step`) to estimate the parameters of an equivalent second-order
mass-spring-damper system.

# Arguments
- `sam::SymbolicAWEModel`: The model from which to take physical parameters (mass).
- `tether_lens::Matrix{Float64}`: A matrix of tether length time series data.
- `F_step::Float64`: The magnitude of the applied step force.

# Keywords
- `p::Int=5`: The percentage band used to determine the settling time.
- `prn::Bool=false`: If true, enables printing of detailed calculations.

# Returns
- `Tuple{Vector{Float64}, Vector{Float64}}`: A tuple containing two vectors:
    1.  `k_values` (spring constants [N/m])
    2.  `c_values` (damping coefficients [Ns/m])
"""
function calc_spring_props(sam::SymbolicAWEModel, tether_lens, F_step; p=5, prn=false)
    (; tethers, segments) = sam.sys_struct
    set = sam.set
    dt = 1/set.sample_freq

    k_values = zeros(4)
    c_values = zeros(4)

    first_segments = [segments[tether.segment_idxs[1]] for tether in tethers]
    mass_per_meter = [seg.density * π * (seg.diameter / 2)^2
        for seg in first_segments]

    for j in eachindex(tethers)
        tether_len_series = tether_lens[j, :]
        initial_len = tether_len_series[1]
        final_len = tether_len_series[end]
        delta_x_ss = final_len - initial_len

        m = mass_per_meter[j] * 0.5 * tethers[j].stretched_len
        @assert m > 0

        if abs(delta_x_ss) < 1e-6
            @warn "Steady-state change too small for Tether $j; skipping."
            k_values[j] = NaN; c_values[j] = NaN
            continue
        end

        # Spring stiffness
        k = F_step / delta_x_ss
        k_values[j] = k
        ω_n = sqrt(k/m)

        # Find settling time index according to your in_percent_band function:
        T_s_index = -1
        for i in 1:length(tether_len_series)
            if in_percent_band(tether_len_series, final_len, delta_x_ss, i, p)
                T_s_index = i
                break
            end
        end

        if T_s_index == -1
            @warn "Could not find settling time ($p% criterion) for Tether $j; using fallback."
            # Fallback: use time constant method to find tau
            target_len = initial_len + (1 - 1/ℯ) * delta_x_ss
            tau_idx = findfirst(i ->
                (delta_x_ss > 0 && tether_len_series[i] >= target_len) ||
                (delta_x_ss < 0 && tether_len_series[i] <= target_len),
                1:length(tether_len_series))
            if tau_idx === nothing
                @warn "Cannot determine tau for Tether $j"
                c_values[j] = NaN
            else
                tau = (tau_idx-1)*dt
                c = k * tau
                c_values[j] = c
                println("Tether $j fallback c=", c)
            end
            continue
        end

        T_s = (T_s_index - 1) * dt
        # Calculate damping ratio based on variable percentage settling criterion:
        X = -log(p / 100)
        ζ = X / (ω_n * T_s)
        c = 2 * ζ * sqrt(k * m)
        c_values[j] = c
    end

    prn && for j in eachindex(tethers)
        println("Tether $(j): k = $(k_values[j]) N/m, c = $(c_values[j]) Ns/m")
    end
    return k_values, c_values
end

"""
    step(sam, steps, F_step, F_0; abs_tol, consecutive_steps_needed, prn) -> Matrix

Apply a step force to a model and simulate its dynamic response.

This function records the length of each tether over a specified number of simulation
steps. It includes an early exit condition if the system's state settles.

# Arguments
- `sam::SymbolicAWEModel`: The model to be simulated.
- `steps::Int`: The total number of simulation steps.
- `F_step::Float64`: The magnitude of the step force to apply.
- `F_0::Vector{KVec3}`: The initial force vector for each tether attachment point.

# Keywords
- `abs_tol::Float64=1e-6`: Absolute tolerance for the settling check.
- `consecutive_steps_needed::Int=10`: Number of consecutive steps required to be
  within tolerance to be considered settled.
- `prn::Bool=false`: If true, enables printing of status messages.

# Returns
- `Matrix{Float64}`: A matrix where each row corresponds to a tether and each
  column to a time step, containing the tether lengths.
"""
function step(sam::SymbolicAWEModel, steps, F_step, F_0;
              abs_tol=1e-6,
              consecutive_steps_needed=10,
              prn=false)

    (; points, tethers) = sam.sys_struct

    initial_tether_lens = [norm(points[i].pos_w) for i in eachindex(tethers)]
    [points[i].disturb .= F_0[i] .+ F_step * normalize(points[i].pos_w) for i in eachindex(tethers)]

    tether_lens = zeros(length(tethers), steps+1)
    tether_lens[:, 1] .= initial_tether_lens # Store the initial lengths
    settled_steps = 0
    for step in 1:steps
        next_step!(sam; vsm_interval=0)
        for j in eachindex(tethers)
            tether_lens[j, step+1] = norm(points[j].pos_w)
        end
        # Check absolute delta for all tethers
        step_deltas = abs.(tether_lens[:, step+1] .- tether_lens[:, step])
        max_delta = maximum(step_deltas)
        if max_delta < abs_tol
            settled_steps += 1
        else
            settled_steps = 0
        end
        if settled_steps >= consecutive_steps_needed
            prn && println("Stopped at step $step: all tethers within $abs_tol for $consecutive_steps_needed steps.")
            tether_lens[:, step+2:end] .= tether_lens[:, step+1]
            break
        end
    end
    if settled_steps < consecutive_steps_needed
        @warn "Stepping simulation did not settle within the given steps."
    end
    return tether_lens
end

"""
    update_segment_forces!(sys_struct::SystemStructure)

Calculate and update spring forces for all segments in-place.

This function computes the spring-damper forces for each segment using
the same formulas as in `generate_system.jl`. It updates the `len` and
`force` fields of each segment in the SystemStructure.

The spring-damper force follows Hooke's law with damping:

```math
F = k(l - l_0) - c\\dot{l}
```

where:
- `k = unit_stiffness / l` (tension) or
  `k = compression_frac * unit_stiffness / l` (compression)
- `l` is current length, `l_0` is unstretched length
- `c = unit_damping / l` is damping coefficient
- `\\dot{l} = (v₁ - v₂) ⋅ û` is extension rate

# Arguments
- `sys_struct::SystemStructure`: The system structure containing points
  and segments.

# Returns
- `nothing`: Modifies segment.len and segment.force in-place.

# Example
```julia
update_segment_forces!(sam.sys_struct)
for segment in sam.sys_struct.segments
    println("Segment \$(segment.idx): force = \$(segment.force) N")
end
```
"""
function update_segment_forces!(sys_struct::SystemStructure)
    (; points, segments) = sys_struct

    for segment in segments
        p1, p2 = segment.point_idxs

        # Segment vector and length
        segment_vec = points[p2].pos_w - points[p1].pos_w
        len = norm(segment_vec)
        unit_vec = segment_vec / len

        # Relative velocity along segment
        rel_vel = points[p1].vel_w - points[p2].vel_w
        spring_vel = rel_vel ⋅ unit_vec

        # Stiffness (handles compression)
        if len > segment.l0
            stiffness = segment.unit_stiffness / len
        else
            stiffness = segment.compression_frac *
                        segment.unit_stiffness / len
        end

        # Damping
        damping = segment.unit_damping / len

        # Update segment fields in-place
        segment.len = len
        segment.force = stiffness * (len - segment.l0) -
                        damping * spring_vel
    end

    return nothing
end
