# Copyright (c) 2025 Bart van de Lint, Jelle Poland
# SPDX-License-Identifier: MPL-2.0

"""
Example: 2D trajectory with cosine steering input where steering = turn rate

This plots the shape traced when:
- Turn rate ψ̇ = A·cos(ω·t)  (positive = turn right)
- Constant forward speed v
- Heading ψ(t) = ψ₀ - (A/ω)·sin(ω·t)

The trajectory closes when A/ω equals a zero of J₀ (Bessel function).

Run from REPL:
    include("examples/cosine_steering_trajectory.jl")

    # Closed trajectory using 1st zero of J₀
    plot_trajectory(j0_zero=1)

    # Slightly open trajectory
    plot_trajectory(j0_zero=1, δ=0.01)

    # More complex closed pattern (2nd zero)
    plot_trajectory(j0_zero=2)

    # Animated version - watch the simulation evolve
    animate_trajectory(j0_zero=1)
    animate_trajectory(j0_zero=1, δ=0.05, speed=0.5)  # slower
    animate_trajectory(j0_zero=2, n_periods=2)        # more loops
"""

using GLMakie
using Interpolations
using NLsolve
using Statistics

# First 10 zeros of J₀ (precomputed for efficiency)
const J0_ZEROS = [2.4048255576957728, 5.5200781102863106, 8.653727912911012,
                  11.791534439014281, 14.930917708487785, 18.071063967910922,
                  21.211636629879258, 24.352471530749302, 27.493479132040254,
                  30.634606468431975]

"""
    get_closed_amplitude(; ω=1.0, j0_zero=1, δ=0.0)

Calculate the turn rate amplitude A that gives a closed trajectory.

- `j0_zero`: which zero of J₀ to use (1, 2, 3, ...)
- `δ`: relative perturbation (e.g., 0.01 shifts A by 1%)
"""
function get_closed_amplitude(; ω = 1.0, j0_zero = 1, δ = 0.0)
    B = J0_ZEROS[j0_zero]  # A/ω for closed trajectory
    B_perturbed = B * (1 + δ)
    return B_perturbed * ω
end

"""
    compute_steering_vs_x(; j0_zero=1, n_points=1000)

Compute the steering(x) relationship for a closed trajectory.
Returns x values and corresponding steering values.
"""
function compute_steering_vs_x(; j0_zero = 1, n_points = 1000)
    B = J0_ZEROS[j0_zero]  # A/ω for closed trajectory
    ψ₀ = π/2

    # Parametric: x(θ) and steering(θ) where θ = ωt ∈ [0, π] (half cycle)
    θ = range(0, π, length = n_points)

    # Heading: ψ = ψ₀ - B·sin(θ)
    ψ = ψ₀ .- B .* sin.(θ)

    # Steering (normalized): A·cos(θ)/A = cos(θ)
    steering_normalized = cos.(θ)

    # x position and arc length distance
    dθ = θ[2] - θ[1]
    x = zeros(n_points)
    dist = zeros(n_points)
    for i in 2:n_points
        dx = cos(ψ[i-1]) * dθ
        dy = sin(ψ[i-1]) * dθ
        x[i] = x[i-1] + dx
        dist[i] = dist[i-1] + sqrt(dx^2 + dy^2)
    end

    # Normalize distance to cycle fraction (0 to 1)
    cycle_length = dist[end]
    cycle_frac = dist ./ cycle_length

    return (; x, cycle_frac, steering_normalized, θ, ψ, B, cycle_length)
end

"""
    create_steering_interp(; j0_zero=1)

Create interpolation functions for steering(x) and steering(dist).
Returns interpolators and data.
"""
function create_steering_interp(; j0_zero = 1, n_points = 1000)
    data = compute_steering_vs_x(; j0_zero, n_points)

    # Linear interpolation with flat extrapolation using Interpolations.jl
    itp_x = linear_interpolation(data.x, data.steering_normalized, extrapolation_bc=Flat())
    # cycle fraction based (0 to 1)
    itp_frac = linear_interpolation(data.cycle_frac, data.steering_normalized, extrapolation_bc=Flat())

    return (; itp_x, itp_frac, data)
end

"""
    simulate_from_steering_vs_x(; ...)

Simulate trajectory using blended steering from x-position and distance lookups.
- x_weight=1.0: pure x-based steering
- x_weight=0.0: pure distance-based steering
- x_weight=0.5: 50/50 blend (default)
"""
function simulate_from_steering_vs_x(;
    v = 1.0,           # forward speed
    ω = 1.0,           # angular frequency (determines A from A/ω)
    j0_zero = 1,       # which zero of J₀ to use
    δ = 0.0,           # relative perturbation of steering amplitude
    ψ₀ = π/2,          # initial heading (rad)
    T = 2π/ω,          # simulation time
    dt = 0.01,         # time step
    x_weight = 0.5,    # blend: 1.0 = pure x, 0.0 = pure distance, 0.5 = 50/50
    v_mult = 1.0,      # velocity multiplier (applied to position, not distance)
    v_offset = 0.0,    # velocity offset (applied to position, not distance)
    ψ_dot_mult = 1.0,  # turn rate multiplier (applied to heading, not stored in ψ_dot)
    ψ_dot_offset = 0.0 # turn rate offset (applied to heading, not stored in ψ_dot)
)
    # Get steering interpolations
    interp = create_steering_interp(; j0_zero)
    A = J0_ZEROS[j0_zero] * ω * (1 + δ)  # steering amplitude
    cycle_length = interp.data.cycle_length

    t = 0:dt:T
    n = length(t)

    x = zeros(n)
    y = zeros(n)
    ψ = zeros(n)
    ψ_dot = zeros(n)
    cycle_dist = zeros(n)  # distance within current half-cycle (resets on direction change)
    cycle_frac = zeros(n)  # cycle fraction (0 to 1)

    ψ[1] = ψ₀
    prev_vx = cos(ψ₀)  # initial x-velocity direction

    moving_right = true  # start moving right (vx > 0 initially since ψ₀ = π/2)

    for i in 1:n
        # Compute cycle fraction (0 to 1)
        raw_frac = clamp(cycle_dist[i] / cycle_length, 0.0, 1.0)
        # When moving right: 0→1, when moving left: 1→0
        cycle_frac[i] = moving_right ? raw_frac : (1.0 - raw_frac)

        # Blend steering from x and cycle fraction lookups
        steering_from_x = interp.itp_x(x[i])
        steering_from_frac = interp.itp_frac(cycle_frac[i])
        steering_normalized = x_weight * steering_from_x + (1 - x_weight) * steering_from_frac

        ψ_dot[i] = A * steering_normalized * (1 + δ)

        if i < n
            # Update heading: dψ/dt = -steering (positive steering = turn right = decrease ψ)
            # Apply turn rate disturbance (mult/offset) but ψ_dot stores clean control law
            ψ_dot_actual = ψ_dot[i] * ψ_dot_mult + ψ_dot_offset
            ψ[i + 1] = ψ[i] - ψ_dot_actual * dt

            # Compute clean position update (for distance tracking)
            vx = cos(ψ[i])
            vy = sin(ψ[i])
            dx_clean = v * vx * dt
            dy_clean = v * vy * dt

            # Apply velocity disturbance (mult/offset) to position (but NOT to distance calculation)
            v_actual = v * v_mult + v_offset
            dx = v_actual * vx * dt
            dy = v_actual * vy * dt
            x[i + 1] = x[i] + dx
            y[i + 1] = y[i] + dy

            # Detect direction changes and reset cycle (use clean distance)
            if prev_vx >= 0 && vx < 0
                # Transition right→left: start counting down from 1
                cycle_dist[i + 1] = 0.0
                moving_right = false
            elseif prev_vx < 0 && vx >= 0
                # Transition left→right: start counting up from 0
                cycle_dist[i + 1] = 0.0
                moving_right = true
            else
                cycle_dist[i + 1] = cycle_dist[i] + sqrt(dx_clean^2 + dy_clean^2)
            end
            prev_vx = vx
        end
    end

    return (; t, ψ_dot, ψ, x, y, cycle_frac, steering_data=interp.data)
end

"""
    simulate_cosine_steering(; ...)

Time-based simulation with cosine steering.

Control parameters for trajectory adjustment:
- `x_shift`: Sine term at 2ω → shifts x (integrates to zero, no tilt)
- `y_shift`: Amplitude perturbation → shifts y extent
- `tilt`: Constant turn rate offset → causes progressive tilt
- `φ`: Phase offset in steering cosine (rad)
"""
function simulate_cosine_steering(;
    v = 1.0,           # forward speed
    ω = 1.0,           # angular frequency of steering input (rad/s)
    j0_zero = 1,       # which zero of J₀ to use (1, 2, 3, ...)
    x_shift = 0.0,     # sine term at 2ω → x shift
    y_shift = 0.0,     # amplitude perturbation → y shift
    tilt = 0.0,        # constant turn rate offset → progressive tilt
    φ = 0.0,           # phase offset in steering (rad)
    ψ₀ = π/2,          # initial heading (rad), π/2 = +y direction
    T = 2π/ω,          # simulation time (one full period by default)
    dt = 0.01          # time step
)
    A = get_closed_amplitude(; ω, j0_zero, δ=y_shift)
    t = 0:dt:T
    n = length(t)

    # Turn rate = cosine + sine(2ω) + offset
    # - cos term: base figure-8
    # - sin(2ω) term: x shift (integrates to zero)
    # - tilt: progressive rotation
    ψ_dot = A .* cos.(ω .* t .+ φ) .+ x_shift .* sin.(ω .* t * 2) .+ tilt

    # Heading from integrating turn rate numerically (analytical is complex with modulation)
    ψ = zeros(n)
    ψ[1] = ψ₀
    for i in 2:n
        ψ[i] = ψ[i-1] - ψ_dot[i-1] * dt
    end

    # Integrate position using heading
    x = zeros(n)
    y = zeros(n)
    for i in 2:n
        x[i] = x[i-1] + v * cos(ψ[i-1]) * dt
        y[i] = y[i-1] + v * sin(ψ[i-1]) * dt
    end

    return (; t, ψ_dot, ψ, x, y)
end

"""
    compute_trajectory_metrics(x, y)

Compute centroid and tilt of a trajectory.
Returns (centroid_x, centroid_y, tilt) where tilt is in radians.
Tilt is measured as the angle of the principal axis (via PCA).
"""
function compute_trajectory_metrics(x, y)
    # Centroid
    cx = mean(x)
    cy = mean(y)

    # Tilt via PCA: find principal axis direction
    xc = x .- cx
    yc = y .- cy

    # Covariance matrix
    cov_xx = mean(xc .* xc)
    cov_yy = mean(yc .* yc)
    cov_xy = mean(xc .* yc)

    # Principal axis angle
    tilt = 0.5 * atan(2 * cov_xy, cov_xx - cov_yy)

    return (; centroid_x=cx, centroid_y=cy, tilt)
end

"""
    solve_for_target(target_x, target_y, target_tilt; kwargs...)

Solve for (x_shift, y_shift, tilt) parameters that place the trajectory
centroid at (target_x, target_y) with the given tilt.

Returns named tuple with x_shift, y_shift, tilt and convergence info.
"""
function solve_for_target(target_x, target_y, target_tilt;
    v = 1.0, ω = 1.0, j0_zero = 1, φ = 0.0, dt = 0.01
)
    T = 2π / ω  # one full period

    function residual!(F, params)
        xs, ys, tp = params
        result = simulate_cosine_steering(; v, ω, j0_zero, x_shift=xs, y_shift=ys, tilt=tp, φ, T, dt)
        metrics = compute_trajectory_metrics(result.x, result.y)

        F[1] = metrics.centroid_x - target_x
        F[2] = metrics.centroid_y - target_y
        # Normalize tilt difference to [-π, π]
        tilt_diff = metrics.tilt - target_tilt
        F[3] = atan(sin(tilt_diff), cos(tilt_diff))
    end

    # Initial guess
    initial_guess = [0.0, 0.0, 0.0]

    sol = nlsolve(residual!, initial_guess; ftol=1e-10, iterations=100)

    x_shift, y_shift, tilt_param = sol.zero
    return (; x_shift, y_shift, tilt=tilt_param, converged=converged(sol), sol)
end

"""
    plot_steering_law(; j0_zero=1)

Plot the steering(x) control law - the key relationship for position-based control.
"""
function plot_steering_law(; j0_zero = 1)
    data = compute_steering_vs_x(; j0_zero)

    fig = Figure(size = (800, 400))
    ax = Axis(fig[1, 1],
        xlabel = "x position",
        ylabel = "Steering (normalized)",
        title = "Steering vs x Control Law (J₀ zero #$j0_zero, B = $(round(data.B, digits=3)))")

    lines!(ax, data.x, data.steering_normalized, linewidth = 2, color = :coral)
    scatter!(ax, [data.x[1]], [data.steering_normalized[1]], color = :green, markersize = 15, label = "Start")

    display(fig)
    fig
end

"""
    animate_x_based(; ...)

Animate using x-position-based steering control.

Disturbance parameters (applied to simulation, not affecting control law or distance calc):
- `v_mult`: Velocity multiplier (default 1.0)
- `v_offset`: Velocity offset (default 0.0)
- `ψ_dot_mult`: Turn rate multiplier (default 1.0)
- `ψ_dot_offset`: Turn rate offset (default 0.0)
"""
function animate_x_based(; ω = 1.0, j0_zero = 1, δ = 0.0, n_periods = 1, fps = 30,
                         speed = 1.0, dt = 0.01, x_weight = 0.5,
                         v_mult = 1.0, v_offset = 0.0,
                         ψ_dot_mult = 1.0, ψ_dot_offset = 0.0)
    T = n_periods * 2π / ω
    result = simulate_from_steering_vs_x(; ω, j0_zero, δ, T, dt, x_weight,
                                          v_mult, v_offset, ψ_dot_mult, ψ_dot_offset)
    A = J0_ZEROS[j0_zero] * ω * (1 + δ)
    n = length(result.t)

    # Observables
    idx = Observable(1)
    x_trace = @lift result.x[1:$idx]
    y_trace = @lift result.y[1:$idx]
    cycle_frac_trace = @lift result.cycle_frac[1:$idx]
    ψ_dot_trace = @lift result.ψ_dot[1:$idx]
    x_cur = @lift [result.x[$idx]]
    y_cur = @lift [result.y[$idx]]
    cycle_frac_cur = @lift [result.cycle_frac[$idx]]
    ψ_dot_cur = @lift [result.ψ_dot[$idx]]

    fig = Figure(size = (1400, 500))

    # Margins
    margin = 0.1
    x_range = maximum(result.x) - minimum(result.x)
    y_range = maximum(result.y) - minimum(result.y)
    ψ_dot_margin = 0.1 * (maximum(result.ψ_dot) - minimum(result.ψ_dot))

    # Steering vs x
    ax1 = Axis(fig[1, 1],
        xlabel = "x", ylabel = "Steering (rad/s)",
        title = "Steering vs x")
    xlims!(ax1, minimum(result.x) - margin * x_range, maximum(result.x) + margin * x_range)
    ylims!(ax1, minimum(result.ψ_dot) - ψ_dot_margin, maximum(result.ψ_dot) + ψ_dot_margin)
    # Show the control law curve
    lines!(ax1, result.steering_data.x, result.steering_data.steering_normalized .* A,
           linewidth = 1, color = (:coral, 0.3))
    # Animated trace
    lines!(ax1, x_trace, ψ_dot_trace, linewidth = 2, color = :coral)
    scatter!(ax1, [result.x[1]], [result.ψ_dot[1]], color = :green, markersize = 12)
    scatter!(ax1, x_cur, ψ_dot_cur, color = :red, markersize = 12)

    # Steering vs cycle fraction
    ax2 = Axis(fig[1, 2],
        xlabel = "cycle fraction", ylabel = "Steering (rad/s)",
        title = "Steering vs cycle fraction")
    xlims!(ax2, 0, 1.1)
    ylims!(ax2, minimum(result.ψ_dot) - ψ_dot_margin, maximum(result.ψ_dot) + ψ_dot_margin)
    # Show the control law curve
    lines!(ax2, result.steering_data.cycle_frac, result.steering_data.steering_normalized .* A,
           linewidth = 1, color = (:teal, 0.3))
    # Animated trace
    lines!(ax2, cycle_frac_trace, ψ_dot_trace, linewidth = 2, color = :teal)
    scatter!(ax2, [0.0], [result.ψ_dot[1]], color = :green, markersize = 12)
    scatter!(ax2, cycle_frac_cur, ψ_dot_cur, color = :red, markersize = 12)

    # 2D trajectory
    has_disturb = (v_mult != 1.0 || v_offset != 0.0 || ψ_dot_mult != 1.0 || ψ_dot_offset != 0.0)
    disturb_str = has_disturb ? ", v×$(round(v_mult,digits=2))+$(round(v_offset,digits=2)), ψ̇×$(round(ψ_dot_mult,digits=2))+$(round(ψ_dot_offset,digits=2))" : ""
    ax3 = Axis(fig[1, 3],
        xlabel = "x", ylabel = "y",
        title = "2D Trajectory (x_weight = $x_weight$disturb_str)",
        aspect = DataAspect())
    xlims!(ax3, minimum(result.x) - margin * x_range, maximum(result.x) + margin * x_range)
    ylims!(ax3, minimum(result.y) - margin * y_range, maximum(result.y) + margin * y_range)
    lines!(ax3, result.x, result.y, linewidth = 1, color = (:blue, 0.2))
    lines!(ax3, x_trace, y_trace, linewidth = 2, color = :blue)
    scatter!(ax3, [result.x[1]], [result.y[1]], color = :green, markersize = 15)
    scatter!(ax3, x_cur, y_cur, color = :red, markersize = 12)

    display(fig)

    # Animation loop
    running = Ref(true)
    on(events(fig).keyboardbutton) do event
        if event.action == Keyboard.press
            if event.key == Keyboard.escape || event.key == Keyboard.q
                running[] = false
            end
        end
    end

    dt_sim = result.t[2] - result.t[1]
    frame_dt = dt_sim / speed
    skip = max(1, round(Int, 1 / (fps * frame_dt)))

    println("Animation running... Press Escape, Q, close window, or Ctrl+C to stop")
    is_running() = running[] && !isempty(fig.scene.current_screens)

    try
        while is_running()
            for i in 1:skip:n
                !is_running() && break
                idx[] = i
                sleep(1 / fps)
            end
            if is_running()
                idx[] = n
                sleep(0.3)
            end
        end
    catch e
        e isa InterruptException || rethrow(e)
    end
    println("Animation stopped")
    fig
end

function plot_trajectory(; ω = 1.0, j0_zero = 1, x_shift = 0.0, y_shift = 0.0, tilt = 0.0, φ = 0.0, n_periods = 1)
    T = n_periods * 2π / ω
    result = simulate_cosine_steering(; ω, j0_zero, x_shift, y_shift, tilt, φ, T)

    fig = Figure(size = (1200, 800))

    # 2D trajectory
    params_str = "x=$(round(x_shift, digits=2)), y=$(round(y_shift, digits=3)), tilt=$(round(tilt, digits=3))"
    ax1 = Axis(fig[1:2, 1],
        xlabel = "x", ylabel = "y",
        title = "2D Trajectory ($params_str)",
        aspect = DataAspect())
    lines!(ax1, result.x, result.y, linewidth = 2, color = :blue)
    scatter!(ax1, [result.x[1]], [result.y[1]], color = :green, markersize = 15, label = "Start")
    scatter!(ax1, [result.x[end]], [result.y[end]], color = :red, markersize = 15, label = "End")
    axislegend(ax1, position = :lt)

    # Steering vs x
    ax2 = Axis(fig[1, 2],
        xlabel = "x", ylabel = "ψ̇ (steering, rad/s)",
        title = "Steering vs x")
    lines!(ax2, result.x, result.ψ_dot, linewidth = 2, color = :coral)
    scatter!(ax2, [result.x[1]], [result.ψ_dot[1]], color = :green, markersize = 12)
    scatter!(ax2, [result.x[end]], [result.ψ_dot[end]], color = :red, markersize = 12)

    # Steering vs y
    ax3 = Axis(fig[2, 2],
        xlabel = "y", ylabel = "ψ̇ (steering, rad/s)",
        title = "Steering vs y")
    lines!(ax3, result.y, result.ψ_dot, linewidth = 2, color = :mediumpurple)
    scatter!(ax3, [result.y[1]], [result.ψ_dot[1]], color = :green, markersize = 12)
    scatter!(ax3, [result.y[end]], [result.ψ_dot[end]], color = :red, markersize = 12)

    # Steering vs time
    ax4 = Axis(fig[3, 1], xlabel = "t", ylabel = "ψ̇ (rad/s)", title = "Steering")
    lines!(ax4, result.t, result.ψ_dot, linewidth = 2, color = :orange)

    # Heading vs time
    ax5 = Axis(fig[3, 2], xlabel = "t", ylabel = "ψ (rad)", title = "Heading")
    lines!(ax5, result.t, result.ψ, linewidth = 2, color = :purple)

    display(fig)
    fig
end

function animate_trajectory(; ω = 1.0, j0_zero = 1, x_shift = 0.0, y_shift = 0.0, tilt = 0.0, φ = 0.0,
                            n_periods = 1, fps = 30, speed = 1.0)
    T = n_periods * 2π / ω
    result = simulate_cosine_steering(; ω, j0_zero, x_shift, y_shift, tilt, φ, T)
    n = length(result.t)

    # Observables for current index
    idx = Observable(1)

    # Observables for traces (up to current point)
    x_trace = @lift result.x[1:$idx]
    y_trace = @lift result.y[1:$idx]
    ψ_trace = @lift result.ψ[1:$idx]
    ψ_dot_trace = @lift result.ψ_dot[1:$idx]
    t_trace = @lift result.t[1:$idx]

    # Observables for current point
    x_cur = @lift [result.x[$idx]]
    y_cur = @lift [result.y[$idx]]
    ψ_cur = @lift [result.ψ[$idx]]
    ψ_dot_cur = @lift [result.ψ_dot[$idx]]
    t_cur = @lift [result.t[$idx]]

    fig = Figure(size = (1200, 800))

    # Margins for axis limits
    margin = 0.1
    x_range = maximum(result.x) - minimum(result.x)
    y_range = maximum(result.y) - minimum(result.y)
    ψ_dot_margin = 0.1 * (maximum(result.ψ_dot) - minimum(result.ψ_dot))

    # 2D trajectory
    params_str = "x=$(round(x_shift, digits=2)), y=$(round(y_shift, digits=3)), tilt=$(round(tilt, digits=3))"
    ax1 = Axis(fig[1:2, 1],
        xlabel = "x", ylabel = "y",
        title = "2D Trajectory ($params_str)",
        aspect = DataAspect())
    xlims!(ax1, minimum(result.x) - margin * x_range, maximum(result.x) + margin * x_range)
    ylims!(ax1, minimum(result.y) - margin * y_range, maximum(result.y) + margin * y_range)
    lines!(ax1, result.x, result.y, linewidth = 1, color = (:blue, 0.2))
    lines!(ax1, x_trace, y_trace, linewidth = 2, color = :blue)
    scatter!(ax1, [result.x[1]], [result.y[1]], color = :green, markersize = 15)
    scatter!(ax1, x_cur, y_cur, color = :red, markersize = 12)

    # Steering vs x
    ax2 = Axis(fig[1, 2],
        xlabel = "x", ylabel = "ψ̇ (steering, rad/s)",
        title = "Steering vs x")
    xlims!(ax2, minimum(result.x) - margin * x_range, maximum(result.x) + margin * x_range)
    ylims!(ax2, minimum(result.ψ_dot) - ψ_dot_margin, maximum(result.ψ_dot) + ψ_dot_margin)
    lines!(ax2, result.x, result.ψ_dot, linewidth = 1, color = (:coral, 0.2))
    lines!(ax2, x_trace, ψ_dot_trace, linewidth = 2, color = :coral)
    scatter!(ax2, [result.x[1]], [result.ψ_dot[1]], color = :green, markersize = 12)
    scatter!(ax2, x_cur, ψ_dot_cur, color = :red, markersize = 12)

    # Steering vs y
    ax3 = Axis(fig[2, 2],
        xlabel = "y", ylabel = "ψ̇ (steering, rad/s)",
        title = "Steering vs y")
    xlims!(ax3, minimum(result.y) - margin * y_range, maximum(result.y) + margin * y_range)
    ylims!(ax3, minimum(result.ψ_dot) - ψ_dot_margin, maximum(result.ψ_dot) + ψ_dot_margin)
    lines!(ax3, result.y, result.ψ_dot, linewidth = 1, color = (:mediumpurple, 0.2))
    lines!(ax3, y_trace, ψ_dot_trace, linewidth = 2, color = :mediumpurple)
    scatter!(ax3, [result.y[1]], [result.ψ_dot[1]], color = :green, markersize = 12)
    scatter!(ax3, y_cur, ψ_dot_cur, color = :red, markersize = 12)

    # Steering vs time
    ax4 = Axis(fig[3, 1], xlabel = "t", ylabel = "ψ̇ (rad/s)", title = "Steering")
    xlims!(ax4, 0, T)
    lines!(ax4, result.t, result.ψ_dot, linewidth = 1, color = (:orange, 0.2))
    lines!(ax4, t_trace, ψ_dot_trace, linewidth = 2, color = :orange)
    scatter!(ax4, t_cur, ψ_dot_cur, color = :red, markersize = 12)

    # Heading vs time
    ax5 = Axis(fig[3, 2], xlabel = "t", ylabel = "ψ (rad)", title = "Heading")
    xlims!(ax5, 0, T)
    lines!(ax5, result.t, result.ψ, linewidth = 1, color = (:purple, 0.2))
    lines!(ax5, t_trace, ψ_trace, linewidth = 2, color = :purple)
    scatter!(ax5, t_cur, ψ_cur, color = :red, markersize = 12)

    display(fig)

    # Track if we should stop
    running = Ref(true)

    # Listen for keyboard events (Escape or Q to quit)
    on(events(fig).keyboardbutton) do event
        if event.action == Keyboard.press
            if event.key == Keyboard.escape || event.key == Keyboard.q
                running[] = false
            end
        end
    end

    # Animate
    dt_sim = result.t[2] - result.t[1]
    frame_dt = dt_sim / speed
    skip = max(1, round(Int, 1 / (fps * frame_dt)))

    println("Animation running... Press Escape, Q, close window, or Ctrl+C to stop")

    # Helper to check if we should keep running
    is_running() = running[] && !isempty(fig.scene.current_screens)

    try
        while is_running()
            for i in 1:skip:n
                !is_running() && break
                idx[] = i
                sleep(1 / fps)
            end
            if is_running()
                idx[] = n
                sleep(0.3)  # Brief pause at end before looping
            end
        end
    catch e
        if e isa InterruptException
            println("\nAnimation stopped")
        else
            rethrow(e)
        end
    end
    println("Animation stopped")

    fig
end

function compare_j0_zeros()
    fig = Figure(size = (1200, 400))

    for i in 1:3
        result = simulate_cosine_steering(; j0_zero = i, T = 2π)
        B = J0_ZEROS[i]
        ax = Axis(fig[1, i],
            xlabel = "x", ylabel = "y",
            title = "J₀ zero #$i (A/ω = $(round(B, digits=2)))",
            aspect = DataAspect())
        lines!(ax, result.x, result.y, linewidth = 2)
        scatter!(ax, [result.x[1]], [result.y[1]], color = :green, markersize = 10)
    end

    Label(fig[0, :], "Closed Trajectories at Different J₀ Zeros", fontsize = 20)

    display(fig)
    fig
end

function compare_deltas(; j0_zero = 1)
    fig = Figure(size = (1200, 400))

    deltas = [-0.05, 0.0, 0.05]

    for (i, δ) in enumerate(deltas)
        result = simulate_cosine_steering(; j0_zero, δ, T = 4π)
        B = J0_ZEROS[j0_zero] * (1 + δ)
        ax = Axis(fig[1, i],
            xlabel = "x", ylabel = "y",
            title = "δ = $δ (A/ω = $(round(B, digits=2)))",
            aspect = DataAspect())
        lines!(ax, result.x, result.y, linewidth = 2)
        scatter!(ax, [result.x[1]], [result.y[1]], color = :green, markersize = 10)
    end

    Label(fig[0, :], "Effect of δ Perturbation (J₀ zero #$j0_zero)", fontsize = 20)

    display(fig)
    fig
end

# Animate closed trajectory by default when included
animate_x_based(;n_periods=1, fps=1000, dt=0.001)
