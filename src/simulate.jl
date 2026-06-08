# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

"""
    sim!(sam, set_values; dt, total_time, vsm_interval,
         prn, lin_model, y_op)

Run a generic simulation with a matrix of control inputs.
Optionally, also simulate a provided linear model.

# Arguments
- `sam::SymbolicAWEModel`: Initialized AWE model.
- `set_values::Matrix{Float64}`: Control torque offsets [Nm]
  per step. Rows = steps, columns = winches.

# Keywords
- `dt`: Time step [s]. Default `1/sam.set.sample_freq`.
- `total_time`: Simulation duration [s]. Default 10.0.
- `vsm_interval`: Steps between VSM updates. Default 3.
- `prn`: Print performance summary. Default true.
- `lin_model`: Optional `StateSpace` for linear comparison.
- `y_op`: Operating point output vector. Required when
  `lin_model` is provided.

# Returns
- `(SysLog, Nothing)` or `(SysLog, SysLog)` when `lin_model`
  is provided.
"""
function sim!(
    sam::SymbolicAWEModel,
    set_values::Matrix{Float64};
    dt=1/sam.set.sample_freq,
    total_time=10.0,
    vsm_interval=3,
    prn=true,
    lin_model::Union{Nothing, <:NamedTuple, StateSpace}=nothing,
    y_op::Union{Nothing, AbstractVector}=nothing,
    torque_damp=0.9,
)
    steps = Int(round(total_time / dt))
    if size(set_values, 1) != steps
        error("The number of rows in set_values ($(size(set_values, 1))) must match the number of simulation steps ($steps).")
    end
    if lin_model isa NamedTuple
        lin_model = ss(lin_model...)
    end

    logger = Logger(sam, steps+1)
    sys_state = SysState(sam)
    sys_state.time = 0.0
    # log!(logger, sys_state)

    if prn
        @info "Starting nonlinear simulation..."
    end
    step_time = 0.0
    vsm_time = 0.0
    integ_time = 0.0
    set_torques = similar(set_values)

    steady_torque = calc_steady_torque(sam)

    # --- Nonlinear Simulation Loop ---
    elapsed = @elapsed for step in 1:steps
        t = step * dt

        steady_torque = torque_damp * steady_torque + (1-torque_damp) * calc_steady_torque(sam)
        set_torques[step, :] = steady_torque .+ set_values[step, :]

        try
            step_time += @elapsed next_step!(sam;
                                             set_values=set_torques[step, :],
                                             dt, vsm_interval=vsm_interval)
            integ_time += sam.t_step
            vsm_time += sam.t_vsm
            if step < steps ÷ 2
                step_time, integ_time, vsm_time = 0.0, 0.0, 0.0
            end
        catch exception
            if exception isa AssertionError
                if prn
                    @warn "Crashed at t=$t"
                end
                break
            else
                rethrow(exception)
            end
        end
        update_sys_state!(sys_state, sam)
        sys_state.time = t
        log!(logger, sys_state)
    end

    # --- Linear Simulation ---
    lin_sim_log = nothing
    if !isnothing(lin_model)
        isnothing(y_op) && error(
            "y_op (operating point output) is required " *
            "when lin_model is provided")
        t_vec = 0:dt:(total_time - dt)
        ΔU = permutedims(set_torques) .- set_torques[1, :]
        lin_res = lsim(lin_model, ΔU, t_vec)

        # Reconstruct full output from deviation output:
        # y = y_op + Δy
        ΔY = lin_res.y
        lin_y_full = ΔY .+ y_op

        # Log the complete linear simulation result
        lin_logger = Logger(sam, steps)
        lin_sys_state = SysState(sam)
        update_sys_state!(lin_sys_state, collect(y_op), sam, t_vec[1])
        for step in 1:steps
            y_k = lin_y_full[:, step]
            update_sys_state!(
                lin_sys_state, y_k, sam, t_vec[step])
            log!(lin_logger, lin_sys_state)
        end
        save_log(lin_logger, "tmp_run_lin")
        lin_sim_log = load_log("tmp_run_lin")
    end

    mkpath(get_data_path())
    save_log(logger, "tmp_run")
    sim_log = load_log("tmp_run")
    
    if prn
        lines = [
            "Performance Summary:",
            @sprintf("%-12s | %-12s | %-10s", "Component", "Speedup (×)", "Total Time"),
            "-------------|--------------|------------",
            @sprintf("%-12s | %12.2f | %10.2f", "Simulation", total_time / elapsed / 2, elapsed),
            @sprintf("%-12s | %12.2f | %10.2f", "Step",       total_time / step_time / 2, step_time),
            @sprintf("%-12s | %12.2f | %10.2f", "Integrator", total_time / integ_time / 2, integ_time),
            @sprintf("%-12s | %12.2f | %10.2f", "VSM",        total_time / vsm_time / 2, vsm_time)
        ]
        @info join(lines, "\n")
    end

    return (sim_log, lin_sim_log)
end

"""
    sim_reposition!(sam; dt, total_time, reposition_interval_s, target_elevation_deg,
                    target_azimuth_deg, prn)

Run a simulation that periodically resets the kite's elevation and azimuth.

This function simulates the AWE model and, at a specified time interval, calls
`reposition!` to reposition the kite to a target elevation and azimuth. It logs
the entire simulation and returns a `SysLog`.

# Arguments
- `sam::SymbolicAWEModel`: Initialized AWE model.

# Keywords
- `dt::Float64`: Time step [s]. Defaults to `1/sam.set.sample_freq`.
- `total_time::Float64`: Total simulation duration [s]. Defaults to 20.0.
- `reposition_interval_s::Float64`: The interval in seconds at which to reset the pose. Defaults to 5.0.
- `target_elevation::Float64`: The target elevation in rad for repositioning. Defaults to deg2rad(45.0).
- `target_azimuth::Float64`: The target azimuth in rad for repositioning. Defaults to 0.0.
- `target_heading::Float64`: The target heading in rad for repositioning. Defaults to 0.0.
- `prn::Bool`: If true, prints status messages during the simulation. Defaults to true.

# Returns
- `SysLog`: A log of the complete simulation.
"""
function sim_reposition!(
    sam::SymbolicAWEModel;
    dt=1/sam.set.sample_freq,
    total_time=20.0,
    reposition_interval_s=5.0,
    target_elevation=deg2rad(45.0),
    target_azimuth=0.0,
    target_heading=0.0,
    prn=true
)
    # 1. --- Initialization ---
    sys_struct = sam.sys_struct
    steps = Int(round(total_time / dt))
    reposition_interval_steps = Int(round(reposition_interval_s / dt))
    set_values = zeros(Float64, steps, length(sys_struct.winches))
    vsm_interval = 1 ÷ dt
    
    logger = Logger(sam, steps+1)
    sys_state = SysState(sam)
    sys_state.time = 0.0
    log!(logger, sys_state)

    if prn
        println("--- Starting simulation with periodic repositioning ---")
        println("Total time: $(total_time)s, Reposition interval: $(reposition_interval_s)s")
    end

    # 2. --- Simulation Loop ---
    time = @elapsed for step in 1:steps
        t = step * dt
        
        # Hold the kite in place by countering the tether forces with winch torques
        set_values[step, :] = -sam.set.drum_radius .* [norm(winch.force) for winch in sys_struct.winches]
        
        # --- Repositioning Logic ---
        if step > 1 && (step - 1) % reposition_interval_steps == 0
            if prn
                println("\n>>> Time: $(round(t, digits=2))s. Repositioning kite...")
                println(">>> Target Elevation: $(rad2deg(target_elevation))°,"*
                        " Target Azimuth: $(rad2deg(target_azimuth))°")
            end
            
            # Update the transform with the new target pose
            sys_struct.transforms[1].elevation = target_elevation
            sys_struct.transforms[1].azimuth   = target_azimuth
            sys_struct.transforms[1].heading   = target_heading
            
            # Apply the transformation without changing velocities
            SymbolicAWEModels.reposition!(sys_struct.transforms, sys_struct)
            
            # Reinitialize the solver to handle the state discontinuity
            local_prob = sam.prob
            if local_prob isa ProbWithAttributes
                SymbolicAWEModels.reinit!(sam, local_prob, FBDF())
            end

            if prn
                # Verify the new pose after one step
                next_step!(sam; dt=dt, set_values=set_values[step, :], vsm_interval)
                updated_elevation_deg = rad2deg(sys_struct.wings[1].elevation)
                updated_azimuth_deg = rad2deg(sys_struct.wings[1].azimuth)
                updated_heading_deg = rad2deg(sys_struct.wings[1].heading)
                println(">>> Pose updated." *
                        " Elevation: $(round(updated_elevation_deg, digits=2))°, " *
                        " Azimuth: $(round(updated_azimuth_deg, digits=2))°, " *
                        " Heading: $(round(updated_heading_deg, digits=2))°.\n")
            end
        else
            # --- Normal simulation step ---
            next_step!(sam; dt=dt, set_values=set_values[step, :], vsm_interval)
        end
        
        # Log the state at the current time step
        update_sys_state!(sys_state, sam)
        sys_state.time = t
        log!(logger, sys_state)
    end

    if prn
        println("--- Simulation Finished ---")
        println("Runtime: $time")
        println("Times realtime: $(dt*steps/time)")
    end

    # Save and return the log
    mkpath(get_data_path())
    save_log(logger, "tmp_reposition_run")
    return load_log("tmp_reposition_run")
end


"""
    make_lin_sys_state(y::AbstractVector, sam::SymbolicAWEModel, t::Real)

Construct a SysState for logging linear state-space simulation output y (ordered as
sam.outputs).
"""
function make_lin_sys_state(y::AbstractVector, sam::SymbolicAWEModel, t::Real)
    sys_state = SysState(sam)
    update_sys_state!(sys_state, y, sam, t)
    return sys_state
end

"""
    update_sys_state!(sys_state, y::AbstractVector, sam::SymbolicAWEModel, t::Real)

Update a SysState for a linear state-space simulation, using output y and model sam.
"""
function update_sys_state!(sys_state, y::AbstractVector, sam::SymbolicAWEModel, t::Real)
    sys = sam.prob.sys
    outputs = sam.outputs
    for (i, sym) in enumerate(outputs)
        if isequal(sym, sys.heading[1])
            sys_state.heading = y[i]
        elseif isequal(sym, sys.turn_rate[3, 1])
            sys_state.turn_rates[1,3] = y[i]
        elseif isequal(sym, sys.tether_len[1])
            sys_state.l_tether[1] = y[i]
        elseif isequal(sym, sys.tether_len[2])
            sys_state.l_tether[2] = y[i]
        elseif isequal(sym, sys.tether_len[3])
            sys_state.l_tether[3] = y[i]
        elseif isequal(sym, sys.winch_vel[1])
            sys_state.v_reelout[1] = y[i]
        elseif isequal(sym, sys.winch_vel[2])
            sys_state.v_reelout[2] = y[i]
        elseif isequal(sym, sys.winch_vel[3])
            sys_state.v_reelout[3] = y[i]
        elseif isequal(sym, sys.winch_force[1])
            sys_state.winch_force[1] = y[i]
        elseif isequal(sym, sys.winch_force[2])
            sys_state.winch_force[2] = y[i]
        elseif isequal(sym, sys.winch_force[3])
            sys_state.winch_force[3] = y[i]
        elseif isequal(sym, sys.angle_of_attack[1])
            sys_state.AoA = mod(y[i] + π, 2π) - π  # Wrap to [-π, π]
        elseif isequal(sym, sys.elevation[1])
            sys_state.elevation = y[i]
        elseif isequal(sym, sys.azimuth[1])
            sys_state.azimuth = y[i]
        elseif isequal(sym, sys.course[1])
            sys_state.course = y[i]
        end
    end
    sys_state.time = t
    return sys_state
end
