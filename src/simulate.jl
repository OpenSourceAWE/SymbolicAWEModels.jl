# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

"""
    sim!(sam, set_values; dt, total_time, vsm_interval, prn, lin_model)

Run a generic simulation for a given AWE model and a matrix of control inputs.
Optionally, also simulate a provided linear model, returning both logs.

# Arguments
- `sam::SymbolicAWEModel`: Initialized AWE model.
- `set_values::Matrix{Float64}`: A matrix of external control torques [Nm] to be
  applied at each time step. The number of rows must equal the number of
  simulation steps, and the number of columns must equal the number of winches.

# Keywords
- `dt::Float64`: Time step [s]. Defaults to `1/sam.set.sample_freq`.
- `total_time::Float64`: Total simulation duration [s]. Defaults to 10.0.
- `vsm_interval::Int`: Interval for the value state machine updates. Defaults to 3.
- `prn::Bool`: If true, prints a performance summary upon completion. Defaults to true.
- `lin_model`: (optional) a continuous-time `StateSpace` object from `ControlSystemsBase`.
    If provided, the linear model is simulated in parallel and a second log is returned.

# Returns
- If `lin_model` is not provided: `(SysLog, Nothing)` (nonlinear log, nothing)
- If `lin_model` is provided: `(SysLog, SysLog)` (nonlinear, linear logs)
"""
function sim!(
    sam::SymbolicAWEModel,
    set_values::Matrix{Float64};
    dt=1/sam.set.sample_freq,
    total_time=10.0,
    vsm_interval=3,
    prn=true,
    lin_model::Union{Nothing, <:NamedTuple, StateSpace}=nothing
)
    steps = Int(round(total_time / dt))
    sys_struct = sam.sys_struct
    if size(set_values, 1) != steps
        error("The number of rows in set_values ($(size(set_values, 1))) must match the number of simulation steps ($steps).")
    end
    if lin_model isa NamedTuple
        lin_model = ss(lin_model...)
    end

    logger = Logger(length(sys_struct.points), steps)
    sys_state = SysState(sam)

    if prn
        @info "Starting nonlinear simulation..."
    end
    step_time = 0.0
    vsm_time = 0.0
    integ_time = 0.0
    set_torques = similar(set_values)
    y_op = sam.simple_lin_model.get_y(sam.integrator)

    # --- Nonlinear Simulation Loop ---
    elapsed = @elapsed for step in 1:steps
        t = (step-1) * dt

        set_torques[step, :] = -sam.set.drum_radius .* [norm(winch.force) for winch in sys_struct.winches]
        set_torques[step, :] .+= set_values[step, :]

        try
            step_time += @elapsed next_step!(sam;
                                             set_values=set_torques[step, :],
                                             dt, vsm_interval=vsm_interval)
            integ_time += sam.t_step
            vsm_time += sam.t_vsm
            if step < steps ÷ 2
                step_time, integ_time, vsm_time = 0.0, 0.0, 0.0
            end
        catch e
            if e isa AssertionError
                if prn
                    @warn "Crashed at t=$t"
                end
                break
            else
                rethrow(e)
            end
        end
        update_sys_state!(sys_state, sam)
        sys_state.time = t
        log!(logger, sys_state)
    end

    # --- Linear Simulation ---
    lin_lg = nothing
    if !isnothing(lin_model)
        t_vec = 0:dt:(total_time - dt)
        ΔU = permutedims(set_torques) .- set_torques[1, :]
        lin_res = lsim(lin_model, ΔU, t_vec)
        
        # Reconstruct full output from deviation output: y = y_op + Δy
        ΔY = lin_res.y
        lin_y_full = ΔY .+ y_op
        
        # Log the complete linear simulation result
        n_points = length(sys_struct.points)
        lin_logger = Logger(n_points, steps)
        lin_sys_state = SysState(y_op, sam, t_vec[1])
        for step in 1:steps
            y_k = lin_y_full[:, step]
            update_sys_state!(lin_sys_state, y_k, sam, t_vec[step])
            log!(lin_logger, lin_sys_state)
        end
        save_log(lin_logger, "tmp_run_lin")
        lin_lg = load_log("tmp_run_lin")
    end

    mkpath(get_data_path())
    save_log(logger, "tmp_run")
    lg = load_log("tmp_run")
    
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

    return (lg, lin_lg)
end

"""
    sim_oscillate!(sam; dt, total_time, steering_freq, steering_magnitude, vsm_interval,
                   bias, prn, lin_model)

Run a simulation with oscillating steering input on the given AWE model.
Optionally also simulate a provided linear model.

# Keywords (see sim!)
- `lin_model`: (optional) a continuous-time `StateSpace` object from `ControlSystemsBase`.

# Returns
- If `lin_model` is not provided: `(SysLog, Nothing)` (nonlinear log, nothing)
- If `lin_model` is provided: `(SysLog, SysLog)` (nonlinear, linear logs)
"""
function sim_oscillate!(
    sam::SymbolicAWEModel;
    dt=1/sam.set.sample_freq,
    total_time=10.0,
    steering_freq=0.5,
    steering_magnitude=10.0,
    vsm_interval=3,
    bias = 0.0,
    prn=false,
    lin_model=nothing
)
    sys_struct = sam.sys_struct
    steps = Int(round(total_time / dt))
    num_winches = length(sys_struct.winches)
    @assert num_winches == 3
    set_values = zeros(Float64, steps, num_winches)

    if prn
        @info "Simulating using oscillating steering inputs\n" *
              "\twith frequency = $steering_freq Hz\n" *
              "\tand magnitude = $steering_magnitude N."
    end

    for step in 1:steps
        t = (step-1) * dt
        steering = steering_magnitude * cos(2π * steering_freq * t + bias)
        set_values[step, :] = [0.0, steering, -steering]
    end

    return sim!(sam, set_values; dt=dt, total_time=total_time, vsm_interval=vsm_interval,
                prn=prn, lin_model=lin_model)
end

"""
    sim_turn!(sam; dt, total_time, steering_time, steering_magnitude, vsm_interval, prn,
              lin_model)

Run a simulation with a constant steering input for a specified duration.
Optionally also simulate a provided linear model.

# Keywords (see sim!)
- `lin_model`: (optional) a continuous-time `StateSpace` object from `ControlSystemsBase`.

# Returns
- If `lin_model` is not provided: `(SysLog, Nothing)` (nonlinear log, nothing)
- If `lin_model` is provided: `(SysLog, SysLog)` (nonlinear, linear logs)
"""
function sim_turn!(
    sam::SymbolicAWEModel;
    dt=1/sam.set.sample_freq,
    total_time=10.0,
    steering_time=2.0,
    steering_magnitude=10.0,
    vsm_interval=3,
    prn=false,
    lin_model=nothing
)
    sys_struct = sam.sys_struct
    steps = Int(round(total_time / dt))
    steering_steps = Int(round(steering_time / dt))
    num_winches = length(sys_struct.winches)
    @assert num_winches == 3
    set_values = zeros(Float64, steps, num_winches)

    if prn
        @info "Generating turn commands..."
    end

    for step in 1:steps
        if step <= steering_steps
            set_values[step, :] = [0.0, steering_magnitude, -steering_magnitude]
        else
            set_values[step, :] = zeros(num_winches)
        end
    end

    return sim!(sam, set_values; dt=dt, total_time=total_time, vsm_interval=vsm_interval,
                prn=prn, lin_model=lin_model)
end


"""
    SysState(y::AbstractVector, sam::SymbolicAWEModel, t::Real; zoom=1.0)

Construct a SysState for logging linear state-space simulation output y (ordered as
sam.outputs).
"""
function SysState(y::AbstractVector, sam::SymbolicAWEModel, t::Real; zoom=1.0)
    P = length(sam.sys_struct.points)
    ss = SysState{P}()
    update_sys_state!(ss, y, sam, t; zoom)
    return ss
end

"""
    update_sys_state!(ss::SysState, y::AbstractVector, sam::SymbolicAWEModel, t::Real;
                      zoom=1.0)

Update a SysState for a linear state-space simulation, using output y and model sam.
"""
function update_sys_state!(ss::SysState, y::AbstractVector, sam::SymbolicAWEModel, t::Real;
                           zoom=1.0)
    sys = sam.sys
    outputs = sam.outputs
    for (i, sym) in enumerate(outputs)
        if isequal(sym, sys.heading[1])
            ss.heading = y[i]
        elseif isequal(sym, sys.turn_rate[1,3])
            ss.turn_rates[1,3] = y[i]
        elseif isequal(sym, sys.tether_len[1])
            ss.l_tether[1] = y[i]
        elseif isequal(sym, sys.tether_len[2])
            ss.l_tether[2] = y[i]
        elseif isequal(sym, sys.tether_len[3])
            ss.l_tether[3] = y[i]
        elseif isequal(sym, sys.tether_vel[1])
            ss.v_reelout[1] = y[i]
        elseif isequal(sym, sys.tether_vel[2])
            ss.v_reelout[2] = y[i]
        elseif isequal(sym, sys.tether_vel[3])
            ss.v_reelout[3] = y[i]
        elseif isequal(sym, sys.winch_force[1])
            ss.force[1] = y[i]
        elseif isequal(sym, sys.winch_force[2])
            ss.force[2] = y[i]
        elseif isequal(sym, sys.winch_force[3])
            ss.force[3] = y[i]
        elseif isequal(sym, sys.angle_of_attack[1])
            ss.AoA = y[i]
        elseif isequal(sym, sys.elevation[1])
            ss.elevation = y[i]
        elseif isequal(sym, sys.azimuth[1])
            ss.azimuth = y[i]
        elseif isequal(sym, sys.course[1])
            ss.course = y[i]
        end
    end
    ss.time = t
    return ss
end
