"""
CSV Replay with Equilibrium Turn Rate Analysis

Based on csv_replay.jl. At the end of simulation, uses NonlinearSolve to find
the turn rate that would cause VSM moment[3] to be zero, and compares to
the actual coupled simulation turn rate and CSV data.

Usage:
    julia --project=examples examples/coupled/csv_replay_equilibrium.jl
"""

using SymbolicAWEModels, VortexStepMethod, KiteUtils
using SymbolicAWEModels: reposition!, rotate_around_z, rotate_around_y, calc_steady_torque
using CSV, DataFrames
using GLMakie
using Statistics
using Rotations
using UnPack
using LinearAlgebra
using OrdinaryDiffEqBDF
using KiteUtils: calc_elevation, azimuth_east
using NonlinearSolve, ADTypes

# Configuration parameters
CSV_PATH = "data/v3/2025-10-09_16-58-33_ProtoLogger_lidar.csv"
START_FRAME = 22068 + 120 # First frame to replay
END_FRAME = START_FRAME + 50 # Last frame to replay (use nothing for all frames)
REMAKE_CACHE = false
WING_CD = 0.0  # Drag coefficient for WING type points

# V3 Kite steering/depower calibration (from KCU documentation)
# Steering calibration
STEERING_L0 = 1.6  # Neutral steering tape length (m)
STEERING_GAIN = 1.2  # Maximum differential (m) at |u_s| = 1
STEERING_MULTIPLIER = 1.0

# Depower calibration
DEPOWER_L0 = 0.2 # SUPPOSED TO BE 0.2
DEPOWER_GAIN = 5.0
DEPOWER_OFFSET = -15.0

INITIAL_DAMPING = [0.0, 300.0, 600.0]
DECAY_TIME = 1.0
MIN_DAMPING = [0.0, 30, 60]



"""
    steering_percentage_to_lengths(percentage)

Convert CSV steering percentage to left/right tape lengths (m).
Percentage convention: negative = left turn, positive = right turn.
"""
function steering_percentage_to_lengths(percentage)
    u_s = percentage / 100.0  # Convert percentage to [-1, 1]
    L_left = STEERING_L0 + STEERING_GAIN * STEERING_MULTIPLIER * u_s
    L_right = STEERING_L0 - STEERING_GAIN * STEERING_MULTIPLIER * u_s
    return L_left, L_right
end

"""
    depower_percentage_to_length(percentage)

Convert CSV depower percentage to tape length (m).
"""
function depower_percentage_to_length(percentage)
    u_p = percentage / 100.0  # Convert percentage to [0, 1]
    L_depower = DEPOWER_L0 + DEPOWER_GAIN * u_p
    return L_depower
end

"""
    steering_length_to_percentage(L_right)

Convert right steering tape length (m) back to steering percentage.
"""
function steering_length_to_percentage(L_right)
    u_s = (STEERING_L0 - L_right) / (STEERING_GAIN * STEERING_MULTIPLIER)
    return u_s * 100.0
end

"""
    depower_length_to_percentage(L_depower)

Convert depower tape length (m) back to depower percentage.
"""
function depower_length_to_percentage(L_depower)
    u_p = (L_depower - DEPOWER_L0) / DEPOWER_GAIN
    return u_p * 100.0 - DEPOWER_OFFSET
end

"""
    load_flight_data(csv_path::String)

Load and parse CSV flight data using space/multi-space delimiter.
Returns a named tuple with all CSV columns as time series data.
"""
function load_flight_data(csv_path::String)
    @info "Loading CSV data from: $csv_path"
    df = CSV.read(csv_path, DataFrame;
                  delim=' ',
                  silencewarnings=true,
                  normalizenames=true,
                  types=Float64,
                  strict=false)
    t0 = df.time[START_FRAME]
    df.time .= df.time .- t0
    col_names = Tuple(Symbol(name) for name in names(df))
    data = NamedTuple{col_names}(Tuple(eachcol(df)))
    return data
end

"""
    limit_frames(data; start_frame=1, end_frame=nothing)

Limit data to frame range [start_frame, end_frame].
Automatically slices all fields in the named tuple.
"""
function limit_frames(data; start_frame=1, end_frame=nothing)
    n_total = length(data.time)
    start_idx = max(1, min(start_frame, n_total))
    if isnothing(end_frame)
        end_idx = n_total
    else
        end_idx = max(start_idx, min(end_frame, n_total))
    end
    limited = NamedTuple{keys(data)}(
        Tuple(field[start_idx:end_idx] for field in data)
    )
    return limited
end

"""
    add_distance_column(data)

Add distance and cumulative_distance columns to CSV data.
Calculates 3D Euclidean distance between consecutive kite positions.
"""
function add_distance_column(data)
    n = length(data.time)
    distances = zeros(Float64, n)
    cumulative_distances = zeros(Float64, n)

    for i in 2:n
        dx = data.kite_pos_east[i] - data.kite_pos_east[i-1]
        dy = data.kite_pos_north[i] - data.kite_pos_north[i-1]
        dz = data.kite_height[i] - data.kite_height[i-1]
        distances[i] = sqrt(dx^2 + dy^2 + dz^2)
        cumulative_distances[i] = cumulative_distances[i-1] + distances[i]
    end

    # Add new fields to the named tuple
    return merge(data, (distance=distances, cumulative_distance=cumulative_distances))
end

"""
    interpolate_csv_data(data, n_substeps)

Linearly interpolate CSV data to create n_substeps points between each pair
of original data points. Handles missing values by propagating them.
"""
function interpolate_csv_data(data, n_substeps)
    n_original = length(data.time)
    n_interp = (n_original - 1) * n_substeps + 1

    # Create interpolated arrays for each field
    interp_data = Dict{Symbol, Vector}()

    for field in keys(data)
        # Determine element type (handle Union{Float64, Missing})
        eltype_field = eltype(data[field])
        interp_values = Vector{eltype_field}(undef, n_interp)

        for i in 1:(n_original-1)
            # Starting index in interpolated array
            start_idx = (i-1) * n_substeps + 1

            val_i = data[field][i]
            val_next = data[field][i+1]

            # Interpolate between data[i] and data[i+1]
            for j in 0:(n_substeps-1)
                idx = start_idx + j
                alpha = j / n_substeps

                # Handle missing values
                if ismissing(val_i) || ismissing(val_next)
                    interp_values[idx] = missing
                else
                    interp_values[idx] = (1 - alpha) * val_i + alpha * val_next
                end
            end
        end

        # Add the last point
        interp_values[end] = data[field][end]
        interp_data[field] = interp_values
    end

    # Convert back to named tuple
    col_names = Tuple(keys(data))
    return NamedTuple{col_names}(Tuple(interp_data[k] for k in col_names))
end

"""
    euler_to_quaternion(roll_deg, pitch_deg, yaw_deg)

Convert Euler angles (in degrees) from NED to ENU quaternion.
CSV data is in NED (North East Down) frame, but Q_b_w requires ENU (East North Up).

NED to ENU transformation:
  X_ENU = Y_NED (East)
  Y_ENU = X_NED (North)
  Z_ENU = -Z_NED (Up = -Down)
"""
function euler_to_quaternion(roll_deg, pitch_deg, yaw_deg)
    # Convert degrees to radians
    roll_rad = deg2rad(roll_deg)
    pitch_rad = deg2rad(pitch_deg)
    yaw_rad = deg2rad(yaw_deg)
    # Create rotation in NED frame using Rotations.jl (ZYX convention)
    rot_ned = RotZYX(yaw_rad, pitch_rad, roll_rad)
    R_ned_to_enu = [0.0 1.0 0.0;    # X_ENU = Y_NED (East)
                    1.0 0.0 0.0;    # Y_ENU = X_NED (North)
                    0.0 0.0 -1.0]   # Z_ENU = -Z_NED (Up)
    rot_enu = R_ned_to_enu * Matrix(rot_ned)
    q = SymbolicAWEModels.rotation_matrix_to_quaternion(rot_enu)
    return q
end

function calc_R_b_w(sys_struct::SystemStructure)
    @unpack points, wings, wind_vec_gnd = sys_struct
    wing = wings[1]
    R_b_w, origin = SymbolicAWEModels.calc_refine_wing_frame(
        points,
        wing.z_ref_points,
        wing.y_ref_points,
        wing.origin_idx
    )
    return R_b_w
end

"""
    wrap_to_pi(angle)

Wrap angle to [-π, π] range.
"""
function wrap_to_pi(angle)
    return mod(angle + π, 2π) - π
end

"""
    calc_feedforward_torque(tether_force_n, winch)

Calculate feed-forward torque from tether force.
Similar to calc_steady_torque but uses measured force.
"""
function calc_feedforward_torque(tether_force_n, winch)
    torque = -winch.drum_radius / winch.gear_ratio * tether_force_n + winch.friction
    return torque
end

"""
    calc_heading(sys_struct::SystemStructure, R_b_w)

Calculate heading angle from rotation matrix, wrapped to [-π, π].
"""
function calc_heading(sys_struct::SystemStructure, R_b_w)
    e_x = R_b_w[:, 1]
    wind_norm = [1,0,0]
    # Project -e_x onto plane perpendicular to wind
    minus_e_x = -e_x
    proj_on_wind = dot(minus_e_x, wind_norm) * wind_norm
    e_x_perp = minus_e_x - proj_on_wind
    # Heading components in wind-perpendicular plane
    wind_cross_z = [wind_norm[2], -wind_norm[1], 0]
    heading_x = dot(e_x_perp, wind_cross_z)
    heading_z = e_x_perp[3]
    heading = atan(heading_x, heading_z)
    return wrap_to_pi(heading)
end

"""
    calc_csv_heading(roll_deg, pitch_deg, yaw_deg, sys_struct)

Calculate heading from CSV Euler angles, wrapped to [-π, π].
"""
function calc_csv_heading(roll_deg, pitch_deg, yaw_deg, sys_struct)
    quat = euler_to_quaternion(roll_deg, pitch_deg, yaw_deg)
    R = SymbolicAWEModels.quaternion_to_rotation_matrix(quat)
    heading = calc_heading(sys_struct, R)
    return wrap_to_pi(heading + π)
end

# Mutable state for heading spike filter
const PREV_CSV_HEADING = Ref{Float64}(NaN)
const MAX_HEADING_CHANGE = deg2rad(20.0)  # Max allowed change per step

"""
    filter_csv_heading(heading)

Filter out spikes in CSV heading data. Rejects changes > MAX_HEADING_CHANGE.
"""
function filter_csv_heading(heading)
    if isnan(PREV_CSV_HEADING[])
        PREV_CSV_HEADING[] = heading
        return heading
    end
    delta = wrap_to_pi(heading - PREV_CSV_HEADING[])
    if abs(delta) > MAX_HEADING_CHANGE
        # Spike detected, keep previous value
        return PREV_CSV_HEADING[]
    end
    PREV_CSV_HEADING[] = heading
    return heading
end

function apply_force!(sys, control)
    wing = sys.wings[1]
    R_b_w = wing.R_b_w
    for point in sys.points
        distance_frac = point.pos_w ⋅ normalize(wing.pos_w) / norm(wing.pos_w)
        point.disturb .= -R_b_w[:, 1] * control * distance_frac
    end
end

# Mutable state for simulation cumulative distance
const SIM_PREV_POS = Ref{Vector{Float64}}(zeros(3))
const SIM_CUMULATIVE_DIST = Ref{Float64}(0.0)

function reset_distance_tracker!()
    SIM_PREV_POS[] = zeros(3)
    SIM_CUMULATIVE_DIST[] = 0.0
end

function update_sim_distance!(wing_pos)
    if SIM_PREV_POS[] == zeros(3)
        SIM_PREV_POS[] = copy(wing_pos)
        return 0.0
    end
    dist = norm(wing_pos - SIM_PREV_POS[])
    SIM_CUMULATIVE_DIST[] += dist
    SIM_PREV_POS[] = copy(wing_pos)
    return SIM_CUMULATIVE_DIST[]
end

function update_vel_from_csv!(sys, row, brake)
    @unpack wings, points, winches, segments = sys
    wing = wings[1]

    # Set wind speed from CSV
    sys.set.v_wind = row.wind_at_kite

    # Apply steering from CSV directly (no controller)
    steering = clamp(row.steering, -100.0, 100.0)
    L_left, L_right = steering_percentage_to_lengths(steering)
    segments[87].l0 = L_left
    segments[89].l0 = L_right

    # Winch with feed-forward torque from CSV tether force
    winch = winches[1]
    ff_torque = calc_feedforward_torque(row.tether_force, winch)
    winch.brake = brake
    winch.set_value = ff_torque

    # Update depower from CSV
    L_depower = depower_percentage_to_length(row.depower + DEPOWER_OFFSET)
    segments[88].l0 = L_depower

    return winch.set_value
end

"""
    update_sys_struct_from_csv!(sys_struct, csv_row_data)

Update system structure from a single CSV row.
Updates wing orientation from Euler angles and position via transform system.
Uses nonlinear solve to find rotation angle that achieves desired heading.
"""
function update_sys_struct_from_csv!(sys, row)
    @unpack wings, points, winches, segments, transforms = sys
    wing = wings[1]
    transform = transforms[1]

    # calc target heading from CSV
    quat = euler_to_quaternion(row.roll, row.pitch, row.yaw)
    csv_heading = calc_heading(sys,
        SymbolicAWEModels.quaternion_to_rotation_matrix(quat)) + pi
    wing.R_b_w = calc_R_b_w(sys)
    curr_heading = calc_heading(sys, wing.R_b_w)

    # calc needed transform
    csv_pos = [row.x, row.y, row.z]
    curr_pos = wing.pos_w
    delta_pos = csv_pos - curr_pos
    # apply transform
    for (n, point_idx) in enumerate(39:44)
        points[point_idx].pos_b .= [0.0, 0.0, -n * row.tether_len / 6 * 1.01]
    end
    transform.elevation = KiteUtils.calc_elevation(csv_pos)
    transform.azimuth = KiteUtils.azimuth_east(csv_pos)
    transform.heading = csv_heading
    SymbolicAWEModels.reinit!([transform], sys)

    # apply vel
    csv_vel = [row.vx, row.vy, row.vz]
    wing.vel_w .= csv_vel
    for point in points
        transform_frac = point.pos_w ⋅ normalize(wing.pos_w) / norm(wing.pos_w)
        point.vel_w .= transform_frac * csv_vel
    end

    # update tether length and velocity
    winches[1].brake = true

    # Convert CSV percentages to tape lengths
    L_left, L_right = steering_percentage_to_lengths(row.steering)
    L_depower = depower_percentage_to_length(row.depower + DEPOWER_OFFSET)

    segments[87].l0 = L_left   # Left steering tape
    segments[89].l0 = L_right  # Right steering tape
    segments[88].l0 = L_depower  # Depower tape
end

"""
    find_equilibrium_turn_rate(wing)

Find the turn rate that would cause VSM yaw moment (moment[3]) to be zero.
Uses NonlinearSolve to find the equilibrium turn rate.

Returns the equilibrium turn rate in rad/s.
"""
function find_equilibrium_turn_rate(wing)
    body_aero = wing.vsm_aero
    solver = wing.vsm_solver
    va_b = wing.va_b

    v_a = norm(va_b)
    if v_a < 1e-6
        @warn "Apparent velocity too small: $v_a m/s"
        return NaN
    end

    α = atan(va_b[3], va_b[1])
    β = asin(clamp(va_b[2] / v_a, -1.0, 1.0))

    gamma_prev = copy(solver.sol.gamma_distribution)

    function residual!(F, turn_rate, p)
        omega = [0.0, 0.0, turn_rate[1]]
        va_vec = [cos(α) * cos(β), sin(β), sin(α)] * v_a
        VortexStepMethod.set_va!(body_aero, va_vec, omega)
        VortexStepMethod.solve!(solver, body_aero, gamma_prev; log=false)
        F[1] = solver.sol.moment[3]
        return nothing
    end

    # Initial guess: current coupled sim turn rate
    u0 = [wing.turn_rate[3]]

    prob = NonlinearProblem(residual!, u0)
    sol = solve(prob, NewtonRaphson(; autodiff=AutoFiniteDiff());
                abstol=1e-8, reltol=1e-8)

    if sol.retcode == ReturnCode.Success
        return sol.u[1]
    else
        @warn "NonlinearSolve failed with retcode: $(sol.retcode)"
        return NaN
    end
end

"""
    run_physics_replay(csv_path; kwargs...)

Replay CSV data using physics simulation with CSV inputs.
Updates tether length and steering/depower segments from CSV data at each step.
"""
function run_physics_replay(csv_path::String;
                        start_frame=START_FRAME,
                        end_frame=END_FRAME,
                        n_substeps=5)

    raw_data = load_flight_data(csv_path)
    limited_data = limit_frames(raw_data; start_frame, end_frame)
    limited_data = add_distance_column(limited_data)

    # Interpolate CSV data for smoother control
    @info "Interpolating CSV data with $n_substeps substeps"
    csv_data = interpolate_csv_data(limited_data, n_substeps)

    @info "Loading v3 kite system structure from YAML"
    set_data_path("data/v3")
    set = Settings("system.yaml")
    set.g_earth = 9.81
    set.l_tether = 212.68
    vsm_set_path = joinpath(get_data_path(), "vsm_settings_reduced_for_coupling.yaml")
    vsm_set = VortexStepMethod.VSMSettings(vsm_set_path; data_prefix=false)
    vsm_set.wings[1].geometry_file = "data/v3/aero_geometry_stable.yaml"
    sys_struct = load_sys_struct_from_yaml("data/v3/struc_geometry_stable.yaml";
        system_name="v3", set, wing_type=SymbolicAWEModels.REFINE, vsm_set)
    csv_sys_struct = load_sys_struct_from_yaml("data/v3/struc_geometry_stable.yaml";
        system_name="v3", set, wing_type=SymbolicAWEModels.REFINE, vsm_set)
    num_wing_points = count(p -> p.type == SymbolicAWEModels.WING, sys_struct.points)
    wing_area = sys_struct.wings[1].vsm_aero.projected_area / num_wing_points
    for point in sys_struct.points
        if point.type == SymbolicAWEModels.WING
            point.area = wing_area
            point.drag_coeff = WING_CD
        end
    end
    sam = SymbolicAWEModel(set, sys_struct)
    init!(sam)

    n_steps = length(csv_data.time)
    @info "Creating log with $n_steps timesteps"
    sys_state = SysState(sam)
    logger = Logger(sam, n_steps)

    # Create CSV reference model for visualization
    csv_sam = SymbolicAWEModel(set, csv_sys_struct)
    init!(csv_sam)
    csv_state = SysState(csv_sam)
    csv_logger = Logger(csv_sam, n_steps)

    # Storage for tape percentages (for plotting)
    csv_tape_times = Float64[]
    csv_tape_steering_pct = Float64[]
    csv_tape_depower_pct = Float64[]
    phys_tape_times = Float64[]
    phys_tape_steering_pct = Float64[]
    phys_tape_depower_pct = Float64[]

    # Loop through CSV data and update sys_struct
    @info "Replaying CSV data..."
    replay_start = time()
    sys = sam.sys_struct
    SymbolicAWEModels.set_body_frame_damping(sys, INITIAL_DAMPING)

    # Calculate dt from interpolated CSV timesteps
    dt = csv_data.time[2] - csv_data.time[1]
    @info "Using timestep dt = $dt s"

    # Calculate initial delta between set.l_tether and CSV tether length
    first_csv_tether_len = csv_data.ground_tether_length[1]
    tether_len_delta = set.l_tether - first_csv_tether_len
    @info "Tether length delta" set_l_tether=set.l_tether csv_tether=first_csv_tether_len delta=tether_len_delta

    # Reset heading spike filter and distance tracker
    PREV_CSV_HEADING[] = NaN
    reset_distance_tracker!()

    function get_row(csv_data, step)
        csv_row = (
            time = csv_data.time[step],
            roll = csv_data.kite_0_roll[step],
            pitch = csv_data.kite_0_pitch[step],
            yaw = csv_data.kite_0_yaw[step],
            x = csv_data.kite_pos_east[step],
            y = csv_data.kite_pos_north[step],
            z = csv_data.kite_height[step],
            vx = csv_data.kite_est_vx[step],
            vy = csv_data.kite_est_vy[step],
            vz = csv_data.kite_est_vz[step],
            tether_len = csv_data.ground_tether_length[step],
            tether_vel = csv_data.ground_tether_reelout_speed[step],
            tether_force = csv_data.ground_tether_force[step] * 9.81,  # Convert kg to N
            steering = csv_data.kite_actual_steering[step],
            depower = csv_data.kite_actual_depower[step],
            distance = csv_data.distance[step],
            cumulative_distance = csv_data.cumulative_distance[step],
            wind_at_kite = coalesce(csv_data.lidar_wind_velocity_at_kite_mps[step], 10.0),
            angle_of_attack = deg2rad(coalesce(csv_data.airspeed_angle_of_attack[step], 0.0))
        )
    end

    try
        for step in 1:n_steps-1
            # Create row data structure
            csv_row = get_row(csv_data, step)

            # Update CSV reference model and log
            update_sys_struct_from_csv!(csv_sam.sys_struct, csv_row)
            SymbolicAWEModels.reinit!(csv_sam, csv_sam.prob, FBDF())
            update_sys_state!(csv_state, csv_sam)
            csv_state.winch_force[1] = csv_row.tether_force
            csv_state.AoA = csv_row.angle_of_attack
            csv_state.time = csv_row.time
            csv_state.l_tether[1] = csv_row.tether_len
            csv_state.v_reelout[1] = csv_row.tether_vel
            csv_state.v_wind_gnd[1] = csv_row.wind_at_kite
            log!(csv_logger, csv_state)

            # Store CSV tape percentages for plotting
            push!(csv_tape_times, csv_row.time)
            push!(csv_tape_steering_pct, csv_row.steering)
            push!(csv_tape_depower_pct, csv_row.depower)

            # Update system structure from CSV on first step
            if step == 1
                update_sys_struct_from_csv!(sam.sys_struct, csv_row)
                SymbolicAWEModels.reinit!(sam, sam.prob, FBDF())
            end

            # Update damping
            t = csv_row.time
            if t <= DECAY_TIME
                current_damping = (INITIAL_DAMPING - MIN_DAMPING) *
                                  (1.0 - t / DECAY_TIME) + MIN_DAMPING
                SymbolicAWEModels.set_body_frame_damping(sam.sys_struct, current_damping)
            end

            # Apply control and step
            brake = true
            set_value = update_vel_from_csv!(sam.sys_struct, csv_row, brake)

            # Update winch tether length and velocity from CSV
            sam.sys_struct.winches[1].tether_len = csv_row.tether_len + tether_len_delta
            sam.sys_struct.winches[1].tether_vel = csv_row.tether_vel
            SymbolicAWEModels.reinit!(sam, sam.prob, FBDF())

            next_step!(sam; dt=dt, set_values=[set_value])

            # Log state
            update_sys_state!(sys_state, sam)
            sys_state.time = t
            log!(logger, sys_state)

            # Store physics tape percentages for plotting
            push!(phys_tape_times, t)
            push!(phys_tape_steering_pct,
                steering_length_to_percentage(sam.sys_struct.segments[89].l0))
            push!(phys_tape_depower_pct,
                depower_length_to_percentage(sam.sys_struct.segments[88].l0))

            # Progress reporting
            if step % max(1, div(n_steps, 10)) == 0 || step == n_steps
                elapsed = time() - replay_start
                @info "  Step $step/$n_steps " *
                      "(t = $(round(csv_row.time, digits=2)) s)"
            end
        end
    catch err
        if err isa AssertionError
            @warn "Still plotting"
        else
            rethrow(err)
        end
    end

    replay_elapsed = time() - replay_start
    @info "CSV replay logged in $(round(replay_elapsed, digits=2)) s"

    @info "Saving logs..."
    save_log(logger, "csv_replay")
    save_log(csv_logger, "csv_reference")
    syslog = load_log("csv_replay")
    csvlog = load_log("csv_reference")

    # Create tape percentages data for plotting
    phys_tape_pct = (
        time = phys_tape_times,
        steering = phys_tape_steering_pct,
        depower = phys_tape_depower_pct
    )
    csv_tape_pct = (
        time = csv_tape_times,
        steering = csv_tape_steering_pct,
        depower = csv_tape_depower_pct
    )

    return sam, syslog, csv_sam, csvlog, csv_data, raw_data, phys_tape_pct, csv_tape_pct
end

# Main execution
sam, syslog, csv_sam, csvlog, csv_data, raw_data, phys_tape_pct, csv_tape_pct =
    run_physics_replay(CSV_PATH)

# Find equilibrium turn rate at end of simulation
@info "Finding equilibrium turn rate using VSM..."
wing = sam.sys_struct.wings[1]

# Get coupled simulation turn rate from heading derivative
dt_log = syslog.syslog.time[end] - syslog.syslog.time[end-1]
coupled_turn_rate = (syslog.syslog.heading[end] - syslog.syslog.heading[end-1]) / dt_log
coupled_turn_rate_deg = rad2deg(coupled_turn_rate)

# Get CSV turn rate from heading derivative
csv_turn_rate = (csvlog.syslog.heading[end] - csvlog.syslog.heading[end-1]) / dt_log
csv_turn_rate_deg = rad2deg(csv_turn_rate)

# Get current yaw moment before solving
body_aero = wing.vsm_aero
solver = wing.vsm_solver
current_moment_z = solver.sol.moment[3]

# Find equilibrium turn rate
equilibrium_turn_rate = find_equilibrium_turn_rate(wing)
equilibrium_turn_rate_deg = rad2deg(equilibrium_turn_rate)

# Print comparison
println("\n" * "="^60)
println("TURN RATE COMPARISON")
println("="^60)
println("Coupled sim turn rate:     $(round(coupled_turn_rate_deg, digits=4))°/s " *
        "($(round(coupled_turn_rate, digits=6)) rad/s)")
println("CSV turn rate:             $(round(csv_turn_rate_deg, digits=4))°/s " *
        "($(round(csv_turn_rate, digits=6)) rad/s)")
println("Equilibrium turn rate:     $(round(equilibrium_turn_rate_deg, digits=4))°/s " *
        "($(round(equilibrium_turn_rate, digits=6)) rad/s)")
println("Difference (eq - coupled): $(round(equilibrium_turn_rate_deg - coupled_turn_rate_deg, digits=4))°/s")
println("Current yaw moment (M_z):  $(round(current_moment_z, digits=4)) Nm")
println("="^60)

# Plot results
fig = plot([sam.sys_struct, csv_sam.sys_struct], [syslog, csvlog];
     plot_tether=true, plot_aero_force=false, plot_kite_vel=true,
     plot_wind=true, plot_reelout=false, plot_v_app=false, plot_turn_rates=true,
     tape_lengths=[phys_tape_pct, csv_tape_pct],
     suffixes=["phys", "csv"])
display(fig)
sphere = plot_sphere_trajectory([syslog,csvlog])
