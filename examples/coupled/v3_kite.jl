# Copyright (c) 2025 Jelle Poland, Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

"""
TU Delft V3 Kite Comparison: REFINE vs QUATERNION Wing Types

This example compares REFINE and QUATERNION wing models:
- REFINE: Panel forces applied to structural points, deformable wing
- QUATERNION: Rigid body dynamics with group twist DOFs

The script runs both simulations and plots them for direct comparison.
"""

using SymbolicAWEModels
using VortexStepMethod
using LinearAlgebra
using Statistics
using GLMakie
using KiteUtils
using DiscretePIDs

REFINE = true
QUAT = false

# Heading PID controller parameters
MAX_HEADING = 10.0  # Maximum heading amplitude [degrees]
PERIOD = 60.0       # Oscillation period [seconds]
HEADING_P = 0.0     # Proportional gain
HEADING_I = 0.1     # Integral gain
HEADING_D = 0.0     # Derivative gain

# Tether winch PID controller parameters
WINCH_P = 1000.0    # Proportional gain [N/m]
WINCH_I = 100.0     # Integral gain [N/(m·s)]
WINCH_D = 50.0      # Derivative gain [N·s/m]

# V3 Kite steering/depower calibration (from KCU documentation)
STEERING_L0 = 1.6  # Neutral steering tape length (m)
STEERING_GAIN = 1.2  # Maximum differential (m) at |u_s| = 1

# Depower calibration
DEPOWER_L0 = 0.2
DEPOWER_GAIN = 5.0

"""
    steering_percentage_to_lengths(percentage)

Convert steering percentage to left/right tape lengths (m).
Percentage convention: negative = left turn, positive = right turn.
"""
function steering_percentage_to_lengths(percentage)
    u_s = percentage / 100.0  # Convert percentage to [-1, 1]
    L_left = STEERING_L0 - STEERING_GAIN * u_s
    L_right = STEERING_L0 + STEERING_GAIN * u_s
    return L_left, L_right
end

"""
    depower_percentage_to_length(percentage)

Convert depower percentage to tape length (m).
"""
function depower_percentage_to_length(percentage)
    u_p = percentage / 100.0  # Convert percentage to [0, 1]
    L_depower = DEPOWER_L0 + DEPOWER_GAIN * u_p
    return L_depower
end

"""
    run_v3_kite(wing_type::WingType; kwargs...)

Run a v3 kite simulation with the specified wing type.

# Arguments
- `wing_type::WingType`: Either `REFINE` or `QUATERNION`

# Keyword Arguments
- `sim_time::Float64=300.0`: Simulation duration [s]
- `fps::Int=1`: Frames per second for logging
- `remake_cache::Bool=false`: Force rebuild of cached model
- `initial_damping::Float64=10.0`: Initial world frame damping [N·s/m]
- `decay_time::Float64=2.0`: Time for damping decay [s]
- `max_steering::Float64=0.2`: Maximum steering line length change [m]
- `show_plots::Bool=false`: Display 3D plots during simulation
- `v_wind::Float64=15.4`: Wind speed [m/s]
- `upwind_dir::Float64=-90.0`: Wind direction [°]
- `max_heading::Float64=50.0`: Maximum heading amplitude for sine wave [°]
- `period::Float64=20.0`: Oscillation period for sine wave [s]
- `heading_p::Float64=0.0`: Proportional gain for heading controller
- `heading_i::Float64=0.1`: Integral gain for heading controller
- `heading_d::Float64=0.0`: Derivative gain for heading controller
- `winch_p::Float64=1000.0`: Proportional gain for winch controller [N/m]
- `winch_i::Float64=100.0`: Integral gain for winch controller [N/(m·s)]
- `winch_d::Float64=50.0`: Derivative gain for winch controller [N·s/m]

# Returns
- `SysLog`: The simulation log containing time history data
"""
function run_v3_kite(wing_type::WingType;
                     sim_time=300.0,
                     fps=1,
                     remake_cache=false,
                     initial_damping=100.0,
                     decay_time=2.0,
                     max_steering=0.1,
                     show_plots=false,
                     v_wind=15.4,
                     upwind_dir=-90.0,
                     max_heading=50.0,
                     period=20.0,
                     heading_p=0.0,
                     heading_i=0.1,
                     heading_d=0.0,
                     winch_p=1000.0,
                     winch_i=100.0,
                     winch_d=50.0)

    wing_type_str = wing_type == SymbolicAWEModels.REFINE ? "REFINE" : "QUATERNION"
    @info "Running v3 kite simulation with $wing_type_str wing type..."

    # Load settings
    set_data_path("data/v3")
    set = Settings("system.yaml")
    set.v_wind = v_wind
    set.upwind_dir = upwind_dir

    # Load YAML structure path
    model_name = wing_type == QUATERNION ? "v3_quat" : "v3"
    struc_yaml_path = joinpath("data", "v3", "struc_geometry_stable.yaml")

    # Load VSMSettings
    vsm_set_path = joinpath(get_data_path(), "vsm_settings_reduced_for_coupling.yaml")
    vsm_set = VortexStepMethod.VSMSettings(vsm_set_path; data_prefix=false)

    # Use 36 panels for both wing types (matches vsm_settings.yaml default)
    vsm_set.wings[1].n_panels = 36
    # Note: n_unrefined_sections is automatically inferred from YAML geometry

    # Load system structure with wing_type and vsm_set parameters
    sys = load_sys_struct_from_yaml(struc_yaml_path;
        system_name=model_name, set, wing_type, vsm_set)


    # Initialize damping
    SymbolicAWEModels.set_world_frame_damping(sys, initial_damping)

    wing_points = [p for p in sys.points if p.type == WING]
    n_unrefined = sys.wings[1].vsm_wing.n_unrefined_sections
    @info "$wing_type_str wing setup:" n_wing_points=length(wing_points) n_groups=length(sys.groups) n_unrefined=n_unrefined n_panels=length(sys.wings[1].vsm_aero.panels) n_segments=length(sys.segments)

    # Create symbolic model
    sam = SymbolicAWEModel(set, sys)

    # Apply steering
    # sys.segments[87].l0 += max_steering
    # sys.segments[89].l0 -= max_steering

    # Initialize model
    n_unrefined = sys.wings[1].vsm_wing.n_unrefined_sections
    SymbolicAWEModels.init!(sam; remake=remake_cache, ignore_l0=false, remake_vsm=true)

    # # Stabilization phase
    # @info "Stabilizing system..."
    # [point.fix_static = true for point in sys.points if point.type == WING]
    # if wing_type == QUATERNION
    #     sys.wings[1].fix_sphere = true
    # end
    # @time next_step!(sam; dt=10.0)
    # [point.fix_static = false for point in sys.points if point.type == WING]
    # if wing_type == QUATERNION
    #     sys.wings[1].fix_sphere = false
    # end

    # Create logger
    n_steps = Int(round(fps * sim_time))
    Δt = sim_time / n_steps
    logger = Logger(sam, n_steps + 1)
    sys_state = SysState(sam)
    sys_state.time = 0.0
    log!(logger, sys_state)

    # Store nominal segment lengths for PID control
    nominal_l0_87 = sys.segments[87].l0
    nominal_l0_89 = sys.segments[89].l0

    # Storage for heading setpoint
    heading_setpoint = Float64[]
    push!(heading_setpoint, 0.0)  # Initial setpoint

    # Create PID controller for heading
    max_heading_rad = deg2rad(max_heading)
    angular_freq = 2π / period  # rad/s
    heading_pid = DiscretePID(;
        K = heading_p > 0 ? heading_p : 1.0,
        Ti = heading_i > 0 ? 1.0 / heading_i : false,
        Td = heading_d > 0 ? heading_d : false,
        Ts = Δt,
        umin = -abs(max_steering),
        umax = abs(max_steering)
    )
    @info "Heading PID controller initialized" max_heading period heading_p heading_i heading_d
    @info "  Sine wave: ±$(max_heading)°, period=$(period)s"

    # Store nominal tether length for winch controller
    nominal_tether_length = sys.winches[1].tether_len

    # Initialize winch torque based on initial tether force
    winch = sys.winches[1]
    initial_force = norm(winch.force)
    initial_torque = -winch.drum_radius / winch.gear_ratio * initial_force + winch.friction
    winch.set_value = initial_torque
    @info "Winch initialized" initial_force initial_torque drum_radius=winch.drum_radius gear_ratio=winch.gear_ratio friction=winch.friction

    # Create PID controller for tether winch (maintain constant length)
    # Output is force [N], will be converted to torque in control loop
    max_force = 50000.0  # N
    winch_pid = DiscretePID(;
        K = winch_p,
        Ti = winch_i > 0 ? winch_p / winch_i : false,
        Td = winch_d > 0 ? winch_d / winch_p : false,
        Ts = Δt,
        umin = -max_force,
        umax = max_force
    )
    @info "Winch PID controller initialized" nominal_tether_length winch_p winch_i winch_d

    # Optional initial plot
    if show_plots
        scene = plot(sam.sys_struct)
        display(scene)
    end

    # Time-marching loop
    @info "Starting simulation: $n_steps steps, Δt = $(round(Δt, digits=4)) s"
    sim_start_time = time()

    for step in 1:n_steps
        t = step * Δt

        # Update damping
        if t <= decay_time
            current_damping = initial_damping * (1.0 - t / decay_time)
            SymbolicAWEModels.set_world_frame_damping(sam.sys_struct, current_damping)
        else
            SymbolicAWEModels.set_world_frame_damping(sam.sys_struct, 0.0)
        end

        # PID heading control with sine wave setpoint
        target_heading_rad = max_heading_rad * sin(angular_freq * t)
        current_heading = sam.sys_struct.wings[1].heading
        steering_control = heading_pid(target_heading_rad, current_heading, 0.0)

        # Store setpoint for plotting
        push!(heading_setpoint, target_heading_rad)

        # Apply differential steering (opposite signs for turning moment)
        sys.segments[87].l0 = nominal_l0_87 + steering_control
        sys.segments[89].l0 = nominal_l0_89 - steering_control

        # PID winch control to maintain constant tether length
        current_tether_length = sys.winches[1].tether_len
        winch_force_control = winch_pid(nominal_tether_length, current_tether_length, 0.0)

        # Convert force to torque: τ = -r/G * F + friction
        winch_torque = -winch.drum_radius / winch.gear_ratio * winch_force_control + winch.friction
        sys.winches[1].set_value = -winch_torque

        # Advance simulation
        try
            next_step!(sam; set_values=[-winch_torque], dt=Δt, vsm_interval=1)
        catch err
            if err isa AssertionError
                @error "next_step! failed" step
                break
            end
            rethrow(err)
        end

        # Log state
        update_sys_state!(sys_state, sam)
        sys_state.time = t
        log!(logger, sys_state)

        # Progress updates
        if step % max(1, div(n_steps, 10)) == 0 || step == n_steps
            elapsed = time() - sim_start_time
            times_realtime = t / elapsed
            @info "  Step $step/$n_steps (t = $(round(t, digits=2)) s)" times_realtime=round(times_realtime, digits=2)
        end
    end

    # Calculate performance
    total_wall_time = time() - sim_start_time
    final_times_realtime = sim_time / total_wall_time
    @info "Simulation completed: $wing_type_str" wall_time=round(total_wall_time, digits=2) times_realtime=round(final_times_realtime, digits=2)

    # Save and load log
    log_name = "tmp_run_$(lowercase(wing_type_str))"
    save_log(logger, log_name)
    syslog = load_log(log_name)

    return syslog, sam, heading_setpoint
end

# ============= Main Execution =============

@info "V3 Kite Comparison: Running REFINE and QUATERNION simulations..."

# Run both simulations
heading_setpoint_refine = nothing
heading_setpoint_quat = nothing

if REFINE
    syslog_refine, sam_refine, heading_setpoint_refine = run_v3_kite(SymbolicAWEModels.REFINE;
        sim_time=200.0, fps=60, show_plots=false,
        max_heading=MAX_HEADING, period=PERIOD,
        heading_p=HEADING_P, heading_i=HEADING_I, heading_d=HEADING_D,
        winch_p=WINCH_P, winch_i=WINCH_I, winch_d=WINCH_D)
end
if QUAT
    syslog_quat, sam_quat, heading_setpoint_quat = run_v3_kite(QUATERNION;
        sim_time=50.0, fps=24, show_plots=false,
        max_heading=MAX_HEADING, period=PERIOD,
        heading_p=HEADING_P, heading_i=HEADING_I, heading_d=HEADING_D,
        winch_p=WINCH_P, winch_i=WINCH_I, winch_d=WINCH_D)
end

@info "Both simulations complete. Creating comparison plots..."

fig = nothing
# Create comparison plot
if REFINE && QUAT
    fig = plot([sam_refine.sys_struct, sam_quat.sys_struct], [syslog_refine, syslog_quat];
               plot_turn_rates=true,
               plot_azimuth=true,
               plot_elevation=true,
               plot_aoa=true,
               plot_heading=true,
               plot_default=false,
               plot_aero_force=true,
               plot_tether=true,
               heading_setpoint=[heading_setpoint_refine, heading_setpoint_quat])
end

if QUAT && !REFINE
    fig = plot(sam_quat.sys_struct, syslog_quat;
               plot_aero_moment=true,
               plot_tether=true,
               heading_setpoint=heading_setpoint_quat)
end
if !QUAT && REFINE
    fig = plot(sam_refine.sys_struct, syslog_refine;
               plot_tether=true,
               heading_setpoint=heading_setpoint_refine)
end
display(fig)

@info "Comparison plot created!"

nothing
