# Copyright (c) 2025 Jelle Poland, Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

"""
TU Delft V3 Kite: Zenith hold followed by circular flight

This example keeps the kite at a target azimuth (zenith) and then transitions
into a circular flight pattern using the same simulated system.
"""

using SymbolicAWEModels
using VortexStepMethod
using LinearAlgebra
using Statistics
using GLMakie
using KiteUtils
using DiscretePIDs
using Dates
using StaticArrays

"""
    run_v3_kite(wing_type::WingType; kwargs...)

Run a two-phase v3 kite simulation with the specified wing type
(zenith azimuth hold, then circular flight). The two phases can have
independent durations/FPS.

# Arguments
- `wing_type::WingType`: Either `REFINE` or `QUATERNION`

# Keyword Arguments
- `sim_time_zenith::Float64=300.0`: Duration of zenith phase [s]
- `fps_zenith::Int=1`: Frames per second for zenith logging
- `sim_time_circles::Float64=300.0`: Duration of circular flight phase [s]
- `fps_circles::Int=1`: Frames per second for circular flight logging
- `remake_cache::Bool=false`: Force rebuild of cached model
- `initial_damping::Float64=10.0`: Initial world frame damping [N·s/m]
- `decay_time::Float64=2.0`: Time for damping decay [s] (applied in both phases)
- `up::Float64=0.4`: Power-tape input (0–1) mapped to segment 88 length
- `ramp_time_up::Float64=25.0`: Power-tape ramp time during zenith phase [s]
- `ramp_time_us::Float64=25.0`: Steering/power ramp time during circular phase [s]
- `start_ramp_time::Float64=0.0`: Time offset before starting power-tape ramp [s] (zenith phase)
- `max_us_zenith::Float64=0.1`: Max steering tape change for azimuth PID [m/1.4]
- `us::Float64=0.1`: Steering tape command for the circular phase
- `show_plots::Bool=false`: Display 3D plots during simulation
- `v_wind::Float64=15.4`: Wind speed [m/s]
- `v_wind_base::Float64=15.0`: Baseline wind used for circular phase ramp [m/s]
- `upwind_dir::Float64=-90.0`: Wind direction [°]
- `heading_p::Float64=0.0`: Proportional gain for azimuth controller
- `heading_i::Float64=0.1`: Integral gain for azimuth controller
- `heading_d::Float64=0.0`: Derivative gain for azimuth controller
- `target_azimuth::Float64=0.0`: Azimuth setpoint during zenith phase [rad]
- `winch_p::Float64=1000.0`: Proportional gain for winch controller [N/m]
- `winch_i::Float64=100.0`: Integral gain for winch controller [N/(m·s)]
- `winch_d::Float64=50.0`: Derivative gain for winch controller [N·s/m]
- `tube_bending_resistance::Float64=0.0`: Outward body-frame force magnitude applied to tube nodes during circular flight

# Returns
- `SysLog`: The simulation log containing time history data
"""
function run_v3_kite(wing_type::WingType;
                     sim_time_zenith=300.0,
                     fps_zenith=1,
                     sim_time_circles=300.0,
                     fps_circles=1,
                     remake_cache=false,
                     initial_damping=100.0,
                     decay_time=2.0,
                     up=0.4,
                     ramp_time_up=25.0,
                     ramp_time_us=25.0,
                     start_ramp_time=0.0,
                     max_us_zenith=0.1,
                     show_plots=false,
                     v_wind=15.4,
                     upwind_dir=-90.0,
                     heading_p=0.0, # Proportional gain
                     heading_i=0.1, # Integral gain
                     heading_d=0.0, # Derivative gain
                     winch_p=1000.0, # Proportional gain [N/m]
                     winch_i=100.0, # Integral gain [N/(m·s)]
                     winch_d=50.0, # Derivative gain [N·s/m]
                     REFINE = true,
                     target_azimuth=0.0,
                     us=0.1,
                     v_wind_base=15.0,
                     tube_bending_resistance=0.0
                     )

    wing_type_str = wing_type == SymbolicAWEModels.REFINE ? "REFINE" : "QUATERNION"
    @info "Running v3 kite simulation with $wing_type_str wing type (zenith -> circular)..."

    # Load settings
    set_data_path("data/v3")
    set = Settings("system.yaml")
    set.v_wind = v_wind
    set.upwind_dir = upwind_dir

    # Load YAML structure path
    model_name = wing_type == QUATERNION ? "v3_quat" : "v3"
    struc_yaml_path = joinpath("data", "v3", "CORRECT_struc_geometry.yaml")

    # Load VSMSettings
    vsm_set_path = joinpath(get_data_path(), "CORRECT_vsm_settings.yaml")
    vsm_set = VortexStepMethod.VSMSettings(vsm_set_path; data_prefix=false)
    vsm_set.wings[1].geometry_file = "data/v3/CORRECT_aero_geometry.yaml"

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

    # Initialize model
    n_unrefined = sys.wings[1].vsm_wing.n_unrefined_sections
    SymbolicAWEModels.init!(sam; remake=remake_cache, ignore_l0=false, remake_vsm=true)


    # Create logger (two phases with independent settings)
    n_steps_zenith = max(1, Int(round(fps_zenith * sim_time_zenith)))
    Δt_zenith = sim_time_zenith / n_steps_zenith
    n_steps_circles = max(1, Int(round(fps_circles * sim_time_circles)))
    Δt_circles = sim_time_circles / n_steps_circles
    total_steps = n_steps_zenith + n_steps_circles
    logger = Logger(sam, total_steps + 1)
    sys_state = SysState(sam)
    sys_state.time = 0.0
    log!(logger, sys_state)

    # Store nominal segment lengths for PID/control
    nominal_l0_87 = sys.segments[87].l0
    nominal_l0_88 = sys.segments[88].l0
    nominal_l0_89 = sys.segments[89].l0
    power_tape_change = ((200 + 5000 * up) / 1000) - nominal_l0_88  # Convert mm to m

    # Storage for azimuth setpoint (targeting azimuth=0)
    azimuth_setpoint = Float64[]
    push!(azimuth_setpoint, target_azimuth)

    # Create PID controller for azimuth (reuse heading gains/limits)
    max_steering = max_us_zenith * 1.4  # Max steering line length change [m]

    azimuth_pid = DiscretePID(;
        K = heading_p > 0 ? heading_p : 1.0,
        Ti = heading_i > 0 ? 1.0 / heading_i : false,
        Td = heading_d > 0 ? heading_d : false,
        Ts = Δt_zenith,
        umin = -abs(max_steering),
        umax = abs(max_steering)
    )
    @info "Azimuth PID controller initialized" target_azimuth target_azimuth_deg=rad2deg(target_azimuth) heading_p heading_i heading_d

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
        Ts = Δt_zenith,
        umin = -max_force,
        umax = max_force
    )
    @info "Winch PID controller initialized" nominal_tether_length winch_p winch_i winch_d

    # Optional initial plot
    if show_plots
        scene = plot(sam.sys_struct)
        display(scene)
    end

    # Phase 1: Zenith/azimuth hold
    @info "Starting zenith phase: $n_steps_zenith steps, Δt = $(round(Δt_zenith, digits=4)) s"
    @info "Power-tape ramp" up=round(up, digits=3) target_l0_88=round(nominal_l0_88 + power_tape_change, digits=4) ramp_time=ramp_time_up start_ramp_time=start_ramp_time
    sim_start_time = time()

    for step in 1:n_steps_zenith
        t = step * Δt_zenith

        # Update damping (decays to zero by decay_time)
        if t <= decay_time
            current_damping = initial_damping * (1.0 - t / decay_time)
            SymbolicAWEModels.set_world_frame_damping(sam.sys_struct, current_damping)
        else
            SymbolicAWEModels.set_world_frame_damping(sam.sys_struct, 0.0)
        end

        # PID azimuth control: drive azimuth to target
        current_azimuth = sam.sys_struct.wings[1].azimuth
        steering_control = azimuth_pid(target_azimuth, current_azimuth, 0.0)

        push!(azimuth_setpoint, target_azimuth)

        # Ramp power-tape change on segment 88
        ramp_factor = ramp_time_up > 0 ? clamp((t - start_ramp_time) / ramp_time_up, 0.0, 1.0) : 1.0
        power_control = power_tape_change * ramp_factor
        sys.segments[88].l0 = nominal_l0_88 + power_control

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
            next_step!(sam; set_values=[-winch_torque], dt=Δt_zenith, vsm_interval=1)
        catch err
            if err isa AssertionError
                @error "next_step! failed during zenith phase" step
                break
            end
            rethrow(err)
        end

        # Log state
        update_sys_state!(sys_state, sam)
        sys_state.time = t
        log!(logger, sys_state)

        # Progress updates
        if step % max(1, div(n_steps_zenith, 10)) == 0 || step == n_steps_zenith
            elapsed = time() - sim_start_time
            times_realtime = t / elapsed
            @info "  Zenith step $step/$n_steps_zenith (t = $(round(t, digits=2)) s)" times_realtime=round(times_realtime, digits=2)
        end
    end

    # Phase 2: Circular flight
    @info "Switching to circular flight phase" phase_time=round(sys_state.time, digits=2)
    SymbolicAWEModels.set_world_frame_damping(sam.sys_struct, initial_damping)

    # Fix tether length for circular flight
    winch.brake = true
    winch.set_value = 0.0

    steering_tape_change = 1400 * us / 1000  # Convert mm to m
    vw_change = v_wind - v_wind_base

    power_target = nominal_l0_88 + power_tape_change
    steer_target_87 = nominal_l0_87 + steering_tape_change
    steer_target_89 = nominal_l0_89 - steering_tape_change
    power_start = sys.segments[88].l0
    steer_start_87 = sys.segments[87].l0
    steer_start_89 = sys.segments[89].l0

    for step in 1:n_steps_circles
        t_stage = step * Δt_circles
        t_total = sim_time_zenith + t_stage

        # Update damping (reset for this phase)
        if t_stage <= decay_time
            current_damping = initial_damping * (1.0 - t_stage / decay_time)
            SymbolicAWEModels.set_world_frame_damping(sam.sys_struct, current_damping)
        else
            SymbolicAWEModels.set_world_frame_damping(sam.sys_struct, 0.0)
        end

        ramp_factor = min(t_stage / ramp_time_us, 1.0)

        # Power and steering ramp toward circular flight targets
        sys.segments[88].l0 = power_start + (power_target - power_start) * ramp_factor
        sys.segments[87].l0 = steer_start_87 + (steer_target_87 - steer_start_87) * ramp_factor
        sys.segments[89].l0 = steer_start_89 + (steer_target_89 - steer_start_89) * ramp_factor

        # Update wind speed (kept constant here, ramp hook retained)
        sam.sys_struct.set.v_wind = v_wind_base + vw_change

        # Apply outward tube bending resistance forces in body-frame ±y directions
        if tube_bending_resistance != 0
            R_b_w = sam.sys_struct.wings[1].R_b_w
            force_pos_y = R_b_w * [0.0, tube_bending_resistance, 0.0]
            force_neg_y = R_b_w * [0.0, -tube_bending_resistance, 0.0]
            sam.sys_struct.points[2].disturb .= force_pos_y
            sam.sys_struct.points[3].disturb .= force_pos_y
            sam.sys_struct.points[20].disturb .= force_neg_y
            sam.sys_struct.points[21].disturb .= force_neg_y
        end

        # Fixed tether length: brake engaged; no winch torque
        winch_torque = 0.0
        sys.winches[1].set_value = -winch_torque

        # Advance simulation
        try
            next_step!(sam; set_values=[-winch_torque], dt=Δt_circles, vsm_interval=1)
        catch err
            if err isa AssertionError
                @error "next_step! failed during circular phase" step
                break
            end
            rethrow(err)
        end

        # Log state
        update_sys_state!(sys_state, sam)
        sys_state.time = t_total
        log!(logger, sys_state)

        # Progress updates
        if step % max(1, div(n_steps_circles, 10)) == 0 || step == n_steps_circles
            elapsed = time() - sim_start_time
            times_realtime = t_total / elapsed
            @info "  Circle step $step/$n_steps_circles (t = $(round(t_total, digits=2)) s)" times_realtime=round(times_realtime, digits=2)
        end
    end

    # Calculate performance
    total_wall_time = time() - sim_start_time
    total_sim_time = sim_time_zenith + sim_time_circles
    final_times_realtime = total_sim_time / total_wall_time
    @info "Simulation completed: $wing_type_str (zenith + circular)" wall_time=round(total_wall_time, digits=2) times_realtime=round(final_times_realtime, digits=2)

    # Save and load log
    log_name = "tmp_run_$(lowercase(wing_type_str))"
    save_log(logger, log_name)
    syslog = load_log(log_name)

    # Permanent save
    save_dir = joinpath("processed_data", "v3_kite")
    isdir(save_dir) || mkpath(save_dir)
    timestamp = Dates.format(Dates.now(), "yyyy_mm_dd_HH_MM")
    up_tag = Int(round(up*100))
    us_tag = Int(round(us*100))
    v_wind_tag = Int(round(v_wind))
    log_name = "zenith_circle__up_$(up_tag)_us_$(us_tag)_vw_$(v_wind_tag)_date_$(timestamp)"
    save_log(logger, log_name; path=save_dir)

    return syslog, sam, azimuth_setpoint
end

# ============= Main Execution =============

# Run both simulations
syslog_refine, sam_refine, azimuth_setpoint_refine = run_v3_kite(SymbolicAWEModels.REFINE;
    # general settings
    v_wind=8,
    v_wind_base=15,
    up=0.2,
    # settings zenith initialisation flight
    sim_time_zenith=100, 
    fps_zenith=60,
    start_ramp_time=0.1,
    ramp_time_up=10.0,
    initial_damping=100.0,
    decay_time=2.0,
    max_us_zenith = 0.02,
    target_azimuth = 0.0,
    # settings circular flight
    sim_time_circles=300,
    fps_circles=60,
    ramp_time_us = 15.0,
    us=0.2,
)


#### Plot results#
fig = plot(sam_refine.sys_struct, syslog_refine;
    plot_turn_rates=true,
    plot_reelout=false,
    plot_gk=true,
    plot_aoa=true,
    plot_heading=false,
    plot_elevation=true,
    plot_azimuth=true,
    plot_winch_force=false,
    plot_set_values=false,
    plot_us=true,)

scene = replay(syslog_refine, sam_refine.sys_struct)

scr1 = display(fig)
wait(scr1)
scr2 = display(scene)
wait(scr2)

nothing
