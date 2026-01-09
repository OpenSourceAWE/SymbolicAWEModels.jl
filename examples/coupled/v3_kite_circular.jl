# Copyright (c) 2025 Jelle Poland, Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

"""
TU Delft V3 Kite: REFINE Wing

This example runs a REFINE wing model and plots the results.
"""

using SymbolicAWEModels
using VortexStepMethod
using LinearAlgebra
using Statistics
using GLMakie
using KiteUtils
using DiscretePIDs
using Dates
using StaticArrays: SVector

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

"""
    run_v3_kite(; kwargs...)

Run a v3 kite simulation using the REFINE wing type.

# Keyword Arguments
- `sim_time::Float64=300.0`: Simulation duration [s]
- `fps::Int=1`: Frames per second for logging
- `remake_cache::Bool=false`: Force rebuild of cached model
- `initial_damping::Float64=10.0`: Initial body frame damping [N·s/m]
- `damping_pattern::Vector{Float64}=[0.0, 1.0, 1.0]`: Per-axis damping pattern [x, y, z]
- `decay_time::Float64=2.0`: Time for damping decay [s]
- `min_damping::Float64=0.0`: Minimum damping after decay [N·s/m]
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
- `tube_bending_resistance::Float64=0.0`: Outward body-frame force magnitude applied to nodes 2, 3 (+y) and 20, 21 (-y)
- `save_subdir::AbstractString=""`: Subfolder under `processed_data/v3_kite` for permanent saves
- `run_tag::AbstractString=""`: Extra tag appended to the log name

# Returns
- `SysLog`: The simulation log containing time history data
"""
function run_v3_kite(;
                     sim_time=300.0,
                     fps=4,
                     remake_cache=false,
                     initial_damping=100.0,
                     damping_pattern=[0.0, 1.0, 1.0],
                     decay_time=10.0,
                     min_damping=0.0,
                     up = 0.4,
                     us=0.1,
                     show_plots=false,
                     v_wind=15.4,
                     upwind_dir=-90.0,
                     ramp_time=25.0,
                     v_wind_base=15,
                     tether_length=150.0,
                     max_heading=50.0,
                     period=20.0,
                     heading_p=0.0,
                     heading_i=0.1,
                     heading_d=0.0,
                     winch_p=1000.0,
                     winch_i=100.0,
                     winch_d=50.0,
                     tube_bending_resistance=0.0,
                     save_subdir="",
                     run_tag="")

    wing_type = SymbolicAWEModels.REFINE
    wing_type_str = "REFINE"
    @info "Running v3 kite simulation in n_steps: $(Int(round(fps * sim_time)))"

    # Load settings
    set_data_path("data/v3")
    set = Settings("system.yaml")
    set.v_wind = v_wind
    set.upwind_dir = upwind_dir
    # set.l_tethers[1] = tether_length
    set.v_reel_outs[1] = 0.0

    # Load YAML structure path based on depower and tether length
    model_name = "v3"
    depower_int = Int(round(up * 100))
    tether_int = Int(round(tether_length))
    struc_yaml_path = joinpath("data", "v3",
        "struc_geometry_depower$(depower_int)_tether$(tether_int).yaml")
    aero_yaml_path = "data/v3/aero_geometry_depower$(depower_int)_tether$(tether_int).yaml"

    # Load VSMSettings
    vsm_set_path = joinpath(get_data_path(), "CORRECT_vsm_settings.yaml")
    vsm_set = VortexStepMethod.VSMSettings(vsm_set_path; data_prefix=false)
    vsm_set.wings[1].geometry_file = aero_yaml_path

    # Use 36 panels for both wing types (matches vsm_settings.yaml default)
    vsm_set.wings[1].n_panels = 36
    # Note: n_unrefined_sections is automatically inferred from YAML geometry

    # Load system structure with wing_type and vsm_set parameters
    sys = load_sys_struct_from_yaml(struc_yaml_path;
        system_name=model_name, set, wing_type, vsm_set)

    # Initialize damping with per-axis values [x, y, z]
    SymbolicAWEModels.set_body_frame_damping(sys, damping_pattern * initial_damping)

    wing_points = [p for p in sys.points if p.type == WING]
    n_unrefined = sys.wings[1].vsm_wing.n_unrefined_sections
    @info "REFINE wing setup:" n_wing_points=length(wing_points) n_groups=length(sys.groups) n_unrefined=n_unrefined n_panels=length(sys.wings[1].vsm_aero.panels) n_segments=length(sys.segments)

    # Create symbolic model
    sam = SymbolicAWEModel(set, sys)

    # Apply steering
    # sys.segments[87].l0 += max_steering
    # sys.segments[89].l0 -= max_steering

    # Initialize model
    n_unrefined = sys.wings[1].vsm_wing.n_unrefined_sections
    sam.set.l_tether = tether_length
    sam.sys_struct.winches[1].tether_len = tether_length
    sam.sys_struct.winches[1].brake = true
    SymbolicAWEModels.init!(sam; remake=remake_cache, ignore_l0=false, remake_vsm=true)
    # SymbolicAWEModels.reinit!(sam, sam.prob, SymbolicAWEModels.FBDF())

    # Create logger
    n_steps = Int(round(fps * sim_time))
    Δt = sim_time / n_steps
    logger = Logger(sam, n_steps + 1)
    sys_state = SysState(sam)
    sys_state.time = 0.0
    log!(logger, sys_state)

    # Store nominal segment lengths for PID control
    nominal_l0_87 = sys.segments[87].l0
    nominal_l0_88 = sys.segments[88].l0
    nominal_l0_89 = sys.segments[89].l0

    # Storage for heading setpoint
    heading_setpoint = Float64[]
    push!(heading_setpoint, 0.0)  # Fixed steering: keep heading setpoint at zero
    # # Create PID controller for heading
    # max_heading_rad = deg2rad(max_heading)
    # angular_freq = 2π / period  # rad/s
    # heading_pid = DiscretePID(;
    #     K = heading_p > 0 ? heading_p : 1.0,
    #     Ti = heading_i > 0 ? 1.0 / heading_i : false,
    #     Td = heading_d > 0 ? heading_d : false,
    #     Ts = Δt,
    #     umin = -abs(max_steering),
    #     umax = abs(max_steering)
    # )
    # @info "Heading PID controller initialized" max_heading period heading_p heading_i heading_d
    # @info "  Sine wave: ±$(max_heading)°, period=$(period)s"

    # Keep winch torque at zero; brake is engaged to fix tether length
    winch = sys.winches[1]
    winch.set_value = 0.0

    # Optional initial plot
    if show_plots
        scene = plot(sam.sys_struct)
        display(scene)
    end

    ## ACTUATION
    #len_power-tape = 200 + 5000 * up
    #len_steering_tape = 1600 + 1400 * us
    #us < 0 shortens right tape, and lengths left tape, causing a right turn
    #RIGHT len_steering_tape (us = -1) = 1600 - 1400 = 200
    #LEFT len_steering_tape (us = 1) = 1600 + 1400 = 3000

    steering_tape_change = 1400 * us / 1000  # Convert mm to m
    power_tape_change = ((200 + 5000 * up) / 1000) - nominal_l0_88  # Convert mm to m
    vw_change = v_wind - v_wind_base


    # Time-marching loop
    @info "Starting simulation: $n_steps steps, Δt = $(round(Δt, digits=4)) s"
    @info " Initial lengths (m): segment 87: $(round(nominal_l0_87, digits=4)), segment 88: $(round(nominal_l0_88, digits=4)), segment 89: $(round(nominal_l0_89, digits=4))"
    @info " Steering tape change (m): $(round(steering_tape_change, digits=4)), Power tape change (m): $(round(power_tape_change, digits=4))"
    sim_start_time = time()
    aoa_log_interval_steps = max(1, Int(round(3.0 / Δt)))  # roughly every 3 seconds

    # Storage for segment stretch statistics (after t > 1.0)
    max_stretch_samples = Float64[]
    mean_stretch_samples = Float64[]
    max_idx_samples = Int[]

    wings = sam.sys_struct.wings
    wing = wings[1]

    for step in 1:n_steps
        t = step * Δt

        # Update body frame damping: decay to min_damping
        current_damping = max(initial_damping * (1.0 - t / decay_time), min_damping)
        SymbolicAWEModels.set_body_frame_damping(sam.sys_struct, damping_pattern * current_damping)

        # World frame damping: decay from 1000 to 0 over first 1.0s
        world_damping = t <= 1.0 ? 1000.0 * (1.0 - t) : 0.0
        SymbolicAWEModels.set_world_frame_damping(sam.sys_struct, world_damping)

        # Fixed tether length: brake engaged; only steering ramp is applied
        ramp_factor = min(t / ramp_time, 1.0)
        steering_control = steering_tape_change * ramp_factor
        power_control = power_tape_change * ramp_factor
        push!(heading_setpoint, 0.0)  # Keep heading setpoint flat for plotting

        # Apply power control
        sys.segments[88].l0 = nominal_l0_88 + power_control
        # print every x second & time below 30 sec
        if step % Int(round(3.0 / Δt)) == 0 && t <= 30.0 && power_tape_change > 1e-4
            @info "power-tape = $(round(sys.segments[88].l0, digits=4)) m at t=$(round(t, digits=2)) s"
        end

        # Apply differential steering (opposite signs for turning moment)
        sys.segments[87].l0 = nominal_l0_87 + steering_control
        sys.segments[89].l0 = nominal_l0_89 - steering_control

        # Update wind speed linearly
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

        # Convert force to torque: τ = -r/G * F + friction
        winch_torque = 0.0
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

        # Collect segment stretch statistics after t > 1.0
        if t > 1.0
            max_stretch, mean_stretch, max_idx = segment_stretch_stats(sam.sys_struct)
            push!(max_stretch_samples, max_stretch)
            push!(mean_stretch_samples, mean_stretch)
            push!(max_idx_samples, max_idx)
        end

        # Log AoA periodically (value stored in sys_state by update_sys_state!)
        # if step % aoa_log_interval_steps == 0
        #     alpha = sys_state.AoA
        #     alpha = atan(wing.va_b[3], wing.va_b[1])
        #     vsm_alpha = wing.vsm_solver.sol.alpha_dist[length(wing.vsm_solver.sol.alpha_dist) ÷ 2 + (length(wing.vsm_solver.sol.alpha_dist) % 2)]
        #     @info "---> Angle of attack" t=round(t, digits=2) sys_state.AoA=round(rad2deg(alpha), digits=2) vsm_alpha=round(rad2deg(vsm_alpha), digits=2)

        # end

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

    # Report segment stretch statistics (for t > 1.0)
    if !isempty(max_stretch_samples)
        overall_max_stretch = maximum(max_stretch_samples)
        overall_mean_stretch = mean(mean_stretch_samples)
        max_stretch_idx = max_idx_samples[argmax(max_stretch_samples)]
        @info "Segment stretch statistics (t > 1.0):" max_relative=round(overall_max_stretch, digits=6) max_percentage=round(overall_max_stretch*100, digits=4) mean_relative=round(overall_mean_stretch, digits=6) mean_percentage=round(overall_mean_stretch*100, digits=4) max_segment_idx=max_stretch_idx
    end

    lt_tag = Int(round(tether_length))

    # Save and load log
    log_name = "tmp_run_$(lowercase(wing_type_str))_lt_$(lt_tag)"
    save_log(logger, log_name)
    syslog = load_log(log_name)

    # Permanent save
    save_root = joinpath("processed_data", "v3_kite")
    save_dir = isempty(save_subdir) ? save_root : joinpath(save_root, save_subdir)
    isdir(save_dir) || mkpath(save_dir)
    timestamp = Dates.format(Dates.now(), "yyyy_mm_dd_HH_MM_SS")
    up_tag = Int(round(up*100))
    us_tag = Int(round(us*100))
    v_wind_tag = Int(round(v_wind))
    log_name = "circle__up_$(up_tag)" * "_" * "us_$(us_tag)" * "_" * "vw_$(v_wind_tag)" * "_" * "lt_$(lt_tag)"
    if !isempty(run_tag)
        log_name *= "_" * run_tag
    end
    log_name *= "_date_" * timestamp
    save_log(logger, log_name; path=save_dir)


    return syslog, sam, heading_setpoint
end

# ==========================================
# ============= Main Execution =============
# ==========================================
    us = 0.15  # {{{ 0.0  <> 0.30 }}} suitable range ~kite-as-a-sensor
    up = 0.4 # 22  # {{{ 0.4 <> 0.5 }}} 0.5858 is baseline ~PIM's thesis
    #0.4151powered and #0.5012depowered #0.39 during turns
    vw = 11.1  # {{{ 10.  <> 15.0 }}} suitable range?
    lt = 225  # problems when changing...

    sim_time = 5.0
    decay_time = 10.0 #2secs works better than 3 somehow
    ramp_time = 10.0 #2sec
    fps = 600
    initial_damping = 10.0
    damping_pattern = [100.0, 100.0, 100.0]
    min_damping = 1.0
    tube_bending_resistance = 0  # N


    syslog_refine, sam_refine, heading_setpoint_refine = run_v3_kite(
        sim_time=sim_time, fps=fps,
        up=up, us=us, v_wind=vw, tether_length=lt,
        decay_time=decay_time, ramp_time=ramp_time,
        initial_damping=initial_damping, damping_pattern=damping_pattern, min_damping=min_damping,
        max_heading=MAX_HEADING, period=PERIOD,
        tube_bending_resistance=tube_bending_resistance,
        heading_p=HEADING_P, heading_i=HEADING_I, heading_d=HEADING_D,
        winch_p=WINCH_P, winch_i=WINCH_I, winch_d=WINCH_D)


    fig = plot(sam_refine.sys_struct,
               syslog_refine;
               plot_turn_rates=false,
               plot_reelout=false,
               plot_twist=false,
               plot_yaw_rate_paper=true,
               plot_v_app=true,
               plot_kite_vel=true,
               plot_gk=true,
               plot_aoa=true, 
               plot_heading=false, 
               plot_elevation=true,
               plot_azimuth=true, 
               plot_winch_force=false, 
               plot_set_values=false,
               gk_ylims=(0.0, 15.0), 
               aoa_ylims=(0.0, 15.0),
               yaw_rate_paper_ylims=(0.0, 50.0),
               plot_tether_actual=true,
               plot_us=true)

scene = replay(syslog_refine, sam_refine.sys_struct)

scr1 = display(fig)
wait(scr1)
scr2 = display(scene)
wait(scr2)



# # Report final geometric AoA using hardcoded mid-panel corners (world frame)
# last_state = syslog_refine.syslog[end]
# X = last_state.X; Y = last_state.Y; Z = last_state.Z
# # Mid-panel corners: 10,11,12,13 (11/13 front; 10/12 back)
# back = 0.5 .* ([X[10], Y[10], Z[10]] .+ [X[12], Y[12], Z[12]])
# front = 0.5 .* ([X[11], Y[11], Z[11]] .+ [X[13], Y[13], Z[13]])

# delta_z = front[3] - back[3]
# delta_x = front[1] - back[1]
# aoa_wrt_horizontal = -rad2deg(atan(delta_z, delta_x))
# # @info "alpha wrt horizontal $(round(aoa_wrt_horizontal, digits=2))"

# chord_w = front .- back
# wing = sam_refine.sys_struct.wings[1]
# v_app_w = wing.R_b_w * wing.va_b
# @info "v_app" v_app_w=round.(v_app_w, digits=2)
# aoa_geom_deg = rad2deg(acos(clamp(dot(chord_w, v_app_w) / (norm(chord_w) * norm(v_app_w) + 1e-12), -1.0, 1.0)))
# @info "alpha wrt v_app $(round(aoa_geom_deg, digits=2))"


# ##########################
# ### compute L/D system ###
# ##########################
# sl = syslog_refine.syslog
# last_state = sl[end]
# prev_state = sl[end - 1]
# dt = (last_state.time - prev_state.time) + 1e-12

# # compute wing v_a
# min1 = sl[end - 1]
# last_state = sl[end]

# X_last = last_state.X; Y_last = last_state.Y; Z_last = last_state.Z
# X_min1 = min1.X; Y_min1 = min1.Y; Z_min1 = min1.Z

# X_last_back = 0.5 * (X_last[10] + X_last[12])
# Y_last_back = 0.5 * (Y_last[10] + Y_last[12])
# Z_last_back = 0.5 * (Z_last[10] + Z_last[12])

# X_last_front = 0.5 * (X_last[11] + X_last[13])
# Y_last_front = 0.5 * (Y_last[11] + Y_last[13])
# Z_last_front = 0.5 * (Z_last[11] + Z_last[13])

# X_min1_back = 0.5 * (X_min1[10] + X_min1[12])
# Y_min1_back = 0.5 * (Y_min1[10] + Y_min1[12])
# Z_min1_back = 0.5 * (Z_min1[10] + Z_min1[12])

# X_min1_front = 0.5 * (X_min1[11] + X_min1[13])
# Y_min1_front = 0.5 * (Y_min1[11] + Y_min1[13])
# Z_min1_front = 0.5 * (Z_min1[11] + Z_min1[13])

# va_mid_panel_front = SVector{3,Float64}(
#     (X_last_front - X_min1_front) / (dt) - last_state.v_wind_kite[1],
#     (Y_last_front - Y_min1_front) / (dt) - last_state.v_wind_kite[2],
#     (Z_last_front - Z_min1_front) / (dt) - last_state.v_wind_kite[3],
# )
# va_mid_panel_back = SVector{3,Float64}(
#     (X_last_back - X_min1_back) / (dt) - last_state.v_wind_kite[1],
#     (Y_last_back - Y_min1_back) / (dt) - last_state.v_wind_kite[2],
#     (Z_last_back - Z_min1_back) / (dt) - last_state.v_wind_kite[3],
# )
# va_mid_panel = -0.5 .* (va_mid_panel_front .+ va_mid_panel_back)
# va_mid_panel_unit = va_mid_panel / (norm(va_mid_panel) + 1e-12)
# # @info "v_app mid-panel $(round.(va_mid_panel, digits=5))"

# # Aero forces in world frame
# R_b_w = SymbolicAWEModels.quaternion_to_rotation_matrix(last_state.orient)
# F_aero_b = last_state.aero_force_b
# F_aero_world = R_b_w * F_aero_b

# drag_dir = va_mid_panel_unit              # drag is positive aligned with v_a
# lift_dir = cross(drag_dir, SVector(0.0, 1.0, 0.0))  # lift is positive perpendicular to drag and upwards
# # @info "checking lift dir" lift_dir=round.(lift_dir, digits=4) drag_dir=round.(drag_dir, digits=4)

# drag_wing = dot(F_aero_world, drag_dir)
# lift_wing = dot(F_aero_world, lift_dir)

# # Recompute tether drag from segment states (using the last two snapshots for velocity)
# segments = sam_refine.sys_struct.segments
# points = sam_refine.sys_struct.points
# n_points = length(last_state.X)
# pos_last = [SVector(last_state.X[i], last_state.Y[i], last_state.Z[i]) for i in 1:n_points]
# pos_prev = [SVector(prev_state.X[i], prev_state.Y[i], prev_state.Z[i]) for i in 1:n_points]
# vel_est = [(pos_last[i] - pos_prev[i]) ./ dt for i in 1:n_points]

# wind_vec_gnd = last_state.v_wind_gnd
# cd_tether = SymbolicAWEModels.get_cd_tether(sam_refine.set)
# # @info "cd_tether $(cd_tether)"

# let drag_bridles = 0.0, drag_tether = 0.0, lift_bridles = 0.0, lift_tether = 0.0
#     for segment in segments
#         p1, p2 = segment.point_idxs
#         p1_pos = pos_last[p1]; p2_pos = pos_last[p2]
#         # Skip structural wing segments (drag handled by VSM)
#         if points[p1].type == SymbolicAWEModels.WING && points[p2].type == SymbolicAWEModels.WING
#             continue
#         end
#         seg_vec = p2_pos - p1_pos
#         seg_len = norm(seg_vec) + 1e-12
#         seg_dir = seg_vec / seg_len

#         seg_vel = 0.5 .* (vel_est[p1] + vel_est[p2])
#         seg_height = max(0.0, 0.5 * (p1_pos[3] + p2_pos[3]))
#         wind_factor = SymbolicAWEModels.calc_wind_factor(sam_refine.am, max(seg_height, 1.0), sam_refine.set)
#         wind_vel = wind_factor .* wind_vec_gnd
#         va_seg = wind_vel - seg_vel
#         app_perp = va_seg .- dot(va_seg, seg_dir) * seg_dir

#         area = seg_len * segment.diameter
#         rho = SymbolicAWEModels.calc_rho(sam_refine.am, seg_height)
#         v_perp_mag = norm(app_perp)
#         Tether_force = 0.5 * rho * cd_tether * area * v_perp_mag .* app_perp

#         drag_scalar = dot(Tether_force, drag_dir)
#         lift_scalar = dot(Tether_force, lift_dir)
        
#         if 47 <= segment.idx <= 89
#             drag_bridles += drag_scalar
#             lift_bridles += lift_scalar
#         elseif 90 <= segment.idx <= 95
#             drag_tether += drag_scalar
#             lift_tether += lift_scalar
#         end
        
#     end

#     # Total aero (wing + tether drag approximation)
#     total_drag = drag_wing + drag_bridles + drag_tether
#     total_lift = lift_wing + lift_tether + lift_bridles
#     @info "L/D wing" lift_wing = round(lift_wing, digits=2) drag_wing = round(drag_wing, digits=2) L_over_D = round(lift_wing / (drag_wing + 1e-12), digits=2)
#     @info "Bridle aero" drag_bridles = round(drag_bridles, digits=2) lift_bridles = round(lift_bridles, digits=2)
#     @info "Tether aero" drag_tether = round(drag_tether, digits=2) lift_tether = round(lift_tether, digits=2)
#     @info "L/D system" lift_total = round(total_lift, digits=2) drag_total = round(total_drag, digits=2) L_over_D = round(total_lift / (total_drag + 1e-12), digits=2)
# end

nothing
