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
using Dates
using StaticArrays

# --- Helper: geometric AoA of mid panel using panel axes (robust, no corner ordering) ---
mid_panel_index(n::Integer) = cld(n, 2)  # middle panel index

# AoA from panel body-frame axes: project apparent wind into chord-normal plane
function aoa_of_panel_body(panel, va_b::AbstractVector{<:Real})
    v = collect(va_b)
    vnorm = norm(v)
    vnorm > 1e-9 || return 0.0
    v̂ = v ./ vnorm
    x̂ = panel.x_airf              # chord direction (body frame)
    ŷ = panel.y_airf              # span direction (body frame)
    ẑ = panel.z_airf              # normal direction (body frame)
    v2 = v̂ .- dot(v̂, ŷ) .* ŷ   # remove span component
    v2_norm = norm(v2)
    v2_norm > 1e-9 || return 0.0
    v2̂ = v2 ./ v2_norm
    return atan(dot(v2̂, ẑ), dot(v2̂, x̂))
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
- `up::Float64=0.4`: Power-tape input (0–1) mapped to segment 88 length
- `ramp_time::Float64=25.0`: Time to ramp power-tape change [s]
- `start_ramp_time::Float64=0.0`: Time offset before starting power-tape ramp [s]
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
                     up=0.4,
                     ramp_time=25.0,
                     start_ramp_time=0.0,
                     max_us=0.1,
                     show_plots=false,
                     v_wind=15.4,
                     upwind_dir=-90.0,
                     max_heading=0.0, # Maximum heading amplitude degrees
                     period=60.0, # Oscillation period seconds 
                     heading_p=0.0, # Proportional gain
                     heading_i=0.1, # Integral gain
                     heading_d=0.0, # Derivative gain
                     winch_p=1000.0, # Proportional gain [N/m]
                     winch_i=100.0, # Integral gain [N/(m·s)]
                     winch_d=50.0, # Derivative gain [N·s/m]
                     REFINE = true,
                     target_azimuth=-90
                     )

    wing_type_str = wing_type == SymbolicAWEModels.REFINE ? "REFINE" : "QUATERNION"
    @info "Running v3 kite simulation with $wing_type_str wing type..."

    # Load settings
    set_data_path("data/v3")
    set = Settings("system.yaml")
    set.v_wind = v_wind
    set.upwind_dir = upwind_dir

    # Load YAML structure path
    model_name = wing_type == QUATERNION ? "v3_quat" : "v3"
    struc_yaml_path = joinpath("data", "v3", "struc_geometry.yaml")

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
    nominal_l0_88 = sys.segments[88].l0
    nominal_l0_89 = sys.segments[89].l0
    power_tape_change = ((200 + 5000 * up) / 1000) - nominal_l0_88  # Convert mm to m

    # Storage for azimuth setpoint (targeting azimuth=0)
    azimuth_setpoint = Float64[]
    push!(azimuth_setpoint, 0.0)

    # Create PID controller for azimuth (reuse heading gains/limits)
    max_steering = max_us * 1.4  # Max steering line length change [m]

    azimuth_pid = DiscretePID(;
        K = heading_p > 0 ? heading_p : 1.0,
        Ti = heading_i > 0 ? 1.0 / heading_i : false,
        Td = heading_d > 0 ? heading_d : false,
        Ts = Δt,
        umin = -abs(max_steering),
        umax = abs(max_steering)
    )
    @info "Azimuth PID controller initialized" target_azimuth_deg=0.0 heading_p heading_i heading_d

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
    @info "Power-tape ramp" up=round(up, digits=3) target_l0_88=round(nominal_l0_88 + power_tape_change, digits=4) ramp_time=ramp_time start_ramp_time=start_ramp_time
    sim_start_time = time()

    for step in 1:n_steps
        t = step * Δt

        # Update damping
        if t <= decay_time
            # Linear decay: reaches zero exactly at decay_time
            current_damping = initial_damping * (1.0 - t / decay_time)
            SymbolicAWEModels.set_world_frame_damping(sam.sys_struct, current_damping)
            # @info "  Damping update" step t=round(t, digits=2) current_damping=round(current_damping, digits=2)
        else
            SymbolicAWEModels.set_world_frame_damping(sam.sys_struct, 0.0)
        end

        # PID azimuth control: drive azimuth to 0
        current_azimuth = sam.sys_struct.wings[1].azimuth
        steering_control = azimuth_pid(target_azimuth, current_azimuth, 0.0)

        # Store setpoint for plotting/inspection
        push!(azimuth_setpoint, target_azimuth)

        # Ramp power-tape change on segment 88
        ramp_factor = ramp_time > 0 ? clamp((t - start_ramp_time) / ramp_time, 0.0, 1.0) : 1.0
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

    # Permanent save
    save_dir = joinpath("processed_data", "v3_kite")
    isdir(save_dir) || mkpath(save_dir)
    timestamp = Dates.format(Dates.now(), "yyyy_mm_dd_HH_MM")
    up_tag = Int(round(up*100))
    us = (nominal_l0_87 - sys.segments[87].l0)/1.4  # Steering input derived from segment length change
    us_tag = round(Int, us*100)
    if length(string(us_tag)) == 1
        us_tag = "0$(us_tag)"
    end
    v_wind_tag = Int(round(v_wind))
    log_name = "zenith__up_$(up_tag)_us_$(us_tag)_vw_$(v_wind_tag)_date_$(timestamp)"
    save_log(logger, log_name; path=save_dir)

    return syslog, sam, azimuth_setpoint
end

# ============= Main Execution =============

# Run both simulations
syslog_refine, sam_refine, azimuth_setpoint_refine = run_v3_kite(SymbolicAWEModels.REFINE;
    sim_time=50, 
    fps=240,
    v_wind=15,
    up=0.4, 
    start_ramp_time=0.1,
    ramp_time=10.0,
    initial_damping=100.0,
    decay_time=2.0,
    max_us = 0.02,
    target_azimuth = 0
)

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


#### Plot results#
fig = plot(sam_refine.sys_struct, syslog_refine;
    plot_turn_rates=true, plot_reelout=false, plot_gk=true,
    plot_aoa=true, plot_heading=false, plot_elevation=true,
    plot_azimuth=true, plot_winch_force=false, plot_set_values=false,
    plot_us=true,)

scene = replay(syslog_refine, sam_refine.sys_struct)

scr1 = display(fig)
wait(scr1)
scr2 = display(scene)
wait(scr2)

nothing
