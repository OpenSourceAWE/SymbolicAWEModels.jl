#!/usr/bin/env julia
# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

"""
Settle REFINE wing with world frame damping

This script loads the V3 kite with REFINE wing type, runs it with constant
world frame damping while maintaining alignment through repositioning, and
updates the YAML point positions to the settled equilibrium state.

The script creates a backup of the original YAML file before modifying it.
"""

using SymbolicAWEModels
using VortexStepMethod
using LinearAlgebra
using KiteUtils
using OrdinaryDiffEqCore
using CairoMakie, GLMakie
using CSV, DataFrames
using UnPack
using Dates
using DiscretePIDs

include("../examples/coupled/utils.jl")

# Configuration - Simulation
WORLD_DAMPING = 1000.0  # Ns/m
NUM_STEPS = 4000
MIN_DAMPING = 10.0  # Ns/m
DT = 0.1  # seconds
STEERING_PERCENTAGE = 0.0  # steering [-100, 100]
DEPOWER_PERCENTAGE = 40   # depower [0, 100]
WIND_VEL = 10.72
ELEVATION = 60
TETHER_LENGTH = 240  # Total tether length (m), 6 segments
EXTRA_POINTS_CSV = "data/v3/straight_flight_reelout_frame_7182.csv"
EXTRA_POINTS_FRAME = 7182
FLIGHT_CSV = "data/v3/v3_2025-10-09-ekf.csv"

# Video frame mapping (video_frame 7182 = UTC 15:36:31.0)
VIDEO_FRAME_REF = 7182
UTC_REF_SECONDS = 15*3600 + 36*60 + 31.0
VIDEO_FPS = 29.97

# Geometry reductions
TE_FRAC = 0.95       # Factor for TE wires (segments 20-28), 1.0 = no change
TIP_REDUCTION = 0.4  # Tip LE reduction (m), 0.0 = no change
L0_REDUCTION = 0.2   # Reduction for steering/depower L0 (m), 0.0 = no change

# Build destination filename suffix
DEST_SUFFIX = "l0$(L0_REDUCTION)_tip$(TIP_REDUCTION)_te$(TE_FRAC)"

SOURCE_STRUC_PATH = "data/v3/CORRECT_struc_geometry.yaml"
DEST_STRUC_PATH = "data/v3/struc_geometry_$(DEST_SUFFIX).yaml"
SOURCE_AERO_PATH = "data/v3/CORRECT_aero_geometry.yaml"
DEST_AERO_PATH = "data/v3/aero_geometry_$(DEST_SUFFIX).yaml"

# V3 Kite steering/depower calibration (from KCU documentation)
STEERING_L0 = 1.6 - L0_REDUCTION  # Neutral steering tape length (m)
STEERING_GAIN = 1.4  # Maximum differential (m) at |u_s| = 1
DEPOWER_L0 = 0.2 - L0_REDUCTION   # Neutral depower tape length (m)
DEPOWER_GAIN = 5.0

# Heading controller parameters
HEADING_KP = 0.5
HEADING_TAU_I = false  # No integral
HEADING_SETPOINT = -1.562  # Target heading in radians

"""
    steering_percentage_to_lengths(percentage)

Convert steering percentage to left/right tape lengths (m).
Percentage convention: negative = left turn, positive = right turn.
"""
function steering_percentage_to_lengths(percentage)
    u_s = percentage / 100.0  # Convert percentage to [-1, 1]
    L_left = STEERING_L0 - (STEERING_GAIN / 2.0) * u_s
    L_right = STEERING_L0 + (STEERING_GAIN / 2.0) * u_s
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

"""Convert unix timestamp to UTC seconds since midnight."""
function unix_to_utc_seconds(unix_timestamp::Float64)
    dt = Dates.unix2datetime(unix_timestamp)
    return Dates.hour(dt)*3600 + Dates.minute(dt)*60 +
           Dates.second(dt) + Dates.millisecond(dt)/1000
end

"""Convert UTC seconds to video frame number."""
function utc_to_video_frame(utc_seconds::Float64)
    delta_seconds = utc_seconds - UTC_REF_SECONDS
    return round(Int, VIDEO_FRAME_REF + delta_seconds * VIDEO_FPS)
end

@info "Settling REFINE wing with world frame damping..."
@info "Configuration" WORLD_DAMPING NUM_STEPS DT total_time=NUM_STEPS*DT

# Load settings
set_data_path("data/v3")
set = Settings("system.yaml")
set.g_earth = 9.81
set.v_wind = WIND_VEL
set.l_tether = TETHER_LENGTH

# Load VSMSettings
vsm_set_path = joinpath(get_data_path(), "vsm_settings_reduced_for_coupling.yaml")
vsm_set = VortexStepMethod.VSMSettings(vsm_set_path; data_prefix=false)
vsm_set.wings[1].geometry_file = SOURCE_AERO_PATH
vsm_set.wings[1].n_panels = 36

# Load system structure with REFINE wing type
struc_yaml_path = SOURCE_STRUC_PATH
sys = load_sys_struct_from_yaml(struc_yaml_path;
    system_name="v3", set,
    wing_type=SymbolicAWEModels.REFINE, vsm_set)
sys.transforms[1].elevation = deg2rad(ELEVATION * cos(HEADING_SETPOINT))
sys.transforms[1].azimuth = deg2rad(ELEVATION * sin(HEADING_SETPOINT))
sys.transforms[1].heading = HEADING_SETPOINT

# Update tether length: points 39-44 and segments 90-95
segment_len = TETHER_LENGTH / 6.0
for (i, point_idx) in enumerate(39:44)
    sys.points[point_idx].pos_cad .= [0.0, 0.0, -i * segment_len]
end
for seg_idx in 90:95
    sys.segments[seg_idx].l0 = segment_len
end
@info "Tether configured" TETHER_LENGTH segment_len

# Apply geometry reductions
sys.segments[47].l0 -= TIP_REDUCTION
sys.segments[48].l0 -= TIP_REDUCTION
sys.segments[57].l0 -= TIP_REDUCTION
sys.segments[58].l0 -= TIP_REDUCTION
for seg_idx in 20:28
    sys.segments[seg_idx].l0 *= TE_FRAC
end
@info "Reductions applied" TIP_REDUCTION TE_FRAC

# Set initial world frame damping (decays linearly to MIN_DAMPING)
SymbolicAWEModels.set_world_frame_damping(sys, WORLD_DAMPING)
SymbolicAWEModels.set_body_frame_damping(sys, 300.0)

wing_points = [p for p in sys.points if p.type == WING]
@info "System setup" n_wing_points=length(wing_points) n_points=length(sys.points) n_segments=length(sys.segments)

# Create symbolic model
sam = SymbolicAWEModel(set, sys)

# Initialize model
SymbolicAWEModels.init!(sam; remake=false, ignore_l0=false, remake_vsm=true)

# Enable winch brakes to lock tether length
for winch in sys.winches
    winch.brake = true
end
@info "Winch brakes enabled"

# Initialize steering and depower from configuration
L_left, L_right = steering_percentage_to_lengths(STEERING_PERCENTAGE)
L_depower = depower_percentage_to_length(DEPOWER_PERCENTAGE)

sys.segments[87].l0 = L_left   # Left steering tape
sys.segments[89].l0 = L_right  # Right steering tape
sys.segments[88].l0 = L_depower  # Depower tape

@info "Steering/depower initialized" steering=STEERING_PERCENTAGE depower=DEPOWER_PERCENTAGE L_left L_right L_depower

# Initialize heading PID controller
heading_pid = DiscretePID(;
    K = HEADING_KP,
    Ti = HEADING_TAU_I,
    Td = false,
    Ts = DT,
    umin = -1.0,
    umax = 1.0)

# Create logger
logger = Logger(sam, NUM_STEPS + 1)
sys_state = SysState(sam)
sys_state.time = 0.0
log!(logger, sys_state)

# Simulation loop with repositioning every step
@info "Starting settling simulation..."
for step in 1:NUM_STEPS
    t = step * DT

    # Linear damping decrease over entire simulation
    damping = max(WORLD_DAMPING * (1.0 - step / NUM_STEPS), MIN_DAMPING)
    SymbolicAWEModels.set_world_frame_damping(sys, damping)

    # Heading control
    wing = sys.wings[1]
    wing.R_b_w = SymbolicAWEModels.calc_refine_wing_frame(
        sys.points, wing.z_ref_points, wing.y_ref_points, wing.origin_idx)[1]
    curr_heading = calc_heading(sys, wing.R_b_w)
    delta_heading = -wrap_to_pi(HEADING_SETPOINT - curr_heading)
    steering_control = DiscretePIDs.calculate_control!(
        heading_pid, 0.0, delta_heading, 0.0)
    L_left, L_right = steering_percentage_to_lengths(STEERING_PERCENTAGE)
    sys.segments[87].l0 = L_left + STEERING_GAIN * steering_control
    sys.segments[89].l0 = L_right - STEERING_GAIN * steering_control

    # Advance one timestep
    try
        next_step!(sam; dt=DT, vsm_interval=1)
    catch err
        if err isa AssertionError
            @error "Simulation failed" step t
            break
        end
        rethrow(err)
    end

    # Log state
    update_sys_state!(sys_state, sam)
    sys_state.time = t
    log!(logger, sys_state)

    # Progress updates
    if step % 20 == 0 || step == NUM_STEPS
        wing = sys.wings[1]
        @info "Step $step/$NUM_STEPS (t = $(round(t, digits=2)) s)" damping=round(damping, digits=1) elevation=round(rad2deg(wing.elevation), digits=2) azimuth=round(rad2deg(wing.azimuth), digits=2) heading=round(rad2deg(wing.heading), digits=2)
    end
end

@info "Simulation completed. Updating YAML with settled positions..."

# Update YAML files with settled positions
update_yaml_from_sys_struct!(sys, SOURCE_STRUC_PATH, DEST_STRUC_PATH,
                            SOURCE_AERO_PATH, DEST_AERO_PATH)

# Save and load the log
log_name = "settle_refine_wing"
save_log(logger, log_name)
syslog = load_log(log_name)

@info "Settling complete." DEST_STRUC_PATH DEST_AERO_PATH

# Load extra points for interactive use
extra_pts, extra_groups = load_extra_points(EXTRA_POINTS_CSV, sam.sys_struct)

# Load flight CSV and find setpoint at EXTRA_POINTS_FRAME
@info "Loading flight CSV for setpoints..."
flight_df = CSV.read(FLIGHT_CSV, DataFrame; delim=' ', silencewarnings=true,
                     normalizenames=true, types=Float64, strict=false)

# Find row closest to EXTRA_POINTS_FRAME
frame_idx = nothing
min_diff = Inf
for (i, unix_t) in enumerate(flight_df.unix_time)
    ismissing(unix_t) && continue
    frame = utc_to_video_frame(unix_to_utc_seconds(unix_t))
    diff = abs(frame - EXTRA_POINTS_FRAME)
    if diff < min_diff
        global min_diff = diff
        global frame_idx = i
    end
end
if min_diff > 1
    @warn "Closest frame is $min_diff frames away from target $EXTRA_POINTS_FRAME"
end

# Extract setpoints from the matching row
if !isnothing(frame_idx)
    n_steps = length(syslog.syslog.time)
    # Compute heading from EKF Euler angles
    roll = flight_df.ekf_kite_roll[frame_idx]
    pitch = flight_df.ekf_kite_pitch[frame_idx]
    yaw = flight_df.ekf_kite_yaw[frame_idx]
    heading_val = calc_csv_heading(roll, pitch, yaw, sam.sys_struct)
    # Compute v_kite from EKF velocity components
    vx = flight_df.ekf_kite_velocity_x[frame_idx]
    vy = flight_df.ekf_kite_velocity_y[frame_idx]
    vz = flight_df.ekf_kite_velocity_z[frame_idx]
    v_kite_val = sqrt(vx^2 + vy^2 + vz^2)
    # Get v_app from EKF
    v_app_val = flight_df.ekf_kite_apparent_windspeed[frame_idx]
    # Get angle of attack from EKF
    aoa_val = deg2rad(flight_df.ekf_wing_angle_of_attack[frame_idx])
    # Get tether length and force
    l_tether_val = flight_df.ekf_tether_length[frame_idx]
    winch_force_val = flight_df.ground_tether_force[frame_idx]
    setpoints = Dict(
        :heading => fill(heading_val, n_steps),
        :l_tether => fill(l_tether_val, n_steps),
        :winch_force => fill(winch_force_val, n_steps),
        :v_kite => fill(v_kite_val, n_steps),
        :v_app => fill(v_app_val, n_steps),
        :aoa => fill(aoa_val, n_steps),
    )
    @info "Setpoints from frame $EXTRA_POINTS_FRAME" heading=rad2deg(heading_val) l_tether=l_tether_val winch_force=winch_force_val v_kite=v_kite_val v_app=v_app_val AoA=rad2deg(aoa_val)
else
    @warn "Could not find frame $EXTRA_POINTS_FRAME in flight CSV"
    setpoints = nothing
end

# Save plots as PDFs
CairoMakie.activate!()
for dir in (:front, :side, :top)
    scene = plot_body_frame_local(sam.sys_struct;
        extra_points=extra_pts, extra_groups=extra_groups, dir)
    pdf_filename = "data/v3/body_frame_$(dir)_settle_frame$(EXTRA_POINTS_FRAME)_$(DEST_SUFFIX).pdf"
    save(pdf_filename, scene)
    @info "Plot saved" pdf_filename
end

# Define ylims for plots
plot_ylims = Dict(
    :winch_force => (0, 2000),
    :heading => (-180, 180),
    :aoa => (0, 25),
)

# Save 2D plot
fig = plot([sam.sys_struct], [syslog];
     plot_tether=false, plot_aero_force=false, plot_kite_vel=true,
     plot_wind=false, plot_reelout=false, plot_v_app=true, plot_turn_rates=false,
     plot_winch_force=true, plot_heading=true, plot_course=false, plot_aoa=true,
     setpoints, ylims=plot_ylims, suffixes=["sim"])
pdf_2d = "data/v3/settle_2d_frame$(EXTRA_POINTS_FRAME)_$(DEST_SUFFIX).pdf"
save(pdf_2d, fig)
@info "Plot saved" pdf_2d

# Display with GLMakie
GLMakie.activate!()
fig = plot([sam.sys_struct], [syslog];
     plot_tether=false, plot_aero_force=false, plot_kite_vel=true,
     plot_wind=false, plot_reelout=false, plot_v_app=true, plot_turn_rates=false,
     plot_winch_force=true, plot_heading=true, plot_course=false, plot_aoa=true,
     setpoints, ylims=plot_ylims, suffixes=["sim"])
