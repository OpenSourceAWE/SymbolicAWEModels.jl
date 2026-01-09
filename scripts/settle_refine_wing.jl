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
using GLMakie
using CSV, DataFrames

# Configuration
WORLD_DAMPING = 1000.0  # Ns/m
DECAY_STEPS = 5000     # Steps over which damping decays to zero
NUM_STEPS = 2000
DT = 0.1  # seconds
STEERING_PERCENTAGE = 0.0  # steering [-100, 100]
DEPOWER_PERCENTAGE = 40   # depower [0, 100]
SOURCE_STRUC_PATH = "data/v3/CORRECT_struc_geometry.yaml"
DEST_STRUC_PATH = "data/v3/struc_geometry_stable_$(DEPOWER_PERCENTAGE).yaml"
SOURCE_AERO_PATH = "data/v3/CORRECT_aero_geometry.yaml"
DEST_AERO_PATH = "data/v3/aero_geometry_stable_$(DEPOWER_PERCENTAGE).yaml"
WIND_VEL = 20.0
ELEVATION = 75
TETHER_LENGTH = 212.68  # Total tether length (m), 6 segments
EXTRA_POINTS_CSV = "data/v3/straight_flight_reelout_frame_7182.csv"
LE_FRAC = 0.95  # Factor to reduce l0 of LE struts (segments 20-28)
TIP_REDUCTION = 0.4

# V3 Kite steering/depower calibration (from KCU documentation)
STEERING_L0 = 0.8  # Neutral steering tape length (m)
STEERING_GAIN = 1.4  # Maximum differential (m) at |u_s| = 1
DEPOWER_L0 = 0.0
DEPOWER_GAIN = 5.0

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

"""
    load_extra_points(csv_path, sys_struct)

Load extra points from CSV and transform from camera frame to simulation frame.
"""
function load_extra_points(csv_path::String, sys_struct; body_offset=[0.3, 0.0, 0.2])
    df = CSV.read(csv_path, DataFrame)

    # CSV reference: LE[10], LE[11] (0-indexed in CSV, so Julia indices 11, 12)
    le_pts = [[r.x, r.y, r.z] for r in eachrow(df) if r.group == "LE"]
    csv_le10, csv_le11 = le_pts[11], le_pts[12]
    csv_le_center = (csv_le10 + csv_le11) / 2

    # CSV strut centers: strut3[1]/strut4[1] are at TE, [end] are at LE
    strut3 = [[r.x, r.y, r.z] for r in eachrow(df) if r.group == "strut3"]
    strut4 = [[r.x, r.y, r.z] for r in eachrow(df) if r.group == "strut4"]
    csv_le_center = (strut3[end] + strut4[end]) / 2
    csv_te_center = (strut3[1] + strut4[1]) / 2

    # Sim reference: points 10, 12 (center LE), point 27 (bridle)
    sim_p10 = collect(sys_struct.points[10].pos_w)
    sim_p12 = collect(sys_struct.points[12].pos_w)
    sim_le_center = (sim_p10 + sim_p12) / 2
    point_27 = collect(sys_struct.points[27].pos_w)
    cam_pos = point_27 + sys_struct.wings[1].R_b_w * [0, 0.2, 0]

    # Direction vectors
    csv_span = normalize(strut4[end] - strut3[end])

    # CSV basis: y=spanwise, z from wing center geometry, x from cross
    csv_y = csv_span
    csv_wing_center = (csv_le_center + csv_te_center) / 2
    @show csv_wing_center
    csv_z = normalize(csv_wing_center - csv_y * 0.84/2)
    csv_x = cross(csv_y, csv_z)

    # Sim basis: directly from wing rotation matrix
    R_b_w = sys_struct.wings[1].R_b_w
    sim_x = R_b_w[:, 1]
    sim_y = R_b_w[:, 2]
    sim_z = R_b_w[:, 3]

    # Rotation: R * csv_basis = sim_basis
    csv_basis = hcat(csv_x, csv_y, csv_z)
    sim_basis = hcat(sim_x, sim_y, sim_z)
    R = sim_basis * csv_basis'

    # Translation: align LE centers
    T = sim_le_center - R * csv_le_center + R_b_w * body_offset

    # Transform all points (including camera origin marker)
    all_pts = [[row.x, row.y, row.z] for row in eachrow(df)]
    push!(all_pts, zeros(3))
    transformed = [Tuple(R * p + T) for p in all_pts]

    return transformed
end

@info "Settling REFINE wing with world frame damping..."
@info "Configuration" WORLD_DAMPING DECAY_STEPS NUM_STEPS DT total_time=NUM_STEPS*DT

# Load settings
set_data_path("data/v3")
set = Settings("system.yaml")
set.g_earth = 9.81
set.v_wind = WIND_VEL
set.l_tether = TETHER_LENGTH

# Load VSMSettings
vsm_set_path = joinpath(get_data_path(), "vsm_settings_reduced_for_coupling.yaml")
vsm_set = VortexStepMethod.VSMSettings(vsm_set_path; data_prefix=false)
vsm_set.wings[1].geometry_file = "data/v3/aero_geometry.yaml"
vsm_set.wings[1].n_panels = 36

# Load system structure with REFINE wing type
struc_yaml_path = joinpath("data", "v3", "CORRECT_struc_geometry.yaml")
sys = load_sys_struct_from_yaml(struc_yaml_path;
    system_name="v3", set,
    wing_type=SymbolicAWEModels.REFINE, vsm_set)
sys.transforms[1].elevation = deg2rad(ELEVATION)

# Update tether length: points 39-44 and segments 90-95
segment_len = TETHER_LENGTH / 6.0
for (i, point_idx) in enumerate(39:44)
    sys.points[point_idx].pos_cad .= [0.0, 0.0, -i * segment_len]
end
for seg_idx in 90:95
    sys.segments[seg_idx].l0 = segment_len
end
@info "Tether configured" TETHER_LENGTH segment_len

# Apply LE strut l0 reduction factor (segments 20-28)
if LE_FRAC != 1.0
    for seg_idx in 20:28
        sys.segments[seg_idx].l0 *= LE_FRAC
    end
    @info "LE struts l0 reduced" LE_FRAC
end

sys.segments[47].l0 -= TIP_REDUCTION
sys.segments[48].l0 -= TIP_REDUCTION
sys.segments[57].l0 -= TIP_REDUCTION
sys.segments[58].l0 -= TIP_REDUCTION

# Set initial world frame damping (will decay over DECAY_STEPS)
SymbolicAWEModels.set_world_frame_damping(sys, WORLD_DAMPING)

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

# Create logger
logger = Logger(sam, NUM_STEPS + 1)
sys_state = SysState(sam)
sys_state.time = 0.0
log!(logger, sys_state)

# Simulation loop with repositioning every step
@info "Starting settling simulation..."
for step in 1:NUM_STEPS
    t = step * DT

    # Decay world damping linearly over DECAY_STEPS
    if step <= DECAY_STEPS
        damping = WORLD_DAMPING * (1.0 - step / DECAY_STEPS)
        SymbolicAWEModels.set_world_frame_damping(sys, damping)
    end

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
        current_damping = step <= DECAY_STEPS ?
            WORLD_DAMPING * (1.0 - step / DECAY_STEPS) : 0.0
        @info "Step $step/$NUM_STEPS (t = $(round(t, digits=2)) s)" damping=round(current_damping, digits=1) elevation=round(rad2deg(wing.elevation), digits=2) azimuth=round(rad2deg(wing.azimuth), digits=2) heading=round(rad2deg(wing.heading), digits=2)
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

@info "Creating plots..."

# Load extra points for comparison
extra_pts = load_extra_points(EXTRA_POINTS_CSV, sam.sys_struct)

# Create 3D plot with extra points
scene = plot(sam.sys_struct; extra_points=extra_pts, body_frame=true)
display(scene)

@info "Plot created! Settling complete."
