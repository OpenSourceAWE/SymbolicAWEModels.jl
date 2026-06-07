# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# KPS4 parking flight comparison: SymbolicAWEModels vs KiteModels
#
# Runs the same 30s parking flight with both packages and plots
# a side-by-side comparison using the Makie extension.

using Pkg
if Base.active_project() != joinpath(@__DIR__, "Project.toml")
    Pkg.activate(@__DIR__)
end

using GLMakie
using SymbolicAWEModels
using SymbolicAWEModels: Point
using KiteModels
using KitePodModels
using KiteUtils: init!, next_step!, update_sys_state!
using LinearAlgebra

pkg_root = dirname(@__DIR__)
set_data_path(joinpath(pkg_root, "data", "kps4"))

# ==================== SHARED PARAMETERS =================== #
SIM_TIME = 30.0
dt = 0.05
N_STEPS = round(Int, SIM_TIME / dt)
UPWIND_DIR = -pi/2 + deg2rad(10)
PRE_STRESS = 0.9975

# ==================== KITEMODELS SIMULATION ================ #
println("=" ^ 60)
println("Running KiteModels simulation...")
println("=" ^ 60)

km_set = Settings("system.yaml")
km_set.abs_tol = 0.0006
km_set.rel_tol = 0.00001
km_set.solver = "DFBDF"
km_set.upwind_dir = rad2deg(UPWIND_DIR)

kcu_km = KCU(km_set)
kps4 = KPS4(kcu_km)

km_logger = Logger(km_set.segments + 5, N_STEPS)

integrator = KiteModels.init!(
    kps4; delta=0.001, stiffness_factor=0.1,
    prn=false, steady_state=false)
next_step!(kps4, integrator;
    set_speed=0, upwind_dir=UPWIND_DIR, dt=1e-10)

println("KiteModels initial state:")
println("  aoa=[$(round(kps4.alpha_2, digits=2)), " *
    "$(round(kps4.alpha_3, digits=2)), " *
    "$(round(kps4.alpha_4, digits=2))]")

km_elapsed = 0.0
for step in 1:N_STEPS
    global km_elapsed += @elapsed next_step!(
        kps4, integrator;
        set_speed=0, upwind_dir=UPWIND_DIR, dt=dt)
    local sys_state = SysState(kps4)
    sys_state.var_01 = kps4.pitch
    sys_state.var_02 = kps4.pitch_rate
    log!(km_logger, sys_state)
end
km_rt = (N_STEPS * dt) / km_elapsed
println("KiteModels realtime factor: $(round(km_rt, digits=2))")

lift, drag = KiteModels.lift_drag(kps4)
println("KiteModels final lift, drag [N]: " *
    "$(round(lift, digits=2)), $(round(drag, digits=2))")

save_log(km_logger, "kps4_km")
km_syslog = load_log("kps4_km")

# ==================== SYMBOLICAWEMODELS SIMULATION ========= #
println()
println("=" ^ 60)
println("Running SymbolicAWEModels simulation...")
println("=" ^ 60)

set = Settings("system.yaml")
set.upwind_dir = rad2deg(UPWIND_DIR)

# Geometry from KiteUtils
particles = KiteUtils.get_particles(
    set.height_k, set.h_bridle,
    set.width, set.m_k)
pos_kcu = particles[2]
pos_nose = particles[3]
pos_top = particles[4]
pos_right = particles[5]
pos_left = particles[6]

# Mass distribution
kite_mass = set.mass
k_nose = set.rel_nose_mass * kite_mass
k_top = set.rel_top_mass *
        (1.0 - set.rel_nose_mass) * kite_mass
k_side = 0.5 * (1.0 - set.rel_top_mass) *
         (1.0 - set.rel_nose_mass) * kite_mass
set.mass = 0.0

points = Point[]
push!(points, Point(:ground, zeros(3), STATIC))
push!(points, Point(:kcu, pos_kcu, DYNAMIC;
    extra_mass=set.kcu_mass,
    transform=:main_tf))
push!(points, Point(:nose, pos_nose, DYNAMIC;
    extra_mass=k_nose, transform=:main_tf))
push!(points, Point(:top, pos_top, WING;
    extra_mass=k_top, wing=:plate_wing,
    transform=:kite_tilt))
push!(points, Point(:right, pos_right, WING;
    extra_mass=k_side, wing=:plate_wing,
    transform=:kite_tilt))
push!(points, Point(:left, pos_left, WING;
    extra_mass=k_side, wing=:plate_wing,
    transform=:kite_tilt))

pos_map = Dict(:kcu => pos_kcu, :nose => pos_nose,
               :top => pos_top, :right => pos_right,
               :left => pos_left)
bridle_l0(a, b) = norm(pos_map[b] - pos_map[a]) * PRE_STRESS

segments = [
    Segment(:kcu_nose,   set, :kcu,   :nose;
        l0=bridle_l0(:kcu, :nose),
        diameter_mm=set.d_line),
    Segment(:right_nose, set, :right, :nose;
        l0=bridle_l0(:right, :nose),
        diameter_mm=set.d_line),
    Segment(:right_left, set, :right, :left;
        l0=bridle_l0(:right, :left),
        diameter_mm=set.d_line),
    Segment(:top_right,  set, :top,   :right;
        l0=bridle_l0(:top, :right),
        diameter_mm=set.d_line),
    Segment(:left_kcu,   set, :left,  :kcu;
        l0=bridle_l0(:left, :kcu),
        diameter_mm=set.d_line),
    Segment(:right_kcu,  set, :right, :kcu;
        l0=bridle_l0(:right, :kcu),
        diameter_mm=set.d_line),
    Segment(:top_left,   set, :top,   :left;
        l0=bridle_l0(:top, :left),
        diameter_mm=set.d_line),
    Segment(:left_nose,  set, :left,  :nose;
        l0=bridle_l0(:left, :nose),
        diameter_mm=set.d_line),
    Segment(:nose_top,   set, :nose,  :top;
        l0=bridle_l0(:nose, :top),
        diameter_mm=set.d_line),
]

tethers = [Tether(:main_tether, set.l_tethers[1];
    start_point=:ground, end_point=:kcu,
    n_segments=set.segments)]

winches = [Winch(:winch, set, [:main_tether];
                 winch_point=:ground)]

# Plate surfaces and wing
rel_side_area = set.rel_side_area / 100.0
K = 1.0 - rel_side_area

surfaces = [
    PlateSurface(:main, [1,0,0], [0,1,0],
        set.area, :top;
        twist=deg2rad(set.alpha_zero)),
    PlateSurface(:right_tip, [1,0,0], [0,0,-1],
        set.area * rel_side_area, :right;
        twist=deg2rad(set.alpha_ztip)),
    PlateSurface(:left_tip, [1,0,0], [0,0,1],
        set.area * rel_side_area, :left;
        twist=deg2rad(set.alpha_ztip)),
]

cl_interp, cd_interp = create_plate_interpolations(
    set.alpha_cl, set.cl_list, set.cd_list;
    alpha_cd=set.alpha_cd)

plate_wing = PlateWing(
    :plate_wing, surfaces, cl_interp, cd_interp;
    dynamics_type=PARTICLE_DYNAMICS,
    z_ref_points=([:right, :left], :top),
    y_ref_points=(:left, :right),
    origin=:kcu,
    drag_corr=0.93 * K,
    cmq=set.cmq, smc=set.smc,
    cord_length=set.cord_length)

alpha_depower = calc_alpha_depower(KCU(set), 0.25)
plate_wing.surfaces[1].twist =
    deg2rad(set.alpha_zero) - alpha_depower

KITE_ANGLE = 3.83
transforms = [
    Transform(:main_tf,
        deg2rad(set.elevation), 0.0, 0.0;
        base_pos=zeros(3), base_point=:ground,
        wing=:plate_wing),
    Transform(:kite_tilt,
        deg2rad(set.elevation - KITE_ANGLE),
        0.0, 0.0;
        base_transform=:main_tf,
        rot_point=:top),
]

sys = SystemStructure("kps4", set;
    points, segments, tethers, winches,
    wings=[plate_wing], transforms)
sys.winches[1].brake = true

sam = SymbolicAWEModel(set, sys)
init!(sam; remake=false, prn=true)

w = sam.sys_struct.wings[1]
aoas = [round(s.aoa, digits=2) for s in w.surfaces]
println("SymAWE initial aoa=$(aoas)")

sam_logger = Logger(sam, N_STEPS + 1)
sys_state = SysState(sam)
sys_state.time = 0.0
log!(sam_logger, sys_state)

# Warmup step
next_step!(sam; dt, set_values=[0.0])
update_sys_state!(sys_state, sam)
sys_state.time = dt
log!(sam_logger, sys_state)

sam_elapsed = 0.0
for step in 2:N_STEPS
    global sam_elapsed += @elapsed next_step!(
        sam; dt, set_values=[0.0])
    update_sys_state!(sys_state, sam)
    sys_state.time = step * dt
    log!(sam_logger, sys_state)
end
sam_rt = ((N_STEPS - 1) * dt) / sam_elapsed
println("SymAWE realtime factor: $(round(sam_rt, digits=2))")

w = sam.sys_struct.wings[1]
va_dir = normalize(w.va_b)
drag_val = w.aero_force_b ⋅ va_dir
lift_val = norm(w.aero_force_b - drag_val * va_dir)
println("SymAWE final lift, drag [N]: " *
    "$(round(lift_val, digits=2)), " *
    "$(round(drag_val, digits=2))")

save_log(sam_logger, "kps4_sam")
sam_syslog = load_log("kps4_sam")

# ==================== COMPARISON PLOTS ==================== #
println()
println("=" ^ 60)
println("Plotting comparison...")
println("=" ^ 60)

fig = plot(
    [sam.sys_struct, sam.sys_struct],
    [sam_syslog, km_syslog];
    suffixes=["SymAWE", "KiteModels"],
    size=(1200, 800),
    plot_default=false,
    plot_elevation=true,
    plot_azimuth=true)
display(fig)
