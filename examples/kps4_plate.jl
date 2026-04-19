# KPS4 flat-plate kite model using SymbolicAWEModels
#
# Replicates KiteModels.jl's simulate_simple.jl example:
# - 4 kite particles (KCU, nose, right tip, left tip)
# - 6-segment main tether
# - 9 bridle springs
# - Flat-plate CL/CD aerodynamics
#
# Uses the same system.yaml settings as KiteModels.jl.

using Pkg
if Base.active_project() != joinpath(@__DIR__, "Project.toml")
    Pkg.activate(@__DIR__)
end

using GLMakie
using SymbolicAWEModels
using SymbolicAWEModels: Point
using KitePodModels
using KiteUtils: init!, next_step!, update_sys_state!
using LinearAlgebra

pkg_root = dirname(@__DIR__)
set_data_path(joinpath(pkg_root, "data", "kps4"))
set = deepcopy(load_settings("system.yaml"))
UPWIND_DIR = -pi/2 + deg2rad(10)
set.upwind_dir = rad2deg(UPWIND_DIR)

# Simulation parameters
SIM_TIME = 30.0
dt = 0.05
N_STEPS = round(Int, SIM_TIME / dt)

# KPS4 constants
PRE_STRESS = 0.9975

# ==================== GEOMETRY ============================ #
particles = KiteUtils.get_particles(
    set.height_k, set.h_bridle, set.width, set.m_k)
# particles = [origin, kcu, nose(A), top(B), right(C), left(D)]
pos_kcu = particles[2]
pos_nose = particles[3]
pos_top = particles[4]
pos_right = particles[5]
pos_left = particles[6]

# ==================== POINTS ============================== #
# Mass distribution from KPS4.init_masses!
kite_mass = set.mass
k_nose = set.rel_nose_mass * kite_mass
k_top = set.rel_top_mass *
        (1.0 - set.rel_nose_mass) * kite_mass
k_side = 0.5 * (1.0 - set.rel_top_mass) *
         (1.0 - set.rel_nose_mass) * kite_mass
set.mass = 0.0  # masses set explicitly on points

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

# ==================== SEGMENTS ============================ #
pos_map = Dict(:kcu => pos_kcu, :nose => pos_nose,
               :top => pos_top, :right => pos_right,
               :left => pos_left)
bridle_l0(a, b) = norm(pos_map[b] - pos_map[a]) * PRE_STRESS

segments = [
    Segment(:kcu_nose,   set, :kcu,   :nose;
        l0=bridle_l0(:kcu, :nose), diameter_mm=set.d_line),
    Segment(:right_nose, set, :right, :nose;
        l0=bridle_l0(:right, :nose), diameter_mm=set.d_line),
    Segment(:right_left, set, :right, :left;
        l0=bridle_l0(:right, :left), diameter_mm=set.d_line),
    Segment(:top_right,  set, :top,   :right;
        l0=bridle_l0(:top, :right), diameter_mm=set.d_line),
    Segment(:left_kcu,   set, :left,  :kcu;
        l0=bridle_l0(:left, :kcu), diameter_mm=set.d_line),
    Segment(:right_kcu,  set, :right, :kcu;
        l0=bridle_l0(:right, :kcu), diameter_mm=set.d_line),
    Segment(:top_left,   set, :top,   :left;
        l0=bridle_l0(:top, :left), diameter_mm=set.d_line),
    Segment(:left_nose,  set, :left,  :nose;
        l0=bridle_l0(:left, :nose), diameter_mm=set.d_line),
    Segment(:nose_top,   set, :nose,  :top;
        l0=bridle_l0(:nose, :top), diameter_mm=set.d_line),
]

# ==================== TETHER ============================== #
tethers = [Tether(:main_tether, set.l_tethers[1];
    start_point=:ground, end_point=:kcu,
    n_segments=set.segments)]

# ==================== WINCH =============================== #
winches = [Winch(:winch, set, [:main_tether];
                 winch_point=:ground)]

# ==================== PLATE SURFACES ====================== #
# Body frame (same as VSMWing): x=chord, y=span right, z=down
# Surface axes are in body frame, transformed to world by R_b_w.
rel_side_area = set.rel_side_area / 100.0
K = 1.0 - rel_side_area  # drag area correction (0.694)

surfaces = [
    # Main surface: full area (KPS4 convention)
    PlateSurface(:main, [1,0,0], [0,1,0],
        set.area, :top;
        twist=deg2rad(set.alpha_zero)),
    # Right tip: full side area fraction
    PlateSurface(:right_tip, [1,0,0], [0,0,-1],
        set.area * rel_side_area, :right;
        twist=deg2rad(set.alpha_ztip)),
    # Left tip: mirrored, full side area fraction
    PlateSurface(:left_tip, [1,0,0], [0,0,1],
        set.area * rel_side_area, :left;
        twist=deg2rad(set.alpha_ztip)),
]

# ==================== CL/CD INTERPOLATIONS ================ #
cl_interp, cd_interp = create_plate_interpolations(
    set.alpha_cl, set.cl_list, set.cd_list;
    alpha_cd=set.alpha_cd)

# ==================== PLATE WING ========================== #
wing = PlateWing(:plate_wing, surfaces, cl_interp, cd_interp;
    wing_type=REFINE,
    z_ref_points=([:right, :left], :top),
    y_ref_points=(:left, :right),
    origin=:kcu,
    drag_corr=0.93 * K,
    cmq=set.cmq, smc=set.smc,
    cord_length=set.cord_length)

# Set initial twist: depower on main surface (KPS4 convention)
alpha_depower = calc_alpha_depower(KCU(set), 0.25)
wing.surfaces[1].twist =
    deg2rad(set.alpha_zero) - alpha_depower

# ==================== TRANSFORMS =========================== #
KITE_ANGLE = 3.83
transforms = [
    # Position KCU/nose at elevation + azimuth
    Transform(:main_tf,
        deg2rad(set.elevation), deg2rad(10.0), 0.0;
        base_pos=zeros(3), base_point=:ground,
        wing=:plate_wing),
    # Tilt kite body by KITE_ANGLE on top of main
    Transform(:kite_tilt,
        deg2rad(set.elevation - KITE_ANGLE),
        deg2rad(10.0), 0.0;
        base_transform=:main_tf,
        rot_point=:top),
]

# ==================== SYSTEM STRUCTURE ==================== #
sys = SystemStructure("kps4", set;
    points, segments, tethers, winches,
    wings=[wing], transforms)
sys.winches[1].brake = true

# ==================== MODEL + INIT ======================== #
sam = SymbolicAWEModel(set, sys)
init!(sam; remake=false, prn=true)
# find_steady_state!(sam)
w = sam.sys_struct.wings[1]
aoas = [round(s.aoa, digits=2) for s in w.surfaces]
println("Initial: aoa=$(aoas)°, elevation=$(round(
    rad2deg(w.elevation), digits=2))°, " *
    "va_b=$(round.(w.va_b, digits=2)), " *
    "v_wind=$(round.(w.v_wind, digits=2))")
top_pt = sam.sys_struct.points[:top]
wf = SymbolicAWEModels.calc_wind_factor(
    sam.sys_struct.am, 0.0, 0.0, top_pt.pos_w[3],
    sam.sys_struct)
println("height=$(round(top_pt.pos_w[3],
    digits=2))m, " *
    "wind_factor=$(round(wf, digits=4)), " *
    "v_wind_gnd=$(round.(sam.set.wind_vec, digits=2))")
Rbw = w.R_b_to_w
println("R_b_to_w:\n  x=$(round.(Rbw[:,1], digits=3))" *
    "\n  y=$(round.(Rbw[:,2], digits=3))" *
    "\n  z=$(round.(Rbw[:,3], digits=3))")
va_dir = normalize(w.va_b)
init_drag = w.aero_force_b ⋅ va_dir
init_lift = norm(w.aero_force_b - init_drag * va_dir)
println("Initial lift, drag  [N]: " *
    "$(round(init_lift, digits=2)), " *
    "$(round(init_drag, digits=2))")

# ==================== SIMULATE + LOG ====================== #
logger = Logger(sam, N_STEPS + 1)
sys_state = SysState(sam)
sys_state.time = 0.0
log!(logger, sys_state)

println("Simulating $(SIM_TIME)s parking flight...")
# Warmup
next_step!(sam; dt, set_values=[0.0])
update_sys_state!(sys_state, sam)
sys_state.time = dt
log!(logger, sys_state)

sim_elapsed = 0.0
integ_elapsed = 0.0
for step in 2:N_STEPS
    global sim_elapsed += @elapsed next_step!(
        sam; dt, set_values=[0.0])
    global integ_elapsed += sam.t_step

    update_sys_state!(sys_state, sam)
    sys_state.time = step * dt
    log!(logger, sys_state)
end
sim_time = (N_STEPS - 1) * dt
avg_rt = sim_time / sim_elapsed
integ_rt = sim_time / integ_elapsed
println("Avg realtime factor: $(round(avg_rt, digits=2))")
println("Integrator realtime factor: " *
    "$(round(integ_rt, digits=2))")

# Lift and drag in body frame
w = sam.sys_struct.wings[1]
va_dir = normalize(w.va_b)
drag_val = w.aero_force_b ⋅ va_dir
lift_val = norm(w.aero_force_b - drag_val * va_dir)
println("Ground wind speed: $(round(set.v_wind, digits=2)) m/s")
println("lift, drag  [N]: $(round(lift_val, digits=2)), " *
        "$(round(drag_val, digits=2))")

# ==================== SAVE + REPLAY ======================= #
save_log(logger, "kps4_run")
syslog = load_log("kps4_run")
scene = replay(syslog, sam.sys_struct;
               autoplay=false, loop=true)
display(scene)
