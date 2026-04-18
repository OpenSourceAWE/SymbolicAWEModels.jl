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
using KiteUtils
using KiteUtils: init!, next_step!, update_sys_state!
using LinearAlgebra

pkg_root = dirname(@__DIR__)
set_data_path(joinpath(pkg_root, "data", "kps4"))
set = deepcopy(load_settings("system.yaml"))

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
    extra_mass=set.kcu_mass))
push!(points, Point(:nose, pos_nose, DYNAMIC;
    extra_mass=k_nose))
push!(points, Point(:top, pos_top, WING;
    extra_mass=k_top, wing=:plate_wing))
push!(points, Point(:right, pos_right, WING;
    extra_mass=k_side, wing=:plate_wing))
push!(points, Point(:left, pos_left, WING;
    extra_mass=k_side, wing=:plate_wing))

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

surfaces = [
    # Main surface: chord=x, span=y (xz-plane projection)
    PlateSurface(:main, [1,0,0], [0,1,0],
        set.area * (1.0 - rel_side_area), :top;
        twist_b=-1.0, alpha_offset=set.alpha_zero),
    # Right tip: chord=x, span=z (xy-plane projection)
    PlateSurface(:right_tip, [1,0,0], [0,0,1],
        set.area * rel_side_area / 2, :right;
        twist_a=1.0, alpha_offset=set.alpha_ztip),
    # Left tip: chord=x, span=-z (mirrored)
    PlateSurface(:left_tip, [1,0,0], [0,0,-1],
        set.area * rel_side_area / 2, :left;
        twist_a=-1.0, alpha_offset=set.alpha_ztip),
]

# ==================== CL/CD INTERPOLATIONS ================ #
cl_interp, cd_interp = create_plate_interpolations(
    set.alpha_cl, set.cl_list, set.cd_list;
    alpha_cd=set.alpha_cd)

# ==================== PLATE WING ========================== #
wing = PlateWing(:plate_wing, surfaces, cl_interp, cd_interp;
    wing_type=REFINE,
    z_ref_points=(:kcu, :top),
    y_ref_points=(:right, :left),
    origin=:kcu,
    drag_corr=0.93,
    cmq=set.cmq, smc=set.smc,
    cord_length=set.cord_length)

# ==================== SYSTEM STRUCTURE ==================== #
sys = SystemStructure("kps4", set;
    points, segments, tethers, winches, wings=[wing])

# ==================== MODEL + INIT ======================== #
sam = SymbolicAWEModel(set, sys)
init!(sam; prn=true)

# ==================== SIMULATE + LOG ====================== #
logger = Logger(sam, N_STEPS + 1)
sys_state = SysState(sam)
sys_state.time = 0.0
log!(logger, sys_state)

println("Simulating $(SIM_TIME)s parking flight...")
total_elapsed = @elapsed for step in 1:N_STEPS
    next_step!(sam; dt, set_values=[0.0])

    update_sys_state!(sys_state, sam)
    sys_state.time = step * dt
    log!(logger, sys_state)
end
avg_rt = SIM_TIME / total_elapsed
h = round(sam.sys_struct.points[:top].pos_w[3]; digits=1)
println("Done. h=$(h)m, avg $(round(avg_rt; digits=1))x RT")

# ==================== SAVE + REPLAY ======================= #
save_log(logger, "kps4_run")
syslog = load_log("kps4_run")
scene = replay(syslog, sam.sys_struct;
               autoplay=false, loop=true)
display(scene)
