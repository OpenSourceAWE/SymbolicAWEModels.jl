# Copyright (c) 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Horizontal beam built from a chain of rigid-body segments connected by 6-DOF
# elastic joints — the multi-rigid-body model (e.g. an inflatable leading-edge
# tube). No wind. The root segment is fixed; the beam sags under gravity and
# settles into a static cantilever deflection set by the joints' bending
# stiffness EI. The motion is logged to a SysLog and shown with `replay`.
# Set `g_earth = 0.0` below for a free, gravity-free beam.

using Pkg
Pkg.activate(@__DIR__)
using SymbolicAWEModels
using KiteUtils
using GLMakie

# ----- beam parameters -----
n_segments = 10
seg_length = 0.5                 # [m] per segment
seg_mass = 0.5                   # [kg] per segment
seg_radius = 0.02                # [m] equivalent rod radius (axial inertia)
g_earth = 9.81                   # set to 0.0 for a gravity-free beam

# Slender-rod principal inertia: small about the long (x) axis, m·L²/12 about
# the transverse axes.
inertia = [0.5 * seg_mass * seg_radius^2,
           seg_mass * seg_length^2 / 12,
           seg_mass * seg_length^2 / 12]

# Joint stiffness: stiff in stretch/shear, finite bending/torsion.
stiffness_axial = 1.0e5
stiffness_shear = 1.0e5
stiffness_torsion = 5.0e3
stiffness_bending = 5.0e3        # ↓ for a floppier beam, ↑ for a stiffer one
damping_trans = 50.0
damping_rot = 20.0

# ----- settings (no wind; gravity as configured) -----
settings_yaml = """
system:
    log_file: "data/beam"
    g_earth: $g_earth
solver:
    solver: "FBDF"
    abs_tol: 1.0e-7
    rel_tol: 1.0e-7
kite:
    model: ""
    foil_file: "ram_air_kite/ram_air_kite_foil.dat"
    physical_model: "beam_6dof"
    mass: 0.0
    quasi_static: false
tether:
    cd_tether: 0.958
    unit_damping: 0.0
    unit_stiffness: 0.0
    rho_tether: 724.0
    e_tether: 5.5e10
winch:
    winch_model: "TorqueControlledMachine"
    drum_radius: 0.110
    gear_ratio: 1.0
    inertia_total: 0.024
    f_coulomb: 122.0
    c_vf: 30.6
environment:
    rho_0: 1.225
    v_wind: 0.0
    upwind_dir: -90.0
    upwind_elevation: 0.0
    wind_vec: [0.0, 0.0, 0.0]
    profile_law: 0
"""

pkg_root = dirname(@__DIR__)
tmpdir = mktempdir()
data_path = joinpath(tmpdir, "2plate_kite")
cp(joinpath(pkg_root, "data", "2plate_kite"), data_path; force=true)
write(joinpath(data_path, "settings.yaml"), settings_yaml)
write(joinpath(data_path, "system.yaml"),
    "system:\n  sim_settings: settings.yaml\n")
set_data_path(data_path)
set = Settings("system.yaml")

# ----- build the beam -----
bodies = RigidBody[]
for i in 1:n_segments
    origin_x = (i - 0.5) * seg_length        # segment centered at its origin
    push!(bodies, RigidBody(Symbol("seg_$i");
        mass=seg_mass, inertia_principal=inertia,
        pos=[origin_x, 0.0, 0.0],
        fixed=(i == 1)))                       # root segment clamped
end

joints = ElasticJoint[]
for i in 1:(n_segments - 1)
    push!(joints, ElasticJoint(Symbol("joint_$i"),
        Symbol("seg_$i"), Symbol("seg_$(i+1)");
        anchor_a=[seg_length / 2, 0.0, 0.0],   # right end of segment i
        anchor_b=[-seg_length / 2, 0.0, 0.0],  # left end of segment i+1
        stiffness_axial, stiffness_shear,
        stiffness_torsion, stiffness_bending,
        damping_trans, damping_rot))
end

sys = SystemStructure("beam_6dof", set;
    rigid_bodies=bodies, elastic_joints=joints)
sam = SymbolicAWEModel(set, sys)
init!(sam)

# ----- simulate + log -----
dt = 0.02
t_end = 5.0
n_steps = round(Int, t_end / dt)
logger = Logger(sam, n_steps + 1)
sys_state = SysState(sam)
tip = sam.sys_struct.rigid_bodies[Symbol("seg_$n_segments")]

for step in 0:n_steps
    step > 0 && next_step!(sam; dt, vsm_interval=0)
    update_sys_state!(sys_state, sam)
    sys_state.time = step * dt
    log!(logger, sys_state)
end
println("Tip deflection after $t_end s: $(round(tip.pos_w[3]; digits=4)) m")

# ----- save, reload, replay -----
save_log(logger, "beam_6dof")
syslog = load_log("beam_6dof")
# vector_scale shrinks the per-body RGB frame arrows to suit 0.5 m segments.
scene = replay(syslog, sam.sys_struct; vector_scale=0.3)
display(scene)

nothing
