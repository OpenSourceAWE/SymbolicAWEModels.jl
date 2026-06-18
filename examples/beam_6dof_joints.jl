# Copyright (c) 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Horizontal beam built from a chain of rigid-body segments connected by 6-DOF
# elastic joints — the multi-rigid-body model (e.g. an inflatable leading-edge
# tube). No wind. The root segment is fixed; the beam sags under gravity and
# settles into a static cantilever deflection set by the joints' bending
# stiffness EI. The motion is logged to a SysLog and shown with `replay`.
# Set `g_earth: 0.0` in data/beam/settings.yaml for a free, gravity-free beam.

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

# ----- settings -----
# Environment settings (no wind; gravity) live in the committed data dir
# `data/beam/settings.yaml`; edit `g_earth` there for a gravity-free beam. Using
# the package data dir keeps the compiled model cached and reused across runs.
set_data_path(joinpath(dirname(@__DIR__), "data", "beam"))
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
