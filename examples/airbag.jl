# Copyright (c) 2025 Bart van de Lint, Jelle Poland
# SPDX-License-Identifier: LGPL-3.0-only

"""
Airbag simulation: a pressurized square membrane (N×N particle grid)
inflated by a uniform internal gauge pressure. The boundary nodes are
fixed (STATIC); interior nodes are dynamic. At each time step the panel
normals are recomputed from the current geometry and the resulting
pressure forces (P × triangle_area / 3 per corner) are fed in as
`disturb` forces on the dynamic nodes.

The membrane starts flat; the pressure inflates it into a pillow shape.

Reference: Thedens, P. Dissertation (2022), p.45
           https://github.com/awegroup/Particle_System_Simulator
"""

using Pkg
if Base.active_project() != joinpath(@__DIR__, "Project.toml")
    Pkg.activate(@__DIR__)
end

using GLMakie
using KiteUtils: init!, next_step!, update_sys_state!
using SymbolicAWEModels
import SymbolicAWEModels: Point   # resolve ambiguity with GLMakie
using KiteUtils
using LinearAlgebra

# ── Parameters ─────────────────────────────────────────────────────────
const N         = 4        # grid divisions per side (N+1 nodes per edge)
const L         = 1.0      # side length [m]
const Z0        = 0.0      # membrane height above ground [m]
const P_GAUGE   = 100.0    # internal gauge pressure [Pa]

# Membrane material
const SEG_STIFF = 2_000.0  # unit stiffness EA [N]
const SEG_DAMP  = 50.0     # unit damping [N·s]
const SEG_DIA   = 0.002    # segment diameter [m]
const NODE_MASS = 0.01     # lumped particle mass [kg]
const WF_DAMP   = 2.0      # world-frame velocity damping [N·s/m]
const Z_PERTURB = 0.002    # tiny initial z-offset for interior nodes [m]
                            # (breaks flat-plane symmetry to aid solver)

# ── Helper functions ────────────────────────────────────────────────────
pname(i, j) = Symbol("p_$(i)_$(j)")

# 1-based linear index: outer loop j=0:N, inner loop i=0:N
pidx(i, j) = j * (N + 1) + i + 1

# ── Settings ────────────────────────────────────────────────────────────
set_data_path(joinpath(dirname(@__DIR__), "data"))
set = Settings("base/system.yaml")
set.v_wind = 0.0

dx = L / N

# ── Build point grid ────────────────────────────────────────────────────
points = Point[]
for j in 0:N
    for i in 0:N
        is_boundary = (i == 0 || i == N || j == 0 || j == N)
        type = is_boundary ? STATIC : DYNAMIC
        # Small upward perturbation on interior nodes to break symmetry
        z = is_boundary ? Z0 : Z0 + Z_PERTURB
        kw = is_boundary ? (;) :
            (; extra_mass=NODE_MASS, world_frame_damping=WF_DAMP)
        push!(points, Point(pname(i, j),
            [i * dx, j * dx, z], type; kw...))
    end
end

# ── Build segments ──────────────────────────────────────────────────────
segments = Segment[]
seg_id   = 1

# Structural edges: horizontal (x-direction) and vertical (y-direction)
for j in 0:N, i in 0:N
    global seg_id
    if i < N
        push!(segments, Segment(seg_id, pname(i, j), pname(i + 1, j),
              SEG_STIFF, SEG_DAMP, SEG_DIA; compression_frac=0.0))
        seg_id += 1
    end
    if j < N
        push!(segments, Segment(seg_id, pname(i, j), pname(i, j + 1),
              SEG_STIFF, SEG_DAMP, SEG_DIA; compression_frac=0.0))
        seg_id += 1
    end
end

# Diagonal shear segments (both diagonals per quad cell)
for j in 0:N-1, i in 0:N-1
    global seg_id
    push!(segments, Segment(seg_id, pname(i, j), pname(i + 1, j + 1),
          SEG_STIFF, SEG_DAMP, SEG_DIA; compression_frac=0.0))
    seg_id += 1
    push!(segments, Segment(seg_id, pname(i + 1, j), pname(i, j + 1),
          SEG_STIFF, SEG_DAMP, SEG_DIA; compression_frac=0.0))
    seg_id += 1
end
# Rest lengths (l0) are automatically set to initial CAD distances by SystemStructure.

# ── Build model ─────────────────────────────────────────────────────────
sys = SystemStructure("airbag", set; points, segments)
sam = SymbolicAWEModel(set, sys)
init!(sam; remake=false)

# ── Pressure force computation ───────────────────────────────────────────
"""
    apply_pressure!(sam, pressure)

Compute outward-normal pressure forces from the current panel geometry
and inject them as `disturb` forces on the dynamic interior nodes.

Each quad cell (i,j) is split into two triangles. The force on each
triangle is `P × area_triangle`. This force is divided equally among
the three corners, and only DYNAMIC corners accumulate the contribution.
The integrator holds a reference to the mutable `sys_struct`, so
mutating `disturb` fields takes effect at the next integration step.
"""
function apply_pressure!(sam, pressure)
    pts = sam.sys_struct.points

    # Reset all disturb forces to zero
    for pt in pts
        pt.disturb .= 0.0
    end

    for j in 0:N-1, i in 0:N-1
        ia = pidx(i,     j    )
        ib = pidx(i + 1, j    )
        ic = pidx(i + 1, j + 1)
        id = pidx(i,     j + 1)

        pa = pts[ia].pos_w
        pb = pts[ib].pos_w
        pc = pts[ic].pos_w
        pd = pts[id].pos_w

        # Triangle 1: a – b – c  (area-weighted outward normal via cross product)
        n1   = (pb .- pa) × (pc .- pa)
        mag1 = norm(n1)
        if mag1 > 1e-12
            # P × triangle_area / 3 per corner; triangle_area = |n|/2
            fvec1 = (pressure / 6) .* n1
            for idx in (ia, ib, ic)
                pts[idx].type == DYNAMIC && (pts[idx].disturb .+= fvec1)
            end
        end

        # Triangle 2: a – c – d
        n2   = (pc .- pa) × (pd .- pa)
        mag2 = norm(n2)
        if mag2 > 1e-12
            fvec2 = (pressure / 6) .* n2
            for idx in (ia, ic, id)
                pts[idx].type == DYNAMIC && (pts[idx].disturb .+= fvec2)
            end
        end
    end
    return nothing
end

# ── Simulation loop ──────────────────────────────────────────────────────
n_steps  = 50
logger   = Logger(sam, n_steps)
sys_state = SysState(sam)

for i in 1:n_steps
    apply_pressure!(sam, P_GAUGE)
    next_step!(sam; dt=0.001)
    update_sys_state!(sys_state, sam)
    sys_state.time = i / set.sample_freq
    log!(logger, sys_state)
end

save_log(logger, "airbag")
syslog = load_log("airbag")
scene  = replay(syslog, sam.sys_struct)
display(scene)
