# NetworkDynamics backend — status report

**Branch:** `nd-bodies-autonomous` (off `networkdynamics-backend`)
**North star:** a *complete* alternative backend with *all* monolith features;
acceptance = every monolith test passes unchanged except swapping `backend=`, results
matching to solver tolerance (~1e-10 instantaneous RHS parity).

This branch adds **rigid bodies, BODY_STATIC ride points, Timoshenko/Elastic beam joints,
rigid-wing VSM aero, a winch + tether attaching to the wing body, and the DYNAMIC
twist-surface DOF** to the NetworkBackend, all validated by exact backend-swap parity
against the monolith. The final structural gap — the combined multi-segment edge for kite
steering pulleys (ND's SimpleGraph forbids the parallel edges they create) — is now
**closed**: the full SK100 structural model assembles and evaluates its RHS on the ND
backend (state 67, warm compile ~7.9 s, first-RHS f_ip JIT ~5.3 s); full-run parity is the
next step. See §3.

**Two full monolith test suites now pass on the NetworkBackend by backend swap only:**
`test_rigid_body.jl` **19/19** and `test_timoshenko_joint.jl` **20/20** (the latter checks
analytical beam physics — Timoshenko-vs-Euler deflection, axial, torsion, and two
nonlinear callable-rigidity softening laws — so this is *correctness*, not just parity).

---

## 1. Delivered and validated on this branch

Each item is a shared-component + ND-edge/vertex pair validated by building the *same*
`SystemStructure` on both backends and comparing after N steps.

| Feature | Parity vs monolith | Model |
|---|---|---|
| **Singleton rigid body** (6-DOF vertex) | exact | free-fall (a=g), spin-up (α=τ/I), COM-offset |
| **BODY_STATIC ride point + wrench edge** | Δ=0.0 | body + offset ride point + spring (ω_b=9.18) |
| **Timoshenko joint edge** (body↔body beam) | Δ~1e-12 | two bodies + joint, A pitched by ext_moment |
| **Elastic joint edge** (body↔body lumped) | Δ~2.6e-12 | two bodies + joint, A pitched by ext_moment |
| **Body↔body segment** (dual wrench) | Δ~4e-12 | bridle line between two plates via ride points |
| Point-only (SK100-style) regression | Δ=0.0 | mass-on-spring, unchanged narrow path |

### The design (pure ODE, no DAE — the DAE approach was vetoed)
ND's coreloop forbids general feedthrough vertices, so a ride point that forwards its
aggregated segment force back to a body can't be a plain vertex. Resolution
(**cluster-vertex**, agreed design):

- A **rigid body is one vertex** carrying the 13-state principal pose
  (`com_w`,`com_vel`,`Q_p_to_w`,`ω_p`); its pose is an *output* (function of its own
  states → non-feedthrough).
- A **BODY_STATIC ride point is absorbed** into its body vertex — it is not a graph
  vertex. Its position is reconstructed inside the incident **wrench edge** from the
  body pose + the point's `anchor_b`, and the edge delivers force + `arm×force` moment
  back to the body. Pure ODE, no recompute.
- A **TimoshenkoJoint between two bodies is a wide edge**: it reads both bodies' poses,
  evaluates the corotational element wrench, and delivers the equal-and-opposite
  restoring wrench to each. Joints are stateless — the DOF live in the node-bodies.

### Wide interface (§8.5) — per-network, zero cost to SK100
ND fixes one vertex-output width and one edge-output width per graph. Body-containing
models use a **wide superset**: vertex output `[pos, vel, pulley_len, pose_R(9),
pose_com, pose_com_vel, pose_omega]` (points zero-fill the pose slots); vertex input
`[force, mass, tension, moment]`; edge I/O widened to match with a `moment` output.
A `wide` flag threads through the point/segment builders; point-only models keep the
**narrow** interface, so SK100 pays nothing and does not regress.

### No duplicated physics
The 6-DOF body math (`rigid_body_pose_expressions`) and the Timoshenko element wrench
(`timoshenko_element_wrench`) are backend-agnostic helpers in `components.jl`; the
monolith `rigid_body_eqs!` and `timoshenko_joint_eqs!` delegate to them. Monolith
`test_rigid_body` (19/19) and `test_timoshenko_joint` (20/20) pass unchanged after the
refactors.

---

## 2. Key lessons (recorded so they don't recur)

- **Edgeless ND networks** need a fork patch: with no edges, `edepth` must equal the
  vertices' input width (not 0), else the zero-filled aggregation view is 0-wide and
  the vertex reads its inputs out of bounds → NaN. (Fork `nonnumeric-params`,
  `construction.jl`.)
- **No runtime `ifelse` on a build-time constant.** The wrench edge first selected the
  body endpoint with `ifelse(body_at_src, src_pose, dst_pose)`; that `ifelse`
  propagated through the matrix/cross/spring expressions and **hung `mtkcompile`
  indefinitely**. `body_at_src` is known at assembly time → bake it into the kernel
  (one specialized kernel per orientation). Compiles in seconds.
- The fork's `MTKExt_utils.jl:218` "dropped metadata … could NOT be fixed" `@warn` is
  **benign** for these kernels (outputs and derivatives are numerically exact).
- **Tear nested pose→geometry→force expressions** into edge-internal variables, and use
  **scalar (not vector) `make_param`** for explicit edge params — both otherwise stall
  the compile / trip the fork's I/O scalarization (see the dual-wrench edge).
- **Edge nonnumeric (callable) params work** (not just vertex): joint rigidity laws
  `EIy(κ)` route through the fork's nonnumeric SII path exactly like the aero polars.

## Hermite ride-point design (two-edge, no clustering) — DONE

**Implemented and parity-validated** (commit `37301974`; Δω=3.9e-7). A point riding a
`TimoshenkoJoint`'s deformed centerline depends on **both** end bodies (feedthrough on two
vertices). The plan's "merge into one cluster vertex" would force doubling the uniform
vertex interface (2 poses + 2 wrench channels on *every* vertex) — large churn. The
cleaner equivalent shipped here: a bridle knot `K` attached to a Hermite point on
joint(a,b) becomes **two wrench edges**, `K↔a` and `K↔b`. Each reads the *other* body's
pose via **`extin`** (edge extin is proven — the winched-tether segment uses it), so both
compute the *identical* Hermite ride position + spring force `F`; edge `K↔a` applies the
`(1−sfrac)` share to `a` and `−(1−sfrac)F` to `K`, edge `K↔b` the `sfrac` share to `b`
and `−sfrac F` to `K`, so `K` gets the full reaction `−F` and each body its beam-fraction
split. No vertex widening, no cluster. Shipped as: shared `beam_hermite_ride_expressions`
(components.jl, with torn frame/theta/pos + a `ride_velocity` accessor so velocity uses
the torn position), a `network_hermite_wrench_segment` kernel (torn intermediates + extin
to the other body), and assembly that emits the two edges + wires the extin. **Required a
fork fix** (`external_inputs.jl:45`, `463affe`): the Hermite edge is the first to read a
component *output* (body pose) via `extin`, which indexed `sym(cm)` (states) → "invalid
index: nothing"; fixed to `outsym_flat(cm)`. This was the last structural gap before
rigid-VSM aero.

---

## 3. Remaining to the north star

| Layer | State |
|---|---|
| Singleton rigid body | ✅ exact |
| BODY_STATIC ride points + wrench edges | ✅ Δ=0.0 |
| TimoshenkoJoint body↔body edge | ✅ Δ~1e-12 |
| ElasticJoint edge | ✅ Δ~2.6e-12 |
| Body↔body segments (bridle between plates) | ✅ Δ~4e-12 |
| STATIC (clamped) body vertex | ✅ Δ~9e-12 |
| Full `test_yaml_bodies` model by backend swap | ✅ Δ~9e-12 (STATIC+Timo+BODY_STATIC) |
| Multi-body Timoshenko beam (3-body chain + 2 joint edges) | ✅ Δ~6e-10 (no new code — bodies=vertices, joints=edges) |
| Callable (nonlinear) joint rigidities `EIy(κ)` | ✅ Δ~5e-11 (edge nonnumeric params via the fork) |
| Full `test_timoshenko_joint.jl` on NetworkBackend | ✅ **20/20** by backend swap (analytical physics: Timoshenko-vs-Euler, axial, torsion, 2× nonlinear) |
| **Hermite ride point on a beam joint** (two-edge + extin) | ✅ **PARITY PASS** — Δω=3.9e-7, Δpos=4.9e-9, ΔQ=7.7e-9 (with segment damping 0). High segment damping numerically amplifies a ~1e-7 backend expression-order difference to ~6.5e-6 (not a physics bug). |
| **→ the whole structural layer (bodies/joints/beams/ride points) is DONE** | ✅ |
| **RIGID wing VSM aero on a body vertex** (AeroDirect + AeroLinearized) | ✅ **PARITY PASS** — Δvel~1e-10 after 10 steps, both modes exact. `WingBodyVertex` = `BodyVertex` + the rigid `aero_component` nested, connectors driven from the body's own pose (no extin), force/moment transported to COM and folded into the wrench. Getter scatters `va_b`/`R_b_to_w`/`v_wind` so `refresh_rigid_aero!` sees fresh wind. STATIC twist surfaces. |
| **BODY_STATIC ride-point aerodynamic drag on the body** | ✅ **PARITY PASS** — Δvel=1.5e-10 with a draggy kcu (area=0.1, cd=1.0). `body_ride_drag_wrench`; per-point arm/va/drag torn (else mtkcompile stalls), distinct per-point params. Wired into `WingBodyVertex`; shared helper ready for plain `BodyVertex`. |
| **Winch + tether coexisting with a rigid wing body** | ✅ **PARITY PASS** — wing Δcom=4.5e-12, tether Δlen=3e-13, winch Δvel=2.6e-11 (`14a93869`). Winch/pulley/tether kernels widened to the body superset (`wide_pose_appendix`); getter `dyn_idxs` carries `(vertex,point)` pairs (they diverge once ride points compact vertex indices). |
| **Body-touching winched tether** (tether attaches at a BODY_STATIC point) | ✅ **PARITY PASS** — wing Δcom=3.8e-12, tether Δlen=4.2e-13 (`566faa1c`). `wrench_tether` edge (l0 = tether_len/n_segs) sharing one wrench core with the plain/pulley wrench edges; intra-body segments skipped (inert on a rigid body). |
| **DYNAMIC twist-surface DOF** | ✅ **PARITY PASS** — twist Δ=1.4e-16, twist_ω Δ=1.2e-14 (`8d9879db`). `WingBodyVertex` integrates `free_twist`/`twist_ω` states; `aero_moment` = `subsys.twist_moment[k]`, `inertia = mass·\|chord\|²/3` build constant. **`tether_moment=0`** (exact only where surface points carry no structural force). |
| Pulley steering on a body (both legs + kcu_steering to one body vertex) | ✅ **COMPILES** — `network_combined_wrench_segment` sums the parallel steering edges into one `EdgeModel` (task #28). Full SK100 structural (VSM rigid wing + twist + winch + 6-seg tether + 2 pulleys) **assembles and evaluates its RHS** (state 67, warm compile ~7.9 s, first-RHS f_ip JIT ~5.3 s, du finite). **Not yet run to parity** (next: full step vs monolith). |
| Twist `tether_moment` channel (per-point structural force → twist axis) | ⬜ (needed for DYNAMIC twist where surface points carry bridle/steering force; SK100 uses kinematic twist so lower priority) |
| Intra-body drag bridle (`le_*→kcu`, kcu not a twist member) | ⬜ (item 4 — needs a body-internal drag wrench mirroring `body_ride_drag_wrench`) |
| Flight-mechanics diagnostics (observed vars) | ⬜ |
| Knot-position getter for Hermite ride points | ⬜ (dynamics are correct; the diagnostic knot position isn't scattered back to the struct yet) |
| SK_coupled on ND (milestone) | ⬜ (structural layer now compiles; needs full-run parity + aero coupling) |

Everything not yet supported errors with a clear message in `build_body_mixed_network`
rather than silently mis-simulating.

### Roadmap to SK100 — structural layer now compiles

Pieces A (DYNAMIC twist), B (winch/tether + body), and the **combined multi-segment edge**
all landed this session. **The full SK100 structural model now assembles and evaluates its
RHS on the ND backend** — the last structural blocker is cleared.

**The combined multi-segment edge (ND parallel-edge limit) — DONE.** ND's `Network` uses
`Graphs.SimpleGraph`, which cannot hold two edges between the same vertex pair. Because all
of a rigid body's BODY_STATIC ride points collapse to a *single* body vertex, kite steering
creates parallel edges: pulley "left" has two legs (`le_steering_left`, `te_steering_left`)
**plus** `kcu_steering_left`, all three joining the body vertex to the one `steering_left`
point. Fixed: `build_body_mixed_network` accumulates `seg_pairs[key]`, and any pair with
>1 wrench-family segment becomes one `network_combined_wrench_segment` `EdgeModel`
(`build_combined_edges`, compiled per pair) that sums each member's wrench and emits once,
summing per-leg pulley tension to the pulley endpoint. Two subtleties that cost real time:
(1) **must tear the geometry, not just the summed outputs** — the first cut teared only the
force sums, leaving each member's nested pose→√→force expression inside the sum, and
mtkcompile ballooned to **29 GB RSS** and hung; the fix mirrors `network_dual_wrench_segment`
(tear anchor world pos `cxr_i` and forces `cff_i`/`cfr_i` per member). (2) A **one-line ND
fork bug** then surfaced at `Network()` assembly — `isdense` sorted `outidxs` twice instead
of `extidxs`, so the winch-tether `extin` indices (now interleaved by the combined edges)
failed the density assertion; fixed in `~/Code/SciML/NetworkDynamics/src/network_structure.jl`.

**Benchmark (rigid_structural_geometry.yaml, warm session):** state dim 67, compile
(`build_prob!`) ~7.9 s, first RHS eval (f_ip JIT) ~5.3 s, warm RHS ~52 µs, `du` finite. A
cold first-ever build adds ~500 s of one-time JIT of the mtkcompile/scalarization pipeline
(amortized across a session). **Next: run it to full-step parity against the monolith.**

Two smaller follow-ups for *exact* 2plate: **(item 4)** intra-body bridle segments touching
a non-twist point (`le_*→kcu`) carry monolith drag → need a body-internal drag wrench
(mirror `body_ride_drag_wrench`, `unit_stiffness=0` so only drag survives); and the twist
**`tether_moment` channel** (per-point structural force projected onto the twist axis) for
DYNAMIC twist where the surface points carry bridle/steering force — SK100 uses kinematic/
derived twist, so this is lower priority.

---

## 4. Files touched on this branch

- `src/components/components.jl` — wide `body_io`/`point_io`/`segment_io`,
  `finish_vertex`, `vertex_pose_io`, `segment_io_wide`, `BodyVertex`, `StaticBody`, the
  shared `timoshenko_element_wrench`/`elastic_joint_wrench`/`beam_hermite_ride_expressions`;
  `WingBodyVertex` (rigid VSM aero + DYNAMIC twist-surface DOF on a body) and
  `body_ride_drag_wrench` (BODY_STATIC ride-point drag, torn); `wide` flag +
  `wide_pose_appendix` on the point/winch/pulley builders (`WinchPoint`/`PulleyPoint`
  widened to the body superset).
- `src/generate_system/timoshenko_joint_eqs.jl`, `joint_eqs.jl` — delegate to the shared
  wrench helpers (monolith `test_timoshenko_joint` 20/20, `test_joint` unchanged).
- `ext/SymbolicAWEModelsNetworkDynamicsExt.jl` — `build_body_mixed_network` (BODY_STATIC
  absorption + graph remap; DYNAMIC + STATIC body vertices; **winch/pulley vertices +
  tether/pulley + wrench-tether/wrench-pulley edges; intra-body self-loop skip; `seg_pairs`
  accumulation → single vs combined edges**), `network_wrench_segment`,
  `network_wrench_tether_segment`, `network_wrench_pulley_segment` (shared
  `wrench_edge_loads`/`wrench_edge_emit_eqs`/`wrench_loads_from`), **`network_combined_wrench_segment`
  + `build_combined_edges` + `record_combined_edge_params!`** (sum parallel steering edges
  into one `EdgeModel`, geometry+force torn per member, pulley sign baked),
  `network_dual_wrench_segment`, `network_timoshenko_joint_edge`, `network_elastic_joint_edge`,
  wide edge kernels, per-orientation joint kernels, edge callables, `set_wing_twist_states!`,
  BODY_STATIC + body-pose + DYNAMIC-twist getter (`dyn_idxs` as `(vertex,point)` pairs).
- `~/Code/SciML/NetworkDynamics` (fork) — edgeless-network support in `construction.jl`; and
  the `isdense` `sort!(extidxs)` fix in `network_structure.jl:325` (was a dup `sort!(outidxs)`)
  that let the combined-edge + winch-tether `extin` model pass `Network()` assembly.
- `docs/src/private_functions.md` — new docstrings registered.

Validated by backend-swap scratchpad tests (see the session scratchpad): singleton body,
BODY_STATIC ride+wrench, Timoshenko/Elastic joints, body↔body segment, STATIC body, full
`test_yaml_bodies`, 3-body beam, nonlinear rigidity, and `test_timoshenko_joint` 20/20.
