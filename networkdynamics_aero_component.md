# Aero as a NetworkDynamics component

**Status:** design (not yet implemented). Working doc, not committed — like
`networkdynamics_implementation_plan.md`. This is the trickiest step of the
network backend, so it gets its own note. Everything below is grounded in how the
monolith already assembles aero (`src/generate_system/aero_eqs.jl`, the
`VSMEngine` in `src/system_structure/types.jl`, `refresh_aero!`) and in what the
structural network backend already does and validates today.

---

## 1. Why aero looks hard, and why it isn't as hard as it looks

The Vortex Step Method solves a circulation distribution whose aerodynamic
influence coefficient (AIC) matrix is **dense**: every panel's induced velocity
depends on every panel's circulation. Taken literally, coupling that into a
nearest-neighbour graph integrator (NetworkDynamics couples only along edges)
would need an all-to-all block — the `O(N²)` edges of open question **O1** in the
implementation plan.

The escape hatch already exists in this codebase: **the VSM circulation is not
re-solved inside the ODE RHS.** It is *frozen and linearised* around an operating
point and refreshed between steps. Concretely, `VSMEngine`
(`types.jl:151`) carries

- `aero_y` — the operating-point inputs `[alpha, beta, ω₁, ω₂, ω₃, twist…]`,
- `aero_x` — the baseline wind-axis coefficients `[CL, CD, CS, CM…]`,
- `aero_jac` — the dense Jacobian `d(aero_x)/d(aero_y)`.

Inside the RHS the aero force is the **affine** map
`aero_x ≈ aero_x₀ + aero_jac · (aero_y − aero_y₀)`, evaluated symbolically, then
mapped to per-panel/per-point forces. The expensive dense solve happens in
`refresh_aero!` every `vsm_interval` steps and only updates `aero_x₀`, `aero_y₀`
and `aero_jac`.

**Consequence for the network backend:** the AIC never becomes graph edges. It
lives entirely inside the frozen `aero_jac`, which is a **parameter** — exactly
the kind of value the backend already refreshes between steps (`sync_params!` on
the monolith, `NWParameter` mutation on the network). O1 resolves in favour of
the "global-operator parameter" option, and it does so *for free* because the
monolith is already structured this way. No `O(N²)` edges.

So the aero component is: **a per-wing affine force law over a small, fixed
operating point, with a dense Jacobian parameter refreshed between steps.**

---

## 2. What the aero subsystem reads and writes (from `aero_eqs.jl`)

The monolith already wires each wing's aero as an isolated subsystem
`aero_component(wing.aero, wing, sys_struct)` with a clean connector interface —
this is the same `aero_component` dispatch the plan's D1b relies on. Its I/O is
exactly what a network component needs:

### Particle wing (`PARTICLE_DYNAMICS`, e.g. the 2plate/SK100 kite)
Per aero-surface node `k` (`aero_eqs.jl:45-56`):

- **in:** `point_pos[:,k]`, `point_vel[:,k]` (body frame, relative to the wing
  origin), `va[:,k]` (apparent wind, body frame), `rho[k]`; plus per-panel flap
  `delta[i]`.
- **out:** `point_force[:,k]` — the aero force lumped onto that structural node.

So the subsystem is a map **N node-states → N node-forces**. That is a
one-component-per-wing hub reading its nodes and writing back per-node forces.

### Rigid wing (`RIGID_DYNAMICS`)
Whole-wing (`aero_eqs.jl:70-90`):

- **in:** `va` (body-frame apparent wind), `rho`, `R_b_w`, `omega`,
  `twist`/`twist_vel` per twist-surface.
- **out:** `force`, `moment` (the wing wrench, body frame) and a per-twist-surface
  aero moment.

So the subsystem is a map **wing pose/rate + twist state → one wrench (+ twist
moments)**.

Either way the operating point is *small and fixed-size* (a handful of scalars
plus one entry per twist-surface), independent of node count. This is why the
frozen-Jacobian parameter is tiny and the component compiles flat in N, like the
point/segment kernels already do.

---

## 3. Graph topology: a wing "aero hub" vertex + N spoke edges

The structural backend today makes every point a vertex and every segment an
edge, with the aggregation sign convention already validated (D3). Aero adds one
construct per wing:

```
            (wing-node vertices)                      structural
   pt_a ────────┐                                     segments stay
   pt_b ────────┤   spoke edges (state → hub,         ordinary edges
   pt_c ────────┼──▶  force ← hub)          ┌───────────────────────┐
   pt_d ────────┘                           │  aero hub vertex      │
                                            │  • fits wing frame    │
                                            │  • forms aero_y       │
                                            │  • affine aero_x      │
                                            │  • frozen aero_jac(p) │
                                            └───────────────────────┘
```

- **Aero hub vertex** (one per wing). Owns the frozen VSM operator as parameters
  (`aero_x₀`, `aero_y₀`, `aero_jac`, panel geometry, `point_to_vsm_point`,
  `wing_segments`). For a rigid wing it *is* the wing body vertex (Phase 4); for a
  particle wing it is a stateless hub whose only job is the aero force law.
- **Spoke edges** (one per aero-surface node): each carries the node's
  position/velocity **to** the hub (needed to fit the particle-wing frame and
  build `aero_y`) and carries that node's aero force **back** to the node vertex's
  force input. This is a *single ordinary ND edge* — src→dst and dst→src outputs
  are the standard `EdgeModel` signature the segment already uses. **N spokes per
  wing, not N².**

The dense AIC coupling between nodes is *not* represented as node-to-node edges;
it is summarised in the hub's `aero_jac` parameter. The spokes only move
per-node kinematics in and per-node forces out.

Twist-surfaces (D1b) become their own vertices carrying the twist DOF
(`twist_angle`, `twist_ω`, or a KINEMATIC flap δ) and connect to the hub by the
same spoke pattern, so the hub reads twist state instead of re-deriving it
through the current global-array handoff.

### Why the hub can fit the wing frame locally
The particle-wing frame is fitted from its node positions (`R_v_to_w`, `wing_pos`
in `create_sys.jl`). Because every spoke delivers its node's position to the hub,
the hub has all the inputs it needs to fit the frame and form `alpha/beta/ω` —
this is a many→one aggregation, exactly what an ND vertex does with its incident
edges. No global array, no non-local read.

---

## 4. Between-step refresh (the `vsm_interval` mechanism)

`next_step!` already re-linearises the VSM every `vsm_interval` steps via
`refresh_aero!`, then re-syncs parameters. On the network backend this becomes:

1. step the ND integrator `dt`;
2. scatter node/wing state back to `sys_struct` (the network `get_all_state`
   getter already does this for points);
3. if `iter % vsm_interval == 0`: `refresh_aero!` re-solves the VSM at the new
   operating point and updates the `VSMEngine` fields **in place** (mutable
   struct);
4. push the updated `aero_x₀/aero_y₀/aero_jac` into the hub vertex's `NWParameter`
   entries (the network analogue of `sync_params!`).

Nothing here is new physics — it is the existing refresh, retargeted from the
monolith's flat parameter buffer to the hub vertex's parameters. The frozen model
between refreshes is what keeps the RHS an affine, edge-local computation.

---

## 5. Shared force law (D2), same as the structural kernels

The per-node / per-wing force computation is factored into a shared helper the
same way `point_acceleration` / `segment_endpoint_loads` are shared today
(`src/components/`). Both backends call it:

- monolith: inside the `aero_component` subsystem wired by `aero_eqs!`;
- network: inside the aero hub vertex.

The environment (`rho`, apparent wind) is injected by the caller, never read from
the atmosphere inside the kernel — identical to the structural kernels, so there
is one source of truth for the aero force and no drift between backends. The only
backend-specific piece is aggregation sign on the spoke edges (D3), which is
already handled for the structural edges.

---

## 6. Validation plan (mirrors what the structural network already passed)

The structural backend was validated by exact instantaneous parity against the
monolith on the real 2plate ram model (tether points to 0.0, wing points to
~1e-3; the only divergence is the deliberately-frozen pulley/winch length). Aero
gets the same treatment:

1. **Frozen-operating-point RHS parity.** With the VSM linearisation held fixed
   (same `aero_jac`, `aero_x₀`), assert the network's per-node aero force equals
   the monolith's `aero_force_point_b` (particle) / wrench (rigid) to ~1e-10 at
   t=0, on the 2plate model with real aero instead of `AeroNone`.
2. **Refresh parity.** Run several `vsm_interval` refreshes and confirm the hub's
   parameter update reproduces the monolith's `refresh_aero!` result.
3. **Trajectory parity.** With pulley/winch also handled (or frozen on both
   sides), integrate and confirm the wing pose tracks the monolith.

Because the AIC is a parameter and the operating point is small and fixed-size,
each per-wing aero kernel compiles once and runs flat in node count — the whole
point of the backend.

---

## 6c. IMPLEMENTED design — live particle aero in the vertex f-pass (supersedes §6b)

**§6b's aero-hub-vertex design is invalid** and was NOT used. It assumed a stateless
hub vertex whose *output* `g` reads the wing points via `extin`. The ND coreloop
(`coreloop.jl`) collects `extin` (`collect_externals!`) **after** every vertex-`g` pass
(lines 39/55) and before the edge-`f`/vertex-`f` passes (76/97). So **`extin` is legal
only in an f-pass** — a vertex `g` cannot read it. A hub-`g` reading wing points is
impossible.

**What shipped instead:** the aero force is computed in the **wing node's own
f-pass** (`SAM.LiveAeroWingNodePoint`, `components.jl`). Each PARTICLE wing whose mode
computes live force (`needs_live_aero` = `is_wing && !provides_aero_override`; true for
ContinuousAero/AeroPressure/AeroPlate, false for AeroNone/AeroDirect) gets **one shared
kernel per wing**:

- It reads via `extin`: the ref-point frame + origin velocity (`body_damp_extin`), the
  origin position, and **every** wing point's position/velocity (`live_aero_extin`).
- It fits `Rbw`/`wing_pos` symbolically and instantiates the **unchanged** shared
  `aero_component` as the `aero` subsystem, wiring its `point_pos`/`point_vel`/`va`/`rho`
  connectors exactly as `aero_eqs!` does (`live_aero_connector_eqs`).
- It selects **its own** point force from `aero_component`'s output by the per-instance
  `aero_slot` parameter, rotates to world, and adds `force/mass` to `D(vel)` on top of
  the shared point acceleration + body-frame damping — matching the monolith's
  `aero_force_w = R_b_to_w·aero_force_point_b` (`point_eqs.jl`).

Every live-aero node of a wing shares this kernel and the same `extin` (they differ only
in `aero_slot`), so it compiles **once per wing** — aligned with the time-to-first-solve
goal. The whole wing's aero is evaluated in each node (`O(npts)` per node, `O(npts²)`
per step); the expensive VSM solve stays out of the RHS (frozen `v_ind`, refreshed every
`vsm_interval`). Trade-off accepted for one-compile flatness; can be optimised later.

**Params (the callable work paid off):** the aero subsystem uses a **full-path** param
registry (`ParamView` default / monolith-style names) so a mode reading several
twist-surfaces in one kernel (AeroPlate's `chord`/`area`) does not collide under the
network's per-field memoisation. After nesting they are `aero₊p_wings_i_aero_v_ind_c_j`
(matrix, scalarised column-major), `aero₊p_wings_i_aero_cl/cd/cm` (callables). The
scalar/matrix params sync through the `Float64` `ParamGroup`; the callables through the
ND fork's nonnumeric SII route (`record_aero_params!`, a `MultiCallableSetter`). The ND
fork's `_scal_flat` was fixed to `vec`-flatten matrix params (`v_ind` is 3×n_panels).

**Validated (t=0 instantaneous acceleration parity, 2plate BILLOWING):**
ContinuousAero max |Δacc| = **1.04e-11** across all 6 wing points vs the monolith
(compared at the post-init state, where pulley/winch states are identical). AeroDirect
(frozen) = 1.05e-11 and AeroNone finite — no regression. AeroPressure routes through the
identical kernel but its mode is env-blocked here (VSM `Panel` has no `section_aero` in
the installed VortexStepMethod, fails in the monolith too). AeroPlate needs a
per-point twist-surface fixture (no test coverage exists) but is structurally handled by
the full-path param naming. AeroLinearized is RIGID-only (out of scope until rigid wings
land on the network).

## 6b. Original (SUPERSEDED) spec — aero-hub vertex + extin (design only, not built)

Grounded in verified code: the shared `aero_component` I/O (`aero_eqs.jl:45-66`),
`build_panel_force_eqs` (frame-**covariant** — `common.jl:280`), the callable
`ContinuousPolar` params, and the **already-working** ref-point frame fit the
network does for body damping (`body_damp_extin`, ext:477). This is the D5-faithful
path: reuse the *unchanged* `aero_component`; only the assembly differs.

### Architecture: one stateless aero-hub vertex per wing + pure `extin` (no spokes)
The ND aggregator sums incident-edge contributions, so it can't hand a vertex its
neighbours *individually* — but the aero law needs each strut point separately.
Solution (matches the existing backend, which is `extin`-heavy, not edge-heavy):

- **Aero-hub vertex** (one per wing, `dim=0`, isolated in the graph — a
  `StaticVertex` whose `g` is the aero force law):
  - `extin` reads, by network index:
    - the wing's z/y ref points' `pos_k` and origin `vel_k` (same list as
      `body_damp_extin`) → build `Rbw`, `wing_pos`, `wing_vel` **symbolically**,
      reusing the frame-fit already in `wing_eqs.jl`/body-damp;
    - every wing point's `pos_k`, `vel_k` (the strut LE/TE points).
  - Per wing point `k` it forms the shared connectors exactly as `aero_eqs.jl` does:
    `point_pos_b = Rbw'*(pos_k − wing_pos)`, `point_vel_b = Rbw'*vel_k`,
    `va_b = Rbw'*(wind(z_k) − vel_k)`, `rho_k = calc_rho(z_k)`  (`z_k` = `pos_k[3]`;
    `wind_vec` passed verbatim as today).
  - Instantiates `aero_component(wing.aero, wing, ss; params)` **unchanged**, wires
    the connectors, and the hub **output** (`outsym`, dim `3·n_points`) is each
    point's force rotated back to world: `Rbw * point_force_b[:,k]`.
- **Wing-node vertices** gain one `extin` bundle: `VIndex(hub, :aero_force_k_c)`
  (c=1..3), added into their force balance. ND computes the hub's `g` (output) in
  the g-pass, then the wing nodes read it in the f-pass — no algebraic loop
  (both are functions of states), no new edges.

### Params on the hub (this is what the nonnumeric work unblocked)
`aero_component(::ContinuousAero)` declares `[vind_p, cl, cd, cm]`:
- `vind_p` — `make_array_param` (3×n_panels frozen induced velocity), body frame;
- `cl,cd,cm` — the `ContinuousPolar` **callables**, now carried by ND and set
  through the uniform SII route:
  `record`: `add_callable_param!(builder, VIndex(hub,:cl), reader=…wing.aero.cl)`.
  `sync`: a callable-typed group calling `setp(nw, VIndex(hub,:cl))(prob, polar)`
  (NOT the Float64 `ParamGroup` buffer — callables don't change per step; re-point
  only on refresh/deserialize). Add a `:callable` kind to `record_*_params!` /
  `build_network_param_sync`, mirroring `:scalar`/`:array`.

### Between-step refresh (hub params)
Same as monolith `refresh_particle_aero!`, retargeted: after scatter, if
`iter%vsm_interval==0` run the VSM solve, freeze `v_ind` (array param resync) and
(if the polars' `body_aero` object identity changed) re-point `cl/cd/cm` via
`set_nonnumeric!`. If `ContinuousPolar` wraps a *mutated-in-place* `body_aero`, the
object identity is stable and no callable resync is needed mid-run — verify.

### AeroPressure = ContinuousAero + δ
`aero_component(::AeroPressure)` adds `delta[1:n_panels]` (per-panel flap) and a
traction scatter. Network variant is identical to the hub above plus: the hub reads
each panel's twist-surface δ via `extin` `VIndex(twist_surface_vertex, :twist_delta)`
(twist surfaces become their own vertices, doc §3), wiring `subsys.delta[i]`. The
`build_panel_force_eqs` δ path (`common.jl:316`) and 3-arg polars are already shared.

### t=0 parity gate (before trajectory)
Build 2plate (BILLOWING) on both backends, freeze the same `v_ind`, assert the
network per-point `aero_force_point` == monolith `aero_force_point_b` to ~1e-9 at
t=0 (post-refresh Rbw identical on both sides). Then refresh parity, then trajectory.

**Status:** spec complete + verified against source; ND callable infra committed
(`bcf03fc`). Coding blocked only on memory for the build/test loop — implement the
hub builder + wing-node aero `extin` + `:callable` param sync against the t=0 gate.

## 7. Open sub-questions

- **Frame-fit inside a vertex.** Fitting `R_v_to_w` from node positions inside the
  hub vertex must stay allocation-light; the monolith already does this per step,
  so it is a port, not new work.
- **Particle vs rigid share of the hub.** For rigid wings the hub coincides with
  the wing body vertex (Phase 4); for particle wings it is a standalone hub. Both
  use the same aero kernel; only the output wiring differs (wrench to body vs
  per-node forces to nodes). Keep one `aero_component`-style dispatch, as the
  monolith does.
- **δ / flap.** KINEMATIC twist-surfaces feed `delta` into the panel map
  (`aero_eqs.jl:38-43`); the twist-surface vertices supply this over their spokes,
  no special path.
