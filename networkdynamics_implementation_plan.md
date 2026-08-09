# NetworkDynamics backend for SymbolicAWEModels — implementation plan

**Branch:** `networkdynamics-backend` (off `continuous-pressure-aero`)
**Worktree:** `/home/newbart/Code/Kite/SymbolicAWEModels-nd`
**Status:** working doc, not committed.

---

## ★ NORTH STAR

**The goal is NOT just to run SK100 — it is to make a COMPLETE alternative
backend with ALL of the features of the monolith.** `NetworkBackend` must reach
full feature parity with `MonolithBackend`.

**Acceptance test:** every monolith test passes unchanged except swapping the
`backend=` input (`SymbolicAWEModel(set, sys; backend=NetworkBackend())`). Run
the *same* existing tests on both backends; results must match to solver
tolerance (~1e-10 instantaneous RHS parity, per-point).

**Working rule — fix, don't block.** When a feature is missing or a structure is
unsupported on ND, implement it. Do not revert and declare a limitation. The
motivating win (compile time flat in N) is a bonus; parity is the requirement.

---

## 0. HARD REQUIREMENTS (non-negotiable)

Gate acceptance. Everything below conforms.

1. **No manual `@parameters` in the ext.** Every parameter flows through
   `ParamView`/`ParamRegistry` (`flat_params.jl`), struct-backed, live-synced into
   `NWParameter` each step (the ND analog of `sync_params!`). Mutating a struct
   field must take effect without a rebuild, exactly as on the monolith. ✅ met.
2. **No `pos1/pos2/pos3` scalarization.** Kernels use array vars (`pos(t)[1:3]`).
   Requires a **fork of NetworkDynamics** that scalarizes array-valued I/O across
   the whole codegen pipeline (`NetworkDynamicsMTKExt.jl`), not just input
   resolution. Dev fork at `~/.julia/dev/NetworkDynamics`. ✅ met.
3. **Nonnumeric (callable) params in codegen.** Fork threads callable params
   (aero `cl`/`cd`/`cm` polars) through f/g **and observed-function** codegen via
   a typed compartment (`pargs = (numparams, nonnumeric_params)`), so live-aero
   and observed formulas referencing polars resolve. ✅ met (branch
   `nonnumeric-params`: `ObsfunWrapper.pnn`, `VertexObsfunBatch.pnns`,
   `apply_compobs` dispatch).

Benign/narrow, noted: `total_mass` written into the struct as a build side-effect;
winch force assumes single tether segment at the winch point.

### Live-parameter design (req #1)

Reuses the monolith's `SAM.ParamGroup` + `sync_group!` — no new sync machinery.
Network-specific pieces: (a) setter is `setp(nw, …)` over composite
`VIndex`/`EIndex` (bare param symbols are ambiguous — one per instance);
(b) readers bound **per-instance** while walking points/segments
(`record_*_params!` + `ParamBuilder`), since one kernel compiles per *type* but
instance `i` reads `ss.points[i]`. Non-field-read values use small readers
(`GroundWindReader`, `TetherDragReader`, `PulleyLineMassReader`, `ConstReader`).
**GOTCHA:** `set_value` is excluded from the sync — owned by `set/get_set_values`,
synced-after-control would clobber the setpoint (cost a full-system divergence).

---

## 1. Why NetworkDynamics

`mtkcompile` + codegen builds **one** flattened RHS for all points (SK100: simplify
~63 s, ODEProblem/codegen ~141 s) — scales with total node count. ND compiles
**one kernel per component type** and iterates instances at runtime → compile time
flat in N. This is the performance motivation; per the north star, feature parity
is the actual requirement.

---

## 2. Progress at a glance

| Area | Status |
|---|---|
| Point/Segment components + shared kernels (`src/components/`) | ✅ landed & validated |
| Backend dispatch seam (`src/backends.jl`, `208d4ad`) | ✅ `backend::ModelBackend` field, `build_prob!` dispatch |
| ND ext build + full `init!`/`next_step!`/getter wiring | ✅ real model builds/inits/steps/integrates e2e |
| Point/segment parity | ✅ `max|Δpos|=0.0` over 4000 steps; full 2plate ram structure exact |
| Pulley (local vertex) + winch reeling (`extin`) | ✅ exact 0.0 RHS parity |
| Array-var + autonomous-param rewrite | ✅ req #1+#2; winch 3.2e-11, pulley 9.4e-12 |
| Body-frame damping | ✅ 2plate AeroNone RHS exact 1.6e-10 |
| Live particle aero (ContinuousAero/AeroPressure/AeroPlate) | ✅ t=0 parity ~1e-11 |
| Wing output observables (aero_force_b/moment_b via SII) | 🔨 fork+ext done, validating |
| BODY_STATIC/KINEMATIC ride-along points | ⬜ |
| ElasticJoints | ⬜ |
| Rigid bodies / beam (Timoshenko) dynamics | ⬜ |
| Rigid/pressure aero on beam-coupled wings (VSM) | ⬜ |
| Flight-mechanics diagnostics (elevation/azimuth/heading/course/aoa/turn) | ⬜ |
| SciMLBase-3 / MTK-11 migration | ✅ done (see §5) |

Legend: ✅ done · 🔨 in progress · ⬜ not started.

**Live-aero design** in `networkdynamics_aero_component.md`; **param unification**
in `networkdynamics_param_unification_plan.md`; **component unification** in
`component_unification_plan.md`.

### Current coverage

Runs the full ram-kite particle structure (points, segments, pulleys, reeling
winch, body-frame damping, live particle aero) end-to-end with exact instantaneous
parity to the monolith. Remaining monolith-only (the parity gap to close):
`BODY_STATIC`/`KINEMATIC` points, ElasticJoints, rigid bodies / beam dynamics,
VSM/rigid + pressure aero on beam-coupled wings, flight-mechanics diagnostics.
These are exactly the features SK_coupled (SK100 beam+pressure, in BeyondTheSim)
needs — so "SK_coupled runs on ND" is the concrete milestone that certifies the
structural half of the north star.

---

## 3. Design decisions (locked)

- **D1 — Components define boundaries; the graph is read, not discovered.** Each
  physical type is a small MTK `System`. Topology lives in `SystemStructure`
  (segments carry endpoint indices); assembly is a translator, not discovery.
- **D1a — Every distinct behavioral variant is its own component type.** Split
  axis = **state + connection signature**, not YAML section (ND fixes a kernel's
  state shape). Points split by dynamical variant (§7). Guard against
  over-splitting: same-state variants differing only by an additive term are one
  component with an optional input (e.g. `DYNAMIC` vs `WING PARTICLE_DYNAMICS`
  differ by an aero input → one component).
- **D1b — Aero consumes `TwistSurface` as a subcomponent** (twist angle, ω,
  geometry, flap δ), instead of the global-array handoff.
- **D2 — Single source of truth: both backends assemble from the same
  components/kernels.** The backend toggle selects only the assembly layer.
- **D3 — Shared physics, backend-specific aggregation sign.** MTK `Flow` sums to
  zero at a connect (monolith gathers `−Σ`); ND aggregation sums edge outputs
  into the input (`+Σ`). Kernels are identical; the sign is the one place the two
  wrappers differ (monolith edge emits negated force/half-mass; ND emits positive).
- **D4 — The frozen global-array monolith is the migration oracle** (parity to
  ~1e-10).
- **D5 — Migration is incremental.** Convert one component type at a time; the
  network refuses to switch on until the types it needs are componentized.
- **D6 — Keep full ModelingToolkit** (user call). ND components use
  `mtkcompile=true` on each per-type system (cheap, cached, flat in N), not ND's
  own reducer. ND stays an optional weakdep + extension.

---

## 4. Open questions

- **O1 — VSM AIC coupling:** N² edges (idiomatic) vs a single global-operator
  parameter refreshed per RHS eval (faster). Decide by benchmarking assembly cost
  in the aero phase. The live-particle-aero component already resolves this for
  the particle case (frozen-circulation, one shared kernel per wing, extin reads).
- **O2 — Constraint-index inventory:** does SK100 have any constraint above index
  1 that `mtkcompile` reduces? Determines whether ND's DAE-index-1 limit bites.
  Inventory during the beam/joint phase.

---

## 5. SciMLBase-3 migration (ND needs SciMLBase 3) — ✅ DONE

ND requires **SciMLBase 3**; the stack was on 2. Fix (maintainer approved) keeps
**MTK at v11** (has SciMLBase-3 builds, 11.35+) and moved the solver stack up:
SciMLBase 2.155 → 3.34, MTK 11.26 → 11.37.1, OrdinaryDiffEqCore 3.33 → 4.6,
DiffEqBase 6 → 7.6.1; ND 1.1.0 + Graphs 1.14 added; VSM updated to v4 (adds
`section_aero`, `ObjAdapter`, `NeuralFoil`). The monolith still works on the new
stack (identical correct RHS). **KiteModels** caps DiffEqBase at 6.x and was
removed from `examples/Project.toml` — any example that `using KiteModels` is
broken until KiteModels is updated upstream.

---

## 6. Component model — type → ND construct

| SymAWE type | ND construct | Notes |
|---|---|---|
| `Point` | `VertexModel` (several types, §7) | state pos,vel; out pos,vel; in ΣF, Σmass |
| `Segment` | `EdgeModel`, one type | reads both endpoints, writes ±F + half-mass; stateless |
| `Pulley` | stateful `VertexModel` | rope across ≥2 segments → vertex owning `pulley_len` |
| `Winch` | `VertexModel` | tether-length ODE; segments read `tether_len` via `extin` |
| `Body`/rigid `Wing` | heterogeneous `VertexModel` | pose state; drives attach points; force/moment accumulators → edge aggregation |
| `TwistSurface` | `VertexModel` | deformable section; twist dynamics; aero subcomponent (D1b) |
| VSM aero | dense `EdgeModel` **or** global-operator param | O1 |

### Connector & aggregation (validated, §3-D3)

- Point outputs both `pos` and `vel` (segments need endpoint velocity).
- Mass is a second aggregation channel: segment writes half-mass into a Flow var;
  point translational mass = `extra_mass + Σ half_mass` (dynamic `l0` → dynamic
  half-mass cleanly).
- Shared kernels `src/components/kernels.jl`: `segment_geometry`,
  `segment_spring_force`, `segment_perp_drag`, `point_drag_force`,
  `segment_half_mass` — geometry+material in, forces out; environment injected by
  the caller. Both backends call these via `point_acceleration` /
  `segment_endpoint_loads` in `components.jl`.

---

## 7. Point vertex taxonomy (refines D1a)

The one-big-`point_eqs!` splits into vertex types by state + connection signature:

| Component | states | slaved to | force it returns | status |
|---|---|---|---|---|
| `StaticPoint` | 0 | — (`pos ~ pos_w`) | — | ✅ |
| `DynamicPoint` | 6 | — | aero = optional input (0 for non-wing) | ✅ |
| `BodyRidePoint` (BODY_STATIC) | 0 | one body pose | force + COM moment → body | ⬜ |
| `BeamPoint` (KINEMATIC) | 0 | a Timoshenko joint (2 bodies) | axial-split force + moment → both ends | ⬜ |
| `RigidWingPoint` | 0 | one body + its twist-surface | force + moment → body | ⬜ |

**Slave points are NOT top-level vertices (§8.2).** `BodyRidePoint`/`BeamPoint`/
`RigidWingPoint` have 0 own states and forward an aggregated force — that is
feedthrough, forbidden as a plain ND vertex. They are absorbed into their **cluster
vertex** as outputs (positions) + internal wrench accumulation; only `StaticPoint`
and `DynamicPoint` are stand-alone point vertices.

- **Merge** `DYNAMIC` + `WING-particle` → one `DynamicPoint` (identical 6-state
  integrator; aero is an optional input). Done.
- **Do not merge** slave types with `DynamicPoint` (0 vs 6 states → different kernel).
- **Edges:** `SpringDamperSegment` is one type. Wing-structural special cases
  (rigid→zero spring, particle→no drag) fold into params (`unit_stiffness=0`, a
  drag flag). `l0` (pulley/tether/fixed) is an **input** from the winch/pulley
  vertex, not a branch inside the segment.

---

## 8. ND coreloop constraint → the component-granularity design (pure ODE, no DAE)

### 8.1 The constraint

The ND coreloop (`NetworkDynamics/src/coreloop.jl`) runs **non-feedthrough vertex-g
first** (line 39, `(nothing,nothing)`), then edge-g (41), then feedthrough vertex-g
(55, with `aggbuf`), then `collect_externals!` (61), then the f-passes (76/97, with
`aggbuf`+`extbuf`). Therefore:

- A **non-ff vertex output depends only on its own states** — not on aggregated edge
  inputs and not on `extin`. **`extin` is available only in the f-pass.**
- `construction.jl:54-80`: **general feedthrough vertices are forbidden** (allowed
  only as single-neighbour leaf + one `LoopbackConnection`); `ff_to_constraint(vf)`
  turns feedforward outputs into **algebraic DAE states** — **rejected** (DAE hurts
  perf via a bigger/worse-conditioned Newton solve + consistent init; user vetoed it).

This is why the monolith's `body_force`/`body_moment` `.+=` accumulators do **not**
translate to a plain vertex sum: a ride point forwarding its *aggregated* segment
force back to a body is inherently feedthrough.

### 8.2 The resolution — vertex granularity is a design choice

The monolith is also one flat symbolic system, but `mtkcompile` does **whole-system
tearing** — it topologically sorts one big DAG, so a *stacked* relation ("this point
rides that body and pushes force back") schedules fine no matter how many struct
boundaries the chain crosses. ND has **no tearing and no hierarchy**: a flat graph
with a hard-coded schedule that knows exactly one dependency shape (vertex ← its
edges ← neighbour vertices). A point that is really a *sub-part of a body* has no
honest home as its own top-level vertex.

So don't give it one. **A vertex = a maximal cluster of state-carriers whose ride
points would otherwise be feedthrough.** Absorb the feedthrough into the vertex:

- the cluster's ride-point **positions/velocities are outputs** — pure functions of
  the cluster's own states (non-ff ✓);
- **segments stay ordinary edges** (cluster ↔ free knot); at edge-g each reads the
  ride position off the cluster output and delivers a **full wrench** — force **and**
  moment-about-COM (it knows the arm: ride-pos and COM are both cluster outputs) —
  into the cluster's aggregated input;
- the cluster **f-pass** sums the aggregated wrenches + internal elastic + aero +
  gravity and integrates its states.

Feedthrough-free, pure ODE — no DAE, no recompute, no separate ride vertices. The
force arrives as a *normal aggregated edge input*, exactly what ND is built for.

**Partitioning rule (union-find over bodies) — merge only where forced.** Merge two
bodies iff a ride point depends on **both** — i.e. a Hermite point *rides a
`TimoshenkoJoint` between them*. The trigger is a point on the joint's centerline,
nothing else: a **bare joint** (couples two poses, no point riding it) stays an
**edge**, and its two bodies stay **separate singleton vertices**. The merge is
transitive (a node shared by two ride-point relations can't be split across vertices
without duplicating its states + an equality constraint = DAE), so each *connected
group* = one cluster vertex — but that is the minimal correct grouping, not
over-unification. Common case (plates on bare joints, bridle lines to free knots):
**every body is a singleton vertex, every joint an edge.** Only a true bridle-on-beam
collapses a chain.

### 8.3 The three cases

- **Rigid body + its `BODY_STATIC` points** → singleton cluster vertex (13 states:
  `com_w`,`com_vel`,`Q_p_to_w`,`ω_p`). Ride positions = body-pose functions of the 13
  states (`body_ride_eqs`, non-ff). Segments = edges; wrench back.
- **Flexible Timoshenko beam** (chain of *N* node-bodies + *M* internal elements +
  Hermite ride points) → **one** cluster vertex (13·*N* states). The internal
  Timoshenko elements are **internal equations** (not edges): stateless corotational
  wrenches between consecutive node states (`timoshenko_joint_eqs!`). Hermite ride
  positions are functions of the node states (`beam_hermite_ride_eqs`, non-ff). A
  beam's flexibility **is** its node-body DOF (see 8.4).
- **Inter-cluster joint** (`ElasticJoint` / `TimoshenkoJoint` between two *separate*
  clusters, no shared ride point) → an **EDGE**. Stateless: reads both endpoint
  poses/vels (cluster outputs), computes the equal-and-opposite wrench, delivers into
  both clusters' aggregated inputs. Pure ODE.

**Parallel cluster↔cluster segments** (two bridle lines between the same two plates)
→ merge into one multi-segment edge kernel — fatter edge, **no multigraph needed**.

### 8.4 Timoshenko joints have no DOF of their own

`timoshenko_joint_eqs!` and `joint_eqs!` are **stateless** — pure feedforward from
the two nodes' `(pos, R, com, vel, ω)` to a restoring wrench (corotational frame →
small per-node deformations → consistent `EA/GJ/EIy/EIz` + shear-`Φ` stiffness →
wrench transported to each COM). Nothing is integrated. **The DOF live in the
node-bodies**; the joint only *couples* them. That is exactly why an inter-cluster
element becomes an **edge** (stateless coupling of two vertices' outputs) while an
intra-beam element becomes an **internal equation** of the beam's single vertex.

Shared math is not duplicated: `rigid_body_pose_expressions` (`components.jl`)
sources both backends (`rigid_body_eqs!` delegates — monolith 19/19 green); the
`joint_eqs!` / `timoshenko_joint_eqs!` wrench bodies become the edge-kernel /
internal-equation bodies verbatim.

### 8.5 ND uniform-output-depth constraint (widen + zero-fill, per-model)

`construction.jl:99-105` enforces **one uniform vertex-output width (`vdepth`) and
one uniform edge-output width (`edepth`)** across the whole `Network` — "All vertex
models must have the same output dimension." *States* may differ per vertex
(`dynstates = Σ dim`, so a body's 13 and a point's 6 coexist), but the *interface* is
uniform. A body cannot simply output a wider pose than a point.

Resolution — **widen the common interface to a superset and zero-fill, chosen
per-model**:
- Vertex output superset = `[pos(3), vel(3), pose-extras(R 9, com 3, com_vel 3, ω 3),
  pulley_len(1)]`; a point fills `pos`,`vel` (+`pulley_len`) and zeros the pose slots;
  a body fills the pose slots.
- Edge output superset (wrench) = `[force(3), mass(1), tension(1), moment(3)]`; a
  point↔point segment fills `force`/`mass`/`tension`, zeros `moment`; a point↔body
  segment fills `force`/`moment`, zeros `mass`/`tension`.
- Edge *inputs* (what a kernel reads) may differ per edge model — point↔point reads
  pos/vel both ends; point↔body reads pos/vel from the point and the pose slots from
  the body. Only the *output* widths are constrained.

**Per-model depth (no cost on SK100):** `vdepth`/`edepth` are build-time choices.
Pure-particle models (SK100, thousands of points, no bodies) keep the current narrow
depth — zero regression. Only body-containing models (RAM kite, ~15–100 points, a few
bodies) pay the wider interface, which is negligible at that scale. The **bare body**
(`test_rigid_body.jl`, no points/segments) needs no widening at all — first target.

**Validation ladder:** bare body (`test_rigid_body.jl`: free-fall, spin-up,
COM-offset) → body + one `BODY_STATIC` ride point + one segment (first use of
widen+zero-fill + heterogeneous edge) → Timoshenko beam cluster → inter-cluster
joints → rigid VSM aero.

## 9. Roadmap to the north star

Structural-parity work, in dependency order. Each item is validated by
backend-swap parity against the monolith on the existing tests (§NORTH STAR).

1. **Wing output observables** ✅ — `aero_force_b`/`aero_moment_b` on ND via SII
   observed vars. Fork threads nonnumeric params into observed codegen (req #3); ext
   adds aggregate observed vars + per-wing readers. `test_continuous_aero.jl` **19/19
   on NetworkBackend** (backend-swap only).
2. **Cluster-vertex builder + partitioning** 🔨 — union-find bodies over shared
   Hermite ride points (8.2); emit one `VertexModel` per cluster with wrench-input
   aggregation. Rigid-body case first (singleton cluster), validated against
   `test_wing_dynamics.jl` / `test_rigid_body.jl` by backend-swap.
3. **`BODY_STATIC` ride points as cluster outputs + wrench edges** ⬜ —
   `body_ride_eqs` positions as vertex outputs; segments as wrench-carrying edges;
   merge parallel cluster↔cluster segments.
4. **Timoshenko beam as one cluster vertex** ⬜ — internal `timoshenko_joint_eqs!`
   equations + `beam_hermite_ride_eqs` Hermite outputs inside the cluster; validated
   on the beam tests.
5. **Inter-cluster joints as edges** ⬜ — `joint_eqs!` / `timoshenko_joint_eqs!`
   wrench as `EdgeModel` between two cluster vertices.
6. **Rigid/pressure aero on beam-coupled wings (VSM)** ⬜ — `TwistSurface` as its
   own component (D1b); resolve O1; per-section Γ solve as index-1 rows or an inner
   solve injected as a param.
7. **Flight-mechanics diagnostics** ⬜ — elevation/azimuth/heading/course/aoa/
   turn_rate as observed vars (KiteUtils.calc_* in `scalar_eqs.jl`), materialized by
   the getter each step, matching `model_management.jl`.
8. **Milestone — SK_coupled on ND** ⬜ — SK100 beam+pressure builds/inits/steps on
   `NetworkBackend`, certifying the structural half.
9. **Consolidate** ⬜ — unify the getter/observable API across backends (SII).

---

## 9. Constraints to keep in view

- Re-authoring the monolith via components does **not** speed it up (`mtkcompile`
  flattens back to a monolith). The flat-in-N win is **network-only**; the monolith
  rewrite buys zero drift, not speed.
- ND is ODE + index-1 mass-matrix DAE only — no global Pantelides/dummy-derivative.
  The particle-spring net needs no reduction, but any higher-index constraint
  SymAWE relies on `mtkcompile` to reduce must be reformulated by hand (O2).

---

## 10. References

- ND MTK integration: https://juliadynamics.github.io/NetworkDynamics.jl/stable/mtk_integration/
- `VertexModel(sys, inputs, outputs; vidx)`, `EdgeModel(sys, srcin, dstin, srcout,
  dstout; src, dst)`, `with_mtk_model_cache() do … end`.
- SymAWE seam: `create_sys!` (`src/generate_system/create_sys.jl`),
  `maybe_create_prob!`/`build_prob!` (`src/model_management.jl`).
