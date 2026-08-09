# Component unification (D5): one shared component source, both backends

**Branch:** `networkdynamics-backend`. Working doc, not committed. Builds on
[[networkdynamics_param_unification]], [[networkdynamics_backend_status]].

## The one rule

**Every piece of physics is a shared component in `src/components/`. The ONLY thing
that differs between backends is how those components are assembled.** No physics
lives in `point_eqs.jl`/`segment_eqs.jl` at the end — those files are retired branch
by branch. The network assembles the components as an ND graph; the monolith
assembles the identical components with explicit aggregation equations and
`mtkcompile`. Same components, same params, same defaults.

## Why everything must be a component (the graph benefit)

The NetworkDynamics backend's win — SK100 time-to-first-solve — comes from the system
being a **graph**: vertices with *bounded* coupling, so ND exploits sparsity/locality.
If any physics stays as monolithic global equations:
- on the network it can't be represented graph-locally at all — it either doesn't
  exist in the backend, or breaks the sparsity that makes ND fast;
- and the two backends stay duplicated for exactly the hard parts.

So keeping `point_eqs.jl` for wing/twist/aero is at most a temporary scaffold, never
the endpoint. A hybrid monolith is a transition state, not a design.

## The nuance: "component" means "graph-LOCAL component"

Making something a component only helps if its coupling is **bounded**. A rigid wing
that needed every point to read every other point is not local and buys no speedup.
Here the couplings ARE bounded, so this is achievable:
- **wing frame** → each wing node reads its 4 ref points via `extin` (network already
  does this);
- **rigid body** → ONE vertex owning com pos/vel/quat/ω, points derived from it
  (task #12) — not N mutually-coupled vertices;
- **aero** → the frozen-circulation scheme keeps each point's force local (refreshed
  between steps).

If a coupling truly can't be made local, forcing it into a component is theater. The
design work per subsystem IS finding its local decomposition.

## The reframe: promote network kernels → shared components

The network's `network_*` kernels ALREADY are the graph-local decomposition
(vertex/edge with bounded `extin`). So the move is not "monolith → new components";
it is: **promote each network kernel into a shared `src` component, then have the
monolith assemble from it**, retiring the matching `*_eqs.jl` branch. The graph
structure is the design; the monolith just picks a non-graph assembly of the same
pieces.

## Shared I/O is the coupling interface, used by BOTH backends

The segment carries the full graph I/O (`pulley_len`, `tension`, endpoint loads) —
NOT a burden on the monolith but the shared coupling interface: the network wires it
through the graph, the monolith wires it with an explicit equation
(`seg.pulley_len_in ~ pulley_len`). Points already carry the superset via `point_io`;
segments get the same treatment.

## Approach B (settled): explicit force-input, supersede `Node`/`Flow`/`connect`

A `Flow` connector can only be consumed by MTK `connect`; the network never uses
`connect` (it sums edge outputs into vertex inputs over the ND graph). The only
interface both an ND-graph and an MTK-`connect`-free assembly can share is a plain
**input variable**, which also preserves array I/O (`pos(t)[1:3]`) that `Flow` forces
to scalarize. `Node`/`Flow` superseded. Reverses committed `33bb307e`.

## Backend assembly (the ONE place they differ)

- **Network** (ext): wrap each component as `VertexModel`/`EdgeModel`; ND sums edge
  `*_force`/`*_mass` outputs into the vertex `force_in`/`mass_in`. Naming: one generic
  symbol per component TYPE (`:extra_mass`), addressed by `VIndex`/`EIndex`.
- **Monolith**: `@named point_i = DynamicPoint(s, params, i)` etc.; add explicit
  aggregation equations `point_i.force_in ~ Σ incident seg *_force`, one
  `System(..., systems=[all components])`, `mtkcompile`. Naming: generic `:extra_mass`
  → MTK-namespaced `point_i₊extra_mass`.

## Sequence — REVISED (2026-08-03, user decision: defer the monolith to last)

The monolith gets **no** perf benefit (still one dense ODE) and a phased monolith
migration needs a fragile alias bridge (`pos[:,i] ~ point_i₊pos`) coexisting with the
flat path on the *reference oracle*. Also the shared `DynamicPoint` is a strict subset
of the monolith DYNAMIC point (missing `fix_static`, `fix_point_sphere`, non-wing
body-frame damping — the network lacks these too). So: **build/extend the NETWORK
component-by-component first; switch the monolith to the components ONCE at the end**
as a single `create_sys!` rewrite from the complete, already-parity-proven set. Same
end state (all shared components, assembly the only difference), lower risk, and every
step advances the actual TTFS goal.

1. **Points + plain segments** — ✅ DONE (network side). Shared `SAM.DynamicPoint`/
   `StaticPoint`/`SpringDamperSegment` with full graph I/O (`segment_io` widened to
   carry `pulley_len`/`tension`); network `network_*` kernels delegate to them; the
   pulley/winch/wing kernels build on the shared `SAM.*` helpers. Validated exact
   (combined `Δaccel=0.0`; 2plate full-accel bit-identical).
2. **Wing-frame, twist, rigid body, live aero** — extend the NETWORK, each as a shared
   `src/components/` component + thin ext wrapper (extin/ND I/O only). Grows the graph
   to cover the full kite (task #12, #13). The `DynamicPoint` superset gaps
   (`fix_static`/`fix_sphere`/body-damp) land here as the models that need them arrive.
3. **Monolith switch (LAST)** — rewrite `create_sys!` to `@named`-assemble the shared
   components + explicit aggregation; delete `point_eqs.jl`/`segment_eqs.jl`. Needs the
   namespaced-sync plumbing (below). One coherent cut, no bridge coexistence.

Validate each network extension vs the frozen monolith oracle to ~1e-10.

## Enabling plumbing: per-instance namespaced sync (OQ2)

Monolith components are `@named`, so a generic leaf `:extra_mass` becomes
`point_i₊extra_mass`. Per-step sync must address the namespaced name. The registry
entry `path` already carries the instance index; map it to the `@named` prefix in
`build_param_sync`/`param_address`. `survivor_index` already maps the `₊`-leaf.

## Step loop (task #16) — status

- **Structural stepping VALIDATED EXACT.** `init!` + `next_step!` run end-to-end on the
  network (the plumbing — `build_prob!`/`init_backend!`/`NetworkStateGetter`/
  `NetworkControlSetter`/`param_sync` — already existed). Combined winch+pulley
  trajectory over 40 real steps (dt=0.02, torque winch): max |Δpos|=1.8e-10,
  |Δvel|=1.3e-9, |Δtether_len|=3.9e-12 vs the monolith. Points/segments/pulleys/winch
  all correct over a trajectory.
- **Aero-refresh gap (the remaining work).** `refresh_aero!` is struct-level and
  backend-agnostic, but it reads `point.va_b`/`wing.va_b`/`wing.R_b_to_w`/`wing.pos_w`
  each step. The monolith's `get_all_state` scatters `sys.va_point_b → point.va_b`
  from the integrator; the network's `NetworkStateGetter` does not, and a KINEMATIC
  wing is not an integrator entity (it's fitted from ref points). **Fix:** a
  network wing-state reconstruction getter — numeric `R_b_to_w` from the wing's
  ref points (`wing.z_ref_points`/`y_ref_points`, `wing_frame_columns` formula),
  `wing.pos_w/vel_w` from the origin point, `wing.va_b` and per-point `va_b =
  R_b_to_w'·(wind(pos_z)·wind_gnd − vel)` — run inside/after `NetworkStateGetter` so
  `update_sys_struct!` refreshes it before `refresh_aero!`. Needs a PARTICLE_DYNAMICS
  VSM trajectory test (adapt one of test/test_continuous_aero.jl etc.) to validate.

## Live aero on the network (ContinuousAero / AeroPressure / twist) — design

- **AeroDirect (frozen per-point force) DONE + validated** (commit c6b558c7). The
  network reconstructs each PARTICLE wing's apparent wind in the state getter
  (`wing_kinematics_from_points!`, shared, built on `wing_frame_columns`); `lin_vsm`
  is now a respected kwarg (was hardcoded false → zero aero); aero points selected by
  `is_wing_node`. Aero RHS bit-exact (`max|Δacc|=2.7e-9`).
- **ContinuousAero/pressure need the LIVE symbolic force** (not frozen
  `point.aero_force_b`). Graph-locality holds because the circulation is frozen: a
  panel force = `f(section geom, va, rho, frozen v_ind, δ)`, no panel↔panel coupling
  at eval. A wing node's force = Σ over its adjacent panels, needing neighboring
  struts' LE/TE positions (via `extin`) + frozen `v_ind`/`cl`/`cd`/`cm` (params). So
  decompose `build_panel_force_eqs` **per wing node** (its adjacent-panel slice) in
  the wing-node kernel — NOT a single aero vertex (feed-forward `extin` rule forbids
  wing nodes reading its outputs). Monolith ref: `aero_component(::ContinuousAero)`
  (continuous.jl:216) — `sec_le/te` interpolated from strut points, `panel_force` →
  `point_force[k]` split 0.75/0.25 LE/TE with billow offsets.
- **Twist surface component** feeds `δ` (flap) into the aero. `apply_flap_delta!` →
  `twist_surface_deltas(ss)` = `flap_delta(ts, R_main, R_flap)` from the two hinge
  bodies' orientations. For a first ContinuousAero cut use `δ=0`, then add the twist
  component. Monolith: `twist_surface_eqs!` (twist DOF) + `twist_surface_delta_eqs!`.
- **AeroPressure** (`aero_component(::AeroPressure)`, pressure.jl:102) shares the
  per-panel→per-point structure; reuse the same node decomposition + `δ` wiring.

## Status / done

- **P0 done:** `src/components/components.jl` has explicit-I/O `point_io`/`segment_io`,
  `DynamicPoint`/`StaticPoint`/`SpringDamperSegment` reading `params.<kind>[idx].field`
  and calling shared `point_acceleration`/`segment_endpoint_loads`; `param_unknowns`
  + `ground_wind_vec` in src. Proven to `mtkcompile` to 6 unknowns (`test_chainB.jl`).
- **Reader elimination done + committed** (2b522708): no bespoke readers on the network
  path; cd_tether via struct read + `:structural` drag-free segment type; topology
  constants set once via `set_const_params!`. Only src-monolith `WindFactorReader`
  remains — dies when the monolith adopts the segment component (step 1).
- Each component gets a FRESH `ParamView`/registry (not one global), else
  `param_unknowns` returns cumulative entries; sync collects per-component registries.
