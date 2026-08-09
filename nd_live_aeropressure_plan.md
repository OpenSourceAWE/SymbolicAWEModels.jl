# Wire live AeroPressure onto beam-ride wing nodes (ND backend)

Branch `networkdynamics-backend`. Goal: the full SK_coupled.jl pressure kite runs WITH
aero on `NetworkBackend`. The structural DAE already builds/compiles/param-syncs/getter-
constructs (3 blockers fixed 2026-08-07: parallel-spring combiner, `is_parameter` prune
filter, weighted-ref in the getter). The remaining gap is the **live aero force term**:
the RHS currently carries no AeroPressure force.

## Why the current path doesn't cover it

- AeroPressure force is **live symbolic** (not a frozen per-point read): `aero_component`
  (`src/aero_modes/pressure.jl:109-179`) builds each wing point's force from live apparent
  wind `va` + frozen pattern params (`v_ind`/`traction`/`traction_net`/`cl`/`cd`/`cm`),
  via shared `build_panel_force_eqs` (`common.jl:280-348`). `provides_aero_override` is
  false → ND already routes it through the LIVE path, not the frozen `AeroDirect` path.
- The existing ND live path `LiveAeroWingNodePoint` (`components.jl:1254-1292`) only serves
  **DYNAMIC** wing nodes (it integrates pos/vel). The pressure kite's wing nodes are
  **BODY_STATIC ride points** on beam bodies (pos/vel algebraic from the body pose), built
  by `network_rigid_ride_point`/`network_hermite_ride_point`. They currently receive NO
  aero.
- Monolith reference for the physics: aero force is folded into `point_force` BEFORE the
  ride split (`point_eqs.jl:283-307`); `aero_force_w = R_b_to_w * aero_force_point_b`, then
  `body_ride_eqs`/`beam_hermite_ride_eqs` transport force+moment to the body/beam by axial
  fraction. The ND analog: add `aero_force_w` to the ride's `force_out`; the existing
  `network_mount_edge` already delivers `mount_frac·force_out` + `(pos−com)×that` with the
  correct beam fraction — **no mount-edge change needed**.

## Injection point

`ride_force_out_eq` (ext `:408-415`): `force_out ~ force_in + drag [+ gravity]`. Add an
optional `extra_force` (world-frame `aero_force_w`) into `load`. That is the whole coupling
on the delivery side.

## Newly-discovered sub-requirement: weighted-ref frame IN the kernel

The live-aero kernel fits the wing frame from **single** ref points
(`wing_node_inputs()` → 4 scalar ref points, `wing_frame_rotation`), and `live_aero_extin`
(ext `:1314-1327`) calls `ref_single_id(wing.origin)`. The pressure wing's `origin` and
`z_ref[1]` are the **confluence centroid** — multi-point `WeightedRefPoints`. So the
in-kernel frame + the extin must support weighted refs:
- Kernel: declare per-ref extin reading EVERY id of each `WeightedRefPoints`, compute the
  weighted centroid in-equation, then feed `wing_frame_rotation`. (Getter already does this
  via `get_ref_position_from_points`; the kernel needs the symbolic analog.)
- `live_aero_extin`: bind each ref id (not just `ref_single_id`) and the origin ids.

Single-point wings (weight 1, one id) must stay byte-identical (existing ContinuousAero
live test).

## Implementation steps

1. **`components.jl` — factor `live_aero_force_terms(s, params, aero_params, wing, slot,
   num_points)`** returning `(ext_vars, subsys, aero_slot, conn_eqs, rot, ovel,
   aero_force_w, wing_force, wing_moment, agg_force, agg_moment)`. Rewrite
   `LiveAeroWingNodePoint` to consume it (behavior unchanged; `aero = aero_force_w/mass`).
   This de-dups the aero wiring so the ride kernels reuse it.

2. **Weighted-ref frame** — generalize `wing_node_inputs` (or add
   `weighted_wing_frame_inputs(z_ref_points, y_ref_points)`) to declare extin for every ref
   id + weights, and compute weighted centroids feeding `wing_frame_rotation`. Update
   `live_aero_node_inputs` origin to weighted (multi-id). Single-id path identical.

3. **Ride aero kernels (ext)** — give `network_rigid_ride_point`/`network_hermite_ride_point`
   an optional `aero=(wing, slot, num_points, aero_params)`; when present call
   `SAM.live_aero_force_terms`, append `ext_vars` + `wing_force`/`wing_moment`, nest
   `subsys`, add `aero_force_w` via `ride_force_out_eq(...; extra_force=aero_force_w)`, and
   append `aero_slot` to the param unknowns. `ride_force_out_eq` gains `extra_force=nothing`.

4. **`build_body_mixed_network` wiring** — for each PARTICLE wing with `needs_live_aero`
   whose wing nodes are rides:
   - `point_idxs = wing_node_points(ss, wing)` (slot order). `npts = length`.
   - Per wing build a rigid aero base and/or hermite aero base (one `network_view` aero
     registry `areg` per wing) with `ff_to_constraint=true`; rebind each aero-ride point via
     `VertexModel(base; extin = [ride_pose_extin(point); weighted_live_aero_extin(ss, wing,
     point_idxs)])`.
   - `record_aero_params!(builder, callables, ss, ride_vertex, wing, areg)` once per wing
     (params are per-wing frozen buffers; same on every node). Slots via
     `set_live_aero_slots!`-style write of `aero_slot` per aero-ride vertex.
   - Route these points to the aero base in `vmodels[i]` instead of the plain ride base.

5. **Getter** — `network_wing_geoms` `live_vertex` must be an aero-RIDE vertex (already
   `aero_pts[1]`); it now exposes `wing_aero_force_b_1..3`, so `wing_aero_reader` works
   unchanged. Confirm `aero_pts[1]` is one of the aero rides.

6. **Param prune** — the `is_parameter` filter (already added) covers any aero param
   `mtkcompile` prunes.

## Verification

2plate ContinuousAero live particle test must stay green (single-ref path unchanged);
then SK_coupled.jl build+init+step on ND, t=0 aero force parity vs monolith (~1e-8), then
a few steps. Flap `delta` coupling stays out (both backends' ND live path is flap-free;
`components.jl:1251`).

## Watch-outs

- Whole-wing aero evaluated per ride kernel, slot-selected → O(N²); accepted for the
  DYNAMIC path, same here.
- Each aero ride reads ALL wing points' pos/vel via extin (`wpos_k_c`/`wvel_k_c`) — large
  extin count; mirrors DYNAMIC live nodes.
- Don't double-count: BODY_STATIC aero rides must NOT also be built as DYNAMIC live-aero
  nodes (`build_live_aero_vmodels` only wires DYNAMIC; the body-mixed path builds rides).
- Mass/gravity unchanged: rigid aero ride still adds no gravity (body carries it); Hermite
  aero ride keeps its `extra_mass+mass_in` gravity.
