# Param unification: one `params.X[idx].field` registry for both backends

**Branch:** `networkdynamics-backend`. Working doc, not committed (like the other
`networkdynamics_*.md`). Goal set by the user: **no hand-listed `@parameters` /
`make_param` blocks in the ND components.** The component equations must read
`params.points[idx].extra_mass` directly — the same custom-getproperty view the
monolith uses — and the backend difference becomes a *dispatch*, not a parallel
parameter list.

---

## 1. What we have today (the duplication to kill)

Two parallel, hand-maintained descriptions of the same parameter set:

- **Monolith** (`flat_params.jl`): equations read `params.points[i].extra_mass`.
  `ParamView`/`PathView` getproperty walks the path, reads the struct default, and
  in **one** call both (a) mints the symbol and (b) records a `ParamEntry(param,
  PathReader(path), kind)` into the `ParamRegistry`. Single source of truth.
- **Network** (`ext`): the kernel hand-lists symbols
  (`particle_params()` → `make_param(:extra_mass, 0.0)`, …) **and** a second
  function (`record_particle_params!` / `record_edge_params!`) re-lists the *same*
  fields as `add_param!(builder, VIndex(i, :extra_mass), PathReader((:points, i,
  :extra_mass)))`. The kernel symbol list and the reader list must be kept in
  lockstep by hand; adding `aero_force_b` meant editing both.

This also breeds **name/field aliases** that can silently drift:

| kernel param symbol | struct field path |
|---|---|
| `world_damping` | `world_frame_damping` |
| `body_damp` | `body_frame_damping` |
| `wind_gnd` | *computed* (`ground_wind`) |
| `extra_mass`,`drag_coeff`,`area`,`aero_force_b` | same name ✓ |

The unified design **mandates param-name == struct-field-name** (it comes from the
path), deleting the alias table; the two genuine computed values (`wind_gnd`,
`cd_tether`) keep the existing `param_computed!` escape hatch.

## 2. The one real constraint

ND compiles **one kernel per component *type***, not per instance. So the *symbol*
in the kernel equation must be **generic** (`:extra_mass`, identical for every
vertex); the instance index lives only in the runtime **address**
(`VIndex(i, :extra_mass)`) and the **reader** (`PathReader((:points, i, …))`).

Three layers, only the middle two carry the index:

| layer | monolith | network |
|---|---|---|
| symbol (compile once) | `p_points_5_extra_mass` (index baked in name) | `:extra_mass` (generic) |
| value slot (runtime) | that unique symbol | `VIndex(5, :extra_mass)` |
| reader (per step) | `PathReader((:points,5,:extra_mass))` | **identical** |

The reader is already backend-agnostic; sync already shared (`ParamGroup` +
`sync_params!`). **Only the symbol-naming and the address-resolution differ** —
exactly the two things to put behind dispatch.

## 3. Design

### 3.1 Backend-tagged view

Parameterise the view on the backend so getproperty dispatches:

```julia
struct ParamView{B<:ModelBackend}; reg::ParamRegistry; end
struct PathView{B<:ModelBackend}; reg::ParamRegistry; path::Tuple; end
```

`params = ParamView{typeof(sam.backend)}(reg)`. Component code is **identical**
across backends:

```julia
accel = point_acceleration(s, pos, vel, force,
    params.points[idx].extra_mass + mass_in,
    params.points[idx].drag_coeff, params.points[idx].area,
    params.points[idx].world_frame_damping, wind_gnd, …)
```

### 3.2 Symbol policy — dispatch in `leaf_param!`

- `MonolithBackend`: `make_param(param_name(path), value)` — current behaviour,
  memoised per full path (distinct symbol per instance).
- `NetworkBackend`: `make_param(leaf_name(path), value)` — generic, memoised per
  **(component-kind, leaf_name)** so the kernel reuses one `Num` and ND sees one
  symbol per type.

### 3.3 Address + entry — store the scope, not a baked name

Generalise the entry to carry the *path* (it already has the reader):

```julia
struct ParamEntry; symbol; reader; kind; path::Tuple; end
```

`param_address(::MonolithBackend, entry) = entry.symbol`
`param_address(::NetworkBackend, entry)  = component_index(entry.path)` →
`VIndex(i, leaf)` for `(:points, i, …)` / `(:winches, …)`, `EIndex(j, leaf)` for
`(:segments, j, …)`. `build_param_sync` then sets the setter target via
`param_address(backend, entry)`; the monolith path is unchanged.

### 3.4 Field-set discovery replaces `record_*_params!`

The kernel is built once, so its getproperty calls only register the *representative*
instance. To bind every instance without re-listing fields:

1. **Kernel build** evaluates equations with a representative `idx`, registering the
   set of leaf paths touched under each component kind (`:points`, `:segments`, …).
2. **Replay**: for each real instance `i` of that kind, re-touch the same leaf names
   (`params.points[i].field`) — a mechanical walk over the discovered field set, not
   a hand list. Each call is idempotent per address and yields the
   `VIndex(i, field) => reader` binding.

Net effect: `particle_params()`, `record_particle_params!`, `record_edge_params!`,
`record_wing_node_params!`, `record_winch_params!`, `record_pos_w_params!`, and the
ext `ParamBuilder` all **collapse** into "component uses `params.…` + generic replay
walk." Computed leaves (`wind_gnd`, `cd_tether`, pulley line-mass, const pulley/tether
signs) stay as `param_computed!` readers, addressed the same way.

### 3.5 Worked example — a component after the refactor

`idx` is a representative instance; `params::ParamView{NetworkBackend}`. The same
source compiles on `MonolithBackend` (namespaced/baked symbols) with no edits:

```julia
function network_wing_node_point(s, params, idx; name)
    vars, pos, vel, force, mass_in, _, pulley_len_out = vertex_io()
    ext, zp1, zp2, yp1, yp2, ovel = body_damp_inputs()
    append!(vars, ext)
    p    = params.points[idx]                     # lazy PathView (:points, idx)
    mass = p.extra_mass + mass_in                 # leaf getproperty registers here
    accel = SAM.point_acceleration(s, collect(pos), collect(vel), collect(force),
        mass, p.drag_coeff, p.area,
        collect(p.world_frame_damping), collect(ground_wind_param(p)))
    rot  = wing_frame_rotation(zp1, zp2, yp1, yp2)
    damp = body_frame_damp_accel(vel, collect(p.body_frame_damping), rot, ovel)
    aero = (rot * collect(p.aero_force_b)) ./ mass
    eqs = [ D.(collect(pos)) .~ collect(vel);
            D.(collect(vel)) .~ accel .+ aero .- damp;
            pulley_len_out ~ 0.0 ]
    return System(eqs, t, vars, param_unknowns(params); name)
end
```

`p = params.points[idx]` binds an **immutable partial PathView** — it registers
nothing; each leaf (`p.extra_mass`, `p.aero_force_b`) is the getproperty that mints
the symbol + `PathReader` + scope. So splitting the prefix into a local is safe and
idiomatic (identical path to `params.points[idx].field`). `particle_params()` and
the whole `record_*_params!` family are deleted.

### 3.6 Extensibility — new fields & custom component models come for free

Because the view has no hardcoded field list, a **custom winch model** (the deferred
`winchmodels_ext_idea`) with new struct fields just writes
`params.winches[idx].inertia`, `params.winches[idx].r_drum`, … in its
`winch_component` equations; symbol, default, reader, per-instance address, and sync
are all generated for both backends. Bounded by four conditions:

1. The field is a **leaf** (`Real` → scalar, numeric array → array, function →
   callable) — or a container type added to `param_descend` (custom winch models are
   already `<: AbstractWinchModel`, which is whitelisted).
2. **Computed** (non-field) values keep the `param_computed!` escape hatch + a named
   reader (e.g. `wind_gnd`, `cd_tether`).
3. **Control inputs** must be carved out of sync (like `set_value`) — the model
   author marks them; not auto-derivable.
4. A new component *type* gets its **own ND kernel** (one per type, automatic), so
   generic symbol names never collide across types.

Not solved by this refactor: making a model's *dynamics* symbolic (e.g. WinchModels'
`calc_acceleration` must be MTK-expressible). The view handles the component's
**data**, not its RHS form.

## 4. Migration steps

1. Land the correctness fix first (see §6) so we refactor known-good code.
2. Add the backend type param to `ParamView`/`PathView`; thread `sam.backend` in.
   Monolith path must stay byte-identical (regression-guard: monolith parity + the
   existing suite).
3. Dispatch `leaf_param!` symbol policy; add `path` to `ParamEntry`; add
   `param_address`; route `build_param_sync` / `build_network_param_sync` through it.
4. Rewrite ND kernels to read `params.<kind>[idx].<field>` (drop `particle_params`,
   `make_array_param(:aero_force_b …)`, `make_param(:pulley_mass …)`, …). Rename
   kernel-local field names to match struct fields (`world_frame_damping`,
   `body_frame_damping`).
5. Replace the `record_*_params!` family with the field-set replay walk.
6. Re-validate: winch/pulley/combined exact parity (0.0 / e-11) and 2plate wing
   parity must be unchanged.

## 5. Open questions

- **O-P1 — representative index & extin.** Kernels rebound per instance via `extin`
  (wing frame, winch `tether_len`) already vary by instance without recompiling; the
  field-set replay must not double-count these non-param inputs.
- **O-P2 — edge vs vertex kind from path.** `component_index` must map `(:segments,
  j, …)`→`EIndex(j)` and everything else →`VIndex(i)`; pulley/winch live on their
  vertex. Confirm the graph-edge ordering `j` matches `edgelist` (already handled in
  `record_edge_params!`).
- **O-P3 — `set_value` exclusion.** Still must be excluded from sync (owned by the
  control setter). The replay walk must skip it — keep the current explicit carve-out.
- **O-P4 — D2 endgame.** Once the monolith is itself assembled from `@named`
  components, its params get MTK-namespaced (`point_5₊extra_mass`) from the same
  generic `:extra_mass`; then even the symbol policy converges and this dispatch
  reduces to address resolution only. Out of scope here, but this design is the step
  toward it.

## 6. Prerequisite — pt7 body-damp under real aero ✅ FIXED (2026-08)

Was: 2plate + `AeroDirect` RHS parity showed steering points off by ~8.16. Root
cause: the monolith's `point_damping_accel` damps **any** DYNAMIC point with nonzero
`body_frame_damping` (generic DYNAMIC branch `point_eqs.jl:424`, not only wing nodes);
under AeroNone the steering bfd was zero so the earlier `is_wing_node` gate coincided,
under AeroDirect it's `[10,10,10]` so the network under-damped. Fixed:
`network_wing_node` now gates on `is_wing_node` **or** nonzero `body_frame_damping`.
Whole-model matched-state RHS parity with a nonzero frozen aero force = MAX |Δaccel|
1.1e-10. Refactor now sits on correct code.
