# Plan: Flatten struct reads into MTK parameters in the ODE RHS

## Goal

The 2plate_kite ODE RHS spends a large fraction of its time reading mutable
struct fields through `@register_symbolic` getters rather than on physics.
Eliminate that by flattening every getter-read value into MTK parameters and
syncing them once per `next_step!`, without changing results or reintroducing
allocations.

## Baseline (measured)

Script: `test/bench_rhs_profile.jl` (benchmark + thread-filtered CPU profile;
idle worker threads sleeping in libc are filtered out — they otherwise swamp
~97% of samples).

- Model: 2plate_kite, RIGID_DYNAMICS, winch braked.
  15 points, 27 segments, 1 wing, 1 winch, 3 twist_surfaces, state length 67.
- **RHS: 3.79 µs median, 0 allocations.**

Self-time breakdown (14.5k active samples):

| bucket | share | nature |
| --- | --- | --- |
| `datatype.c` (runtime field-offset lookup) | 13.1% | abstract-field tax — removable |
| `julia.h` (layout/type introspection) | 10.5% | abstract-field tax — mostly removable |
| `gf.c` (dynamic dispatch) | 12.9% | partly removable |
| `types.jl` `getproperty` forwarding | 5.8% | removable |
| `accessors.jl` getter bodies | 5.0% | the access path itself |
| `float`/`int`/`math` (arithmetic) | ~20% | floor |
| RGF glue + `getindex` (u/du/param indexing) | ~19% | floor |

Top single self-time leaves: `*` (float.jl, 6.8%),
`macro expansion` (RuntimeGeneratedFunctions.jl:200, 5.1%),
`ijl_field_index` (datatype.c, 4.7%), `get_aero_jac` (accessors.jl:368, 4.4%),
`getproperty` (types.jl:125-126, 3.2% + 2.0%).

Note: this overhead is pure CPU (method lookup + offset computation). The
getters return `isbits` scalars, so dynamic dispatch does **not** box — which is
why the benchmark is correctly 0 allocs and yet the cost exists.
Allocation-free ≠ dispatch-free.

## Upside

Flattening removes the `@register_symbolic` boundary (`gf.c`) *and* the
abstract-field tax (`datatype.c` + `julia.h` + `types.jl` ≈ 29%) together; the
RHS read becomes a direct load from a flat buffer. Remaining cost is then real
arithmetic (~20%) + u/du/param indexing (~19%). Estimate from self-time shares,
not yet measured — the `Segment`-only prototype (below) converts it to a number.

## Result: Segment prototype (DONE) — the win is not in concrete fields

The Segment family (all-scalar) was converted end-to-end (`flat_params.jl`:
registry + `params` view + `setp`-based `sync_params!`). Verified:

- **Correct & live:** perturbing `segments[i].l0` flips exactly one
  `tunable` slot and shifts `du` (proves the param is compiled into the RHS and
  `sync_params!` propagates). `test_segment.jl` 61/61, `test_bench.jl` 15/15,
  RHS still **0 allocations**, "no package allocations in RHS" still passes.
- **RHS time unchanged: 3.79 → 3.79 µs (0% change).**

**Why no win:** segment fields are *concrete scalar reads* through `psys`
(`SystemStructure` is concretely typed) — a field-offset on a concrete `Segment`
is cheap. The 29% "abstract-field tax" in the baseline is **not** concrete field
reads; it is the **abstract `aero::AbstractAeroModel` field inside `VSMWing`** —
every `wing.aero.engine.…` read does a *runtime* offset lookup (`datatype.c`).
Flattening cheap getters removes cheap work.

**Measured (DONE):** two concrete-typing passes, both 0 allocations:

| stage | RHS median | vs baseline |
| --- | --- | --- |
| baseline | 3.79 µs | — |
| `Wing{A}` (concrete `aero`, abstract `engine`) | 3.469 µs | 8.5% |
| + concrete `engine` (`AeroLinearized{E}` etc.) | **1.703 µs** | **55% (2.23×)** |

`Wing{A}` removed the `datatype.c` offset lookup on `wing.aero.engine` (small).
The big win was the **engine pass**: the force reconstruction reads the whole
Jacobian element-wise (`aero_x + aero_jac·Δ`), and with `engine::Union{Nothing,
VSMEngine}` (bare `VSMEngine` is a UnionAll) every coefficient access was a
separate dynamic dispatch. Parameterizing the aero modes by the concrete engine
type (`AeroLinearized{E}`/`AeroDirect{E}`/`ContinuousAero{E}`, engine built and
attached via `attach_engine!` *before* the wing is constructed) collapsed all of
them. (Initial estimate for this pass was ~2-4% — wrong by an order of magnitude;
it was the dominant cost.)

**Pivot:** the lever is the abstract `.aero` field, not flattening per se.
Concretely typing the wing (`VSMWing{A}`, `aero::A`) turns the runtime offset
lookup into a direct load, so the existing `@register_symbolic` aero getters
become cheap *without* flattening aero. Decision (with Bart): **do both** — keep
the flat-param machinery for genuinely-mutable scalar/array struct data, but for
the hot path **concretely type the wing and keep `@register_symbolic` on aero**.
The Segment flattening stays (correct, harmless, exercises the machinery) but is
not where time is spent.

## Interpolations: callable parameters, not `@register_symbolic`

A parameter can carry any object (`@parameters p::MyType`, symtype must match the
value). To *call* it inside equations, declare it **callable** and invoke it
directly — no registration:

```julia
@parameters (interp::typeof(itp))(..) = itp   # (..) = callable parameter
eqs = [f ~ interp(Δ)]                          # symbolic call, no @register_symbolic
```

This is exactly how MTKStandardLibrary's `Interpolation` block works (its older
`get_sampled_data(t, buffer::Parameter)` path *does* use `@register_symbolic` +
a hand-written derivative rule — we copy the *callable-parameter* block instead).
So we implement it ourselves in one line; **no standard-library dependency.**

For elastic joints (`get_joint_force(psys, j, kind, Δ)`):
- linear stiffness → plain numeric param, `k·Δ` (already symbolic, no register)
- interpolation stiffness → callable param, `k(Δ)` (no register, no `psys` read)

Both also remove the `psys` read. **Deserialization caveat (Bart):** a callable
param holds the object, so after loading a cached model it must be re-pointed to
the live `sys_struct` object via `setp` (callable params live in the nonnumeric
buffer, not `tunable`) — same repoint `psys` already gets in `init!`. Δ flows
through as a ForwardDiff `Dual`, so the interpolation must be AD-compatible
(DataInterpolations.jl is).

**DONE.** `joint_eqs.jl` now uses `joint_stiffness_term` (numeric `k·Δ` or
callable `k(Δ)`, branched on the build-time field type) and reads
`damping_trans`/`damping_rot` as numeric flat params. `flat_params.jl` gained
`make_callable_param` + a `CallableParamSync` group bundled with the numeric group
in a `ParamSync`; callables are re-pointed in `init!`/`reinit!`. Anchors
(array-valued) are **left as registered getters** (array fields aren't supported
through the view; joints aren't a hot path).

## Profiled the 1.703 µs RHS (settles what's left)

After the engine pass, profiling the RHS (3M calls, thread-filtered):
- **Aero coefficient reads ≈ 25%** — `get_aero_jac` alone 15.5% (self 6.4%), plus
  `get_aero_x`/`get_aero_y` + the `getindex` chains. Still concrete (fast per
  read) but there are *many* (the Jacobian element-wise) and each pays the
  `@register_symbolic`/`getindex`-through-`psys` boundary. **Flattening these into
  flat array params is a real win** — contradicts the earlier "concrete ⇒ 0%"
  guess (segments were few; aero is many). TODO.
- `calc_wind_factor` ≈ 1.8% and it is **arithmetic** (profile-law math); as a
  callable param the math stays → ~0%. Same for `calc_rho`/polars. Convert for
  *cleanliness* (zero `@register_symbolic`), not speed.

Decision (Bart): **flatten everything, eliminate `@register_symbolic`.** Frozen
`(psys,idx)` reads → flat params (numeric/array); state-dependent calls
(atmosphere, polars, joint interp) → **callable params** (called with live state,
math runs inside, drops the macro); pure functions of state with no data wrap as
stateless callable params or get inlined symbolically.

### Quaternion cleanup (DONE)

`rotation_matrix_to_quaternion` rewritten from `if`/`elseif` to `ifelse` + nested
selection so it evaluates symbolically on `Num` (and still on `Float64`), removing
the 4 `@register_symbolic` component getters (`rotation_matrix_to_quaternion_w/x/y/z`)
that each recomputed the *whole* quaternion. Now one CSE-shared symbolic
expression. Round-trip verified to 2e-15 over 2000 random rotations + edge cases.
It's `getsym`-only (wing orientation output, once per `next_step!`), so this is
**cleanup, not a RHS win** — but −4 registers.

## Validated milestone

| stage | RHS | vs baseline | tests |
| --- | --- | --- | --- |
| baseline | 3.79 µs | — | — |
| `Wing{A}` + concrete `engine` | 1.703 µs | 55% | 1739/1739 |
| + aero coeffs flattened (`aero_y/x/jac`, `drag_frac`) | **1.375 µs** | **64% (2.76×)** | **full suite 1739/1739** |

All 0 allocations. The aero flatten validated the hard machinery: path-view
(`params.wings[w].aero_jac`), array flat params (`make_array_param`), and the
**subsystem name-matching sync** (`survivor_index` matches the registry's bare
param to the namespaced `aero_1₊…` survivor). Shipped also: joints
(numeric+callable params), segment flat params, quaternion symbolic rewrite (−4
registers), parametric-constructor fix.

## Remaining "flatten all" — scope + the one real obstacle

Goal (Bart): zero `@register_symbolic`. Two tiers:
- **Top-level accessors** (points/segments/winches/wings/twist/pulleys/bodies):
  mechanical — path-view + array flat-param support, convert reads, delete
  getters. Perf ~0% (concrete reads), pure cleanliness. No subsystem issue.
- **Subsystem-entangled** (aero coeffs `aero_jac/x/y`; continuous/plate polars):
  these reads live inside `aero_component`'s child `System(eqs,t,vars,[psys])`.
  Flat params declared there get **namespaced** (`aero_1₊…`) at flattening, so a
  registry's bare param object won't match `setp`. Fix: match registry params to
  `parameters(sys)` **by name across the namespace boundary** (the same kind of
  machinery `make_psys_setter` uses, which filters psys by type, not object).
  Aero coeffs are the measured **~25%** win; polars are arithmetic (~0%).
- **Atmosphere** (`calc_rho`/`calc_wind_factor`): top-level callable params,
  arithmetic (~0%), cleanliness.

## Flatten-all phase (data → params, functions → register)

Decision refined (Bart): **a callable param wrapping a function is no improvement
over `@register_symbolic`** — keep functions registered, flatten only *data*.
Reverted the callable-param mechanism (atmosphere `calc_wind_factor`, joint
interpolation stay registered; `rotation_matrix_to_quaternion` stays the inlined
*symbolic* version — that one's a real win, not a callable param).

Converted the remaining **data** getters to flat params (points, segments, wings,
twist surfaces, pulleys, rigid bodies, settings, winch top-level): **RHS 1.375 →
1.197 µs (now 3.17× over baseline), 0 allocations; full suite 1739/1739.**

Still registered, by design: IC getters (defaults/guesses — see rule below),
genuine functions (`calc_wind_factor`, joint interpolation, polars), `psys` anchor
(`get_body_R_b_to_p` for point-less models), winch-component getters (custom
interface), and the aero-mode subsystem data getters (direct overrides,
continuous `v_ind`, plate data) — not yet converted.

Two hard rules learned (each cost a rebuild):
- **Equations → flat params; defaults/guesses → keep registered getters.** A flat
  param referenced *only* in a default/guess (initial conditions: `pos_w`/`vel_w`
  of dynamic points, `com_w_0`/`Q_p_to_w_0`, `twist` of dynamic surfaces,
  `set_value`, `tether_len`, `pulley_len/vel`, `winch_vel`) is **pruned** by
  `mtkcompile`, then the default references a missing param → build error. ICs
  aren't in the RHS anyway, so registered getters there cost nothing.
- **`psys` must stay anchored in an equation.** When every equation-level getter
  in a model is flattened, `psys` becomes orphaned (only in IC defaults) →
  `mtkcompile`: "psys is present but not an unknown". Point models keep it alive
  via `calc_wind_factor(…, psys)`; rigid-body/joint/beam models (no points) need
  `get_body_R_b_to_p` left registered to re-anchor it.
- **Custom winch components can't take a `params` kwarg** (user-defined
  `winch.model`), so the default winch component's getters stay registered (1
  winch, negligible) rather than break the extension interface.

## Status: data flattening COMPLETE (1739/1739, 3.17×)

Every equation-level data getter is now a flat MTK param — top-level *and*
aero-mode subsystems (direct overrides, continuous `v_ind`, plate geometry). RHS
**3.79 → 1.197 µs (3.17×), 0 allocations; full suite 1739/1739** (validated after
each phase). Three commits on `rhs-perf`: `01b8c1c5` (Wing/engine concrete-typing),
`dc1fb82a` (top-level data), `e7e6c769` (aero-mode subsystems).

`@register_symbolic` deliberately retained (not removable as flat params):
- **Functions of state:** `calc_wind_factor`, joint interpolation `get_joint_force`,
  polars `get_plate_cl/cd`, `get_continuous_cl/cd/cm`.
- **Initial-condition reads** (defaults/guesses — pruned if flattened): `get_pos_w`,
  `get_vel_w`, `get_com_w/com_vel/Q_p_to_w` (+ body variants), `get_pulley_len/vel`,
  `get_winch_vel`, `get_set_value`, `get_tether_len`, `get_twist`.
- **`psys` anchor:** `get_body_R_b_to_p` (point-less models).
- **Winch component:** `get_winch_*` (custom `winch.model` interface).

**Dead-code deletion (DONE).** Removed ~65 unused registered getters from
`accessors.jl` (537→148 lines) + `engine_aero_*` barriers + dead plate/continuous
data getters; `test_bench.jl` per-getter allocation tests dropped. Full suite
**1727/1727** (12 fewer than 1739 = the removed getter-allocation tests).
Subtlety: 5 getters my first ASCII/`(psys`-only scan misclassified as dead were
actually live and had to be restored — `get_l0`, `get_twist_ω`, `get_ω_p`,
`get_body_ω_p` (Greek names + a pulley guess the regex missed) and `get_wind_vec`
(used as a bare function reference in `param_computed!`, kept as a plain reader).

**Cache caveat:** deleting the getters means old cached `.bin` models that
reference them fail to deserialize. The load path's try/catch rebuilds, but stale
`v0.12.0` caches in `data/` had to be purged here. Users rebuild on first run
(or `remake=true`).

**ICs could be uniform too (not done).** MTK's `Initial(x)` parameters are
settable via `setp` — `initial.points[i].pos_w → Initial(pos[:,i])` would give a
uniform `initial`-view mirroring `params`, syncing ICs the same way instead of
registered getters in defaults. Feasible follow-up; current ICs use registered
getters.

## Core idea

Every `@register_symbolic get_field(psys, idx)` call site already encodes the
full spec to flatten: a component family (in the name), a field, an index,
scalar-vs-array, and the size. The registered getters in `accessors.jl` *are*
the schema — we reinterpret them as MTK parameters instead of runtime reads.

A getter qualifies iff **all its arguments are `(psys, constant_idx)`** — its
value is constant across a `step!`. Getters that take a symbolic `@variable`
(e.g. `get_continuous_cl(..., alpha)`, `calc_wind_factor(..., z, ...)`) are
genuine nonlinear functions of state and **must stay registered**. The
classifier is mechanical: *does the call site pass any unknown?*

**Param vs state is per-instance, not per-field.** A field like `pos_w` is a
param for STATIC points (`pos ~ params.points[i].pos_w`) but ODE state for
DYNAMIC points (`D(pos) ~ vel`, where the param is simply never referenced). So
there is no global "exclude these fields" list — the equation generator already
branches on the component's role and only reaches for the param where the value
*is* a param. This is exactly the branch the registered getters already live
under (`get_pos_w` is only called on the static branch).

At RHS time a parameter read compiles to a **direct load from the flat
`MTKParameters.tunable` buffer** — none of the `datatype.c` / `gf.c` /
`getproperty` tax. The struct is read only once per step, by a sync function,
not thousands of times by the RHS.

## Why this is exact, not an approximation

Nothing mutates `sys_struct` *inside* the RHS — it is only written by
`update_sys_struct!` / `refresh_aero!`, which run *after* `step!`. So the "live"
reads never actually change during a solve. Freezing them per `next_step!` (and
re-syncing the aero params after each `refresh_aero!`) gives **bitwise-identical**
dynamics. The only obligation is *completeness*: sync everything the RHS reads.
(FBDF's internal sub-steps + Newton iterations within one `step!` see constant
params — which is exactly what params are for.)

## MTK facts this relies on (verified, MTK 11.26.8 / SII 0.3.48)

- `MTKParameters` (in `ModelingToolkitBase`) stores tunable params in one flat
  `Vector{Float64}` (`.tunable`). Scalar and array-valued params
  (`@parameters p[1:n]`, N-D `p[1:a,1:b,1:n]`) are laid out **contiguously**.
- `setp(sys, syms)` builds an **allocation-free in-place setter** once; calling
  it `copyto!`s into the buffer. `parameter_index` exposes offsets for a
  direct-write fast path later.

## Measured: MTK prunes unreferenced params — so register, don't enumerate

`mtkcompile` **deletes parameters that no equation references.** This rules out
"just declare a param for every field / every instance and sync them all," and
dictates the registration design below. Verified with `test/mwe_hole_array.jl`
(MTK 11.26.8); three cases, `n=4` "points", only the two static ones referenced:

| case | setup | result |
| --- | --- | --- |
| **A** — array `pos[1:4]`, only `pos[1],pos[2]` used | hole in the array | **fails at problem build** — `Invalid parameter pos`; the array param shrinks to its referenced subset, a full-width value is rejected |
| **B** — scalar param per point, **only static declared** | no unreferenced params | ✅ `parameters=[p1,p2]`, tunable `[10,20]`, one-call `setp([p1,p2])` works |
| **C** — scalar param for **every** point, dynamic unreferenced | "param for everything" | unreferenced `q3,q4` **silently pruned** from `parameters(sys)`; buffer comes back `[10,20]`; `setp([q1,q2,q3,q4])` throws `MethodError` |

Takeaways:
- A param must exist iff some equation references it. Declaring extras is not
  free bloat — they get pruned and then `setp`/`sync_params!` on them breaks (C).
- A partially-referenced **array** param is worse than per-instance scalars: it
  doesn't just prune, it makes the whole array unusable at build time (A).
- The robust pattern (B) is **size each field's param to exactly the referenced
  instances**. Since referencing is decided per-instance by the generator, the
  clean way to get there is lazy registration (next section), not a static enum.

## The `params` mirror

`params` is a **build-time symbolic mirror** of `sys_struct`, not a runtime
object. It exists only while equations are constructed; its leaves are MTK array
parameters, and indexing one returns a symbolic term baked into the RHS. At
runtime there is no `params` object — it compiled to flat-buffer loads.

| | what it is | exists | leaf type |
| --- | --- | --- | --- |
| `sys_struct.points[i].extra_mass` | mutable data, source of truth | runtime, between steps | `Float64` |
| `params.points[i].extra_mass` | symbolic array param (via view) | build time only | `Symbolic` |

`sync_params!(integ, sys_struct)` is the one bridge: copies struct fields into
the flat buffer. Called at the top of `next_step!`; aero params re-synced after
`refresh_aero!`.

## The spec is lazy registration, not a static enum

The pruning result (above) means we can't pre-declare params from a field list —
we must declare exactly the `(field, instance)` slots the equations actually
reference. The generator already touches each one, at the point where it today
calls a registered getter. So make a registry entry **at the call site**:

```julia
# registry entry = (symbolic param, value thunk)
param!(reg, family, field, idx, value_thunk) -> symbolic_param
```

Wherever the generator today writes `get_extra_mass(psys, i)`, it instead writes
`param!(reg, :points, :extra_mass, i, () -> sys_struct.points[i].extra_mass)`,
which returns the symbolic param *and* records `(symbol, thunk)`. After
generation the registry holds exactly the referenced slots — no holes, nothing
to prune, nothing to enumerate by hand.

The registry is the single source of truth, and it drives everything:

- **param construction** — group entries by `(family, field)`; build each array
  param sized to the recorded instance count (with the instance→row map). Fields
  read for every instance (`l0`, `extra_mass`) get full-width arrays; fields read
  only on one branch (`pos_w` → static points) get subset-sized arrays.
- **the `params` mirror** — the view layer (below) reads from these arrays.
- **`sync_params!`** — for each entry write `thunk()` into the param's buffer
  slot. The thunk is the *only* thing that knows the struct field, so the same
  mechanism handles generic components and custom aero structs alike (see Aero
  modes).

"Less manual work than register_symbolic" still holds: a getter call site goes
from a `@register_symbolic` declaration + a getter body to a single `param!`
call that carries the field read inline as the thunk.

## Access ergonomics

Storage is field-major underneath (one MTK array param per `(family, field)`,
component index last). A thin build-time view layer makes call sites read
**character-for-character identical to `sys_struct`** — same plural collection
names — so migrating the `*_eqs.jl` files is a pure root-variable rename
(`sys_struct` / `psys` → `params`), nothing else.

### Field-major storage

`params._points` is a `NamedTuple` of MTK array params; the component index is
the last index.

```julia
function build_family(prefix, collection, fields)
    n = length(collection); sample = first(collection)
    pairs = map(fields) do f
        dims = size(getfield(sample, f))           # () scalar, (3,) vec, (3,3) matrix
        name = Symbol(prefix, :_, f)               # :point_extra_mass
        f => only(@parameters $name[dims..., 1:n]) # scalar→p[1:n], vec→p[1:3,1:n]
    end
    NamedTuple(pairs)
end
```

Raw storage access (if ever needed directly):
- scalar: `params._points.extra_mass[i]`
- vector: `params._points.body_frame_damping[:, i]`
- matrix: `params._wings.R_b_to_p[:, :, i]`

### Build-time view → exact `sys_struct` mirror

`params.points` returns an indexable family; `params.points[i]` returns a tiny
**build-time** accessor whose `getproperty(:extra_mass)` indexes the underlying
field-major flat param. Pure symbolic sugar — runs only during equation
construction, **zero runtime cost** (it yields the same symbolic term as the raw
storage access above, just nicer to write).

```julia
struct ParamFamilyView; family; idx; end           # build-time only
Base.getindex(fam::ParamFamily, i) = ParamFamilyView(fam, i)
function Base.getproperty(v::ParamFamilyView, f::Symbol)
    p = getfield(getfield(v, :family), f)           # the field-major array param
    nd = ndims(p)
    return nd == 1 ? p[v.idx] : @view p[ntuple(_->(:), nd-1)..., v.idx]
end
```

Call sites are then identical to `sys_struct`, scalars and arrays alike:
- `params.points[i].extra_mass`         ≡ `sys_struct.points[i].extra_mass`
- `params.points[i].body_frame_damping` ≡ `sys_struct.points[i].body_frame_damping`
- `params.wings[i].R_b_to_p`            ≡ `sys_struct.wings[i].R_b_to_p`

The array-index dimension is hidden inside the view, so there is no `[:, i]`
wart at the call site. The only cost is the small `getindex`/`getproperty`
overload (build-time only — it never touches RHS perf).

## Codegen the sync function — two implementations

Both are driven by the same registry of `(param, thunk)` entries and both are
allocation-free. They differ in *where the value goes*: **let SII own the buffer
layout (`setp`)** vs **bake the buffer offsets into the function ourselves
(`@generated`)**.

### Implementation 1 — `setp` setters (SII owns the layout)

Group the registry by `(family, field)`; build one setter per group once at
model-build time:

```julia
setters = (extra_mass = setp(sys, points_extra_mass), ...)  # built once, from registry
staging = (extra_mass = Vector{Float64}(undef, n_referenced), ...)  # preallocated once
function sync_params!(integ, sys_struct)
    for (k, thunk) in enumerate(group_thunks.extra_mass)
        staging.extra_mass[k] = thunk()            # thunk reads sys_struct.points[i].extra_mass
    end
    setters.extra_mass(integ, staging.extra_mass)  # SII copyto!s into .tunable
    # ... one such block per (family, field) group
end
```

- `setp` returns a closure that knows the parameter's slice of `.tunable` and
  does the `copyto!` for us. We never compute an offset.
- **Pros:** robust and future-proof — if MTK changes the buffer layout, scrambles
  ordering, or splits a param across buffers, `setp` still does the right thing.
  Easy to verify against `getp`. This is the documented, supported path.
- **Cons:** one indirection per group (the setter is a runtime callable, not
  inlined), plus we gather into a staging array and then `setp` copies it in —
  effectively two passes over each group's data. Both are cheap, but not nothing.
- The gather loop is generated from the registry — no `~80` hand-written copies.

### Implementation 2 — `@generated` direct buffer write (we own the offsets)

Capture each parameter's offset range in `.tunable` once via `parameter_index`,
then emit straight-line stores into the flat buffer:

```julia
@generated function sync_params!(integ, sys_struct)
    # at compile time, walk the registry and emit one unrolled store per entry:
    #   buf = parameter_values(integ).tunable
    #   buf[off + k] = sys_struct.points[i].extra_mass        # scalar
    #   copyto!(view(buf, rng), sys_struct.points[i].body_frame_damping)  # array
end
```

- One pass, writes directly into `.tunable`, no staging array, no per-field
  closure call — the whole thing inlines to a flat sequence of stores.
- **Pros:** the lowest-overhead form; truly zero dispatch, zero indirection.
- **Cons:** we hard-code offsets, so it is **coupled to MTK's buffer layout**. If
  a future MTK reorders or re-buckets parameters (tunable vs initials vs
  constant), the offsets silently go wrong. Needs the offsets recaptured whenever
  the system is rebuilt, and a guard test comparing against `getp`.

### Which, when

Start with **Implementation 1** — it is correct by construction and the only
thing to code-gen is the trivial gather loop. Move to **Implementation 2** only
if the `Segment` prototype shows `sync_params!` itself on the profile. It almost
certainly won't: sync runs **once per `next_step!`**, the RHS runs **many times
per step**, so even a 2-pass `setp` sync is negligible against the RHS savings.
The win we are chasing is in the RHS (flat loads), which is **identical for both**
— the sync implementation only affects the once-per-step overhead.

## Custom components (aero, and winch)

**One mechanism, not two.** The `params` view auto-registers by recursively
descending the *concrete* `sys_struct` at build time (we have the instance, so
`getfield`/`fieldnames` work even for abstractly-typed fields). Rule at each
`.field`: numeric leaf → register a param; struct → descend and extend the path.
This reaches custom structs too — `params.wings[w].aero.engine.aero_jac`
auto-registers exactly like `params.points[i].extra_mass`. Cadence is inferred
from the path (under `.aero` → refresh-cadence; elsewhere → step-cadence).
Explicit `param!` shrinks to the **escape hatch** for values that aren't a plain
reachable field read (computed/derived, or not on a stable `sys_struct` path).

The author safeguard is unchanged: auto-registration only fires on fields read
*through the view*. Live-state-dependent values (a polar at the live `alpha`) are
built from connector `@variables`, never read through the view, so they're never
frozen — they stay registered symbolic functions.

**Make winch match aero.** Today `winch.model::Function` selects a builder, with
winch data on the core `Winch` struct. Refactor it to the aero shape: a
`winch.model` **struct instance** dispatching `winch_component(model, sys_struct,
winch_idx; name)` (mirroring `wing.aero` / `aero_component`). This (a) gives
custom winches a per-instance data slot the recursive view can descend into —
same auto-param path as aero — and (b) makes "custom component = struct +
dispatch" one consistent extension pattern across the codebase.

**Superseded by the Segment measurement.** The original claim here — that
flattening aero dissolves the `datatype.c`/`gf.c` tax "without separate
concrete-typing surgery" — is **wrong**: the Segment prototype showed concrete
field flattening buys nothing, because the tax is the *abstract `.aero` field*,
not the `psys` read. The fix is concrete-typing the wing (next section), which
makes the existing `@register_symbolic` aero getters cheap. Aero is **not**
flattened.

## Follow-up: concretely type the wing (enforce homogeneity)

Separate from flattening (the RHS already doesn't touch `.aero`), but composes
with it and worth doing: replace `VSMWing.aero::AbstractAeroModel` with a type
parameter `VSMWing{A<:AbstractAeroModel}`, propagated through
`SystemStructure{W}`. Decision: **wings are homogeneous in aero mode**, and the
concrete type *enforces* it — a mixed-mode model won't type-check into one
`NamedCollection{VSMWing{A}}`, turning the pathological case into a compile error
instead of a silent abstract-collection perf cliff.

- **Payoff:** type-stable per-step coupling loops (`refresh_aero!`,
  `update_sys_struct!`, VSM glue) — removes the dispatch/allocs flattening leaves
  in the cold path. Not a RHS win; a per-step-latency + cleanliness win.
- **Main cost — construction ordering:** a type-param field can't be mutated
  across types, so the current "build wing, then `setup_aero!` swaps in the
  engine" flow must flip to **resolve the aero mode first, then build the wing
  with its final `A`**. Everything else is propagating the type param up.
- **Assumption baked in (confirm):** no single model needs mixed aero (e.g. main
  wing REFINE + tail in a cheaper mode). Taken as out of scope for these models.

## Adding a param (author's view)

Which package owns each piece — the recorder is **SymbolicAWEModels**, not MTK:

| symbol | package | role |
| --- | --- | --- |
| `param!` | **SymbolicAWEModels** (new) | records `(param, thunk)` in the registry; returns the param |
| `@parameters` / `Symbolics.variable` | MTK / Symbolics | constructs the symbolic param |
| `@register_symbolic` | MTK | the mechanism being **removed** |
| `setp` / `getp` | SymbolicIndexingInterface | write/read the flat buffer in `sync_params!` |

Name the recorder `param!` (not `register`) — `register` collides with MTK's
`@register_symbolic`, the very thing we're replacing.

### Default: just read the field through the view (auto)

Add the field to its struct (core *or* custom), then read it through `params`;
the recursive view registers param + thunk + sync, with cadence inferred from the
path:

```julia
ramp = params.winches[w].ramp_time            # field on the core Winch struct
jac  = params.wings[w].aero.engine.aero_jac   # field inside the custom aero struct
```

No `@register_symbolic`, no getter body, no sync code, no `@parameters` — strictly
less than today (which also needed a getter in `accessors.jl`).

### Escape hatch: explicit `param!` for non-field-read values

Only when the value isn't a plain reachable field — a computed/derived value, or
one not on a stable `sys_struct` path:

```julia
p = param!(reg, :stem, () -> compute(...); dims=(nx, ny), cadence=:refresh)
```

`param!` builds the symbol via `Symbolics.variable` (functional `@parameters`) and
records the thunk. Same rule as everywhere: a value that depends on a live
connector input (a polar at the live `alpha`) is *not* a param — keep it a
registered symbolic function.

## Risks

- **Compile-time / scalarization blowup.** The referenced params scalarize into
  many scalar params in the index cache. RHS runtime is ideal (flat loads) but
  `mtkcompile` time and codegen size may grow. **Prototype one family
  (`Segment`, all-scalar, full-width) end-to-end and measure RHS µs *and* build
  time before converting the rest.**
- **Pruning is a sharp edge, not just a caveat (measured above).** Register a
  param only where an equation references it; never pre-declare from a field
  list. Prefer per-instance / subset-sized arrays over full-width arrays with
  holes (case A fails at build; case C silently prunes and breaks `setp`).
- **Don't sweep up the must-stay-registered getters** (continuous-aero polars,
  wind factor) — they depend on `u`. Applies inside aero modes too.
- **Structural flags** (`fix_static`, `speed_controlled`) are build-time-known —
  bake into equation generation as Julia `if`, don't make them params.

## Rollout

1. Build the registry + `param!` + `setp`-based `sync_params!` (Implementation 1)
   for `Segment` only (all-scalar, full-width — no per-instance subset, no aero).
   Wire through `bench_rhs_profile.jl`; compare RHS µs + build time vs baseline.
2. If the win tracks the estimate, roll out family by family, deleting the
   corresponding registered getters from `accessors.jl` as each is converted.
   `pos_w`/`vel_w` introduce the per-instance subset-sizing path.
3. Convert aero last: replace `get_*` reads inside `aero_component` with
   `params.wings[w].aero…` view reads (auto-register, refresh cadence). Keep the
   live-`alpha` polar getters registered. Refactor `winch.model` to the
   struct+`winch_component` dispatch shape so it rides the same auto path.
4. Only if `sync_params!` shows on the profile, switch it to the `@generated`
   direct-buffer form (Implementation 2).

Background: `test/mwe_hole_array.jl` is the standalone reproducer for the pruning
behaviour — keep or delete once the registry design is locked in.

## Verification

- `test/bench_rhs_profile.jl` — RHS median µs before/after; profile shifts.
- `test/test_bench.jl` — must stay **0 allocations**, and the "no package
  allocations in RHS" test must still pass.
- `test/test_aero_modes.jl` and the rigid/particle suites — results unchanged
  (the params hold the same values; only the access path changes).
