# `ScheduledBackend` — autonomous experiment brief

> **READ THIS SECTION BEFORE ANYTHING ELSE.**
>
> **This is a self-contained experiment, to be executed end-to-end by Claude,
> independently.** It lives in its own worktree (`SymbolicAWEModels-scheduled`,
> branch `scheduled-backend`, forked off `nd-bodies-autonomous` at `68670a7f`).
> Nothing here should be merged, and no other worktree should be touched.
>
> **You are expected to run this to completion without asking for direction.**
> Every decision that could otherwise become a question is pinned in §3. If you hit
> something genuinely unpinned, choose the option most consistent with §3 and §7,
> write down what you chose and why in `DECISIONS.md`, and keep going. Do not stop
> to confirm.
>
> **Standing instruction (2026-08-07, from the user): do not ask questions, and
> commit without asking.** This overrides the repo's usual "ask before committing"
> convention *for this worktree*, and it overrides §6's parenthetical. Commit at
> every green boundary and keep going until §2 is satisfied. Plan files
> (`DECISIONS.md`, `own_backend_plan.md`, this file) still stay uncommitted.
>
> **This experiment is not finished until §2 is fully satisfied.** A backend that
> runs points and segments is not a result. A backend that runs everything the
> monolith runs, with code that reads well, is the result.
>
> This document supersedes `own_backend_plan.md`. That plan's compile-time
> justification was **measured and refuted** (§4) — do not re-derive it, and do not
> re-open the "should we do this at all" question. The decision to run the
> experiment has been taken on other grounds.

---

## 1. What we are building and why

One more `ModelBackend` alongside `MonolithBackend` and `NetworkBackend`: a
hand-rolled, batched, **depth-N** runtime over per-type compiled MTK kernels.

The motivation is **not speed** (§4). It is that NetworkDynamics forces two
structural constraints on us — a fixed depth-2 schedule and one uniform I/O width per
network — and every ugly thing in `ext/SymbolicAWEModelsNetworkDynamicsExt.jl`
descends from them. We own the runtime, so we simply do not have those constraints.

**The experiment succeeds if the resulting code is something a new reader
understands quickly.** If it ends up as complicated as the ND extension, the
experiment has failed even if every test passes.

---

## 2. Definition of done

Three conditions, all required.

### 2.1 Feature parity with the monolith

Every component the monolith supports, working under `ScheduledBackend`. The
monolith's equation generators in `src/generate_system/` are the checklist:

| generator | must work under `ScheduledBackend` |
|---|---|
| `point_eqs.jl` | DYNAMIC, STATIC, BODY_STATIC points |
| `segment_eqs.jl` | spring–damper segments, drag, nonlinear stiffness |
| `pulley_eqs.jl` | pulleys, incl. wing-node pulleys |
| `winch_eqs.jl` | torque and cascaded-length winches |
| `tether_eqs.jl` | tethers, tether init |
| `body_eqs.jl`, `rigid_body_eqs.jl` | rigid bodies, principal frame, ride points |
| `timoshenko_joint_eqs.jl` | Timoshenko beam joints |
| `joint_eqs.jl` | elastic joints |
| `twist_surface_eqs.jl` | twist surfaces, FIXED/KINEMATIC/DYNAMIC twist |
| `wing_eqs.jl`, `aero_eqs.jl` | all aero modes: VSM direct/linearized, ContinuousAero, AeroPressure, AeroPlate, AeroNone |
| `scalar_eqs.jl` | elevation, azimuth, heading, aoa and the rest of the diagnostics |
| `initial_conditions.jl` | `init!` including `reset_vel`, custom ICs |

### 2.2 The test suite passes by backend swap

The acceptance oracle already exists. These are the test files on `origin/main` at
`v0.13.0`; each must pass with only the backend swapped, except where a test is
inherently monolith-only (linearization, control functions — note them explicitly in
`DECISIONS.md` rather than silently skipping):

```
test_aero_modes          test_beam_replay        test_continuous_aero
test_deserialization     test_flap_aero          test_flap_beam
test_getter_allocations  test_heading_calculation test_helpers
test_joint               test_match_aero_sections test_multi_section_group
test_point               test_pressure_aero      test_principal_body_frame
test_principal_frame_invariance                  test_profile_law
test_pulley              test_quaternion_auto_groups
test_quaternion_conversions                      test_rigid_body
test_section_alignment   test_segment            test_segment_nonlinear
test_static_twist        test_tether_init        test_tether_winch
test_timoshenko_joint    test_transform          test_twist_alignment
test_weighted_ref_points test_wing               test_wing_dynamics
test_yaml_bodies         test_yaml_weighted_ref
```

Parity target: match the monolith to solver tolerance at `t = 0` and after stepping.
The ND backend already achieves ~1e-11 on most of these; treat that as the bar.

### 2.3 The code is good

Non-negotiable, and the actual point of the experiment. See §7.

---

## 3. Pinned decisions

These are settled. Do not relitigate them.

### D1 — Vendor NetworkDynamics' codegen; do not depend on NetworkDynamics

Copy the MTK→kernel compiler into `src/scheduled/codegen/`. It is **MIT licensed**
(`Copyright (c) 2019-2025 Frank Hellmann, Michael Lindner and Hans Würfel and
Contributors`) — retain the copyright notice and state clearly in a header comment
which upstream file each part came from and at which commit.

Measured scope: the transitive closure from `generate_io_function` on the
`mtkcompile=true` path is **43 functions, ~884 lines** (an undercount — the scanner
missed functions passed as values, so budget ~1000–1100) out of 3252:

```
371  ext/MTKExt_utils.jl            metadata / alias / defaults repair
279  ext/NetworkDynamicsMTKExt.jl   generate_io_function itself
234  ext/MTKExt_simplification.jl   simplify_with_mtkcompile
```

Why vendor rather than depend: it takes the two invasive fork patches (array-valued
I/O scalarization, nonnumeric/callable params) with it, so **the fork disappears and
the package becomes registerable**. Depending on ND would mean carrying 18k lines and
an unregisterable fork for one internal function.

Cost, accepted knowingly: we own MTK-version churn in that ~371 lines of metadata
repair. Mitigation: keep the vendored files as close to upstream as possible so
`diff` against a fresh ND checkout stays cheap, and record the upstream commit.

**Strip while vendoring** — these exist only to serve ND's depth-2 coreloop:
- the `ff_to_constraint` promotion of feedforward outputs to algebraic states
  (we allow feedthrough, so outputs stay observed)
- `fftype` as a *scheduling* concept; keep it only insofar as it determines `g`'s
  argument list, and prefer making that signature uniform
- `assume_io_coupling` and its `implicit_output` machinery
- the `Network`/graph/`VIndex`/`EIndex` surface — never construct a `Network`

**Keep `cse = true`.** Measured: 167× compile and 3.8–9.8× runtime on our kernels,
agreeing to ≤8e-15. ND's `cse=false` default is tuned for 2–4 state kernels.

### D2 — Parameters are ours, addressed from `SystemStructure`

We have the topology in `SystemStructure`, so we do not need ND's global symbol
namespace. This kills two documented problems outright:

- `nd_param_disambiguation_issue.md` — ND names params by bare field name, so
  `points[src].anchor_b` and `points[dst].anchor_b` collide in one kernel. With an
  `(instance, slot)` address derived from topology, the collision cannot occur.
- callable params needed a whole fork branch because ND's `p` must be flat numeric.

Reuse verbatim: `ParamView`, `ParamRegistry`, `PathReader`, `ParamGroup`,
`sync_params!`. They are already backend-agnostic. Only the *address* changes.

```julia
struct ScheduledParams{C}
    numeric::Vector{Float64}   # layout owned by us, derived from the topology
    callables::C               # per-instance polars / rigidity laws, typed
end
```

**Be precise about what "not flat" means.** The generated kernels come from
`build_function`, so the numeric argument must stay an indexable buffer — we cannot
hand a struct to generated code. What we escape is ND's *naming and addressing*, not
the buffer:

- layout is computed from `SystemStructure` at assembly time, `(instance, slot)` →
  flat index, which is exactly what `ParamRegistry` already does;
- callables live in a separate typed field, not smuggled through the numeric buffer;
- runtime updates keep working the way they do today — mutate the struct field, then
  `sync_params!` writes the buffer before the step. Where a value must be live
  *during* RHS evaluation, the existing registered-symbolic-function pattern (as used
  by the winch) still applies and is unchanged.

### D3 — Components stay equation-generating functions

`src/components/components.jl` already sources both existing backends, and it works.
Keep that style. Everything that is currently a `wide`-flagged variant collapses into
one function, because there is no widening any more (§5).

Real MTK components (`@mtkmodel`, connectors) are **an optional investigation only,
after §2 is met**, and must not block it. The known blocker is that array-valued
variables misbehave in places; if you look at it, timebox it, characterise the bug
precisely, and write it up in `DECISIONS.md`. Do not start there.

### D4 — Depth-N schedule, feedthrough allowed

Sort the output-dependency graph at build time, layer it, batch each layer by kernel.
Our real chain is depth 4: body states → body pose → ride-point position → segment
force → body `f`. Error on a cycle, naming the components involved.

### D5 — Keep the other backends working

`MonolithBackend` is the parity oracle and must not regress. `NetworkBackend` is a
second oracle; leave it alone. Nothing in `src/components/components.jl` may change in
a way that breaks either, except removing the `wide` machinery once nothing uses it.

---

## 4. Already measured — do not redo

Gate 0 ran on 2026-08-06. Its conclusion inverted the original plan's premise.

ND kernel compiles that looked catastrophic (Timoshenko edge 578 s, `build_prob!`
854 s) were **two ND defaults**, not ND's architecture:

1. `chk_component` (`ND/src/doctor.jl:89`, default on) re-runs each fresh kernel
   through an `AccessTracker`, forcing Julia to codegen the whole body specialised for
   a debug wrapper type and then discard it — **86 %** of wall clock (9 500 of 11 007
   profile samples). `mtkcompile` was 3.4 %.
2. `cse=false` emits the hashconsed DAG as a tree — 276 036 characters of source for a
   model with **535 unique DAG nodes**.

With both fixed: Timoshenko edge **0.40 s**, `build_prob!` **20 s**.

**Consequences for this experiment:**

- **There is no compile-time argument for this backend.** Do not claim one, do not
  benchmark for one, do not let it justify a design choice. We are here for the code.
- Kernel compile costs ~0.1–0.4 s each after a one-off ~60 s warm-up of the
  MTK/Symbolics machinery (paid once per session, visible as whichever kernel is built
  first). Budget accordingly; do not panic at the first kernel.
- Symbolics hashconsing makes a second identical build ~13× cheaper. **Only cold
  first calls are meaningful timings.**
- Count symbolic size on the **DAG** (memoised walk). A naive tree walk reports 45 436
  nodes for a 535-node model.
- **Profile before extrapolating.** Three plausible hypotheses were all wrong; one
  `Profile.@profile` settled it in ten minutes.

---

## 5. What must not exist in the result

The experiment is judged partly on absence. None of the following may appear:

- **Widening.** No `wide` kwarg, no `WIDE_*` constant, no zero-fill superset, no
  uniform I/O width. Each kernel declares its own input and output width; wiring is by
  slot. (`wide` appears 85× in the ND ext and 32× in `components.jl` today.)
- **Feedthrough workarounds.** BODY_STATIC ride points are first-class components with
  their own kernel, not absorbed into their body and re-derived inside every incident
  edge. No cluster-vertex partitioning, no union-find over bodies, no
  `build_body_mixed_network` (268 lines).
- **Per-orientation baked kernels.** Order each edge's endpoints when building its
  wiring; there is no `body_at_src` variant.
- **Combined multi-segment edges.** Parallel connections are free — nothing is keyed
  by vertex pair. No `SimpleGraph`, and no graph at all unless it earns its place.
- **The two-edge Hermite trick** and `extin` plumbing.
- **Hand-minted param names** to dodge symbol collisions (D2 removes the cause).

If you find yourself reintroducing any of these, stop and reconsider the design — it
means a constraint has crept back in that we do not actually have.

---

## 6. Phases

Each phase ends with tests green. Commit at each boundary (ask before committing, per
repo convention). Keep `DECISIONS.md` current as you go.

- **P0 — vendor + spike.** Copy the codegen per D1, strip the depth-2 parts, get
  `DynamicPoint` and `SpringDamperSegment` compiling through it. Hand-wire a
  mass-on-spring and match the monolith. Proves array I/O and struct-backed params
  with no ND dependency. *This is the stop/go: if the vendored codegen cannot be made
  to stand alone cleanly, say so and stop.*
- **P1 — runtime core.** Layout, wiring (gather/scatter), schedule, RHS, buffer cache
  keyed on `eltype(u)`, param sync, state getter. Points, segments, pulleys, winches,
  tethers green.
- **P2 — structural layer.** Bodies, joints, ride points, Hermite points, wrench
  segments, twist surfaces. Mostly *deleting* the ND workarounds against existing
  oracles. `test_rigid_body`, `test_timoshenko_joint`, `test_joint` green.
- **P3 — aero.** All modes per §2.1. `test_aero_modes`, `test_continuous_aero`,
  `test_pressure_aero`, `test_flap_aero`, `test_wing*` green.
- **P4 — completeness.** The full §2.2 list. Diagnostics/getters, deserialization,
  allocation tests.
- **P5 — read it back.** Re-read every file you wrote as if you had never seen it.
  Delete what is not needed. This phase is not optional; §2.3 is a completion
  condition, not a nicety.

---

## 7. Code quality bar

The repo conventions apply in full (`CLAUDE.md`), in particular:

- Lines within 92 characters.
- **No duplicate code.** If two kernels share physics, they share a function.
- Docstrings (`"""`) on every type and function, explaining *how it works*, not the
  story behind it. Not verbose — nobody reads a long docstring.
- Inline comments only for a genuinely non-obvious fact, one line maximum. Prefer
  moving the explanation into the docstring. No session-relevant commentary.
- No `_` name prefixes. No short abbreviations or acronyms — `ref_point`, not `wrp`.
- No `const` outside modules in scripts. No `deepcopy`.
- Everything with a docstring goes in the docs, or the docs build fails.

Specific to this experiment:

- The runtime should be roughly four concepts and no more: **kernel**, **instance**,
  **wiring**, **schedule**. If a fifth appears, justify it in `DECISIONS.md`.
- Nothing below the assembler may mention `SystemStructure`, `Point`, `Body` or any
  other domain type. The assembler translates topology into instances and wiring; the
  runtime sees integers and buffers. **This boundary is what would let the runtime
  become its own package later** — enforce it even though we are the only user.
- Vendored code stays in `src/scheduled/codegen/` and is not restyled to repo
  conventions, so it can be diffed against upstream. Everything we write does follow
  them. Do not blur that line.

---

## 8. Reference material

**In this worktree:** `own_backend_plan.md`, the predecessor to this brief. Its §3–§6
(backend seam, architecture sketch, worked example, deletion list) are still the best
detailed design reference — read them. Its §7 compile estimate is refuted by §4 here,
and its "should we do this" framing is settled. Treat it as an appendix, not a plan.

`~/Code/Kite/SymbolicAWEModels-nd` (branch `nd-bodies-autonomous`) has:

- `nd_param_disambiguation_issue.md` — the param collision D2 removes.
- `nd_cse_issue.md`, `nd_depth_n_issue.md` — upstream issue drafts, context only.
- `networkdynamics_status_report.md`, `networkdynamics_aero_component.md` — what the
  ND backend achieved and how the aero decomposition was designed. Reuse the physics
  decomposition; discard the ND-shaped plumbing.
- `ext/SymbolicAWEModelsNetworkDynamicsExt.jl` — 2900+ lines. Read it as a **catalogue
  of what not to do**, and as the source of the enumeration/classification logic
  (`classify_segments`, `record_*_params!`, geometry setters), which is sound and
  worth carrying over.

Upstream codegen source: `~/Code/SciML/NetworkDynamics`, branch `nonnumeric-params`.
Record the exact commit you vendor from.
