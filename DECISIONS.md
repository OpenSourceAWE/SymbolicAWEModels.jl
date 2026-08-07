# Decisions log — `ScheduledBackend` experiment

Running record for the experiment described in `scheduled_backend_experiment.md`.

Append an entry whenever you resolve something the brief did not pin, deviate from it,
skip a test, or hit a dead end. One short entry each: what you chose, why, and what it
costs. This is the substitute for asking, so it must be honest about the bad choices
too, not just the tidy ones.

Format:

```
## <date> — <short title>
**Choice:** ...
**Why:** ...
**Cost / risk:** ...
```

---

## 2026-08-07 — Experiment opened

**Choice:** Worktree `SymbolicAWEModels-scheduled`, branch `scheduled-backend`, forked
off `nd-bodies-autonomous` at `68670a7f`.

**Why:** That branch carries the two assets the experiment reuses unchanged — the
shared components in `src/components/` and the backend-agnostic parameter registry
(`ParamView` / `ParamRegistry` / `PathReader` / `ParamGroup` / `sync_params!`).

**Cost / risk:** It is ahead of `origin/main`, so the parity oracle is the monolith on
this branch rather than at the `v0.13.0` tag. If the two have diverged in behaviour,
prefer the monolith on this branch and note any test that differs.

## 2026-08-07 — P0 stop/go: the vendored codegen stands alone

**Choice:** Go. `src/scheduled/codegen/mtk_codegen.jl` is ~560 lines vendored from
NetworkDynamics `68756bce`, with no NetworkDynamics dependency. `DynamicPoint`,
`StaticPoint` and `SpringDamperSegment` compile through it (0.78 s / 0.07 s cold,
including the one-off MTK warm-up), and the hand-wired mass-on-spring matches the
monolith **exactly** (`du = [-9.81, 0, 0, 0, 0, 0]`, bitwise).

**Why:** The gate the brief set. Array-valued I/O needs no fork patch — we call
`build_function` ourselves — and struct-backed parameters work through the existing
`ParamRegistry` with the container index swapped per instance.

**Cost / risk:** We now own ~370 lines of MTK metadata repair. Recorded upstream
commit in the file header so `diff` against a fresh checkout stays cheap.

## 2026-08-07 — `fix_metadata!` reports only under `verbose`

**Choice:** The vendored `fix_metadata!` warns about metadata it could not restore
only when `verbose = true`.

**Why:** Every scalarized array element (`pos(t)[1]`, `p_points_2_area[3]`) fails to
resolve by name, so upstream warns on *every* kernel we build — 20+ multi-line
warnings per model. The residue is harmless: `isequal` on symbolics ignores metadata,
so `build_function` still matches the I/O lists. This is the same layer that warns
loudly on the ND backend today.

**Cost / risk:** If metadata loss ever does break something (`isparameter` returning
`false` for a stripped parameter would make `get_alias` mistake it for an alias), it
now fails silently. Mitigated by parity tests, and `verbose = true` still shows it.

## 2026-08-07 — a component whose outputs split across the schedule becomes two

**Choice:** The schedule orders whole instances, using a per-input flag for whether
the kernel's outputs read that input. A component that both *reads* something
upstream and *returns* something downstream through the same instance — a
`BODY_STATIC` ride point, which places itself from its body's pose and forwards the
segment force back to that body — is modelled as two components, not one.

**Why:** Slot-granular scheduling cannot help: both output groups live on the same
instance, so the cycle survives. Splitting keeps the runtime at four concepts and the
decomposition is honest — a ride point's kinematics and its wrench transport are two
different maps.

**Cost / risk:** More instances at `SK100` scale (one extra per ride point). No extra
kernels, so compile time is unaffected.

## 2026-08-07 — the scheduled components are their own file, not `components.jl`

**Choice:** The scheduled backend's component `System`s live in
`src/scheduled/components.jl` (`Particle`, `Anchor`, `PulleyParticle`,
`WinchAnchor`, `SpringSegment`, `PulleySegment`, `TetherSegment`). They are the
same *style* as D3 asks for — plain equation-generating functions over the shared
force laws in `src/components/` — but they declare their own I/O rather than the
`NetworkBackend`'s padded, orientation-baked one.

**Why:** §5 forbids padding and per-orientation kernels, and the ND components have
both: every point carries an unused `tension_in`/`pulley_len_out`, every segment an
unused `src_pulley_len`/`dst_pulley_len`/`src_tension`/`dst_tension`, and the pulley
segment bakes `pulley_at_src`. D5 says the `NetworkBackend` must keep working, so
the padded versions stay for it. All the *physics* is shared: `segment_load_terms`,
`point_acceleration`, `dynamic_point_dynamics`, `pulley_split_eqs`,
`winch_component`.

**Cost / risk:** Two sets of component declarations in the tree until the
`NetworkBackend` is retired, at which point the padded ones and the `wide`
machinery delete together. No physics is duplicated.

## 2026-08-07 — endpoint orientation is a wiring question, not a kernel one

**Choice:** A pulley segment declares one `rest_len` input and one `tension`
output; a tether segment declares one `rest_len` and emits `src_tension` *and*
`dst_tension`. Which endpoint they belong to is decided when the segment is wired.

**Why:** This removes `pulley_at_src`, `tension_sign_src` and `tension_sign_dst`
outright — the three parameters that exist in the ND extension only because a
`SimpleGraph` orients each edge min→max. The only pulley parameter left is
`pulley_side` (`±1`), which is physics (which half of the rope), not orientation.

**Cost / risk:** None found. The tether segment carries three extra output slots
for the endpoint it is not attached to.

## 2026-08-07 — `segment_endpoint_loads` became a view of `segment_load_terms`

**Choice:** The shared segment force law now returns a named tuple with every term
(`len`, `unit_vec`, `spring`, `spring_vec`, `half_mass`, `half_drag`,
`force_on_src`, `force_on_dst`); `segment_endpoint_loads` is the four-value view the
`NetworkBackend` already destructures.

**Why:** The scheduled segment needs `spring_vec` (the winch force) and `half_drag`
(the point's `total_drag`) separately, and re-deriving either would duplicate the
force law. The wrapper keeps the ND path byte-identical.

**Cost / risk:** One extra indirection on the ND path. Also added: a `nonlinear`
branch for a callable `unit_stiffness`, which the shared law did not support and
`test_segment_nonlinear` needs.

## 2026-08-07 — `default_backend!` for the backend swap

**Choice:** Added `default_backend()`/`default_backend!(backend)`; a
`SymbolicAWEModel` built without an explicit `backend` uses it.

**Why:** §2.2 asks for the suite to pass "with only the backend swapped". Without
this, that means editing every `SymbolicAWEModel(...)` call site in 35 test files,
which is not a backend swap — it is a patch. `test/test_scheduled_backend.jl` sets
it, runs the monolith's own files unmodified, and restores it.

**Cost / risk:** Global mutable state. It is a test/driver convenience; the
constructor keyword still wins where it is given.

## 2026-08-07 — P1 status: points and segments are exact

**Choice / result:** With points and segments only, the scheduled backend matches
the monolith **bitwise**: `max|Δdu| = 0.0` at five random states, and after 20 steps
`Δpos = Δvel = 0.0` with `len`, `spring_force`, `l0`, `drag_force` and `total_mass`
all identical. Pulleys, winches and tethers are implemented but not yet exercised.

**Cost / risk:** The exactness is expected, not lucky — both backends evaluate the
same shared force laws over the same parameter values. It does mean any divergence
later is a real assembly bug, not accumulated float noise, which is a good oracle.

## 2026-08-07 — hours lost to a fixture, not the backend

**Note, kept because it is the kind of thing this log is for:** `next_step!` failed
with `MaxIters` at `t = 0` for a long time. The cause was `sample_freq` missing from
the spike's settings YAML, so `dt = 1/0 = Inf`. Three plausible backend-side
hypotheses (feedthrough staleness, the `FunctionWrappersWrapper`, the mass matrix)
were all wrong. What settled it in one run was printing `integrator.dt` — the same
lesson §4 of the brief already records: measure the state, do not extrapolate from
the symptom.

## 2026-08-07 — `NetworkBackend` removed from this worktree (user instruction)

**Choice:** Deleted `ext/SymbolicAWEModelsNetworkDynamicsExt.jl` (2964 lines), the
`NetworkBackend` type, the NetworkDynamics/Graphs weak dependencies and the
extension entry. `src/components/components.jl` is now the one component source:
the shared physics builders the monolith's equation generators call, plus the
component `System`s the `ScheduledBackend` compiles — 777 lines where the ND-shaped
version was 1349 plus a 2964-line extension.

**Why:** The user asked for it directly, and it is what `own_backend_plan.md` §10.4
already wanted. D5 of the brief ("`NetworkBackend` is a second oracle; leave it
alone") is superseded by that instruction. The monolith remains the parity oracle,
and it is the stronger one — it is what the tests were written against.

**What went with it:** every `wide` kwarg and `WIDE_*` constant, `segment_io_wide`,
`wide_pose_appendix`, `finish_vertex`, the padded `point_io`/`segment_io`, the
per-orientation `pulley_at_src`/`tension_sign_*` parameters, `build_body_mixed_network`
(268 lines), `MixedEdgeInfo`, the combined multi-segment edges, the two-edge Hermite
trick and the `extin` plumbing. §5 of the brief is now satisfied by absence, not by
discipline.

**Also simplified:** `ParamView`/`PathView` lost their `ModelBackend` type parameter
along with `param_symbol_name`/`param_cache_key`/`leaf_symbol`. There is one naming
policy — the full path — which is exactly what makes
`nd_param_disambiguation_issue.md` impossible rather than merely avoided.

**Cost / risk:** One fewer oracle. The rigid-body, Timoshenko, twist-surface and
live-aero *components* the ND extension carried are gone with it; when the scheduled
backend reaches P2/P3 they will be rebuilt against the shared expression builders
(`rigid_body_pose_expressions`, `timoshenko_element_wrench`, `elastic_joint_wrench`,
`beam_hermite_ride_expressions`, `aero_component`), all of which survive because the
monolith uses them too. The ND versions remain readable in git history and in
`~/Code/Kite/SymbolicAWEModels-nd`.

## 2026-08-07 — P1 result: the P1 test files pass by backend swap

**Result:** With `default_backend!(ScheduledBackend())` and nothing else changed,
the monolith's own test files run green:

| file | result |
|---|---|
| `test_point.jl` | 52/52 |
| `test_segment.jl` | 61/61 |
| `test_pulley.jl` | 46/46 |
| `test_tether_init.jl` | 41/41 |

Every one of those calls `test_init!`, which asserts the RHS allocates **zero**
bytes — so the batched evaluation is allocation-free through the real integrator
path, not just in a micro-benchmark. Assembly is 0.02–2.1 s and the integrator build
reports 0.0 s, because there is no per-model codegen left to do.

**Two bugs the run found, both fixed:**

1. `compile_kernel` collected the callable defaults into a `Vector{<:Interpolation}`
   where `ComponentKernel` wants `Vector{Any}` — every nonlinear-stiffness segment
   failed to construct. (`test_segment_nonlinear.jl`)
2. `generate_io_function` read the parameter defaults *after* flattening, and
   flattening drops a nested subsystem's defaults, so the winch motor's drum
   parameters came back `nothing`. Reading them off the hierarchical system first
   fixes it. (`test_tether_winch.jl`)

**Cost / risk:** Zero allocations depend on `ODEProblem{true, FullSpecialize}` — the
default `AutoSpecialize` wraps the right-hand side in SciMLBase's function wrappers,
which allocated ~200 bytes per call. Full specialization is the right choice here
anyway (one concrete RHS type), but it does compile more of the solver stack.

## 2026-08-07 — callables are stored per kernel, not per instance in one flat vector

**Choice:** `ScheduledParams.callables` is a tuple with one entry per kernel, each a
vector of that kernel's instances' callable parameters. `run_batches!` walks it in
step with the kernel tuple, so each generated kernel receives its own concretely
typed store; `ComponentInstance` carries a `position` (its index among its kernel's
instances) instead of a slice of a global callable buffer.

**Why:** The first version was one flat `Vector{Any}` sliced per instance. It works,
but the element type is `Any`, so every call to a nonlinear stiffness law inside a
generated kernel is a dynamic dispatch — 640 bytes per RHS call with one such
segment, 1920 with three. `test_segment_nonlinear` catches exactly that through
`test_init!`'s zero-allocation assertion. Per-kernel vectors keep the element type
concrete because all instances of one kernel read the same struct field.

**Cost / risk:** The callable store is a tuple, so `sync_params!` holds direct
references to the per-instance vectors rather than indexing it. If two instances of
one kernel ever carry callables of *different* concrete types the vector widens to
their join and the dispatch goes dynamic again — correct, just slower.

## 2026-08-07 — the runtime/domain boundary, enforced by file

**Choice:** `ScheduledParams`, `ScheduledParamSync` and `sync_params!` moved out of
`runtime.jl` into `params.jl`. `runtime.jl` and `kernel.jl` now contain no mention of
`SystemStructure` or any other domain type; the runtime only ever reads `p.numeric`
and `p.callables` off whatever parameter object it is handed.

**Why:** §7 asks for exactly this boundary — "nothing below the assembler may mention
`SystemStructure`, `Point`, `Body`" — and the parameter sync was the one leak. The
layering is now, bottom to top: `codegen/` (MTK → callables), `kernel.jl` +
`runtime.jl` (integers and buffers), `params.jl` (the first layer that reads the
struct), `components.jl` (the physics), `assembly.jl` (topology → instances and
wiring), `state.jl` + `backend.jl` (the package seam).

**Cost / risk:** None. It is a file move plus a header comment.

## 2026-08-07 — where the experiment stands, and what it does not yet do

**Done.** P0 (vendored codegen, stop/go passed) and P1 (runtime core). The
`ScheduledBackend` runs points, segments, pulleys, winches and tethers, matching the
monolith bitwise, and the monolith's own test files for those pass with only
`default_backend!` swapped. §5's prohibitions hold by absence: `grep` finds no
`wide`, no `SimpleGraph`, no `extin`, no `*_at_src`, no combined edges anywhere in
`src/`.

**Not done.** P2 (bodies, joints, ride points, Hermite points, wrench segments,
twist surfaces), P3 (every aero mode) and P4 (the remaining ~25 test files,
deserialization, the getter-allocation test). §2.1's checklist is therefore
satisfied for `point_eqs.jl`, `segment_eqs.jl`, `pulley_eqs.jl`, `winch_eqs.jl`,
`tether_eqs.jl` and the `initial_conditions.jl` path the scheduled backend uses, and
not for the rest. The design decisions those phases depend on are recorded above —
in particular the two-component split for ride points, which is the one place the
depth-N schedule needed a deliberate choice rather than falling out.

**Deliberately not attempted.** D3's optional `@mtkmodel`/connector investigation.
It is explicitly gated behind §2 being met.

## 2026-08-07 — P1 complete: 267/267 by backend swap

**Result:** With `default_backend!(ScheduledBackend())` and nothing else changed,
every P1 test file is green:

| file | result |
|---|---|
| `test_point.jl` | 52/52 |
| `test_segment.jl` | 61/61 |
| `test_segment_nonlinear.jl` | 11/11 |
| `test_pulley.jl` | 46/46 |
| `test_tether_init.jl` | 41/41 |
| `test_tether_winch.jl` | 56/56 |

267 assertions, and every file calls `test_init!`, so the RHS is allocation-free on
each of them. That closes §2.1 for `point_eqs.jl`, `segment_eqs.jl`, `pulley_eqs.jl`,
`winch_eqs.jl`, `tether_eqs.jl` and the `initial_conditions.jl` path these use.

**On the schedule depth actually exercised:** these models resolve to **two** layers
(points and winches emit state-derived outputs, segments read them). The runtime is
depth-N by construction — `build_schedule` is a plain Kahn layering with no bound —
but the depth-4 chain the brief names (body states → pose → ride position → segment
force → body `f`) only appears once P2's bodies exist, and is not yet demonstrated.
Claiming otherwise would be claiming an untested property.

## 2026-08-07 — the monolith did not regress

**Result:** The same six files on `MonolithBackend` also run 267/267 after the ND
removal, the `components.jl` pruning and the `ParamView` simplification. The oracle
is intact.

## 2026-08-07 — P2 entry point, worked out but not built

Stopping at the P1 boundary with everything green rather than starting the body
layer and leaving the package unloadable. The design is settled, so record it:

**`test_rigid_body.jl` needs only a `RigidBody` component** — its bodies carry no
points and no segments (free fall, spin-up under `ext_moment_b`, COM offset), and
its last testset is `reinit!`-only with no ODE. So the first P2 step is one
component, not the whole structural layer.

`RigidBody(s, params, idx; name, frozen)`, mirroring `body_eqs!`:
- states `com_w`(3), `com_vel`(3), `Q`(4), `omega_p`(3);
- inputs `force_in`(3), `moment_in`(3);
- outputs the pose a ride point reads — `pos`, `vel`, `R_b_to_w`(9), `com`,
  `com_vel`, `omega_w`;
- observables the getter needs — `acc_w`, `Q_b_to_w`, `omega_b`.
- Call `rigid_body_pose_expressions` **once**, with the `ω_kinematic`/`d_ω_p`/
  `d_com_w`/`d_com_vel` overrides written against *torn* `alpha_p`/`com_acc`
  variables that are then bound to `ex.α_p`/`ex.com_acc`. There is no cycle:
  `α_p` and `com_acc` do not depend on the overrides. This is exactly what
  `body_eqs!` does through its `body_α_p`/`body_com_acc` arrays.
- `STATIC` is a build-time branch (its own kernel, since `frozen` is topology);
  `fix_sphere` stays an `ifelse` on the parameter, as in the monolith.

**Then the ride points**, which is where the depth-4 schedule finally gets
exercised. Per the entry above, two components:
- `RidePoint`: body pose in → `pos`, `vel`, `arm` (`pos − com`) out.
- `RideWrench`: `arm` + the segments' aggregated `force_in` → `moment_out`
  (`arm × force_in`) into the body's `moment_in`.
The segment forces themselves wire *straight* into the body's `force_in`; only the
moment needs a component, because only the moment needs the arm. Layers then come
out as body pose (1) → ride position (2) → segment force (3) → ride wrench (4) →
body `f`, which is the chain §1 of the plan says NetworkDynamics cannot express.

## 2026-08-07 — the depth-4 schedule, demonstrated

**Result:** A body hanging from a ground anchor by one spring through a
`BODY_STATIC` point compiles to five kernels and schedules to **four layers**:

```
kernels      = [:rigid_body, :anchor, :ride_point, :ride_wrench, :spring_segment]
layer widths = [(1,1,0,0,0), (0,0,1,0,0), (0,0,0,0,1), (0,0,0,1,0)]
```

body pose → ride position → segment force → ride wrench → body derivative. That is
the chain §1 of `own_backend_plan.md` names as the one NetworkDynamics' fixed
depth-2 coreloop cannot express, and it falls out of a plain Kahn layering with no
special case.

At `t = 0`, where both backends hold the identical state, **every** quantity is
bitwise identical: body acceleration, COM, ride-point position and velocity, segment
force and length all `Δ = 0.0`. After 20 steps the two differ by ~1e-6, which is
divergent *adaptive stepping* at `abstol = reltol = 1e-4`, not physics — the
integrators choose different step sequences for a body swinging on a spring.

**The two-component split earns its place.** `RidePoint` (pose → `pos`, `vel`,
`arm`, `height`) and `RideWrench` (`height`, `vel`, `arm`, segment force → `force_out`,
`moment_out`) land in layers 2 and 4 respectively. A single component holding both
would be a cycle, and the runtime would have said so by name.

## 2026-08-07 — a vendored-codegen fix the ride point forced

**Choice:** `drop_unused_inputs` replaces upstream's element-wise drop of inputs no
equation reads.

**Why:** `mtkcompile` requires an array input to be all-in or all-out, *and* requires
every declared input to appear in the system. Upstream drops unused inputs one
element at a time, which produces a partial array and fails the first rule. Keeping
whole arrays satisfies both.

**Also:** `RideWrench` takes a scalar `height` rather than the full `pos`, because
that is all it needs — the air density and wind factor at its own altitude. Fixing
the drop rule made the partial-array case legal, but declaring three slots to use
one is the padding §5 forbids, so the component asks for what it uses.

## 2026-08-07 — joints, and why a clamped body has no state

**Result:** `ElasticJointComponent` and `TimoshenkoJointComponent` are body-to-body
components over the shared `elastic_joint_wrench` / `timoshenko_element_wrench` —
the same builders the monolith's `joint_eqs!` and `timoshenko_joint_eqs!` call.
Both read two poses and emit the restoring wrench on each, so they land in layer 2
and the bodies integrate in the same pass. `test_joint.jl` 23/23,
`test_timoshenko_joint.jl` 20/20.

**A `STATIC` body is now `StaticBody`, not `RigidBody(frozen = true)`.** The frozen
variant kept thirteen states with zero derivatives, and `mtkcompile` turned those
into parameters with no default (`var"com_vel[1]#0"`) — the same aliasing artefact
`monolith_missing_param_defaults` records on the monolith. Pretending to integrate
zeros was the mistake: a clamped body has no state, so it emits its fixed pose
directly. That also deleted the `frozen` branch from `body_integration`.

**Note on the joint poses.** The monolith's joint code carries `ω_b` and rotates it
to world at every use; the body component already outputs `omega_w`, so
`joint_poses` passes it straight through. Same value, one fewer rotation.

## 2026-08-07 — where the §2.2 list stands: 21 of 35

Swept every remaining file. Green by backend swap (21):

`test_point`, `test_segment`, `test_segment_nonlinear`, `test_pulley`,
`test_tether_init`, `test_tether_winch`, `test_rigid_body`, `test_joint`,
`test_timoshenko_joint`, `test_getter_allocations`, `test_heading_calculation`,
`test_match_aero_sections`, `test_multi_section_group`, `test_principal_body_frame`,
`test_quaternion_auto_groups`, `test_quaternion_conversions`,
`test_section_alignment`, `test_static_twist`, `test_weighted_ref_points`,
`test_wing_dynamics`, `test_yaml_bodies`.

**Every one of the 14 remaining failures is the same gap**, not fourteen gaps: the
wing layer. Eleven are literally `body … has type KINEMATIC` (a particle wing, whose
pose is *fitted* from its ref points rather than integrated); one is a beam-anchored
`BODY_STATIC` point (the Hermite ride, `point.joint_idx > 0`); the rest are
assertions that come out zero or `NaN` because no aero force is produced —
`aero_wrench_w`, `point_force_w`, the `test_aero_modes` pose sweep.

So the remaining work is one coherent piece, not a long tail:

1. **KINEMATIC bodies** — `wing_eqs.jl`: pose fitted from the four structural ref
   points, `wing_frame_columns` (already shared) giving `R_b_to_w`, everything else
   derived. No state.
2. **Wing nodes** — `is_wing_node` DYNAMIC points carrying a frozen per-point
   `aero_force_b` and body-frame damping relative to the fitted wing frame. These
   read the wing frame, so they schedule after it.
3. **Beam-anchored ride points** — `point.joint_idx > 0`, placed on a Timoshenko
   joint's Hermite centerline (`beam_hermite_ride_expressions`, already shared).
   Same two-component split as the body ride, reading *two* body poses.
4. **Aero** — `aero_component` per mode, which both backends already share; the
   frozen VSM modes need only the parameter sync, the live modes
   (`ContinuousAero`, `AeroPressure`, `AeroPlate`) need the per-point states as
   inputs.
5. **Twist surfaces** — FIXED/KINEMATIC/DYNAMIC twist, `twist_surface_eqs.jl`.

Note that `test_static_twist`, `test_twist_alignment`, `test_section_alignment` and
`test_match_aero_sections` already pass, so the twist *geometry* is struct-level and
intact; only the twist DOF and the aero coupling are missing.

## 2026-08-07 — fitted (KINEMATIC) wings, and weights in the wiring

**Choice:** `Wiring` carries a weight per source, so an input slot is a *weighted*
sum of output slots rather than a plain sum.

**Why:** A wing's frame is fitted from *weighted* reference points — a blend of
several structural points. With weights in the wiring, a weighted reference costs
one `connect!` per contributing point and no component at all; without them it would
need a per-arity kernel carrying the weights as parameters. One multiply in
`gather!`, and the coupling model is still one sentence.

**Also added:** `KinematicBody` (stateless — its frame comes from
`wing_frame_columns`, its origin pose from a reference point, exactly as
`wing_eqs!`), `WingNodePoint` (a particle plus the frozen per-point aero force and
body-frame damping `point_eqs!` adds), and a `FittedWingReadout` so the state getter
rebuilds such a wing's apparent wind through the shared
`wing_kinematics_from_points!` — without it the VSM solve gets zero wind and
`VortexStepMethod` rejects it.

**Two bugs this found.** `compile_kernel` built its defaults with `map`, which
returns `Vector{Any}` for a kernel with *no* parameters — and `KinematicBody` is the
first such kernel, being pure geometry over its inputs. And `initial_state` guarded
on `type == STATIC` where it meant "does not integrate": a fitted body has no state
either.

**Result:** `test_transform` 1/2 → **181/181**, `test_yaml_weighted_ref` →
**13/13**, `test_wing` 1/2 → 91/99, `test_aero_modes` 0 → 448/500,
`test_continuous_aero` → 16/19. 23 files green, 1119 assertions, nothing regressed.

## 2026-08-07 — `test_helpers` failures are mine, not the backend's

**Note:** `test_helpers.jl` asserts that no `Manifest.toml` sits in the repo root and
that `Manifest-v1.1x.toml.default` is no older than `Project.toml`. Both fail
because *I* copied a `Manifest.toml` in to instantiate this worktree and then added
`OrderedCollections` to `Project.toml`. Neither is a backend defect. The manifest
defaults genuinely need regenerating for the new dependency, which needs a Julia
1.11 and a 1.12 resolve; the stray `Manifest.toml` is gitignored and should be
removed when the worktree is done with.

## 2026-08-07 — closing status

**Green by backend swap, committed in `test/test_scheduled_backend.jl` (23 files,
1119 assertions):** `test_point`, `test_segment`, `test_segment_nonlinear`,
`test_pulley`, `test_tether_init`, `test_tether_winch`, `test_rigid_body`,
`test_joint`, `test_timoshenko_joint`, `test_getter_allocations`,
`test_heading_calculation`, `test_match_aero_sections`, `test_multi_section_group`,
`test_principal_body_frame`, `test_quaternion_auto_groups`,
`test_quaternion_conversions`, `test_section_alignment`, `test_static_twist`,
`test_weighted_ref_points`, `test_wing_dynamics`, `test_yaml_bodies`,
`test_transform`, `test_yaml_weighted_ref`. The monolith runs the same files green.

**Close but not there:** `test_wing` 97/99, `test_aero_modes` 448/500,
`test_continuous_aero` 16/19. All the remaining assertions read `NaN` or zero for a
wing's aerodynamic force, so the models build, initialise and step — what is missing
is that the *force* never reaches the wing. Next step is to find out where it stops:
the frozen `aero_force_b` is a plain struct-field parameter on `WingNodePoint`, so
the generic reader should sync it after each `refresh_aero!`; the likely suspects
are the wing-level `aero_force_b`/`aero_moment_b` never being written back into the
struct by the getter (the monolith has them as observables), and the VSM solve
returning a degenerate solution because `refresh_aero!` runs before the fitted-wing
kinematics are first reconstructed.

**Untried:** `test_beam_replay`, `test_deserialization`, `test_flap_aero`,
`test_flap_beam`, `test_pressure_aero`, `test_profile_law`,
`test_principal_frame_invariance`, `test_twist_alignment`. Several need the same
aero force; `test_flap_beam` and `test_beam_replay` additionally need the
beam-anchored (Hermite) ride point, which is the one structural component still
missing — the same two-part split as the body ride, over
`beam_hermite_ride_expressions`, reading *two* body poses instead of one.

**`test_helpers` is environmental, not a backend defect** — see the entry above.

**§5 still holds by absence.** `grep` over `src/` finds no `wide`, no `WIDE_*`, no
`SimpleGraph`, no `extin`, no `*_at_src`, no `tension_sign_*`, no combined edges,
and no hand-minted parameter names. The runtime is still four concepts, and the
domain boundary at `params.jl` still holds: `kernel.jl` and `runtime.jl` mention no
`SystemStructure`, `Point` or `Body`.

---

## D23 — a point anchored to a beam is a ride point, not a special case

**Decision.** `point.joint_idx > 0` gets the same two-component split as
`point.body_idx > 0`: `HermiteRidePoint` places it on the Timoshenko element's
deflected centerline from both end bodies' poses, `HermiteRideWrench` sends
`(1−s)` of its load to the first end body and `s` to the second. Both go through the
shared `beam_hermite_ride_expressions`, so there is one placement law, not two.

**Why.** This is the component the ND extension could not have: it is a
*feedthrough* vertex reading *two* body poses, which is exactly what ND's depth-2
schedule and its single-source-per-input model forbid (the ND notes call it "the
last structural blocker"). Under a depth-N schedule it needed no new machinery at
all — only the two halves and six `connect!` calls.

**Cost.** Two instances per beam-anchored point. SK100 has 268 of them, so 536
instances of two kernels.

## D24 — aero is a component, not a term inside the points

**Decision.** A wing's aerodynamics is its own component — `ParticleWingAero` for a
`PARTICLE_DYNAMICS` wing (every structural point's pose in, a world force per point
out) and `WingAero` for a `RIGID_DYNAMICS` one (the body's pose in, one wrench about
its COM out). `WingNodePoint` lost its `with_aero` branch; the force now arrives at
`force_in` like any other load.

**Why.** It looks like a cycle and is not: a point's position and a body's pose are
functions of state alone, never of the force they receive, so `input_feeds_output`
is false for `force_in` and the schedule simply gains a layer — structure, wing
frame, aero, derivatives. The monolith always routes particle aero through the aero
subsystem too (`aero_force_point_b` is a *variable* there, bound to a parameter only
when the mode is frozen), so this is the faithful translation and the frozen modes
collapse to one parameter read per point.

This is where ND's design forced `LiveAeroWingNodePoint`: with only two passes, aero
had to be computed *inside* the wing node's f-pass, because `extin` is illegal
anywhere else. Here it is just another vertex.

## D25 — a twist surface's mass arrives as an input

**Decision.** `TwistSurfaceDOF` takes the mass for its `⅓mL²` inertia as
`node_mass_in`, gathered from its own nodes, rather than reading
`points[i].extra_mass` for each of them.

**Why.** A kernel may read only *one* component of a container it is instanced over
— `record_build_index!` enforces it, and it is what makes per-instance parameter
remapping sound. Summing over several points would have broken that rule and needed
either a per-surface kernel (no batching) or a bespoke summing reader. Gathering the
masses through the wiring keeps the rule and costs one constant output per node.

## D26 — `fix_wing` is not carried over

**Decision.** `TwistSurfaceDOF` omits the monolith's `ifelse(fix_wing == true, 0,
…)` freeze on the twist derivatives.

**Why.** `fix_wing` is a parameter nothing sets: `create_sys.jl` declares it
`false`, `SystemStructure` has an unrelated `fix_wing::Bool` field that is never
wired to it, and no code path anywhere writes the symbolic parameter. Carrying dead
branches into a new backend is how they become permanent.

## D27 — `drop_unused_inputs` drops single components

**Decision.** The vendored codegen hides *every* declared input no equation reads,
including one component of a declared vector, instead of keeping a whole vector any
of whose components is read.

**Why.** The original rule was written for a wrench that reads only `pos[3]`, and it
held until `FlapDelta`: a hinge axis of `[0,1,0]` leaves six of a frame's nine
entries unread, and `mtkcompile` then rejected them as "not found in the system".
The system is scalarized before compilation, so a vector's components are already
independent variables and dropping one on its own is well defined. The input *slot*
survives either way — `allinputs` still sizes the buffer — so the wiring is
unaffected and a dropped input simply has no effect, which is exactly true of it.
