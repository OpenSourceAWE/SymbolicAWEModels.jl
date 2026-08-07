# Own batched backend (`ScheduledBackend`) — plan

**Branch:** would fork off `nd-bodies-autonomous`. Working doc, not committed.
Supersedes the runtime half of [[networkdynamics_implementation_plan]]; keeps
[[component_unification_plan]] and [[networkdynamics_param_unification_plan]]
unchanged — the shared components and the parameter registry are the assets we
carry over.

**Thesis.** Replace *NetworkDynamics as the runtime substrate* with ~1k lines we own:
a build-time schedule over per-type compiled MTK kernels plus gather/scatter buffers.
Keep MTK authoring, keep the shared `src/components/`, keep the parameter registry,
keep compile-time flat in N. Drop the fork, the uniform-width layer, the
feedthrough ban, and the parallel-edge blocker.

> **STATUS 2026-08-06 — gate 0 ran, and it refutes the compile-time half of this
> document.** ND's compile cost was never ND's architecture; it was two ND *defaults*,
> both fixable from outside. With both fixed, a Timoshenko edge kernel compiles in
> **0.40 s** (was 578 s) and runs **3.8× faster**. See §1.1. Everything in §7 that
> projects a compile-time win over ND is dead; the surviving case for this backend is
> **complexity only**, and must be argued on those terms. Read §1.1 and §7.1 before §7.

---

## 1. Why — the ND tax, measured

ND's coreloop hardcodes a **depth-2 dependency graph** (vertex-g from own states →
edge-g → vertex-f) and **one uniform I/O width per `Network`**. Our physics is
depth-3/4 feedforward and heterogeneous. Every workaround on this branch traces back
to that one mismatch:

| ND constraint | What it cost us |
|---|---|
| no general feedthrough vertices | cluster-vertex absorption (§8.2–8.4), BODY_STATIC points are not components, ride positions recomputed inside *every* incident edge, the two-edge Hermite trick |
| one `vdepth`/`edepth` per network | the whole `wide` / `wide_pose_appendix` / zero-fill layer threaded through point, winch, pulley and segment builders |
| `SimpleGraph`, no parallel edges | task #28 — **the last structural blocker to SK100** |
| `extin` legal only in the f-pass | the aero-hub design was impossible; live aero forced into the wing-node f-pass |
| array-valued I/O | fork patch |
| callable (nonnumeric) params | fork patch |
| edgeless networks | fork patch |
| params addressed by bare field name | `nd_param_disambiguation_issue.md` |

Two of these are worse than "annoying":

**(a) The fork is a release blocker.** `SymbolicAWEModels` cannot be registered against
`~/Code/SciML/NetworkDynamics#nonnumeric-params`. ND hit v1.0 on 2026-07-22 and v1.1.0
on 2026-07-31 — the API is stabilising, which makes upstreaming three invasive codegen
patches *harder*, not easier.

**(b) ND-as-built is probably not flat in N at SK100 scale.** Because ride points
collapse into their body vertex and `SimpleGraph` forbids parallel edges, every vertex
pair carrying more than one segment needs a *combined* multi-segment edge — and
`build_combined_edges` compiles those **individually** ("each combined edge has a
distinct member set … compiled individually rather than reused"). Counting on the real
`data/sk100/struc_geometry_pressure.yaml`, with ride points mapped to their owning body
and Hermite points fanned to both joint bodies:

```
355 segments      → 788 graph edges after ride-point collapse + Hermite fan-out
370 distinct vertex pairs
168 pairs carry >1 segment   (multiplicities 2,3,4,6,7)
```

That is ~168 per-instance `mtkcompile` + codegen passes. Even grouping by member
signature leaves tens of kernels. The property we went to ND *for* is being eaten by
the limitation we cannot fix from outside ND.

**Both (a) and (b) are fixable upstream and neither survives as a rewrite argument.**
(a) is three PRs; (b) is ND issue #194, which is genuinely small. Only the fixed depth-2
schedule and the uniform I/O width are unfixable from outside — and both cost
*complexity*, not speed. That is the honest ranking; the rest of this section is
supporting detail, not independent reasons.

### 1.1 What gate 0 actually found — ND is not slow

Gate 0 timed per-kernel construction on the current branch. First reading, on a
representative 2-body + Timoshenko + ride-point + segment fixture:

| kernel | as shipped | root cause removed |
|---|---|---|
| TimoshenkoJoint edge (wide) | 578.57 s | **0.40 s** |
| wrench segment edge (wide) | 274.63 s | **0.16 s** |
| BodyVertex (wide) | 76.20 s | **0.18 s** |
| SpringDamper edge (narrow) | 7.05 s | **0.13 s** |
| DynamicPoint (narrow) | 0.15 s | **0.19 s** |
| `build_prob!` (whole fixture) | 854.47 s | — |

The naive reading of the first column — ~16 kernels × the mean ≈ 25–50 min — says the
batched backend is *no better than the monolith's ~30 min*, which would have killed this
plan. That reading is wrong. Profiling the 578 s call (`Profile.@profile`, 11 007
main-thread samples) found the time was not in MTK at all:

| frame | samples | share |
|---|---|---|
| `chk_component` (ND `src/doctor.jl:109`) | 9 500 (9 425 self) | **86 %** |
| `simplify_with_mtkcompile` | 487 | 4.4 % |
| `mtkcompile` proper | 370 | 3.4 % |

Two ND defaults, compounding:

1. **`cse=false`** at all three `build_function` sites in `NetworkDynamicsMTKExt.jl`
   (deliberate — PR #231, replaced by the observed-`Let` sharing of PR #234). The
   Timoshenko edge's model is **23 equations / 535 unique DAG nodes**; emitting that DAG
   as a tree yields **276 036 characters** of source. With `cse=true`: 44 127.
2. **`chk_component`**, on by default, then makes Julia type-infer and codegen that
   276 KB body *specialised for the `AccessTracker` debug wrapper* instead of
   `Vector{Float64}`, calls it once to check for out-of-bounds access, and discards the
   compiled result. The production specialisation is compiled again afterwards.

Both are tuned for ND's design point and neither is a bug for its intended users. A
swing equation is ~50 nodes: `cse` has nothing to share, and `chk_component` costs
milliseconds while catching real dim/sym mistakes. This is the **type-to-instance ratio
mismatch** in its sharpest form — ND assumes many instances of few simple types; we are
few instances of many complex types, so an O(kernel size) default that is free for them
costs us ten minutes.

**`cse=true` is also 1.6–9.8× faster at runtime** (Timoshenko warm call 4.69e-7 →
1.24e-7 s), machine-epsilon identical (max rel Δ ≤ 8e-15). LLVM does not recover the
sharing on its own from the larger IR. So ND's default was costing solve time too.

Status: both are applied. `CHECK_COMPONENT[] = false` is set during the network build
(commit `6d2c4c7a`). An opt-in `cse` kwarg plus a `NetworkDynamics.set_cse!` global
default is patched into the fork (`VertexModel`/`EdgeModel` → `generate_io_function` →
the three `build_function` calls), mirroring the existing `mtkcompile` / `set_mtkcompile!`
pattern with `nothing` inheriting the global. Defaults are unchanged, verified: an
explicit `cse=false` reproduces the old output bitwise (max|Δ| = 0.0). It is the
cleanest of the upstream PRs to offer. **Not yet done:** setting `set_cse!(true)` in our
own build path, and running ND's test suite against the patch.

**Consequence for this document: ND costs us complexity, not compile time.** Any
argument for the rewrite has to be made on §1's constraint table alone.

---

## 2. Prior art — does this exist already?

Checked before proposing to build it.

| Option | Verdict |
|---|---|
| **NetworkDynamics.jl** (v1.1.0, active) | The right idea, the wrong constraints. The three gaps are **open upstream issues, unimplemented**: [#194 Support Multigraphs](https://github.com/JuliaDynamics/NetworkDynamics.jl/issues/194) (open since Jan 2025, no comments), [#383 native support around vertex clusters](https://github.com/JuliaDynamics/NetworkDynamics.jl/issues/383) (opened 2026-08-03), [#149 Multilevel networks](https://github.com/JuliaDynamics/NetworkDynamics.jl/issues/149) (open since Oct 2024). Plus our three fork patches. |
| **JuliaSimCompiler.jl** | Does have loop rerolling for repeated components — but it is **proprietary** (JuliaHub registry, not open source), and mentions of it were removed from open MTK. Not an option for a released LGPL package. |
| **GraphDynamics.jl** (Neuroblox) | Closest *runtime* architecture: batched subsystems + connection rules over an ODEProblem, and we already measured ~170× faster TTFS at N=40 with it. But subsystems are **hand-written Julia, not MTK** — as you said, it is not symbolic, so it does not work for us. Our physics is authored symbolically and our parameters are struct-backed symbolic reads. Rejected as a dependency; useful as a reference for the runtime layer. |
| **Native MTK** | Flattens by design: `mtkcompile` builds one `TearingState` over the fully expanded system. There is no per-type / compile-once switch in MTK 11.x. Confirmed again — unchanged since the original evaluation. |

### Answering the earlier "don't roll your own" verdict

`../BeyondTheSim.jl/networkdynamics_backend_plan.md` concluded: *"Roll your own network
compiler on MTK components: that's reinventing JuliaSimCompiler … a large, low-level
compiler project against MTK internals."* That verdict was right about what it was
judging and does not apply here, because **we are not writing a compiler**.

The compiler stays MTK's, called through public API on one small system per type. What
we write is the *runtime*: layout, schedule, buffers, gather/scatter. In ND's own source
that is `coreloop.jl` (270) + `aggregators.jl` (298) + the scheduling half of
`construction.jl` (~300) ≈ **900 lines of ND's 18 700**. Everything else in ND —
initialization (1644), symbolic indexing (2100), callbacks (894), linear analysis
(1500), metadata (1368), show/doctor/spinners — we either do not use or already bypass.

### 2.1 The real competitor: fix ND upstream

After §1.1 this is the option to beat, not a footnote. Full scope:

| item | size | status |
|---|---|---|
| `cse` opt-in kwarg + `set_cse!` global | ~30 lines | **done in fork**, verified: 167× compile, 3.8× runtime, ≤8e-15 Δ, defaults unchanged |
| `CHECK_COMPONENT[] = false` at build | 1 line | **done** (`6d2c4c7a`) |
| #194 multigraph / parallel edges | small | unblocks task #28, kills combined-edge machinery |
| array-valued I/O scalarization | fork patch | needs upstreaming |
| nonnumeric (callable) params | fork patch | needs upstreaming |
| edgeless networks | fork patch | needs upstreaming |

The first three are clean, self-contained PRs a maintainer has no reason to reject; the
last three are the invasive ones and are why the fork exists. **None of it removes the
depth-2 coreloop or the uniform I/O width** — so the `wide` layer, the 268-line
`build_body_mixed_network`, the per-orientation wrench kernels and the BODY_STATIC
absorption all survive this route. That is the actual trade: upstreaming is far cheaper
and keeps a maintained dependency, at the price of permanently carrying the complexity
catalogued in §1's table.

---

## 3. What the backend seam actually requires

Narrow, and already proven by the ND ext:

```julia
function SAM.build_prob!(::SAM.ScheduledBackend, sam; prn = true)
    # → ProbWithAttributes(; prob, param_sync, set_set_values, get_all_state)
end
```

`init!`, `next_step!` and `refresh_aero!` all work on the `SystemStructure`. **Nothing
in the package requires SymbolicIndexingInterface.** `NetworkStateGetter` uses `getu` /
`NWState` purely as a name→index lookup into `u`; when we own the layout that becomes
direct indexing — simpler and faster. This is the confirmation of your instinct: each
component type knows its own symbol→slot map, and cross-component diagnostics
(elevation, azimuth, heading, aoa) become plain Julia functions of the struct, which is
what `scalar_eqs.jl` computes anyway.

---

## 4. Architecture

Four concepts. Nothing more.

```
ComponentKernel   one compiled MTK type: f!, g!, obs!, symbol→slot maps
ComponentInstance one occurrence: state range, param range, output range, input wiring
Schedule          build-time topological layering of instances, batched by kernel
Wiring            precomputed Int vectors: which output slots sum into which input slot
```

### 4.1 Kernel — one `mtkcompile` per type

A kernel is built from an existing shared component unchanged. `inputs` become
parameters that `mtkcompile` allows to be unbound; `outputs` are read back out of the
solved system.

```julia
"""
    ComponentKernel(system, input_syms, output_syms; name)

Compile one component *type* into callable kernels. `rhs!(dstate, state, input,
param, t)` integrates the component's own states; `out!(output, state, input,
param, t)` evaluates its declared outputs. `state_slot`/`output_slot` map a
variable name to its index inside the component's slice.
"""
struct ComponentKernel{F, G, O}
    name::Symbol
    rhs!::F                       # (dstate, state, input, param, t)
    out!::G                       # (output, state, input, param, t)
    obs!::O                       # (buffer, state, input, param, t) — diagnostics
    n_states::Int
    n_inputs::Int
    n_outputs::Int
    n_params::Int
    state_slot::Dict{Symbol, Int}
    output_slot::Dict{Symbol, Int}
    param_registry::SAM.ParamRegistry
end
```

The build is MTK public API — this is the piece to spike first:

```julia
function compile_kernel(system, inputs, outputs; name)
    compiled = mtkcompile(system; inputs, outputs)
    (_, rhs!), state_syms, param_syms, _ =
        generate_control_function(compiled, inputs; simplify = false)
    out! = build_output_function(compiled, outputs, state_syms, inputs, param_syms)
    return ComponentKernel(name, rhs!, out!, ...)
end
```

`build_output_function` substitutes the compiled system's observed equations into the
output expressions and calls `Symbolics.build_function` — the same ~30 lines ND does at
`NetworkDynamicsMTKExt.jl:586-616`. Because we own that call, array-valued I/O
(`pos(t)[1:3]`) needs no scalarization patch, and callable parameters are just an extra
argument instead of a fork.

### 4.2 Instance and layout

```julia
struct ComponentInstance
    kernel::Int                   # index into the kernel table
    states::UnitRange{Int}        # slice of u
    params::UnitRange{Int}        # slice of p
    outputs::UnitRange{Int}       # slice of the output buffer
    inputs::UnitRange{Int}        # slice of the input buffer
end
```

One kernel, many instances. Instances of the same kernel are stored contiguously so a
batch is a strided loop over views — the same trick that makes ND flat in N.

### 4.3 Wiring — gather/scatter instead of uniform-width aggregation

Every input slot is the sum of a set of output slots. That is the whole coupling model,
and it is where the widening layer disappears: **each kernel declares its own input and
output width**, because we index by slot, never by a global depth.

```julia
"""
    Wiring(sources, offsets)

Flattened compressed-row map: input slot `k` is the sum of the output-buffer slots
`sources[offsets[k]:offsets[k+1]-1]`. Built once at assembly time.
"""
struct Wiring
    sources::Vector{Int}
    offsets::Vector{Int}
end

function gather!(input, output, wiring)
    @inbounds for k in eachindex(input)
        total = zero(eltype(input))
        for j in wiring.offsets[k]:(wiring.offsets[k + 1] - 1)
            total += output[wiring.sources[j]]
        end
        input[k] = total
    end
    return nothing
end
```

A segment writing `force` into a point, a wrench edge writing `force` **and** `moment`
into a body, a ride point forwarding an aggregated force into its body — all the same
mechanism, no zero-fill, no superset.

### 4.4 Schedule — the thing ND cannot express

Outputs depend on states and (sometimes) on other components' outputs. Sort that
dependency graph at build time, layer it, and batch each layer by kernel:

```julia
"""
    build_schedule(instances, wiring) -> Vector{Vector{Batch}}

Topologically layer the instances by output dependency, then group each layer by
kernel so one compiled function still runs over many instances. Errors on a cycle
(a genuine algebraic loop), naming the components involved.
"""
```

Our dependency chain is `body states → body pose → ride position → segment force →
body f`: depth 4, acyclic, and *not* representable in ND's fixed depth-2 loop. With a
schedule it is ordinary. Consequences, all of them deletions:

- **BODY_STATIC ride points become first-class components again** — one kernel, ~300
  instances at SK100 — instead of being absorbed and re-derived inside every incident
  edge. Cluster-vertex partitioning, union-find over bodies, §8.2–8.4: gone.
- **Hermite ride points stop needing the two-edge/`extin` trick** — one component that
  reads two body poses.
- **Parallel edges are free.** Edges are a list; nothing keys them by vertex pair.
  Task #28 disappears, and with it the 168 per-instance combined-edge compiles.
- **Per-orientation baked kernels go away** — we order each edge's endpoints when we
  build its wiring, so `body_at_src` is always true.

### 4.5 The RHS

```julia
function (rhs::ScheduledRHS)(du, u, p, t)
    buffers = get_buffers(rhs, eltype(u))
    for layer in rhs.schedule
        for batch in layer
            kernel = rhs.kernels[batch.kernel]
            for inst in batch.instances
                gather!(view(buffers.input, inst.inputs), buffers.output, rhs.wiring)
                kernel.out!(view(buffers.output, inst.outputs),
                            view(u, inst.states), view(buffers.input, inst.inputs),
                            view(p.numeric, inst.params), t)
            end
        end
    end
    for batch in rhs.state_batches
        kernel = rhs.kernels[batch.kernel]
        for inst in batch.instances
            kernel.rhs!(view(du, inst.states), view(u, inst.states),
                        view(buffers.input, inst.inputs),
                        view(p.numeric, inst.params), t)
        end
    end
    return nothing
end
```

`get_buffers` keys a cache on `eltype(u)` so ForwardDiff duals work if we ever want an
AD Jacobian. Today's solver is `FBDF(; autodiff = AutoFiniteDiff())`, so plain `Float64`
is enough for v1.

### 4.6 Parameters — reuse, do not rewrite

`ParamView` / `ParamRegistry` / `PathReader` / `ParamGroup` / `sync_params!` are already
backend-agnostic and stay verbatim. Only the *address* changes, from `VIndex(i, :sym)`
to `(instance, slot)` → a flat index. Two ND-specific pains vanish:

- **callable params**: `p` is our own struct, so `p.callables[i]` holds the aero polars
  with no codegen patch;
- **name aliasing** (`nd_param_disambiguation_issue.md`): a param address is
  `(instance, slot)`, so two reads of `anchor_b` on different instances can no longer
  collapse onto one symbol.

```julia
struct ScheduledParams{C}
    numeric::Vector{Float64}
    callables::C                  # per-instance NamedTuple of polars / rigidity laws
end
```

### 4.7 Getters

Per-type, as you proposed. `kernel.obs!` gives the component-local diagnostics;
cross-component quantities are plain functions of the struct.

```julia
function (getter::ScheduledStateGetter)(integ, ss)
    for (instance, point) in getter.dynamic_points
        base = first(instance.states) - 1
        point.pos_w .= @view integ.u[(base + 1):(base + 3)]
        point.vel_w .= @view integ.u[(base + 4):(base + 6)]
    end
    ...
end
```

Observables that depend on aggregated inputs (`wing_aero_force_b`, `pose_R`) need the
output buffers current — run one RHS evaluation after the step, or reuse the last one.

---

## 5. Worked mini-example (the spike target)

Mass on a spring, from the shared components unchanged:

```julia
point_kernel = compile_kernel(SAM.DynamicPoint(sam, params, 1; name = :point),
                              [:force_in, :mass_in], [:pos, :vel]; name = :point)
segment_kernel = compile_kernel(SAM.SpringDamperSegment(sam, params, 1; name = :seg),
                                [:src_pos, :src_vel, :dst_pos, :dst_vel],
                                [:src_force, :src_mass, :dst_force, :dst_mass];
                                name = :segment)

assembly = ScheduledAssembly()
anchor = add_instance!(assembly, static_kernel)
mass   = add_instance!(assembly, point_kernel)
spring = add_instance!(assembly, segment_kernel)
connect_output!(assembly, anchor => :pos,  spring => :src_pos)
connect_output!(assembly, mass   => :pos,  spring => :dst_pos)
accumulate!(assembly, spring => :dst_force, mass => :force_in)
accumulate!(assembly, spring => :dst_mass,  mass => :mass_in)

prob = ODEProblem(finalize(assembly), u0, tspan, p0)
```

Acceptance for the spike: `max|Δaccel|` vs the monolith on the existing mass-on-spring
parity test, and no fork of anything.

---

## 6. What we delete

Of the 2928-line ND ext, roughly:

| Survives with edits (~1200) | Deleted (~700) | New (~900–1200) |
|---|---|---|
| `classify_segments`, edge/vertex enumeration, `record_*_params!`, `set_const_params!`, state setters, wing geometry, control setter | `WIDE_*` constants + `wide` flag threading, BODY_STATIC absorption, per-orientation kernels, combined-edge machinery, `extin` plumbing, `MixedEdgeInfo` | kernel compilation, layout, schedule, wiring, RHS, buffer cache, getter |

`src/components/components.jl` (1349 lines of shared physics) is **untouched**, which is
why every feature keeps its existing monolith parity oracle.

---

## 7. Honest compile-time estimate for SK100

> **§7.1 — superseded by measurement.** This section was written before gate 0 and
> argues a 10–30× compile win over the monolith. The *monolith* comparison still stands.
> The implied win over **ND** does not: §1.1 measures ND kernels at 0.13–0.40 s each
> once `cse=true` and `CHECK_COMPONENT[]=false`, so a ~16-kernel SK100 build is
> **seconds of codegen plus a one-off ~60 s warm-up of the MTK/Symbolics machinery**
> (paid once per session, not per kernel — visible as the first kernel built always
> costing ~60 s regardless of which one it is).
>
> An own backend calls the same `mtkcompile` and the same `build_function`, so it lands
> in the same place. **Between ND-with-both-fixes and an own backend, expect no
> compile-time difference at all.** The estimate below should be read as
> *batched-vs-monolith*, which is the architecture both options share — it is not a
> reason to prefer one of them over the other.
>
> Also corrected here: an earlier draft cited 45 436 expression nodes for the Timoshenko
> edge. That was a tree walk over a hashconsed DAG with no memoisation. The real model
> is **535 unique nodes** — our physics is not bloated. The 45 k figure only describes
> what `cse=false` *emits*.

### The model (counted from `data/sk100/struc_geometry_pressure.yaml`)

| | count |
|---|---|
| points | 384 — 292 `BODY_STATIC` (24 on bodies, 268 Hermite receivers on joints), 88 `DYNAMIC`, 4 `STATIC` |
| segments | 355 (271 touch a body or joint-riding point) |
| bodies | 50 |
| Timoshenko joints | 49 |
| twist surfaces | 12 |
| tethers / winches | 4 / 3 |
| wing | 1, `n_panels = 44` |

Runtime cache name for the built model confirms the same order:
`404pnt_379seg_12grp_1wng_3wch_51bdy`.

### Current monolith cost (measured, `../BeyondTheSim.jl/progress.md`)

| stage | time |
|---|---|
| `mtkcompile` | ~63 s |
| `ODEProblem` / codegen | ~141 s |
| getters | ~7 s |
| **RHS JIT** (LLVM -O2, one fused RuntimeGeneratedFunction) | **~17 min** |
| **total first build** | **~30 min** |

~970 states, ~32k inlined observed equations, and a **2.6 GB** serialized model cache.

### The estimate

The point is the ratio of **instances to types**. SK100 is close to the ideal case:
50 bodies are one type, 268 Hermite receivers are one type, 355 segments are two or
three types, 88 dynamic points are one type. A batched build compiles roughly:

| kernel | count |
|---|---|
| dynamic / static / body / static-body vertex | 4 |
| BODY_STATIC ride point, Hermite ride point | 2 |
| plain / tether / pulley segment | 3 |
| wrench segment (body-touching) | 1 |
| Timoshenko joint, elastic joint | 2 |
| winch, pulley | 2 |
| twist surface, aero receiver / wing aero | 2–3 |
| **total** | **~16–17** |

Each is a small system (tens of equations) compiled in isolation. Against ~32k fused
equations today, and given LLVM's super-linear cost in single-function size, the RHS
build should drop by **more** than the equation-count ratio.

| | monolith today | batched (estimate) |
|---|---|---|
| symbolic (`mtkcompile` + codegen) | ~210 s | ~20–60 s (16 small systems) |
| RHS JIT | ~17 min | ~30–120 s (16 small functions) |
| assembly / wiring / schedule | — | seconds (no codegen) |
| **RHS build total** | **~20 min** | **~1–3 min → 10–30×** |
| **full TTFS** | **~30 min** | **~3–6 min → 5–10×** |

**The bigger prize is the scaling, not the constant.** Refining the beam (more chord
elements, more panels, more receivers) currently costs super-linear compile time; on a
batched backend it costs ~zero. That is what makes SK100 iterable.

### What does *not* speed up — be honest about this

- Package load, VSM / NeuralFoil polar setup, struct construction, `init!` solve. These
  are a fixed several minutes and dominate the "5–10×" row above.
- **Per-eval runtime may get slightly worse.** The monolith is one fused function with
  global CSE; batched execution pays buffer traffic and loses cross-component CSE.
  Expect ~1–3× slower per RHS call, partly recovered at large N by cache locality.
- **The Jacobian, unless we do something about it.** FBDF + `AutoFiniteDiff` with no
  `jac_prototype` anywhere means ~970 RHS evaluations per Jacobian. This dominates warm
  solve time and is untouched by any compile-time work.

### Two secondary wins worth more than they look

1. **Sparsity.** We know the coupling graph exactly, so a `jac_prototype` + colorvec is
   nearly free to emit. ~970 dense columns → plausibly 20–40 colors, i.e. an order of
   magnitude fewer RHS evaluations per Jacobian. This may beat the compile win in
   day-to-day use. Not in v1, but only this backend makes it cheap.
2. **The 2.6 GB model cache.** That is one giant serialized generated function. A
   per-type build should serialize to megabytes, making the cache actually cheap to
   load, version and ship.

### Confidence

Low-to-medium on the absolute numbers, high on the direction, and **only against the
monolith** — see §7.1. Nobody has benchmarked a batched backend at SK100 scale;
`progress.md` still lists it as unconfirmed. The supporting evidence is the
GraphDynamics prototype (~170× TTFS at N=40, trivial kernels), the structural argument,
and now §1.1's per-kernel measurements, which are the first real data point: small
kernels do compile in fractions of a second, so the batched *architecture* delivers what
this section claims. ND already implements that architecture.

---

## 8. Risks

1. **MTK codegen for input-carrying subsystems** is the one genuine unknown. ND's
   `MTKExt_utils.jl` is 1264 lines of defaults/metadata/aliasing pain, some of which is
   MTK's fault and would be inherited. Mitigation: gate 0 + gate 1 before committing.
   *Partly answered by §1.1:* `mtkcompile(sys; inputs, outputs)` on the wide Timoshenko
   edge is 3.4 % of the build and reduces cleanly (23 eqs → 0 states + 54 observed), so
   the MTK half is sound. The 1264 lines are ND's metadata/aliasing repair work, and we
   would have to re-derive our own version of it — the `fix_metadata!` "could NOT be
   fixed" warnings on every wide kernel are that layer failing loudly today.
2. **Flatness leaks.** Any per-instance kernel silently restores O(N) compiles — exactly
   what happened to ND via combined edges. Mitigation: assert a kernel-count ceiling in
   the build and log the kernel table.
3. **The aero kernel may not be small.** If a wing-level kernel carries all 44 panels it
   becomes the new bottleneck. Mitigation: keep the per-node slice decomposition already
   designed in `networkdynamics_aero_component.md`.
4. **We own the scheduler.** ~1k lines of our own bugs, in the layer where bugs are
   subtle. Mitigation: the parity oracle already exists for every feature.
5. **Two throwaway weeks** if gate 1 fails. Bounded, and #28 still needs doing either way.

---

## 9. Phased plan

Each gate is a stop/go with a measurement, not a vibe.

- **Gate 0 — DONE (2026-08-06), result in §1.1.** It did its job and inverted the
  premise: the kernels are fast, ND's defaults were slow, and the compile-time
  justification for this project is gone. Two process lessons worth keeping — the
  "repeat call" column was meaningless (Symbolics hashconsing makes a second identical
  build ~13× cheaper, so only cold first calls count), and **profile before
  extrapolating**: one `Profile.@profile` answered in ten minutes what three rounds of
  plausible hypotheses (`assume_io_coupling` inflation, `_get_formulas` inlining,
  `mtkcompile` cost) all got wrong.
- **Gate 0b — NEW, do this before gate 1.** Wire `cse=true` + `CHECK_COMPONENT[]=false`
  into the ND ext, rebuild SK100 end-to-end, and measure real TTFS. If ND-with-fixes
  builds SK100 in minutes, the remaining case for this backend is complexity alone, and
  the honest comparison is *this plan* vs *upstreaming #194 + the fork patches*. Cheap,
  and it decides the project.
- **Gate 1 — the codegen spike (≈1 session).** `compile_kernel` for `DynamicPoint` and
  `SpringDamperSegment`, hand-wired, mass-on-spring parity vs the monolith. Proves array
  I/O and struct-backed params with no fork. **Stop/go on the whole project.** Note it
  can no longer be justified on compile time — it must be justified by the §1 constraint
  table.
- **Phase 2 — runtime core.** Layout, wiring, schedule, RHS, buffer cache, param sync,
  state getter. Re-run the point/segment/pulley/winch parity suite.
- **Phase 3 — port the structural layer.** Bodies, joints, ride points, Hermite points,
  wrench segments. Mostly *deleting* ND workarounds against existing oracles.
  Acceptance: `test_rigid_body.jl` 19/19 and `test_timoshenko_joint.jl` 20/20 by backend
  swap, matching what ND already achieves.
- **Phase 4 — aero + twist**, then SK100 end-to-end. Acceptance: the north-star test —
  every monolith test passing with only `backend=` swapped.
- **Phase 5 — the payoff.** Benchmark TTFS vs the monolith at several N. Then graph
  sparsity for the Jacobian.

Keep the ND backend alive until phase 3 passes; it is a second oracle.

---

## 10. Open questions

0. **The question this document now turns on:** with ND's two defaults fixed, is the
   remaining complexity in §1's table worth a ~1k-line rewrite plus owning the
   scheduler? Gate 0b prices the left-hand side; §2.1 prices the alternative. Do not
   start gate 1 before answering it — the compile-time argument that motivated this
   plan no longer exists.
1. Does `generate_control_function` keep array-valued inputs intact, or does it
   scalarize the way ND's unpatched path did? (Gate 1 answers this.)
2. Callable params: extra `build_function` argument, or a typed `p` struct the generated
   code indexes? Depends on what MTK's codegen will accept as a parameter object.
3. Do we keep a `Graphs.jl` graph at all? Probably not — instances plus wiring are
   enough, and dropping it is what makes parallel edges free.
4. Retire the ND backend after phase 4, or keep it? Keeping it means maintaining the
   fork; retiring it means one fewer oracle.
    Retire ND immediately, strip out the ugly ND stuff immediately, make the code
    readable and nice again.

## References

- ND open gaps: [#194 multigraphs](https://github.com/JuliaDynamics/NetworkDynamics.jl/issues/194),
  [#383 vertex clusters](https://github.com/JuliaDynamics/NetworkDynamics.jl/issues/383),
  [#149 multilevel networks](https://github.com/JuliaDynamics/NetworkDynamics.jl/issues/149)
- [GraphDynamics.jl](https://github.com/Neuroblox/GraphDynamics.jl) — batched runtime,
  non-symbolic
- [JuliaSimCompiler loop rerolling](https://docs.sciml.ai/SciMLBenchmarksOutput/dev/ModelingToolkit/RCCircuit/)
  — proprietary
- [MTK compile time for repeated components](https://discourse.julialang.org/t/modelingtoolkit-jl-performance-for-large-models-with-similar-components/82442)
- MTK API used: `mtkcompile(sys; inputs, outputs)`
  (`ModelingToolkitBase/src/systems/systems.jl:84`), `generate_control_function`
  (`ModelingToolkitBase/src/inputoutput.jl:241`)
- §1.1 sources: ND `src/doctor.jl:89` (`chk_component`), `src/NetworkDynamics.jl:137`
  (`CHECK_COMPONENT`), `ext/NetworkDynamicsMTKExt.jl` `build_function` sites,
  `simplify_with_mtkcompile` (`ext/MTKExt_simplification.jl:887`), `_get_formulas`.
  ND history: PR #231 `hw/disable_cse` (27 Mar 2025), PR #234 `hw/cse_observed`
  (28 Mar 2025, "keep observed separate for faster codegen and execution").
