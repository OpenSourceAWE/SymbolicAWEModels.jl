# Issue: ND backend can't disambiguate two reads of the same field on different instances within one kernel

## Summary

The network backend names each flat parameter by its **bare field name**
(`leaf_symbol(path)` in `flat_params.jl` — trailing indices stripped) so that one
generic kernel is shared across all instances of a component type, with per-instance
values supplied at sync time. This breaks down when a **single kernel reads the same
struct field on two different container instances**: both reads memoize to the one
symbol, and the second silently aliases onto the first.

## Where it bites today

The **dual-wrench edge** (`network_dual_wrench_segment`, `ext/…Ext.jl`) is a body↔body
spring whose *both* endpoints are `BODY_STATIC` ride points. It needs each endpoint's
`anchor_b` (the fixed body-frame offset) as a parameter. Reading both through ParamView:

```julia
params.points[ride_src].anchor_b   # → network symbol :anchor_b
params.points[ride_dst].anchor_b   # → network symbol :anchor_b  (SAME — aliases!)
```

Both collapse to `:anchor_b` because `leaf_symbol((:points, ride_src, :anchor_b))` and
`leaf_symbol((:points, ride_dst, :anchor_b))` are both `:anchor_b`. The two physically
distinct anchors would become one parameter.

## Current workaround (accepted, kept deliberately)

The dual-wrench kernel mints explicitly-named, pre-scalarized params instead of routing
through ParamView:

```julia
anchor_src = [SAM.make_param(Symbol(:anchor_b_src_, k), 0.0) for k in 1:3]
anchor_dst = [SAM.make_param(Symbol(:anchor_b_dst_, k), 0.0) for k in 1:3]
```

and binds them by hand in `record_mixed_edge_params!` with the correct per-endpoint
`PathReader((:points, info.ride_src/ride_dst, :anchor_b, k))`. Two lines, correct,
documented at the call site. Validated: dual-wrench parity Δ~4e-12.

## What this is NOT

It is **not** the "raw `make_array_param` trips the fork's vector-param scalarizer"
problem. Arrays read through ParamView work: the *single*-endpoint wrench edge reads
`params.points[ride_idx].anchor_b` directly and the fork scalarizes the `anchor_b[1:3]`
registry param to `anchor_b_1/_2/_3` cleanly. The dual-wrench issue is purely the
**same-field-two-instances name collision**, orthogonal to arrays. Scalarizing array
leaves on read would not fix it (both endpoints would still scalarize to the same names).

## Proper fix (deferred — cost > benefit right now)

Give the network param path a **role disambiguator**: when one kernel reads the same
`(container, field)` more than once, suffix the symbol with a stable-across-instances
role tag (e.g. `anchor_b_src` / `anchor_b_dst`), not the raw container index (which
varies per edge and would break kernel sharing). This touches two places that must agree:

- `ParamView`/`PathView` `getproperty` in `flat_params.jl` — detect the repeat read and
  emit the suffixed symbol + a role-aware cache key.
- `replay_fields!` in the ext — bind the matching suffixed symbol per endpoint.

This adds permanent complexity to a param path shared by **both** backends for the
benefit of exactly one kernel. Until a second consumer needs it, the two-line explicit
workaround is the better trade. Revisit if/when more edges read multiple same-typed
instances (e.g. multi-point wrench aggregators).

## Related

- `flat_params.jl`: `leaf_symbol`, `param_symbol_name(NetworkBackend, …)`,
  `param_cache_key(NetworkBackend, …)` — the bare-field-name policy.
- `ext/…Ext.jl`: `network_dual_wrench_segment`, `record_mixed_edge_params!`.
- Memory `networkdynamics_param_unification` — the backend-tagged ParamView design.
