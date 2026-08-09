# Issue draft — NetworkDynamics.jl

**Title:** `cse=false` in MTK codegen is a large regression for big component models

---

`generate_io_function` passes `cse=false` to all three `build_function` calls (`f`, `g`, `obsf`). For small component models that's clearly the right call — the observed equations already act as shared named intermediates, and #231's benchmarks show the win. For large component models it inverts sharply.

Our components are structural-mechanics kernels (corotational Timoshenko beam elements, rigid-body wrenches). A representative edge model:

- 23 equations, **535 unique nodes** in the hashconsed DAG
- emitted with `cse=false`: **276,036 characters** of source for `g`
- emitted with `cse=true`: **44,127 characters**

Without CSE each shared subexpression is re-emitted once per reference, so a 535-node DAG expands to a ~45,000-node tree. Measured on ND v1.1.0 (`mtkcompile=true`, `CHECK_COMPONENT[]=false` so this isolates codegen):

| | `cse=false` | `cse=true` |
|---|---|---|
| `build_function` | 0.60 s | 0.09 s |
| eval + first call (LLVM) | 3.14 s | 0.11 s |

End to end via `EdgeModel`/`VertexModel`, `cse=true` compiles 2–28× faster across our kernel set. Runtime also improves 1.6–9.8× (crude `@elapsed` timings on sub-microsecond calls, so treat the factors as approximate but the direction is clear) — LLVM does not recover the sharing from the larger IR on its own. Outputs agree to ≤8e-15 relative.

Would you take a PR adding an opt-in `cse` kwarg on `VertexModel`/`EdgeModel`, threaded to `generate_io_function`, plus a `NetworkDynamics.set_cse!` global mirroring `set_mtkcompile!` (`cse=nothing` inherits)? Defaults unchanged — verified that an explicit `cse=false` reproduces current output bitwise. I have this working locally and can open it straight away.
