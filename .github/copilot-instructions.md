# Copilot Instructions for SymbolicAWEModels.jl

## Project Context
- This repository is a Julia package for symbolic airborne wind energy models.
- Prefer correctness and clarity over clever optimizations.
- Keep changes minimal and focused on the requested task.

## Language and Style
- Use Julia 1.12-compatible syntax unless the surrounding file clearly targets another version.
- Follow existing code style in the touched file (naming, spacing, imports, and structure).
- Do not rename public APIs unless explicitly requested.
- Add short comments only when logic is non-obvious.

## Dependencies and Architecture
- Reuse existing modules and utilities before adding new abstractions.
- Avoid adding new package dependencies unless necessary and justified.
- Keep `src/SymbolicAWEModels.jl` as the package entry point and update exports intentionally.

## Testing and Validation
- When changing behavior, add or update tests under `test/`.
- Run targeted tests first (closest affected test files), then broader tests if needed.
- Do not ignore failing tests; either fix them or explain why they are unrelated.
- To run a single test, use `include("test/<test_file>.jl")` from the Julia REPL, or use `./bin/run_julia test/<test_file>.jl` from the command line.

## Documentation
- Update docs in `docs/src/` when user-facing behavior, APIs, or workflows change.
- Keep examples runnable and consistent with current APIs.

## Safety for Edits
- Preserve existing unrelated changes in the working tree.
- Avoid large refactors for small requests.
- If requirements are ambiguous, state assumptions briefly in the final response.
