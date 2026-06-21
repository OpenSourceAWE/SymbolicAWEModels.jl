#!/usr/bin/env bash
# Benchmark next_step! against a specific git ref and append the result to
# bench/results.md (in the main checkout, not the worktree).
#
#   bench/run_on_ref.sh <git-ref> [julia-args...]
#
# The ref is checked out into a throwaway git worktree, the current bench
# script + 2plate_kite data are copied in (old refs may lack them), and the benchmark
# runs against THAT ref's SymbolicAWEModels source (examples env points at the
# worktree root via [sources]). First run on a fresh ref precompiles — minutes.
#
# Note: this script and the model-construction API must be compatible with the
# ref. Refs predating the current public API will fail to build; that is
# reported, not silently skipped.
set -euo pipefail

ref="${1:?usage: run_on_ref.sh <git-ref>}"
shift || true

repo_root="$(git rev-parse --show-toplevel)"
results="$repo_root/bench/results.md"
worktree="$(mktemp -d)/symawe-${ref//\//_}"

cleanup() { git -C "$repo_root" worktree remove --force "$worktree" 2>/dev/null || true; }
trap cleanup EXIT

git -C "$repo_root" worktree add --detach "$worktree" "$ref"

mkdir -p "$worktree/bench"
cp "$repo_root/bench/bench_next_step.jl" "$worktree/bench/"
if [ ! -d "$worktree/data/2plate_kite" ]; then
    mkdir -p "$worktree/data"
    cp -r "$repo_root/data/2plate_kite" "$worktree/data/"
fi

cd "$worktree"
BENCH_REF="$ref" BENCH_RESULTS="$results" BENCH_DATE="$(git -C "$repo_root" log -1 --format=%cs "$ref")" \
    julia --project=examples "$@" bench/bench_next_step.jl
