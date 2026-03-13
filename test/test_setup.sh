#!/bin/bash

# SPDX-FileCopyrightText: 2026 Bart van de Lint
#
# SPDX-License-Identifier: MPL-2.0

# Integration test for end-user setup.
# Run from the repo root: bash test/test_setup.sh
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PASS=0
FAIL=0

pass() { echo "  PASS: $1"; ((++PASS)); }
fail() { echo "  FAIL: $1"; ((++FAIL)); }

cleanup() {
    [[ -n "${TMPDIR_USER:-}" ]] && rm -rf "$TMPDIR_USER"
}
trap cleanup EXIT

# Use xvfb-run if available (headless CI), otherwise run directly
if command -v xvfb-run &>/dev/null; then
    JULIA="xvfb-run -a julia"
else
    JULIA="julia"
fi

echo "=== End-user setup ==="
TMPDIR_USER=$(mktemp -d)
echo "  tmpdir: $TMPDIR_USER"
cd "$TMPDIR_USER"

# Simulate: mkdir my_project && cd my_project && julia --project=.
# Then: pkg> add SymbolicAWEModels (use dev for unreleased)
# Then: using SymbolicAWEModels; SymbolicAWEModels.init_module()
$JULIA --project=. -e '
    using Pkg
    Pkg.develop(path="'"$REPO_ROOT"'")
    using SymbolicAWEModels
    SymbolicAWEModels.init_module()
' 2>&1 && pass "init_module()" || fail "init_module()"

# Verify example files were copied
for f in menu.jl hanging_mass.jl catenary_line.jl \
         pulley.jl saddle_form.jl \
         coupled_2plate_kite.jl coupled_2plate_kite_linear_vsm.jl \
         coupled_tether_deflection.jl coupled_realtime_visualization.jl \
         coupled_linearize.jl coupled_simple_lin_model.jl \
         static_load_2plate_kite.jl sam_tutorial.jl; do
    [[ -f "examples/$f" ]] && pass "copied examples/$f" \
                           || fail "copied examples/$f"
done

# Verify data directories were copied
for d in data/2plate_kite data/base data/saddle_form; do
    [[ -d "$d" ]] && pass "copied $d/" || fail "copied $d/"
done

# Verify all examples use GLMakie
for f in menu.jl hanging_mass.jl catenary_line.jl \
         pulley.jl saddle_form.jl coupled_2plate_kite.jl; do
    if grep -q "using GLMakie" "examples/$f"; then
        pass "$f uses GLMakie"
    else
        fail "$f uses GLMakie"
    fi
done

# Run examples (ordered from simple to complex)
for f in hanging_mass.jl catenary_line.jl \
         pulley.jl saddle_form.jl sam_tutorial.jl \
         coupled_tether_deflection.jl \
         coupled_2plate_kite.jl \
         coupled_2plate_kite_linear_vsm.jl \
         coupled_linearize.jl \
         coupled_simple_lin_model.jl \
         static_load_2plate_kite.jl; do
    echo "  Running $f..."
    $JULIA --project=. -e '
        include("examples/'"$f"'")
    ' 2>&1 && pass "run $f" || fail "run $f"
done

echo ""
echo "=== Results: $PASS passed, $FAIL failed ==="
[[ $FAIL -eq 0 ]] && exit 0 || exit 1
