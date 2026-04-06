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
         coupled_linearize.jl \
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
         static_load_2plate_kite.jl; do
    echo "  Running $f..."
    $JULIA --project=. -e '
        include("examples/'"$f"'")
    ' 2>&1 && pass "run $f" || fail "run $f"
done

# Test README minimal pendulum example
echo "  Running README pendulum example..."
$JULIA --project=. -e '
    using SymbolicAWEModels
    SymbolicAWEModels.init_module(; force=false)
    set_data_path("data/base")

    set = Settings("system.yaml")
    set.v_wind = 0.0

    points = [
        Point(:anchor, [0, 0, 0], STATIC),
        Point(:mass, [0, 0, -50], DYNAMIC; extra_mass=1.0),
    ]
    segments = [Segment(:spring, set, :anchor, :mass, BRIDLE)]
    transforms = [Transform(:tf, deg2rad(-80), 0.0, 0.0;
        base_pos=[0, 0, 50], base_point=:anchor, rot_point=:mass)]

    sys = SystemStructure("pendulum", set; points, segments, transforms)
    sam = SymbolicAWEModel(set, sys)
    init!(sam)

    for _ in 1:100
        next_step!(sam)
    end
' 2>&1 && pass "README pendulum example" || fail "README pendulum example"

# Test README 2plate kite example
echo "  Running README 2plate kite example..."
$JULIA --project=. -e '
    using SymbolicAWEModels, VortexStepMethod
    SymbolicAWEModels.init_module(; force=false)

    set_data_path("data/2plate_kite")

    struc_yaml = joinpath(get_data_path(), "quat_struc_geometry.yaml")
    aero_yaml = joinpath(get_data_path(), "aero_geometry.yaml")
    update_aero_yaml_from_struc_yaml!(struc_yaml, aero_yaml)

    set = Settings("system.yaml")
    vsm_set = VortexStepMethod.VSMSettings(
        joinpath(get_data_path(), "vsm_settings.yaml"); data_prefix=false)

    sys = load_sys_struct_from_yaml(struc_yaml;
        system_name="2plate_kite", set, vsm_set)

    sam = SymbolicAWEModel(set, sys)
    init!(sam)

    l0_left = sam.sys_struct.segments[:kcu_steering_left].l0
    l0_right = sam.sys_struct.segments[:kcu_steering_right].l0

    for step in 1:600
        t = step * (10.0 / 600)
        ramp = clamp(t / 2.0, 0.0, 1.0)
        sam.sys_struct.segments[:kcu_steering_left].l0 = l0_left - 0.1 * ramp
        sam.sys_struct.segments[:kcu_steering_right].l0 = l0_right + 0.1 * ramp
        next_step!(sam; dt=10.0/600, vsm_interval=1)
    end
' 2>&1 && pass "README 2plate kite example" || fail "README 2plate kite example"

echo ""
echo "=== Results: $PASS passed, $FAIL failed ==="
[[ $FAIL -eq 0 ]] && exit 0 || exit 1
