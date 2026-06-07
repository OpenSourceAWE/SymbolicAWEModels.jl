# Copyright (c) 2025 Bart van de Lint, Jelle Poland
# SPDX-License-Identifier: LGPL-3.0-only

# Generate all documentation figures.
# Run with: julia --project=docs docs/generate_figures.jl
#
# This script produces PNGs in docs/src/assets/ by:
#   1. include()-ing the Literate tutorial (hidden save() lines run)
#   2. generating standalone figures for non-Literate pages

import GLMakie
GLMakie.activate!(; visible=false)

using KiteUtils: init!, next_step!, update_sys_state!
using SymbolicAWEModels: Point

# --- Section 1: Literate tutorial figures ---
include(joinpath(@__DIR__, "src", "literate", "tutorial_julia.jl"))

# --- Section 2: Standalone figures (non-Literate pages) ---
# SymbolicAWEModels already loaded by the tutorial include above

ASSETS = joinpath(@__DIR__, "src", "assets")

# 2-plate kite from YAML (for examples.md)
pkg_root = dirname(@__DIR__)
set_data_path(joinpath(pkg_root, "data", "2plate_kite"))

struc_yaml = joinpath(get_data_path(), "rigid_structural_geometry.yaml")

set_2p = Settings("system.yaml")
set_2p.solver = "FBDF"
vsm_set_2p = VortexStepMethod.VSMSettings(
    joinpath(get_data_path(), "vsm_settings.yaml"); data_prefix=false)

sys_2p = load_sys_struct_from_yaml(struc_yaml;
    system_name="2plate_kite", set=set_2p, vsm_set=vsm_set_2p)
scene = plot(sys_2p)
GLMakie.save(joinpath(ASSETS, "2plate_kite_structure.png"), scene)

println("All figures generated in $ASSETS")
