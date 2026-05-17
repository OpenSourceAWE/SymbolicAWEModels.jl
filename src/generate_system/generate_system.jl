# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Implementation of the ram air wing model using ModelingToolkit.jl
# This directory contains the symbolic equation generation for the AWE model.

# Utilities and accessors first (no dependencies)
include("helpers.jl")
include("accessors.jl")

# Component equations (can be in any order - no inter-dependencies)
include("point_eqs.jl")
include("group_eqs.jl")
include("segment_eqs.jl")
include("pulley_eqs.jl")
include("winch_eqs.jl")
include("tether_eqs.jl")

# Higher-level equations
include("wing_eqs.jl")
include("vsm_eqs.jl")
include("plate_eqs.jl")
include("scalar_eqs.jl")

# Main entry point last (depends on all above)
include("create_sys.jl")
