# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Include files in dependency order
include("definition_models.jl")  # KiteGeometry models of SAM's own columns
include("types.jl")         # Enums, aero modes, dynamics types, segment springs
include("tube_laws.jl")     # Breukels/Comer-Levy tube rigidities for Timoshenko tubes
include("rigid_body.jl")    # Body models and inertia, a tube's element frame, flaps
include("wing.jl")
include("run_values.jl")  # The quantities of a run and how it reads them
include("system_structure_core.jl")  # system_definition and station mapping
