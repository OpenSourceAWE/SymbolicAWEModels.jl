# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Include files in dependency order
include("definition_models.jl")  # KiteGeometry models of SAM's own columns
include("types.jl")         # Enums, aero modes, dynamics types, segment springs
include("tube_laws.jl")     # Breukels/Comer-Levy tube rigidities for TimoshenkoTube
include("rigid_body.jl")    # Body models, tube elements, flaps
include("wing.jl")
include("system_state.jl")  # SystemState: what a run holds of each component
include("system_structure_core.jl")  # system_definition and station mapping
include("utilities.jl")     # Helpers, validation, state management
