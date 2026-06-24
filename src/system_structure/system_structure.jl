# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

"""
System Structure

Defines the physical structure of the kite system, including:
- Basic types: Point, TwistSurface, Segment, Pulley, Tether, Winch, Transform
- Wing types: BaseWing, VSMWing with VSM aerodynamics coupling
- SystemStructure container and initialization

Files are organized as:
- types.jl: Enums and struct definitions (Point, Segment, Transform, etc.)
- wing.jl: Wing types and VSM-related code
- system_structure_core.jl: SystemStructure type and constructor
- transforms.jl: Heading/rotation functions, reinit!/reposition!
- utilities.jl: Helper functions and state management
"""

# Include files in dependency order
include("types.jl")         # Enums, Point, TwistSurface, Segment, Pulley, Tether, Winch, Transform
include("rigid_body.jl")    # Body (embedded in Wing and standalone)
include("wing.jl")          # AbstractWing, Wing (embeds a Body)
include("named_collection.jl")  # NamedCollection wrapper for symbol-based indexing
include("system_structure_core.jl")  # SystemStructure type and constructor (uses Transform)
include("transforms.jl")    # Heading calculations, reinit!/reposition! (uses SystemStructure)
include("utilities.jl")     # Helpers, validation, state management
