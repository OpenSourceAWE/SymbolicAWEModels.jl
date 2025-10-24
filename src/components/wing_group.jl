# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

"""
    @mtkmodel WingGroup

A deformable wing section with twist dynamics.

This component models a spanwise section of a wing that can deform through twisting
about its chord line. Multiple groups compose a complete flexible wing, with each
group having independent twist dynamics driven by:
- Bridle line forces (creating moments about the twist axis)
- Aerodynamic moments
- Structural restoring torque
- Damping

# Parameters
- `moment_frac = 0.4`: Chordwise position of twist axis (0=LE, 1=TE)
- `damping = 50.0`: Twist damping coefficient [Ns/rad]
- `inertia = 1.0`: Rotational inertia for twist [kg⋅m²]
- `max_twist = π/2`: Maximum twist angle [rad]
- `dynamics_type = 1`: 1=DYNAMIC, 2=QUASI_STATIC

# States
- `twist(t)`: Twist angle [rad]
- `twist_ω(t)`: Twist angular velocity [rad/s]

# Variables
- `free_twist(t)`: Unclamped twist angle [rad]
- `twist_α(t)`: Twist angular acceleration [rad/s²]
- `tether_moment(t)`: Moment from bridle forces [Nm]
- `aero_moment(t)`: Aerodynamic moment [Nm]
- `total_moment(t)`: Total moment on group [Nm]

# External Inputs (Parameters)
- `tether_moment_ext = 0.0`: External tether moment [Nm]
- `aero_moment_ext = 0.0`: External aerodynamic moment [Nm]

# Equations

## Twist Dynamics (DYNAMIC mode)
The twist angle evolves according to a damped rotational equation:

```math
\\begin{aligned}
I \\ddot{\\theta} &= M_{total} - c\\dot{\\theta} \\\\
M_{total} &= M_{tether} + M_{aero} \\\\
\\theta &= \\text{clamp}(\\theta_{free}, -\\theta_{max}, \\theta_{max})
\\end{aligned}
```

where:
- ``I`` is the rotational inertia of the wing section
- ``\\theta`` is the twist angle
- ``M_{tether}`` is the moment from bridle line forces
- ``M_{aero}`` is the aerodynamic restoring/destabilizing moment
- ``c`` is the damping coefficient
- Clamping prevents unrealistic deformations

## Quasi-Static Mode
For QUASI_STATIC, the twist angle equilibrates instantaneously:

```math
\\begin{aligned}
M_{total} &= 0 \\\\
\\dot{\\theta} &= 0
\\end{aligned}
```

The twist angle is determined algebraically by moment balance.

## Moment Contributions

### Tether Moment
Bridle line forces create moments about the twist axis:
```math
M_{tether} = \\sum_{i=1}^{n} r_i \\times (\\mathbf{F}_i \\cdot \\hat{\\mathbf{z}})
```

where ``r_i`` is the moment arm from the twist axis to attachment point ``i``.

### Aerodynamic Moment
Depends on angle of attack, twist, and flow conditions (computed by VSM).

# Usage Example
```julia
using ModelingToolkit, SymbolicAWEModels

@named group1 = WingGroup(
    moment_frac = 0.4,
    damping = 50.0,
    inertia = 2.0,
    dynamics_type = 1  # DYNAMIC
)

@named group2 = WingGroup(
    dynamics_type = 2  # QUASI_STATIC
)

# Set moments (would come from bridle forces and aero in full system)
eqs = [
    group1.tether_moment_ext ~ 5.0  # [Nm]
    group1.aero_moment_ext ~ -3.0   # [Nm]
]
```

# Notes
- The twist is clamped to prevent numerical instability
- Damping prevents high-frequency oscillations
- QUASI_STATIC mode is faster but less realistic for dynamic maneuvers
- In a complete system, this would be coupled to:
  - A `RigidWing` component (for the overall wing orientation)
  - Multiple `PointMass` components (bridle attachment points)
  - Aerodynamic calculation (VortexStepMethod)
"""
@mtkmodel WingGroup begin
    @parameters begin
        moment_frac = 0.4, [
            description = "Chordwise position of twist axis (0=LE, 1=TE)"
        ]
        damping = 50.0, [description = "Twist damping coefficient [Ns/rad]"]
        inertia = 1.0, [description = "Rotational inertia for twist [kg⋅m²]"]
        max_twist = pi / 2, [description = "Maximum twist angle [rad]"]
        dynamics_type = 1, [description = "1=DYNAMIC, 2=QUASI_STATIC"]
        # External moment inputs (to be connected to tethers/aero)
        tether_moment_ext = 0.0, [description = "External tether moment [Nm]"]
        aero_moment_ext = 0.0, [description = "External aerodynamic moment [Nm]"]
    end

    @variables begin
        # States
        twist(t), [description = "Twist angle [rad]"]
        twist_ω(t), [description = "Twist angular velocity [rad/s]"]

        # Derived quantities
        free_twist(t), [description = "Unclamped twist angle [rad]"]
        twist_α(t), [description = "Twist angular acceleration [rad/s²]"]

        # Moments
        tether_moment(t), [description = "Moment from bridle forces [Nm]"]
        aero_moment(t), [description = "Aerodynamic moment [Nm]"]
        total_moment(t), [description = "Total moment [Nm]"]
    end

    @equations begin
        # External moment interface
        tether_moment ~ tether_moment_ext
        aero_moment ~ aero_moment_ext

        # Total moment
        total_moment ~ tether_moment + aero_moment

        # Twist angle with clamping
        twist ~ clamp(free_twist, -max_twist, max_twist)

        # Dynamics based on type
        if dynamics_type == 1  # DYNAMIC
            # Angular acceleration
            twist_α ~ total_moment / inertia - damping * twist_ω

            # State evolution
            D(free_twist) ~ twist_ω
            D(twist_ω) ~ twist_α

        elseif dynamics_type == 2  # QUASI_STATIC
            # Moment balance (algebraic)
            total_moment ~ 0.0

            # Zero velocity and acceleration
            twist_ω ~ 0.0
            twist_α ~ 0.0

            # Free twist equals constrained twist (no clamping dynamics)
            free_twist ~ twist
        end
    end
end
