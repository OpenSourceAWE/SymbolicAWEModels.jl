# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

"""
    @mtkmodel PointMass

A 3D point mass with configurable dynamics type.

**Note:** Named `PointMass` (not `MassPoint`) to avoid conflict with legacy `Point` struct.

This component models a particle in 3D space that can be:
- **DYNAMIC**: Moves according to Newton's second law (``\\mathbf{F} = m\\mathbf{a}``)
- **STATIC**: Position is fixed in space
- **QUASI_STATIC**: Forces are balanced instantaneously (``\\mathbf{F} = 0``)
- **WING**: Attached to a wing rigid body (position determined externally)

# Parameters
- `mass = 1.0`: Mass of the point [kg]
- `dynamics_type = 1`: Type of dynamics (1=DYNAMIC, 2=QUASI_STATIC, 3=WING, 4=STATIC)
- `fixed_pos[1:3] = [0.0, 0.0, 0.0]`: Fixed position for STATIC type [m]
- `fix_sphere = false`: If true, constrains motion to a sphere (radial only)
- `bridle_damping = 0.0`: Damping coefficient for bridle points [Ns/m]

# Variables (Connection Interface)
- `pos(t)[1:3]`: Position in world frame [m]
- `vel(t)[1:3]`: Velocity in world frame [m/s]
- `force(t)[1:3]`: Force on the point [N] (sum of all external forces)
- `wing_vel(t)[1:3]`: Wing velocity for bridle damping [m/s] (connect to zeros if no wing)

# Variables (Internal)
- `acc(t)[1:3]`: Acceleration vector [m/s²]
- `disturb_force(t)[1:3]`: External disturbance force [N]
- `net_force(t)[1:3]`: Net force on the point [N]

# Equations
The component implements different equation sets based on `dynamics_type`:

## DYNAMIC (type = 1)
```math
\\begin{aligned}
\\frac{d\\mathbf{r}}{dt} &= \\mathbf{v} \\\\
\\frac{d\\mathbf{v}}{dt} &= \\mathbf{a} \\\\
\\mathbf{a} &= \\frac{\\mathbf{F}_{net}}{m} \\\\
\\mathbf{F}_{net} &= \\mathbf{F}_{external} + \\mathbf{F}_{disturb}
\\end{aligned}
```

With optional sphere constraint:
```math
\\frac{d\\mathbf{r}}{dt} = \\left(\\mathbf{v} \\cdot \\hat{\\mathbf{r}}\\right) \\hat{\\mathbf{r}}
```

## STATIC (type = 4)
```math
\\begin{aligned}
\\mathbf{r} &= \\mathbf{r}_{fixed} \\\\
\\mathbf{v} &= \\mathbf{0} \\\\
\\mathbf{a} &= \\mathbf{0}
\\end{aligned}
```

## QUASI_STATIC (type = 2)
```math
\\begin{aligned}
\\mathbf{v} &= \\mathbf{0} \\\\
\\mathbf{a} &= \\mathbf{0} \\\\
\\mathbf{F}_{net} &= \\mathbf{0}
\\end{aligned}
```

Position is determined algebraically by force balance.

## WING (type = 3)
```math
\\begin{aligned}
\\mathbf{r} &= \\mathbf{r}_{wing\\_COM} + \\mathbf{R}_{b \\to w} \\mathbf{r}_b \\\\
\\mathbf{v} &= \\mathbf{0} \\quad \\text{(determined by wing)} \\\\
\\mathbf{a} &= \\mathbf{0} \\quad \\text{(determined by wing)}
\\end{aligned}
```

# Usage Example
```julia
using ModelingToolkit, SymbolicAWEModels

@named ground = PointMass(dynamics_type=4)
@named kite_point = PointMass(mass=0.5, dynamics_type=1, bridle_damping=2.0)
@named quasi_point = PointMass(dynamics_type=2)  # Force-balanced point

# Connect to other components via direct variable equations
@named seg = TetherSegment()
eqs = [
    # Connect ground to segment end 1
    seg.pos1 ~ ground.pos
    seg.vel1 ~ ground.vel
    ground.force ~ -seg.force1  # Newton's 3rd law
    ground.wing_vel ~ zeros(3)  # No wing attached

    # Connect kite to segment end 2
    seg.pos2 ~ kite_point.pos
    seg.vel2 ~ kite_point.vel
    kite_point.force ~ -seg.force2
    kite_point.wing_vel ~ zeros(3)  # No wing attached
]
```

# Notes
- Forces are connected directly between components (no automatic accumulation)
- Gravity is included automatically via `g_earth` parameter
- For QUASI_STATIC, the solver finds the position where all forces balance
- Use Newton's 3rd law when connecting: component A's force on B = -(component B's force on A)
"""
@mtkmodel PointMass begin
    @parameters begin
        mass = 1.0, [description = "Mass of the point [kg]"]
        dynamics_type = 1, [
            description = "Dynamics type: 1=DYNAMIC, 2=QUASI_STATIC, 3=WING, 4=STATIC"
        ]
        fix_sphere = false, [description = "Constrain to sphere (radial motion only)"]
        bridle_damping = 0.0, [description = "Damping coefficient for bridle points [Ns/m]"]
        g_earth = 9.81, [description = "Gravitational acceleration [m/s²]"]
        fixed_pos[1:3] = [0.0, 0.0, 0.0], [description = "Fixed position for STATIC type [m]"]
    end

    @variables begin
        # Connection interface (exposed to other components)
        pos(t)[1:3], [description = "Position in world frame [m]"]
        vel(t)[1:3], [description = "Velocity in world frame [m/s]"]
        force(t)[1:3], [description = "Force on the point [N]"]
        wing_vel(t)[1:3], [description = "Wing velocity for bridle damping (connect to zeros if no wing) [m/s]"]

        # Internal variables
        acc(t)[1:3], [description = "Acceleration vector [m/s²]"]
        disturb_force(t)[1:3] = zeros(3), [description = "External disturbance force [N]"]
        net_force(t)[1:3], [description = "Net force on the point [N]"]
        axis(t)[1:3], [description = "Normalized radial direction for sphere constraint"]
        bridle_damp_vec(t)[1:3], [
            description = "Bridle damping force vector [N]"
        ]
    end

    @equations begin
        # ===== Always-present equations =====
        # Compute normalized radial axis for sphere constraint
        axis ~ pos / max(1e-6, sqrt(sum(pos .^ 2)))

        # Bridle damping (velocity relative to wing)
        bridle_damp_vec ~ bridle_damping * (vel - wing_vel)

        # Net force calculation (external force + gravity + disturbances)
        net_force ~ force + [0, 0, -mass * g_earth] + disturb_force
    end

    # ===== Conditional equations for different dynamics types =====
    @equations begin
        # Position derivative (depends on dynamics_type and fix_sphere)
        # STATIC (4): no motion
        # DYNAMIC (1): either sphere-constrained or unconstrained
        # QUASI_STATIC (2): position is algebraic (no D(pos) equation)
        # WING (3): position set externally (no D(pos) equation)

        D.(pos) .~ ifelse(
            dynamics_type == 4,  # STATIC
            zeros(3),
            ifelse(
                dynamics_type == 1,  # DYNAMIC
                ifelse(
                    fix_sphere,
                    (sum(vel .* axis)) * axis,  # Sphere constraint
                    vel  # Unconstrained
                ),
                zeros(3)  # QUASI_STATIC or WING: handled algebraically
            )
        )
    end

    @equations begin
        # Velocity derivative (depends on dynamics_type and fix_sphere)
        D.(vel) .~ ifelse(
            dynamics_type == 4,  # STATIC
            zeros(3),
            ifelse(
                dynamics_type == 1,  # DYNAMIC
                ifelse(
                    fix_sphere,
                    (sum(acc .* axis)) * axis,  # Sphere constraint
                    acc  # Unconstrained
                ),
                zeros(3)  # QUASI_STATIC or WING
            )
        )
    end

    @equations begin
        # Position constraints for STATIC type
        pos ~ ifelse(
            dynamics_type == 4,  # STATIC
            fixed_pos,
            pos  # Otherwise position evolves from ODE or algebraically
        )
    end

    @equations begin
        # Velocity constraints
        vel ~ ifelse(
            dynamics_type == 4,  # STATIC
            zeros(3),
            ifelse(
                (dynamics_type == 2) | (dynamics_type == 3),  # QUASI_STATIC or WING
                zeros(3),
                vel  # DYNAMIC: velocity evolves from ODE
            )
        )
    end

    @equations begin
        # Acceleration equation
        acc ~ ifelse(
            dynamics_type == 1,  # DYNAMIC
            net_force / mass - bridle_damp_vec,  # Newton's 2nd law with damping
            zeros(3)  # STATIC, QUASI_STATIC, or WING
        )
    end

    @equations begin
        # Force balance constraint for QUASI_STATIC
        net_force .~ ifelse(
            dynamics_type == 2,  # QUASI_STATIC
            zeros(3),  # Algebraic constraint: forces must balance
            net_force  # Otherwise net_force is just a computed variable
        )
    end
end
