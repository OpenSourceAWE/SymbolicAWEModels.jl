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
- `fixed_pos = [0.0, 0.0, 0.0]`: Fixed position for STATIC type [m]
- `fix_sphere = false`: If true, constrains motion to a sphere (radial only)
- `bridle_damping = 0.0`: Damping coefficient for bridle points [Ns/m]

# Connectors
- `node::MechanicalNode`: Mechanical connection point

# Variables
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

@named ground = PointMass(dynamics_type=4, fixed_pos=[0, 0, 0])
@named kite_point = PointMass(mass=0.5, dynamics_type=1, bridle_damping=2.0)
@named quasi_point = PointMass(dynamics_type=2)  # Force-balanced point

# Connect to other components
@named seg = TetherSegment()
eqs = [
    connect(ground.node, seg.node1)
    connect(kite_point.node, seg.node2)
]
```

# Notes
- The `node.force` connector accumulates forces from all connected components
- Gravity must be added externally via `disturb_force` or through a separate component
- For QUASI_STATIC, the solver finds the position where all forces balance
"""
@mtkmodel PointMass begin
    @structural_parameters begin
        fixed_pos = [0.0, 0.0, 0.0]
        wing_vel = [0.0, 0.0, 0.0]
    end

    @parameters begin
        mass = 1.0, [description = "Mass of the point [kg]"]
        dynamics_type = 1, [
            description = "Dynamics type: 1=DYNAMIC, 2=QUASI_STATIC, 3=WING, 4=STATIC"
        ]
        fix_sphere = false, [description = "Constrain to sphere (radial motion only)"]
        bridle_damping = 0.0, [description = "Damping coefficient for bridle points [Ns/m]"]
        g_earth = 9.81, [description = "Gravitational acceleration [m/s²]"]
    end

    @variables begin
        acc(t)[1:3], [description = "Acceleration vector [m/s²]"]
        disturb_force(t)[1:3] = zeros(3), [description = "External disturbance force [N]"]
        net_force(t)[1:3], [description = "Net force on the point [N]"]
        axis(t)[1:3], [description = "Normalized radial direction for sphere constraint"]
        bridle_damp_vec(t)[1:3], [
            description = "Bridle damping force vector [N]"
        ]
    end

    @components begin
        node = MechanicalNode()
    end

    @equations begin
        # Compute normalized radial axis for sphere constraint
        axis ~ node.pos / max(1e-6, sqrt(sum(node.pos .^ 2)))

        # Bridle damping (velocity relative to wing)
        bridle_damp_vec ~ bridle_damping * (node.vel - wing_vel)

        # Net force calculation (connector force + gravity + disturbances)
        net_force ~ node.force + [0, 0, -mass * g_earth] + disturb_force

        # Equations based on dynamics type
        if dynamics_type == 4  # STATIC
            D.(node.pos) .~ zeros(3)
            D.(node.vel) .~ zeros(3)
            node.pos ~ fixed_pos
            node.vel ~ zeros(3)
            acc ~ zeros(3)

        elseif dynamics_type == 1  # DYNAMIC
            # Apply sphere constraint if requested
            if fix_sphere
                D.(node.pos) .~ (sum(node.vel .* axis)) * axis
                D.(node.vel) .~ (sum(acc .* axis)) * axis
            else
                D.(node.pos) .~ node.vel
                D.(node.vel) .~ acc
            end
            # Newton's second law with bridle damping
            acc ~ net_force / mass - bridle_damp_vec

        elseif dynamics_type == 2  # QUASI_STATIC
            node.vel ~ zeros(3)
            acc ~ zeros(3)
            # Force balance constraint (algebraic)
            net_force .~ zeros(3)

        elseif dynamics_type == 3  # WING (placeholder - actual position set externally)
            # Velocity and acceleration determined by wing dynamics
            node.vel ~ zeros(3)
            acc ~ zeros(3)
            # Position is set via external equation (wing transformation)
            # D.(node.pos) .~ zeros(3)  # Will be overridden by connect
        end
    end
end
