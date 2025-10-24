# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

"""
    @mtkmodel PulleyComponent

A pulley component that redistributes length between two tether segments.

**Note:** Named `PulleyComponent` (not `Pulley`) to avoid conflict with legacy `Pulley` struct.

This component enforces a constant total length constraint while allowing dynamic
redistribution of length between two connected segments. It models the inertia
of the redistributed tether mass.

# Parameters
- `sum_len = 10.0`: Total constant length (l1 + l2) [m]
- `dynamics_type = 1`: 1=DYNAMIC (with inertia), 2=QUASI_STATIC (instant balance)
- `damping = 5.0`: Damping coefficient for length redistribution [1/s]
- `rho_tether = 724.0`: Tether material density [kg/m³]
- `diameter = 0.004`: Tether diameter [m]

# States
- `len(t)`: Length of first segment [m]
- `vel(t)`: Rate of length change [m/s]

# Variables
- `len2(t)`: Length of second segment [m]
- `force(t)`: Net force causing length redistribution [N]
- `mass(t)`: Total mass of both segments [kg]
- `acc(t)`: Acceleration of length change [m/s²]

# Connectors
- `seg1_force::RealInput`: Force in first segment (from external binding)
- `seg2_force::RealInput`: Force in second segment (from external binding)

# Equations

## Length Constraint
```math
l_1 + l_2 = L_{total} = \\text{const}
```

## Force Balance and Dynamics
For DYNAMIC type:
```math
\\begin{aligned}
F_{pulley} &= F_1 - F_2 \\\\
m \\ddot{l}_1 &= F_{pulley} - c\\dot{l}_1 \\\\
m &= \\rho_{tether} \\pi (d/2)^2 L_{total}
\\end{aligned}
```

For QUASI_STATIC type:
```math
\\begin{aligned}
F_1 &= F_2 \\\\
\\dot{l}_1 &= 0 \\\\
\\ddot{l}_1 &= 0
\\end{aligned}
```

# Usage Example
```julia
using ModelingToolkit, SymbolicAWEModels

# Create pulley
@named pulley = PulleyComponent(sum_len=20.0, dynamics_type=1, damping=5.0)

# Note: This component typically requires manual force binding from segments
# In the system builder, you would connect segment forces like:
# pulley_eqs = [
#     pulley.seg1_force ~ seg1.spring_force
#     pulley.seg2_force ~ seg2.spring_force
# ]
```

# Notes
- The pulley does not have MechanicalNode connectors because it doesn't directly
  connect to points - it manages segment lengths
- Forces must be provided via parameters or external equations
- Total mass is distributed across both segments
- This is a simplified model - real pulleys have friction and rotational inertia
"""
@mtkmodel PulleyComponent begin
    @parameters begin
        sum_len = 10.0, [description = "Total constant length [m]"]
        dynamics_type = 1, [description = "1=DYNAMIC, 2=QUASI_STATIC"]
        damping = 5.0, [description = "Damping coefficient [1/s]"]
        rho_tether = 724.0, [description = "Tether density [kg/m³]"]
        diameter = 0.004, [description = "Tether diameter [m]"]
        seg1_force = 0.0, [description = "Force in first segment [N] (external input)"]
        seg2_force = 0.0, [description = "Force in second segment [N] (external input)"]
    end

    @variables begin
        len(t), [description = "Length of first segment [m]"]
        vel(t), [description = "Velocity of length change [m/s]"]
        len2(t), [description = "Length of second segment [m]"]
        force(t), [description = "Net force on pulley [N]"]
        mass(t), [description = "Total mass of segments [kg]"]
        acc(t), [description = "Acceleration of length change [m/s²]"]
    end

    @equations begin
        # Length constraint
        len2 ~ sum_len - len

        # Mass calculation
        mass ~ rho_tether * pi * (diameter / 2)^2 * sum_len

        # Force balance
        force ~ seg1_force - seg2_force

        # Dynamics based on type
        if dynamics_type == 1  # DYNAMIC
            D(len) ~ vel
            D(vel) ~ acc
            acc ~ force / mass - damping * vel

        elseif dynamics_type == 2  # QUASI_STATIC
            vel ~ 0.0
            acc ~ 0.0
            # Algebraic constraint: forces must balance
            force ~ 0.0
        end
    end
end
