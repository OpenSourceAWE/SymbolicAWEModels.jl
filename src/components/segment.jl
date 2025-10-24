# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

"""
    @mtkmodel TetherSegment

A spring-damper element with aerodynamic drag, connecting two mechanical nodes.

**Note:** Named `TetherSegment` (not `Segment`) to avoid conflict with legacy `Segment` struct.

This component models a flexible tether segment with:
- Axial spring-damper mechanics (Hooke's law with damping)
- Asymmetric compression stiffness
- Aerodynamic drag perpendicular to the segment
- Wind shear effects based on height

# Parameters
- `l0 = 1.0`: Unstretched length [m]
- `axial_stiffness = 1e6`: Axial stiffness [N]
- `axial_damping = 1e3`: Axial damping [Ns]
- `diameter = 0.004`: Tether diameter [m]
- `compression_frac = 0.1`: Stiffness reduction factor in compression
- `cd_tether = 1.1`: Drag coefficient for the tether
- `rho_tether = 724.0`: Tether material density [kg/m³]
- `rho_air = 1.225`: Air density at ground level [kg/m³]
- `v_wind_gnd[1:3] = [10.0, 0, 0]`: Ground wind vector [m/s]
- `wind_shear_exp = 0.0`: Wind profile exponent (0=constant, 0.14=typical)

# Connectors
- `node1::MechanicalNode`: First connection point
- `node2::MechanicalNode`: Second connection point

# Variables
- `len(t)`: Current length of segment [m]
- `unit_vec(t)[1:3]`: Unit vector from node1 to node2
- `spring_force(t)`: Axial spring-damper force magnitude [N]
- `spring_force_vec(t)[1:3]`: Axial spring-damper force vector [N]
- `drag_force(t)[1:3]`: Aerodynamic drag force vector [N]
- `rel_vel(t)[1:3]`: Relative velocity between nodes [m/s]
- `spring_vel(t)`: Rate of extension [m/s]
- `stiffness(t)`: Effective axial stiffness [N/m]
- `damping(t)`: Effective axial damping [Ns/m]
- `segment_mass(t)`: Mass of the segment [kg]
- `height(t)`: Average height of segment for wind calculation [m]
- `segment_vel(t)[1:3]`: Average velocity of segment [m/s]
- `segment_rho(t)`: Air density at segment height [kg/m³]
- `wind_vel(t)[1:3]`: Wind velocity at segment height [m/s]
- `va(t)[1:3]`: Apparent wind vector [m/s]
- `area(t)`: Drag area of segment [m²]
- `app_perp_vel(t)[1:3]`: Perpendicular component of apparent wind [m/s]

# Equations

## Spring-Damper Force
The axial force follows Hooke's law with velocity damping:

```math
\\begin{aligned}
\\mathbf{F}_{spring} &= F_{axial} \\hat{\\mathbf{u}} \\\\
F_{axial} &= k_{eff}(l - l_0) - c_{eff}\\dot{l} \\\\
\\hat{\\mathbf{u}} &= \\frac{\\mathbf{r}_2 - \\mathbf{r}_1}{l} \\\\
\\dot{l} &= (\\mathbf{v}_1 - \\mathbf{v}_2) \\cdot \\hat{\\mathbf{u}}
\\end{aligned}
```

Where effective stiffness is:
```math
k_{eff} = \\begin{cases}
k / l & \\text{if } l > l_0 \\text{ (tension)} \\\\
\\alpha k / l & \\text{if } l \\leq l_0 \\text{ (compression)}
\\end{cases}
```

And effective damping:
```math
c_{eff} = \\frac{c}{l}
```

## Aerodynamic Drag
The drag force acts perpendicular to the segment:

```math
\\begin{aligned}
\\mathbf{F}_{drag} &= \\frac{1}{2}\\rho_{air}(h) C_d A |\\mathbf{v}_a| \\mathbf{v}_{a,\\perp} \\\\
\\mathbf{v}_a &= \\mathbf{v}_{wind}(h) - \\mathbf{v}_{segment} \\\\
\\mathbf{v}_{a,\\perp} &= \\mathbf{v}_a - (\\mathbf{v}_a \\cdot \\hat{\\mathbf{u}})\\hat{\\mathbf{u}} \\\\
A &= l \\cdot d \\\\
h &= \\frac{z_1 + z_2}{2}
\\end{aligned}
```

Wind profile (power law):
```math
\\mathbf{v}_{wind}(h) = \\mathbf{v}_{wind,gnd} \\left(\\frac{h}{h_{ref}}\\right)^\\alpha
```

## Force Application
Forces are distributed to the nodes:
```math
\\begin{aligned}
\\mathbf{F}_{node1} &= \\mathbf{F}_{spring} + 0.5\\mathbf{F}_{drag} + 0.5 m g \\hat{\\mathbf{z}} \\\\
\\mathbf{F}_{node2} &= -\\mathbf{F}_{spring} + 0.5\\mathbf{F}_{drag} + 0.5 m g \\hat{\\mathbf{z}}
\\end{aligned}
```

# Usage Example
```julia
using ModelingToolkit, SymbolicAWEModels

# Create a tether segment
@named seg = TetherSegment(
    l0 = 10.0,
    axial_stiffness = 1.2e6,
    axial_damping = 500.0,
    diameter = 0.004,
    compression_frac = 0.1,
    cd_tether = 1.1
)

# Connect to points
@named p1 = PointMass(mass=0.1)
@named p2 = PointMass(mass=0.1)

eqs = [
    connect(p1.node, seg.node1)
    connect(p2.node, seg.node2)
]
```

# Notes
- Half the segment mass and drag are applied to each node
- The compression stiffness factor prevents numerical instability when tethers slacken
- Wind shear can be disabled by setting `wind_shear_exp = 0`
- Drag force is applied to the average position between nodes
"""
@mtkmodel TetherSegment begin
    @parameters begin
        l0 = 1.0, [description = "Unstretched length [m]"]
        axial_stiffness = 1e6, [description = "Axial stiffness [N]"]
        axial_damping = 1e3, [description = "Axial damping [Ns]"]
        diameter = 0.004, [description = "Tether diameter [m]"]
        compression_frac = 0.1, [
            description = "Stiffness reduction factor in compression"
        ]
        cd_tether = 1.1, [description = "Drag coefficient"]
        rho_tether = 724.0, [description = "Tether material density [kg/m³]"]
        rho_air = 1.225, [description = "Air density at ground [kg/m³]"]
        v_wind_gnd[1:3] = [10.0, 0, 0], [description = "Ground wind vector [m/s]"]
        wind_shear_exp = 0.0, [
            description = "Wind profile exponent (0=constant, 0.14=typical)"
        ]
        g_earth = 9.81, [description = "Gravitational acceleration [m/s²]"]
        h_ref = 1.0, [description = "Reference height for wind profile [m]"]
    end

    @variables begin
        # Geometric quantities
        len(t), [description = "Current length [m]"]
        unit_vec(t)[1:3], [description = "Unit vector from node1 to node2"]
        segment_vec(t)[1:3], [description = "Position difference vector [m]"]

        # Spring-damper dynamics
        spring_force(t), [description = "Axial spring-damper force [N]"]
        spring_force_vec(t)[1:3], [description = "Axial force vector [N]"]
        rel_vel(t)[1:3], [description = "Relative velocity [m/s]"]
        spring_vel(t), [description = "Extension rate [m/s]"]
        stiffness(t), [description = "Effective stiffness [N/m]"]
        damping(t), [description = "Effective damping [Ns/m]"]

        # Mass properties
        segment_mass(t), [description = "Mass of segment [kg]"]

        # Aerodynamic quantities
        height(t), [description = "Average height [m]"]
        segment_vel(t)[1:3], [description = "Average velocity [m/s]"]
        segment_rho(t), [description = "Air density at height [kg/m³]"]
        wind_vel(t)[1:3], [description = "Wind velocity at height [m/s]"]
        va(t)[1:3], [description = "Apparent wind [m/s]"]
        area(t), [description = "Drag area [m²]"]
        app_perp_vel(t)[1:3], [description = "Perpendicular apparent wind [m/s]"]
        drag_force(t)[1:3], [description = "Drag force vector [N]"]

        # Wind factor
        wind_factor(t), [description = "Wind speed factor at height"]
    end

    @components begin
        node1 = MechanicalNode()
        node2 = MechanicalNode()
    end

    @equations begin
        # ==================== GEOMETRY ==================== #
        segment_vec ~ node2.pos - node1.pos
        len ~ max(1e-6, sqrt(sum(segment_vec .^ 2)))  # Prevent division by zero
        unit_vec ~ segment_vec / len

        # ==================== SPRING-DAMPER FORCE ==================== #
        rel_vel ~ node1.vel - node2.vel
        spring_vel ~ sum(rel_vel .* unit_vec)

        # Effective properties (length-dependent)
        damping ~ axial_damping / len

        # Asymmetric stiffness (softer in compression)
        stiffness ~ ifelse(
            len > l0,
            axial_stiffness / len,  # Tension
            compression_frac * axial_stiffness / len   # Compression
        )

        # Axial force (Hooke's law + damping)
        spring_force ~ stiffness * (len - l0) - damping * spring_vel
        spring_force_vec ~ spring_force * unit_vec

        # ==================== MASS ==================== #
        segment_mass ~ rho_tether * pi * (diameter / 2)^2 * l0

        # ==================== AERODYNAMIC DRAG ==================== #
        height ~ max(0.0, 0.5 * (node1.pos[3] + node2.pos[3]))
        segment_vel ~ 0.5 * (node1.vel + node2.vel)

        # Wind profile (power law or constant)
        wind_factor ~ ifelse(
            wind_shear_exp == 0.0,
            1.0,
            (max(height, h_ref) / h_ref)^wind_shear_exp
        )
        wind_vel ~ wind_factor * v_wind_gnd

        # Air density (simplified - could use AtmosphericModels)
        segment_rho ~ rho_air

        # Apparent wind and drag
        va ~ wind_vel - segment_vel
        area ~ len * diameter
        app_perp_vel ~ va - (sum(va .* unit_vec)) * unit_vec
        drag_force ~ (0.5 * segment_rho * cd_tether * sqrt(sum(va .^ 2)) * area) * app_perp_vel

        # ==================== FORCE APPLICATION ==================== #
        # Spring force: +F on node2, -F on node1
        # Drag: half to each node
        # Gravity: half segment mass to each node
        node1.force ~ -spring_force_vec + 0.5 * drag_force + [0, 0, -0.5 * segment_mass * g_earth]
        node2.force ~ spring_force_vec + 0.5 * drag_force + [0, 0, -0.5 * segment_mass * g_earth]
    end
end
