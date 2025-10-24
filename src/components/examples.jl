# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

"""
# Component-Based Modeling Examples

This file contains example systems built using the modern MTK component approach.
These examples demonstrate how to compose complex AWE systems from modular components.

## Example Systems

1. `simple_pendulum_system`: Single point mass on a tether
2. `double_tether_system`: Two segments with dynamic pulley
3. `winch_controlled_system`: Ground winch controlling tether length

Each example can be used as a template for building custom systems.
"""

"""
    simple_pendulum_system()

Create a simple pendulum: a kite point mass suspended from a fixed ground point
by a single elastic tether segment.

# System Topology
```
Ground (fixed)
   |
   | Tether (elastic)
   |
Kite (dynamic mass)
```

# Returns
- `sys::ODESystem`: The complete system ready for `structural_simplify`

# Usage
```julia
using ModelingToolkit, OrdinaryDiffEq

sys = simple_pendulum_system()
simple_sys = structural_simplify(sys)
prob = ODEProblem(simple_sys, [], (0.0, 10.0))
sol = solve(prob, Tsit5())
```
"""
function simple_pendulum_system()
    t = ModelingToolkit.t_nounits

    # Components
    @named ground = PointMass(
        dynamics_type = 4,  # STATIC
        fixed_pos = [0.0, 0.0, 0.0],
        mass = 1.0
    )

    @named kite = PointMass(
        dynamics_type = 1,  # DYNAMIC
        mass = 5.0,
        fix_sphere = true,  # Constrain to sphere (like a physical pendulum)
        bridle_damping = 0.0
    )

    @named tether = TetherSegment(
        l0 = 100.0,
        axial_stiffness = 1.2e6,
        axial_damping = 500.0,
        diameter = 0.004,
        compression_frac = 0.1,
        cd_tether = 1.1,
        v_wind_gnd = [10.0, 0.0, 0.0],
        wind_shear_exp = 0.14
    )

    # Connections
    eqs = [
        connect(ground.node, tether.node1)
        connect(tether.node2, kite.node)
    ]

    # Build system
    @named sys = ODESystem(eqs, t;
        systems = [ground, kite, tether]
    )

    return sys
end

"""
    double_tether_system()

Create a system with two tether segments connected by a dynamic pulley,
suspended from a fixed point.

# System Topology
```
Ground (fixed)
   |
   | Segment 1
   |
Pulley (dynamic redistribution)
   |
   | Segment 2
   |
Kite (dynamic mass)
```

The pulley allows length to redistribute between the two segments while
maintaining constant total length.

# Returns
- `sys::ODESystem`: The complete system

# Usage
```julia
sys = double_tether_system()
simple_sys = structural_simplify(sys)
# ... solve as usual
```
"""
function double_tether_system()
    t = ModelingToolkit.t_nounits

    # Points
    @named ground = PointMass(dynamics_type=4, fixed_pos=[0,0,0])
    @named pulley_point = PointMass(dynamics_type=1, mass=0.5)
    @named kite = PointMass(dynamics_type=1, mass=5.0)

    # Segments
    @named seg1 = TetherSegment(
        l0 = 50.0,
        axial_stiffness = 1.2e6,
        diameter = 0.004
    )
    @named seg2 = TetherSegment(
        l0 = 50.0,
        axial_stiffness = 1.2e6,
        diameter = 0.004
    )

    # Pulley dynamics
    @named pulley = PulleyComponent(
        sum_len = 100.0,
        dynamics_type = 1,  # DYNAMIC
        damping = 5.0
    )

    # Connections
    eqs = [
        # Mechanical connections
        connect(ground.node, seg1.node1)
        connect(seg1.node2, pulley_point.node)
        connect(pulley_point.node, seg2.node1)
        connect(seg2.node2, kite.node)

        # Pulley force coupling
        pulley.seg1_force ~ seg1.spring_force
        pulley.seg2_force ~ seg2.spring_force

        # Pulley length control (this would need integration with seg1/seg2.l0)
        # For now, this is a simplified coupling - full implementation would
        # require more sophisticated length redistribution
    ]

    @named sys = ODESystem(eqs, t;
        systems = [ground, pulley_point, kite, seg1, seg2, pulley]
    )

    return sys
end

"""
    winch_controlled_system()

Create a system with a ground-based winch controlling tether length.

# System Topology
```
Winch (motor control)
   |
   | Tether (variable length)
   |
Kite (dynamic mass)
```

The winch applies motor torque to reel the tether in or out, with friction
and inertia effects.

# Returns
- `sys::ODESystem`: The complete system

# Parameters
- Motor torque can be set via `winch.ctrl.value`

# Usage
```julia
sys = winch_controlled_system()
simple_sys = structural_simplify(sys)

# Set initial conditions and parameters
# u0 = [...] # Initial states
# p = [...]  # Parameters (can set motor torque here)
# prob = ODEProblem(simple_sys, u0, (0.0, 100.0), p)
# sol = solve(prob, Rosenbrock23())  # Stiff solver for winch dynamics
```
"""
function winch_controlled_system()
    t = ModelingToolkit.t_nounits

    # Components
    @named winch = WinchComponent(
        gear_ratio = 110.0,
        drum_radius = 0.1632,
        f_coulomb = 133.0,
        c_vf = 7.0,
        inertia_total = 1.47
    )

    @named ground = PointMass(dynamics_type=4, fixed_pos=[0,0,0])
    @named kite = PointMass(dynamics_type=1, mass=5.0)

    @named tether = TetherSegment(
        l0 = 100.0,  # This will be controlled by winch
        axial_stiffness = 1.2e6,
        axial_damping = 500.0,
        diameter = 0.004
    )

    # Connections
    eqs = [
        # Mechanical
        connect(ground.node, tether.node1)
        connect(tether.node2, kite.node)

        # Winch control
        winch.ctrl.value ~ 50.0  # 50 Nm motor torque (constant for now)
        winch.tether_force ~ tether.spring_force

        # Length coupling (simplified - full version would modify tether.l0 dynamically)
        # tether.l0 ~ winch.tether_len  # This coupling needs careful implementation
    ]

    @named sys = ODESystem(eqs, t;
        systems = [winch, ground, kite, tether]
    )

    return sys
end

"""
    build_custom_system(builder::Function)

Helper function for users to build custom systems with a clean syntax.

# Usage
```julia
sys = build_custom_system() do
    @named my_ground = PointMass(dynamics_type=4, fixed_pos=[0,0,0])
    @named my_kite = PointMass(dynamics_type=1, mass=10.0)
    @named my_tether = TetherSegment(l0=200.0, ...)

    eqs = [
        connect(my_ground.node, my_tether.node1)
        connect(my_tether.node2, my_kite.node)
    ]

    return (eqs, [my_ground, my_kite, my_tether])
end
```
"""
function build_custom_system(builder::Function)
    t = ModelingToolkit.t_nounits

    eqs, components = builder()

    @named sys = ODESystem(eqs, t; systems = components)
    return sys
end

# Export example functions
export simple_pendulum_system, double_tether_system, winch_controlled_system, build_custom_system
