# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

"""
    @connector MechanicalNode

A mechanical connector representing a 3D point with position, velocity, and force.

This connector follows Kirchhoff's flow law: the sum of all forces at a connection
point must equal zero. This is enforced by marking `force` as a `Flow` variable.

# Variables
- `pos(t)[1:3]`: Position vector in world frame [m]
- `vel(t)[1:3]`: Velocity vector in world frame [m/s]
- `force(t)[1:3]`: Force vector in world frame [N], marked as Flow variable

# Usage
```julia
@named node = MechanicalNode()
```

# Notes
When multiple components are connected via `connect(node1, node2)`, ModelingToolkit
automatically enforces:
- Position continuity: `node1.pos ~ node2.pos`
- Velocity continuity: `node1.vel ~ node2.vel`
- Force balance: `node1.force + node2.force ~ 0` (Flow variable property)
"""
@connector MechanicalNode begin
    @variables begin
        pos(t)[1:3], [description = "Position in world frame [m]"]
        vel(t)[1:3], [description = "Velocity in world frame [m/s]"]
        force(t)[1:3], [connect = Flow, description = "Force in world frame [N]"]
    end
end

"""
    @connector ControlSignal

A scalar control signal connector for actuator inputs.

# Variables
- `value(t)`: Control signal value (e.g., torque [Nm], force [N])

# Usage
```julia
@named ctrl = ControlSignal()
```
"""
@connector ControlSignal begin
    @variables begin
        value(t), [description = "Control signal value"]
    end
end
