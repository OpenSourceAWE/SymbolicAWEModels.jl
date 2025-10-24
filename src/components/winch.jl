# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

"""
    @mtkmodel WinchComponent

A winch component with motor, drum, gearbox, and friction dynamics.

**Note:** Named `WinchComponent` (not `Winch`) to avoid conflict with legacy `Winch` struct.

This component models a ground-based winch that controls tether length through
torque or speed regulation. It includes:
- Motor torque input
- Gear reduction
- Drum radius conversion (rotation → linear motion)
- Coulomb and viscous friction
- Rotational inertia
- Optional brake

# Parameters
- `gear_ratio = 110.0`: Gear reduction ratio
- `drum_radius = 0.1632`: Drum radius [m]
- `f_coulomb = 133.0`: Coulomb friction torque [Nm]
- `c_vf = 7.0`: Viscous friction coefficient [Ns/m]
- `inertia_total = 1.47`: Total rotational inertia [kg⋅m²]
- `brake = false`: Brake engaged (locks tether length)

# States
- `tether_len(t)`: Tether length [m]
- `tether_vel(t)`: Tether velocity (reel-out rate) [m/s]

# Variables
- `tether_acc(t)`: Tether acceleration [m/s²]
- `set_value(t)`: Control input (motor torque) [Nm]
- `tether_force(t)`: Force in the tether [N] (external input)
- `ω_motor(t)`: Motor angular velocity [rad/s]
- `τ_friction(t)`: Friction torque [Nm]
- `τ_motor(t)`: Motor torque [Nm]
- `τ_total(t)`: Total torque on drum [Nm]
- `α_motor(t)`: Motor angular acceleration [rad/s²]

# Connectors
- `ctrl::ControlSignal`: Control input connector

# Equations

## Kinematic Conversion
```math
\\begin{aligned}
\\omega_{motor} &= \\frac{n}{r_{drum}} v_{tether} \\\\
\\alpha_{motor} &= \\frac{n}{r_{drum}} a_{tether}
\\end{aligned}
```

## Friction Model
```math
\\begin{aligned}
\\tau_{friction} &= \\tau_{Coulomb} + \\tau_{viscous} \\\\
\\tau_{Coulomb} &= \\text{sign}(\\omega_{motor}) \\frac{F_c r_{drum}}{n} \\\\
\\tau_{viscous} &= c_{vf} \\omega_{motor} \\frac{r_{drum}^2}{n^2}
\\end{aligned}
```

Where a smooth sign function is used to avoid discontinuity at zero:
```math
\\text{sign}(x) \\approx \\frac{x}{\\sqrt{x^2 + \\epsilon^2}}
```

## Torque Balance
```math
\\begin{aligned}
\\tau_{total} &= \\tau_{motor} + \\frac{r_{drum}}{n} F_{tether} - \\tau_{friction} \\\\
I_{total} \\alpha_{motor} &= \\tau_{total}
\\end{aligned}
```

## Tether Dynamics
```math
\\begin{aligned}
\\frac{dl}{dt} &= v_{tether} \\\\
\\frac{dv}{dt} &= a_{tether} = \\frac{r_{drum}}{n} \\alpha_{motor}
\\end{aligned}
```

With brake engaged:
```math
\\frac{dl}{dt} = 0, \\quad \\frac{dv}{dt} = 0
```

# Usage Example
```julia
using ModelingToolkit, SymbolicAWEModels

# Create winch with control input
@named winch = WinchComponent(
    gear_ratio = 110.0,
    drum_radius = 0.1632,
    inertia_total = 1.47
)

# Set control torque
eqs = [
    winch.ctrl.value ~ 50.0  # 50 Nm motor torque
    winch.tether_force ~ 1000.0  # External: tether tension
]
```

# Notes
- The `tether_force` parameter must be connected to the actual tether tension
- Positive velocity = reel out, negative = reel in
- The smooth sign function prevents solver issues near zero velocity
- Brake overrides all dynamics and locks the tether length
"""
@mtkmodel WinchComponent begin
    @parameters begin
        gear_ratio = 110.0, [description = "Gear reduction ratio"]
        drum_radius = 0.1632, [description = "Drum radius [m]"]
        f_coulomb = 133.0, [description = "Coulomb friction torque [Nm]"]
        c_vf = 7.0, [description = "Viscous friction coefficient [Ns/m]"]
        inertia_total = 1.47, [description = "Total rotational inertia [kg⋅m²]"]
        brake = false, [description = "Brake engaged"]
        tether_force = 0.0, [description = "Tether tension [N] (external input)"]
    end

    @variables begin
        tether_len(t), [description = "Tether length [m]"]
        tether_vel(t), [description = "Tether velocity [m/s]"]
        tether_acc(t), [description = "Tether acceleration [m/s²]"]
        set_value(t), [description = "Control input (motor torque) [Nm]"]
        ω_motor(t), [description = "Motor angular velocity [rad/s]"]
        τ_friction(t), [description = "Friction torque [Nm]"]
        τ_motor(t), [description = "Motor torque [Nm]"]
        τ_total(t), [description = "Total torque [Nm]"]
        α_motor(t), [description = "Motor angular acceleration [rad/s²]"]
        smooth_sign_val(t), [description = "Smooth sign of motor velocity"]
    end

    @components begin
        ctrl = ControlSignal()
    end

    @equations begin
        # Control input
        set_value ~ ctrl.value

        # Kinematic conversion
        ω_motor ~ gear_ratio / drum_radius * tether_vel

        # Smooth sign function (EPSILON = 6 for smoothness)
        smooth_sign_val ~ ω_motor / sqrt(ω_motor^2 + 36.0)

        # Friction model
        τ_friction ~ (
            smooth_sign_val * f_coulomb * drum_radius / gear_ratio +
            c_vf * ω_motor * drum_radius^2 / gear_ratio^2
        )

        # Motor and tether force torque
        τ_motor ~ set_value

        # Total torque balance
        τ_total ~ τ_motor + drum_radius / gear_ratio * tether_force - τ_friction

        # Rotational dynamics
        α_motor ~ τ_total / inertia_total

        # Tether dynamics
        tether_acc ~ drum_radius / gear_ratio * α_motor

        # State equations (with brake override)
        if brake
            D(tether_len) ~ 0.0
            D(tether_vel) ~ 0.0
        else
            D(tether_len) ~ tether_vel
            D(tether_vel) ~ tether_acc
        end
    end
end
