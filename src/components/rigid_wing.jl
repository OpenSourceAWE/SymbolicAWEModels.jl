# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

"""
    @mtkmodel RigidWing

A rigid body wing with 6-DOF dynamics using quaternion orientation.

This component models a wing as a rigid body with:
- 3-DOF translational motion (Newton's 2nd law)
- 3-DOF rotational motion (Euler equations)
- Quaternion-based orientation (avoids singularities)
- Multiple attachment points for bridle connections
- Aerodynamic force/moment interfaces

# Parameters
- `mass = 25.0`: Wing mass [kg]
- `I_xx = 10.0`: Moment of inertia about x-axis [kg⋅m²]
- `I_yy = 15.0`: Moment of inertia about y-axis [kg⋅m²]
- `I_zz = 20.0`: Moment of inertia about z-axis [kg⋅m²]
- `y_damping = 150.0`: Lateral (y-axis) angular damping [Ns]
- `z_disturb = 0.0`: Disturbance torque about z-axis [Nm]
- `fix_sphere = false`: Constrain to sphere (radial motion only)
- `g_earth = 9.81`: Gravitational acceleration [m/s²]

# States
- `pos(t)[1:3]`: Position of center of mass in world frame [m]
- `vel(t)[1:3]`: Velocity in world frame [m/s]
- `Q_b_w(t)[1:4]`: Quaternion (body to world), [w,x,y,z] format
- `ω_b(t)[1:3]`: Angular velocity in body frame [rad/s]

# Variables
- `acc(t)[1:3]`: Linear acceleration [m/s²]
- `α_b(t)[1:3]`: Angular acceleration in body frame [rad/s²]
- `R_b_w(t)[1:3, 1:3]`: Rotation matrix (body to world)
- `Q_vel(t)[1:4]`: Quaternion derivative
- `aero_force_b(t)[1:3]`: Aerodynamic force in body frame [N]
- `aero_moment_b(t)[1:3]`: Aerodynamic moment in body frame [Nm]
- `tether_force_w(t)[1:3]`: Total tether force in world frame [N]
- `tether_moment_w(t)[1:3]`: Total tether moment in world frame [Nm]

# Connectors
Currently no connectors - attachment points to be added in future version.
Forces/moments are applied via parameters for now.

# Equations

## Quaternion Kinematics
The quaternion derivative is computed using the skew-symmetric matrix formulation:

```math
\\frac{dQ}{dt} = \\frac{1}{2}\\Omega(\\boldsymbol{\\omega})Q
```

where

```math
\\Omega(\\boldsymbol{\\omega}) = \\begin{bmatrix}
0 & -\\omega_x & -\\omega_y & -\\omega_z \\\\
\\omega_x & 0 & \\omega_z & -\\omega_y \\\\
\\omega_y & -\\omega_z & 0 & \\omega_x \\\\
\\omega_z & \\omega_y & -\\omega_x & 0
\\end{bmatrix}
```

## Rotation Matrix from Quaternion
```math
\\mathbf{R} = \\begin{bmatrix}
1-2(q_y^2+q_z^2) & 2(q_xq_y-q_zq_w) & 2(q_xq_z+q_yq_w) \\\\
2(q_xq_y+q_zq_w) & 1-2(q_x^2+q_z^2) & 2(q_yq_z-q_xq_w) \\\\
2(q_xq_z-q_yq_w) & 2(q_yq_z+q_xq_w) & 1-2(q_x^2+q_y^2)
\\end{bmatrix}
```

## Euler Rotation Equations
```math
\\begin{aligned}
\\alpha_x &= \\frac{M_x + (I_y - I_z)\\omega_y\\omega_z}{I_x} \\\\
\\alpha_y &= \\frac{M_y - d_y\\omega_y + (I_z - I_x)\\omega_z\\omega_x}{I_y} \\\\
\\alpha_z &= \\frac{M_z + M_{disturb} + (I_x - I_y)\\omega_x\\omega_y}{I_z}
\\end{aligned}
```

## Translational Dynamics (Newton's 2nd Law)
```math
m\\mathbf{a} = \\mathbf{F}_{aero,w} + \\mathbf{F}_{tether,w} + \\mathbf{F}_{gravity}
```

where aerodynamic force is transformed from body to world frame:
```math
\\mathbf{F}_{aero,w} = \\mathbf{R}_{b\\to w} \\mathbf{F}_{aero,b}
```

## Optional Sphere Constraint
When `fix_sphere=true`, motion is constrained to radial direction:
```math
\\begin{aligned}
\\frac{d\\mathbf{r}}{dt} &= (\\mathbf{v} \\cdot \\hat{\\mathbf{r}})\\hat{\\mathbf{r}} \\\\
\\frac{d\\mathbf{v}}{dt} &= (\\mathbf{a} \\cdot \\hat{\\mathbf{r}})\\hat{\\mathbf{r}}
\\end{aligned}
```

# Usage Example
```julia
using ModelingToolkit, SymbolicAWEModels

@named wing = RigidWing(
    mass = 25.0,
    I_xx = 10.0,
    I_yy = 15.0,
    I_zz = 20.0,
    y_damping = 150.0
)

# Set aerodynamic forces (would come from VSM in full implementation)
eqs = [
    wing.aero_force_b ~ [100.0, 0.0, 50.0]  # [N] in body frame
    wing.aero_moment_b ~ [0.0, 5.0, 0.0]    # [Nm] in body frame
    wing.tether_force_w ~ [0.0, 0.0, -500.0]  # [N] in world frame
]
```

# Notes
- Quaternion normalization is enforced through the differential equations
- The y-axis damping prevents unrealistic lateral oscillations
- Disturbance torque can be used for testing or turbulence modeling
- Future versions will include attachment point connectors for bridle forces
"""
@mtkmodel RigidWing begin
    @parameters begin
        mass = 25.0, [description = "Wing mass [kg]"]
        I_xx = 10.0, [description = "Moment of inertia about x-axis [kg⋅m²]"]
        I_yy = 15.0, [description = "Moment of inertia about y-axis [kg⋅m²]"]
        I_zz = 20.0, [description = "Moment of inertia about z-axis [kg⋅m²]"]
        y_damping = 150.0, [description = "Lateral angular damping [Ns]"]
        z_disturb = 0.0, [description = "Disturbance torque about z [Nm]"]
        fix_sphere = false, [description = "Constrain to sphere (radial only)"]
        g_earth = 9.81, [description = "Gravitational acceleration [m/s²]"]
        # External force/moment inputs (to be connected to aerodynamics/tethers)
        aero_force_b_ext[1:3] = zeros(3), [
            description = "External aero force in body frame [N]"
        ]
        aero_moment_b_ext[1:3] = zeros(3), [
            description = "External aero moment in body frame [Nm]"
        ]
        tether_force_w_ext[1:3] = zeros(3), [
            description = "External tether force in world frame [N]"
        ]
        tether_moment_w_ext[1:3] = zeros(3), [
            description = "External tether moment in world frame [Nm]"
        ]
    end

    @variables begin
        # Translational states
        pos(t)[1:3], [description = "Position in world frame [m]"]
        vel(t)[1:3], [description = "Velocity in world frame [m/s]"]
        acc(t)[1:3], [description = "Acceleration [m/s²]"]

        # Rotational states
        Q_b_w(t)[1:4], [description = "Quaternion (body to world) [w,x,y,z]"]
        ω_b(t)[1:3], [description = "Angular velocity in body frame [rad/s]"]
        α_b(t)[1:3], [description = "Angular acceleration in body frame [rad/s²]"]

        # Derived quantities
        R_b_w(t)[1:3, 1:3], [description = "Rotation matrix (body to world)"]
        Q_vel(t)[1:4], [description = "Quaternion derivative"]
        ω_b_stable(t)[1:3], [description = "Constrained angular velocity"]
        α_b_damped(t)[1:3], [description = "Damped angular acceleration"]

        # Forces and moments
        aero_force_b(t)[1:3], [description = "Aero force in body frame [N]"]
        aero_moment_b(t)[1:3], [description = "Aero moment in body frame [Nm]"]
        tether_force_w(t)[1:3], [description = "Tether force in world frame [N]"]
        tether_moment_w(t)[1:3], [description = "Tether moment in world frame [Nm]"]
        tether_moment_b(t)[1:3], [description = "Tether moment in body frame [Nm]"]

        # Constraint helpers
        axis(t)[1:3], [description = "Radial axis for sphere constraint"]
        axis_b(t)[1:3], [description = "Radial axis in body frame"]
    end

    @equations begin
        # ==================== EXTERNAL FORCE/MOMENT INTERFACE ==================== #
        aero_force_b ~ aero_force_b_ext
        aero_moment_b ~ aero_moment_b_ext
        tether_force_w ~ tether_force_w_ext
        tether_moment_w ~ tether_moment_w_ext

        # ==================== QUATERNION TO ROTATION MATRIX ==================== #
        # R = quaternion_to_rotation_matrix(Q)
        R_b_w[1, 1] ~ 1 - 2 * (Q_b_w[3]^2 + Q_b_w[4]^2)
        R_b_w[1, 2] ~ 2 * (Q_b_w[2] * Q_b_w[3] - Q_b_w[4] * Q_b_w[1])
        R_b_w[1, 3] ~ 2 * (Q_b_w[2] * Q_b_w[4] + Q_b_w[3] * Q_b_w[1])
        R_b_w[2, 1] ~ 2 * (Q_b_w[2] * Q_b_w[3] + Q_b_w[4] * Q_b_w[1])
        R_b_w[2, 2] ~ 1 - 2 * (Q_b_w[2]^2 + Q_b_w[4]^2)
        R_b_w[2, 3] ~ 2 * (Q_b_w[3] * Q_b_w[4] - Q_b_w[2] * Q_b_w[1])
        R_b_w[3, 1] ~ 2 * (Q_b_w[2] * Q_b_w[4] - Q_b_w[3] * Q_b_w[1])
        R_b_w[3, 2] ~ 2 * (Q_b_w[3] * Q_b_w[4] + Q_b_w[2] * Q_b_w[1])
        R_b_w[3, 3] ~ 1 - 2 * (Q_b_w[2]^2 + Q_b_w[3]^2)

        # ==================== QUATERNION KINEMATICS ==================== #
        # dQ/dt = 0.5 * Ω(ω) * Q, where Ω is the skew-symmetric matrix
        Q_vel[1] ~ 0.5 * (
            -ω_b_stable[1] * Q_b_w[2] - ω_b_stable[2] * Q_b_w[3] -
            ω_b_stable[3] * Q_b_w[4]
        )
        Q_vel[2] ~ 0.5 * (
            ω_b_stable[1] * Q_b_w[1] + ω_b_stable[3] * Q_b_w[3] -
            ω_b_stable[2] * Q_b_w[4]
        )
        Q_vel[3] ~ 0.5 * (
            ω_b_stable[2] * Q_b_w[1] - ω_b_stable[3] * Q_b_w[2] +
            ω_b_stable[1] * Q_b_w[4]
        )
        Q_vel[4] ~ 0.5 * (
            ω_b_stable[3] * Q_b_w[1] + ω_b_stable[2] * Q_b_w[2] -
            ω_b_stable[1] * Q_b_w[3]
        )

        D.(Q_b_w) .~ Q_vel

        # ==================== SPHERE CONSTRAINT ==================== #
        axis ~ pos / max(1e-6, sqrt(sum(pos .^ 2)))
        axis_b ~ [
            sum(R_b_w[1, :] .* axis)
            sum(R_b_w[2, :] .* axis)
            sum(R_b_w[3, :] .* axis)
        ]

        # Constrain angular velocity if sphere constraint is active
        if fix_sphere
            ω_b_stable ~ ω_b - (sum(ω_b .* axis_b)) * axis_b
        else
            ω_b_stable ~ ω_b
        end

        # ==================== EULER ROTATION EQUATIONS ==================== #
        # Apply damping and disturbances
        α_b_damped[1] ~ α_b[1]
        α_b_damped[2] ~ α_b[2] - y_damping * ω_b[2]
        α_b_damped[3] ~ α_b[3] + z_disturb

        # Transform tether moment to body frame
        tether_moment_b ~ [
            sum(R_b_w[1, :] .* tether_moment_w)
            sum(R_b_w[2, :] .* tether_moment_w)
            sum(R_b_w[3, :] .* tether_moment_w)
        ]

        # Euler equations with cross-product terms
        α_b[1] ~ (
            aero_moment_b[1] + tether_moment_b[1] +
            (I_yy - I_zz) * ω_b[2] * ω_b[3]
        ) / I_xx
        α_b[2] ~ (
            aero_moment_b[2] + tether_moment_b[2] +
            (I_zz - I_xx) * ω_b[3] * ω_b[1]
        ) / I_yy
        α_b[3] ~ (
            aero_moment_b[3] + tether_moment_b[3] +
            (I_xx - I_yy) * ω_b[1] * ω_b[2]
        ) / I_zz

        # Angular velocity evolution (with sphere constraint)
        if fix_sphere
            D.(ω_b) .~ α_b_damped - (sum(α_b_damped .* axis_b)) * axis_b
        else
            D.(ω_b) .~ α_b_damped
        end

        # ==================== TRANSLATIONAL DYNAMICS ==================== #
        # Total force: aerodynamic (body frame) + tether (world frame) + gravity
        acc ~ (
            [sum(R_b_w[1, :] .* aero_force_b);
             sum(R_b_w[2, :] .* aero_force_b);
             sum(R_b_w[3, :] .* aero_force_b)] +
            tether_force_w + [0, 0, -mass * g_earth]
        ) / mass

        # Position/velocity evolution (with sphere constraint)
        if fix_sphere
            D.(pos) .~ (sum(vel .* axis)) * axis
            D.(vel) .~ (sum(acc .* axis)) * axis
        else
            D.(pos) .~ vel
            D.(vel) .~ acc
        end
    end
end
