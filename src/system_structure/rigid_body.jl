# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Standalone rigid body component. Shares the 6-DOF dynamics generator
# `rigid_body_eqs!` with RIGID_DYNAMICS wings, but carries no aero, tether, or
# transform machinery — it is specified directly by mass, principal inertia, and
# initial conditions, and driven by gravity plus a settable external wrench.

"""
    RigidBody

A free 6-DOF rigid body integrated in the principal frame. Use it on its own or
as a link in a multi-body chain (connected by elastic joints).

The body is specified directly: `mass`, `inertia_principal` (diagonal, principal
frame), initial pose/twist (`pos_w`, `vel_w`, `Q_b_to_w`, `ω_b`), and optional
`com_offset_b`/`R_b_to_p` relating the body frame to the principal frame. Loads
are gravity plus the settable external wrench (`ext_force_w`, `ext_moment_b`),
which `@register_symbolic` getters read live each RHS evaluation.

$(TYPEDFIELDS)
"""
mutable struct RigidBody
    "Index in the rigid_bodies vector (assigned by SystemStructure)."
    idx::Int64
    "Name used for lookup by other components' `_ref` fields."
    const name::Union{Int, Symbol, Nothing}

    "Total mass [kg]."
    mass::SimFloat
    "Principal moments of inertia `[Ixx, Iyy, Izz]` [kg·m²]."
    const inertia_principal::KVec3
    "Constant body→principal rotation."
    const R_b_to_p::Matrix{SimFloat}
    "Offset from body origin to COM, body frame [m]."
    const com_offset_b::KVec3

    "External force applied at the COM, world frame [N] (settable)."
    const ext_force_w::KVec3
    "External moment applied about the COM, body frame [N·m] (settable)."
    const ext_moment_b::KVec3
    "Isotropic angular damping coefficient [N·m·s] (principal frame)."
    angular_damping::SimFloat
    "If true, the body is held at its initial pose (all DOF frozen)."
    fixed::Bool

    # Body frame state: initial conditions in, live output out.
    "Body→world quaternion. Initial condition in; algebraic output out."
    const Q_b_to_w::Vector{SimFloat}
    "Angular velocity, body frame [rad/s]. Initial condition in; output out."
    const ω_b::KVec3
    "Body origin position, world frame [m]. Initial condition in; output out."
    const pos_w::KVec3
    "Body origin velocity, world frame [m/s]. Initial condition in; output out."
    const vel_w::KVec3
    "Body origin acceleration, world frame [m/s²] (output)."
    const acc_w::KVec3

    # Principal frame ODE state (derived at init, then integrated).
    "COM position, world frame [m]."
    const com_w::KVec3
    "COM velocity, world frame [m/s]."
    const com_vel::KVec3
    "Principal→world quaternion."
    const Q_p_to_w::Vector{SimFloat}
    "Angular velocity, principal frame [rad/s]."
    const ω_p::KVec3
end

"""
    RigidBody(name; mass, inertia_principal, pos, vel=zeros, Q_b_to_w=[1,0,0,0],
              ω_b=zeros, com_offset_b=zeros, R_b_to_p=I, angular_damping=0,
              ext_force_w=zeros, ext_moment_b=zeros)

Construct a standalone rigid body. `pos`/`vel` are the body origin's initial
world-frame position/velocity; `Q_b_to_w`/`ω_b` its initial orientation/spin.
The principal-frame ODE state is derived from these by `init_rigid_body!`.
"""
function RigidBody(name;
        mass::Real,
        inertia_principal,
        pos,
        vel = zeros(SimFloat, 3),
        Q_b_to_w = SimFloat[1, 0, 0, 0],
        ω_b = zeros(SimFloat, 3),
        com_offset_b = zeros(SimFloat, 3),
        R_b_to_p = Matrix{SimFloat}(I, 3, 3),
        angular_damping::Real = 0.0,
        fixed::Bool = false,
        ext_force_w = zeros(SimFloat, 3),
        ext_moment_b = zeros(SimFloat, 3),
    )
    return RigidBody(0, name,
        SimFloat(mass), KVec3(inertia_principal),
        Matrix{SimFloat}(R_b_to_p), KVec3(com_offset_b),
        KVec3(ext_force_w), KVec3(ext_moment_b), SimFloat(angular_damping), fixed,
        Vector{SimFloat}(Q_b_to_w), KVec3(ω_b),
        KVec3(pos), KVec3(vel), zeros(KVec3),
        zeros(KVec3), zeros(KVec3), zeros(SimFloat, 4), zeros(KVec3))
end

"""
    init_rigid_body!(body::RigidBody)

Derive the principal-frame ODE state (`com_w`, `com_vel`, `Q_p_to_w`, `ω_p`) from
the body-frame initial conditions (`pos_w`, `vel_w`, `Q_b_to_w`, `ω_b`). Mirrors
`init_principal_frame!` for wings, with `R_p_to_w = R_b_to_w * R_b_to_p'`.
"""
function init_rigid_body!(body::RigidBody)
    R_b_to_w = quaternion_to_rotation_matrix(body.Q_b_to_w)
    body.com_w .= body.pos_w .+ R_b_to_w * body.com_offset_b
    R_p_to_w = R_b_to_w * body.R_b_to_p'
    body.Q_p_to_w .= rotation_matrix_to_quaternion(R_p_to_w)
    ω_w = R_b_to_w * body.ω_b
    r_com_w = R_b_to_w * body.com_offset_b
    body.com_vel .= body.vel_w .+ cross(ω_w, r_com_w)
    body.ω_p .= R_p_to_w' * ω_w
    return nothing
end

"""
    ElasticJoint

A 6-DOF elastic connection between two `RigidBody`s. Anchored at a body-frame
offset on each body, it applies a restoring wrench proportional to the relative
pose of the anchors, decomposed in body A's frame into axial (`EA`), shear
(`GA`, both transverse axes), torsion (`GJ`), and bending (`EI`, both transverse
axes), with optional translational/rotational damping. The equal-and-opposite
wrench is added to both bodies' load accumulators.

The stiffness law is linear here; a single tunable function of the relative DOF
keeps a nonlinear/pressure-dependent law (e.g. `EI(p)`) a drop-in later.

$(TYPEDFIELDS)
"""
mutable struct ElasticJoint
    "Index in the elastic_joints vector (assigned by SystemStructure)."
    idx::Int64
    "Name used for lookup."
    const name::Union{Int, Symbol, Nothing}

    "Resolved index of body A (filled by SystemStructure)."
    body_a_idx::Int64
    "Resolved index of body B (filled by SystemStructure)."
    body_b_idx::Int64
    "Raw reference (name or index) of body A."
    const body_a_ref::NameRef
    "Raw reference (name or index) of body B."
    const body_b_ref::NameRef

    "Anchor offset from body A origin, body A frame [m]."
    const anchor_a_b::KVec3
    "Anchor offset from body B origin, body B frame [m]."
    const anchor_b_b::KVec3

    "Axial stiffness EA [N/m] (body A x-axis)."
    stiffness_axial::SimFloat
    "Shear stiffness GA [N/m] (both transverse axes)."
    stiffness_shear::SimFloat
    "Torsional stiffness GJ [N·m/rad] (about body A x-axis)."
    stiffness_torsion::SimFloat
    "Bending stiffness EI [N·m/rad] (both transverse axes)."
    stiffness_bending::SimFloat
    "Translational damping [N·s/m]."
    damping_trans::SimFloat
    "Rotational damping [N·m·s/rad]."
    damping_rot::SimFloat
end

"""
    ElasticJoint(name, body_a, body_b; anchor_a=zeros, anchor_b=zeros,
                 stiffness_axial, stiffness_shear, stiffness_torsion,
                 stiffness_bending, damping_trans=0, damping_rot=0)

Connect `body_a` to `body_b` (names or indices) with a 6-DOF elastic joint.
`anchor_a`/`anchor_b` are the connection points in each body's frame.
"""
function ElasticJoint(name, body_a, body_b;
        anchor_a = zeros(SimFloat, 3),
        anchor_b = zeros(SimFloat, 3),
        stiffness_axial::Real,
        stiffness_shear::Real,
        stiffness_torsion::Real,
        stiffness_bending::Real,
        damping_trans::Real = 0.0,
        damping_rot::Real = 0.0,
    )
    body_a_ref = body_a isa Integer ? Int(body_a) : Symbol(body_a)
    body_b_ref = body_b isa Integer ? Int(body_b) : Symbol(body_b)
    return ElasticJoint(0, name, 0, 0, body_a_ref, body_b_ref,
        KVec3(anchor_a), KVec3(anchor_b),
        SimFloat(stiffness_axial), SimFloat(stiffness_shear),
        SimFloat(stiffness_torsion), SimFloat(stiffness_bending),
        SimFloat(damping_trans), SimFloat(damping_rot))
end
