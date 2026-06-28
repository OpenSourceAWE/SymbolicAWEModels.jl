# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

"""
    WeightedRefPoints

Weighted combination of reference points for body-frame definition. Supports
single points, equal-weight averaging, and arbitrary weight combinations.

# Fields
- `refs`: Unresolved names/indices (filled at construction)
- `ids`: Resolved point indices (filled by `resolve!`)
- `weights`: Normalized weights (sum to 1.0)
"""
mutable struct WeightedRefPoints
    const refs::Vector{NameRef}
    ids::Vector{Int64}
    const weights::Vector{Float64}
end

"""
    Body{A<:AbstractAeroModel, D<:WingDynamics}

A 6-DOF body integrated (or fitted) in the principal frame, optionally carrying
aerodynamics. A plain body has `aero = AeroNone()`; a "wing" carries a real aero
mode (see [`sys_struct.wings`]). `type` is DYNAMIC (free 6-DOF), KINEMATIC (pose
fitted from structural points — particle wings) or STATIC (clamped). `D` mirrors
the rigid/particle distinction into the type domain for aero dispatch. All bodies
share the 6-DOF generator `rigid_body_eqs!`.

The rigid-body core (`mass`, `inertia_principal`, frames, 6-DOF state) is set
directly or derived from the body's points; aero/wing fields are inert when
`aero` is [`AeroNone`](@ref). Loads are gravity (`-g·mass` for DYNAMIC bodies),
the settable external wrench, joint wrenches, and aerodynamics.

$(TYPEDFIELDS)
"""
mutable struct Body{A<:AbstractAeroModel, D<:WingDynamics}
    "Index in the bodies vector (assigned by SystemStructure)."
    idx::Int64
    "Name used for lookup by other components' `_ref` fields."
    const name::Union{Int, Symbol, Nothing}
    "Resolved transform index (filled by SystemStructure). 0 = no transform."
    transform_idx::Int64
    "Raw transform reference (name or idx). 0 = no transform."
    const transform_ref::Union{Int, Symbol}

    # ---- rigid-body core ----
    "Total mass [kg]."
    mass::SimFloat
    "Principal moments of inertia `[Ixx, Iyy, Izz]` [kg·m²]."
    const inertia_principal::KVec3
    "Constant body→principal rotation."
    const R_b_to_p::Matrix{SimFloat}
    "Principal frame → CAD (from inertia diagonalisation)."
    const R_p_to_c::Matrix{SimFloat}
    "Offset from body origin to COM, body frame [m]."
    const com_offset_b::KVec3
    "External force applied at the COM, world frame [N] (settable)."
    const ext_force_w::KVec3
    "External force applied at the COM, body frame [N] (settable)."
    const ext_force_b::KVec3
    "External moment applied about the COM, body frame [N·m] (settable)."
    const ext_moment_b::KVec3
    "Angular damping `[d_x, d_y, d_z]` [N·m·s] along the body-frame axes."
    damping::KVec3
    "If true, COM motion is confined to a sphere about the world origin."
    fix_sphere::Bool
    "Dynamics type: DYNAMIC (free 6-DOF), KINEMATIC (fitted) or STATIC (frozen)."
    type::DynamicsType
    "Initial body-origin position [m]; `pos_w` is reset to this each `reinit!`."
    const pos_cad::KVec3
    "Initial body→world orientation; `Q_b_to_w` is reset from this each `reinit!`."
    const R_b_to_c::Matrix{SimFloat}

    # ---- body frame state (ICs in, live output out) ----
    const Q_b_to_w::Vector{SimFloat}
    const ω_b::KVec3
    const pos_w::KVec3
    const vel_w::KVec3
    const acc_w::KVec3

    # ---- principal frame ODE state ----
    const com_w::KVec3
    const com_vel::KVec3
    const Q_p_to_w::Vector{SimFloat}
    const ω_p::KVec3

    # ---- aerodynamics + wing machinery (inert when aero === AeroNone) ----
    aero::A
    "Resolved twist_surface indices (filled by SystemStructure)."
    twist_surface_idxs::Vector{Int64}
    "Raw twist_surface references."
    const twist_surface_refs::Vector{NameRef}
    const wind_disturb::KVec3
    drag_frac::SimFloat
    const va_b::KVec3
    const v_wind::KVec3
    const aero_force_b::KVec3
    const aero_moment_b::KVec3
    const tether_moment::KVec3
    const tether_force::KVec3
    elevation::SimFloat
    elevation_vel::SimFloat
    elevation_acc::SimFloat
    azimuth::SimFloat
    azimuth_vel::SimFloat
    azimuth_acc::SimFloat
    heading::SimFloat
    const turn_rate::KVec3
    const turn_acc::KVec3
    course::SimFloat
    aoa::SimFloat
    "Whether in-group (twist_surface) points contribute their moment to the body."
    group_points_moment::Bool
    # Body-frame reference points (define R_b_to_w / pos_w from structural points)
    z_ref_points::Union{Nothing, Tuple{WeightedRefPoints, WeightedRefPoints}}
    y_ref_points::Union{Nothing, Tuple{WeightedRefPoints, WeightedRefPoints}}
    origin::Union{Nothing, WeightedRefPoints}
end

"""`Body` with rigid (DYNAMIC/STATIC) 6-DOF dynamics."""
const RigidWing{A} = Body{A, RigidDynamics}
"""`Body` whose pose is fitted from structural points (particle/KINEMATIC)."""
const ParticleWing{A} = Body{A, ParticleDynamics}

"""
    principal_frame(inertia) -> (inertia_principal, R_to_principal)

Diagonalise a 3×3 symmetric inertia tensor. Returns the principal moments
`inertia_principal` and the rotation `R_to_principal` mapping the input frame to
the principal frame, so that `R · inertia · R' = Diagonal(inertia_principal)`.

The principal axes are permuted and signed to align as closely as possible with
the input-frame axes: a near-diagonal tensor gives `R ≈ I`, and a body symmetric
about a coordinate plane gives a pure rotation about the normal of that plane.
`R` is always a proper rotation (`det = +1`).
"""
function principal_frame(inertia::AbstractMatrix)
    decomposition = eigen(Symmetric(Matrix{SimFloat}(inertia)))
    moments = decomposition.values
    axes = Matrix(decomposition.vectors)   # columns: principal axes, input frame
    # Permute columns so principal axis i aligns with input axis i (R ≈ I).
    best_perm, best_score = (1, 2, 3), -Inf
    for perm in ((1,2,3), (1,3,2), (2,1,3), (2,3,1), (3,1,2), (3,2,1))
        score = abs(axes[1, perm[1]]) + abs(axes[2, perm[2]]) +
                abs(axes[3, perm[3]])
        score > best_score && ((best_perm, best_score) = (perm, score))
    end
    order = collect(best_perm)
    axes, moments = axes[:, order], moments[order]
    # Sign each axis to point along its positive input axis.
    for i in 1:3
        axes[i, i] < 0 && (axes[:, i] .*= -1)
    end
    # Guarantee a proper rotation by flipping the least-aligned axis if reflected.
    if det(axes) < 0
        flip = argmin([abs(axes[i, i]) for i in 1:3])
        axes[:, flip] .*= -1
    end
    return Vector{SimFloat}(moments), Matrix{SimFloat}(axes')
end

"""
    Body(name; mass, inertia_principal | inertia, pos, vel=zeros,
              Q_b_to_w=[1,0,0,0], ω_b=zeros, com_offset_b=zeros, R_b_to_p=I,
              angular_damping=0, fix_sphere=false, type=DYNAMIC, transform=nothing,
              ext_force_w=zeros, ext_moment_b=zeros)

Construct a standalone rigid body. `pos`/`vel` are the body origin's initial
world-frame position/velocity; `Q_b_to_w`/`ω_b` its initial orientation/spin.
The principal-frame ODE state is derived from these by `init_rigid_body!`.

`type` is `DYNAMIC` (free 6-DOF, default) or `STATIC` (clamped to its initial
pose — e.g. a cantilever root). `transform` optionally references a
[`Transform`](@ref) that repositions/rotates the body's initial pose (azimuth,
elevation, heading), like a wing.

`angular_damping` is a scalar (isotropic) or length-3 per-axis principal-frame
damping. `fix_sphere=true` confines COM motion to a sphere about the world
origin (radial DOF frozen), like a `RIGID_DYNAMICS` wing's `fix_sphere`.

Supply the inertia in one of two ways: `inertia_principal` (a length-3 diagonal
principal inertia, with `R_b_to_p` giving the body→principal rotation), or
`inertia` (a full 3×3 body-frame tensor), in which case both `inertia_principal`
and `R_b_to_p` are derived via [`principal_frame`](@ref). Give one, not both.
"""
function Body(name;
        mass::Real,
        inertia_principal = nothing,
        inertia = nothing,
        pos,
        vel = zeros(SimFloat, 3),
        Q_b_to_w = SimFloat[1, 0, 0, 0],
        ω_b = zeros(SimFloat, 3),
        com_offset_b = zeros(SimFloat, 3),
        R_b_to_p = Matrix{SimFloat}(I, 3, 3),
        damping = 0.0,
        fix_sphere::Bool = false,
        type::DynamicsType = DYNAMIC,
        transform = nothing,
        ext_force_w = zeros(SimFloat, 3),
        ext_force_b = zeros(SimFloat, 3),
        ext_moment_b = zeros(SimFloat, 3),
    )
    # Scalar damping broadcasts to an isotropic per-axis vector.
    damping_vec = damping isa Real ? KVec3(damping, damping, damping) : KVec3(damping)
    type in (DYNAMIC, STATIC) || error(
        "Body $name: type must be DYNAMIC or STATIC, got $type.")
    transform_ref = isnothing(transform) ? 0 : transform
    if !isnothing(inertia)
        isnothing(inertia_principal) || error(
            "Body $name: give `inertia` or `inertia_principal`, not both.")
        inertia_principal, R_b_to_p = principal_frame(inertia)
    elseif isnothing(inertia_principal)
        error("Body $name: provide `inertia_principal` or `inertia`.")
    end
    R_b_to_c = quaternion_to_rotation_matrix(Vector{SimFloat}(Q_b_to_w))
    # Plain body: no aero (AeroNone), rigid dynamics, inert aero/wing fields.
    return Body{AeroNone, RigidDynamics}(
        0, name, 0, transform_ref,
        SimFloat(mass), KVec3(inertia_principal), Matrix{SimFloat}(R_b_to_p),
        Matrix{SimFloat}(I, 3, 3), KVec3(com_offset_b),
        KVec3(ext_force_w), KVec3(ext_force_b), KVec3(ext_moment_b),
        damping_vec, fix_sphere, type,
        KVec3(pos), Matrix{SimFloat}(R_b_to_c),
        Vector{SimFloat}(Q_b_to_w), KVec3(ω_b),
        KVec3(pos), KVec3(vel), zeros(KVec3),
        zeros(KVec3), zeros(KVec3), zeros(SimFloat, 4), zeros(KVec3),
        # aero/wing fields (inert for a plain body)
        AeroNone(), Int64[], NameRef[],
        zeros(KVec3), one(SimFloat),
        zeros(KVec3), zeros(KVec3), zeros(KVec3), zeros(KVec3), zeros(KVec3), zeros(KVec3),
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        zeros(KVec3), zeros(KVec3), 0.0, 0.0,
        true, nothing, nothing, nothing)
end

"""
    init_principal_state!(obj)

Derive the principal-frame ODE state (`com_w`, `com_vel`, `Q_p_to_w`, `ω_p`) from
the body-frame initial conditions (`pos_w`, `vel_w`, `Q_b_to_w`, `ω_b`) of a rigid
body or `RIGID_DYNAMICS` wing, with `R_p_to_w = R_b_to_w * R_b_to_p'`. Shared by
[`init_rigid_body!`](@ref) and the rigid branch of `init_principal_frame!` (for a
`RIGID_DYNAMICS` wing `R_b_to_p = R_p_to_c' * R_b_to_c`, so the two agree).
"""
function init_principal_state!(obj)
    R_b_to_w = quaternion_to_rotation_matrix(obj.Q_b_to_w)
    obj.com_w .= obj.pos_w .+ R_b_to_w * obj.com_offset_b
    R_p_to_w = R_b_to_w * obj.R_b_to_p'
    obj.Q_p_to_w .= rotation_matrix_to_quaternion(R_p_to_w)
    ω_w = R_b_to_w * obj.ω_b
    obj.com_vel .= obj.vel_w .+ cross(ω_w, R_b_to_w * obj.com_offset_b)
    obj.ω_p .= R_p_to_w' * ω_w
    return nothing
end

"""
    init_rigid_body!(body::Body)

Derive the body's principal-frame ODE state from its body-frame initial
conditions via [`init_principal_state!`](@ref).
"""
init_rigid_body!(body::Body) = init_principal_state!(body)

"""
    ElasticJoint

A 6-DOF elastic connection between two `Body`s. Anchored at a body-frame
offset on each body, it applies a restoring wrench from the relative pose of the
anchors, decomposed in body A's frame into axial (`EA`), shear (`GA`, both
transverse axes), torsion (`GJ`), and bending (`EI`, both transverse axes), with
optional translational/rotational damping. The equal-and-opposite wrench is added
to both bodies' load accumulators.

Each stiffness is either a `Real` (linear law, force `= k·Δ`) or a callable
interpolation `f` (nonlinear law, force `= f(Δ)`, e.g. a wrinkling/saturating
inflatable beam). They may be mixed per DOF; the type parameter `S` keeps the
fields concrete (`Real`s share one type, interpolations share one type) so the
ODE right-hand side stays allocation-free. The symbolic equations are identical
for `Real` or interpolation stiffnesses (the choice is resolved at runtime in the
registered force function), but the stiffness types are part of the
`SystemStructure` type parameter, so a float vs an interpolation is a distinct
compiled model with its own cache entry.

$(TYPEDFIELDS)
"""
mutable struct ElasticJoint{S}
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

    "Axial stiffness: `Real` EA [N/m], or interpolation force(Δx) (body A x-axis)."
    stiffness_axial::S
    "Shear stiffness: `Real` GA [N/m], or interpolation force(Δ) (both transverse)."
    stiffness_shear::S
    "Torsional stiffness: `Real` GJ [N·m/rad], or interpolation moment(Δθ) (x-axis)."
    stiffness_torsion::S
    "Bending stiffness: `Real` EI [N·m/rad], or interpolation moment(Δθ) (transverse)."
    stiffness_bending::S
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
`anchor_a`/`anchor_b` are the connection points in each body's frame. Each
stiffness is a `Real` (linear) or a callable interpolation of the relative DOF
(nonlinear); interpolations must all be the same type.
"""
function ElasticJoint(name, body_a, body_b;
        anchor_a = zeros(SimFloat, 3),
        anchor_b = zeros(SimFloat, 3),
        stiffness_axial,
        stiffness_shear,
        stiffness_torsion,
        stiffness_bending,
        damping_trans::Real = 0.0,
        damping_rot::Real = 0.0,
    )
    # Reals → SimFloat; interpolations kept as-is. All interps must share a type.
    conv(s) = s isa Real ? SimFloat(s) : s
    stiffs = map(conv, (stiffness_axial, stiffness_shear,
                        stiffness_torsion, stiffness_bending))
    interp_types = unique(typeof(s) for s in stiffs if !(s isa Real))
    length(interp_types) > 1 && error(
        "ElasticJoint: all interpolation stiffnesses must be the same type, " *
        "got $(interp_types). Mix only `Real`s and a single interpolation type.")
    S = Union{map(typeof, stiffs)...}
    body_a_ref = body_a isa Integer ? Int(body_a) : Symbol(body_a)
    body_b_ref = body_b isa Integer ? Int(body_b) : Symbol(body_b)
    return ElasticJoint{S}(0, name, 0, 0, body_a_ref, body_b_ref,
        KVec3(anchor_a), KVec3(anchor_b), stiffs...,
        SimFloat(damping_trans), SimFloat(damping_rot))
end
