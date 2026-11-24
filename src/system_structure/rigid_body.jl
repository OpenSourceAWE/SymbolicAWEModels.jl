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
    BodyModel

What SymbolicAWEModels holds of a body beyond its KiteGeometry `bodies` row: the columns
only it reads, each `nothing` where the row leaves it unset.

$(TYPEDFIELDS)
"""
@kwdef mutable struct BodyModel <: KiteGeometry.AbstractModel
    "The wing whose body frame the body's damping acts in, by name or index."
    wing::Union{Nothing, KiteGeometry.NameRef} = nothing
    "Spin damping about the body axes [1/s], one value or one per axis."
    angular_damping::Union{Nothing, Float64, Vector{Float64}} = nothing
    "Damping per unit mass on the world axes [1/s], one value or one per axis."
    world_frame_damping::Union{Nothing, Float64, Vector{Float64}} = nothing
    "Damping per unit mass on the body axes [1/s], one value or one per axis."
    body_frame_damping::Union{Nothing, Float64, Vector{Float64}} = nothing
    "Whether the body is held on a sphere."
    fix_sphere::Union{Nothing, Bool} = nothing
    "Whether the body is frozen in place."
    fix_static::Union{Nothing, Bool} = nothing
    "External force, world frame [N]."
    ext_force_w::Union{Nothing, Vector{Float64}} = nothing
    "External force, body frame [N]."
    ext_force_b::Union{Nothing, Vector{Float64}} = nothing
    "External moment, body frame [N·m]."
    ext_moment_b::Union{Nothing, Vector{Float64}} = nothing
    "How the principal frame is found: `EIGEN_DECOMP` or `Y_ROTATION`."
    principal_frame_method::Union{Nothing, String} = nothing
    "Entrained-air mass resisting acceleration without adding weight [kg]."
    apparent_mass::Union{Nothing, Float64} = nothing
    "Centre of the body's own mass from the origin of its frame, body frame [m]."
    com_offset_b::Union{Nothing, Vector{Float64}} = nothing
end

KiteGeometry.model_columns(::Type{BodyModel}) = fieldnames(BodyModel)
KiteGeometry.model_references(::Type{BodyModel}) = (wing=:wings,)

"""
    Body

A 6-DOF body integrated (or fitted) in the principal frame: a KiteGeometry `bodies` row
whose model is a [`BodyModel`](@ref). The body of a wing is the one its `wings` row names.
Its dynamics type is DYNAMIC (free 6-DOF), KINEMATIC (pose fitted from its points, the
body of a particle wing) or STATIC (held by its `fixed` column).

The body's own mass properties are its `extra_mass` and `extra_inertia_KA`. A DYNAMIC or
STATIC body adds the points it carries, its `BODY_STATIC` riders and a `RIGID_DYNAMICS`
wing's nodes, as point masses of their `total_mass` at their anchors: a run's
`total_mass`, `com_offset_b` and principal inertia describe the body with them
([`update_mass_properties!`](@ref)). The row's `pos_ENU` is the centre of the body's own
mass; the origin of its frame is `origin_ENU`.
"""
const Body = KiteGeometry.Body{BodyModel}

"""The `wings` row of the wing whose body body record `r` is, `nothing` for another."""
function body_wing(r)
    wing = body_geometry(r).wing_row
    return wing == 0 ? nothing : getfield(run_state(r), :structure).wings[wing]
end

"""The dynamics type of body record `body`: `STATIC` where it is held, `KINEMATIC` for
the body of a particle wing and `DYNAMIC` for another."""
function body_type(body)
    body.fixed && return STATIC
    return body_dynamics_type(body) == PARTICLE_DYNAMICS ? KINEMATIC : DYNAMIC
end

"""The [`WingType`](@ref) of body record `body`: its wing's, `RIGID_DYNAMICS` for a body
that is no wing."""
body_dynamics_type(body) = body_geometry(body).dynamics_type

"""Column `column` of the model of body record `r`, else the same column of its wing's
`wings` row, `nothing` where neither gives it."""
function wing_column(r, column)
    model = getfield(r, :row).model
    value = hasproperty(model, column) ? getproperty(model, column) : nothing
    wing = body_wing(r)
    isnothing(value) && !isnothing(wing) && return getproperty(wing.model, column)
    return value
end

"""How the principal frame of body record `r` is found: its own column, else its wing's,
else `EIGEN_DECOMP`."""
body_principal_frame_method(r) = parse_principal_frame_method(
    something(wing_column(r, :principal_frame_method), "EIGEN_DECOMP"))

"""The spin damping of body record `r` per axis [1/s]: its own column, else its wing's,
else none on a body that is no wing and 150 about y on a wing."""
body_angular_damping(r) = per_axis(something(wing_column(r, :angular_damping),
    isnothing(body_wing(r)) ? 0.0 : [0.0, 150.0, 0.0]))

"""The centre of `body`'s own mass from the origin of its frame, body frame [m]: a row's
`com_offset_b` column, what the run works out for a body record."""
own_com_offset(body::KiteGeometry.Body) =
    something(body.model.com_offset_b, zero(SVector{3, Float64}))
own_com_offset(body) = body_geometry(body).extra_com_offset_b

"""
    origin_ENU(body)

The origin of `body`'s frame in the world [m] in its initial pose: its `pos_ENU`, the
centre of its own mass, less that centre's offset from the origin.
"""
origin_ENU(body) = KVec3(body.pos_ENU) -
    quaternion_to_rotation_matrix(body.Q_KA2ENU) * own_com_offset(body)

"""Put the origin of body row `body`'s frame at `origin` [m] in the world, keeping its
rotation."""
set_origin_ENU!(body, origin) = body.pos_ENU =
    origin + quaternion_to_rotation_matrix(body.Q_KA2ENU) * own_com_offset(body)

"""
    mesh_frame(r) -> (origin, R_mesh_to_w)

The frame the mesh of body record `r`'s wing is authored in: the KA frame of its wing's
`wings` row, origin [m] and rotation into the world, else the world frame.
"""
function mesh_frame(r)
    wing = body_wing(r)
    isnothing(wing) && return zeros(KVec3), Matrix{SimFloat}(I, 3, 3)
    return KVec3(wing.pos_ENU), quaternion_to_rotation_matrix(wing.Q_KA2ENU)
end

"""The vector column `field` of the model of `row`, a zero vector where it is unset."""
column_vector(row, field) =
    something(getproperty(row.model, field), zero(SVector{3, Float64}))

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

"""The principal moments of `inertia` and the rotation into them, found by `method`
(see [`PrincipalFrameMethod`](@ref))."""
principal_axes(inertia, method::PrincipalFrameMethod) =
    method == Y_ROTATION ? calc_inertia_y_rotation(inertia) : principal_frame(inertia)

"""
    point_mass_inertia(mass, r) -> Matrix

Inertia tensor [kg·m²] of a point `mass` [kg] at offset `r` [m] from the reference.
"""
point_mass_inertia(mass, r) = mass * (dot(r, r) * I(3) - r * r')

"""
    calc_inertia_y_rotation(I_tensor) -> (inertia_principal, R_to_principal)

Diagonalize a 3×3 inertia tensor via a closed-form rotation about the Y axis,
returning `(moments, R)` in the same format as [`principal_frame`](@ref).

The rotation angle is `θ = atan(2·I₁₃, I₁₁ − I₃₃) / 2`, zeroing out the
`I[1,3]` / `I[3,1]` cross terms while leaving the Y axis unchanged.

Unlike [`principal_frame`](@ref) (full 3-axis eigendecomposition + permutation
search), this is a unique, closed-form solution — the right choice for a wing
symmetric about the XZ-plane (no Y products of inertia), where the generic
permutation search is ambiguous when two principal moments are close.
"""
function calc_inertia_y_rotation(I_tensor)
    θ = atan(2 * I_tensor[1, 3], I_tensor[1, 1] - I_tensor[3, 3]) / 2
    cθ, sθ = cos(θ), sin(θ)
    Ry = SimFloat[cθ 0 sθ; 0 1 0; -sθ 0 cθ]
    I_diag = Ry * I_tensor * Ry'
    return Vector{SimFloat}(diag(I_diag)), Ry
end

"""
    broadcast_damping(damping) -> KVec3

Widen a scalar damping coefficient to an isotropic per-axis vector, or pass a
length-3 one through.
"""
broadcast_damping(damping) = damping isa Real ?
    KVec3(damping, damping, damping) : KVec3(damping)

"""
    Body(name; extra_mass, inertia_principal | inertia, pos,
              Q_b_to_w=[1,0,0,0], com_offset_b=zeros, R_b_to_p=I,
              angular_damping=0, world_frame_damping=0, body_frame_damping=0,
              fix_sphere=false, fix_static=false, type=DYNAMIC,
              wing=nothing, ext_force_w=zeros, ext_moment_b=zeros)

Construct a standalone rigid body. `pos` is the body origin's initial world-frame
position and `Q_b_to_w` its initial orientation.

`type` is `DYNAMIC` (free 6-DOF, default) or `STATIC` (clamped to its initial
pose — e.g. a cantilever root). `wing` names the parent wing whose frame and
motion the body-frame damping resolves against, as a `Point`'s `wing` does.

`angular_damping` is a scalar (isotropic) or length-3 per-axis spin damping [1/s].
`world_frame_damping` and `body_frame_damping` damp the COM velocity resolved on
the world and body axes respectively; both are per-mass [1/s], like a `Point`'s.
`fix_sphere=true` confines COM motion to a sphere about the world origin (radial
DOF frozen), like a `RIGID_DYNAMICS` wing's `fix_sphere`. `fix_static=true` holds
the body where it is; unlike `type=STATIC` it is a parameter, so it can be
toggled on a built model.

`extra_mass`, `com_offset_b` and the inertia are the body's own, without the
`BODY_STATIC` points riding it: `SystemDefinition` adds those. Supply the
inertia about the body's own COM in one of two ways: `inertia_principal` (a length-3
diagonal principal inertia, with `R_b_to_p` giving the body→principal rotation), or
`inertia` (a full 3×3 body-frame tensor). Give one, not both. Either is stored as a
tensor; a run's `inertia_principal` and `R_b_to_p` diagonalise the total with the riders,
so its `R_b_to_p` can differ from the one passed.
"""
function Body(name;
        extra_mass::Real,
        inertia_principal = nothing,
        inertia = nothing,
        pos,
        Q_b_to_w = SimFloat[1, 0, 0, 0],
        com_offset_b = zeros(SimFloat, 3),
        R_b_to_p = Matrix{SimFloat}(I, 3, 3),
        angular_damping = 0.0,
        world_frame_damping = 0.0,
        body_frame_damping = 0.0,
        fix_sphere::Bool = false,
        fix_static::Bool = false,
        type::DynamicsType = DYNAMIC,
        wing = nothing,
        ext_force_w = zeros(SimFloat, 3),
        ext_force_b = zeros(SimFloat, 3),
        ext_moment_b = zeros(SimFloat, 3),
        principal_frame_method::PrincipalFrameMethod = EIGEN_DECOMP,
    )
    damping_vec = broadcast_damping(angular_damping)
    world_damping_vec = broadcast_damping(world_frame_damping)
    body_damping_vec = broadcast_damping(body_frame_damping)
    type in (DYNAMIC, STATIC) || error(
        "Body $name: type must be DYNAMIC or STATIC, got $type.")
    if !isnothing(inertia)
        isnothing(inertia_principal) || error(
            "Body $name: give `inertia` or `inertia_principal`, not both.")
        inertia_principal, R_b_to_p = principal_axes(inertia, principal_frame_method)
    elseif isnothing(inertia_principal)
        error("Body $name: provide `inertia_principal` or `inertia`.")
    end
    extra_inertia_KA = R_b_to_p' * Diagonal(inertia_principal) * R_b_to_p
    model = BodyModel(; wing=cell_ref(wing), angular_damping=Vector(damping_vec),
        world_frame_damping=Vector(world_damping_vec),
        body_frame_damping=Vector(body_damping_vec), fix_sphere, fix_static,
        ext_force_w=Vector{Float64}(ext_force_w), ext_force_b=Vector{Float64}(ext_force_b),
        ext_moment_b=Vector{Float64}(ext_moment_b),
        principal_frame_method=string(principal_frame_method),
        com_offset_b=Vector{Float64}(com_offset_b))
    return KiteGeometry.Body(; name=Symbol(name),
        pos_ENU=pos + quaternion_to_rotation_matrix(Q_b_to_w) * com_offset_b,
        Q_KA2ENU=Q_b_to_w, extra_mass, extra_inertia_KA, fixed=type == STATIC, model)
end

"""
    initial_rotation(body) -> Matrix

The body→world rotation matrix of `body`'s initial pose `Q_KA2ENU`.
"""
initial_rotation(body) = quaternion_to_rotation_matrix(body.Q_KA2ENU)

"""
    initial_body_position(body, point) -> Vector

Where `point` sits in `body`'s frame, from the body origin, in the initial pose.
"""
initial_body_position(body, point) =
    initial_rotation(body)' * (point.pos_ENU - origin_ENU(body))

"""
    initial_chord_position(wing, point) -> SimFloat

Body-frame x of `point` on `wing` in the initial pose: smaller towards the leading edge.
"""
initial_chord_position(wing, point) = initial_body_position(wing, point)[1]

"""
    AbstractTubeModel

Structural element a [`Tube`](@ref) is simulated as, the model of its KiteGeometry `tubes`
row. Each concrete model carries its rigidities, its damping, the rest state it captures at
`place!`, the tube's index and where its ends sit on its bodies.
"""
abstract type AbstractTubeModel <: KiteGeometry.AbstractModel end

"""
    ElasticTube

Lumped 6-DOF spring between a tube's two anchors: the relative pose of the anchors,
decomposed in body A's frame into axial, shear (both transverse axes), torsion and
bending (both transverse axes) springs, with Rayleigh damping.

Each stiffness is either a `Real` (linear law, force `= k·Δ`) or a callable
interpolation `f` (nonlinear law, force `= f(Δ)`), mixable per DOF; callables must
share one type, which keeps the fields concrete.

$(TYPEDFIELDS)
"""
mutable struct ElasticTube{S} <: AbstractTubeModel
    "Axial stiffness: `Real` [N/m], or interpolation force(Δx) (body A x-axis)."
    stiffness_axial::S
    "Shear stiffness: `Real` [N/m], or interpolation force(Δ) (both transverse)."
    stiffness_shear::S
    "Torsional stiffness: `Real` [N·m/rad], or interpolation moment(Δθ) (x-axis)."
    stiffness_torsion::S
    "Bending stiffness: `Real` [N·m/rad], or interpolation moment(Δθ) (transverse)."
    stiffness_bending::S
    "Rayleigh damping β [s]: each DOF is damped in proportion to its own
    stiffness (C = βK), so rigid motion stays undamped by construction."
    damping::SimFloat
    "Rest anchor offset (body A frame) [m], the placed geometry, set by `place!`."
    const rest_offset_a::KVec3
    "Rest relative rotation `R_a' R_b`, the placed orientation, set by `place!`."
    const R_rel0::Matrix{SimFloat}
    "Index in `SystemDefinition.tubes` (assigned by SystemDefinition)."
    idx::Int64
    "Tube end on body A, as an offset from its origin in its frame [m]."
    const anchor_a_b::KVec3
    "Tube end on body B, as an offset from its origin in its frame [m]."
    const anchor_b_b::KVec3
end

"""
    ElasticTube(; stiffness_axial, stiffness_shear, stiffness_torsion,
                stiffness_bending, damping=0)

A lumped [`ElasticTube`](@ref) with the four given stiffnesses, each a `Real` or a
callable of its relative DOF.
"""
function ElasticTube(; stiffness_axial, stiffness_shear, stiffness_torsion,
                     stiffness_bending, damping::Real = 0.0)
    stiffnesses = concrete_rigidities("ElasticTube", (stiffness_axial,
        stiffness_shear, stiffness_torsion, stiffness_bending))
    return ElasticTube{Union{map(typeof, stiffnesses)...}}(stiffnesses...,
        SimFloat(damping), zeros(KVec3), Matrix{SimFloat}(I, 3, 3), 0, zeros(KVec3),
        zeros(KVec3))
end

"""
    TimoshenkoTube

A 2-node corotational Timoshenko beam element; a chain of tubes forms a beam. The
stiffness couples each node's transverse displacement to its rotation, so
transverse shear is represented.

Each rigidity (`EA` axial, `GA` shear with correction factor `shear_coeff`,
`Φ = 12·EI/(k·GA·L²)`, `GJ` torsion, `EIy`/`EIz` bending about the two transverse
axes) is either a `Real` (linear) or a callable of that mode's strain/curvature
returning the effective rigidity there. All callables must share one type.

$(TYPEDFIELDS)
"""
mutable struct TimoshenkoTube{S} <: AbstractTubeModel
    "Axial rigidity: `Real` EA [N], or callable EA(ε) of axial strain ε=δ/L₀."
    EA::S
    "Shear rigidity (before `shear_coeff`): `Real` GA [N], or callable GA(γ) of shear angle."
    GA::S
    "Torsional rigidity: `Real` GJ [N·m²], or callable GJ(κ) of twist rate κ=Δθ/L₀."
    GJ::S
    "Bending rigidity about y: `Real` EIy [N·m²], or callable EIy(κ) of curvature κ=Δθ/L₀."
    EIy::S
    "Bending rigidity about z: `Real` EIz [N·m²], or callable EIz(κ) of curvature κ=Δθ/L₀."
    EIz::S
    "Shear correction factor k (e.g. 5/6 solid, 8/9 inflated tube)."
    shear_coeff::SimFloat
    "Rayleigh damping β [s]: the element wrench is also evaluated at β times the
    deformation rate (C = βK), so every mode is damped in proportion to its
    rigidity and rigid motion stays undamped by construction."
    damping::SimFloat
    "Rest chord length [m], taken from the placed geometry at the first `place!`."
    rest_length::SimFloat
    "Rest orientation of node A in the element frame (set at reinit!)."
    const R_a_rel0::Matrix{SimFloat}
    "Rest orientation of node B in the element frame (set at reinit!)."
    const R_b_rel0::Matrix{SimFloat}
    "Index in `SystemDefinition.tubes` (assigned by SystemDefinition)."
    idx::Int64
    "Tube end on body A, as an offset from its origin in its frame [m]."
    const anchor_a_b::KVec3
    "Tube end on body B, as an offset from its origin in its frame [m]."
    const anchor_b_b::KVec3
end

"""
    TimoshenkoTube(; EA=nothing, GA=nothing, GJ=nothing, EIy=nothing, EIz=nothing,
                   shear_coeff=5/6, damping=0)

A [`TimoshenkoTube`](@ref). Each rigidity is a `Real` or a callable of its
strain/curvature; the ones left `nothing` are derived from the tube's law when the
[`Tube`](@ref) is built.
"""
TimoshenkoTube(; EA = nothing, GA = nothing, GJ = nothing, EIy = nothing,
               EIz = nothing, shear_coeff::Real = 5 / 6, damping::Real = 0.0) =
    timoshenko_tube((EA, GA, GJ, EIy, EIz), shear_coeff, damping)

function timoshenko_tube(rigidities, shear_coeff, damping)
    rigidities = concrete_rigidities("TimoshenkoTube", rigidities)
    return TimoshenkoTube{Union{map(typeof, rigidities)...}}(rigidities...,
        SimFloat(shear_coeff), SimFloat(damping), 0.0,
        Matrix{SimFloat}(I, 3, 3), Matrix{SimFloat}(I, 3, 3), 0, zeros(KVec3), zeros(KVec3))
end

"""
    rigidity_fields(model::AbstractTubeModel) -> Tuple

The fields of `model` holding a rigidity, each a `Real` or a callable.
"""
rigidity_fields(::ElasticTube) =
    (:stiffness_axial, :stiffness_shear, :stiffness_torsion, :stiffness_bending)
rigidity_fields(::TimoshenkoTube) = (:EA, :GA, :GJ, :EIy, :EIz)

"""
    concrete_rigidities(model, rigidities) -> Tuple

`rigidities` with every `Real` as a `SimFloat` and callables kept, erroring when
the callables are of more than one type.
"""
function concrete_rigidities(model, rigidities)
    converted = map(value -> value isa Real ? SimFloat(value) : value, rigidities)
    callable_types = unique(typeof(value) for value in converted
                            if !(value isa Union{Real, Nothing}))
    length(callable_types) > 1 && error(
        "$model: all callable rigidities must be the same type, got " *
        "$(callable_types). Mix only `Real`s and a single callable type.")
    return converted
end

"""
    with_law_rigidities(model, material, diameter, pressure) -> AbstractTubeModel

A fresh copy of `model` with every rigidity it leaves `nothing` derived from KiteGeometry
`material` at the tube's `diameter` [m] and `pressure` [Pa] ([`tube_law_rigidities`](@ref)).
"""
function with_law_rigidities(model::TimoshenkoTube, material, diameter, pressure)
    fields = rigidity_fields(model)
    rigidities = map(field -> getfield(model, field), fields)
    if any(isnothing, rigidities)
        EA, GA, EI, GJ = tube_law_rigidities(material, diameter, pressure)
        derived = (; EA, GA, GJ, EIy = EI, EIz = EI)
        rigidities = map((field, value) -> something(value, derived[field]), fields,
                         rigidities)
    end
    fresh = timoshenko_tube(rigidities, model.shear_coeff, model.damping)
    fresh.rest_length = model.rest_length
    return fresh
end

function with_law_rigidities(model::ElasticTube, material, diameter, pressure)
    stiffness_axial, stiffness_shear, stiffness_torsion, stiffness_bending =
        map(field -> getfield(model, field), rigidity_fields(model))
    return ElasticTube(; stiffness_axial, stiffness_shear, stiffness_torsion,
                       stiffness_bending, model.damping)
end

"""
    Tube{M<:AbstractTubeModel}

An inflated tube between two `Body`s: a KiteGeometry `tubes` row whose model is the
structural element it is simulated as, a [`TimoshenkoTube`](@ref) beam element or a lumped
[`ElasticTube`](@ref). Its ends sit at an anchor on each body, and the placed bodies fix its
rest length. Its `material` names its stiffness law.
"""
const Tube{M <: AbstractTubeModel} = KiteGeometry.Tube{M}


"""
    Tube(name, body_a, body_b; diameter, pressure, law=:breukels2011,
         model=TimoshenkoTube(), anchor_a=zeros, anchor_b=zeros)

A tube from `body_a` to `body_b` (names or indices), its ends at `anchor_a` and
`anchor_b` in each body's frame. The tube holds a copy of `model`, the rigidities it
leaves `nothing` derived from `law` at `diameter` [m] and `pressure` [Pa].
"""
function Tube(name, body_a, body_b; diameter::Real, pressure::Real,
              law = :breukels2011, model::AbstractTubeModel = TimoshenkoTube(),
              anchor_a = zeros(SimFloat, 3), anchor_b = zeros(SimFloat, 3))
    law = Symbol(law)
    material = KiteGeometry.Material(; name=law, law=String(law),
                                     parameters=KiteGeometry.tube_law_parameters(law))
    model = with_law_rigidities(model, material, diameter, pressure)
    model.anchor_a_b .= anchor_a
    model.anchor_b_b .= anchor_b
    return KiteGeometry.Tube(; name=Symbol(name), bodies=cell_ref.((body_a, body_b)),
        diameter=SimFloat(diameter), pressure=SimFloat(pressure), material=law, model)
end

"""
    timoshenko_element_frame(x_a, x_b, R_a) -> (e1, e2, e3, len)

Orthonormal corotational element frame: `e1` along the chord from node A to node
B, `e2`/`e3` from node A's y-axis projected transverse to the chord. `len` is the
current chord length. Generic over numeric and symbolic inputs.
"""
function timoshenko_element_frame(x_a, x_b, R_a)
    d = x_b .- x_a
    len = sqrt(sum(abs2, d))
    e1 = d ./ len
    y_a = R_a[:, 2]
    e2_raw = y_a .- sum(y_a .* e1) .* e1
    e2 = e2_raw ./ sqrt(sum(abs2, e2_raw))
    e3 = e1 × e2
    return e1, e2, e3, len
end

"""
    tube_initial_frames(tube, bodies) -> (R_a, R_b, anchor_a_w, anchor_b_w)

World rotations of the two connected bodies and the world positions of the tube's
two anchors, in the bodies' initial poses (`pos_ENU`, `Q_KA2ENU`).
"""
function tube_initial_frames(tube::Tube, bodies)
    body_a = bodies[tube.bodies[1]]
    body_b = bodies[tube.bodies[2]]
    R_a = initial_rotation(body_a)
    R_b = initial_rotation(body_b)
    anchor_a_w = origin_ENU(body_a) .+ R_a * tube.anchor_a_b
    anchor_b_w = origin_ENU(body_b) .+ R_b * tube.anchor_b_b
    return R_a, R_b, anchor_a_w, anchor_b_w
end

"""
    init_tube_rest!(tube, bodies)

Capture the rest reference of `tube.model` from the bodies' initial poses, so the
initial geometry is unstrained and the tube wrench there is exactly zero. One
`init_rest!` method per tube model:

- [`ElasticTube`](@ref): rest anchor offset (body-A frame) and rest relative
  rotation `R_a' R_b`.
- [`TimoshenkoTube`](@ref): rest length (on the first call) and per-node
  orientations relative to the corotational element frame.
"""
init_tube_rest!(tube::Tube, bodies) =
    init_rest!(tube.model, tube_initial_frames(tube, bodies)...)

function init_rest!(model::ElasticTube, R_a, R_b, anchor_a_w, anchor_b_w)
    model.rest_offset_a .= R_a' * (anchor_b_w .- anchor_a_w)
    model.R_rel0 .= R_a' * R_b
    return nothing
end

function init_rest!(model::TimoshenkoTube, R_a, R_b, x_a, x_b)
    e1, e2, e3, len = timoshenko_element_frame(x_a, x_b, R_a)
    element_frame = [e1[1] e2[1] e3[1];
                     e1[2] e2[2] e3[2];
                     e1[3] e2[3] e3[3]]
    model.rest_length ≈ 0 && (model.rest_length = len)
    model.R_a_rel0 .= element_frame' * R_a
    model.R_b_rel0 .= element_frame' * R_b
    return nothing
end

"""
    beam_anchor(pos_ENU, tube, bodies) -> (beam_frac, beam_offset_b)

Where a point at `pos_ENU` [m] rides `tube`'s beam element in the bodies' initial poses:
its position `s ∈ [0, 1]` along the element and its offset off the centreline in the
rest element frame [m].
"""
function beam_anchor(pos_ENU, tube::Tube, bodies)
    R_a, _, node_a, node_b = tube_initial_frames(tube, bodies)
    e1, e2, e3, len = timoshenko_element_frame(node_a, node_b, R_a)
    rel = pos_ENU .- node_a
    s = clamp(dot(rel, e1) / len, 0.0, 1.0)
    perp = rel .- (s * len) .* e1
    return s, KVec3(dot(perp, e1), dot(perp, e2), dot(perp, e3))
end

"""
    has_body_flap(station) -> Bool

`true` if the station reads its deflection off two ordered flap bodies.
"""
has_body_flap(station) = length(something(station.flap_bodies, ())) == 2

"""
    has_point_flap(station) -> Bool

`true` if the station reads its deflection off three structural points
`[fore, hinge, aft]` ([`point_flap_delta`](@ref)).
"""
has_point_flap(station) = length(something(station.flap_points, ())) == 3

"""
    has_flap(station) -> Bool

`true` if the station carries a flap hinge either way, so it drives a live
deflection δ ([`flap_delta`](@ref) / [`point_flap_delta`](@ref) /
`station_delta_eqs!`). A station given both reads the points.
"""
has_flap(station) =
    has_point_flap(station) || has_body_flap(station)

"""
    flap_delta(station, R_main, R_flap) -> SimFloat

Signed flap deflection δ [rad] of `station` from the two flap bodies'
world orientations `R_main`, `R_flap`: project each body's reference chord
direction onto the plane normal to the world hinge axis and take the signed angle
about that axis, minus the rest deflection `flap_rest_delta`.
"""
flap_delta(station, R_main, R_flap) =
    flap_delta_expression(station, R_main, R_flap)

"""
    point_flap_delta(station, points, R_wing) -> SimFloat

Deflection δ [rad] of a point flap from the current positions of its three points,
with `R_wing` the owning wing's orientation. The same
[`point_flap_delta_expression`](@ref) the symbolic RHS evaluates, so it is the ground
truth for tests.
"""
function point_flap_delta(station, points, R_wing)
    fore, hinge, aft = station.flap_points
    return point_flap_delta_expression(station, points[fore].pos_w,
        points[hinge].pos_w, points[aft].pos_w, R_wing)
end

"""
    init_station_flap!(station, sys_struct)

Capture a flapped `KINEMATIC` station's rest geometry from the initial pose, so the
undeformed configuration reads δ = 0. A body flap defaults its reference chords to
each body's x-axis and measures the two bodies' `Q_KA2ENU`; a point flap measures
its three points' `pos_ENU`. No-op for surfaces without a flap.
"""
function init_station_flap!(station, sys_struct)
    if has_point_flap(station)
        points = sys_struct.points
        fore, hinge, aft = station.flap_points
        station.flap_rest_delta = 0.0
        station.flap_rest_delta = point_flap_delta_expression(station,
            points[fore].pos_ENU, points[hinge].pos_ENU, points[aft].pos_ENU,
            initial_rotation(sys_struct.wings[station.wing]))
        return nothing
    end
    has_body_flap(station) || return nothing
    bodies = sys_struct.bodies
    isempty(station.flap_chord_refs) &&
        (station.flap_chord_refs = [KVec3(1.0, 0.0, 0.0), KVec3(1.0, 0.0, 0.0)])
    R_main = initial_rotation(bodies[station.flap_bodies[1]])
    R_flap = initial_rotation(bodies[station.flap_bodies[2]])
    station.flap_rest_delta = 0.0
    station.flap_rest_delta = flap_delta(station, R_main, R_flap)
    return nothing
end

