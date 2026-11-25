# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

"""
    BodyModel

What SymbolicAWEModels holds of a body beyond its KiteGeometry `bodies` row: the columns
only it reads, each `nothing` where the row leaves it unset.

A body row is a 6-DOF body integrated (or fitted) in the principal frame. The body of a
wing is the one its `wings` row names. Its dynamics type is DYNAMIC (free 6-DOF),
KINEMATIC (pose fitted from its points, the body of a particle wing) or STATIC (held by its
`fixed` column).

The body's own mass properties are its `extra_mass` and `extra_inertia_KA`. A DYNAMIC or
STATIC body adds the points it carries, its `BODY_STATIC` riders and a `RIGID_DYNAMICS`
wing's nodes, as point masses of their `total_mass` at their anchors: a run's
`total_mass`, `total_com_offset_b` and principal inertia describe the body with them
([`body_mass_properties`](@ref)). The row's `pos_ENU` is the centre of the body's own
mass; the origin of its frame is a run's `origin_ENU`.

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

"""The centre of `body`'s own mass from the origin of its frame, body frame [m]: its
row's `com_offset_b` column."""
own_com_offset(body::KiteGeometry.Body) =
    something(body.model.com_offset_b, zero(SVector{3, Float64}))

"""Put the origin of body row `body`'s frame at `origin` [m] in the world, keeping its
rotation."""
set_origin_ENU!(body, origin) = body.pos_ENU =
    origin + quaternion_to_rotation_matrix(body.Q_KA2ENU) * own_com_offset(body)

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
    initial_rotation(body) -> Matrix

The body→world rotation matrix of `body`'s initial pose `Q_KA2ENU`.
"""
initial_rotation(body) = quaternion_to_rotation_matrix(body.Q_KA2ENU)

"""
    initial_chord_position(sam, b, pos_ENU) -> SimFloat

Body-frame x of the position `pos_ENU` on body `b` of `sam` in the initial pose:
smaller towards the leading edge.
"""
initial_chord_position(sam, b, pos_ENU) = initial_body_position(sam, b, pos_ENU)[1]

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
    flap_delta(sam, s, R_main, R_flap) -> SimFloat

Signed flap deflection δ [rad] of station `s` of `sam` from the two flap bodies'
world orientations `R_main`, `R_flap` ([`flap_delta_expression`](@ref)).
"""
flap_delta(sam, s, R_main, R_flap) = flap_delta_expression(sam, s, R_main, R_flap)

"""
    point_flap_delta(sam, s, R_wing) -> SimFloat

Deflection δ [rad] of the point flap of station `s` of `sam` from the current positions
of its three points, with `R_wing` the owning wing's orientation. The same
[`point_flap_delta_expression`](@ref) the symbolic RHS evaluates, so it is the ground
truth for tests.
"""
function point_flap_delta(sam, s, R_wing)
    fore, hinge, aft = sam.definition.spanwise_stations[s].flap_points
    position(i) = run_value(sam, Val(:points), i, Val(:pos_w))
    return point_flap_delta_expression(sam, s, position(fore), position(hinge),
                                       position(aft), R_wing)
end
