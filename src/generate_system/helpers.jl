# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Helper functions for symbolic equation generation

"""
    WindFactor(am, profile_law)

Callable wind-shear factor, used as a callable flat parameter `w(pos_z)`: the
ratio of wind speed at height `pos_z` to the ground value, from atmospheric model
`am` under `profile_law` (1.0 when `profile_law == 0`). `ForwardDiff.Dual`-safe in
`pos_z`. Read live from `sys_struct` via [`WindFactorReader`](@ref).
"""
struct WindFactor
    am::AtmosphericModel
    profile_law::Int
end
function (w::WindFactor)(pos_z)
    w.profile_law == 0 && return 1.0
    return AtmosphericModels.calc_wind_factor(w.am, max(1.0, pos_z), w.profile_law)
end

"""
    WindFactorReader()

Serializable flat-param reader producing a [`WindFactor`](@ref) from the live
`sys_struct`'s atmospheric model and wind profile law.
"""
struct WindFactorReader end
(::WindFactorReader)(sys_struct) =
    WindFactor(sys_struct.am, sys_struct.set.profile_law)

"""
    calc_angle_of_attack(va_wing_b)

Calculate the angle of attack [rad] from the apparent wind vector `va_wing_b` in the
body frame.
"""
function calc_angle_of_attack(va_wing_b)
    return atan(va_wing_b[3], va_wing_b[1])
end

"""
    smooth_normalize(vec)

Differentiable normalization: `vec / smooth_norm(vec)`.
"""
smooth_normalize(vec) = vec / smooth_norm(vec)

"""
    smooth_norm(v, eps=1e-12)

Differentiable norm: `sqrt(sum(abs2, v) + eps^2)`.
"""
smooth_norm(v, eps=1e-12) = sqrt(sum(abs2, v) + eps^2)

"""
    quaternion_to_rotation_matrix(q)

Convert a quaternion `q` (scalar-first format [w, x, y, z]) to a 3x3 rotation
matrix.
"""
function quaternion_to_rotation_matrix(q::AbstractVector)
    w, x, y, z = q[1], q[2], q[3], q[4]

    return [
        1-2*(y*y+z*z) 2*(x*y-z*w) 2*(x*z+y*w)
        2*(x*y+z*w) 1-2*(x*x+z*z) 2*(y*z-x*w)
        2*(x*z-y*w) 2*(y*z+x*w) 1-2*(x*x+y*y)
    ]
end

"""
    rotation_matrix_to_quaternion(R)

Convert a 3x3 rotation matrix `R` to a quaternion (scalar-first format [w, x, y, z]),
selecting the numerically stable branch (largest of the trace / diagonal). Written
with `ifelse` and `&` (not `if`/`&&`) so it evaluates symbolically on `Num` as well
as on concrete `Float64` — a single symbolic expression, computed once and
common-subexpression-shared across components, with no `@register_symbolic`. The
unselected branches' `sqrt` arguments are clamped to ≥0 so they never error; the
selected branch always has a positive radicand.
"""
function rotation_matrix_to_quaternion(R::AbstractMatrix)
    Rc = collect(R)
    r11, r22, r33 = Rc[1, 1], Rc[2, 2], Rc[3, 3]
    tr_ = r11 + r22 + r33

    S1 = sqrt(max(tr_ + 1.0, 0.0)) * 2
    w1 = 0.25 * S1
    x1 = (Rc[3, 2] - Rc[2, 3]) / S1
    y1 = (Rc[1, 3] - Rc[3, 1]) / S1
    z1 = (Rc[2, 1] - Rc[1, 2]) / S1

    S2 = sqrt(max(1.0 + r11 - r22 - r33, 0.0)) * 2
    w2 = (Rc[3, 2] - Rc[2, 3]) / S2
    x2 = 0.25 * S2
    y2 = (Rc[1, 2] + Rc[2, 1]) / S2
    z2 = (Rc[1, 3] + Rc[3, 1]) / S2

    S3 = sqrt(max(1.0 + r22 - r11 - r33, 0.0)) * 2
    w3 = (Rc[1, 3] - Rc[3, 1]) / S3
    x3 = (Rc[1, 2] + Rc[2, 1]) / S3
    y3 = 0.25 * S3
    z3 = (Rc[2, 3] + Rc[3, 2]) / S3

    S4 = sqrt(max(1.0 + r33 - r11 - r22, 0.0)) * 2
    w4 = (Rc[2, 1] - Rc[1, 2]) / S4
    x4 = (Rc[1, 3] + Rc[3, 1]) / S4
    y4 = (Rc[2, 3] + Rc[3, 2]) / S4
    z4 = 0.25 * S4

    # Trace branch if positive, else the largest-diagonal branch (nested ifelse, no `&`).
    pick(b1, b2, b3, b4) = ifelse(tr_ > 0, b1,
        ifelse(r11 >= r22,
            ifelse(r11 >= r33, b2, b4),
            ifelse(r22 >= r33, b3, b4)))

    return [pick(w1, w2, w3, w4), pick(x1, x2, x3, x4),
            pick(y1, y2, y3, y4), pick(z1, z2, z3, z4)]
end

"""
    rotate_v_around_k(v, k, θ)

Rotate vector `v` around axis `k` by angle `θ` using Rodrigues' rotation formula.
"""
function rotate_v_around_k(v, k, θ)
    k = smooth_normalize(k)
    v_rot = v * cos(θ) + (k × v) * sin(θ) + k * (k ⋅ v) * (1 - cos(θ))
    return v_rot
end

"""
    calc_R_v_to_w(wing_pos, e_x)

Calculate the rotation matrix from the view frame (`_v`) to the world frame (`_w`).

The view frame is defined with its z-axis pointing from the origin to the wing,
and its x-axis aligned with the wing's x-axis projected onto the view plane.

"""
function calc_R_v_to_w(wing_pos, e_x)
    wp1, wp2, wp3 = wing_pos[1], wing_pos[2], wing_pos[3]
    ex1, ex2, ex3 = e_x[1], e_x[2], e_x[3]

    # z = normalize(wing_pos)
    wp_norm = smooth_norm((wp1, wp2, wp3))
    z1, z2, z3 = wp1 / wp_norm, wp2 / wp_norm, wp3 / wp_norm

    # y = normalize(z × e_x)
    zxe1 = z2 * ex3 - z3 * ex2
    zxe2 = z3 * ex1 - z1 * ex3
    zxe3 = z1 * ex2 - z2 * ex1
    zxe_norm = smooth_norm((zxe1, zxe2, zxe3))
    y1, y2, y3 = zxe1 / zxe_norm, zxe2 / zxe_norm, zxe3 / zxe_norm

    # x = y × z
    x1 = y2 * z3 - y3 * z2
    x2 = y3 * z1 - y1 * z3
    x3 = y1 * z2 - y2 * z1

    return [x1 y1 z1; x2 y2 z2; x3 y3 z3]
end

"""
    calc_R_t_to_w(wing_pos)

Calculate the rotation matrix from the local tether frame (`_t`) to the world
frame (`_w`).

The tether frame is a local spherical coordinate system:
- **z-axis**: Aligned with the tether (radial direction).
- **y-axis**: Azimuthal direction, parallel to the XY plane.
- **x-axis**: Elevation direction, tangent to the sphere (`y × z`).
"""
function calc_R_t_to_w(wing_pos)
    z = smooth_normalize(wing_pos)
    if wing_pos[2] ≈ 0.0 && wing_pos[1] ≈ 0.0
        y = [0, 1, 0]
    else
        y = smooth_normalize([-wing_pos[2], wing_pos[1], 0])
    end
    x = y × z
    return [x[1] y[1] z[1]; x[2] y[2] z[2]; x[3] y[3] z[3]]
end

"""
    sym_calc_R_t_to_w(wing_pos)

Symbolic version of `calc_R_t_to_w` that uses explicit element access
to avoid slice scalarization issues.
"""
function sym_calc_R_t_to_w(wing_pos)
    wp1, wp2, wp3 = wing_pos[1], wing_pos[2], wing_pos[3]

    # z = normalize(wing_pos)
    wp_norm = smooth_norm((wp1, wp2, wp3))
    z1, z2, z3 = wp1 / wp_norm, wp2 / wp_norm, wp3 / wp_norm

    # y = normalize([-wp2, wp1, 0])
    yu1, yu2, yu3 = -wp2, wp1, 0
    y_norm = smooth_norm((yu1, yu2, yu3))
    y1, y2, y3 = yu1 / y_norm, yu2 / y_norm, yu3 / y_norm

    # x = y × z
    x1 = y2 * z3 - y3 * z2
    x2 = y3 * z1 - y1 * z3
    x3 = y1 * z2 - y2 * z1

    return [x1 y1 z1; x2 y2 z2; x3 y3 z3]
end

"""
    Base.getindex(x::ModelingToolkit.Symbolics.Arr, idxs::Vector{Int64})

Extend `Base.getindex` to allow indexing a symbolic array with a vector of
integer indices, which is not natively supported by ModelingToolkit.
"""
function Base.getindex(x::ModelingToolkit.Symbolics.Arr, idxs::Vector{Int64})
    Num[x[idx] for idx in idxs]
end
