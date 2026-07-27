# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

"""
    ObjAdapter

Mesh mass-property helpers for wings authored from a triangulated `.obj`: center
of mass, surface inertia tensor, and the per-unit-mass `(com, inertia)` used to
seed a wing's mass properties. Mesh IO is delegated to
`VortexStepMethod.ObjAdapter.read_faces`; this module owns only the mass-property
maths, kept out of VortexStepMethod (which is aero-only).
"""
module ObjAdapter

using LinearAlgebra
import VortexStepMethod

export center_of_mass, calculate_inertia_tensor, unit_inertia_from_obj,
    unit_inertia_matrix, unit_inertia_vector

"""
    center_of_mass(vertices, faces) -> Vector{Float64}

Area-weighted center of mass of a triangulated surface mesh, in the mesh's own
CAD frame. Non-mutating; requires triangular faces. Unlike VortexStepMethod's old
`center_to_com!`, this neither shifts the vertices nor forces the COM onto the
xz-plane.
"""
function center_of_mass(vertices, faces)
    area_total = 0.0
    com = zeros(3)
    for face in faces
        length(face) == 3 ||
            throw(ArgumentError("Triangulate faces in a CAD program first"))
        v1 = vertices[face[1]]; v2 = vertices[face[2]]; v3 = vertices[face[3]]
        area = norm(cross(v2 - v1, v3 - v1)) / 2
        area_total += area
        com += area * (v1 + v2 + v3) / 3
    end
    return com / area_total
end

"""
    calculate_inertia_tensor(vertices, faces, mass, com) -> Matrix{Float64}

Surface-area-weighted 3×3 inertia tensor about `com`, for total `mass` spread
uniformly over the mesh surface. Pass `mass = 1` for the per-unit-mass tensor
[m²]; multiply by the physical mass for [kg·m²].
"""
function calculate_inertia_tensor(vertices, faces, mass, com)
    inertia = zeros(3, 3)
    total_area = 0.0
    for face in faces
        v1 = vertices[face[1]] .- com
        v2 = vertices[face[2]] .- com
        v3 = vertices[face[3]] .- com
        area = norm(cross(v2 - v1, v3 - v1)) / 2
        total_area += area
        for p in (v1, v2, v3), i in 1:3, j in 1:3
            i == j ? (inertia[i, i] += area * (sum(p .^ 2) - p[i]^2)) :
                (inertia[i, j] -= area * (p[i] * p[j]))
        end
    end
    return (mass / total_area) * inertia / 3
end

"""
    unit_inertia_from_obj(obj_path) -> (com, unit_inertia)

Read a triangulated `.obj` and return its center of mass `com` [m] and its
per-unit-mass inertia tensor `unit_inertia` (3×3, [m²], about `com`), both in the
mesh CAD frame. Multiply `unit_inertia` by the wing mass for the physical tensor.
"""
function unit_inertia_from_obj(obj_path)
    vertices, faces = VortexStepMethod.ObjAdapter.read_faces(obj_path)
    com = center_of_mass(vertices, faces)
    return com, calculate_inertia_tensor(vertices, faces, 1.0, com)
end

"""
    unit_inertia_matrix(unit_inertia) -> Matrix{Float64}

Rebuild the symmetric 3×3 inertia tensor from the 6-vector
`[Ixx, Iyy, Izz, Ixy, Ixz, Iyz]`.
"""
unit_inertia_matrix(unit_inertia) = [
    unit_inertia[1] unit_inertia[4] unit_inertia[5]
    unit_inertia[4] unit_inertia[2] unit_inertia[6]
    unit_inertia[5] unit_inertia[6] unit_inertia[3]
]

"""
    unit_inertia_vector(tensor) -> Vector{Float64}

Pack a symmetric 3×3 inertia tensor into `[Ixx, Iyy, Izz, Ixy, Ixz, Iyz]`.
"""
unit_inertia_vector(tensor) = [
    tensor[1, 1], tensor[2, 2], tensor[3, 3],
    tensor[1, 2], tensor[1, 3], tensor[2, 3],
]

end
