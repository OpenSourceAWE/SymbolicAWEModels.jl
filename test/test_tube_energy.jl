# SPDX-FileCopyrightText: 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_tube_energy.jl - A tube's restoring wrench is minus the gradient of its strain
# energy, for both tube models, at finite rotations. The gradient is taken by
# automatic differentiation of the energy, so the wrench formula is checked against
# its definition rather than against a second formula.

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

using Test
using SymbolicAWEModels
using LinearAlgebra
using SymbolicAWEModels: ForwardDiff, cross_matrix, element_frame_matrix,
    timoshenko_node_rotations, timoshenko_strain_energy, timoshenko_nodal_wrench,
    elastic_tube_deformation, elastic_tube_strain_energy, elastic_tube_kinematic_map,
    tube_stiffness_term

"""A tube of `model` between two bodies placed apart, its rest geometry captured."""
function placed_tube(model)
    inertia = [0.01, 0.01, 0.01]
    body_a = Body(:a; extra_mass=1.0, inertia_principal=inertia, pos=[0.0, 0.0, 0.0])
    body_b = Body(:b; extra_mass=1.0, inertia_principal=inertia, pos=[1.0, 0.2, 0.1])
    tube = Tube(:tube, :a, :b; diameter=0.15, pressure=3e4, model)
    tube.body_a_idx = 1
    tube.body_b_idx = 2
    tube.model.idx = 1
    SymbolicAWEModels.init_tube_rest!(tube, [body_a, body_b])
    return tube
end

"""`R` turned by the spin `φ` to first order, which is all a gradient at `φ = 0` sees."""
turned(φ, R) = (I + cross_matrix(φ)) * R

# A deformed pose: each node displaced and turned by up to 0.3 rad from rest.
ANCHOR_A = [0.01, -0.02, 0.03]
ANCHOR_B = [1.05, 0.1, 0.25]
ROTATION_A = exp(cross_matrix([0.2, -0.1, 0.3]))
ROTATION_B = exp(cross_matrix([-0.3, 0.25, 0.1]))

"""Gradient of `energy(x_a, R_a, x_b, R_b)` in the anchors' `(x_a, φ_a, x_b, φ_b)`."""
pose_gradient(energy) = ForwardDiff.gradient(motion -> energy(
    ANCHOR_A .+ motion[1:3], turned(motion[4:6], ROTATION_A),
    ANCHOR_B .+ motion[7:9], turned(motion[10:12], ROTATION_B)), zeros(12))

@testset "Tube wrenches are minus their energy's gradient" begin
    @testset "TimoshenkoTube with Breukels laws" begin
        bending = tube_bending_law(0.075, 0.3)
        tube = placed_tube(TimoshenkoTube(; EA=1e4, GA=1500.0,
            GJ=tube_torsion_law(0.075, 0.3), EIy=bending, EIz=bending,
            shear_coeff=TUBE_SHEAR_COEFF))
        model = tube.model
        function energy(x_a, R_a, x_b, R_b)
            frame, len = element_frame_matrix(x_a, x_b, R_a)
            θ_a, θ_b = timoshenko_node_rotations(model, frame, R_a, R_b)
            return timoshenko_strain_energy(model, [len - model.rest_length; θ_a; θ_b])
        end
        frame, len = element_frame_matrix(ANCHOR_A, ANCHOR_B, ROTATION_A)
        θ_a, θ_b = timoshenko_node_rotations(model, frame, ROTATION_A, ROTATION_B)
        wrench = timoshenko_nodal_wrench(tube, (; tubes=[tube]), frame, len,
            ROTATION_A, ROTATION_B, θ_a, θ_b, zeros(12))
        gradient = pose_gradient(energy)
        @test maximum(abs, wrench .+ gradient) < 1e-12 * maximum(abs, gradient)
    end

    @testset "ElasticTube" begin
        tube = placed_tube(ElasticTube(; stiffness_axial=1e4, stiffness_shear=3e3,
            stiffness_torsion=50.0, stiffness_bending=80.0))
        model = tube.model
        energy(x_a, R_a, x_b, R_b) = elastic_tube_strain_energy(model,
            vcat(elastic_tube_deformation(model, x_a, R_a, x_b, R_b)...))
        deformation = vcat(elastic_tube_deformation(model, ANCHOR_A, ROTATION_A,
                                                    ANCHOR_B, ROTATION_B)...)
        generalized = [tube_stiffness_term(tube, (; tubes=[tube]), kind, deformation[i])
                       for (i, kind) in enumerate((1, 2, 2, 3, 4, 4))]
        wrench = -(elastic_tube_kinematic_map(ROTATION_A, model.R_rel0, ROTATION_B,
                                              ANCHOR_B .- ANCHOR_A)' * generalized)
        gradient = pose_gradient(energy)
        @test maximum(abs, wrench .+ gradient) < 1e-12 * maximum(abs, gradient)
    end
end
nothing
