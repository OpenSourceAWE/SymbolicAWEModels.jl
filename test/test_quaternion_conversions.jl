# SPDX-FileCopyrightText: 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

using Test
using SymbolicAWEModels
using LinearAlgebra
using Rotations

@testset verbose=true "Quaternion conversion functions" begin

    @testset "quaternion_to_rotation_matrix" begin
        # Test identity quaternion
        q_identity = [1.0, 0.0, 0.0, 0.0]  # w, x, y, z
        R = SymbolicAWEModels.quaternion_to_rotation_matrix(q_identity)
        @test R ≈ I(3) atol=1e-10

        # Test 90° rotation around Z axis
        # Quaternion for 90° around Z: [cos(45°), 0, 0, sin(45°)]
        q_z90 = [cos(π/4), 0.0, 0.0, sin(π/4)]
        R_z90 = SymbolicAWEModels.quaternion_to_rotation_matrix(q_z90)
        expected_z90 = [0.0 -1.0 0.0; 1.0 0.0 0.0; 0.0 0.0 1.0]
        @test R_z90 ≈ expected_z90 atol=1e-10

        # Test arbitrary quaternion
        q_arb = [0.9239, 0.3827, 0.0, 0.0]  # Normalize first
        q_arb = q_arb / norm(q_arb)
        R_arb = SymbolicAWEModels.quaternion_to_rotation_matrix(q_arb)

        # Check orthogonality
        @test R_arb' * R_arb ≈ I(3) atol=1e-10
        @test det(R_arb) ≈ 1.0 atol=1e-10
    end

    @testset "rotation_matrix_to_quaternion" begin
        # Test identity matrix
        R_identity = Matrix{Float64}(I(3))
        q = SymbolicAWEModels.rotation_matrix_to_quaternion(R_identity)
        @test q[1] ≈ 1.0 atol=1e-10  # w component
        @test norm(q[2:4]) ≈ 0.0 atol=1e-10  # x, y, z components

        # Test 90° rotation around Z
        R_z90 = [0.0 -1.0 0.0; 1.0 0.0 0.0; 0.0 0.0 1.0]
        q_z90 = SymbolicAWEModels.rotation_matrix_to_quaternion(R_z90)
        expected_q = [cos(π/4), 0.0, 0.0, sin(π/4)]
        @test q_z90 ≈ expected_q atol=1e-10

        # Test 90° rotation around X
        R_x90 = [1.0 0.0 0.0; 0.0 0.0 -1.0; 0.0 1.0 0.0]
        q_x90 = SymbolicAWEModels.rotation_matrix_to_quaternion(R_x90)
        expected_q_x = [cos(π/4), sin(π/4), 0.0, 0.0]
        @test q_x90 ≈ expected_q_x atol=1e-10

        # Test 90° rotation around Y
        R_y90 = [0.0 0.0 1.0; 0.0 1.0 0.0; -1.0 0.0 0.0]
        q_y90 = SymbolicAWEModels.rotation_matrix_to_quaternion(R_y90)
        expected_q_y = [cos(π/4), 0.0, sin(π/4), 0.0]
        @test q_y90 ≈ expected_q_y atol=1e-10

        # Test quaternion normalization
        q_normalized = q_z90 / norm(q_z90)
        @test norm(q_normalized) ≈ 1.0 atol=1e-10
    end

    @testset "Round-trip conversion" begin
        # Test multiple random rotations
        for _ in 1:20
            # Generate random quaternion
            q_rand = randn(4)
            q_rand = q_rand / norm(q_rand)

            # Convert to rotation matrix and back
            R = SymbolicAWEModels.quaternion_to_rotation_matrix(q_rand)
            q_recovered = SymbolicAWEModels.rotation_matrix_to_quaternion(R)

            # Quaternions q and -q represent the same rotation
            # So check if q_recovered ≈ q_rand OR q_recovered ≈ -q_rand
            match_pos = norm(q_recovered - q_rand) < 1e-10
            match_neg = norm(q_recovered + q_rand) < 1e-10
            @test match_pos || match_neg
        end

        # Test round-trip starting from rotation matrix
        for _ in 1:20
            # Generate random rotation matrix using Rotations.jl
            R_rand = rand(RotMatrix{3})

            # Convert to quaternion and back
            q = SymbolicAWEModels.rotation_matrix_to_quaternion(Matrix(R_rand))
            R_recovered = SymbolicAWEModels.quaternion_to_rotation_matrix(q)

            @test R_recovered ≈ R_rand atol=1e-10
        end
    end

    @testset "Edge cases" begin
        # Test rotation matrix with negative trace (edge case branch)
        # 180° rotation around X axis has trace = -1
        R_x180 = [1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 -1.0]
        q_x180 = SymbolicAWEModels.rotation_matrix_to_quaternion(R_x180)
        R_recovered = SymbolicAWEModels.quaternion_to_rotation_matrix(q_x180)
        @test R_recovered ≈ R_x180 atol=1e-10

        # Test rotation matrix with R[2,2] > R[3,3] (different branch)
        # Rotation that emphasizes Y component
        angle = π/3
        R_y = [cos(angle) 0.0 sin(angle); 0.0 1.0 0.0; -sin(angle) 0.0 cos(angle)]
        q_y = SymbolicAWEModels.rotation_matrix_to_quaternion(R_y)
        R_recovered_y = SymbolicAWEModels.quaternion_to_rotation_matrix(q_y)
        @test R_recovered_y ≈ R_y atol=1e-10

        # Test rotation matrix with R[3,3] dominant (final branch)
        angle = π/6
        R_z = [cos(angle) -sin(angle) 0.0; sin(angle) cos(angle) 0.0; 0.0 0.0 1.0]
        q_z = SymbolicAWEModels.rotation_matrix_to_quaternion(R_z)
        R_recovered_z = SymbolicAWEModels.quaternion_to_rotation_matrix(q_z)
        @test R_recovered_z ≈ R_z atol=1e-10
    end

    @testset "Consistency with Rotations.jl" begin
        # Compare results with Rotations.jl package
        for _ in 1:20
            # Generate random quaternion using Rotations.jl
            q_rot = rand(QuatRotation)
            q_array = [q_rot.w, q_rot.x, q_rot.y, q_rot.z]

            # Convert to rotation matrix using both methods
            R_rot = RotMatrix(q_rot)
            R_ours = SymbolicAWEModels.quaternion_to_rotation_matrix(q_array)

            @test R_ours ≈ Matrix(R_rot) atol=1e-10

            # Test reverse conversion
            q_recovered = SymbolicAWEModels.rotation_matrix_to_quaternion(Matrix(R_rot))
            # Account for quaternion sign ambiguity
            match_pos = norm(q_recovered - q_array) < 1e-10
            match_neg = norm(q_recovered + q_array) < 1e-10
            @test match_pos || match_neg
        end
    end
end
nothing
