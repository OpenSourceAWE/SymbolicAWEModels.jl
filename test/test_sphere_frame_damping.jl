# SPDX-FileCopyrightText: 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_sphere_frame_damping.jl
#
# Angular and body-frame damping act on motion relative to the frame that follows a
# wing over the sphere about its tether's base point, so steady flight across that
# sphere is not braked. Tested on `body_integration` directly, then on the
# right-hand side of both backends: the damping a state gets is the rate with the
# damping on minus the rate with it off, at the same state.

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

using Test
using LinearAlgebra
using SymbolicAWEModels
using SymbolicAWEModels: KernelBackend, MonolithBackend, sync_params!, sync_initial!,
                         body_integration, sphere_frame, VortexStepMethod
using KiteUtils

"""Rotation by `angle` [rad] about unit `axis`."""
rotation(axis, angle) = cos(angle) * I + sin(angle) * [0 -axis[3] axis[2];
    axis[3] 0 -axis[1]; -axis[2] axis[1] 0] + (1 - cos(angle)) * axis * axis'

"""
    flying_two_plate(backend, geometry, root, elevation_vel) -> SymbolicAWEModel

The 2plate kite built from `geometry` on `backend`, initialised while its transform
turns it over the sphere about the ground point at `elevation_vel` [rad/s].
"""
function flying_two_plate(backend, geometry, root, elevation_vel)
    data_path = joinpath(root, string(nameof(typeof(backend))), geometry, "2plate_kite")
    mkpath(dirname(data_path))
    cp(joinpath(dirname(@__DIR__), "data", "2plate_kite"), data_path; force = true)
    set_data_path(data_path)
    set = Settings("system.yaml")
    vsm_set = VortexStepMethod.VSMSettings(
        joinpath(data_path, "vsm_settings.yaml"); data_prefix = false)
    sys, transforms = load_sys_struct_from_yaml(
        joinpath(data_path, geometry);
        system_name = "sphere_frame_damping", set, vsm_set)
    transforms[1].elevation_vel = elevation_vel
    sam = SymbolicAWEModel(set, sys; transforms, backend)
    init!(sam; prn = false, remake = true)
    return sam
end

"""Rate of the model's state at its initial state after `change!(sam.state)`."""
function state_rate(sam, change! = identity)
    change!(sam.state)
    sync_initial!(sam.prob.initial_sync, sam.prob.prob, sam.state)
    u = copy(sam.prob.prob.u0)
    du = similar(u)
    sam.prob.prob.f(du, u, sam.integrator.p, sam.integrator.t)
    return du
end

"""
    damping_rate(sam, change! = identity)

What the angular and body-frame damping add to the rate of the model's initial state
after `change!`: the rate with every coefficient as built minus the rate with all of
them zero.
"""
function damping_rate(sam, change! = identity)
    state = sam.state
    angular = [copy(body.angular_damping) for body in state.bodies]
    body_frame = [copy(point.body_frame_damping) for point in state.points]
    damped = state_rate(sam, change!)
    set_angular_damping(state.bodies, 0.0)
    set_body_frame_damping(state.points, 0.0)
    sync_params!(sam.prob.param_sync, sam.integrator, state)
    undamped = state_rate(sam, change!)
    foreach((body, c) -> body.angular_damping = c, state.bodies, angular)
    foreach((point, c) -> point.body_frame_damping = c, state.points, body_frame)
    sync_params!(sam.prob.param_sync, sam.integrator, state)
    return damped .- undamped
end

@testset "Damping against the frame that follows the wing over the sphere" begin
    base = [0.0, 0.0, -20.0]
    elevation = deg2rad(50)
    com = base .+ 200 .* [cos(elevation), 0.0, sin(elevation)]
    sphere_spin = [0.0, -0.15, 0.0]
    com_vel = sphere_spin × (com .- base)
    orientation_p = rotation(normalize([1.0, 2.0, 3.0]), 0.7)
    coefficients = [20.0, 150.0, 5.0]
    params = (; bodies = [(; fix_sphere = false, fix_static = false,
        angular_damping = coefficients, world_frame_damping = zeros(3),
        body_frame_damping = zeros(3), R_b_to_p = Matrix(1.0I, 3, 3))])
    spin_rate(omega_p) = body_integration(params, 1, com, com_vel, omega_p,
        zeros(3), zeros(3), orientation_p;
        frame = sphere_frame(base, com, com_vel)).d_ω_p

    @testset "the frame turns at r × v / |r|² about the base point" begin
        @test sphere_frame(base, com, com_vel).spin ≈ sphere_spin
        @test !(sphere_frame(zeros(3), com, com_vel).spin ≈ sphere_spin)
        @test sphere_frame(nothing, com, com_vel).spin == zeros(3)
    end

    @testset "a body spinning with the sphere feels no damping torque" begin
        @test norm(spin_rate(orientation_p' * sphere_spin)) < 1e-12
    end

    @testset "spin relative to the sphere is damped per body axis" begin
        extra = [0.1, -0.3, 0.2]
        @test spin_rate(orientation_p' * sphere_spin .+ extra) ≈ -coefficients .* extra
    end

    data_path_before = get_data_path()
    root = mktempdir()
    @testset "$(nameof(typeof(backend))) leaves rigid flight over the sphere undamped" for
            backend in (KernelBackend(), MonolithBackend())
        sam = flying_two_plate(backend, "rigid_structural_geometry.yaml", root, 0.4)
        wing = sam.state.wings[1]
        @test norm(wing.ω_b) > 0.3
        @test maximum(abs, damping_rate(sam)) < 1e-9
        flown = copy(wing.ω_b)
        pitched(state) = state.wings[1].ω_b .= flown .+ [0.0, 0.5, 0.0]
        @test maximum(abs, damping_rate(sam, pitched)) > 1.0
    end
    @testset "$(nameof(typeof(backend))) leaves a fitted wing's nodes undamped in rigid \
              flight over the sphere" for backend in (KernelBackend(), MonolithBackend())
        sam = flying_two_plate(backend, "particle_structural_geometry.yaml", root, 0.4)
        @test maximum(norm, point.vel_w for point in sam.state.points) > 1.0
        @test maximum(abs, damping_rate(sam)) < 1e-9
    end
    set_data_path(data_path_before)
end
