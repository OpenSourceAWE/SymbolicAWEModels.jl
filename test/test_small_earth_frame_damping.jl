# SPDX-FileCopyrightText: 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_small_earth_frame_damping.jl
#
# Angular and body-frame damping act on motion relative to the kite's small-earth
# frame: z from the tether's base point to the kite, y the kite's span, stretching
# radially with the tether. Steady flight over the sphere, a turn about the tether
# and reel-out are not braked. Tested on the frame and `body_integration` directly,
# then on the right-hand side of both backends: the damping a state gets is the rate
# with the damping on minus the rate with it off, at the same state.

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

using Test
using LinearAlgebra
using SymbolicAWEModels
using SymbolicAWEModels: KernelBackend, MonolithBackend, sync_initial!, body_integration,
                         small_earth_frame, frame_velocity, damping_base,
                         damping_base_position, parameter_access, wing_body,
                         VortexStepMethod
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
    sys, transforms = load_authoring(joinpath(data_path, geometry);
                                     name = "small_earth_frame_damping", set)
    transforms[1].elevation_vel = elevation_vel
    sam = SymbolicAWEModel(set, sys; transforms, vsm_set, backend)
    init!(sam; prn = false, remake = true)
    return sam
end

"""Rate of the model's state at its initial state after `change!(sam)`."""
function state_rate(sam, change! = identity)
    change!(sam)
    sync_initial!(sam.prob.initial_sync, sam.prob.prob, sam)
    u = copy(sam.prob.prob.u0)
    du = similar(u)
    sam.prob.prob.f(du, u, sam.integrator.p, sam.integrator.t)
    return du
end

"""The `(collection, index, quantity)` of each angular and body-frame damping parameter
the model compiled."""
damping_parameters(sam) = filter(path -> !isnothing(parameter_access(sam, path)),
    [[(:bodies, i, :angular_damping) for i in eachindex(sam.definition.bodies)];
     [(:points, i, :body_frame_damping) for i in eachindex(sam.definition.points)]])

"""
    damping_rate(sam, change! = identity)

What the angular and body-frame damping add to the rate of the model's initial state
after `change!`: the rate with every coefficient as built minus the rate with all of
them zero.
"""
function damping_rate(sam, change! = identity)
    paths = damping_parameters(sam)
    built = [run_value(sam, path...) for path in paths]
    damped = state_rate(sam, change!)
    foreach(path -> set_run_value!(sam, path..., zeros(3)), paths)
    undamped = state_rate(sam, change!)
    foreach((path, c) -> set_run_value!(sam, path..., c), paths, built)
    return damped .- undamped
end

"""The small-earth axes `[x y z]` of a kite at `radius` from its base with span `span`."""
function small_earth_axes(radius, span)
    z = normalize(radius)
    x = normalize(span × z)
    return [x z × x z]
end

"""Position [m] of the base point the model's first wing damps against."""
tether_base(sam) = damping_base_position(damping_base(sam, wing_body(sam, 1)),
                                         point -> run_value(sam, :points, point, :pos_w))

"""Velocity [m/s] of each point of `sam`."""
point_velocities(sam) = [run_value(sam, :points, i, :vel_w)
                         for i in eachindex(sam.definition.points)]

"""Set the velocity of each point `i` of `sam` to `velocity(i, pos_w)`."""
set_point_velocities!(sam, velocity) = foreach(eachindex(sam.definition.points)) do i
    set_run_value!(sam, :points, i, :vel_w,
                   velocity(i, run_value(sam, :points, i, :pos_w)))
end

@testset "Damping against the kite's small-earth frame" begin
    base = [0.0, 0.0, -20.0]
    elevation = deg2rad(50)
    kite = base .+ 200 .* [cos(elevation), 0.0, sin(elevation)]
    kite_vel = [-3.0, 25.0, 4.0]
    kite_spin = [0.05, -0.12, 0.3]
    kite_axes = rotation(normalize([1.0, 2.0, 3.0]), 0.7)
    frame = small_earth_frame(base, kite, kite_vel, kite_axes[:, 2], kite_spin)
    coefficients = [20.0, 150.0, 5.0]

    @testset "the frame turns as its axes do and stretches with the tether" begin
        h = 1e-6
        axes_at(dt) = small_earth_axes(kite .+ kite_vel .* dt .- base,
            rotation(normalize(kite_spin), norm(kite_spin) * dt) * kite_axes[:, 2])
        turn = (axes_at(h) - axes_at(-h)) * axes_at(0)' / (2h)
        @test frame.spin ≈ [turn[3, 2], turn[1, 3], turn[2, 1]] rtol = 1e-7
        @test frame.stretch ≈ (kite .- base) ⋅ kite_vel / norm(kite .- base)^2
        @test frame_velocity(frame, kite) ≈ kite_vel
        @test frame_velocity(frame, base) ≈ zeros(3) atol = 1e-12
        @test !(small_earth_frame(zeros(3), kite, kite_vel, kite_axes[:, 2],
                                  kite_spin).spin ≈ frame.spin)
        still = small_earth_frame(nothing, kite, kite_vel, nothing, nothing)
        @test still.spin == zeros(3) && still.stretch == 0
    end

    # A body hanging below the kite, flying with it on the kite's frame.
    com = kite .+ [0.5, -1.0, -6.0]
    orientation_p = rotation(normalize([3.0, -1.0, 2.0]), 1.1)
    params = (; bodies = [(; fix_sphere = false, fix_static = false,
        angular_damping = coefficients, world_frame_damping = zeros(3),
        body_frame_damping = coefficients, R_b_to_p = Matrix(1.0I, 3, 3))])
    rates(com_vel, omega_p) = body_integration(params, 1, com, com_vel, omega_p,
        zeros(3), zeros(3), orientation_p; frame, wing_frame = kite_axes)
    flying = frame_velocity(frame, com)

    @testset "a body moving and spinning with the frame feels no damping" begin
        still = rates(flying, orientation_p' * frame.spin)
        @test norm(still.d_ω_p) < 1e-12
        @test norm(still.d_com_vel) < 1e-12
    end

    @testset "a turn about the tether is not damped" begin
        along = normalize(kite .- base)
        turning = small_earth_frame(base, kite, kite_vel, kite_axes[:, 2],
                                    kite_spin .+ 0.4 .* along)
        @test turning.spin ≈ frame.spin .+ 0.4 .* along
    end

    @testset "spin and motion relative to the frame are damped per axis" begin
        extra = [0.1, -0.3, 0.2]
        @test rates(flying, orientation_p' * frame.spin .+ extra).d_ω_p ≈
              -coefficients .* extra
        @test rates(flying .+ kite_axes * extra, orientation_p' * frame.spin).d_com_vel ≈
              -kite_axes * (coefficients .* extra)
    end

    data_path_before = get_data_path()
    root = mktempdir()
    @testset "$(nameof(typeof(backend))) leaves rigid flight on the small-earth frame \
              undamped" for backend in (KernelBackend(), MonolithBackend())
        sam = flying_two_plate(backend, "rigid_structural_geometry.yaml", root, 0.4)
        wing(quantity) = run_value(sam, :wings, 1, quantity)
        @test norm(wing(:ω_b)) > 0.3
        @test maximum(abs, damping_rate(sam)) < 1e-9
        flown = wing(:ω_b)
        base = tether_base(sam)
        along = normalize(wing(:com_w) .- base)
        flown_vel = point_velocities(sam)
        flown_wing_vel = wing(:vel_w)
        function turned(sam)
            turn = 0.3 .* along
            set_run_value!(sam, :wings, 1, :ω_b, flown .+ wing(:R_b_to_w)' * turn)
            set_run_value!(sam, :wings, 1, :vel_w,
                           flown_wing_vel .+ turn × (wing(:pos_w) .- base))
            set_point_velocities!(sam, (i, pos) -> flown_vel[i] .+ turn × (pos .- base))
        end
        @test maximum(abs, damping_rate(sam, turned)) < 1e-9
        pitched(sam) = set_run_value!(sam, :wings, 1, :ω_b, flown .+ [0.0, 0.5, 0.0])
        @test maximum(abs, damping_rate(sam, pitched)) > 1.0
    end
    @testset "$(nameof(typeof(backend))) leaves a fitted wing's nodes undamped in rigid \
              flight on the small-earth frame" for
            backend in (KernelBackend(), MonolithBackend())
        sam = flying_two_plate(backend, "particle_structural_geometry.yaml", root, 0.4)
        @test maximum(norm, point_velocities(sam)) > 1.0
        @test maximum(abs, damping_rate(sam)) < 1e-9
        base = tether_base(sam)
        along = normalize(run_value(sam, :wings, 1, :pos_w) .- base)
        flown = point_velocities(sam)
        flown_wing = run_value(sam, :wings, 1, :vel_w)
        field(spin, stretch, pos) = spin × (pos .- base) .+ stretch .* (pos .- base)
        rigidly(spin, stretch) = function (sam)
            set_run_value!(sam, :wings, 1, :vel_w, flown_wing .+
                field(spin, stretch, run_value(sam, :wings, 1, :pos_w)))
            set_point_velocities!(sam, (i, pos) -> flown[i] .+ field(spin, stretch, pos))
        end
        @test maximum(abs, damping_rate(sam, rigidly(0.3 .* along, 0.0))) < 1e-9
        @test maximum(abs, damping_rate(sam, rigidly(zeros(3), 0.02))) < 1e-9
        @test maximum(abs, damping_rate(sam, rigidly(0.3 .* [0.0, 0.0, 1.0] × along,
                                                       0.0))) < 1e-9
        node = findfirst(i -> run_value(sam, :points, i, :is_wing_node),
                         eachindex(sam.definition.points))
        flutter(sam) = set_run_value!(sam, :points, node, :vel_w,
                                      flown[node] .+ [0.0, 0.0, 1.0])
        @test maximum(abs, damping_rate(sam, flutter)) > 1.0
    end
    set_data_path(data_path_before)
end
