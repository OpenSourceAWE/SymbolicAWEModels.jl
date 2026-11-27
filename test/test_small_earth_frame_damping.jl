# SPDX-FileCopyrightText: 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_small_earth_frame_damping.jl
#
# Angular and body-frame damping act on motion relative to the kite's small-earth
# frame: z from the tether's base point to the kite, y the kite's span, stretching
# radially with the tether. Steady flight over the sphere, a turn about the tether
# and reel-out are not braked. A tether's inner points damp against the frame that
# carries the tether's straight line from its anchor to its free end. Tested on the
# frame and `body_integration` directly, then on the right-hand side of both
# backends: the damping a state gets is the rate with the damping on minus the rate
# with it off, at the same state.

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

using Test
using LinearAlgebra
using SymbolicAWEModels
using SymbolicAWEModels: KernelBackend, MonolithBackend, sync_params!, sync_initial!,
                         body_integration, small_earth_frame, frame_velocity,
                         damping_base, base_position, tether_frame,
                         body_frame_damp_accel, damped_tether_ends, VortexStepMethod
using KiteUtils

"""Rotation by `angle` [rad] about unit `axis`."""
rotation(axis, angle) = cos(angle) * I + sin(angle) * [0 -axis[3] axis[2];
    axis[3] 0 -axis[1]; -axis[2] axis[1] 0] + (1 - cos(angle)) * axis * axis'

"""
    flying_two_plate(backend, geometry, root, elevation_vel; tether_damping = nothing)

The 2plate kite built from `geometry` on `backend`, initialised while its transform
turns it over the sphere about the ground point at `elevation_vel` [rad/s].
`tether_damping` [1/s] goes in the tether row's `body_frame_damping` column.
"""
function flying_two_plate(backend, geometry, root, elevation_vel; tether_damping = nothing)
    data_path = joinpath(root, string(nameof(typeof(backend))), geometry, "2plate_kite")
    mkpath(dirname(data_path))
    cp(joinpath(dirname(@__DIR__), "data", "2plate_kite"), data_path; force = true)
    isnothing(tether_damping) || add_tether_damping(joinpath(data_path, geometry),
                                                    tether_damping)
    set_data_path(data_path)
    set = Settings("system.yaml")
    vsm_set = VortexStepMethod.VSMSettings(
        joinpath(data_path, "vsm_settings.yaml"); data_prefix = false)
    sys, transforms = load_sys_struct_from_yaml(
        joinpath(data_path, geometry);
        system_name = "small_earth_frame_damping", set, vsm_set)
    transforms[1].elevation_vel = elevation_vel
    sam = SymbolicAWEModel(set, sys; transforms, backend)
    init!(sam; prn = false, remake = true)
    return sam
end

"""Give the tethers of the structure document at `path` a `body_frame_damping` column."""
function add_tether_damping(path, damping)
    lines = readlines(path)
    header = findfirst(contains("headers: [name, start_point"), lines)
    lines[header] = replace(lines[header], "]" => ", body_frame_damping]")
    row = findfirst(contains("- [main_tether"), lines)
    lines[row] = replace(lines[row], r"\]$" => ", $damping]")
    write(path, join(lines, "\n"))
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

"""The small-earth axes `[x y z]` of a kite at `radius` from its base with span `span`."""
function small_earth_axes(radius, span)
    z = normalize(radius)
    x = normalize(span × z)
    return [x z × x z]
end

"""Position [m] of the base point the model's first wing damps against."""
tether_base(sam) = base_position(damping_base(sam.state, sam.state.wings[1].idx),
                                 point -> sam.state.points[point].pos_w)

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

    @testset "a tether frame carries its line from the anchor to the free end" begin
        anchor, anchor_vel = [1.0, -2.0, 0.5], [0.3, 0.0, -0.1]
        free, free_vel = [150.0, 40.0, 90.0], [-4.0, 25.0, 3.0]
        frame = tether_frame(anchor, anchor_vel, free, free_vel)
        line(s) = anchor .+ s .* (free .- anchor)
        line_vel(s) = anchor_vel .+ s .* (free_vel .- anchor_vel)
        damping(s, vel) = body_frame_damp_accel(vel, fill(3.0, 3), I, frame, line(s))
        @test norm(damping(1.0, free_vel)) < 1e-12
        @test maximum(s -> norm(damping(s, line_vel(s))), 0:0.25:1) < 1e-12
        @test damping(0.4, line_vel(0.4) .+ [0.0, 0.0, 1.0]) ≈ [0.0, 0.0, 3.0]
    end

    @testset "a damped tether needs exactly one world-fixed end" begin
        data_path_before = get_data_path()
        set_data_path(joinpath(dirname(@__DIR__), "data", "2plate_kite"))
        set = Settings("system.yaml")
        set_data_path(data_path_before)
        line(first_end) = SystemStructure("damped_line", set;
            points = [SymbolicAWEModels.Point(:a, [0.0, 0.0, 0.0], first_end),
                      SymbolicAWEModels.Point(:m, [0.0, 0.0, 5.0], DYNAMIC;
                                              body_frame_damping = 1.0),
                      SymbolicAWEModels.Point(:b, [0.0, 0.0, 10.0], DYNAMIC)],
            segments = [Segment(:s1, :a, :m, 1e4, 10.0, 0.002),
                        Segment(:s2, :m, :b, 1e4, 10.0, 0.002)],
            tethers = [Tether(:t, [:s1, :s2])])
        @test damped_tether_ends(line(STATIC), 2) == (1, 3)
        @test damped_tether_ends(line(STATIC), 3) === nothing
        @test_throws ErrorException damped_tether_ends(line(DYNAMIC), 2)
    end

    data_path_before = get_data_path()
    root = mktempdir()
    @testset "$(nameof(typeof(backend))) leaves rigid flight on the small-earth frame \
              undamped" for backend in (KernelBackend(), MonolithBackend())
        sam = flying_two_plate(backend, "rigid_structural_geometry.yaml", root, 0.4)
        wing = sam.state.wings[1]
        @test norm(wing.ω_b) > 0.3
        @test maximum(abs, damping_rate(sam)) < 1e-9
        flown = copy(wing.ω_b)
        base = tether_base(sam)
        along = normalize(wing.com_w .- base)
        flown_vel = [copy(point.vel_w) for point in sam.state.points]
        flown_wing_vel = copy(wing.vel_w)
        function turned(state)
            turn = 0.3 .* along
            state.wings[1].ω_b .= flown .+ wing.R_b_to_w' * turn
            state.wings[1].vel_w .= flown_wing_vel .+ turn × (wing.pos_w .- base)
            foreach((point, vel) -> point.vel_w .= vel .+ turn × (point.pos_w .- base),
                    state.points, flown_vel)
        end
        @test maximum(abs, damping_rate(sam, turned)) < 1e-9
        pitched(state) = state.wings[1].ω_b .= flown .+ [0.0, 0.5, 0.0]
        @test maximum(abs, damping_rate(sam, pitched)) > 1.0
    end
    @testset "$(nameof(typeof(backend))) leaves a fitted wing's nodes and its tether \
              undamped in rigid flight on the small-earth frame" for
            backend in (KernelBackend(), MonolithBackend())
        sam = flying_two_plate(backend, "particle_structural_geometry.yaml", root, 0.4;
                               tether_damping = 2.0)
        tether_points = [point for point in sam.state.points
                         if startswith(string(point.name), "main_tether")]
        @test length(tether_points) == 5
        @test all(point.body_frame_damping == fill(2.0, 3) for point in tether_points)
        @test maximum(norm, point.vel_w for point in sam.state.points) > 1.0
        @test maximum(abs, damping_rate(sam)) < 1e-9
        base = tether_base(sam)
        along = normalize(sam.state.wings[1].pos_w .- base)
        flown = [copy(point.vel_w) for point in sam.state.points]
        flown_wing = copy(sam.state.wings[1].vel_w)
        field(spin, stretch, pos) = spin × (pos .- base) .+ stretch .* (pos .- base)
        rigidly(spin, stretch) = function (state)
            wing = state.wings[1]
            wing.vel_w .= flown_wing .+ field(spin, stretch, wing.pos_w)
            foreach((point, vel) -> point.vel_w .= vel .+ field(spin, stretch, point.pos_w),
                    state.points, flown)
        end
        @test maximum(abs, damping_rate(sam, rigidly(0.3 .* along, 0.0))) < 1e-9
        @test maximum(abs, damping_rate(sam, rigidly(zeros(3), 0.02))) < 1e-9
        @test maximum(abs, damping_rate(sam, rigidly(0.3 .* [0.0, 0.0, 1.0] × along,
                                                       0.0))) < 1e-9
        node = findfirst(point -> point.wing_idx > 0, sam.state.points)
        flutter(state) = state.points[node].vel_w .= flown[node] .+ [0.0, 0.0, 1.0]
        @test maximum(abs, damping_rate(sam, flutter)) > 1.0
        wave = tether_points[3].idx
        function waved(state)
            state.wings[1].vel_w .= flown_wing
            foreach((point, vel) -> point.vel_w .= vel, state.points, flown)
            state.points[wave].vel_w .+= [0.0, 0.5, 0.0]
        end
        rate = damping_rate(sam, waved)
        @test count(>(1e-9), abs.(rate)) == 1
        @test maximum(abs, rate) ≈ 2.0 * 0.5
    end
    set_data_path(data_path_before)
end
