# SPDX-FileCopyrightText: 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

using Test
using SymbolicAWEModels
using SymbolicAWEModels: KernelBackend, MonolithBackend, RigidLayout, VortexStepMethod
using KiteUtils
using LinearAlgebra

"""
    load_2plate(geometry, system_name; backend=KernelBackend()) -> SymbolicAWEModel

A model of the 2plate kite on `geometry`, loaded from a scratch copy of `data/2plate_kite`.
"""
function load_2plate(geometry, system_name; backend = KernelBackend())
    data_path = joinpath(mktempdir(), "2plate_kite")
    cp(joinpath(dirname(@__DIR__), "data", "2plate_kite"), data_path)
    set_data_path(data_path)
    set = Settings("system.yaml")
    vsm_set = VortexStepMethod.VSMSettings(joinpath(data_path, "vsm_settings.yaml");
                                           data_prefix = false)
    sys, transforms = load_sys_struct_from_yaml(joinpath(data_path, geometry);
                                                system_name, set, vsm_set)
    return SymbolicAWEModel(set, sys; transforms, backend)
end

"""Place `sam` at the given pose [deg] and start its run there on braked winches."""
function start_at_pose!(sam; elevation = 50.0, azimuth = 0.0, heading = 0.0)
    transform = sam.state.transforms[1]
    transform.elevation = deg2rad(elevation)
    transform.azimuth = deg2rad(azimuth)
    transform.heading = deg2rad(heading)
    place!(sam)
    init!(sam; prn = false)
    foreach(winch -> winch.brake = 1.0, sam.state.winches)
    return sam
end

"""Every artificial damping value of `sam`, as `==` compares them."""
damping_values(sam) = (
    [(copy(point.body_frame_damping), copy(point.world_frame_damping))
     for point in sam.state.points],
    [(copy(body.body_frame_damping), copy(body.world_frame_damping),
      copy(body.angular_damping)) for body in sam.state.bodies])

"""Positions of every point of `sam`, one vector per point."""
point_positions(sam) = [Vector(point.pos_w) for point in sam.state.points]

"""Largest distance [m] `reposition!` moves a point of `sam`, zero on its held pose."""
function pose_error(sam)
    placed = point_positions(sam)
    reposition!(sam.state)
    return maximum(norm.(point_positions(sam) .- placed))
end

@testset "quasi-steady finder" begin
    particle = load_2plate("particle_structural_geometry.yaml", "quasi_steady_particle")
    solver = (; dt = 0.01, check_every = 25, prn = false)

    @testset "a perturbed kite converges to the state an unperturbed one finds" begin
        start_at_pose!(particle)
        @test find_quasi_steady!(particle; solver...).converged
        settled = point_positions(particle)

        start_at_pose!(particle)
        layout = RigidLayout(particle)
        particle.integrator.u[layout.point_vel] .+= 2.0
        particle.integrator.u[layout.point_pos[3, :]] .+= 0.05
        before = quasi_steady_residual(particle)
        @test before.force_ratio > 0.1

        result = find_quasi_steady!(particle; solver...)

        @test result.converged
        @test result.status === :converged
        @test max(result.force_ratio, result.moment_ratio) <= 0.01
        @test maximum(norm.(point_positions(particle) .- settled)) < 0.05
        @test pose_error(particle) < 1e-6
    end

    @testset "a parked kite converges barely moving, with its damping back" begin
        start_at_pose!(particle; elevation = 55.0)
        particle.state.points[1].body_frame_damping .= 3.0
        damping_before = damping_values(particle)

        result = find_quasi_steady!(particle; solver...)

        @test result.converged
        @test result.speed < 0.5
        @test abs(result.turn_rate) < 0.02
        @test damping_values(particle) == damping_before
        @test pose_error(particle) < 1e-6
        particle.state.points[1].body_frame_damping .= 0.0
    end

    @testset "the residual judges the state the finder leaves" begin
        start_at_pose!(particle)
        result = find_quasi_steady!(particle; solver...)
        residual = quasi_steady_residual(particle)
        @test residual.spin ≈ result.spin atol = 1e-6
        @test max(residual.force_ratio, residual.moment_ratio) <= 0.02
    end

    @testset "a pose with no equilibrium is reported with its residual" begin
        for azimuth in (60.0, 90.0)
            start_at_pose!(particle; azimuth)
            result = find_quasi_steady!(particle; solver...)
            @test !result.converged
            @test result.status in (:aero_failed, :max_steps)
            @test isfinite(result.force_ratio) && isfinite(result.moment_ratio)
            @test max(result.force_ratio, result.moment_ratio) > 0.01
            @test norm(result.force) > 0.0
            @test norm(result.moment) > 0.0
        end
    end

    @testset "the log holds one frame of the state the finder leaves" begin
        start_at_pose!(particle)
        find_quasi_steady!(particle; solver...)
        log = quasi_steady_log(particle)
        @test length(log.syslog) == 1
        frame = log.syslog[1]
        for (idx, point) in enumerate(particle.state.points)
            @test [frame.X[idx], frame.Y[idx], frame.Z[idx]] ≈ point.pos_w rtol = 1e-5
        end
    end

    @testset "a rigid wing's body turns with the pose and its torque is judged" begin
        rigid = load_2plate("rigid_structural_geometry.yaml", "quasi_steady_rigid")
        start_at_pose!(rigid)
        result = find_quasi_steady!(rigid; solver..., damping = 5.0)
        @test result.converged
        @test pose_error(rigid) < 1e-6
    end

    @testset "a tube's damping is off inside the judgement and back after it" begin
        set = Settings("system.yaml")
        bodies = [Body(:clamped; extra_mass = 1.0, inertia_principal = [0.01, 0.01, 0.01],
                       pos = [0.0, 0.0, 0.0], type = STATIC),
                  Body(:free; extra_mass = 1.0, inertia_principal = [0.01, 0.01, 0.01],
                       pos = [1.0, 0.0, 0.0])]
        tube = Tube(:tube, :clamped, :free; diameter = 0.1, pressure = 3e4,
                    model = TimoshenkoTube(; EA = 1e4, GA = 1500.0, GJ = 50.0, EIy = 100.0,
                                           EIz = 100.0, damping = 0.05))
        sys = SystemStructure("quasi_steady_tube", set; bodies, tubes = [tube])
        sam = SymbolicAWEModel(set, sys; backend = KernelBackend())
        init!(sam; prn = false)
        layout = RigidLayout(sam)
        u = copy(sam.integrator.u)
        u[layout.body_vel] .= 1.0
        acceleration() = begin
            du = similar(u)
            sam.integrator.f(du, u, sam.integrator.p, 0.0)
            du[layout.body_vel]
        end
        damped = acceleration()
        @test SymbolicAWEModels.without_damping(acceleration, sam) != damped
        @test sam.state.tubes[:tube].model.damping == 0.05
        @test acceleration() == damped
    end

    @testset "a monolith model is refused by name" begin
        monolith = load_2plate("particle_structural_geometry.yaml", "quasi_steady_monolith";
                               backend = MonolithBackend())
        @test_throws "KernelBackend" find_quasi_steady!(monolith)
    end
end
