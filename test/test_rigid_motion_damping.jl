# SPDX-FileCopyrightText: 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_rigid_motion_damping.jl
#
# Body-frame damping measured against the transform's rigid motion
# (`body_damping_reference = :rigid_motion`), tested on the right-hand side itself:
# the damping a state gets is the rate with the damping on minus the rate with it
# off, at the same state.

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

@isdefined(test_init!) || include(joinpath(@__DIR__, "util.jl"))

using Test
using LinearAlgebra
using SymbolicAWEModels
using SymbolicAWEModels: KernelBackend, MonolithBackend, buffer_slots, sync_params!,
                         sync_initial!, update_sys_struct!, rigid_motion_of_wing,
                         get_model_name, reposition!

"""Velocity slots of every free point in the kernel state, with the point's index."""
function free_point_slots(sam)
    model = sam.prob.initial_sync.model
    slots(instance, name) = collect(buffer_slots(model.system, instance, :states, name))
    return [(slots(model.point_instances[idx], :pos), slots(model.point_instances[idx],
             :vel), idx)
            for (idx, role) in enumerate(model.point_roles)
            if role.kind in (:particle, :wing_node)]
end

"""Rate of state `u` under the model's right-hand side."""
function state_rate(sam, u)
    du = similar(u)
    sam.prob.prob.f(du, u, sam.integrator.p, sam.integrator.t)
    return du
end

"""
    damping_accelerations(sam, u, points) -> Vector of 3-vectors

The acceleration the body-frame damping adds to each free point at state `u`: the
rate with the damping on minus the rate with every coefficient zeroed.
"""
function damping_accelerations(sam, u, points)
    sys = sam.sys_struct
    damped = state_rate(sam, u)
    coefficients = [copy(point.body_frame_damping) for point in sys.points]
    foreach(point -> point.body_frame_damping .= 0.0, sys.points)
    sync_params!(sam.prob.param_sync, sam.integrator, sys)
    undamped = state_rate(sam, u)
    foreach(((point, c),) -> point.body_frame_damping .= c, zip(sys.points, coefficients))
    sync_params!(sam.prob.param_sync, sam.integrator, sys)
    return [damped[vs] .- undamped[vs] for (_, vs, _) in points]
end

"""`u` with every free point's velocity replaced by `velocity_of(pos, vel)`."""
function with_point_velocities(velocity_of, u, points)
    moved = copy(u)
    for (ps, vs, _) in points
        moved[vs] .= velocity_of(u[ps], u[vs])
    end
    return moved
end

"""Velocity at `pos` of the rigid motion `spin`, `stretch` about `base`."""
rigid_velocity(spin, stretch, base, pos) =
    cross(spin, pos .- base) .+ stretch .* (pos .- base)

"""Moment about `origin` of the damping forces `accelerations` puts on the points."""
damping_moment(sam, u, points, accelerations, origin) =
    sum(cross(u[ps] .- origin, sam.sys_struct.points[idx].total_mass .* acceleration)
        for ((ps, _, idx), acceleration) in zip(points, accelerations))

"""
    added_spin_moment(sam, u, points, origin, spin)

How much the damping moment about `origin` changes when every point of state `u`
also spins at `spin` about `origin`.
"""
function added_spin_moment(sam, u, points, origin, spin)
    spun = with_point_velocities((pos, vel) -> vel .+ cross(spin, pos .- origin), u,
                                 points)
    return damping_moment(sam, spun, points, damping_accelerations(sam, spun, points),
                          origin) .-
        damping_moment(sam, u, points, damping_accelerations(sam, u, points), origin)
end

@testset "Body-frame damping against the rigid motion" begin
    data_path_before = get_data_path()
    root = mktempdir()
    sam = two_plate_model(KernelBackend(), RIGID_MOTION_GEOMETRY, root;
                          prepare = write_rigid_motion_geometry,
                          system_name = "rigid_motion_damping")
    sys = sam.sys_struct
    wing = sys.wings[1]
    transform = sys.transforms[1]
    points = free_point_slots(sam)
    u0 = copy(sam.integrator.u)
    base = collect(transform.base_w)
    origin = collect(wing.pos_w)
    R = collect(wing.R_b_to_w)
    radial = normalize(origin .- base)
    span = R[:, 2]
    # Spun about the span through the wing origin: a motion the damping must see.
    pitched = with_point_velocities((pos, vel) -> vel .+ cross(0.5 .* span,
                                    pos .- origin), u0, points)
    scale = maximum(norm, damping_accelerations(sam, pitched, points))

    @testset "only the transform names it in the model" begin
        @test transform.body_damping_reference === :rigid_motion
        @test occursin("_rigiddamp", get_model_name(sam.set, sys))
        @test_throws ErrorException Transform(:t, 0.0, 0.0, 0.0; wing = 1,
            base_pos = zeros(3), base_point = 1, body_damping_reference = :fit)
    end

    @testset "leaves every rigid motion about the base undamped" begin
        @test scale > 1.0
        spin = [0.1, -0.2, 0.3]
        for (label, motion_spin, stretch) in (("turn", spin, 0.0),
                ("turn plus turn rate about the radial", spin .+ 0.5 .* radial, 0.0),
                ("turn plus reel-out", spin, 0.02))
            u = with_point_velocities((pos, vel) ->
                rigid_velocity(motion_spin, stretch, base, pos), u0, points)
            left = maximum(norm, damping_accelerations(sam, u, points))
            @test left < 1e-9 * scale
        end
    end

    @testset "is -R(c .* R'(v - v_ref)) against the wing's own rigid motion" begin
        sam.integrator.u .= pitched
        update_sys_struct!(sam.prob, sam.integrator, sys)
        motion = rigid_motion_of_wing(wing.pos_w, wing.vel_w,
                                      collect(wing.R_b_to_w) * collect(wing.ω_b), base)
        worst = 0.0
        power = 0.0
        for ((ps, vs, idx), acceleration) in zip(points,
                                                damping_accelerations(sam, pitched, points))
            point = sys.points[idx]
            relative = pitched[vs] .- rigid_velocity(motion.spin, motion.stretch, base,
                                                     pitched[ps])
            expected = -(R * (collect(point.body_frame_damping) .* (R' * relative)))
            worst = max(worst, norm(acceleration .- expected))
            power += point.total_mass * dot(acceleration, relative)
        end
        @test worst < 1e-9 * scale
        @test power < 0
        sam.integrator.u .= u0
        update_sys_struct!(sam.prob, sam.integrator, sys)
    end

    @testset "damps a spin about the span, not the turn rate about the radial" begin
        pitch = added_spin_moment(sam, u0, points, origin, 0.5 .* span)
        @test dot(pitch, span) < 0
        @test norm(added_spin_moment(sam, u0, points, origin, 0.5 .* radial)) <
            1e-9 * norm(pitch)
    end

    @testset "reposition! hands back a state it damps the same on the wing axes" begin
        before = [R' * a for a in damping_accelerations(sam, pitched, points)]
        sam.integrator.u .= pitched
        update_sys_struct!(sam.prob, sam.integrator, sys)
        transform.heading += deg2rad(3.0)
        transform.elevation += deg2rad(1.0)
        reposition!(sys.transforms, sys)
        moved = copy(pitched)
        sync_initial!(sam.prob.initial_sync, (; u0 = moved), sys)
        sam.integrator.u .= moved
        update_sys_struct!(sam.prob, sam.integrator, sys)
        R_moved = collect(wing.R_b_to_w)
        turn = R_moved * R'
        after = [R_moved' * a for a in damping_accelerations(sam, moved, points)]
        wrong_way = with_point_velocities((pos, vel) -> turn' * (turn' * vel), moved,
                                          points)
        wrong = [R_moved' * a for a in damping_accelerations(sam, wrong_way, points)]
        error_after = maximum(norm.(after .- before))
        @test error_after < 1e-9 * scale
        @test maximum(norm.(wrong .- before)) > 1e-3 * scale
    end

    set_data_path(data_path_before)
end
