# SPDX-FileCopyrightText: 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_principal_frame_invariance.jl
#
# The principal frame is internal ODE state: any orthonormal frame that
# diagonalizes the wing's inertia tensor must produce the exact same
# simulation, because everything physical goes through
# R_b_to_w = R_p_to_w * R_b_to_p and R_b_to_p absorbs the frame choice.
# Regression guard for the A1-15 failure, where a ~90°-flipped frame
# leaked into the body frame and diverged: fly the wing with the frame each
# principal-frame method picks, ~90° apart about y, and assert world-frame
# forces and a short dynamic run are unchanged.

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

using Test
using SymbolicAWEModels
using SymbolicAWEModels: VortexStepMethod, RIGID_DYNAMICS,
    quaternion_to_rotation_matrix, wing_body
using KiteUtils
using LinearAlgebra

"""
    frame_snapshot(sam)

World-frame observables: wing attitude, aero force/moment, point
positions and forces (tether/bridle loads), segment tensions, twists.
"""
function frame_snapshot(sam)
    sys = sam.definition
    wing(quantity) = run_value(sam, :wings, 1, quantity)
    items(collection) = eachindex(getproperty(sys, collection))
    value(collection, quantity) =
        [run_value(sam, collection, i, quantity) for i in items(collection)]
    R_b_to_w = quaternion_to_rotation_matrix(Vector(wing(:Q_b_to_w)))
    return (;
        R_b_to_w,
        aero_wrench_w = [R_b_to_w * Vector(wing(:aero_force_b));
                         R_b_to_w * Vector(wing(:aero_moment_b))],
        point_pos_w = value(:points, :pos_w),
        point_force_w = value(:points, :force),
        scalars = [value(:segments, :force); value(:spanwise_stations, :twist);
                   norm(run_value(sam, :winches, 1, :force))])
end

"""
    compare_snapshots(snap, ref; rtol, atol_pos)

Assert a variant snapshot matches the baseline: forces relative to the
largest baseline force, positions with `atol_pos`.
"""
function compare_snapshots(snap, ref; rtol, atol_pos)
    force_scale = max(norm(ref.aero_wrench_w),
        maximum(norm, ref.point_force_w), maximum(abs, ref.scalars))
    @test norm(snap.R_b_to_w - ref.R_b_to_w) < rtol
    @test norm(snap.aero_wrench_w - ref.aero_wrench_w) < rtol * force_scale
    @test maximum(norm.(snap.point_pos_w .- ref.point_pos_w)) < atol_pos
    @test maximum(norm.(snap.point_force_w .- ref.point_force_w)) <
        rtol * force_scale
    @test maximum(abs.(snap.scalars .- ref.scalars)) < rtol * force_scale
    return nothing
end

"""
    run_case(sam, method; steps=10, dt=0.05)

Give the wing the principal frame `method` finds, place and init the structure, take one
tiny step so all force outputs are realised, then step the dynamics; returns the snapshots
(the reported failure grew over time, so the run must be long enough to expose a
diverging variant).
"""
function run_case(sam, method; steps=10, dt=0.05)
    sam.definition.bodies[wing_body(sam, 1)].principal_frame_method = method
    place!(sam)
    init!(sam; prn=false)
    next_step!(sam; dt=1e-5, vsm_interval=1)
    snaps = [frame_snapshot(sam)]
    for _ in 1:steps
        next_step!(sam; dt, vsm_interval=1)
        push!(snaps, frame_snapshot(sam))
    end
    return snaps
end

@testset "Principal frame invariance" begin
    pkg_root = dirname(@__DIR__)
    tmpdir = mktempdir()
    data_path = joinpath(tmpdir, "2plate_kite")
    cp(joinpath(pkg_root, "data", "2plate_kite"), data_path; force=true)
    set_data_path(data_path)
    set = Settings("system.yaml")
    # Tight tolerances so frame-choice leaks can't hide in integrator noise.
    set.abs_tol = 1e-7
    set.rel_tol = 1e-7
    set.v_wind = 15.0

    vsm_set = VortexStepMethod.VSMSettings(
        joinpath(data_path, "vsm_settings.yaml"); data_prefix=false)
    sys, transforms = load_authoring(
        joinpath(data_path, "rigid_structural_geometry.yaml");
        name="frame_invariance", set)
    # Yaw inertia on the wing makes its z axis the largest principal one, which
    # EIGEN_DECOMP keeps near the body z axis and Y_ROTATION turns ~90° onto x.
    wing_row = sys.bodies[sys.wings[1].body]
    wing_row.extra_inertia_KA =
        something(wing_row.extra_inertia_KA, zeros(3, 3)) + Diagonal([0.0, 0.0, 3.0])
    aero = Dict(wing.name => AeroLinearized() for wing in sys.wings)
    sam = SymbolicAWEModel(set, sys; transforms, vsm_set, aero)
    @test body_dynamics_type(sam, wing_body(sam, 1)) == RIGID_DYNAMICS
    ref = run_case(sam, "EIGEN_DECOMP")
    @test norm(ref[1].aero_wrench_w) > 1.0
    @test maximum(norm, ref[1].point_force_w) > 1.0
    R_p_to_w_ref = quaternion_to_rotation_matrix(Vector(run_value(sam, :wings, 1,
                                                                  :Q_p_to_w)))

    # The 90°-about-y flip from the A1-15 failure.
    S = [0.0 0 1; 0 1 0; -1 0 0]
    snaps = run_case(sam, "Y_ROTATION")
    # Prove the variant took effect: principal attitude rotated by S.
    R_p_to_w = quaternion_to_rotation_matrix(Vector(run_value(sam, :wings, 1, :Q_p_to_w)))
    @test norm(R_p_to_w - R_p_to_w_ref) > 0.5
    @test norm(R_p_to_w * S' - R_p_to_w_ref) < 0.05

    compare_snapshots(snaps[1], ref[1]; rtol=1e-4, atol_pos=1e-6)
    for idx in 2:lastindex(ref)
        compare_snapshots(snaps[idx], ref[idx]; rtol=1e-2, atol_pos=1e-3)
    end

    rm(tmpdir; recursive=true)
end
nothing
