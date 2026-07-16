# SPDX-FileCopyrightText: 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_principal_frame_invariance.jl
#
# The principal frame is internal ODE state: any orthonormal frame that
# diagonalizes the wing's inertia tensor must produce the exact same
# simulation, because everything physical goes through
# R_b_to_w = R_p_to_w * R_b_to_p and R_b_to_p absorbs the frame choice.
# Regression guard for the A1-15 failure, where a ~90°-flipped frame
# leaked into the body frame and diverged: swap in the flipped frame and
# assert world-frame forces and a short dynamic run are unchanged.

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

using Test
using SymbolicAWEModels
using SymbolicAWEModels: VortexStepMethod, RIGID_DYNAMICS,
    quaternion_to_rotation_matrix
using KiteUtils
using LinearAlgebra

"""
    apply_principal_frame!(wing, S)

Re-express the wing's principal frame through the axis change `S`
(new-principal → old-principal, a proper signed permutation so the
inertia stays diagonal).
"""
function apply_principal_frame!(wing, S)
    @assert S * S' ≈ I(3) && det(S) ≈ 1.0
    inertia_new = S' * Diagonal(collect(wing.inertia_principal)) * S
    @assert norm(inertia_new - Diagonal(diag(inertia_new))) <
        1e-9 * norm(diag(inertia_new))
    wing.inertia_principal .= diag(inertia_new)
    wing.damping .= diag(S' * Diagonal(collect(wing.damping)) * S)
    wing.R_p_to_c .= wing.R_p_to_c * S
    wing.R_b_to_p .= S' * wing.R_b_to_p
    return nothing
end

"""
    frame_snapshot(sys_struct)

World-frame observables: wing attitude, aero force/moment, point
positions and forces (tether/bridle loads), segment tensions, twists.
"""
function frame_snapshot(sys_struct)
    wing = sys_struct.wings[1]
    R_b_to_w = quaternion_to_rotation_matrix(Vector(wing.Q_b_to_w))
    return (;
        R_b_to_w,
        aero_wrench_w = [R_b_to_w * Vector(wing.aero_force_b);
                         R_b_to_w * Vector(wing.aero_moment_b)],
        point_pos_w = [copy(point.pos_w) for point in sys_struct.points],
        point_force_w = [copy(point.force) for point in sys_struct.points],
        scalars = [[segment.force for segment in sys_struct.segments];
                   [surface.twist for surface in sys_struct.twist_surfaces];
                   norm(sys_struct.winches[1].force)])
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
    run_case(sam; steps=10, dt=0.05)

Init, take one tiny step so all force outputs are realised, then step
the dynamics; returns the snapshots (the reported failure grew over
time, so the run must be long enough to expose a diverging variant).
"""
function run_case(sam; steps=10, dt=0.05)
    init!(sam; prn=false)
    next_step!(sam; dt=1e-5, vsm_interval=1)
    snaps = [frame_snapshot(sam.sys_struct)]
    for _ in 1:steps
        next_step!(sam; dt, vsm_interval=1)
        push!(snaps, frame_snapshot(sam.sys_struct))
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
    sys = load_sys_struct_from_yaml(
        joinpath(data_path, "rigid_structural_geometry.yaml");
        system_name="frame_invariance", set, vsm_set,
        aero_mode=AeroLinearized())
    wing = sys.wings[1]
    @test wing.dynamics_type == RIGID_DYNAMICS

    sam = SymbolicAWEModel(set, sys)
    ref = run_case(sam)
    @test norm(ref[1].aero_wrench_w) > 1.0
    @test maximum(norm, ref[1].point_force_w) > 1.0
    R_p_to_w_ref = quaternion_to_rotation_matrix(Vector(wing.Q_p_to_w))

    # The 90°-about-y flip from the A1-15 failure.
    S = [0.0 0 1; 0 1 0; -1 0 0]
    apply_principal_frame!(wing, S)
    @test wing.R_b_to_p ≈ wing.R_p_to_c' * wing.R_b_to_c atol=1e-12

    snaps = run_case(sam)
    # Prove the variant took effect: principal attitude rotated by S.
    R_p_to_w = quaternion_to_rotation_matrix(Vector(wing.Q_p_to_w))
    @test norm(R_p_to_w - R_p_to_w_ref) > 0.5
    @test norm(R_p_to_w * S' - R_p_to_w_ref) < 0.05

    compare_snapshots(snaps[1], ref[1]; rtol=1e-4, atol_pos=1e-6)
    for idx in 2:lastindex(ref)
        compare_snapshots(snaps[idx], ref[idx]; rtol=1e-2, atol_pos=1e-3)
    end

    rm(tmpdir; recursive=true)
end
nothing
