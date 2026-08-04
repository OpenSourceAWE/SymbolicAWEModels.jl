# Copyright (c) 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_flap_beam.jl — the two structure-driven aero-coupling features:
#   A. KINEMATIC flap twist_surface: a live deflection δ = the signed angle between
#      two bodies about a hinge axis (fed to the (α, δ) polars by AeroPressure).
#   B. beam-anchored point: a BODY_STATIC point rides a TimoshenkoJoint's deformed
#      corotational-Hermite centerline instead of one rigid body.
#
# Covered here (no pressure-aero fixture needed):
#   1. δ extraction matches the known hinge rotation (Julia ground truth).
#   4. a KINEMATIC flap leaks no ODE state (no twist DOF).
#   5. the beam-anchored point lies on the bent-beam centerline (Hermite),
#      off the straight LE→TE chord.
# The δ→force response of AeroPressure (tests 2–3) needs a δ-swept cp/cf fixture
# and is exercised in the beam-kite integration (BeyondTheSim).

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

@isdefined(test_init!) || include(joinpath(@__DIR__, "util.jl"))

using Test
using SymbolicAWEModels
const S = SymbolicAWEModels
using SymbolicAWEModels: KVec3
using KiteUtils
using LinearAlgebra

quat_y(phi) = [cos(phi / 2), 0.0, sin(phi / 2), 0.0]  # rotation about body y-axis

SETTINGS_YAML = """
system:
    log_file: "data/flap_beam_test"
    g_earth: 0.0
solver:
    solver: "FBDF"
    abs_tol: 1.0e-8
    rel_tol: 1.0e-8
kite:
    model: ""
    foil_file: "ram_air_kite/ram_air_kite_foil.dat"
    physical_model: "flap_beam_test"
    mass: 0.0
tether:
    cd_tether: 0.958
    unit_damping: 0.0
    unit_stiffness: 0.0
    rho_tether: 724.0
    e_tether: 5.5e10
winch:
    winch_model: "TorqueControlledMachine"
    drum_radius: 0.110
    gear_ratio: 1.0
    inertia_total: 0.024
    f_coulomb: 122.0
    c_vf: 30.6
environment:
    rho_0: 1.225
    v_wind: 0.0
    upwind_dir: -90.0
    upwind_elevation: 0.0
    wind_vec: [0.0, 0.0, 0.0]
    profile_law: 0
"""

"""
Corotational-Hermite world position of a beam-anchored point — the Julia mirror
of the symbolic branch in `point_eqs!`, used as the ground truth.
"""
function hermite_point_w(sam, joint, point)
    bodies = sam.sys_struct.bodies
    body_a = bodies[joint.body_a_idx]
    body_b = bodies[joint.body_b_idx]
    R_a = S.quaternion_to_rotation_matrix(body_a.Q_b_to_w)
    R_b = S.quaternion_to_rotation_matrix(body_b.Q_b_to_w)
    x_a = body_a.pos_w .+ R_a * joint.anchor_a_b
    x_b = body_b.pos_w .+ R_b * joint.anchor_b_b
    e1, e2, e3, len = S.timoshenko_element_frame(x_a, x_b, R_a)
    element_frame = [e1 e2 e3]
    Da = (element_frame' * R_a) * joint.R_a_rel0'
    Db = (element_frame' * R_b) * joint.R_b_rel0'
    θ_a = [0.5 * (Da[3, 2] - Da[2, 3]), 0.5 * (Da[1, 3] - Da[3, 1]),
           0.5 * (Da[2, 1] - Da[1, 2])]
    θ_b = [0.5 * (Db[3, 2] - Db[2, 3]), 0.5 * (Db[1, 3] - Db[3, 1]),
           0.5 * (Db[2, 1] - Db[1, 2])]
    s = point.beam_frac
    N2 = len * (s - 2s^2 + s^3)
    N4 = len * (-s^2 + s^3)
    v = N2 * θ_a[3] + N4 * θ_b[3]
    w = -(N2 * θ_a[2] + N4 * θ_b[2])
    x_center = x_a .+ (s * len) .* e1 .+ v .* e2 .+ w .* e3
    return x_center .+ element_frame * point.beam_offset_b
end

@testset "Flap δ + beam-anchored point" begin

    @testset "flap δ extraction (Julia ground truth)" begin
        # Rest capture at φ=0 → rest_delta 0; then δ equals the hinge rotation.
        for phi in (0.0, 0.3, -0.25, 0.5)
            main = Body(:main; mass=1.0, inertia_principal=[1.0, 1.0, 1.0],
                        pos=[0.0, 0, 0])
            flap = Body(:flap; mass=1.0, inertia_principal=[1.0, 1.0, 1.0],
                        pos=[1.0, 0, 0])
            bodies = [main, flap]
            ts = TwistSurface(:flap_ts, Int[], KINEMATIC, 0.0;
                wing=1, flap_bodies=[1, 2], flap_axis=[0.0, 1.0, 0.0])
            ts.flap_body_idxs = [1, 2]
            S.init_twist_surface_flap!(ts, bodies)
            @test ts.flap_rest_delta ≈ 0.0 atol=1e-12
            @test ts.flap_chord_refs == [KVec3(1, 0, 0), KVec3(1, 0, 0)]
            bodies[2].Q_b_to_w .= quat_y(phi)
            R_main = S.quaternion_to_rotation_matrix(bodies[1].Q_b_to_w)
            R_flap = S.quaternion_to_rotation_matrix(bodies[2].Q_b_to_w)
            @test S.flap_delta(ts, R_main, R_flap) ≈ phi atol=1e-10
        end
        # Nonzero rest: capture at φ0, deflect to φ0+Δ → δ = Δ.
        phi0, delta = 0.2, 0.15
        main = Body(:main; mass=1.0, inertia_principal=[1.0, 1.0, 1.0],
                    pos=[0.0, 0, 0], Q_b_to_w=quat_y(0.0))
        flap = Body(:flap; mass=1.0, inertia_principal=[1.0, 1.0, 1.0],
                    pos=[1.0, 0, 0], Q_b_to_w=quat_y(phi0))
        bodies = [main, flap]
        ts = TwistSurface(:flap_ts, Int[], KINEMATIC, 0.0;
            wing=1, flap_bodies=[1, 2], flap_axis=[0.0, 1.0, 0.0])
        ts.flap_body_idxs = [1, 2]
        S.init_twist_surface_flap!(ts, bodies)
        @test ts.flap_rest_delta ≈ phi0 atol=1e-10
        bodies[2].Q_b_to_w .= quat_y(phi0 + delta)
        R_main = S.quaternion_to_rotation_matrix(bodies[1].Q_b_to_w)
        R_flap = S.quaternion_to_rotation_matrix(bodies[2].Q_b_to_w)
        @test S.flap_delta(ts, R_main, R_flap) ≈ delta atol=1e-10
    end

    tmpdir = mktempdir()
    data_path = joinpath(tmpdir, "2plate_kite")
    cp(joinpath(dirname(@__DIR__), "data", "2plate_kite"), data_path; force=true)
    write(joinpath(data_path, "settings.yaml"), SETTINGS_YAML)
    write(joinpath(data_path, "system.yaml"),
        "system:\n  sim_settings: settings.yaml\n")
    set_data_path(data_path)
    set = Settings("system.yaml")

    beam_length = 1.0
    EI = 100.0; GA = 1500.0; EA = 1.0e4; GJ = 50.0; kshear = 5 / 6
    inertia = [0.01, 0.01, 0.01]

    # Fresh components per call: `SystemStructure` resolves indices in place.
    function build_beam_sys(name; with_flap=true)
        nodeA = Body(:nodeA; mass=1.0, inertia_principal=inertia,
                     pos=[0.0, 0.0, 0.0], type=STATIC)
        nodeB = Body(:nodeB; mass=1.0, inertia_principal=inertia,
                     pos=[beam_length, 0.0, 0.0])
        joint = TimoshenkoJoint(:joint, :nodeA, :nodeB;
            EA, GA, GJ, EIy=EI, EIz=EI, shear_coeff=kshear,
            damping_trans=200.0, damping_rot=3.0)
        # KINEMATIC flap twist_surface: δ = the hinge angle nodeA→nodeB.
        flap = TwistSurface(:flap, Int[], KINEMATIC, 0.0;
            wing=:nodeA, flap_bodies=[:nodeA, :nodeB], flap_axis=[0.0, 1.0, 0.0])
        # Beam-anchored point at midspan, offset 0.1 in +z off the centerline.
        bridle = Point(:bridle, [0.5, 0.0, 0.1], BODY_STATIC; joint=:joint)
        return SystemStructure(name, set;
            points=[bridle], bodies=[nodeA, nodeB], timoshenko_joints=[joint],
            twist_surfaces=with_flap ? [flap] : TwistSurface[])
    end

    sys = build_beam_sys("flap_beam_test")

    @testset "structural resolution" begin
        @test sys.twist_surfaces[:flap].flap_body_idxs == [1, 2]
        @test sys.points[:bridle].joint_idx == 1
        @test sys.points[:bridle].beam_frac ≈ 0.5 atol=1e-6
        @test sys.points[:bridle].beam_offset_b[3] ≈ 0.1 atol=1e-6
    end

    sam = SymbolicAWEModel(set, sys)
    jt = sam.sys_struct.timoshenko_joints[:joint]
    rb = sam.sys_struct.bodies[:nodeB]
    br = sam.sys_struct.points[:bridle]

    rb.ext_force_w .= [0.0, 0.0, 5.0]
    test_init!(sam; prn=false)
    for _ in 1:4000
        next_step!(sam; dt=0.01, vsm_interval=0)
        (norm(rb.vel_w) < 1e-8 && norm(rb.ω_b) < 1e-8) && break
    end

    @testset "KINEMATIC flap adds no ODE state" begin
        # Absolute counts shift with MTK's alias elimination, so compare against
        # the same beam without the twist_surface: a DYNAMIC DOF would add 2.
        bare = SymbolicAWEModel(set, build_beam_sys("flap_beam_bare";
                                                    with_flap=false))
        test_init!(bare; prn=false)
        @test length(sam.integrator.u) == length(bare.integrator.u)
        @test all(isfinite, rb.pos_w)
    end

    @testset "beam-anchored point rides the bent centerline" begin
        expected = hermite_point_w(sam, jt, br)
        @test br.pos_w ≈ expected atol=1e-6
        # It bulges off the straight LE→TE chord midpoint when the beam bends.
        ba = sam.sys_struct.bodies[:nodeA]
        straight_mid = 0.5 .* (ba.pos_w .+ rb.pos_w) .+ [0.0, 0.0, 0.1]
        @test norm(br.pos_w - straight_mid) > 1e-4
    end

    rm(tmpdir; recursive=true)
end
nothing
