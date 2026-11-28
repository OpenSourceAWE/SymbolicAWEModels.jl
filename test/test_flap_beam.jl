# Copyright (c) 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_flap_beam.jl — the two structure-driven aero-coupling features:
#   A. KINEMATIC flap station: a live deflection δ = the signed angle between
#      two bodies about a hinge axis (fed to the (α, δ) polars by AeroPressure).
#   B. beam-anchored point: a BODY_STATIC point rides a Timoshenko tube's deformed
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
Corotational-Hermite world position of point `point` of `sam`, anchored on tube `tube` —
the Julia mirror of the symbolic branch in `point_eqs!`, used as the ground truth.
"""
function hermite_point_w(sam, tube, point)
    body_a, body_b = sam.definition.tubes[tube].bodies
    body(name, quantity) = run_value(sam, :bodies, name, quantity)
    element(quantity) = run_value(sam, :tubes, tube, quantity)
    R_a = S.quaternion_to_rotation_matrix(body(body_a, :Q_b_to_w))
    R_b = S.quaternion_to_rotation_matrix(body(body_b, :Q_b_to_w))
    x_a = body(body_a, :pos_w) .+ R_a * element(:anchor_a_b)
    x_b = body(body_b, :pos_w) .+ R_b * element(:anchor_b_b)
    e1, e2, e3, len = S.timoshenko_element_frame(x_a, x_b, R_a)
    element_frame = [e1 e2 e3]
    Da = (element_frame' * R_a) * element(:R_a_rel0)'
    Db = (element_frame' * R_b) * element(:R_b_rel0)'
    θ_a = [0.5 * (Da[3, 2] - Da[2, 3]), 0.5 * (Da[1, 3] - Da[3, 1]),
           0.5 * (Da[2, 1] - Da[1, 2])]
    θ_b = [0.5 * (Db[3, 2] - Db[2, 3]), 0.5 * (Db[1, 3] - Db[3, 1]),
           0.5 * (Db[2, 1] - Db[1, 2])]
    s = run_value(sam, :points, point, :beam_frac)
    N2 = len * (s - 2s^2 + s^3)
    N4 = len * (-s^2 + s^3)
    v = N2 * θ_a[3] + N4 * θ_b[3]
    w = -(N2 * θ_a[2] + N4 * θ_b[2])
    x_center = x_a .+ (s * len) .* e1 .+ v .* e2 .+ w .* e3
    return x_center .+ element_frame * run_value(sam, :points, point, :beam_offset_b)
end

@testset "Flap δ + beam-anchored point" begin
    tmpdir = mktempdir()
    data_path = joinpath(tmpdir, "2plate_kite")
    cp(joinpath(dirname(@__DIR__), "data", "2plate_kite"), data_path; force=true)
    write(joinpath(data_path, "settings.yaml"), SETTINGS_YAML)
    write(joinpath(data_path, "system.yaml"),
        "system:\n  sim_settings: settings.yaml\n")
    set_data_path(data_path)
    set = Settings("system.yaml")

    # A run of a wing on body `main` whose KINEMATIC station `flap` hinges from `main`
    # to body `flap`, the flap at rest turned by `flap_Q`.
    function flap_run(flap_Q)
        unit_inertia = Matrix(Diagonal([1.0, 1.0, 1.0]))
        main = Body(; name=:main, pos_ENU=[0.0, 0, 0], extra_mass=1.0,
                    extra_inertia_KA=unit_inertia)
        flap = Body(; name=:flap, pos_ENU=[1.0, 0, 0], Q_KA2ENU=flap_Q, extra_mass=1.0,
                    extra_inertia_KA=unit_inertia)
        wing = KiteGeometry.Wing(; name=:main, pos_ENU=zeros(3), body=:main,
            model=WingModel(; aero_mode="AeroNone", stations=[:flap]))
        station = SpanwiseStation(; name=:flap, wing=:main, points=Int[],
            model=StationModel(; type="KINEMATIC", moment_frac=0.0,
                               flap_bodies=[:main, :flap], flap_axis=[0.0, 1.0, 0.0]))
        definition = system_definition("flap_delta_test", set; wings=[wing],
            bodies=[main, flap], spanwise_stations=[station])
        return SymbolicAWEModel(nothing, definition)
    end

    @testset "flap δ extraction (Julia ground truth)" begin
        # Rest capture at φ=0 → rest_delta 0; then δ equals the hinge rotation.
        for phi in (0.0, 0.3, -0.25, 0.5)
            sam = flap_run(quat_y(0.0))
            @test S.flap_rest_delta(sam, 1) ≈ 0.0 atol=1e-12
            @test S.flap_chord_refs(sam, 1) == [KVec3(1, 0, 0), KVec3(1, 0, 0)]
            R_main = S.initial_rotation(sam.definition.bodies[:main])
            R_flap = S.quaternion_to_rotation_matrix(quat_y(phi))
            @test S.flap_delta_expression(sam, 1, R_main, R_flap) ≈ phi atol=1e-10
        end
        # Nonzero rest: capture at φ0, deflect to φ0+Δ → δ = Δ.
        phi0, delta = 0.2, 0.15
        sam = flap_run(quat_y(phi0))
        @test S.flap_rest_delta(sam, 1) ≈ phi0 atol=1e-10
        R_main = S.initial_rotation(sam.definition.bodies[:main])
        R_flap = S.quaternion_to_rotation_matrix(quat_y(phi0 + delta))
        @test S.flap_delta_expression(sam, 1, R_main, R_flap) ≈ delta atol=1e-10
    end

    beam_length = 1.0
    EI = 100.0; GA = 1500.0; EA = 1.0e4; GJ = 50.0; kshear = 5 / 6
    inertia = [0.01, 0.01, 0.01]

    # Fresh components per call: `SystemDefinition` resolves indices in place.
    function build_beam_sys(name; with_flap=true, wing_body_last=false)
        # The wing owning the flap: a rigid wing with no aero, held where it is.
        wing = KiteGeometry.Wing(; name=:nodeA, pos_ENU=zeros(3), body=:nodeA,
            model=WingModel(; aero_mode="AeroNone",
                            stations=with_flap ? [:flap] : nothing))
        nodeA = Body(; name=:nodeA, pos_ENU=zeros(3), extra_mass=1.0,
                     extra_inertia_KA=Matrix(Diagonal(inertia)), fixed=true,
                     model=BodyModel(; angular_damping=0.0))
        nodeB = Body(; name=:nodeB, pos_ENU=[beam_length, 0.0, 0.0], extra_mass=1.0,
                     extra_inertia_KA=Matrix(Diagonal(inertia)))
        materials = [Material(; name=:breukels2011, law="breukels2011",
            parameters=KiteGeometry.tube_law_parameters(:breukels2011))]
        tube = Tube(; name=:tube, bodies=(:nodeA, :nodeB), diameter=0.1, pressure=3e4,
            material=:breukels2011,
            model=TimoshenkoTubeModel(; EA, GA, GJ, EIy=EI, EIz=EI, shear_coeff=kshear,
                                      damping=0.05))
        # KINEMATIC flap station: δ = the hinge angle nodeA→nodeB.
        flap = SpanwiseStation(; name=:flap, wing=:nodeA, points=Int[],
            model=StationModel(; type="KINEMATIC", moment_frac=0.0,
                               flap_bodies=[:nodeA, :nodeB], flap_axis=[0.0, 1.0, 0.0]))
        # Beam-anchored point at midspan, offset 0.1 in +z off the centerline.
        bridle = Point(; name=:bridle, pos_ENU=[0.5, 0.0, 0.1], tube=:tube)
        return system_definition(name, set;
            points=[bridle], wings=[wing],
            bodies=wing_body_last ? [nodeB, nodeA] : [nodeA, nodeB], tubes=[tube],
            materials, spanwise_stations=with_flap ? [flap] : SpanwiseStation[])
    end

    sys = build_beam_sys("flap_beam_test")
    sam = SymbolicAWEModel(set, sys)
    node_b(quantity) = run_value(sam, :bodies, :nodeB, quantity)
    bridle_pos() = run_value(sam, :points, :bridle, :pos_w)

    @testset "structural resolution" begin
        @test sys.spanwise_stations[:flap].flap_bodies == [1, 2]
        @test sys.points[:bridle].tube == 1
        @test run_value(sam, :points, :bridle, :beam_frac) ≈ 0.5 atol=1e-6
        @test run_value(sam, :points, :bridle, :beam_offset_b)[3] ≈ 0.1 atol=1e-6
    end

    @testset "a wing's body need not lead the bodies" begin
        last_sam = SymbolicAWEModel(set, build_beam_sys("flap_beam_test";
                                                        wing_body_last=true))
        @test SymbolicAWEModels.wing_body(last_sam, 1) == 2
        for (collection, item, quantity) in ((:bodies, :nodeA, :total_mass),
                (:bodies, :nodeA, :inertia_principal), (:bodies, :nodeB, :total_mass),
                (:points, :bridle, :beam_frac), (:points, :bridle, :beam_offset_b),
                (:spanwise_stations, :flap, :flap_rest_delta))
            @test run_value(last_sam, collection, item, quantity) ≈
                  run_value(sam, collection, item, quantity)
        end
    end

    test_init!(sam; prn=false)
    set_run_value!(sam, :bodies, :nodeB, :ext_force_w, [0.0, 0.0, 5.0])
    for _ in 1:4000
        next_step!(sam; dt=0.01, vsm_interval=0)
        (norm(node_b(:vel_w)) < 1e-8 && norm(node_b(:ω_b)) < 1e-8) && break
    end

    @testset "KINEMATIC flap adds no ODE state" begin
        # Absolute counts shift with MTK's alias elimination, so compare against
        # the same beam without the station: a DYNAMIC DOF would add 2.
        bare = SymbolicAWEModel(set, build_beam_sys("flap_beam_bare";
                                                    with_flap=false))
        test_init!(bare; prn=false)
        @test length(sam.integrator.u) == length(bare.integrator.u)
        @test all(isfinite, node_b(:pos_w))
    end

    @testset "beam-anchored point rides the bent centerline" begin
        expected = hermite_point_w(sam, :tube, :bridle)
        @test bridle_pos() ≈ expected atol=1e-6
        # It bulges off the straight LE→TE chord midpoint when the beam bends.
        node_a_pos = run_value(sam, :bodies, :nodeA, :pos_w)
        straight_mid = 0.5 .* (node_a_pos .+ node_b(:pos_w)) .+ [0.0, 0.0, 0.1]
        @test norm(bridle_pos() - straight_mid) > 1e-4
    end

    rm(tmpdir; recursive=true)
end
nothing
