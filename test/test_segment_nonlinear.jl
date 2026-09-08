# SPDX-FileCopyrightText: 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_segment_nonlinear.jl - Validate a nonlinear segment spring: `unit_stiffness`
# as a callable force law F(ε) of the axial strain ε = (len − l0)/l0 (returning
# force [N]), as opposed to a constant per-unit stiffness. A mass hangs on one
# vertical segment under gravity; at equilibrium the spring tension balances mg,
# which fixes ε via the supplied law, so the settled length has a closed form.

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

@isdefined(test_init!) || include(joinpath(@__DIR__, "util.jl"))

using Test
using SymbolicAWEModels
using SymbolicAWEModels: KVec3
using KiteUtils
using LinearAlgebra

SETTINGS_YAML = """
system:
    log_file: "data/segment_nl_test"
    g_earth: 9.81
solver:
    solver: "FBDF"
    abs_tol: 1.0e-8
    rel_tol: 1.0e-8
kite:
    model: ""
    foil_file: "ram_air_kite/ram_air_kite_foil.dat"
    physical_model: "segment_nl_test"
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

# Step until the hanging mass is at rest, so the comparison is limited by the
# spring law, not by an unconverged transient.
function settle!(sam, point; dt=0.01, max_steps=8000, vtol=1e-8)
    for _ in 1:max_steps
        next_step!(sam; dt, vsm_interval=0)
        norm(point.vel_w) < vtol && break
    end
end

@testset "Nonlinear segment spring" begin
    pkg_root = dirname(@__DIR__)
    tmpdir = mktempdir()
    data_path = joinpath(tmpdir, "2plate_kite")
    cp(joinpath(pkg_root, "data", "2plate_kite"), data_path; force=true)
    write(joinpath(data_path, "settings.yaml"), SETTINGS_YAML)
    write(joinpath(data_path, "system.yaml"),
        "system:\n  sim_settings: settings.yaml\n")
    set_data_path(data_path)
    set = Settings("system.yaml")

    rest_length = 2.0
    mass = 3.0
    damping = 500.0
    diameter = 0.001
    weight = mass * set.g_earth

    # Force law F(ε) = c1·ε [N]; piecewise-linear in ε so the LinearInterpolation
    # is exact, and the strain argument (not a constant) is what's exercised.
    c1 = 3000.0
    knots = collect(-0.05:0.005:0.05)
    force_law = SymbolicAWEModels.LinearInterpolation(c1 .* knots, knots)

    ground = Point(:ground, KVec3(0.0, 0.0, 0.0), STATIC)
    hang = Point(:hang, KVec3(0.0, 0.0, -rest_length), DYNAMIC; extra_mass=mass)
    spring = Segment(:spring, :ground, :hang, force_law, damping, diameter;
                     l0=rest_length, density=0.0)
    sys = SystemStructure("segment_nl_test", set;
        points=[ground, hang], segments=[spring])

    @testset "Model setup" begin
        @info "Wiring: one segment with a callable (nonlinear) unit_stiffness."
        @test length(sys.segments) == 1
        @test !(sys.segments[:spring].unit_stiffness isa Real)
        @test sys.points[:hang].extra_mass == mass
    end

    sam = SymbolicAWEModel(set, sys)
    node = sam.sys_struct.points[:hang]

    @testset "Hanging-mass equilibrium" begin
        test_init!(sam; prn=false)
        settle!(sam, node)
        # Tension balances gravity: c1·ε = m·g ⇒ ε = mg/c1, len = l0·(1+ε).
        strain = weight / c1
        expected_z = -rest_length * (1 + strain)
        @info "F(ε)=c1·ε balances mg ⇒ len=l0(1+mg/c1)." measured=node.pos_w[3] expected=expected_z
        @test node.pos_w[3] ≈ expected_z rtol=1e-4
        @test abs(node.pos_w[1]) < 1e-6 && abs(node.pos_w[2]) < 1e-6
    end

    @testset "Route 2 tether propagation" begin
        tether_len = 3.0
        end_mass = 3.0
        ground2 = Point(:ground, KVec3(0.0, 0.0, 0.0), STATIC)
        end_pt = Point(:mass, KVec3(0.0, 0.0, -tether_len), DYNAMIC;
                       extra_mass=end_mass)
        line = Tether(:line; start_point=:ground, end_point=:mass,
            n_segments=3, unit_stiffness=force_law, unit_damping=200.0,
            diameter=0.005, density=724.0)
        sys_t = SystemStructure("tether_nl_test", set;
            points=[ground2, end_pt], tethers=[line])

        @info "Route 2 tether: callable unit_stiffness fanned out to every segment."
        @test length(sys_t.segments) == 3
        @test all(!(s.unit_stiffness isa Real) for s in sys_t.segments)
        @test all(s.unit_stiffness === force_law for s in sys_t.segments)

        sam_t = SymbolicAWEModel(set, sys_t)
        node_t = sam_t.sys_struct.points[:mass]
        test_init!(sam_t; prn=false)
        settle!(sam_t, node_t)
        # Uniform tension ≈ end weight: each segment strains ε = mg/c1, so the
        # whole line stretches to total_l0·(1+ε) (tiny segment mass perturbs it).
        strain = end_mass * set.g_earth / c1
        expected_z = -tether_len * (1 + strain)
        @info "Hanging tether: total length = l0(1+mg/c1)." measured=node_t.pos_w[3] expected=expected_z
        @test node_t.pos_w[3] ≈ expected_z rtol=5e-3
    end

    rm(tmpdir; recursive=true)
end
nothing
