# SPDX-FileCopyrightText: 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_timoshenko_joint.jl - Validate the corotational Timoshenko joint (the
# element of a 2-node Timoshenko beam) against closed-form beam theory. One joint
# connects node A (clamped, STATIC) to node B (free, DYNAMIC); a constant external
# load is applied and the settled equilibrium is compared to the analytic result.
#
# 1. Transverse tip load: δ = PL³/3EI + PL/(kGA). The shear term (second) is what
#    distinguishes Timoshenko from Euler-Bernoulli, so we also assert the measured
#    deflection exceeds the Euler-only value.
# 2. Axial load: δ = PL/EA.
# 3. Torsion moment: φ = TL/GJ.
# 4. Nonlinear rigidities: callable EA(ε)/EIy(κ) softening laws, compared to the
#    self-consistent closed-form equilibrium (proves the strain/curvature argument
#    is passed and the effective-rigidity convention is right).

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
    log_file: "data/timoshenko_test"
    g_earth: 0.0
solver:
    solver: "FBDF"
    abs_tol: 1.0e-8
    rel_tol: 1.0e-8
kite:
    model: ""
    foil_file: "ram_air_kite/ram_air_kite_foil.dat"
    physical_model: "timoshenko_test"
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

# Step until the free body is at rest (static equilibrium), so the comparison is
# limited by the element model, not by an unconverged transient.
function settle!(sam, body; dt=0.01, max_steps=6000, vtol=1e-8)
    for _ in 1:max_steps
        next_step!(sam; dt, vsm_interval=0)
        (norm(body.vel_w) < vtol && norm(body.ω_b) < vtol) && break
    end
end

@testset "Timoshenko joint element" begin
    pkg_root = dirname(@__DIR__)
    tmpdir = mktempdir()
    data_path = joinpath(tmpdir, "2plate_kite")
    cp(joinpath(pkg_root, "data", "2plate_kite"), data_path; force=true)
    write(joinpath(data_path, "settings.yaml"), SETTINGS_YAML)
    write(joinpath(data_path, "system.yaml"),
        "system:\n  sim_settings: settings.yaml\n")
    set_data_path(data_path)
    set = Settings("system.yaml")

    beam_length = 1.0
    EI = 100.0; GA = 1500.0; EA = 1.0e4; GJ = 50.0; kshear = 5 / 6
    inertia = [0.01, 0.01, 0.01]

    nodeA = Body(:nodeA; mass=1.0, inertia_principal=inertia,
                 pos=[0.0, 0.0, 0.0], type=STATIC)
    nodeB = Body(:nodeB; mass=1.0, inertia_principal=inertia,
                 pos=[beam_length, 0.0, 0.0])
    joint = TimoshenkoJoint(:joint, :nodeA, :nodeB;
        EA, GA, GJ, EIy=EI, EIz=EI, shear_coeff=kshear,
        damping_trans=200.0, damping_rot=3.0)
    sys = SystemStructure("timoshenko_test", set;
        bodies=[nodeA, nodeB], timoshenko_joints=[joint])

    @testset "Model setup" begin
        @info "Structure wiring: one joint resolved, rest length from geometry."
        @test length(sys.timoshenko_joints) == 1
        @test sys.timoshenko_joints[:joint].body_a_idx == 1
        @test sys.timoshenko_joints[:joint].body_b_idx == 2
        @test sys.timoshenko_joints[:joint].rest_length ≈ beam_length
    end

    sam = SymbolicAWEModel(set, sys)
    rb = sam.sys_struct.bodies[:nodeB]
    ra = sam.sys_struct.bodies[:nodeA]

    @testset "Transverse tip load (Timoshenko vs Euler)" begin
        load = 5.0
        rb.ext_force_w .= [0.0, 0.0, load]
        rb.ext_moment_b .= 0.0
        test_init!(sam; prn=false)
        settle!(sam, rb)
        timoshenko = load * beam_length^3 / (3 * EI) +
                     load * beam_length / (kshear * GA)
        euler = load * beam_length^3 / (3 * EI)
        @info "Cantilever transverse: δ = PL³/3EI + PL/(kGA), shear term active." measured=rb.pos_w[3] expected=timoshenko
        @test rb.pos_w[3] ≈ timoshenko rtol=0.002  # floor: corotational nonlinearity
        @test rb.pos_w[3] > euler + 0.5 * (timoshenko - euler)  # shear active
        @test norm(ra.pos_w) < 1e-9                              # clamped end fixed
    end

    @testset "Axial load" begin
        load = 20.0
        rb.ext_force_w .= [load, 0.0, 0.0]
        rb.ext_moment_b .= 0.0
        test_init!(sam; prn=false)
        settle!(sam, rb)
        expected = load * beam_length / EA
        @info "Axial stretch: δ = PL/EA." measured=(rb.pos_w[1] - beam_length) expected=expected
        @test rb.pos_w[1] - beam_length ≈ expected rtol=1e-4
        @test abs(rb.pos_w[3]) < 1e-6
    end

    @testset "Torsion moment" begin
        torque = 1.0
        rb.ext_force_w .= 0.0
        rb.ext_moment_b .= [torque, 0.0, 0.0]
        test_init!(sam; prn=false)
        settle!(sam, rb)
        R_b_to_w = SymbolicAWEModels.quaternion_to_rotation_matrix(rb.Q_b_to_w)
        twist = asin(clamp(R_b_to_w[3, 2], -1.0, 1.0))
        @info "Torsion: φ = TL/GJ." measured=twist expected=(torque * beam_length / GJ)
        @test twist ≈ torque * beam_length / GJ rtol=5e-4
        rb.ext_moment_b .= 0.0
    end

    # Nonlinear rigidities: each callable returns the effective rigidity at the
    # current strain/curvature (the inflated-tube/Breukels convention). Softening
    # laws EA(ε)=EA0-aε, EIy(κ)=EI0-bκ make the equilibrium self-consistent, which
    # is solvable in closed form so the comparison stays exact.
    EA0 = 1.0e4; axial_slope = 4.0e5
    EI0 = 100.0; bend_slope = 400.0
    eps_knots = collect(-0.02:0.001:0.02)
    kappa_knots = collect(-0.12:0.005:0.12)
    EA_law = SymbolicAWEModels.LinearInterpolation(
        EA0 .- axial_slope .* abs.(eps_knots), eps_knots)
    EIy_law = SymbolicAWEModels.LinearInterpolation(
        EI0 .- bend_slope .* abs.(kappa_knots), kappa_knots)

    nodeA2 = Body(:nodeA; mass=1.0, inertia_principal=inertia,
                  pos=[0.0, 0.0, 0.0], type=STATIC)
    nodeB2 = Body(:nodeB; mass=1.0, inertia_principal=inertia,
                  pos=[beam_length, 0.0, 0.0])
    joint_nl = TimoshenkoJoint(:joint, :nodeA, :nodeB;
        EA=EA_law, GA, GJ, EIy=EIy_law, EIz=EI, shear_coeff=kshear,
        damping_trans=200.0, damping_rot=3.0)
    sys_nl = SystemStructure("timoshenko_test", set;
        bodies=[nodeA2, nodeB2], timoshenko_joints=[joint_nl])
    sam_nl = SymbolicAWEModel(set, sys_nl)
    rb_nl = sam_nl.sys_struct.bodies[:nodeB]

    @testset "Nonlinear axial (callable rigidity)" begin
        @test !(sys_nl.timoshenko_joints[:joint].EA isa Real)  # callable path
        load = 20.0
        rb_nl.ext_force_w .= [load, 0.0, 0.0]
        rb_nl.ext_moment_b .= 0.0
        test_init!(sam_nl; prn=false)
        settle!(sam_nl, rb_nl)
        # P = EA_eff·ε with EA_eff = EA0 - axial_slope·ε ⇒ quadratic, larger root.
        EA_eff = (EA0 + sqrt(EA0^2 - 4 * axial_slope * load)) / 2
        expected = load / EA_eff
        @info "Softening axial: P=EA(ε)·ε self-consistent." measured=(rb_nl.pos_w[1] - beam_length) expected=expected
        @test rb_nl.pos_w[1] - beam_length ≈ expected rtol=1e-3
        @test expected > load / EA0                              # softer than linear
    end

    @testset "Nonlinear bending (callable rigidity)" begin
        load = 5.0
        rb_nl.ext_force_w .= [0.0, 0.0, load]
        rb_nl.ext_moment_b .= 0.0
        test_init!(sam_nl; prn=false)
        settle!(sam_nl, rb_nl)
        # κ = PL/(2·EI_eff), EI_eff = EI0 - bend_slope·κ ⇒ quadratic, larger root.
        EI_eff = (EI0 + sqrt(EI0^2 - 2 * bend_slope * load * beam_length)) / 2
        expected = load * beam_length^3 / (3 * EI_eff) +
                   load * beam_length / (kshear * GA)
        linear = load * beam_length^3 / (3 * EI0) +
                 load * beam_length / (kshear * GA)
        @info "Softening bending: κ=PL/2EI(κ) self-consistent." measured=rb_nl.pos_w[3] expected=expected
        @test rb_nl.pos_w[3] ≈ expected rtol=5e-3
        @test rb_nl.pos_w[3] > linear                            # softer than linear
    end

    rm(tmpdir; recursive=true)
end
nothing
