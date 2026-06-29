# SPDX-FileCopyrightText: 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_joint.jl - 6-DOF elastic joint between two RigidBodies.
#
# Gravity is disabled (g_earth = 0) to isolate the joint dynamics. One model
# (two bodies + one joint) is built once; each test varies stiffness and the
# initial conditions through the live registered accessors.
#
# 1. Axial oscillation: only EA. Relative x oscillates at ω = √(EA·(1/m1+1/m2)).
# 2. Torsional oscillation: only GJ, body A near-fixed. Body B twists at
#    ω = √(GJ/Ixx_B).
# 3. Momentum conservation: COM of the pair stays fixed (no external force).

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

@isdefined(test_init!) || include(joinpath(@__DIR__, "util.jl"))

using Test
using SymbolicAWEModels
using KiteUtils
using LinearAlgebra
using Statistics: mean

SETTINGS_YAML = """
system:
    log_file: "data/joint_test"
    g_earth: 0.0
solver:
    solver: "FBDF"
    abs_tol: 1.0e-8
    rel_tol: 1.0e-8
kite:
    model: ""
    foil_file: "ram_air_kite/ram_air_kite_foil.dat"
    physical_model: "joint_test"
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

"""Period from zero crossings of a signal oscillating about zero."""
function period_from_crossings(times, signal)
    crossings = Int[]
    for i in 2:lastindex(signal)
        signal[i-1] * signal[i] < 0 && push!(crossings, i)
    end
    @assert length(crossings) >= 3 "too few zero crossings"
    half_periods = [times[crossings[i+1]] - times[crossings[i]]
                    for i in 1:(length(crossings)-1)]
    return 2 * mean(half_periods)
end

@testset "ElasticJoint" begin
    pkg_root = dirname(@__DIR__)
    src_data_path = joinpath(pkg_root, "data", "2plate_kite")
    tmpdir = mktempdir()
    data_path = joinpath(tmpdir, "2plate_kite")
    cp(src_data_path, data_path; force=true)
    write(joinpath(data_path, "settings.yaml"), SETTINGS_YAML)
    write(joinpath(data_path, "system.yaml"),
        "system:\n  sim_settings: settings.yaml\n")
    set_data_path(data_path)
    set = Settings("system.yaml")

    inertia = [0.1, 0.2, 0.3]
    body1 = Body(:b1; mass=1.0, inertia_principal=inertia,
                      pos=[0.0, 0.0, 0.0])
    body2 = Body(:b2; mass=1.0, inertia_principal=inertia,
                      pos=[1.0, 0.0, 0.0])
    # Anchors meet at the midpoint [0.5, 0, 0] when relaxed.
    joint = ElasticJoint(:j1, :b1, :b2;
        anchor_a=[0.5, 0.0, 0.0], anchor_b=[-0.5, 0.0, 0.0],
        stiffness_axial=0.0, stiffness_shear=0.0,
        stiffness_torsion=0.0, stiffness_bending=0.0)
    sys = SystemStructure("joint_test", set;
        bodies=[body1, body2], elastic_joints=[joint])

    @testset "Model setup" begin
        @test length(sys.bodies) == 2
        @test length(sys.elastic_joints) == 1
        @test sys.elastic_joints[:j1].body_a_idx == 1
        @test sys.elastic_joints[:j1].body_b_idx == 2
    end

    sam = SymbolicAWEModel(set, sys)
    test_init!(sam)
    b1 = sam.sys_struct.bodies[:b1]
    b2 = sam.sys_struct.bodies[:b2]
    jt = sam.sys_struct.elastic_joints[:j1]

    # Perturb the CAD home (pos_cad / R_b_to_c); init resets pos_w/Q_b_to_w
    # to these, so the stretch/twist survives without a warm-start flag.
    function reset_bodies!()
        for (b, x) in ((b1, 0.0), (b2, 1.0))
            b.pos_cad .= [x, 0.0, 0.0]
            b.R_b_to_c .= [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
            b.vel_w .= 0.0
            b.ω_b .= 0.0
        end
    end

    @testset "Axial oscillation frequency" begin
        jt.stiffness_axial = 100.0
        jt.stiffness_shear = 0.0
        jt.stiffness_torsion = 0.0
        jt.stiffness_bending = 0.0
        reset_bodies!()
        # CAD pose is unstrained, so excite the axial mode with a velocity kick.
        b2.vel_w .= [0.5, 0.0, 0.0]
        test_init!(sam; prn=false, reset_vel=false)

        ω_expected = sqrt(100.0 * (1/1.0 + 1/1.0))
        dt = 0.001
        n_steps = round(Int, 4 * (2pi / ω_expected) / dt)
        times = Float64[]
        rel_x = Float64[]
        for _ in 1:n_steps
            next_step!(sam; dt, vsm_interval=0)
            push!(times, sam.integrator.t)
            push!(rel_x, (b2.pos_w[1] - b1.pos_w[1]) - 1.0)
        end
        T_measured = period_from_crossings(times, rel_x)
        @test 2pi / T_measured ≈ ω_expected rtol=0.02
        # Transverse DOFs stay quiet
        @test abs(b2.pos_w[2]) < 1e-6
        @test abs(b2.pos_w[3]) < 1e-6
    end

    @testset "Momentum conservation (COM constant velocity)" begin
        jt.stiffness_axial = 100.0
        reset_bodies!()
        # No external force: a velocity kick conserves total momentum, so the COM
        # drifts at constant velocity while the bodies oscillate internally.
        b2.vel_w .= [0.3, 0.0, 0.0]
        test_init!(sam; prn=false, reset_vel=false)
        com0 = (b1.pos_w + b2.pos_w) / 2
        vcom = (b1.vel_w + b2.vel_w) / 2
        dt = 0.001
        n_steps = 200
        for _ in 1:n_steps
            next_step!(sam; dt, vsm_interval=0)
        end
        com = (b1.pos_w + b2.pos_w) / 2
        @test com ≈ com0 + vcom * (n_steps * dt) atol=1e-6
    end

    @testset "Torsional oscillation frequency" begin
        # Body A near-fixed via huge inertia: B is a torsional pendulum.
        b1.inertia_principal .= [1.0e4, 1.0e4, 1.0e4]
        Ixx = 0.1
        b2.inertia_principal .= [Ixx, 0.2, 0.3]
        jt.stiffness_axial = 0.0
        jt.stiffness_torsion = 5.0
        jt.stiffness_bending = 0.0
        jt.stiffness_shear = 0.0
        reset_bodies!()
        # CAD orientation is unstrained, so excite the torsional mode with a spin kick.
        b2.ω_b .= [0.5, 0.0, 0.0]
        test_init!(sam; prn=false, reset_vel=false)

        ω_expected = sqrt(5.0 / Ixx)
        dt = 0.001
        n_steps = round(Int, 4 * (2pi / ω_expected) / dt)
        times = Float64[]
        ωx = Float64[]
        for _ in 1:n_steps
            next_step!(sam; dt, vsm_interval=0)
            push!(times, sam.integrator.t)
            push!(ωx, b2.ω_b[1])
        end
        T_measured = period_from_crossings(times, ωx)
        @test 2pi / T_measured ≈ ω_expected rtol=0.03
        b1.inertia_principal .= inertia
        b2.inertia_principal .= inertia
    end

    @testset "Interpolated (nonlinear) stiffness" begin
        # A LinearInterpolation reproducing force = EA·Δ must give the same axial
        # frequency as the float law — exercises the float/interp mix, the
        # zero-alloc function barrier, and the Dual derivative through the interp.
        EA = 100.0
        knots = collect(-0.6:0.05:0.6)
        f_axial = SymbolicAWEModels.LinearInterpolation(EA .* knots, knots)
        b1i = Body(:b1; mass=1.0, inertia_principal=inertia, pos=[0.0, 0.0, 0.0])
        b2i = Body(:b2; mass=1.0, inertia_principal=inertia, pos=[1.0, 0.0, 0.0])
        joint_i = ElasticJoint(:j1, :b1, :b2;
            anchor_a=[0.5, 0.0, 0.0], anchor_b=[-0.5, 0.0, 0.0],
            stiffness_axial=f_axial,       # interpolation ...
            stiffness_shear=0.0, stiffness_torsion=0.0, stiffness_bending=0.0)  # ... mixed with floats
        @test joint_i.stiffness_axial === f_axial
        sys_i = SystemStructure("joint_test", set;
            bodies=[b1i, b2i], elastic_joints=[joint_i])
        sam_i = SymbolicAWEModel(set, sys_i)
        test_init!(sam_i; prn=false)   # zero-alloc RHS with the interpolation

        body1 = sam_i.sys_struct.bodies[:b1]
        body2 = sam_i.sys_struct.bodies[:b2]
        # CAD pose is unstrained, so excite the axial mode with a velocity kick.
        body2.vel_w .= [0.5, 0.0, 0.0]
        test_init!(sam_i; prn=false, reset_vel=false)

        ω_expected = sqrt(EA * (1/1.0 + 1/1.0))
        dt = 0.001
        n_steps = round(Int, 4 * (2pi / ω_expected) / dt)
        times = Float64[]
        rel_x = Float64[]
        for _ in 1:n_steps
            next_step!(sam_i; dt, vsm_interval=0)
            push!(times, sam_i.integrator.t)
            push!(rel_x, (body2.pos_w[1] - body1.pos_w[1]) - 1.0)
        end
        T_measured = period_from_crossings(times, rel_x)
        @test 2pi / T_measured ≈ ω_expected rtol=0.02
    end

    @testset "As-placed geometry is unstrained (zero wrench at init)" begin
        # Body B is both offset (1 m anchor gap) and rotated 30° about z relative
        # to A. With CAD-as-rest capture the joint wrench is zero at init, so B
        # must stay put despite stiff springs; without it the 1 m gap alone would
        # fling it.
        half = π / 12  # half-angle of 30°
        Q_rot = Float64[cos(half), 0.0, 0.0, sin(half)]
        bodyA = Body(:b1; mass=1.0, inertia_principal=inertia, pos=[0.0, 0.0, 0.0])
        bodyB = Body(:b2; mass=1.0, inertia_principal=inertia,
                     pos=[1.0, 0.0, 0.0], Q_b_to_w=Q_rot)
        joint_r = ElasticJoint(:j1, :b1, :b2;
            stiffness_axial=1000.0, stiffness_shear=1000.0,
            stiffness_torsion=1000.0, stiffness_bending=1000.0,
            damping_trans=10.0, damping_rot=10.0)
        sys_r = SystemStructure("joint_test", set;
            bodies=[bodyA, bodyB], elastic_joints=[joint_r])
        sam_r = SymbolicAWEModel(set, sys_r)
        test_init!(sam_r; prn=false)

        jr = sam_r.sys_struct.elastic_joints[:j1]
        Rz30 = [cos(2half) -sin(2half) 0.0; sin(2half) cos(2half) 0.0; 0.0 0.0 1.0]
        @info "Rest captured from CAD: anchor offset and relative rotation."
        @test jr.rest_offset_a ≈ [1.0, 0.0, 0.0] atol=1e-9
        @test norm(jr.R_rel0 - Rz30) < 1e-9

        nodeB = sam_r.sys_struct.bodies[:b2]
        pos0 = copy(nodeB.pos_w)
        for _ in 1:50
            next_step!(sam_r; dt=0.01, vsm_interval=0)
        end
        R_meas = SymbolicAWEModels.quaternion_to_rotation_matrix(nodeB.Q_b_to_w)
        @info "Body B unmoved ⇒ wrench was zero at the placed config." drift=norm(nodeB.pos_w - pos0)
        @test norm(nodeB.pos_w - pos0) < 1e-7
        @test norm(nodeB.vel_w) < 1e-7
        @test norm(R_meas - Rz30) < 1e-6
    end

    rm(tmpdir; recursive=true)
end
nothing
