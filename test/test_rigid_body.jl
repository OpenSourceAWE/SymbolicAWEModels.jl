# SPDX-FileCopyrightText: 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_rigid_body.jl - Standalone RigidBody component tests.
#
# A RigidBody shares the 6-DOF generator (rigid_body_eqs!) with rigid wings but
# carries no aero/tether/transform machinery. Tests:
# 1. Free fall: translational acc = g, no spurious rotation.
# 2. Spin-up under applied body-frame moment: α = τ / I.
# 3. COM-offset body-frame output: com_w = pos_w + R_b_to_w * com_offset_b.

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
    log_file: "data/rigid_body_test"
    g_earth: 9.81
solver:
    solver: "FBDF"
    abs_tol: 0.0001
    rel_tol: 0.0001
kite:
    model: ""
    foil_file: "ram_air_kite/ram_air_kite_foil.dat"
    physical_model: "rigid_body_test"
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

@testset "RigidBody component" begin
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
    body = RigidBody(:body1; mass=2.0, inertia_principal=inertia,
                     pos=[0.0, 0.0, 10.0])
    sys = SystemStructure("rigid_body_test", set; rigid_bodies=[body])

    @testset "Model setup" begin
        @test length(sys.rigid_bodies) == 1
        @test sys.rigid_bodies[:body1].mass ≈ 2.0
    end

    sam = SymbolicAWEModel(set, sys)
    test_init!(sam)

    @testset "Free fall acceleration = g" begin
        rb = sam.sys_struct.rigid_bodies[:body1]
        dt = 0.01
        for _ in 1:5
            next_step!(sam; dt, vsm_interval=0)
        end
        vel_before = copy(rb.vel_w)
        t_before = sam.integrator.t
        for _ in 1:10
            next_step!(sam; dt, vsm_interval=0)
        end
        elapsed = sam.integrator.t - t_before
        acc = (rb.vel_w - vel_before) / elapsed
        @test acc[3] ≈ -9.81 atol=0.1
        @test abs(acc[1]) < 0.1
        @test abs(acc[2]) < 0.1
        @test norm(rb.ω_b) < 1e-6   # no spurious rotation
    end

    @testset "Spin-up under applied moment" begin
        rb = sam.sys_struct.rigid_bodies[:body1]
        rb.pos_w .= [0.0, 0.0, 10.0]
        rb.vel_w .= 0.0
        rb.ω_b .= 0.0
        rb.Q_b_to_w .= [1.0, 0.0, 0.0, 0.0]
        torque = 0.05
        rb.ext_moment_b .= [torque, 0.0, 0.0]
        test_init!(sam; prn=false, reset_vel=false)

        dt = 0.005
        t_before = sam.integrator.t
        ω_before = copy(rb.ω_b)
        for _ in 1:20
            next_step!(sam; dt, vsm_interval=0)
        end
        elapsed = sam.integrator.t - t_before
        α_measured = (rb.ω_b[1] - ω_before[1]) / elapsed
        @test α_measured ≈ torque / inertia[1] atol=0.02
        @test abs(rb.ω_b[2]) < 1e-3
        @test abs(rb.ω_b[3]) < 1e-3
        rb.ext_moment_b .= 0.0
    end

    @testset "COM-offset body-frame output" begin
        offset = [0.5, 0.0, 0.0]
        body2 = RigidBody(:body2; mass=1.0, inertia_principal=inertia,
                          pos=[1.0, 2.0, 5.0], com_offset_b=offset)
        sys2 = SystemStructure("rigid_body_test", set; rigid_bodies=[body2])
        sam2 = SymbolicAWEModel(set, sys2)
        test_init!(sam2)
        rb = sam2.sys_struct.rigid_bodies[:body2]
        # identity orientation: com_w = pos_w + offset
        @test rb.com_w ≈ rb.pos_w .+ offset atol=1e-6
        for _ in 1:10
            next_step!(sam2; dt=0.01, vsm_interval=0)
        end
        R_b_to_w = SymbolicAWEModels.quaternion_to_rotation_matrix(rb.Q_b_to_w)
        @test rb.com_w ≈ rb.pos_w .+ R_b_to_w * offset atol=1e-5
        @test norm(rb.ω_b) < 1e-6
    end

    rm(tmpdir; recursive=true)
end
nothing
