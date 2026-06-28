# SPDX-FileCopyrightText: 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_beam_replay.jl - SysLog round-trip for a rigid-body beam.
#
# Exercises the multi-frame orientation logging added on top of KiteUtils'
# `orients`: a chain of RigidBodies + ElasticJoints is logged via SysState,
# saved/loaded as an Arrow SysLog, and reconstructed with
# `update_from_sysstate!` (the path `replay` uses). No GLMakie rendering here.

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

@isdefined(test_init!) || include(joinpath(@__DIR__, "util.jl"))

using Test
using SymbolicAWEModels
using SymbolicAWEModels: update_from_sysstate!
using KiteUtils
using LinearAlgebra

SETTINGS_YAML = """
system: {log_file: "data/beam", g_earth: 9.81}
solver: {solver: "FBDF", abs_tol: 1.0e-7, rel_tol: 1.0e-7}
kite: {model: "", foil_file: "ram_air_kite/ram_air_kite_foil.dat", physical_model: "beam_replay", mass: 0.0}
tether: {cd_tether: 0.958, unit_damping: 0.0, unit_stiffness: 0.0, rho_tether: 724.0, e_tether: 5.5e10}
winch: {winch_model: "TorqueControlledMachine", drum_radius: 0.110, gear_ratio: 1.0, inertia_total: 0.024, f_coulomb: 122.0, c_vf: 30.6}
environment: {rho_0: 1.225, v_wind: 0.0, upwind_dir: -90.0, upwind_elevation: 0.0, wind_vec: [0.0, 0.0, 0.0], profile_law: 0}
"""

@testset "Beam SysLog replay round-trip" begin
    n = 6; L = 0.5; m = 0.5; r = 0.02
    inertia = [0.5*m*r^2, m*L^2/12, m*L^2/12]

    pkg_root = dirname(@__DIR__)
    tmpdir = mktempdir()
    data_path = joinpath(tmpdir, "2plate_kite")
    cp(joinpath(pkg_root, "data", "2plate_kite"), data_path; force=true)
    write(joinpath(data_path, "settings.yaml"), SETTINGS_YAML)
    write(joinpath(data_path, "system.yaml"),
        "system:\n  sim_settings: settings.yaml\n")
    set_data_path(data_path)
    set = Settings("system.yaml")

    make_bodies() = [Body(Symbol("seg_$i"); mass=m,
        inertia_principal=inertia, pos=[(i-0.5)*L, 0.0, 0.0],
        type=(i==1 ? STATIC : DYNAMIC)) for i in 1:n]
    joints = [ElasticJoint(Symbol("j_$i"), Symbol("seg_$i"), Symbol("seg_$(i+1)");
        anchor_a=[L/2,0.0,0.0], anchor_b=[-L/2,0.0,0.0],
        stiffness_axial=1e5, stiffness_shear=1e5, stiffness_torsion=5e3,
        stiffness_bending=5e3, damping_trans=50.0, damping_rot=20.0)
        for i in 1:(n-1)]

    sys = SystemStructure("beam_replay", set;
        bodies=make_bodies(), elastic_joints=joints)
    sam = SymbolicAWEModel(set, sys)
    init!(sam)

    ss = SysState(sam)
    @test typeof(ss) == SysState{n, n}   # n body slots, n orientation frames

    dt = 0.02; nsteps = 80
    logger = Logger(sam, nsteps + 1)
    for step in 0:nsteps
        step > 0 && next_step!(sam; dt, vsm_interval=0)
        update_sys_state!(ss, sam); ss.time = step*dt; log!(logger, ss)
    end
    tip = sam.sys_struct.bodies[:seg_6]
    @test tip.pos_w[3] < -0.01           # beam sagged under gravity

    save_log(logger, "beam_replay_test")
    lg = load_log("beam_replay_test")
    @test length(lg.syslog) == nsteps + 1
    @test length(lg.syslog.Qw[1]) == n   # O frames preserved through Arrow
    @test norm(collect(lg.syslog.orient[end])) ≈ 1.0 atol=1e-3

    # Reconstruct geometry from the last frame (the path replay uses).
    sys2 = SystemStructure("beam_replay", set;
        bodies=make_bodies(), elastic_joints=joints)
    update_from_sysstate!(sys2, lg.syslog[end])
    @test sys2.bodies[:seg_6].pos_w[3] ≈ tip.pos_w[3] atol=1e-3
    @test sys2.bodies[:seg_1].pos_w[1] ≈ 0.25 atol=1e-4  # fixed root

    # On Windows, load_log keeps an Arrow mmap handle open, so the temp dir
    # may still be locked here; eager cleanup is best-effort.
    try
        rm(tmpdir; recursive=true)
    catch err
        err isa Base.IOError || rethrow()
    end
end
nothing
