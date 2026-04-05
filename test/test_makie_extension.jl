# SPDX-FileCopyrightText: 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

# test_makie_extension.jl - Tests for the Makie extension
#
# Verifies:
# 1. Multi-system plot creates vector-typed colors (no crash)
# 2. Single-system record produces output file
# 3. Multi-system record produces output file
# 4. Replay single system
# 5. Replay multiple systems

using Test

# GLMakie requires OpenGL — skip tests on CI runners without GPU drivers
const GLMAKIE_AVAILABLE = try
    @eval using GLMakie
    GLMakie.activate!(; visible=false)
    true
catch e
    @warn "GLMakie not available, skipping Makie extension tests" exception=e
    false
end

if !GLMAKIE_AVAILABLE
    @testset "Makie Extension" begin
        @test true skip=true  # GLMakie unavailable
    end
else

using SymbolicAWEModels
using SymbolicAWEModels: KVec3
using KiteUtils

# ============================================================================
# Minimal 2-point, 1-segment YAML (same pattern as test_segment.jl)
# ============================================================================
MAKIE_TEST_YAML = """
points:
  headers: [name, pos_cad, type, wing_idx, transform_idx,
            extra_mass, body_frame_damping, world_frame_damping,
            area, drag_coeff]
  data:
    - [anchor, [0.0, 0.0, 0.0], STATIC, nothing, nothing,
       0.0, 0.0, 0.0, 0.0, 0.0]
    - [mass_point, [0.0, 0.0, -10.0], DYNAMIC, nothing,
       nothing, 1.0, 0.0, 0.0, 0.0, 0.0]

segments:
  headers: [name, point_i, point_j, l0, diameter_mm,
            unit_stiffness, unit_damping, compression_frac]
  data:
    - [seg1, anchor, mass_point, 10.0, 5.0,
       1000.0, 10.0, 0.1]
"""

SETTINGS_YAML = """
system:
    log_file: "data/makie_test"
    g_earth: 9.81

solver:
    solver: "FBDF"
    abs_tol: 0.0001
    rel_tol: 0.0001
    relaxation: 0.6

kite:
    model: ""
    foil_file: "ram_air_kite/ram_air_kite_foil.dat"
    physical_model: "2plate"
    struc_geometry_path: "refine_struc_geometry.yaml"
    aero_geometry_path: "aero_geometry.yaml"
    mass: 0.0
    quasi_static: false

tether:
    cd_tether: 0.0
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
    profile_law: 0
"""

function build_test_syslog(sam, sys, n_steps, dt)
    logger = Logger(sam, n_steps)
    sys_state = SysState(sam)
    for i in 1:n_steps
        next_step!(sam; dt=dt, vsm_interval=0)
        update_sys_state!(sys_state, sam)
        sys_state.time = i * dt
        log!(logger, sys_state)
    end
    save_log(logger, "makie_test_log")
    return load_log("makie_test_log")
end

@testset "Makie Extension" begin
    tmpdir = mktempdir()
    yaml_path = joinpath(tmpdir, "refine_struc_geometry.yaml")
    write(yaml_path, MAKIE_TEST_YAML)
    settings_path = joinpath(tmpdir, "settings.yaml")
    write(settings_path, SETTINGS_YAML)
    system_yaml = "system:\n  sim_settings: settings.yaml\n"
    write(joinpath(tmpdir, "system.yaml"), system_yaml)

    set_data_path(tmpdir)
    set = Settings("system.yaml")

    # Build SAM once (expensive) and reuse
    sys1 = load_sys_struct_from_yaml(
        yaml_path; system_name="makie_test_1", set=set)
    sam = SymbolicAWEModel(set, sys1)
    init!(sam; remake=true, prn=false)

    # Build SysLog with a few frames
    lg1 = build_test_syslog(sam, sys1, 5, 0.05)

    # Reset and build second log
    init!(sam; prn=false)
    lg2 = build_test_syslog(sam, sys1, 5, 0.05)

    # Second SystemStructure for multi-system tests
    sys2 = load_sys_struct_from_yaml(
        yaml_path; system_name="makie_test_2", set=set)

    # Reset sys structs to first frame for plotting
    update_from_sysstate!(sys1, lg1.syslog[1])
    update_from_sysstate!(sys2, lg2.syslog[1])

    # ================================================================
    # Test 1: Multi-system plot creates vector-typed colors
    # ================================================================
    @testset "Multi-system plot vector colors" begin
        scene = plot([sys1, sys2]; use_observables=true)
        @test scene isa GLMakie.Scene
    end

    # ================================================================
    # Test 2: Single-system record produces output file
    # ================================================================
    @testset "Single-system record" begin
        outfile = joinpath(tmpdir, "single.mp4")
        scene = SymbolicAWEModels.record(
            lg1, sys1, outfile; framerate=10)
        @test scene isa GLMakie.Scene
        @test isfile(outfile)
        @test filesize(outfile) > 0
    end

    # ================================================================
    # Test 3: Multi-system record produces output file
    # ================================================================
    @testset "Multi-system record" begin
        outfile = joinpath(tmpdir, "multi.mp4")
        scene = SymbolicAWEModels.record(
            [lg1, lg2], [sys1, sys2], outfile; framerate=10)
        @test scene isa GLMakie.Scene
        @test isfile(outfile)
        @test filesize(outfile) > 0
    end

    # ================================================================
    # Test 4: Replay single system
    # ================================================================
    @testset "Single-system replay" begin
        scene = replay(lg1, sys1)
        @test scene isa GLMakie.Scene
    end

    # ================================================================
    # Test 5: Replay multiple systems
    # ================================================================
    @testset "Multi-system replay" begin
        scene = replay([lg1, lg2], [sys1, sys2])
        @test scene isa GLMakie.Scene
    end

    rm(tmpdir; recursive=true)
end

end # if GLMAKIE_AVAILABLE
