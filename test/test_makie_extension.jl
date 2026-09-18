# SPDX-FileCopyrightText: 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_makie_extension.jl - Tests for the Makie extension
#
# Verifies:
# 1. Multi-system plot creates vector-typed colors (no crash)
# 2. Single-system record produces output file
# 3. Multi-system record produces output file
# 4. Replay single system
# 5. Replay multiple systems
# 6. The replay spring-force checkbox recolours the segments

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

@isdefined(test_init!) || include(joinpath(@__DIR__, "util.jl"))

using Test

# GLMakie requires OpenGL — skip tests on CI runners without GPU drivers
const GLMAKIE_AVAILABLE = try
    @eval using GLMakie
    @eval using MakieControlPlots
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
# Minimal 3-point, 2-segment YAML (same pattern as test_segment.jl). Two
# segments, because one carries the lower mass and the other carries both, so
# the spring-force colour ramp has a range to span.
# ============================================================================
MAKIE_TEST_YAML = """
points:
  headers: [name, pos_cad, type, wing_idx, transform_idx,
            extra_mass, body_frame_damping, world_frame_damping,
            area, drag_coeff]
  data:
    - [anchor, [0.0, 0.0, 0.0], STATIC, nothing, nothing,
       0.0, 0.0, 0.0, 0.0, 0.0]
    - [mid_point, [0.0, 0.0, -5.0], DYNAMIC, nothing,
       nothing, 1.0, 0.0, 0.0, 0.0, 0.0]
    - [mass_point, [0.0, 0.0, -10.0], DYNAMIC, nothing,
       nothing, 1.0, 0.0, 0.0, 0.0, 0.0]

segments:
  headers: [name, point_i, point_j, l0, diameter_mm,
            unit_stiffness, unit_damping, compression_frac]
  data:
    - [seg1, anchor, mid_point, 5.0, 5.0,
       1000.0, 10.0, 0.1]
    - [seg2, mid_point, mass_point, 5.0, 5.0,
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
    struc_geometry_path: "particle_structural_geometry.yaml"
    aero_geometry_path: "aero_geometry.yaml"
    mass: 0.0

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
    upwind_elevation: 0.0
    wind_vec: [0.0, 0.0, 0.0]
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
    yaml_path = joinpath(tmpdir, "particle_structural_geometry.yaml")
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
    test_init!(sam)

    # Build SysLog with a few frames
    lg1 = build_test_syslog(sam, sys1, 5, 0.05)

    # Reset and build second log
    test_init!(sam; prn=false)
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
        scene = MakieControlPlots.plot([sys1, sys2]; use_observables=true)
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

    # ================================================================
    # Test 6: The replay checkbox colours the segments by spring force
    # ================================================================
    # `sys2` never ran the simulation, so the only forces it can colour by are
    # the ones the replayed log put there.
    @testset "Spring force checkbox recolours the replayed segments" begin
        ext = Base.get_extension(SymbolicAWEModels, :SymbolicAWEModelsMakieExt)
        replay(lg2, sys2)
        colors = ext.PLOT_SEGMENT_COLORS_OBS[]
        @test allequal(colors[])            # plain segment colour by default

        ext.apply_view_toggle!(:segment_colors_obs, true)
        @test ext.PLOT_FORCE_COLOR[]
        @test !allequal(colors[])           # the two segments differ in tension

        ext.apply_view_toggle!(:segment_colors_obs, false)
        @test !ext.PLOT_FORCE_COLOR[]
        @test allequal(colors[])

        # Starting the replay with the box ticked needs no click.
        replay(lg2, sys2; force_color=true)
        @test !allequal(ext.PLOT_SEGMENT_COLORS_OBS[][])
    end

    # ================================================================
    # Test 7: A bare record filename lands in the output folder
    # ================================================================
    @testset "Bare record filename lands in the output folder" begin
        SymbolicAWEModels.record(lg1, sys1, "bare_name.mp4"; framerate=10)
        @test isfile(joinpath(get_output_path(), "bare_name.mp4"))
        @test !isfile("bare_name.mp4")
    end

    # No teardown: load_log mmaps the Arrow file, and Windows locks a mapped file.
end

end # if GLMAKIE_AVAILABLE
nothing
