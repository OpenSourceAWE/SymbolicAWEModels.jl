# SPDX-FileCopyrightText: 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_profile_law.jl - Wind profile law tests
#
# Tests wind profile calculations at different altitudes using STATIC points as "probes".
# Verifies that wind_at_point matches expected wind speed for each profile law.
#
# Profile laws (see helpers.jl calc_wind_factor):
#   0 = CONST (constant wind, factor = 1.0)
#   1 = EXP (delegated to AtmosphericModels)
#   2 = LOG (delegated to AtmosphericModels)
#   3 = EXPLOG (delegated to AtmosphericModels)

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

using Test
using SymbolicAWEModels
using SymbolicAWEModels: KVec3
using KiteUtils
using LinearAlgebra
using SymbolicIndexingInterface: getu

# ============================================================================
# YAML Configuration - Static probes at different heights
# Need at least one DYNAMIC point with a segment for valid ODE system
# ============================================================================
PROFILE_LAW_TEST_YAML = """
##############################
## Wind Profile Probe System #
##############################

points:
  headers: [name, pos_cad, type, wing_idx, transform_idx, extra_mass, body_frame_damping, world_frame_damping, area, drag_coeff]
  data:
    # Static probes at different heights - measure wind at each location
    - [probe_10m, [0.0, 0.0, 10.0], STATIC, nothing, nothing, 0.0, 0.0, 0.0, 0.0, 0.0]
    - [probe_50m, [0.0, 0.0, 50.0], STATIC, nothing, nothing, 0.0, 0.0, 0.0, 0.0, 0.0]
    - [probe_100m, [0.0, 0.0, 100.0], STATIC, nothing, nothing, 0.0, 0.0, 0.0, 0.0, 0.0]
    - [probe_200m, [0.0, 0.0, 200.0], STATIC, nothing, nothing, 0.0, 0.0, 0.0, 0.0, 0.0]
    - [probe_500m, [0.0, 0.0, 500.0], STATIC, nothing, nothing, 0.0, 0.0, 0.0, 0.0, 0.0]
    # Dynamic point to ensure valid ODE system (connected by segment)
    - [dynamic_point, [0.0, 0.0, 10.5], DYNAMIC, nothing, nothing, 1.0, 0.0, 0.0, 0.0, 0.0]

segments:
  headers: [name, point_i, point_j, l0, diameter_mm, unit_stiffness, unit_damping, compression_frac]
  data:
    # Connect dynamic point to 10m probe with very stiff segment
    - [anchor_segment, probe_10m, dynamic_point, 0.5, 1.0, 1000000.0, 1000.0, 0.1]
"""

@testset "Profile Law Tests" begin
    # Write YAML to temp file
    tmpdir = mktempdir()
    yaml_path = joinpath(tmpdir, "test_profile_law_geometry.yaml")
    write(yaml_path, PROFILE_LAW_TEST_YAML)

    # Create minimal settings file with wind profile parameters
    settings_yaml = """
system:
    log_file: "data/profile_law_test"
    g_earth: 0.0   # No gravity for wind-only tests

solver:
    solver: "FBDF"
    abs_tol: 0.001
    rel_tol: 0.001
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
    v_wind: 10.0           # Reference wind speed [m/s]
    upwind_dir: -90.0      # Wind blows in +x direction
    upwind_elevation: 0.0
    wind_vec: [10.0, 0.0, 0.0]
    profile_law: 0         # Will be overridden per test
    h_ref: 6.0             # Reference height for wind profile [m]
    alpha: 0.08163         # Exponent of wind profile law
    z0: 0.0002             # Surface roughness [m]
"""
    settings_path = joinpath(tmpdir, "settings.yaml")
    write(settings_path, settings_yaml)

    system_yaml = """
system:
  sim_settings: settings.yaml
"""
    system_path = joinpath(tmpdir, "system.yaml")
    write(system_path, system_yaml)

    # Set data path and load settings
    set_data_path(tmpdir)
    set = Settings("system.yaml")

    # Load system structure to verify YAML loading
    sys = load_sys_struct_from_yaml(yaml_path; system_name="profile_law_test", set=set)

    # ========================================================================
    # YAML Loading Verification
    # ========================================================================
    @testset "YAML Loading Verification" begin
        @test length(sys.points) == 6  # 5 probes + 1 dynamic
        @test haskey(sys.points, :probe_10m)
        @test haskey(sys.points, :probe_50m)
        @test haskey(sys.points, :probe_100m)
        @test haskey(sys.points, :probe_200m)
        @test haskey(sys.points, :probe_500m)
        @test haskey(sys.points, :dynamic_point)

        # Verify probe types and positions
        @test sys.points[:probe_10m].type == SymbolicAWEModels.STATIC
        @test sys.points[:probe_10m].pos_cad == KVec3(0.0, 0.0, 10.0)
        @test sys.points[:probe_500m].pos_cad == KVec3(0.0, 0.0, 500.0)
        @test sys.points[:dynamic_point].type == SymbolicAWEModels.DYNAMIC

        println("\n  ====== Loaded 5 static probes at heights: 10m, 50m, 100m, 200m, 500m ======\n")
    end

    # Helper to get wind speed at a point
    function get_wind_speed_at_point(sam, point_name)
        point = sam.sys_struct.points[point_name]
        # Create getter for wind_at_point for this specific point
        wind_getter = getu(sam.prob.sys, sam.prob.sys.wind_at_point[:, point.idx])
        wind_vec = wind_getter(sam.integrator)
        return norm(wind_vec)
    end

    # Helper to get wind vector at a point
    function get_wind_vector_at_point(sam, point_name)
        point = sam.sys_struct.points[point_name]
        wind_getter = getu(sam.prob.sys, sam.prob.sys.wind_at_point[:, point.idx])
        return wind_getter(sam.integrator)
    end

    # ========================================================================
    # Test 1: Profile law 0 (CONST) - uniform wind at all heights
    # ========================================================================
    @testset "Profile law 0 (CONST) - uniform wind" begin
        set.profile_law = 0
        set.v_wind = 10.0

        sys = load_sys_struct_from_yaml(yaml_path; system_name="profile_const", set=set)
        sam = SymbolicAWEModel(set, sys)
        init!(sam; remake=true)

        # All probes should see the same wind speed
        probes = [:probe_10m, :probe_50m, :probe_100m, :probe_200m, :probe_500m]
        for probe in probes
            wind_speed = get_wind_speed_at_point(sam, probe)
            @test wind_speed ≈ 10.0 atol=0.01
        end

        println("\n  ====== Profile law 0 (CONST): All probes see 10.0 m/s wind ======\n")
    end

    # ========================================================================
    # Test 2: Profile law 1 (EXP) - wind increases with height
    # Delegated to AtmosphericModels - verify general behavior
    # ========================================================================
    @testset "Profile law 1 (EXP) - wind increases with height" begin
        set.profile_law = 1
        set.v_wind = 10.0

        sys = load_sys_struct_from_yaml(yaml_path; system_name="profile_exp", set=set)
        sam = SymbolicAWEModel(set, sys)
        init!(sam; remake=true)

        # Wind should increase with height
        probes_ordered = [:probe_10m, :probe_50m, :probe_100m, :probe_200m, :probe_500m]
        wind_speeds = [get_wind_speed_at_point(sam, p) for p in probes_ordered]

        # Verify monotonic increase
        for i in 1:(length(wind_speeds)-1)
            @test wind_speeds[i+1] > wind_speeds[i]
        end

        println("\n  ====== Profile law 1 (EXP): Wind increases with height")
        println("    Heights: 10m, 50m, 100m, 200m, 500m")
        println("    Wind speeds: $(round.(wind_speeds, digits=2)) m/s ======\n")
    end

    # ========================================================================
    # Test 3: Profile law 2 (LOG) - wind increases with height
    # Delegated to AtmosphericModels - verify general behavior
    # ========================================================================
    @testset "Profile law 2 (LOG) - wind increases with height" begin
        set.profile_law = 2
        set.v_wind = 10.0

        sys = load_sys_struct_from_yaml(yaml_path; system_name="profile_log", set=set)
        sam = SymbolicAWEModel(set, sys)
        init!(sam; remake=true)

        # Wind should increase with height
        probes_ordered = [:probe_10m, :probe_50m, :probe_100m, :probe_200m, :probe_500m]
        wind_speeds = [get_wind_speed_at_point(sam, p) for p in probes_ordered]

        # Verify monotonic increase
        for i in 1:(length(wind_speeds)-1)
            @test wind_speeds[i+1] > wind_speeds[i]
        end

        println("\n  ====== Profile law 2 (LOG): Wind increases with height")
        println("    Heights: 10m, 50m, 100m, 200m, 500m")
        println("    Wind speeds: $(round.(wind_speeds, digits=2)) m/s ======\n")
    end

    # NOTE: Profile law 4 (FAST_EXP) is not supported by
    # AtmosphericModels (only 1=EXP, 2=LOG, 3=EXPLOG).

    # ========================================================================
    # Test 5: Wind direction consistency - all probes same direction
    # ========================================================================
    @testset "Wind direction consistency" begin
        set.profile_law = 0  # Use CONST for uniform wind
        set.v_wind = 10.0
        set.upwind_dir = -90.0  # Wind from +x direction

        sys = load_sys_struct_from_yaml(yaml_path; system_name="profile_direction", set=set)
        sam = SymbolicAWEModel(set, sys)
        init!(sam; remake=true)

        probes = [:probe_10m, :probe_50m, :probe_100m, :probe_200m, :probe_500m]

        println("\n  ====== Wind direction test (upwind_dir=-90 -> wind in +x) ======")
        for probe in probes
            wind_vec = get_wind_vector_at_point(sam, probe)
            wind_speed = norm(wind_vec)

            # Wind should be primarily in +x direction
            @test wind_vec[1] > 0  # Positive x component
            @test abs(wind_vec[2]) < 0.01  # Near-zero y
            @test abs(wind_vec[3]) < 0.01  # Near-zero z

            # Direction should be purely +x (normalized)
            if wind_speed > 0.01
                wind_dir = wind_vec ./ wind_speed
                @test wind_dir[1] ≈ 1.0 atol=0.01
            end
        end
        println("    All probes have wind in +x direction\n")
    end

    # ========================================================================
    # Test 6: Zero wind - all probes see zero apparent velocity
    # ========================================================================
    @testset "Zero wind - no apparent velocity" begin
        set.profile_law = 0
        set.v_wind = 0.0

        sys = load_sys_struct_from_yaml(yaml_path; system_name="profile_zero", set=set)
        sam = SymbolicAWEModel(set, sys)
        init!(sam; remake=true)

        probes = [:probe_10m, :probe_50m, :probe_100m, :probe_200m, :probe_500m]

        for probe in probes
            wind_speed = get_wind_speed_at_point(sam, probe)
            @test wind_speed < 0.001
        end

        println("\n  ====== Zero wind: All probes see 0 m/s ======\n")
    end

    # ========================================================================
    # Test 7: Profile law 3 (EXPLOG) - wind increases with height
    # Delegated to AtmosphericModels - verify general behavior
    # ========================================================================
    @testset "Profile law 3 (EXPLOG) - wind increases with height" begin
        set.profile_law = 3
        set.v_wind = 10.0

        sys = load_sys_struct_from_yaml(yaml_path; system_name="profile_explog", set=set)
        sam = SymbolicAWEModel(set, sys)
        init!(sam; remake=true)

        # Wind should increase monotonically with height
        probes_ordered = [:probe_10m, :probe_50m, :probe_100m, :probe_200m, :probe_500m]
        wind_speeds = [get_wind_speed_at_point(sam, p) for p in probes_ordered]

        for i in 1:(length(wind_speeds)-1)
            @test wind_speeds[i+1] > wind_speeds[i]
        end

        println("\n  ====== Profile law 3 (EXPLOG): Wind increases with height")
        println("    Heights: 10m, 50m, 100m, 200m, 500m")
        println("    Wind speeds: $(round.(wind_speeds, digits=2)) m/s ======\n")
    end

    # Cleanup
    rm(tmpdir; recursive=true)
end
nothing
