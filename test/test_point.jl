# SPDX-FileCopyrightText: 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

# test_point.jl - Point mass dynamics tests
#
# Tests point dynamics including drag forces and gravity.
# Verifies:
# 1. No gravity, no wind: point is stationary
# 2. No gravity, with wind: correct drag acceleration
# 3. With gravity, zero area: pure free fall acceleration

using Test
using SymbolicAWEModels
using SymbolicAWEModels: KVec3
using KiteUtils
using LinearAlgebra

# ============================================================================
# YAML Configuration - Point with drag capability
# ============================================================================
const POINT_TEST_YAML = """
##############################
## Point Test System #########
##############################

###########################
## Points #################
###########################
points:
  headers: [name, pos_cad, type, wing_idx, transform_idx, extra_mass, body_frame_damping, world_frame_damping, area, drag_coeff]
  data:
    - [test_point, [0.0, 0.0, 100.0], DYNAMIC, nothing, nothing, 0.69, 0.0, 0.0, 0.1, 0.45]
"""

# YAML for zero-area point (no drag)
const POINT_NO_DRAG_YAML = """
##############################
## Point Test - No Drag ######
##############################

###########################
## Points #################
###########################
points:
  headers: [name, pos_cad, type, wing_idx, transform_idx, extra_mass, body_frame_damping, world_frame_damping, area, drag_coeff]
  data:
    - [test_point, [0.0, 0.0, 100.0], DYNAMIC, nothing, nothing, 0.69, 0.0, 0.0, 0.0, 0.0]
"""

# World frame damping YAML - point with 0.7 N/(m/s) damping
const POINT_WORLD_DAMPING_YAML = """
##############################
## Point Test - World Damping #
##############################

points:
  headers: [name, pos_cad, type, wing_idx, transform_idx, extra_mass, body_frame_damping, world_frame_damping, area, drag_coeff]
  data:
    - [test_point, [0.0, 0.0, 100.0], DYNAMIC, nothing, nothing, 2.3, 0.0, 0.7, 0.0, 0.0]
"""

# Aerodynamic drag YAML - point with area=4.2, cd=0.33
const POINT_AERO_DRAG_YAML = """
##############################
## Point Test - Aero Drag #####
##############################

points:
  headers: [name, pos_cad, type, wing_idx, transform_idx, extra_mass, body_frame_damping, world_frame_damping, area, drag_coeff]
  data:
    - [test_point, [0.0, 0.0, 100.0], DYNAMIC, nothing, nothing, 2.3, 0.0, 0.0, 4.2, 0.33]
"""

# Heavy point YAML for drag test - large mass for slow,
# precise acceleration measurement
const POINT_HEAVY_DRAG_YAML = """
##############################
## Point Test - Heavy Drag ###
##############################

points:
  headers: [name, pos_cad, type, wing_idx, transform_idx,
            extra_mass, body_frame_damping,
            world_frame_damping, area, drag_coeff]
  data:
    - [test_point, [0.0, 0.0, 100.0], DYNAMIC, nothing,
       nothing, 50.0, 0.0, 0.0, 0.5, 1.0]
"""

# High altitude YAML - point at 5000m where air is thinner
const POINT_HIGH_ALTITUDE_YAML = """
##############################
## Point Test - High Altitude #
##############################

points:
  headers: [name, pos_cad, type, wing_idx, transform_idx, extra_mass, body_frame_damping, world_frame_damping, area, drag_coeff]
  data:
    - [test_point, [0.0, 0.0, 5000.0], DYNAMIC, nothing, nothing, 2.3, 0.0, 0.0, 4.2, 0.33]
"""

@testset "Point Tests" begin
    # Write YAML to temp file
    tmpdir = mktempdir()
    yaml_path = joinpath(tmpdir, "test_point_geometry.yaml")
    write(yaml_path, POINT_TEST_YAML)

    yaml_no_drag_path = joinpath(
        tmpdir, "test_point_no_drag_geometry.yaml")
    write(yaml_no_drag_path, POINT_NO_DRAG_YAML)

    yaml_heavy_drag_path = joinpath(
        tmpdir, "test_point_heavy_drag_geometry.yaml")
    write(yaml_heavy_drag_path, POINT_HEAVY_DRAG_YAML)

    # Create minimal settings file
    settings_yaml = """
system:
    log_file: "data/2plate"  # filename without extension  [replay only]
                                   #   use / as path delimiter, even on Windows 
    g_earth:     9.81

solver:
    solver: "FBDF"
    abs_tol: 0.01          # absolute tolerance of the DAE solver [m, m/s]
    rel_tol: 0.01          # relative tolerance of the DAE solver [-]
    relaxation: 0.6        # relaxation factor of inner linear Newton solver, needed for quasi-steady solver

kite:
    model: ""     # 3D model of the kite
    foil_file: "ram_air_kite/ram_air_kite_foil.dat" # filename for the foil shape
    physical_model: "2plate"            # name of the kite model to use (2plate, ram, etc.)
    struc_geometry_path: "refine_struc_geometry.yaml"  # structural YAML
    aero_geometry_path: "aero_geometry.yaml"    # aerodynamic YAML
    mass: 0.0                               # kite mass [kg]
    quasi_static: false                     # whether to use quasi static kite points or not

tether:
    cd_tether: 0.958
    unit_damping: 0.0
    unit_stiffness: 0.0
    rho_tether: 724.0
    e_tether: 5.5e10


winch:
    winch_model: "TorqueControlledMachine" # or AsynchMachine
    drum_radius: 0.110    # radius of the drum                              [m]
    gear_ratio: 1.0        # gear ratio of the winch                         [-]   
    inertia_total: 0.024   # total inertia, as seen from the motor/generator [kgm²]
    f_coulomb: 122.0       # coulomb friction                                [N]
    c_vf: 30.6             # coefficient for the viscous friction            [Ns/m]

environment:
    rho_0: 1.225               # air density at sea level               [kg/m^3]
    v_wind: 0.0              # wind speed at reference height         [m/s]
    upwind_dir: -90.0        # upwind direction                       [deg]
    profile_law: 0           # 1=EXP, 2=LOG, 3=EXPLOG, 4=FAST_EXP, 5=FAST_LOG, 6=FAST_EXPLOG
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

    # Load system structure from YAML
    sys = load_sys_struct_from_yaml(yaml_path; system_name="point_test", set=set)

    # ========================================================================
    # YAML Loading Verification
    # ========================================================================
    @testset "YAML Loading Verification" begin
        # Verify points were loaded correctly
        @test length(sys.points) == 1
        @test haskey(sys.points, :test_point)

        # Verify point properties
        test_point = sys.points[:test_point]
        @test test_point.type == SymbolicAWEModels.DYNAMIC
        @test test_point.pos_cad == KVec3(0.0, 0.0, 100.0)
        @test test_point.extra_mass == 0.69
        @test test_point.area == 0.1
        @test test_point.drag_coeff == 0.45

        println("\n  ====== Loaded point: mass=$(test_point.extra_mass)kg, area=$(test_point.area)m², Cd=$(test_point.drag_coeff) ======\n")
    end

    # ========================================================================
    # Physics Test 1: No gravity, no wind, zero area - stationary
    # ========================================================================
    @testset "No gravity, no wind - stationary" begin
        set.g_earth = 0.0
        set.v_wind = 0.0

        # Use no-drag YAML for cleaner test
        sys = load_sys_struct_from_yaml(yaml_no_drag_path; system_name="point_test_stat", set=set)
        @test sys.points[:test_point].area == 0.0

        sam = SymbolicAWEModel(set, sys)
        init!(sam; remake=true)

        # Record initial state
        initial_pos = copy(sam.sys_struct.points[:test_point].pos_w)
        initial_vel = copy(sam.sys_struct.points[:test_point].vel_w)

        # Run simulation
        dt = 0.01
        for _ in 1:100
            next_step!(sam; dt=dt, vsm_interval=0)
        end

        final_pos = sam.sys_struct.points[:test_point].pos_w
        final_vel = sam.sys_struct.points[:test_point].vel_w

        # Position should be unchanged (within numerical tolerance)
        @test norm(final_pos - initial_pos) ≈ 0.0 atol=1e-8

        # Velocity should remain near zero
        @test norm(final_vel) < 0.001  # Less than 1mm/s

        println("\n  ====== Position drift: $(round(norm(final_pos - initial_pos)*1000, digits=3)) mm (limit: 0 mm) ======\n")
    end

    # ========================================================================
    # Physics Test 2: No gravity, with wind - drag acceleration
    # ========================================================================
    @testset "No gravity, with wind - drag acceleration" begin
        set.g_earth = 0.0
        set.v_wind = 15.0  # 15 m/s wind
        set.profile_law = 0  # Constant wind profile

        # Heavy point (50 kg) with large drag area for
        # slow, precisely measurable acceleration
        sys = load_sys_struct_from_yaml(
            yaml_heavy_drag_path;
            system_name="point_test_drag", set=set)

        @test sys.points[:test_point].area == 0.5
        @test sys.points[:test_point].drag_coeff == 1.0
        @test sys.points[:test_point].extra_mass == 50.0

        sam = SymbolicAWEModel(set, sys)
        init!(sam; remake=true)

        # Physics: F = 0.5 * rho * Cd * A * v^2
        rho = 1.225
        Cd = 1.0
        A = 0.5
        m = 50.0
        v_wind = 15.0

        F_drag = 0.5 * rho * Cd * A * v_wind^2
        a_expected = F_drag / m
        # a = 0.5*1.225*1.0*0.5*225 / 50 = 1.378 m/s^2

        # Free-floating point: wind drag accelerates it.
        # Heavy mass keeps velocity low so drag stays
        # nearly constant over the measurement window.
        dt = 0.001
        n_steps = 10
        initial_vel = copy(
            sam.sys_struct.points[:test_point].vel_w)

        for _ in 1:n_steps
            next_step!(sam; dt=dt, vsm_interval=0)
        end

        final_vel =
            sam.sys_struct.points[:test_point].vel_w
        t_elapsed = n_steps * dt
        speed = norm(final_vel - initial_vel)
        a_measured = speed / t_elapsed

        @test a_measured ≈ a_expected rtol=0.10

        println(
            "\n  ====== Wind drag: " *
            "a_expected=" *
            "$(round(a_expected, digits=2)), " *
            "a_measured=" *
            "$(round(a_measured, digits=2))" *
            " m/s² ======\n"
        )
    end


    # ========================================================================
    # Physics Test 5: Free fall acceleration (no damping, no drag)
    # ========================================================================
    @testset "Free fall - gravity acceleration" begin
        set.g_earth = 9.81
        set.v_wind = 0.0

        # Use no-drag YAML - just a point with mass, no segments
        sys = load_sys_struct_from_yaml(yaml_no_drag_path; system_name="freefall_test", set=set)

        @test sys.points[:test_point].area == 0.0
        @test sys.points[:test_point].world_frame_damping == KVec3(0.0, 0.0, 0.0)
        @test sys.points[:test_point].extra_mass == 0.69

        sam = SymbolicAWEModel(set, sys)
        init!(sam; remake=true)

        # Initial state
        z0 = sam.sys_struct.points[:test_point].pos_w[3]
        vz0 = sam.sys_struct.points[:test_point].vel_w[3]

        # Run for short time
        dt = 0.001
        n_steps = 100  # 0.1 seconds
        total_time = n_steps * dt

        for _ in 1:n_steps
            next_step!(sam; dt=dt, vsm_interval=0)
        end

        vz_final = sam.sys_struct.points[:test_point].vel_w[3]
        z_final = sam.sys_struct.points[:test_point].pos_w[3]

        # Velocity should be v0 - g*t (mass doesn't affect free fall acceleration)
        vz_expected = vz0 - set.g_earth * total_time
        @test vz_final ≈ vz_expected rtol=0.001

        # Position should follow free fall equation: z = z0 + v0*t - 0.5*g*t^2
        z_expected = z0 + vz0 * total_time - 0.5 * set.g_earth * total_time^2
        @test z_final ≈ z_expected atol=0.001

        println("\n  ====== Free fall: v=$(round(vz_final, digits=3)) m/s, expected=$(round(vz_expected, digits=3)) m/s ======\n")
    end

    # ========================================================================
    # Physics Test 6: World frame damping - terminal velocity
    # ========================================================================
    @testset "World frame damping - terminal velocity" begin
        set.g_earth = 9.81
        set.v_wind = 0.0

        yaml_damping_path = joinpath(tmpdir, "world_damping_geometry.yaml")
        write(yaml_damping_path, POINT_WORLD_DAMPING_YAML)

        sys = load_sys_struct_from_yaml(yaml_damping_path; system_name="damping_test", set=set)

        @test sys.points[:test_point].world_frame_damping == KVec3(0.7, 0.7, 0.7)
        @test sys.points[:test_point].area == 0.0
        @test sys.points[:test_point].extra_mass == 2.3

        sam = SymbolicAWEModel(set, sys)
        init!(sam; remake=true)

        # Terminal velocity calculation:
        # world_frame_damping is a per-mass damping coefficient (damping ratio)
        # F_damping = damping * m * v, so at terminal velocity: damping * m * v = m * g
        # v_terminal = g / damping = 9.81 / 0.7 = 14.01 m/s
        damping = 0.7
        v_terminal_expected = set.g_earth / damping

        # Run until terminal velocity is reached (let it settle)
        dt = 0.01
        n_steps = 2000  # 20 seconds should be enough

        for _ in 1:n_steps
            next_step!(sam; dt=dt, vsm_interval=0)
        end

        # Check terminal velocity (downward, so negative z velocity)
        vz_final = sam.sys_struct.points[:test_point].vel_w[3]
        @test abs(vz_final) ≈ v_terminal_expected atol=0.001
        @test vz_final < 0  # Moving downward

        println("\n  ====== World damping terminal velocity: measured=$(round(abs(vz_final), digits=2)) m/s, expected=$(round(v_terminal_expected, digits=2)) m/s ======\n")
    end

    # ========================================================================
    # Physics Test 7: Aerodynamic drag - terminal velocity
    # ========================================================================
    @testset "Aerodynamic drag - terminal velocity" begin
        set.g_earth = 9.81
        set.v_wind = 0.0

        yaml_aero_path = joinpath(tmpdir, "aero_drag_geometry.yaml")
        write(yaml_aero_path, POINT_AERO_DRAG_YAML)

        sys = load_sys_struct_from_yaml(yaml_aero_path; system_name="aero_drag_test", set=set)

        @test sys.points[:test_point].area == 4.2
        @test sys.points[:test_point].drag_coeff == 0.33
        @test sys.points[:test_point].world_frame_damping == KVec3(0.0, 0.0, 0.0)
        @test sys.points[:test_point].extra_mass == 2.3
        @test sys.points[:test_point].total_mass == 0.0 # Not initialized yet

        sam = SymbolicAWEModel(set, sys)
        init!(sam; remake=true)
        @test sys.points[:test_point].total_mass == 2.3

        # Verify point properties are preserved after init
        point = sam.sys_struct.points[:test_point]
        @test point.area == 4.2
        @test point.drag_coeff == 0.33

        m = 2.3
        Cd = 0.33
        A = 4.2

        # Run until terminal velocity is reached at lower altitude
        dt = 0.1
        n_steps = 1000  # 100 seconds to fall to lower altitude with higher air density

        for _ in 1:n_steps
            next_step!(sam; dt=dt, vsm_interval=0)
        end

        # Calculate air density at final height using same model as simulation
        # Note: simulation clamps height to max(0.0, h) for density calculation
        final_height = point.pos_w[3]
        rho = SymbolicAWEModels.calc_rho(sam.am, max(0.0, final_height))

        # Terminal velocity: 0.5 * rho * Cd * A * v^2 = m * g
        # v = sqrt(2 * m * g / (rho * Cd * A))
        v_terminal_expected = sqrt(2 * m * set.g_earth / (rho * Cd * A))

        # Check terminal velocity (downward, so negative z velocity)
        vz_final = point.vel_w[3]
        @test vz_final < 0  # Moving downward

        println("\n  ====== Aero drag terminal velocity: measured=$(round(abs(vz_final), digits=2)) m/s, expected=$(round(v_terminal_expected, digits=2)) m/s (h=$(round(final_height, digits=0))m) ======\n")
        @test abs(vz_final) ≈ v_terminal_expected rtol=0.01
    end

    # ========================================================================
    # Physics Test 8: High altitude point - lower air density = faster fall
    # Demonstrates that rho is NOT clamped for positive heights
    # ========================================================================
    @testset "High altitude point - thin air faster fall" begin
        set.g_earth = 9.81
        set.v_wind = 0.0
        set.profile_law = 0  # Reset to constant profile

        yaml_high_alt_path = joinpath(tmpdir, "high_alt_geometry.yaml")
        write(yaml_high_alt_path, POINT_HIGH_ALTITUDE_YAML)

        sys = load_sys_struct_from_yaml(yaml_high_alt_path; system_name="high_alt_point_test", set=set)

        @test sys.points[:test_point].area == 4.2
        @test sys.points[:test_point].drag_coeff == 0.33
        @test sys.points[:test_point].extra_mass == 2.3

        sam = SymbolicAWEModel(set, sys)
        init!(sam; remake=true)

        point = sam.sys_struct.points[:test_point]
        m = 2.3
        Cd = 0.33
        A = 4.2

        # Run simulation - shorter time to stay at altitude
        dt = 0.1
        n_steps = 200  # 20 seconds

        for _ in 1:n_steps
            next_step!(sam; dt=dt, vsm_interval=0)
        end

        # Get final state
        final_height = point.pos_w[3]

        # Height should still be positive (above sea level)
        @test final_height > 0

        # Calculate air density at altitude (no clamping needed since h > 0)
        rho_at_altitude = SymbolicAWEModels.calc_rho(sam.am, final_height)
        rho_at_sea_level = SymbolicAWEModels.calc_rho(sam.am, 0.0)

        # Air should be thinner at altitude
        @test rho_at_altitude < rho_at_sea_level

        # Terminal velocity at altitude vs sea level
        v_terminal_altitude = sqrt(2 * m * set.g_earth / (rho_at_altitude * Cd * A))
        v_terminal_sea_level = sqrt(2 * m * set.g_earth / (rho_at_sea_level * Cd * A))

        # Should fall faster at altitude due to thinner air
        @test v_terminal_altitude > v_terminal_sea_level

        # Check actual terminal velocity matches expected
        vz_final = point.vel_w[3]

        @test vz_final < 0  # Moving downward
        @test abs(vz_final) ≈ v_terminal_altitude rtol=0.01

        println("\n  ====== High altitude (h=$(round(final_height, digits=0))m): v=$(round(abs(vz_final), digits=2)) m/s")
        println("  ====== Sea level would be: v=$(round(v_terminal_sea_level, digits=2)) m/s ($(round((v_terminal_altitude/v_terminal_sea_level - 1)*100, digits=1))% faster at altitude) ======\n")
    end

    # Cleanup
    rm(tmpdir; recursive=true)
end
