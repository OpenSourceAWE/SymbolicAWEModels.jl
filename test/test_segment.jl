# SPDX-FileCopyrightText: 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_segment.jl - Spring-damper segment dynamics tests
#
# Verifies:
# 1. No gravity, no wind: point stays still (STATIC + DYNAMIC system)
# 2. With gravity: oscillation frequency, equilibrium, damping ratio (STATIC + DYNAMIC)
# 3. Horizontal segment gravity drag: terminal velocity (two DYNAMIC points)
# 4. Vertical segment wind drag: terminal velocity matches wind (two DYNAMIC points)

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
using Statistics

# ============================================================================
# YAML Configuration - Minimal 2-point system with 1 segment
# Note: unit_stiffness and unit_damping are per-unit-length properties:
#   effective k = unit_stiffness / length [N/m]
#   effective c = unit_damping / length [N·s/m]
# ============================================================================
const SEGMENT_TEST_YAML = """
##############################
## Segment Test System #######
##############################

###########################
## Materials ##############
###########################
materials:
  headers: [name, youngs_modulus, density, damping_per_stiffness]
  data:
    - [test_material, 55000000000.0, 724, 0.00077]

###########################
## Points #################
###########################
points:
  headers: [name, pos_cad, type, wing_idx, transform_idx, extra_mass, body_frame_damping, world_frame_damping, area, drag_coeff]
  data:
    - [anchor, [0.0, 0.0, 0.0], STATIC, nothing, nothing, 0.0, 0.0, 0.0, 0.0, 0.0]
    - [mass_point, [0.0, 0.0, -10.0], DYNAMIC, nothing, nothing, 1.0, 0.0, 0.0, 0.0, 0.0]

###########################
## Segments ###############
###########################
segments:
  headers: [name, point_i, point_j, l0, diameter_mm, unit_stiffness, unit_damping, compression_frac]
  data:
    # unit_stiffness=1000, l0=10 -> k=100 N/m; unit_damping=10, l0=10 -> c=1 N·s/m
    - [test_segment, anchor, mass_point, 10.0, 5.0, 1000.0, 10.0, 0.1]
"""

# YAML for oscillation test
# Uses 5.0 kg point mass (>> segment mass ~0.14 kg) for clean spring-damper dynamics
# Mass hangs below anchor (z=-10), so gravity stretches spring
# Note: effective k = unit_stiffness/l0, effective c = unit_damping/l0
const SEGMENT_LOW_DAMP_YAML = """
##############################
## Segment Test - Oscillation #
##############################

materials:
  headers: [name, youngs_modulus, density, damping_per_stiffness]
  data:
    - [test_material, 55000000000.0, 724, 0.00077]

points:
  headers: [name, pos_cad, type, wing_idx, transform_idx, extra_mass, body_frame_damping, world_frame_damping, area, drag_coeff]
  data:
    - [anchor, [0.0, 0.0, 0.0], STATIC, nothing, nothing, 0.0, 0.0, 0.0, 0.0, 0.0]
    - [mass_point, [0.0, 0.0, -10.0], DYNAMIC, nothing, nothing, 5.0, 0.0, 0.0, 0.0, 0.0]

segments:
  headers: [name, point_i, point_j, l0, diameter_mm, unit_stiffness, unit_damping, compression_frac]
  data:
    # unit_stiffness=1000 -> k=100 N/m, unit_damping=100 -> c=10 N*s/m -> zeta=0.224
    - [test_segment, anchor, mass_point, 10.0, 5.0, 1000.0, 100.0, 0.1]
"""

# YAML for horizontal segment drag test - two dynamic points, no extra mass
const SEGMENT_HORIZONTAL_DRAG_YAML = """
##############################
## Horizontal Segment Drag ###
##############################

points:
  headers: [name, pos_cad, type, wing_idx, transform_idx, extra_mass, body_frame_damping, world_frame_damping, area, drag_coeff]
  data:
    - [point_left, [-5.0, 0.0, 50.0], DYNAMIC, nothing, nothing, 0.0, 0.0, 0.0, 0.0, 0.0]
    - [point_right, [5.0, 0.0, 50.0], DYNAMIC, nothing, nothing, 0.0, 0.0, 0.0, 0.0, 0.0]

segments:
  headers: [name, point_i, point_j, l0, diameter_mm, unit_stiffness, unit_damping, compression_frac]
  data:
    - [horiz_segment, point_left, point_right, 10.0, 4.0, 100000.0, 100.0, 0.1]
"""

# YAML for vertical segment wind drag test - two dynamic points, no extra mass
const SEGMENT_VERTICAL_WIND_YAML = """
##############################
## Vertical Segment Wind #####
##############################

points:
  headers: [name, pos_cad, type, wing_idx, transform_idx, extra_mass, body_frame_damping, world_frame_damping, area, drag_coeff]
  data:
    - [point_top, [0.0, 0.0, 60.0], DYNAMIC, nothing, nothing, 0.0, 0.0, 0.0, 0.0, 0.0]
    - [point_bottom, [0.0, 0.0, 50.0], DYNAMIC, nothing, nothing, 0.0, 0.0, 0.0, 0.0, 0.0]

segments:
  headers: [name, point_i, point_j, l0, diameter_mm, unit_stiffness, unit_damping, compression_frac]
  data:
    - [vert_segment, point_top, point_bottom, 10.0, 4.0, 100000.0, 100.0, 0.1]
"""

# YAML for high altitude drag test - starts at 5000m where air is thinner
const SEGMENT_HIGH_ALTITUDE_YAML = """
##############################
## High Altitude Segment #####
##############################

points:
  headers: [name, pos_cad, type, wing_idx, transform_idx, extra_mass, body_frame_damping, world_frame_damping, area, drag_coeff]
  data:
    - [point_left, [-5.0, 0.0, 5000.0], DYNAMIC, nothing, nothing, 0.0, 0.0, 0.0, 0.0, 0.0]
    - [point_right, [5.0, 0.0, 5000.0], DYNAMIC, nothing, nothing, 0.0, 0.0, 0.0, 0.0, 0.0]

segments:
  headers: [name, point_i, point_j, l0, diameter_mm, unit_stiffness, unit_damping, compression_frac]
  data:
    - [horiz_segment, point_left, point_right, 10.0, 4.0, 100000.0, 100.0, 0.1]
"""

@testset "Segment Tests" begin
    # Write YAML to temp files
    tmpdir = mktempdir()
    yaml_path = joinpath(tmpdir, "test_segment_geometry.yaml")
    write(yaml_path, SEGMENT_TEST_YAML)

    yaml_low_damp_path = joinpath(tmpdir, "test_segment_low_damp_geometry.yaml")
    write(yaml_low_damp_path, SEGMENT_LOW_DAMP_YAML)

    yaml_horiz_drag_path = joinpath(tmpdir, "test_horiz_drag_geometry.yaml")
    write(yaml_horiz_drag_path, SEGMENT_HORIZONTAL_DRAG_YAML)

    yaml_vert_wind_path = joinpath(tmpdir, "test_vert_wind_geometry.yaml")
    write(yaml_vert_wind_path, SEGMENT_VERTICAL_WIND_YAML)

    yaml_high_alt_path = joinpath(tmpdir, "test_high_alt_geometry.yaml")
    write(yaml_high_alt_path, SEGMENT_HIGH_ALTITUDE_YAML)

    # Create minimal settings file
    settings_yaml = """
system:
    log_file: "data/segment_test"  # filename without extension  [replay only]
                                   #   use / as path delimiter, even on Windows
    g_earth:     9.81

solver:
    solver: "FBDF"
    abs_tol: 0.0001          # absolute tolerance of the DAE solver [m, m/s]
    rel_tol: 0.0001          # relative tolerance of the DAE solver [-]
    relaxation: 0.6        # relaxation factor of inner linear Newton solver, needed for quasi-steady solver

kite:
    model: ""     # 3D model of the kite
    foil_file: "ram_air_kite/ram_air_kite_foil.dat" # filename for the foil shape
    physical_model: "2plate"            # name of the kite model to use (2plate, ram, etc.)
    struc_geometry_path: "particle_structural_geometry.yaml"  # structural YAML
    aero_geometry_path: "aero_geometry.yaml"    # aerodynamic YAML
    mass: 0.0                               # kite mass [kg]
    quasi_static: false                     # whether to use quasi static kite points or not

tether:
    cd_tether: 0.0             # disable segment aero drag for pure spring-damper test
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
    upwind_elevation: 0.0
    wind_vec: [0.0, 0.0, 0.0]
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
    sys = load_sys_struct_from_yaml(yaml_path; system_name="segment_test", set=set)

    # ========================================================================
    # YAML Loading Verification
    # ========================================================================
    @testset "YAML Loading Verification" begin
        # Verify points were loaded correctly
        @test length(sys.points) == 2
        @test haskey(sys.points, :anchor)
        @test haskey(sys.points, :mass_point)

        # Verify point properties
        anchor = sys.points[:anchor]
        @test anchor.type == SymbolicAWEModels.STATIC
        @test anchor.pos_cad == KVec3(0.0, 0.0, 0.0)
        @test anchor.extra_mass == 0.0

        mass_point = sys.points[:mass_point]
        @test mass_point.type == SymbolicAWEModels.DYNAMIC
        @test mass_point.pos_cad == KVec3(0.0, 0.0, -10.0)
        @test mass_point.extra_mass == 1.0
        @test mass_point.area == 0.0
        @test mass_point.drag_coeff == 0.0

        # Verify segment was loaded correctly
        @test length(sys.segments) == 1
        @test haskey(sys.segments, :test_segment)

        segment = sys.segments[:test_segment]
        @test segment.l0 == 10.0
        @test segment.unit_stiffness == 1000.0
        @test segment.unit_damping == 10.0
        @test segment.diameter == 0.005  # 5mm in meters
        @test segment.compression_frac == 0.1

        println("\n  ====== Loaded segment: l0=$(segment.l0)m, unit_stiffness=$(segment.unit_stiffness)N, unit_damping=$(segment.unit_damping)N·s ======\n")
    end

    # ========================================================================
    # Physics Test 1: No gravity, no wind - point stays still
    # ========================================================================
    @testset "No gravity, no wind - stationary" begin
        set.g_earth = 0.0
        set.v_wind = 0.0

        sys = load_sys_struct_from_yaml(yaml_path; system_name="segment_test", set=set)
        sam = SymbolicAWEModel(set, sys)
        test_init!(sam)

        # Record initial position
        initial_z = sam.sys_struct.points[:mass_point].pos_w[3]

        # Run simulation for 1 second
        dt = 0.01
        n_steps = 100
        for _ in 1:n_steps
            next_step!(sam; dt=dt, vsm_interval=0)
        end

        # Position should be unchanged
        final_z = sam.sys_struct.points[:mass_point].pos_w[3]
        @test abs(final_z - initial_z) < 0.001  # Should not move more than 1mm

        println("\n  ====== Position drift: $(round(abs(final_z - initial_z)*1000, digits=3)) mm (limit: 1 mm) ======\n")
    end

    # ========================================================================
    # Physics Test 2: With gravity - oscillation dynamics
    # ========================================================================
    @testset "With gravity - oscillation dynamics" begin
        set.g_earth = 9.81
        set.v_wind = 0.0

        # Use low damping YAML with 5.0 kg point mass (>> segment mass ~0.14 kg)
        # This ensures clean spring-damper dynamics where point mass dominates
        sys = load_sys_struct_from_yaml(yaml_low_damp_path; system_name="segment_test", set=set)

        # Verify the properties were loaded
        @test sys.segments[:test_segment].unit_damping == 100.0
        @test sys.points[:mass_point].extra_mass == 5.0

        sam = SymbolicAWEModel(set, sys)
        test_init!(sam)

        # Physics parameters
        l0 = 10.0  # rest length [m]
        z0 = -10.0  # initial z position (unstretched spring)
        unit_stiffness = 1000.0
        unit_damping = 100.0

        # Calculate expected total mass: extra_mass + half segment mass
        segment = sys.segments[:test_segment]
        half_segment_mass = 0.5 * set.rho_tether * π * (segment.diameter / 2)^2 * l0
        expected_total_mass = sys.points[:mass_point].extra_mass + half_segment_mass

        # Verify total_mass is correctly computed
        m = sam.sys_struct.points[:mass_point].total_mass
        @test m ≈ expected_total_mass rtol=0.01

        # Expected equilibrium stretch, accounting for k = unit_stiffness / len
        # At equilibrium: k * stretch = m * g
        # (unit_stiffness / (l0 + stretch)) * stretch = m * g
        # Solving: stretch = m * g * l0 / (unit_stiffness - m * g)
        stretch_eq = m * set.g_earth * l0 / (unit_stiffness - m * set.g_earth)
        len_eq = l0 + stretch_eq
        z_eq_expected = z0 - stretch_eq

        # At equilibrium length, effective k and c are:
        k = unit_stiffness / len_eq
        c = unit_damping / len_eq

        # Expected natural frequency: omega_n = sqrt(k/m)
        omega_n_expected = sqrt(k / m)

        # Expected damping ratio: zeta = c / (2 * sqrt(k*m))
        zeta_expected = c / (2 * sqrt(k * m))

        # Expected damped frequency: omega_d = omega_n * sqrt(1 - zeta^2)
        omega_d_expected = omega_n_expected * sqrt(1 - zeta_expected^2)

        # Run simulation for several periods
        dt = 0.01  # Small timestep for accuracy
        total_time = 50.0  # Several oscillation periods
        n_steps = Int(ceil(total_time / dt))

        z_history = Float64[]
        t_history = Float64[]

        for i in 1:n_steps
            next_step!(sam; dt=dt, vsm_interval=0)
            push!(z_history, sam.sys_struct.points[:mass_point].pos_w[3])
            push!(t_history, i * dt)
        end

        # Check equilibrium (final position should converge)
        # Calculate stiffness from measured equilibrium and compare to expected
        z_final = z_history[end]
        stretch_measured = z0 - z_final
        len_measured = l0 + stretch_measured
        k_measured = unit_stiffness / len_measured
        @test k_measured ≈ k rtol=0.001

        # Find peaks for frequency and damping estimation
        peaks_t = Float64[]
        peaks_z = Float64[]
        for i in 2:(length(z_history)-1)
            if z_history[i] > z_history[i-1] && z_history[i] > z_history[i+1]
                push!(peaks_t, t_history[i])
                push!(peaks_z, z_history[i])
            end
        end

        @test length(peaks_t) >= 3  # Should have at least 3 peaks

        # Estimate damped frequency from peak spacing
        periods = diff(peaks_t)
        avg_period = mean(periods)
        omega_d_measured = 2π / avg_period
        @test omega_d_measured ≈ omega_d_expected rtol=0.1  # Within 10%

        # Estimate damping ratio from logarithmic decrement
        # δ = ln(A₁/A₂) = 2πζ/√(1-ζ²)
        # Solving for ζ: ζ = δ/√(4π² + δ²)
        # Use actual equilibrium (z_final) for amplitude calculation, not theoretical
        amps = abs.(peaks_z .- z_final)
        @test amps[end] < amps[1]  # Amplitude should decrease

        # Average multiple log decrements for robustness
        log_decrements = [log(amps[i] / amps[i+1]) for i in 1:min(5, length(amps)-1) if amps[i+1] > 0.001]
        log_decrement = mean(log_decrements)
        zeta_measured = log_decrement / sqrt(4π^2 + log_decrement^2)

        println("\n  ====== Damping ratio: measured=$(round(zeta_measured, digits=3)), expected=$(round(zeta_expected, digits=3))")
        println("  ====== Stiffness: measured=$(round(k_measured, digits=1)) N/m, expected=$(round(k, digits=1)) N/m ======\n")
        @test zeta_measured ≈ zeta_expected rtol=0.2  # Within 20%
    end

    # ========================================================================
    # Physics Test 3: Horizontal segment with gravity - drag terminal velocity
    # ========================================================================
    @testset "Horizontal segment - gravity drag terminal velocity" begin
        set.g_earth = 9.81
        set.v_wind = 0.0
        set.cd_tether = 0.958  # enable segment aero drag for this test

        sys = load_sys_struct_from_yaml(yaml_horiz_drag_path; system_name="segment_test", set=set)
        segment = sys.segments[:horiz_segment]
        L = segment.l0  # 10.0 m
        d = segment.diameter  # 0.004 m (4mm)

        # Segment mass from tether density
        rho_tether = set.rho_tether  # 724 kg/m^3
        segment_mass = rho_tether * π * (d/2)^2 * L
        m_total = segment_mass
        cd = set.cd_tether  # 0.958

        # Verify no extra mass on points (segment mass only)
        @test sys.points[:point_left].extra_mass == 0.0
        @test sys.points[:point_right].extra_mass == 0.0

        sam = SymbolicAWEModel(set, sys)
        test_init!(sam)

        # Run simulation until terminal velocity is reached
        dt = 0.1
        n_steps = 1000

        for _ in 1:n_steps
            next_step!(sam; dt=dt, vsm_interval=0)
        end

        # Check segment stays horizontal (both points at same height)
        final_z_left = sam.sys_struct.points[:point_left].pos_w[3]
        final_z_right = sam.sys_struct.points[:point_right].pos_w[3]
        @test abs(final_z_left - final_z_right) < 0.01  # Within 1cm

        # Calculate air density at final height using same model as simulation
        # Note: simulation clamps height to max(0.0, h) for density calculation
        final_height = (final_z_left + final_z_right) / 2
        rho_air = SymbolicAWEModels.calc_rho(sam.am, max(0.0, final_height))

        # Terminal velocity: v_t = sqrt(2 * m * g / (rho * cd * L * d))
        v_terminal_expected = sqrt(2 * m_total * set.g_earth / (rho_air * cd * L * d))

        # Check terminal velocity (downward)
        vz_left = sam.sys_struct.points[:point_left].vel_w[3]
        vz_right = sam.sys_struct.points[:point_right].vel_w[3]
        avg_vz = (vz_left + vz_right) / 2

        @test avg_vz < 0  # Moving downward
        @test abs(avg_vz) ≈ v_terminal_expected rtol=0.001

        println("\n  ====== Terminal velocity: measured=$(round(abs(avg_vz), digits=2)) m/s, expected=$(round(v_terminal_expected, digits=2)) m/s (h=$(round(final_height, digits=0))m) ======\n")

        # Verify drag_force field: each point gets half the
        # segment drag. At terminal velocity total drag = m*g.
        pl = sam.sys_struct.points[:point_left]
        pr = sam.sys_struct.points[:point_right]
        total_drag = pl.drag_force + pr.drag_force
        @test norm(total_drag) ≈ m_total * set.g_earth rtol=0.01
        # Both halves should be roughly equal
        @test norm(pl.drag_force) ≈ norm(pr.drag_force) rtol=0.05
        # Drag should point upward (opposing downward fall)
        @test total_drag[3] > 0
        # Points have zero area so drag is purely from segment
        @test pl.area == 0.0
        @test pl.drag_coeff == 0.0
    end

    # ========================================================================
    # Physics Test 4: Vertical segment with wind - wind drag terminal velocity
    # ========================================================================
    @testset "Vertical segment - wind drag terminal velocity" begin
        set.g_earth = 0.0  # No gravity
        set.v_wind = 10.0  # 10 m/s wind
        set.cd_tether = 0.958  # enable segment aero drag for this test
        set.profile_law = 0  # Use constant wind profile

        # Verify no extra mass on points (segment mass only)
        sys = load_sys_struct_from_yaml(yaml_vert_wind_path; system_name="segment_test", set=set)
        @test sys.points[:point_top].extra_mass == 0.0
        @test sys.points[:point_bottom].extra_mass == 0.0

        sam = SymbolicAWEModel(set, sys)
        test_init!(sam)

        # Record initial positions
        initial_z_top = sam.sys_struct.points[:point_top].pos_w[3]
        initial_z_bottom = sam.sys_struct.points[:point_bottom].pos_w[3]

        # Run simulation until terminal velocity is reached
        dt = 0.1
        n_steps = 1000

        for _ in 1:n_steps
            next_step!(sam; dt=dt, vsm_interval=0)
        end

        # Check segment stays vertical (points at same height difference)
        final_z_top = sam.sys_struct.points[:point_top].pos_w[3]
        final_z_bottom = sam.sys_struct.points[:point_bottom].pos_w[3]
        height_diff_initial = initial_z_top - initial_z_bottom
        height_diff_final = final_z_top - final_z_bottom
        @test abs(height_diff_final - height_diff_initial) < 0.1  # Within 10cm

        # Check vertical position unchanged (no gravity)
        @test abs(final_z_top - initial_z_top) < 0.5  # Within 50cm (some drift acceptable)
        @test abs(final_z_bottom - initial_z_bottom) < 0.5

        # Check terminal velocity matches wind velocity
        # Wind blows in x direction (upwind_dir = -90 means wind from +x)
        # At steady state, segment velocity = wind velocity
        vx_top = sam.sys_struct.points[:point_top].vel_w[1]
        vx_bottom = sam.sys_struct.points[:point_bottom].vel_w[1]
        avg_vx = (vx_top + vx_bottom) / 2

        # Wind direction: upwind_dir = -90 deg means wind blows in x direction
        @test avg_vx > 0  # Moving in +x direction (with wind)
        @test abs(avg_vx) ≈ set.v_wind rtol=0.01  # Within 1% of wind speed

        println("\n  ====== Wind-driven velocity: measured=$(round(avg_vx, digits=2)) m/s, expected=$(set.v_wind) m/s ======\n")
    end

    # ========================================================================
    # Physics Test 5: High altitude segment - lower air density = faster fall
    # Demonstrates that rho is NOT clamped for positive heights
    # ========================================================================
    @testset "High altitude segment - thin air faster fall" begin
        set.g_earth = 9.81
        set.v_wind = 0.0
        set.cd_tether = 0.958
        set.profile_law = 0  # Reset to constant profile

        # Get segment properties
        sys = load_sys_struct_from_yaml(yaml_high_alt_path; system_name="segment_test", set=set)
        segment = sys.segments[:horiz_segment]
        L = segment.l0  # 10.0 m
        d = segment.diameter  # 0.004 m (4mm)

        # Segment mass from tether density
        rho_tether = set.rho_tether
        segment_mass = rho_tether * π * (d/2)^2 * L
        m_total = segment_mass
        cd = set.cd_tether

        sam = SymbolicAWEModel(set, sys)
        test_init!(sam)

        # Run simulation - shorter time since we want to stay at altitude
        dt = 0.1
        n_steps = 200  # 20 seconds

        for _ in 1:n_steps
            next_step!(sam; dt=dt, vsm_interval=0)
        end

        # Get final state
        final_z_left = sam.sys_struct.points[:point_left].pos_w[3]
        final_z_right = sam.sys_struct.points[:point_right].pos_w[3]
        final_height = (final_z_left + final_z_right) / 2

        # Height should still be positive (above sea level)
        @test final_height > 0

        # Calculate air density at altitude (no clamping needed since h > 0)
        rho_at_altitude = SymbolicAWEModels.calc_rho(sam.am, final_height)
        rho_at_sea_level = SymbolicAWEModels.calc_rho(sam.am, 0.0)

        # Air should be thinner at altitude
        @test rho_at_altitude < rho_at_sea_level

        # Terminal velocity at altitude vs sea level
        v_terminal_altitude = sqrt(2 * m_total * set.g_earth / (rho_at_altitude * cd * L * d))
        v_terminal_sea_level = sqrt(2 * m_total * set.g_earth / (rho_at_sea_level * cd * L * d))

        # Should fall faster at altitude due to thinner air
        @test v_terminal_altitude > v_terminal_sea_level

        # Check actual terminal velocity matches expected
        vz_left = sam.sys_struct.points[:point_left].vel_w[3]
        vz_right = sam.sys_struct.points[:point_right].vel_w[3]
        avg_vz = (vz_left + vz_right) / 2

        @test avg_vz < 0  # Moving downward
        @test abs(avg_vz) ≈ v_terminal_altitude rtol=0.01

        println("\n  ====== High altitude (h=$(round(final_height, digits=0))m): v=$(round(abs(avg_vz), digits=2)) m/s")
        println("  ====== Sea level would be: v=$(round(v_terminal_sea_level, digits=2)) m/s ($(round((v_terminal_altitude/v_terminal_sea_level - 1)*100, digits=1))% faster at altitude) ======\n")
    end

    # Cleanup
    rm(tmpdir; recursive=true)
end
nothing
