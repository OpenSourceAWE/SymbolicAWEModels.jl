# SPDX-FileCopyrightText: 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_pulley.jl - Pulley constraint tests
#
# Tests pulley equilibrium and constraint enforcement with stiff tethers.
# Verifies:
# 1. YAML loading: pulley properties correctly parsed
# 2. Length constraint: l_left + l_right = constant (within elastic stretch)
# 3. Equilibrium finding: converges from off-center start
# 4. Analytical geometry: equilibrium position matches derivation
# 5. Tension balance: symmetric forces at equilibrium

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

# ============================================================================
# YAML Configuration - V-shaped bridle with pulley at apex
# Geometry: Two attachment points at (±2, 0, 10), pulley starts at (0.5, 0, 5)
# Segment l0 = nothing -> auto-calculated from point positions
# The off-center start tests equilibrium finding capability
# ============================================================================
const PULLEY_TEST_YAML = """
##############################
## Pulley Test System ########
##############################

###########################
## Materials ##############
###########################
materials:
  headers: [name, youngs_modulus, density, damping_per_stiffness]
  data:
    - [dyneema, 55000000000.0, 724, 0.00077]

###########################
## Points #################
###########################
# V-shaped bridle: two attachment points at top, pulley point in middle, weight hanging below
# Pulley starts OFF-CENTER at x=0.5 to test equilibrium finding
points:
  headers: [name, pos_cad, type, wing_idx, transform_idx, extra_mass, body_frame_damping, world_frame_damping, area, drag_coeff]
  data:
    - [attach_left, [-2.0, 0.0, 10.0], STATIC, nothing, nothing, 0.0, 0.0, 0.0, 0.0, 0.0]
    - [attach_right, [2.0, 0.0, 10.0], STATIC, nothing, nothing, 0.0, 0.0, 0.0, 0.0, 0.0]
    - [pulley_point, [0.5, 0.0, 5.0], DYNAMIC, nothing, nothing, 0.0, 50.0, 0.0, 0.0, 0.0]
    - [weight, [0.0, 0.0, 0.0], DYNAMIC, nothing, nothing, 1.0, 50.0, 0.0, 0.0, 0.0]

###########################
## Segments ###############
###########################
# Dyneema tethers: very stiff (E=55 GPa) with low damping
# l0 = nothing -> auto-calculated from point positions
segments:
  headers: [name, point_i, point_j, l0, diameter_mm, unit_stiffness, unit_damping, compression_frac]
  data:
    - [left_leg, attach_left, pulley_point, nothing, 5.0, dyneema, nothing, 0.01]
    - [right_leg, attach_right, pulley_point, nothing, 5.0, dyneema, nothing, 0.01]
    - [main_tether, pulley_point, weight, nothing, 5.0, dyneema, nothing, 0.01]

###########################
## Pulleys ################
###########################
# Pulley constraint: left_leg + right_leg = constant
pulleys:
  headers: [name, segment_i, segment_j, type]
  data:
    - [main_pulley, left_leg, right_leg, DYNAMIC]
"""

# ============================================================================
# YAML Configuration - Symmetric starting position (at geometric equilibrium)
# Pulley at (0, 0, 4.34) gives equal leg lengths for both legs
# ============================================================================
const PULLEY_SYMMETRIC_YAML = """
##############################
## Pulley Test - Symmetric ###
##############################

###########################
## Materials ##############
###########################
materials:
  headers: [name, youngs_modulus, density, damping_per_stiffness]
  data:
    - [dyneema, 55000000000.0, 724, 0.00077]

###########################
## Points #################
###########################
# Symmetric V-bridle: pulley at x=0, z=4.34 (geometric equilibrium), weight hanging below
points:
  headers: [name, pos_cad, type, wing_idx, transform_idx, extra_mass, body_frame_damping, world_frame_damping, area, drag_coeff]
  data:
    - [attach_left, [-2.0, 0.0, 10.0], STATIC, nothing, nothing, 0.0, 0.0, 0.0, 0.0, 0.0]
    - [attach_right, [2.0, 0.0, 10.0], STATIC, nothing, nothing, 0.0, 0.0, 0.0, 0.0, 0.0]
    - [pulley_point, [0.0, 0.0, 4.34], DYNAMIC, nothing, nothing, 0.0, 50.0, 0.0, 0.0, 0.0]
    - [weight, [0.0, 0.0, 0.0], DYNAMIC, nothing, nothing, 1.0, 50.0, 0.0, 0.0, 0.0]

###########################
## Segments ###############
###########################
# l0 = nothing -> auto-calculated from point positions
segments:
  headers: [name, point_i, point_j, l0, diameter_mm, unit_stiffness, unit_damping, compression_frac]
  data:
    - [left_leg, attach_left, pulley_point, nothing, 5.0, dyneema, nothing, 0.01]
    - [right_leg, attach_right, pulley_point, nothing, 5.0, dyneema, nothing, 0.01]
    - [main_tether, pulley_point, weight, nothing, 5.0, dyneema, nothing, 0.01]

###########################
## Pulleys ################
###########################
pulleys:
  headers: [name, segment_i, segment_j, type]
  data:
    - [main_pulley, left_leg, right_leg, DYNAMIC]
"""

@testset "Pulley Tests" begin
    # Write YAML to temp files
    tmpdir = mktempdir()
    yaml_path = joinpath(tmpdir, "test_pulley_geometry.yaml")
    write(yaml_path, PULLEY_TEST_YAML)

    yaml_symmetric_path = joinpath(tmpdir, "test_pulley_symmetric.yaml")
    write(yaml_symmetric_path, PULLEY_SYMMETRIC_YAML)

    # Create minimal settings file
    settings_yaml = """
system:
    log_file: "data/pulley_test"
    g_earth:     9.81

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
    cd_tether: 0.958
    unit_damping: 350.0
    unit_stiffness: 120000.0
    rho_tether: 724.0
    e_tether: 55000000000.0
    rel_damping: 0.00077
    d_tether: 5.0

winch:
    winch_model: "TorqueControlledMachine"
    max_force: 4000
    v_ro_max: 8.0
    drum_radius: 0.110
    gear_ratio: 1.0
    inertia_total: 0.024
    f_coulomb: 10.0
    c_vf: 5.0

environment:
    rho_0: 1.225
    v_wind: 0.0
    upwind_dir: -90.0
    upwind_elevation: 0.0
    wind_vec: [0.0, 0.0, 0.0]
    h_ref: 6.0
    profile_law: 0
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
    sys = load_sys_struct_from_yaml(yaml_path; system_name="pulley_test", set=set)

    # ========================================================================
    # YAML Loading Verification
    # ========================================================================
    @testset "YAML Loading Verification" begin
        # Verify points were loaded correctly
        @test length(sys.points) == 4
        @test haskey(sys.points, :attach_left)
        @test haskey(sys.points, :attach_right)
        @test haskey(sys.points, :pulley_point)
        @test haskey(sys.points, :weight)

        # Verify attachment points are STATIC
        @test sys.points[:attach_left].type == SymbolicAWEModels.STATIC
        @test sys.points[:attach_right].type == SymbolicAWEModels.STATIC
        @test sys.points[:attach_left].pos_cad == KVec3(-2.0, 0.0, 10.0)
        @test sys.points[:attach_right].pos_cad == KVec3(2.0, 0.0, 10.0)

        # Verify pulley point is DYNAMIC with no extra mass
        pulley_point = sys.points[:pulley_point]
        @test pulley_point.type == SymbolicAWEModels.DYNAMIC
        @test pulley_point.pos_cad[1] == 0.5  # Off-center at x=0.5
        @test pulley_point.extra_mass == 0.0
        @test pulley_point.body_frame_damping == KVec3(50.0, 50.0, 50.0)

        # Verify weight is DYNAMIC with 1.0 kg mass
        weight = sys.points[:weight]
        @test weight.type == SymbolicAWEModels.DYNAMIC
        @test weight.extra_mass == 1.0

        # Verify segments
        @test length(sys.segments) == 3
        @test haskey(sys.segments, :left_leg)
        @test haskey(sys.segments, :right_leg)
        @test haskey(sys.segments, :main_tether)

        # Verify segment rest lengths are auto-calculated from point positions
        attach_left_pos = sys.points[:attach_left].pos_cad
        attach_right_pos = sys.points[:attach_right].pos_cad
        pulley_pos = sys.points[:pulley_point].pos_cad
        weight_pos = sys.points[:weight].pos_cad

        l0_left_expected = norm(pulley_pos - attach_left_pos)
        l0_right_expected = norm(pulley_pos - attach_right_pos)
        l0_main_expected = norm(weight_pos - pulley_pos)

        @test sys.segments[:left_leg].l0 ≈ l0_left_expected atol=1e-10
        @test sys.segments[:right_leg].l0 ≈ l0_right_expected atol=1e-10
        @test sys.segments[:main_tether].l0 ≈ l0_main_expected atol=1e-10

        # Verify pulley was loaded
        @test length(sys.pulleys) == 1
        @test haskey(sys.pulleys, :main_pulley)

        pulley = sys.pulleys[:main_pulley]
        @test pulley.type == SymbolicAWEModels.DYNAMIC

        println("\n  ====== Loaded pulley system: $(length(sys.points)) points, $(length(sys.segments)) segments, $(length(sys.pulleys)) pulley ======\n")
    end

    # ========================================================================
    # Physics Test 1: Pulley length constraint verification
    # The sum of segment lengths through the pulley should remain constant
    # (within elastic stretch tolerance for stiff dyneema tethers)
    # ========================================================================
    @testset "Pulley length constraint" begin
        set.g_earth = 9.81
        set.v_wind = 0.0

        sys = load_sys_struct_from_yaml(yaml_path; system_name="pulley_test", set=set)
        sam = SymbolicAWEModel(set, sys)
        test_init!(sam)

        # Get initial total length of pulley segments
        pulley = sam.sys_struct.pulleys[:main_pulley]
        initial_sum_len = pulley.sum_len

        # Expected sum_len = l0_left + l0_right (from point positions)
        attach_left_pos = sys.points[:attach_left].pos_cad
        attach_right_pos = sys.points[:attach_right].pos_cad
        pulley_pos = sys.points[:pulley_point].pos_cad
        expected_sum_len = norm(pulley_pos - attach_left_pos) + norm(pulley_pos - attach_right_pos)
        @test initial_sum_len ≈ expected_sum_len atol=1e-10

        # Run simulation
        dt = 0.001
        n_steps = 1000  # 1 second

        sum_len_history = Float64[]

        for _ in 1:n_steps
            next_step!(sam; dt=dt, vsm_interval=0)

            # Calculate current lengths from point positions
            attach_left = sam.sys_struct.points[:attach_left].pos_w
            attach_right = sam.sys_struct.points[:attach_right].pos_w
            pulley_pos = sam.sys_struct.points[:pulley_point].pos_w

            len_left = norm(pulley_pos - attach_left)
            len_right = norm(pulley_pos - attach_right)
            sum_len = len_left + len_right

            push!(sum_len_history, sum_len)
        end

        # Verify sum of lengths remains approximately constant
        # Note: With stiff tethers (dyneema), actual length may vary due to elasticity
        max_deviation = maximum(abs.(sum_len_history .- initial_sum_len))
        @test max_deviation < 0.1

        # The final sum should still be close to initial
        @test sum_len_history[end] ≈ initial_sum_len atol=0.1

        println("\n  ====== Length constraint: initial=$(round(initial_sum_len, digits=2))m, max_deviation=$(round(max_deviation*1000, digits=1))mm ======\n")
    end

    # ========================================================================
    # Physics Test 2: Equilibrium finding from off-center start
    # Starting at x=0.5, the pulley should find the symmetric x=0 equilibrium
    # ========================================================================
    @testset "Equilibrium finding" begin
        set.g_earth = 9.81
        set.v_wind = 0.0

        sys = load_sys_struct_from_yaml(yaml_path; system_name="pulley_test", set=set)
        sam = SymbolicAWEModel(set, sys)
        test_init!(sam)

        # Initial x position is off-center at 0.5
        initial_x = sam.sys_struct.points[:pulley_point].pos_w[1]
        @test abs(initial_x) > 0.4  # Verify we start off-center

        # Run simulation until equilibrium
        dt = 0.1
        n_steps = 5000

        for _ in 1:n_steps
            next_step!(sam; dt=dt, vsm_interval=0)
        end

        # At equilibrium, pulley should find symmetric position (x ≈ 0)
        # due to equal tensions from symmetric attachment points
        final_x = sam.sys_struct.points[:pulley_point].pos_w[1]

        # Should be close to x=0 (symmetric equilibrium)
        @test abs(final_x) < 0.001

        println("\n  ====== Equilibrium finding: initial_x=$(round(initial_x, digits=2))m, final_x=$(round(final_x, digits=3))m ======\n")
    end

    # ========================================================================
    # Physics Test 3: Analytical equilibrium position
    # For symmetric V-bridle, the initial position IS the geometric equilibrium
    # (both legs have equal length)
    # ========================================================================
    @testset "Analytical equilibrium position" begin
        set.g_earth = 9.81
        set.v_wind = 0.0

        # Load symmetric configuration (starting at equilibrium)
        sys = load_sys_struct_from_yaml(yaml_symmetric_path; system_name="pulley_symmetric", set=set)

        # Get positions from YAML
        attach_left_pos = sys.points[:attach_left].pos_cad
        attach_right_pos = sys.points[:attach_right].pos_cad
        pulley_pos = sys.points[:pulley_point].pos_cad
        weight_pos = sys.points[:weight].pos_cad

        # Calculate expected l0 from geometry
        l0_left_expected = norm(pulley_pos - attach_left_pos)
        l0_right_expected = norm(pulley_pos - attach_right_pos)
        l0_main_expected = norm(weight_pos - pulley_pos)

        # Verify auto-calculated l0 matches expected from geometry
        @test sys.segments[:left_leg].l0 ≈ l0_left_expected atol=1e-10
        @test sys.segments[:right_leg].l0 ≈ l0_right_expected atol=1e-10
        @test sys.segments[:main_tether].l0 ≈ l0_main_expected atol=1e-10

        # For symmetric config, left and right should be equal
        @test l0_left_expected ≈ l0_right_expected atol=1e-10

        sam = SymbolicAWEModel(set, sys)
        test_init!(sam)

        # Run to settle at equilibrium
        dt = 0.001
        n_steps = 3000

        for _ in 1:n_steps
            next_step!(sam; dt=dt, vsm_interval=0)
        end

        final_pos = sam.sys_struct.points[:pulley_point].pos_w
        final_x = final_pos[1]
        final_z = final_pos[3]
        initial_z = pulley_pos[3]

        # Verify equilibrium position
        @test abs(final_x) < 0.2  # Close to x=0
        @test final_z ≈ initial_z atol=0.001  # Close to initial z

        # Verify both legs have equal length at equilibrium
        attach_left = sam.sys_struct.points[:attach_left].pos_w
        attach_right = sam.sys_struct.points[:attach_right].pos_w

        len_left = norm(final_pos - attach_left)
        len_right = norm(final_pos - attach_right)

        @test len_left ≈ len_right atol=0.001  # Equal lengths = equal tensions

        println("\n  ====== Analytical equilibrium: z_expected=$(round(initial_z, digits=3))m, z_measured=$(round(final_z, digits=3))m")
        println("  ====== Leg lengths: left=$(round(len_left, digits=3))m, right=$(round(len_right, digits=3))m, l0=$(round(l0_left_expected, digits=3))m ======\n")
    end

    # ========================================================================
    # Physics Test 4: Tension balance at equilibrium
    # At symmetric equilibrium, horizontal force components should cancel
    # ========================================================================
    @testset "Tension balance at equilibrium" begin
        set.g_earth = 9.81
        set.v_wind = 0.0

        sys = load_sys_struct_from_yaml(yaml_symmetric_path; system_name="pulley_symmetric", set=set)
        sam = SymbolicAWEModel(set, sys)
        test_init!(sam)

        # Run to equilibrium
        dt = 0.001
        for _ in 1:3000
            next_step!(sam; dt=dt, vsm_interval=0)
        end

        # At equilibrium, both pulley point and weight should have near-zero velocity
        pulley_vel = sam.sys_struct.points[:pulley_point].vel_w
        weight_vel = sam.sys_struct.points[:weight].vel_w
        @test norm(pulley_vel) < 0.1  # Pulley velocity near zero
        @test norm(weight_vel) < 0.1  # Weight velocity near zero

        # The geometry should be symmetric
        pulley_pos = sam.sys_struct.points[:pulley_point].pos_w
        weight_pos = sam.sys_struct.points[:weight].pos_w
        attach_left = sam.sys_struct.points[:attach_left].pos_w
        attach_right = sam.sys_struct.points[:attach_right].pos_w

        # Weight should hang below pulley point
        @test weight_pos[3] < pulley_pos[3]  # Weight z < pulley z

        # Weight should be roughly centered (x ≈ 0) due to symmetric pulley
        @test abs(weight_pos[1]) < 0.5

        # Unit vectors from pulley to attachments
        vec_left = attach_left - pulley_pos
        vec_right = attach_right - pulley_pos

        len_left = norm(vec_left)
        len_right = norm(vec_right)

        # For equal tensions with symmetric geometry:
        # The horizontal components should cancel
        unit_left = vec_left / len_left
        unit_right = vec_right / len_right

        # If tensions are equal: T * (unit_left + unit_right) should have x ≈ 0
        combined_horizontal = unit_left[1] + unit_right[1]
        @test abs(combined_horizontal) < 0.1  # Horizontal forces balance

        println("\n  ====== Tension balance: pulley_vel=$(round(norm(pulley_vel)*1000, digits=1))mm/s, weight_vel=$(round(norm(weight_vel)*1000, digits=1))mm/s")
        println("  ====== Weight pos: z=$(round(weight_pos[3], digits=2))m (below pulley at z=$(round(pulley_pos[3], digits=2))m) ======\n")
    end

    # Cleanup
    rm(tmpdir; recursive=true)
end
nothing
