# SPDX-FileCopyrightText: 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_per_point_wind.jl - Wind set per point (PerPointWind)
#
# Verifies:
# 1. A uniform per-point wind reproduces the height profile of a constant law
# 2. A segment's tether drag flies in the mean of its two endpoints' winds
# 3. A wind written between steps reaches that step's own point

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

# Two free points on a stiff vertical segment: both carry point drag (area,
# drag_coeff) and the segment carries tether drag, so every wind consumer is live.
PER_POINT_WIND_YAML = """
##############################
## Per-Point Wind Test #######
##############################

points:
  headers: [name, pos_cad, type, wing_idx, transform_idx, extra_mass, body_frame_damping, world_frame_damping, area, drag_coeff]
  data:
    - [point_top, [0.0, 0.0, 60.0], DYNAMIC, nothing, nothing, 0.5, 0.0, 0.0, 0.02, 1.0]
    - [point_bottom, [0.0, 0.0, 50.0], DYNAMIC, nothing, nothing, 0.5, 0.0, 0.0, 0.02, 1.0]

segments:
  headers: [name, point_i, point_j, l0, diameter_mm, unit_stiffness, unit_damping, compression_frac]
  data:
    - [vert_segment, point_top, point_bottom, 10.0, 4.0, 100000.0, 100.0, 0.1]
"""

PER_POINT_WIND_SETTINGS = """
system:
    log_file: "data/per_point_wind_test"
    g_earth:     0.0

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
    v_wind: 10.0
    upwind_dir: -90.0
    upwind_elevation: 0.0
    wind_vec: [10.0, 0.0, 0.0]
    use_wind_vec: true
    profile_law: 0
"""

"""
Step a fresh model on `wind_mode` and return both points' position and velocity.
`winds` gives the per-point wind written before every step, and `point_drag=false`
zeroes the points' own drag area so only the segment's tether drag acts.
"""
function stepped_state(set, yaml_path, wind_mode; steps = 40, dt = 0.05,
                       winds = nothing, point_drag = true)
    sys = load_sys_struct_from_yaml(yaml_path;
        system_name = "per_point_wind_test", set, wind_mode)
    sam = SymbolicAWEModel(set, sys)
    test_init!(sam; prn=false)
    point_drag || foreach(point -> (point.area = 0.0), sam.sys_struct.points)
    for _ in 1:steps
        isnothing(winds) || for (point, wind) in zip(sam.sys_struct.points, winds)
            point.wind_vec .= wind
        end
        next_step!(sam; dt, vsm_interval=0)
    end
    return [copy(point.pos_w) for point in sam.sys_struct.points],
           [copy(point.vel_w) for point in sam.sys_struct.points]
end

@testset "Per-point wind tests" begin
    tmpdir = mktempdir()
    yaml_path = joinpath(tmpdir, "per_point_wind_geometry.yaml")
    write(yaml_path, PER_POINT_WIND_YAML)
    write(joinpath(tmpdir, "settings.yaml"), PER_POINT_WIND_SETTINGS)
    write(joinpath(tmpdir, "system.yaml"), "system:\n  sim_settings: settings.yaml\n")

    data_path_before = get_data_path()
    set_data_path(tmpdir)
    set = Settings("system.yaml")
    uniform = [KVec3(set.wind_vec), KVec3(set.wind_vec)]

    # ========================================================================
    # A uniform per-point wind is the constant profile law
    # ========================================================================
    @testset "Uniform per-point wind matches the constant profile" begin
        profile_pos, profile_vel = stepped_state(set, yaml_path, ProfileWind())
        per_point_pos, per_point_vel =
            stepped_state(set, yaml_path, PerPointWind(); winds = uniform)

        for k in eachindex(profile_pos)
            @test per_point_pos[k] ≈ profile_pos[k] rtol=1e-6
            @test per_point_vel[k] ≈ profile_vel[k] rtol=1e-6
        end
        # The wind actually drove the points downwind, so the parity is not a
        # comparison of two systems that both did nothing.
        @test profile_vel[1][1] > 1.0
    end

    # ========================================================================
    # A segment's drag flies in the mean of its endpoints' winds
    # ========================================================================
    @testset "Tether drag takes the mean of the two endpoints" begin
        split = [KVec3(20.0, 0.0, 0.0), KVec3(0.0, 0.0, 0.0)]
        mean_pos, mean_vel = stepped_state(set, yaml_path, PerPointWind();
            winds = uniform, point_drag = false)
        split_pos, split_vel = stepped_state(set, yaml_path, PerPointWind();
            winds = split, point_drag = false)

        for k in eachindex(mean_pos)
            @test split_pos[k] ≈ mean_pos[k] rtol=1e-6
            @test split_vel[k] ≈ mean_vel[k] rtol=1e-6
        end
        @test mean_vel[1][1] > 0.1
        # Only the segment pushes, and it pushes both endpoints alike.
        @test split_vel[1][1] ≈ split_vel[2][1] rtol=1e-6
    end

    # ========================================================================
    # A wind written per point acts at that point
    # ========================================================================
    @testset "Per-point wind drives its own point" begin
        windward = [KVec3(20.0, 0.0, 0.0), KVec3(0.0, 0.0, 0.0)]
        _, vel = stepped_state(set, yaml_path, PerPointWind();
                               steps = 20, winds = windward)
        @test vel[1][1] > vel[2][1]          # the windward point leads
        @test vel[2][1] > 0.0                # dragged along by the segment
        @test vel[1][1] < windward[1][1]     # but not faster than its own wind

        _, reversed = stepped_state(set, yaml_path, PerPointWind();
                                    steps = 20, winds = reverse(windward))
        @test reversed[2][1] > reversed[1][1]
    end

    set_data_path(data_path_before)
end
