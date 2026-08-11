# SPDX-FileCopyrightText: 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_analytic_jacobian.jl - the KernelBackend's composed Jacobian
#
# The Jacobian the solver is handed is composed from per-kernel local Jacobians
# and the constant wiring rather than differentiated globally, so it has to be
# checked against a global differentiation of the same right-hand side.
# Verifies:
# 1. It matches forward mode entry for entry
# 2. Every entry it writes is inside the declared sparsity pattern
# 3. It refreshes: a second state gives a different, still correct, matrix
# 4. `analytic_jacobian=false` leaves the solver to differentiate itself

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

@isdefined(test_init!) || include(joinpath(@__DIR__, "util.jl"))

using Test
using SymbolicAWEModels
using SymbolicAWEModels: KernelBackend, ForwardDiff
using KiteUtils
using LinearAlgebra

# ============================================================================
# YAML Configuration - a pulley bridle carrying a mass, which gives the wiring
# every shape the composition has to handle: several producers feeding one
# input, a weighted sum, and a component whose derivatives read an input its
# outputs ignore.
# ============================================================================
JACOBIAN_TEST_YAML = """
##################################
## Analytical Jacobian Test ######
##################################

variables:
  dyneema:
    youngs_modulus: 55.0e9
    damping_per_stiffness: 0.00077
    density: 724.0

points:
  headers: [name, pos_cad, type, wing_idx, transform_idx, extra_mass, body_frame_damping, world_frame_damping, area, drag_coeff]
  data:
    - [attach_left, [-2.0, 0.0, 10.0], STATIC, nothing, nothing, 0.0, 0.0, 0.0, 0.0, 0.0]
    - [attach_right, [2.0, 0.0, 10.0], STATIC, nothing, nothing, 0.0, 0.0, 0.0, 0.0, 0.0]
    - [pulley_point, [0.5, 0.0, 5.0], DYNAMIC, nothing, nothing, 0.0, 5.0, 0.0, 0.2, 0.9]
    - [weight, [0.0, 0.0, 0.0], DYNAMIC, nothing, nothing, 1.0, 5.0, 0.0, 0.3, 0.9]

segments:
  headers: [name, point_i, point_j, l0, diameter_mm, youngs_modulus, damping_per_stiffness, density, compression_frac]
  data:
    - [left_leg, attach_left, pulley_point, nothing, 5.0, dyneema, 0.01]
    - [right_leg, attach_right, pulley_point, nothing, 5.0, dyneema, 0.01]
    - [main_tether, pulley_point, weight, nothing, 5.0, dyneema, 0.01]

pulleys:
  headers: [name, segment_i, segment_j, type]
  data:
    - [main_pulley, left_leg, right_leg, DYNAMIC]
"""

JACOBIAN_SETTINGS_YAML = """
system:
    log_file: "data/jacobian_test"
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
    v_wind: 8.0
    upwind_dir: -90.0
    upwind_elevation: 0.0
    wind_vec: nothing
    h_ref: 6.0
    profile_law: 1
"""

"""
    jacobian_test_model(; analytic_jacobian)

The bridle model on the [`KernelBackend`](@ref), initialised and stepped once.

The step is required, not cosmetic. `init!` sets each segment's rest length from its
point positions, so an unstepped model sits exactly where a segment's force changes
between slack and taut and the derivative does not exist there. Two differentiations
of the same right-hand side are then free to pick different one-sided slopes, and
they disagree by `1/compression_frac` — which looks exactly like a broken Jacobian
and is not one.
"""
function jacobian_test_model(; analytic_jacobian)
    tmpdir = mktempdir()
    write(joinpath(tmpdir, "jacobian_geometry.yaml"), JACOBIAN_TEST_YAML)
    write(joinpath(tmpdir, "settings.yaml"), JACOBIAN_SETTINGS_YAML)
    write(joinpath(tmpdir, "system.yaml"), "system:\n  sim_settings: settings.yaml\n")
    set_data_path(tmpdir)
    set = Settings("system.yaml")
    sys = load_sys_struct_from_yaml(joinpath(tmpdir, "jacobian_geometry.yaml");
                                    system_name="jacobian_test", set=set)
    sam = SymbolicAWEModel(set, sys; backend=KernelBackend())
    init!(sam; prn=false, remake=true, analytic_jacobian)
    next_step!(sam; dt=0.05)
    return sam
end

"""The composed Jacobian and a global forward-mode one, at the model's state."""
function jacobian_pair(sam)
    jacobian = sam.prob.prob.f.jac
    rhs = jacobian.rhs
    integrator = sam.integrator
    u, p, t = copy(integrator.u), integrator.p, integrator.t
    composed = zeros(length(u), length(u))
    jacobian(composed, u, p, t)
    reference = ForwardDiff.jacobian((du, x) -> rhs(du, x, p, t), similar(u), copy(u))
    return composed, reference, rhs
end

@testset "Analytical Jacobian" begin
    data_path_before = get_data_path()
    sam = jacobian_test_model(; analytic_jacobian=true)

    @testset "Matches global forward mode" begin
        composed, reference, _ = jacobian_pair(sam)
        scale = maximum(abs, reference)
        @test scale > 1.0
        @test maximum(abs, composed .- reference) < 1e-8 * scale
        @test count(!iszero, composed) == count(!iszero, reference)
    end

    @testset "Stays inside the declared sparsity pattern" begin
        composed, _, rhs = jacobian_pair(sam)
        pattern = rhs.system.sparsity
        @test all(iszero(composed[i, j]) || pattern[i, j]
                  for i in axes(composed, 1), j in axes(composed, 2))
    end

    @testset "Refreshes with the state" begin
        composed, _, rhs = jacobian_pair(sam)
        integrator = sam.integrator
        moved = copy(integrator.u)
        moved .+= 0.05 .* (1:length(moved))
        again = zeros(length(moved), length(moved))
        sam.prob.prob.f.jac(again, moved, integrator.p, integrator.t)
        @test again != composed
        reference = ForwardDiff.jacobian(
            (du, x) -> rhs(du, x, integrator.p, integrator.t), similar(moved), moved)
        @test maximum(abs, again .- reference) < 1e-8 * maximum(abs, reference)
    end

    @testset "Can be turned off" begin
        plain = jacobian_test_model(; analytic_jacobian=false)
        @test isnothing(plain.prob.prob.f.jac)
        @test plain.integrator.t > 0.0
    end

    set_data_path(data_path_before)
end
