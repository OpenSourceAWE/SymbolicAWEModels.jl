# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

# Test that VSM aerodynamic coefficients are independent
# of wind speed magnitude. In incompressible flow, CL, CD,
# and moment coefficients depend only on flow direction
# (alpha, beta) and geometry, not on speed.

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

using Test
using SymbolicAWEModels
using SymbolicAWEModels: VortexStepMethod,
    _vsm_solve_coeffs!, QUATERNION
using KiteUtils
using LinearAlgebra

# ── setup ───────────────────────────────────────────────
pkg_root = dirname(@__DIR__)
src_data = joinpath(pkg_root, "data", "2plate_kite")
tmpdir = mktempdir()
data_path = joinpath(tmpdir, "2plate_kite")
cp(src_data, data_path; force=true)
set_data_path(data_path)

struc_yaml = joinpath(
    data_path, "quat_struc_geometry.yaml")
set = Settings("system.yaml")
vsm_set = VortexStepMethod.VSMSettings(
    joinpath(data_path, "vsm_settings.yaml");
    data_prefix=false)
sys = load_sys_struct_from_yaml(struc_yaml;
    system_name="speed_test", set, vsm_set)
sam = SymbolicAWEModel(set, sys)
init!(sam; remake=false, prn=false)

wing = sam.sys_struct.wings[1]
groups = sam.sys_struct.groups
n_unrefined = wing.vsm_wing.n_unrefined_sections
moment_frac =
    groups[first(wing.group_idxs)].moment_frac

va_b_0 = copy(wing.va_b)
omega_b = copy(wing.ω_b)
theta_0 = zeros(n_unrefined)
for gidx in wing.group_idxs
    g = groups[gidx]
    for ui in g.unrefined_section_idxs
        theta_0[ui] = g.twist
    end
end

# ── tests ───────────────────────────────────────────────
@testset "VSM speed independence" begin
    # Solve at baseline speed
    x_base = _vsm_solve_coeffs!(
        wing, theta_0, va_b_0, omega_b,
        moment_frac, groups)

    # Solve at half speed (same direction)
    va_half = va_b_0 * 0.5
    x_half = _vsm_solve_coeffs!(
        wing, theta_0, va_half, omega_b,
        moment_frac, groups)

    # Solve at double speed (same direction)
    va_double = va_b_0 * 2.0
    x_double = _vsm_solve_coeffs!(
        wing, theta_0, va_double, omega_b,
        moment_frac, groups)

    # All coefficients should be speed-independent
    @testset "half speed" begin
        @test x_half ≈ x_base atol=1e-3
    end
    @testset "double speed" begin
        @test x_double ≈ x_base atol=1e-3
    end
end
nothing
