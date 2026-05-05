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
    _update_quaternion_wing!, QUATERNION
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

"""Drive _update_quaternion_wing! at scaled va, return aero_x copy."""
function coeffs_at(scale)
    wing.va_b .= va_b_0 .* scale
    _update_quaternion_wing!(wing, sam.am, groups)
    return copy(wing.aero_x)
end

# ── tests ───────────────────────────────────────────────
@testset "VSM speed independence" begin
    x_base = coeffs_at(1.0)
    x_half = coeffs_at(0.5)
    x_double = coeffs_at(2.0)

    # All coefficients should be speed-independent
    @testset "half speed" begin
        @test x_half ≈ x_base atol=1e-3
    end
    @testset "double speed" begin
        @test x_double ≈ x_base atol=1e-3
    end
end
nothing
