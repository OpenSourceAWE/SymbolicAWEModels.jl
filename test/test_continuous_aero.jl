# Copyright (c) 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_continuous_aero.jl
# What is ContinuousAero's alone, on the particle-dynamics 2plate kite:
# - the BILLOWING spanwise distribution is required, and its absence is reported
# - the quarter-chord strut couple: each panel's load reaches its bounding
#   struts as 0.75·force + couple at the LE station and 0.25·force − couple at
#   the TE station, so the force splits 3:1 and the couple cancels
#
# The contract ContinuousAero shares with AeroPressure is in
# test_continuous_modes.jl, and the one it shares with every aero mode in
# test_aero_modes.jl.

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

@isdefined(test_init!) || include(joinpath(@__DIR__, "util.jl"))

using Test
using SymbolicAWEModels
using SymbolicAWEModels: VortexStepMethod, aero_scatter_entries,
    aero_section_columns, wing_points
using KiteUtils
using LinearAlgebra

"""
    load_continuous_sys(data_path, set; billowing)

`SystemStructure` for the particle 2plate kite carrying [`ContinuousAero`](@ref),
with the VSM spanwise distribution either left at the file's `SPLIT_PROVIDED` or
switched to `BILLOWING`.
"""
function load_continuous_sys(data_path, set; billowing)
    vsm_set = VortexStepMethod.VSMSettings(
        joinpath(data_path, "vsm_settings.yaml"); data_prefix=false)
    if billowing
        for vsm_wing_settings in vsm_set.wings
            vsm_wing_settings.spanwise_panel_distribution =
                VortexStepMethod.BILLOWING
            vsm_wing_settings.billowing_percentage = 8.0
        end
    end
    return load_sys_struct_from_yaml(
        joinpath(data_path, "particle_structural_geometry.yaml");
        system_name="continuous_test", set, vsm_set, aero_mode=ContinuousAero())
end

@testset "ContinuousAero" begin
    tmpdir = mktempdir()
    data_path = joinpath(tmpdir, "2plate_kite")
    cp(joinpath(dirname(@__DIR__), "data", "2plate_kite"), data_path; force=true)
    set_data_path(data_path)
    set = Settings("system.yaml")

    @testset "BILLOWING is required" begin
        @test_throws ErrorException load_continuous_sys(
            data_path, set; billowing=false)
    end

    sys = load_continuous_sys(data_path, set; billowing=true)
    wing = sys.wings[1]
    mode = wing.aero
    n_panels = length(wing.vsm_aero.panels)

    @testset "quarter-chord strut couple" begin
        nodes = wing_points(sys, wing)
        column = aero_section_columns(wing, nodes)
        leading = Set(k for ((_, station), k) in column if station === :LE)
        trailing = Set(k for ((_, station), k) in column if station === :TE)
        entries = aero_scatter_entries(mode, wing, nodes)
        @test !isempty(entries)
        for panel in 1:n_panels
            reached = [e for e in entries if e[1] == panel]
            @test !isempty(reached)
            @test all(e -> e[2] in leading || e[2] in trailing, reached)
            le_force = sum(e[3] for e in reached if e[2] in leading)
            te_force = sum(e[3] for e in reached if e[2] in trailing)
            @test isapprox(le_force, 0.75; atol=1e-10)
            @test isapprox(te_force, 0.25; atol=1e-10)
            @test isapprox(sum(e[4] for e in reached), 0.0; atol=1e-10)
        end
    end

    rm(tmpdir; recursive=true, force=true)
end
nothing
