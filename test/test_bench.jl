# SPDX-FileCopyrightText: 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

# Benchmark tests for ODE RHS of 2plate_kite model.
# Verifies allocation counts for registered functions and
# ensures no allocations originate from package source code.
#
# Usage: jl test/test_bench.jl

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

using Test
using SymbolicAWEModels
using SymbolicAWEModels: VortexStepMethod, SystemStructure, KVec3
using KiteUtils
using BenchmarkTools
using Statistics
using Printf
using Profile

function setup_bench_sam()
    pkg_root = dirname(@__DIR__)
    src_data = joinpath(pkg_root, "data", "2plate_kite")
    tmpdir = mktempdir()
    dp = joinpath(tmpdir, "2plate_kite")
    cp(src_data, dp; force=true)
    set_data_path(dp)
    set = Settings("system.yaml")
    vsm_set = VortexStepMethod.VSMSettings(
        joinpath(dp, "vsm_settings.yaml");
        data_prefix=false)
    struc_yaml = joinpath(dp, "quat_struc_geometry.yaml")
    sys = load_sys_struct_from_yaml(struc_yaml;
        system_name="bench", set, vsm_set)
    sys.winches[:main_winch].brake = true
    sam = SymbolicAWEModel(set, sys)
    init!(sam; remake=false, prn=false)
    return sam
end

function get_pkg_src_files()
    pkg_root = dirname(@__DIR__)
    src_files = Set{String}()
    for (_, _, files) in walkdir(joinpath(pkg_root, "src"))
        for f in files
            endswith(f, ".jl") && push!(src_files, f)
        end
    end
    return src_files
end

@testset verbose = true "Benchmarks" begin
    sam = setup_bench_sam()
    sys = sam.sys_struct
    idx = Int64(1)

    @testset "ODE RHS" begin
        f = sam.integrator.f
        u = copy(sam.integrator.u)
        p = sam.integrator.p
        t = sam.integrator.t
        du = similar(u)
        f(du, u, p, t)

        rhs = @benchmark $f($du, $u, $p, $t) samples = 10
        med_ms = median(rhs.times) / 1e6
        @printf("  median: %8.3f ms | allocs: %d (%d B)\n",
                med_ms, rhs.allocs, rhs.memory)
        @test rhs.allocs == 0
    end

    # --- warmup all accessors once ---
    SymbolicAWEModels.get_l0(sys, idx)
    SymbolicAWEModels.get_extra_mass(sys, idx)
    SymbolicAWEModels.get_pos_w(sys, idx)
    SymbolicAWEModels.get_Q_b_to_w(sys, idx)
    SymbolicAWEModels.get_com_w(sys, idx)
    SymbolicAWEModels.get_R_b_to_p(sys, idx)
    SymbolicAWEModels.get_inertia_principal(sys, idx)
    SymbolicAWEModels.get_vsm_y(sys, idx, 1)
    SymbolicAWEModels.get_vsm_x(sys, idx, 1)
    SymbolicAWEModels.get_vsm_jac(sys, idx, 1, 1)
    SymbolicAWEModels.get_aero_force_override(sys, idx, 1)
    SymbolicAWEModels.get_aero_moment_override(sys, idx, 1)

    # Julia 1.11 has extra allocations in @register_symbolic
    # accessors that 1.12+ optimizes away.
    v11 = VERSION < v"1.12"

    @testset "@register_symbolic" begin
        a = @allocations SymbolicAWEModels.get_l0(sys, idx)
        @test a <= (v11 ? 2 : 0)
        a = @allocations SymbolicAWEModels.get_extra_mass(
            sys, idx)
        @test a <= (v11 ? 2 : 0)
        a = @allocations SymbolicAWEModels.get_pos_w(
            sys, idx)
        @test a <= 1
        a = @allocations SymbolicAWEModels.get_Q_b_to_w(
            sys, idx)
        @test a <= (v11 ? 1 : 0)
        a = @allocations SymbolicAWEModels.get_com_w(
            sys, idx)
        @test a <= (v11 ? 1 : 0)
        a = @allocations SymbolicAWEModels.get_R_b_to_p(
            sys, idx)
        @test a <= (v11 ? 1 : 0)
        a = @allocations SymbolicAWEModels.get_inertia_principal(
            sys, idx)
        @test a <= (v11 ? 1 : 0)
    end

    @testset "VSM accessors" begin
        a = @allocations SymbolicAWEModels.get_vsm_y(
            sys, idx, 1)
        @test a <= (v11 ? 2 : 0)
        a = @allocations SymbolicAWEModels.get_vsm_x(
            sys, idx, 1)
        @test a <= (v11 ? 2 : 0)
        a = @allocations SymbolicAWEModels.get_vsm_jac(
            sys, idx, 1, 1)
        @test a <= 2
        a = @allocations SymbolicAWEModels.get_aero_force_override(
            sys, idx, 1)
        @test a <= 4
        a = @allocations SymbolicAWEModels.get_aero_moment_override(
            sys, idx, 1)
        @test a <= (v11 ? 2 : 0)
    end

    @testset "No package allocations in RHS" begin
        f = sam.integrator.f
        u = copy(sam.integrator.u)
        p = sam.integrator.p
        t = sam.integrator.t
        du = similar(u)
        f(du, u, p, t)
        f(du, u, p, t)
        GC.gc()

        Profile.Allocs.clear()
        Profile.Allocs.@profile sample_rate = 1.0 begin
            f(du, u, p, t)
        end
        results = Profile.Allocs.fetch()
        src_files = get_pkg_src_files()

        skip_bases = ("int.jl", "float.jl", "promotion.jl",
                      "number.jl", "boot.jl")

        pkg_locs = String[]
        for a in results.allocs
            for frame in a.stacktrace
                frame.line <= 0 && continue
                file = string(frame.file)
                any(p -> contains(file, p),
                    ("gc-", "Profile", "datatype.c")) &&
                    continue
                (endswith(file, ".c") ||
                    endswith(file, ".h")) && continue
                file == ":-1" && continue
                fname = basename(file)
                endswith(fname, ".so") && continue
                fname in skip_bases && continue
                if fname in src_files
                    loc = "$fname:$(frame.line)"
                    loc ∉ pkg_locs && push!(pkg_locs, loc)
                end
                break  # only check first meaningful frame
            end
        end

        if !isempty(pkg_locs)
            println("  Allocating locations in package src:")
            for loc in pkg_locs
                println("    ", loc)
            end
        end
        @test isempty(pkg_locs)
    end
end
nothing
