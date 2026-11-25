# SPDX-FileCopyrightText: 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_backend_parity.jl
#
# The two backends must scatter the same `SystemDefinition` back out of a solved
# model. Every other test reads a handful of named fields, so a quantity the
# monolith writes and the kernel does not stays invisible until someone reads it
# months later and finds a stale value — which is how the kernel came to solve
# its aero on a refitted apparent wind rather than the model's own.
#
# This diffs *every* real-valued field of *every* component, so a new read-back
# on one backend and not the other fails here rather than in a user's log. The
# models are built from the same YAML and only initialised — no settling, no
# state restore — so any difference is the read-back alone.
#
# The comparison is at `init!` and nowhere else on purpose. Stepping first would
# test the integrator, not the read-back: the two backends evaluate the same
# right-hand side through different floating-point association, and a stiff
# mass-spring system amplifies that difference exponentially, so after a few steps
# the two are at genuinely different states and every field differs for a reason
# that is not a bug. Every gap this file has caught so far was visible at `init!`.

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

@isdefined(test_init!) || include(joinpath(@__DIR__, "util.jl"))
@isdefined(write_pressure_fixture) ||
    include(joinpath(@__DIR__, "pressure_fixture.jl"))

using Test
using SymbolicAWEModels
using SymbolicAWEModels: update_run_state!, MonolithBackend, KernelBackend,
                         VortexStepMethod
using KiteUtils

"""
    parity_model(backend, geometry, root)

The 2plate kite on `backend`, built from `geometry` and initialised, with its
run state scattered back out of the solved model. Each backend gets its own copy of
the fixture so neither can read the other's model cache.
"""
function parity_model(backend, geometry, root; aero_mode=nothing, fixture=false)
    data_path = joinpath(root, string(nameof(typeof(backend))), "2plate_kite")
    mkpath(dirname(data_path))
    cp(joinpath(dirname(@__DIR__), "data", "2plate_kite"), data_path; force=true)
    fixture && write_pressure_fixture(data_path)
    set_data_path(data_path)
    set = Settings("system.yaml")
    vsm_set = VortexStepMethod.VSMSettings(
        joinpath(data_path, "vsm_settings.yaml"); data_prefix=false)
    sys, transforms = load_authoring(joinpath(data_path, geometry);
                                     name="backend_parity", set=set)
    aero = isnothing(aero_mode) ? Dict{Symbol, SymbolicAWEModels.AbstractAeroModel}() :
        Dict(wing.name => aero_mode() for wing in sys.wings)
    sam = SymbolicAWEModel(set, sys; transforms, backend, vsm_set, aero)
    init!(sam; prn=false, remake=true)
    update_run_state!(sam.prob, sam.integrator, sam)
    return sam
end

"""Whether `value` is a number or an array of them, and so comparable at all."""
comparable(value) = value isa Real ||
    (value isa AbstractArray && eltype(value) <: Real)

"""
    push_mismatch!(found, tag, idx, name, left, right, atol)

Record `(tag, idx, name, max|Δ|)` in `found` where `left` and `right` are comparable
and disagree by more than `atol`. A value `NaN` on both agrees — that is a quantity
neither backend defines, not a gap between them.
"""
function push_mismatch!(found, tag, idx, name, left, right, atol)
    comparable(left) && comparable(right) || return found
    av, bv = collect(left), collect(right)
    all(isnan, av) && all(isnan, bv) && return found
    delta = maximum(abs.(av .- bv); init=0.0)
    isnan(delta) && (delta = Inf)
    delta > atol && push!(found, (tag, idx, name, delta))
    return found
end

"""
    run_mismatches(kernel, monolith; atol)

Every read-back mismatch between two runs: each quantity of each item of every
collection, and each field of every tube model.
"""
function run_mismatches(kernel, monolith; atol=1e-6)
    found = Tuple{String, Int, Symbol, Float64}[]
    for (collection, quantities) in pairs(SymbolicAWEModels.QUANTITIES)
        for idx in eachindex(getproperty(kernel.definition, collection)),
                name in keys(quantities)
            push_mismatch!(found, String(collection), idx, name,
                           run_value(kernel, collection, idx, name),
                           run_value(monolith, collection, idx, name), atol)
        end
    end
    for idx in eachindex(kernel.definition.tubes)
        a = run_value(kernel, :tubes, idx, :model)
        b = run_value(monolith, :tubes, idx, :model)
        for name in propertynames(a)
            push_mismatch!(found, "tube models", idx, name, getproperty(a, name),
                           getproperty(b, name), atol)
        end
    end
    return found
end

"""One line per mismatch, worst first, for a failing test's message."""
function describe(mismatches)
    sorted = sort(mismatches; by=last, rev=true)
    return join(("$tag[$idx].$name differs by $delta"
                 for (tag, idx, name, delta) in sorted), "\n")
end

@testset "Backend read-back parity" begin
    data_path_before = get_data_path()
    root = mktempdir()
    # Only `AeroPressure` scatters panel loads, exercising the read-back path.
    geometries = (("particle", "particle_structural_geometry.yaml", nothing, false),
                  ("rigid", "rigid_structural_geometry.yaml", nothing, false),
                  ("particle pressure", "particle_structural_geometry.yaml",
                   () -> AeroPressure(), true))

    for (name, geometry, aero_mode, fixture) in geometries
        @testset "$name wing" begin
            kernel = parity_model(KernelBackend(), geometry,
                                  joinpath(root, name); aero_mode, fixture)
            monolith = parity_model(MonolithBackend(), geometry,
                                    joinpath(root, name); aero_mode, fixture)

            mismatches = run_mismatches(kernel, monolith)
            @test isempty(mismatches)
            isempty(mismatches) || println(describe(mismatches))
        end
    end

    set_data_path(data_path_before)
end
