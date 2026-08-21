# SPDX-FileCopyrightText: 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_backend_parity.jl
#
# The two backends must scatter the same `SystemStructure` back out of a solved
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

using Test
using SymbolicAWEModels
using SymbolicAWEModels: update_sys_struct!, MonolithBackend, KernelBackend,
                         VortexStepMethod
using KiteUtils

"""
    parity_model(backend, geometry, root)

The 2plate kite on `backend`, built from `geometry` and initialised, with its
struct scattered back out of the solved model. Each backend gets its own copy of
the fixture so neither can read the other's model cache.
"""
function parity_model(backend, geometry, root)
    data_path = joinpath(root, string(nameof(typeof(backend))), "2plate_kite")
    mkpath(dirname(data_path))
    cp(joinpath(dirname(@__DIR__), "data", "2plate_kite"), data_path; force=true)
    set_data_path(data_path)
    set = Settings("system.yaml")
    vsm_set = VortexStepMethod.VSMSettings(
        joinpath(data_path, "vsm_settings.yaml"); data_prefix=false)
    sys = load_sys_struct_from_yaml(joinpath(data_path, geometry);
                                    system_name="backend_parity", set=set,
                                    vsm_set=vsm_set)
    sam = SymbolicAWEModel(set, sys; backend)
    init!(sam; prn=false, remake=true)
    update_sys_struct!(sam.prob, sam.integrator, sam.sys_struct)
    return sam
end

"""Whether `value` is a number or an array of them, and so comparable at all."""
comparable(value) = value isa Real ||
    (value isa AbstractArray && eltype(value) <: Real)

"""
    field_mismatches(tag, kernel_components, monolith_components, atol)

Every `(tag, index, field, max|Δ|)` where the two backends' components disagree by
more than `atol`. A field left `NaN` on both agrees — that is a quantity neither
backend defines, not a gap between them.
"""
function field_mismatches(tag, kernel_components, monolith_components, atol)
    found = Tuple{String, Int, Symbol, Float64}[]
    for (idx, (a, b)) in enumerate(zip(kernel_components, monolith_components))
        for name in fieldnames(typeof(a))
            left, right = getfield(a, name), getfield(b, name)
            comparable(left) && comparable(right) || continue
            av, bv = collect(left), collect(right)
            all(isnan, av) && all(isnan, bv) && continue
            delta = maximum(abs.(av .- bv); init=0.0)
            isnan(delta) && (delta = Inf)
            delta > atol && push!(found, (tag, idx, name, delta))
        end
    end
    return found
end

"""
    struct_mismatches(kernel, monolith; atol)

Every read-back mismatch between two structs, over all component groups. `wings` is
left out because it is a filtered view of `bodies` — the same `Body` objects, so
including it would only report each wing twice.
"""
function struct_mismatches(kernel, monolith; atol=1e-6)
    groups = (:points, :twist_surfaces, :segments, :pulleys, :tethers, :winches,
              :bodies, :elastic_joints, :timoshenko_joints)
    return reduce(vcat, (field_mismatches(String(group), getfield(kernel, group),
                                          getfield(monolith, group), atol)
                         for group in groups))
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
    geometries = ("particle" => "particle_structural_geometry.yaml",
                  "rigid" => "rigid_structural_geometry.yaml")

    for (name, geometry) in geometries
        @testset "$name wing" begin
            kernel = parity_model(KernelBackend(), geometry,
                                  joinpath(root, name))
            monolith = parity_model(MonolithBackend(), geometry,
                                    joinpath(root, name))

            mismatches = struct_mismatches(kernel.sys_struct,
                                           monolith.sys_struct)
            @test isempty(mismatches)
            isempty(mismatches) || println(describe(mismatches))
        end
    end

    set_data_path(data_path_before)
end
