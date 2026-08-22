# SPDX-FileCopyrightText: 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

"""
    EXCLUDED_FILES

Test files the suite never runs: a precompile workload and an interactive menu.
"""
EXCLUDED_FILES = ["test_for_precompile.jl", "test_menu.jl"]

"""
    KERNEL_SKIP_FILES

Why each test file is left out when `SYMAWE_TEST_BACKEND=kernel`; every other
file in the suite runs. Entries are printed with their reason at startup, so a
gap in the [`KernelBackend`](@ref) has to be named here to be tolerated rather
than going unnoticed.
"""
KERNEL_SKIP_FILES = Dict(
    "test_aqua.jl" => "package-quality checks, backend-independent",
    "test_helpers.jl" => "covers the test helpers themselves, builds no model",
    "test_yaml_variables.jl" => "YAML loader only, builds no model",
    "test_bench.jl" => "benchmarks; every file here already asserts a " *
                       "zero-allocation RHS through test_init!",
    "test_backend_parity.jl" => "builds both backends itself, so the " *
                                "monolith run already covers it",
)

"""
    FILE_DURATIONS

Measured runtime [s] of each test file, taken from a monolith CI run. Only used
to balance the groups, so stale numbers cost balance, never correctness.
"""
FILE_DURATIONS = Dict(
    "test_aero_modes.jl" => 350,
    "test_makie_extension.jl" => 287,
    "test_backend_parity.jl" => 179,
    "test_getter_allocations.jl" => 146,
    "test_bench.jl" => 142,
    "test_flap_aero.jl" => 115,
    "test_analytic_jacobian.jl" => 93,
    "test_beam_replay.jl" => 79,
    "test_joint.jl" => 73,
    "test_continuous_modes.jl" => 72,
    "test_wing.jl" => 68,
    "test_tether_winch.jl" => 66,
    "test_plate_aero.jl" => 62,
    "test_aqua.jl" => 50,
    "test_transform.jl" => 42,
    "test_flap_beam.jl" => 38,
    "test_principal_frame_invariance.jl" => 36,
    "test_deserialization.jl" => 30,
    "test_wing_dynamics.jl" => 24,
    "test_segment_nonlinear.jl" => 23,
    "test_joint_invariance.jl" => 22,
    "test_segment.jl" => 18,
    "test_pulley.jl" => 18,
    "test_static_twist.jl" => 18,
    "test_timoshenko_joint.jl" => 17,
    "test_twist_alignment.jl" => 16,
    "test_profile_law.jl" => 14,
    "test_point.jl" => 13,
    "test_rigid_body.jl" => 10,
    "test_yaml_weighted_ref.jl" => 10,
    "test_match_aero_sections.jl" => 9,
    "test_yaml_bodies.jl" => 8,
    "test_tether_init.jl" => 6,
    "test_pressure_aero.jl" => 4,
    "test_section_alignment.jl" => 2,
    "test_continuous_aero.jl" => 2,
    "test_principal_body_frame.jl" => 1,
    "test_yaml_variables.jl" => 1,
    "test_heading_calculation.jl" => 1,
    "test_quaternion_conversions.jl" => 1,
    "test_tube_laws.jl" => 1,
    "test_weighted_ref_points.jl" => 1,
    "test_multi_section_group.jl" => 1,
    "test_helpers.jl" => 1,
    "test_quaternion_auto_groups.jl" => 1,
)

"""
    UNKNOWN_DURATION

Runtime [s] assumed for a test file that [`FILE_DURATIONS`](@ref) does not list,
high enough that a newly added heavy file does not land on a full group.
"""
UNKNOWN_DURATION = 60

"""
    DEFAULT_GROUP_COUNT

Number of jobs each CI configuration splits its test files over. A runner spends
about 23 minutes precompiling the package before the first test runs, against 36
minutes of tests in total, so finer splits buy less than they cost; five
configurations times three groups also fits the account's 20 concurrent jobs in
a single wave. `SYMAWE_TEST_GROUPS` overrides it.
"""
DEFAULT_GROUP_COUNT = 3

"""
    file_duration(file)

Measured runtime [s] of `file`, or [`UNKNOWN_DURATION`](@ref) if it has none.
"""
file_duration(file) = get(FILE_DURATIONS, file, UNKNOWN_DURATION)

"""
    all_test_files(dir=@__DIR__)

Every test file the suite runs, sorted by name.
"""
function all_test_files(dir=@__DIR__)
    return sort(filter(readdir(dir)) do file
        startswith(file, "test_") && endswith(file, ".jl") && file ∉ EXCLUDED_FILES
    end)
end

"""
    test_backend_name()

Backend the suite runs on, from `SYMAWE_TEST_BACKEND`; empty means the monolith.
"""
function test_backend_name()
    name = lowercase(strip(get(ENV, "SYMAWE_TEST_BACKEND", "")))
    name in ("", "monolith", "kernel") ||
        error("Unknown SYMAWE_TEST_BACKEND=\"$name\"; " *
              "expected \"monolith\" or \"kernel\".")
    return isempty(name) ? "monolith" : name
end

"""
    backend_skip_reasons(backend_name, files=all_test_files())

Files `backend_name` cannot run, mapped to the reason they are skipped.
"""
function backend_skip_reasons(backend_name, files=all_test_files())
    backend_name == "kernel" || return Dict{String, String}()
    unknown = sort(filter(!in(files), collect(keys(KERNEL_SKIP_FILES))))
    isempty(unknown) ||
        error("KERNEL_SKIP_FILES names files the suite does not have: $unknown")
    return KERNEL_SKIP_FILES
end

"""
    backend_test_files(backend_name, files=all_test_files())

The test files `backend_name` runs, sorted by name.
"""
function backend_test_files(backend_name, files=all_test_files())
    return filter(!in(keys(backend_skip_reasons(backend_name, files))), files)
end

"""
    group_count(files)

How many groups `files` is split over. `SYMAWE_TEST_GROUPS` overrides
[`DEFAULT_GROUP_COUNT`](@ref); the value `file` gives one group per file.
"""
function group_count(files)
    setting = lowercase(strip(get(ENV, "SYMAWE_TEST_GROUPS", "")))
    isempty(setting) && return min(DEFAULT_GROUP_COUNT, length(files))
    setting == "file" && return length(files)
    count = tryparse(Int, setting)
    (isnothing(count) || count < 1) &&
        error("SYMAWE_TEST_GROUPS=\"$setting\" is not a positive integer or \"file\".")
    return min(count, length(files))
end

"""
    test_file_groups(files, count=group_count(files))

Split `files` over `count` groups of roughly equal measured runtime, by handing
the longest remaining file to the group that is cheapest so far.
"""
function test_file_groups(files, count=group_count(files))
    groups = [String[] for _ in 1:count]
    loads = zeros(count)
    for file in sort(files; by=file -> -file_duration(file))
        cheapest = argmin(loads)
        push!(groups[cheapest], file)
        loads[cheapest] += file_duration(file)
    end
    return sort.(groups)
end
