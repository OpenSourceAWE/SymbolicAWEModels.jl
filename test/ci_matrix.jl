# SPDX-FileCopyrightText: 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Prints the `include:` list of the CI test matrix as JSON, one entry per runner
# configuration and test group. Runs on a bare Julia with no project, so it has
# to stay free of dependencies.

include("suite_groups.jl")

"""
    RUNNER_CONFIGS

Runner configurations the suite runs on, each crossed with the test groups of
its backend. `coverage` marks the ones that report to Codecov: instrumented code
runs slower, and `Pkg.test` precompiles the whole tree again under the coverage
configuration because Julia keys its caches by that configuration.
"""
RUNNER_CONFIGS = [
    (os="ubuntu-latest", version="1.11", arch="x64", backend="monolith",
     prefix="xvfb-run -a", coverage=false),
    (os="ubuntu-latest", version="1.12", arch="x64", backend="monolith",
     prefix="xvfb-run -a", coverage=true),
    (os="ubuntu-latest", version="1.12", arch="x64", backend="kernel",
     prefix="xvfb-run -a", coverage=true),
    (os="windows-latest", version="1.12", arch="x64", backend="monolith",
     prefix="", coverage=false),
    (os="macOS-latest", version="1.12", arch="aarch64", backend="monolith",
     prefix="", coverage=false),
]

"""
    json_value(value)

Render `value` as a JSON scalar.
"""
json_value(value::Bool) = value ? "true" : "false"
json_value(value::Integer) = string(value)
function json_value(value::AbstractString)
    return '"' * replace(value, "\\" => "\\\\", "\"" => "\\\"") * '"'
end

"""
    json_object(entry)

Render a named tuple as a JSON object.
"""
function json_object(entry)
    fields = ("$(json_value(String(name))):$(json_value(value))"
              for (name, value) in pairs(entry))
    return "{" * join(fields, ",") * "}"
end

"""
    matrix_entries(configs=RUNNER_CONFIGS)

One matrix entry per configuration and test group. `files` holds the group's
files as `Pkg.test` arguments, `group` and `groups` name the job.
"""
function matrix_entries(configs=RUNNER_CONFIGS)
    entries = []
    for config in configs
        groups = test_file_groups(backend_test_files(config.backend))
        for (index, files) in enumerate(groups)
            names = join((replace(file, ".jl" => "") for file in files), " ")
            push!(entries, (; config..., group=index, groups=length(groups),
                            files=names))
        end
    end
    return entries
end

print("[", join((json_object(entry) for entry in matrix_entries()), ","), "]")
