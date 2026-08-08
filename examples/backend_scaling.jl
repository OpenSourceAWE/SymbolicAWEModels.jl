# SPDX-FileCopyrightText: 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

"""
How the two backends scale with model size, measured on a system structure
holding several independent 2-plate kites. The kites share nothing — no segment,
pulley or winch crosses from one to the next — so the state count grows in
proportion to their number and nothing else about the problem changes.

Reports, per kite count and backend, the time to build the model from scratch and
the time for one right-hand-side evaluation.

Build times are only honest in a process that has not built this shape before. The
`KernelBackend` splits its cost in two: assembling the kernels, which the bin saves
and which stays a fraction of a second at every size here, and Julia compiling the
right-hand side and the solver for that particular tuple of kernel types, which is
paid once per process and which the bin does not save. Measure one count and one
backend per process — `SCALING_KITES` and `SCALING_BACKENDS` select them — or the
second measurement reads the first one's warm code and reports the wrong number.
"""

using Pkg
if Base.active_project() != joinpath(@__DIR__, "Project.toml")
    Pkg.activate(@__DIR__)
end

using YAML
using Printf
using LinearAlgebra
using KiteUtils: init!
using SymbolicAWEModels, VortexStepMethod

# Override from the shell to push the sweep further, e.g. SCALING_KITES=8,16.
KITE_COUNTS = haskey(ENV, "SCALING_KITES") ?
              parse.(Int, split(ENV["SCALING_KITES"], ",")) : [1, 2, 4]

# `SCALING_BACKENDS=kernel` or `=monolith` measures one of them, which is what a
# build time needs; the default runs both and is only good for the right-hand side.
BACKENDS = let choice = get(ENV, "SCALING_BACKENDS", "both")
    choice == "kernel" ? (KernelBackend(),) :
    choice == "monolith" ? (MonolithBackend(),) :
    (MonolithBackend(), KernelBackend())
end
KITE_SPACING = 60.0
RHS_CALLS = 2000

# `AeroNone` keeps the measurement about what grows with the kite count — points,
# segments, pulleys and winches. With a VSM mode every copy also brings its own
# nonlinear aero solve, which is per wing rather than per state and drowns the
# structural trend. Set this to `AeroDirect()` to put the aero back in.
AERO_MODE = AeroNone()

pkg_root = dirname(@__DIR__)
source_data = joinpath(pkg_root, "data", "2plate_kite")

"""
    defined_names(geometry) -> Set{String}

Every name the geometry defines, apart from materials. Components refer to each
other by name, so these are exactly the tokens a copy has to rename; materials stay
shared between copies and keep theirs.
"""
function defined_names(geometry)
    names = Set{String}()
    for (section, block) in geometry
        section == "materials" && continue
        block isa Dict || continue
        rows = get(block, "data", nothing)
        rows === nothing && continue
        for row in rows
            if row isa Dict
                haskey(row, "name") && push!(names, string(row["name"]))
            elseif row isa AbstractVector && !isempty(row)
                push!(names, string(row[1]))
            end
        end
    end
    return names
end

"""
    rename_value(value, names, suffix)

`value` with every reference to a name in `names` suffixed, descending into lists.
Numbers, `nothing` and tokens like `DYNAMIC` are left as they are.
"""
function rename_value(value, names, suffix)
    if value isa AbstractString
        return string(value) in names ? string(value) * suffix : value
    elseif value isa AbstractVector
        return [rename_value(element, names, suffix) for element in value]
    else
        return value
    end
end

"""
    copy_section(block, names, suffix, offset)

One geometry section rewritten for a copy: names and name references suffixed, and
a transform's `base_pos` moved sideways by `offset` so the copies do not overlap.
"""
function copy_section(block, names, suffix, offset)
    rows = block["data"]
    copied = map(rows) do row
        if row isa Dict
            entry = Dict{Any, Any}(key => rename_value(value, names, suffix)
                                   for (key, value) in row)
            if haskey(entry, "base_pos")
                base = collect(Float64, row["base_pos"])
                entry["base_pos"] = [base[1], base[2] + offset, base[3]]
            end
            return entry
        end
        return [rename_value(cell, names, suffix) for cell in row]
    end
    fresh = Dict{Any, Any}(key => value for (key, value) in block)
    fresh["data"] = copied
    return fresh
end

"""
    replicate_geometry(path, count, spacing) -> Dict

The geometry at `path` with its components repeated `count` times, each copy
suffixed `_1`, `_2`, … and displaced `spacing` metres along y. Materials are
declared once and shared.
"""
function replicate_geometry(path, count, spacing)
    geometry = YAML.load_file(path)
    count == 1 && return geometry
    names = defined_names(geometry)
    replicated = Dict{Any, Any}()
    for (section, block) in geometry
        if section == "materials" || !(block isa Dict) || !haskey(block, "data")
            replicated[section] = block
            continue
        end
        rows = Any[]
        for k in 1:count
            append!(rows, copy_section(block, names, "_$k", (k - 1) * spacing)["data"])
        end
        fresh = Dict{Any, Any}(key => value for (key, value) in block)
        fresh["data"] = rows
        replicated[section] = fresh
    end
    return replicated
end

"""
    build_model(data_path, geometry_path, count, backend) -> SymbolicAWEModel

A model of `count` kites on `backend`, built from scratch so the reported time is a
build and not a cache read.
"""
function build_model(data_path, geometry_path, count, backend)
    set_data_path(data_path)
    set = Settings("system.yaml")
    set.g_earth = 0.0
    vsm_set = VortexStepMethod.VSMSettings(
        joinpath(data_path, "vsm_settings.yaml"); data_prefix = false)
    sys = load_sys_struct_from_yaml(geometry_path;
                                    system_name = "scaling_$(count)kite", set, vsm_set,
                                    aero_mode = AERO_MODE, prn = false)
    model = SymbolicAWEModel(set, sys; backend)
    init!(model; remake = true, prn = false)
    return model
end

"""
    rhs_seconds(model, calls) -> Float64

Seconds for one right-hand-side evaluation, best of five runs of `calls` calls.
"""
function rhs_seconds(model, calls)
    integrator = model.integrator
    u = copy(integrator.u)
    du = similar(u)
    f, p, t = integrator.f, integrator.p, integrator.t
    for _ in 1:20
        f(du, u, p, t)
    end
    best = Inf
    for _ in 1:5
        elapsed = @elapsed for _ in 1:calls
            f(du, u, p, t)
        end
        best = min(best, elapsed / calls)
    end
    return best
end

"""
    measure(count, backend, work_dir) -> NamedTuple

Build time, right-hand-side time and state count for `count` kites on `backend`.
"""
function measure(count, backend, work_dir)
    data_path = joinpath(work_dir, "kites_$(count)")
    cp(source_data, data_path; force = true)
    geometry_path = joinpath(data_path, "particle_structural_geometry.yaml")
    YAML.write_file(geometry_path,
                    replicate_geometry(joinpath(source_data,
                                                "particle_structural_geometry.yaml"),
                                       count, KITE_SPACING))
    build = @elapsed model = build_model(data_path, geometry_path, count, backend)
    structure = model.sys_struct
    @assert length(structure.bodies) == count "expected $count bodies, got " *
                                              "$(length(structure.bodies))"
    return (; count, backend = nameof(typeof(backend)), build,
            rhs = rhs_seconds(model, RHS_CALLS), states = length(model.integrator.u),
            points = length(structure.points))
end

work_dir = mktempdir()
@printf("%-16s %6s %7s %8s %12s %12s\n", "backend", "kites", "points", "states",
        "build s", "rhs µs")
for count in KITE_COUNTS, backend in BACKENDS
    result = measure(count, backend, work_dir)
    @printf("%-16s %6d %7d %8d %12.1f %12.2f\n", result.backend, result.count,
            result.points, result.states, result.build, result.rhs * 1e6)
    flush(stdout)
end
