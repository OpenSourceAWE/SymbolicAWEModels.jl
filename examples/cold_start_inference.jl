# SPDX-FileCopyrightText: 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

"""
What a cold `KernelBackend` start spends on inference, and what the right-hand
side costs once it is warm.

A cold start is almost entirely Julia compiling code for this model's kernel
types, and the danger is a type whose inference grows faster than the model does.
One such type cost a large kite 853 s of a 943 s `init!` once GLMakie was loaded, because
loading Makie invalidates the cached result and the kernels were held in a
heterogeneous tuple that `Base.deepcopy_internal(::Tuple, …)` types through a
closure capturing the whole tuple. This reproduces that class of problem on the
14-kernel `AeroPressure` 2-plate kite in a few minutes instead of the large kite's hour.

Reports the inference total, then the dearest inferred instances whose signature
names a kernel, then the right-hand side's bytes and time per call. Naming a kernel
is normal — the evaluation path is built on it. What to look for is a single
instance holding a large share of the total, `Core._typeof_captured_variable` above
all, and a right-hand side that has stopped allocating nothing.

`COLD_START_MAKIE=1` loads GLMakie first, which is what makes the difference
visible; without it the same inference is cached and free. The model is built in
a fresh directory every run, so no serialized bin is ever reused.
"""

using Pkg
if Base.active_project() != joinpath(@__DIR__, "Project.toml")
    Pkg.activate(@__DIR__)
end

using Printf
using SnoopCompileCore

get(ENV, "COLD_START_MAKIE", "0") == "1" && @eval using GLMakie

using KiteUtils: init!
using SymbolicAWEModels, VortexStepMethod

package_root = normpath(joinpath(@__DIR__, ".."))
include(joinpath(package_root, "test", "pressure_fixture.jl"))

"""The 2-plate kite on `AeroPressure`, in a fresh data directory."""
function pressure_model()
    data_path = joinpath(mktempdir(), "2plate_kite")
    cp(joinpath(package_root, "data", "2plate_kite"), data_path; force = true)
    write_pressure_fixture(data_path)
    set_data_path(data_path)
    settings = Settings("system.yaml")
    settings.g_earth = 0.0
    vsm_settings = VortexStepMethod.VSMSettings(
        joinpath(data_path, "vsm_settings.yaml"); data_prefix = false)
    structure = load_sys_struct_from_yaml(
        joinpath(data_path, "particle_structural_geometry.yaml");
        system_name = "cold_start", set = settings, vsm_set = vsm_settings,
        aero_mode = AeroPressure())
    haskey(structure.winches, :main_winch) &&
        (structure.winches[:main_winch].brake = true)
    return SymbolicAWEModel(settings, structure; backend = KernelBackend())
end

"""The signature a frame inferred, in full, or a note that it is unavailable."""
function frame_signature(frame)
    return try
        string(frame.ci.def.specTypes)
    catch
        "unavailable"
    end
end

"""Run `count` right-hand side evaluations and return `(seconds, bytes)`."""
function measure_rhs(rhs, du, u, params, count)
    stats = @timed for _ in 1:count
        rhs(du, u, params, 0.0)
    end
    return stats.time, stats.bytes
end

makie = get(ENV, "COLD_START_MAKIE", "0") == "1"
@printf("cold start inference, makie=%d\n", makie)

model = pressure_model()
init_seconds = @elapsed inference = @snoop_inference init!(model; remake = true,
                                                           prn = false)

using SnoopCompile

frames = SnoopCompile.flatten(inference; sortby = SnoopCompile.exclusive)
total = sum(SnoopCompile.exclusive, frames; init = 0.0)
@printf("init! %.2f s, inference %.2f s over %d instances\n", init_seconds, total,
        length(frames))

captures = [frame for frame in frames
            if occursin("ComponentKernel", frame_signature(frame))]
@printf("instances typing a kernel container: %d, %.3f s\n", length(captures),
        sum(SnoopCompile.exclusive, captures; init = 0.0))
for frame in Iterators.reverse(last(captures, 5))
    @printf("  %8.3f s  %s\n", SnoopCompile.exclusive(frame),
            first(frame_signature(frame), 160))
end

problem = model.prob.prob
rhs = problem.f.f
state = copy(problem.u0)
derivative = similar(state)
rhs(derivative, state, problem.p, 0.0)
for count in (1_000, 10_000)
    seconds, bytes = measure_rhs(rhs, derivative, state, problem.p, count)
    @printf("rhs calls=%6d  %6.2f us/call  %5.1f bytes/call\n", count,
            1e6 * seconds / count, bytes / count)
end
