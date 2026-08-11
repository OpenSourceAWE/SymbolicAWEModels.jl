# SPDX-FileCopyrightText: 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

using PackageCompiler, TOML

root = dirname(@__DIR__)
examples = joinpath(root, "examples")

# SymbolicAWEModels stays out so that Revise can still update it
skip = ("SymbolicAWEModels", "SnoopCompile", "SnoopCompileCore")
deps = keys(TOML.parsefile(joinpath(examples, "Project.toml"))["deps"])
packages = Symbol.(sort(collect(setdiff(deps, skip))))

GC.gc(true)
@info "Free memory: $(round(Sys.free_memory() / 1024^2; digits=1)) MB"
@info "Creating sysimage with $(length(packages)) packages: $(join(packages, ", "))"

PackageCompiler.create_sysimage(
    packages;
    sysimage_path=joinpath(root, "kps-image_tmp.so"),
    project=examples,
    precompile_execution_file=joinpath(root, "test", "test_for_precompile.jl"),
    cpu_target=get(ENV, "SYSIMAGE_CPU_TARGET", "native")
)
