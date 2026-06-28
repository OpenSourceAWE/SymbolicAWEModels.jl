# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Install the General registry, cloning directly from GitHub when the Julia
# package server is unavailable (e.g. a transient 404 from pkg.julialang.org).
using Pkg

if any(reg -> reg.name == "General", Pkg.Registry.reachable_registries())
    @info "General registry already installed; skipping"
else
    try
        Pkg.Registry.add("General")
    catch err
        @warn "Pkg server registry install failed; cloning from GitHub" exception = err
        ENV["JULIA_PKG_SERVER"] = ""
        Pkg.Registry.add("General")
    end
end
