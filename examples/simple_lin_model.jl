# Copyright (c) 2024, 2025 Bart van de Lint, Uwe Fechner
# SPDX-License-Identifier: MIT

using Timers
tic()
@info "Loading packages "

PLOT = false
using Pkg
if ! ("LaTeXStrings" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end
using ControlPlots, LaTeXStrings
using SymbolicAWEModels, KiteUtils, LinearAlgebra, Statistics

toc()

# Initialize model
set = Settings("system_ram.yaml")
set.segments = 3
set_values = [-50, 0.0, 0.0]  # Set values of the torques of the three winches. [Nm]
set.quasi_static = false
set.physical_model = "ram"

@info "Creating wing, aero, vsm_solver, sys_struct and symbolic_awe_model:"
sam = SymbolicAWEModel(set)
sam.set.abs_tol = 1e-3
sam.set.rel_tol = 1e-3
toc()

# Initialize at elevation
set.l_tethers[2] += 0.4
set.l_tethers[3] += 0.4
init_sim!(sam; remake=false, reload=false)
sys = sam.sys

@info "System initialized at:"
toc()

# Stabilize system
find_steady_state!(sam)
simple_linearize(sam)


