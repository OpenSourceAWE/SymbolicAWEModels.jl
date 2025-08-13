# Copyright (c) 2024, 2025 Bart van de Lint, Uwe Fechner
# SPDX-License-Identifier: MPL-2.0

using Timers
tic()
@info "Loading packages "
using SymbolicAWEModels

PLOT = false
using Pkg
if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end
using ControlPlots, LaTeXStrings
toc()

# Simulation parameters
dt = 0.05
total_time = 10.0  # Longer simulation to see oscillations
vsm_interval = 3

# Steering parameters
steering_freq = 1/2  # Hz - full left-right cycle frequency
steering_magnitude = 10.0      # Magnitude of steering input [Nm]

# Initialize model
set = Settings("system.yaml")
sam = SymbolicAWEModel(set)
SymbolicAWEModels.init!(sam)

find_steady_state!(sam)
bias = 0.15
if set.physical_model == "4_attach_ram"
    bias = 0.05
end
sl, _ = sim_oscillate!(sam; dt, total_time, vsm_interval, steering_freq, steering_magnitude, 
                         bias, prn=true)
display(plot(sam.sys_struct, sl))
