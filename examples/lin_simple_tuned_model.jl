# SPDX-FileCopyrightText: 2025 Bart van de Lint <bart@vandelint.net>
#
# SPDX-License-Identifier: MPL-2.0

using SymbolicAWEModels
using Pkg
if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end
using ControlPlots

set = Settings("system.yaml")
sam = SymbolicAWEModel(set, "ram")
init!(sam)

tether_set = Settings("system.yaml")
tether_sam = SymbolicAWEModel(tether_set, "tether")
init!(tether_sam)

simple_set = Settings("system.yaml")
simple_sam = SymbolicAWEModel(simple_set, "simple_ram")
init!(simple_sam)

sim_oscillate!(sam; total_time=1.0)
SymbolicAWEModels.copy_to_simple!(sam, tether_sam, simple_sam)

bias = 0.2
# sl, _ = sim_oscillate!(sam; total_time=5.0, prn=true, bias) # TODO: add first frac ram model
# display(plot(sam.sys_struct, sl))

lin_model = linearize!(simple_sam)
sl, lin_sl = sim_oscillate!(simple_sam; total_time=1.0, prn=true, bias, lin_model)
display(plot(simple_sam.sys_struct, sl; plot_all=false, plot_heading=true, suffix=" - simple"))
display(plot(simple_sam.sys_struct, lin_sl; plot_all=false, plot_heading=true, suffix=" - lin"))
