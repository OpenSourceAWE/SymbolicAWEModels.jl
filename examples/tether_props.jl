# SPDX-FileCopyrightText: 2025 Bart van de Lint <bart@vandelint.net>
#
# SPDX-License-Identifier: MPL-2.0

using Pkg
if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end
using SymbolicAWEModels, KiteUtils, Printf, ControlPlots, LaTeXStrings

# --- Setup Models ---
set = Settings("system.yaml")
set.sample_freq = 600
set.abs_tol = 1e-6
set.rel_tol = 1e-6
set.segments = 1
one_seg_sam = SymbolicAWEModel(set, "ram")
init!(one_seg_sam)
one_seg_tether_sam = SymbolicAWEModel(set, "tether")
init!(one_seg_tether_sam)

# --- Calculate Properties and Get Step Response Data ---
find_steady_state!(one_seg_sam; t=10.0, dt=3.0)
axial_stiffness, axial_damping, tether_lens, dt = 
    SymbolicAWEModels.calc_spring_props(one_seg_sam, one_seg_tether_sam; F_step=-0.1)

# --- Print Comparison Table ---
segments = one_seg_sam.sys_struct.segments
tethers = one_seg_sam.sys_struct.tethers
segments = [segments[tether.segment_idxs[1]] for tether in tethers]
real_axial_stiffness = [segment.axial_stiffness for segment in segments]
real_axial_damping = [segment.axial_damping for segment in segments]

println("\n--- Tether Spring Properties ---")
@printf "%-8s | %-15s %-15s %-10s | %-15s %-15s %-10s\n" "Tether" "Calc. Stiffness" "Real Stiffness" "Error (%)" "Calc. Damping" "Real Damping" "Error (%)"
println(repeat("-", 100))
for i in 1:4
    # Calculate relative errors in percent
    stiffness_err = 100 * abs(axial_stiffness[i] - real_axial_stiffness[i]) / real_axial_stiffness[i]
    damping_err   = 100 * abs(axial_damping[i] - real_axial_damping[i]) / real_axial_damping[i]
    # Print data rows
    @printf "%-8d | %-15.2f %-15.2f %-10.2f | %-15.2f %-15.2f %-10.2f\n" i axial_stiffness[i] real_axial_stiffness[i] stiffness_err axial_damping[i] real_axial_damping[i] damping_err
end

# --- Plot Step Response ---
steps = size(tether_lens, 2) - 1
display(plotx(
    dt .* collect(0:steps), 
    tether_lens[1,:] .- tether_lens[1,1], 
    tether_lens[2,:] .- tether_lens[2,1], 
    tether_lens[3,:] .- tether_lens[3,1], 
    tether_lens[4,:] .- tether_lens[4,1];
    title="Force step response",
    ylabels=["Left power [m]", "Right power [m]", "Left steering [m]", "Right steering [m]"],
))

# --- Calculate Properties and Get Step Response Data ---
find_steady_state!(one_seg_sam; t=10.0, dt=3.0)
axial_stiffness, axial_damping, tether_lens, dt = 
    SymbolicAWEModels.calc_spring_props(one_seg_sam, one_seg_tether_sam; F_step=-0.1)

# --- Print Comparison Table ---
segments = one_seg_sam.sys_struct.segments
tethers = one_seg_sam.sys_struct.tethers
segments = [segments[tether.segment_idxs[1]] for tether in tethers]
real_axial_stiffness = [segment.axial_stiffness for segment in segments]
real_axial_damping = [segment.axial_damping for segment in segments]

println("\n--- Tether Spring Properties ---")
@printf "%-8s | %-15s %-15s %-10s | %-15s %-15s %-10s\n" "Tether" "Calc. Stiffness" "Real Stiffness" "Error (%)" "Calc. Damping" "Real Damping" "Error (%)"
println(repeat("-", 100))
for i in 1:4
    # Calculate relative errors in percent
    stiffness_err = 100 * abs(axial_stiffness[i] - real_axial_stiffness[i]) / real_axial_stiffness[i]
    damping_err   = 100 * abs(axial_damping[i] - real_axial_damping[i]) / real_axial_damping[i]
    # Print data rows
    @printf "%-8d | %-15.2f %-15.2f %-10.2f | %-15.2f %-15.2f %-10.2f\n" i axial_stiffness[i] real_axial_stiffness[i] stiffness_err axial_damping[i] real_axial_damping[i] damping_err
end

# --- Plot Step Response ---
steps = size(tether_lens, 2) - 1
display(plotx(
    dt .* collect(0:steps), 
    tether_lens[1,:] .- tether_lens[1,1], 
    tether_lens[2,:] .- tether_lens[2,1], 
    tether_lens[3,:] .- tether_lens[3,1], 
    tether_lens[4,:] .- tether_lens[4,1];
    title="Force step response",
    ylabels=["Left power [m]", "Right power [m]", "Left steering [m]", "Right steering [m]"],
))
