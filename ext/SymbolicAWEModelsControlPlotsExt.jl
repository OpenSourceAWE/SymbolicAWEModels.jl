# SPDX-FileCopyrightText: 2025 Bart van de Lint, Uwe Fechner
#
# SPDX-License-Identifier: MPL-2.0

module SymbolicAWEModelsControlPlotsExt
using ControlPlots, LaTeXStrings, KiteUtils, SymbolicAWEModels

export plot

"""
    ControlPlots.plot(sys::SystemStructure, lg::SysLog; kwargs...)

Create a multi-panel plot of key simulation results from a `SysLog`.

This function visualizes various aspects of the kite's performance and state,
such as turn rates, reel-out speeds, aerodynamic forces, and wing deformation.
Each panel can be individually enabled or disabled via keyword arguments.

# Arguments
- `sys::SystemStructure`: The system structure, used to get component counts (e.g., number of groups).
- `lg::SysLog`: The simulation log data to be plotted.

# Keyword Arguments
- `plot_all::Bool=true`: Defaults to true, enabling all plot panels. If false, all panels are disabled.
- `plot_turn_rates::Bool=plot_all`: Show the panel with the wing's angular velocities (ω_x, ω_y, ω_z).
- `plot_reelout::Bool=plot_all`: Show the panel with the reel-out velocities of the steering winches.
- `plot_aero::Bool=plot_all`: Show the panel with the z-components of aerodynamic force and moment.
- `plot_twist::Bool=plot_all`: Show the panel with the twist angles for each wing group.
- `plot_aoa::Bool=plot_all`: Show the panel with the angle of attack.
- `plot_heading::Bool=plot_all`: Show the panel with the kite's heading angle.

# Example
```julia
# Plot only the angle of attack and heading
plot(model.sys_struct, log, plot_turn_rates=false, plot_reelout=false, plot_aero=false, plot_twist=false)
```
"""
function ControlPlots.plot(sys::SystemStructure, lg::SysLog;
                           plot_all=true,
                           plot_turn_rates=plot_all,
                           plot_reelout=plot_all,
                           plot_aero=plot_all,
                           plot_twist=plot_all,
                           plot_aoa=plot_all,
                           plot_heading=plot_all,
                           plot_force=plot_all,
                           suffix=" - " * sys.name)
    sl = lg.syslog

    # Initialize empty containers for plot data and labels
    plot_data = []
    plot_labels = []
    plot_ylabels = []

    # Conditionally add data for each plot panel
    if plot_turn_rates
        turn_rates_deg = rad2deg.(hcat(sl.turn_rates...))
        push!(plot_data, [turn_rates_deg[1,:], turn_rates_deg[2,:], turn_rates_deg[3,:]])
        push!(plot_labels, [L"\omega_x"*suffix, L"\omega_y"*suffix, L"\omega_z"*suffix])
        push!(plot_ylabels, "turn rates [°/s]")
    end

    if plot_reelout
        v_reelout_23 = ([sl.v_reelout[i][2] for i in eachindex(sl.v_reelout)], [sl.v_reelout[i][3] for i in eachindex(sl.v_reelout)])
        push!(plot_data, v_reelout_23)
        push!(plot_labels, ["v_ro[2]"*suffix, "v_ro[3]"*suffix])
        push!(plot_ylabels, L"v_{ro}~[m/s]")
    end

    if plot_aero
        aero_force_z = [sl.aero_force_b[i][3] for i in eachindex(sl.aero_force_b)]
        aero_moment_z = [sl.aero_moment_b[i][3] for i in eachindex(sl.aero_moment_b)]
        push!(plot_data, [aero_force_z, aero_moment_z])
        push!(plot_labels, [L"F_{aero,z}"*suffix, L"M_{aero,z}"*suffix])
        push!(plot_ylabels, "aero F/M")
    end

    if plot_twist && !isempty(sys.groups)
        twist_angles_deg = rad2deg.(hcat(sl.twist_angles...))[eachindex(sys.groups),:]
        twist_labels = ["twist[$i]"*suffix for i in eachindex(sys.groups)]
        push!(plot_data, [twist_angles_deg[i,:] for i in eachindex(sys.groups)])
        push!(plot_labels, twist_labels)
        push!(plot_ylabels, "twist [°]")
    end

    if plot_aoa
        AoA_deg = rad2deg.(sl.AoA)
        push!(plot_data, [AoA_deg])
        push!(plot_labels, ["AoA"*suffix])
        push!(plot_ylabels, "AoA [°]")
    end

    if plot_heading
        heading_deg = rad2deg.(sl.heading)
        push!(plot_data, [heading_deg])
        push!(plot_labels, ["heading"*suffix])
        push!(plot_ylabels, "heading [°]")
    end

    if plot_force
        winch_force = [[sl.force[i][j] for i in eachindex(sl.force)] for j in 1:3]
        push!(plot_data, winch_force)
        push!(plot_labels, [L"F_{winch,1}"*suffix, L"F_{winch,2}"*suffix, L"F_{winch,3}"*suffix])
        push!(plot_ylabels, "Winch force [N]")
    end

    # Only create a plot if there is data to show
    if isempty(plot_data)
        @warn "No plot sections enabled. Nothing to display."
        return
    end

    # Calculate a dynamic figure height based on the number of subplots
    num_plots = length(plot_ylabels)
    dynamic_ysize = max(5, num_plots * 2.5) # At least 5, plus 2.5 for each subplot

    # Call the plotx function with the dynamically built arguments
    ControlPlots.plotx(sl.time,
        plot_data...; # Splat the data arrays into separate arguments
        ylabels=plot_ylabels,
        labels=plot_labels,
        fig="Oscillating Steering Input Response")
end

function ControlPlots.plot(sys::SystemStructure, reltime::Real;
                           l_tether=50.0, wing_pos=nothing, zoom=false, front=false, xy=nothing)
    pos = [sys.points[i].pos_w for i in eachindex(sys.points)]
    !isnothing(wing_pos) && (pos = [pos..., wing_pos...])
    seg = [[sys.segments[i].point_idxs[1], sys.segments[i].point_idxs[2]] for i in eachindex(sys.segments)]
    if zoom && !front
        xlim = (pos[end][1] - 6, pos[end][1]+6)
        ylim = (pos[end][3] - 10, pos[end][3]+2)
    elseif zoom && front
        xlim = (pos[end][2] - 6, pos[end][2]+6)
        ylim = (pos[end][3] - 10, pos[end][3]+2)
    elseif !zoom && !front
        xlim = (-0.1l_tether, 1.2l_tether)
        ylim = (-0.1l_tether, 1.2l_tether)
    elseif !zoom && front
        xlim = (-0.75l_tether, 0.75l_tether)
        ylim = (-0.1l_tether, 1.2l_tether)
    end
    ControlPlots.plot2d(pos, seg, reltime; zoom, front, xlim, ylim, dz_zoom=0.6, xy)
end

function ControlPlots.plot(sam::SymbolicAWEModel, reltime::Real; kwargs...)
    wings = sam.sys_struct.wings
    if length(wings) > 0
        wing_pos = [wing.pos_w for wing in wings]
    else
        wing_pos = nothing
    end
    plot(sam.sys_struct, reltime; sam.set.l_tether, wing_pos, kwargs...)
end

end
