# Copyright (c) 2025 Bart van de Lint, Jelle Poland
# SPDX-License-Identifier: MPL-2.0

"""
Custom Makie-based polar plotting functions for VortexStepMethod results.
This provides a temporary replacement for PyPlot-based plot_polars during the Makie migration.
"""

using GLMakie
using VortexStepMethod

"""
    generate_polar_data(solver, body_aero, angle_range; 
                        angle_type="angle_of_attack",
                        angle_of_attack=0.0,
                        side_slip=0.0,
                        v_a=10.0)

Generate polar data by sweeping through angle range and solving for each configuration.

# Returns
- `polar_data`: Dictionary with keys :alpha, :CL, :CD, :CS (side force coefficient)
- `rey`: Reynolds number
"""
function generate_polar_data(solver, body_aero, angle_range;
                            angle_type="angle_of_attack",
                            angle_of_attack=0.0,
                            side_slip=0.0,
                            v_a=10.0)
    
    CL_list = Float64[]
    CD_list = Float64[]
    CS_list = Float64[]
    alpha_list = Float64[]
    
    num_angles = length(angle_range)
    for (idx, angle) in enumerate(angle_range)
        # Set angle based on type
        if angle_type == "angle_of_attack"
            alpha = angle
            beta = side_slip
        elseif angle_type == "side_slip"
            alpha = angle_of_attack
            beta = angle
        else
            error("angle_type must be 'angle_of_attack' or 'side_slip'")
        end
        
        # Set flight conditions
        α = deg2rad(alpha)
        β = deg2rad(beta)
        va = v_a * [cos(α)*cos(β), sin(β), sin(α)*cos(β)]
        VortexStepMethod.set_va!(body_aero, va)
        
        @info "Polar sweep $(idx)/$num_angles at $(round(angle, digits=2))°"

        # Solve
        results = VortexStepMethod.solve(solver, body_aero; log=false)
        
        # Extract coefficients
        push!(alpha_list, alpha)
        push!(CL_list, results.CL)
        push!(CD_list, results.CD)
        push!(CS_list, results.CS)
    end
    
    # Calculate Reynolds number (using first wing panel as reference)
    wing = body_aero.bodies[1]
    chord = wing.chord_root  # Use root chord as reference
    rey = v_a * chord / 1.5e-5  # kinematic viscosity of air at sea level
    
    polar_data = Dict(
        :alpha => alpha_list,
        :CL => CL_list,
        :CD => CD_list,
        :CS => CS_list
    )
    
    return polar_data, rey
end

"""
    plot_polars_makie(solver_list, body_aero_list, label_list;
                      angle_range=range(0, 20, 21),
                      angle_type="angle_of_attack",
                      angle_of_attack=0.0,
                      side_slip=0.0,
                      v_a=10.0,
                      title="polar",
                      size=(1200, 800))

Plot polar curves (CL, CD, CS vs angle) using Makie.

# Arguments
- `solver_list`: Vector of VortexStepMethod.Solver objects
- `body_aero_list`: Vector of VortexStepMethod.BodyAerodynamics objects
- `label_list`: Vector of String labels for each case

# Keyword Arguments
- `angle_range`: Range of angles to sweep (default: 0° to 20° in 21 steps)
- `angle_type`: "angle_of_attack" or "side_slip"
- `angle_of_attack`: Fixed angle of attack when sweeping sideslip (default: 0.0°)
- `side_slip`: Fixed sideslip when sweeping angle of attack (default: 0.0°)
- `v_a`: Airspeed in m/s (default: 10.0)
- `title`: Plot title (default: "polar")
- `size`: Figure size in pixels (default: (1200, 800))

# Returns
- `fig`: Makie Figure object
"""
function plot_polars_makie(solver_list, body_aero_list, label_list;
                          angle_range=range(0, 20, 21),
                          angle_type="angle_of_attack",
                          angle_of_attack=0.0,
                          side_slip=0.0,
                          v_a=10.0,
                          title="polar",
                          size=(1200, 800))
    
    # Validate inputs
    if length(solver_list) != length(body_aero_list) || length(solver_list) != length(label_list)
        throw(ArgumentError("Mismatch in number of solvers ($(length(solver_list))), " *
                          "body_aero ($(length(body_aero_list))), and labels ($(length(label_list)))"))
    end
    
    # Generate polar data for all cases
    polar_data_list = []
    labels_with_re = String[]
    
    for (i, (solver, body_aero)) in enumerate(zip(solver_list, body_aero_list))
        polar_data, rey = generate_polar_data(
            solver, body_aero, angle_range;
            angle_type,
            angle_of_attack,
            side_slip,
            v_a
        )
        push!(polar_data_list, polar_data)
        # Update label with Reynolds number
        label_with_re = "$(label_list[i]) Re = $(round(Int64, rey*1e-5))e5"
        push!(labels_with_re, label_with_re)
    end
    
    # Create figure with 3 subplots
    fig = Figure(size=size)
    
    # Determine x-axis label
    x_label = angle_type == "angle_of_attack" ? "Angle of Attack [deg]" : "Sideslip Angle [deg]"
    
    # Plot CL vs angle
    ax1 = Axis(fig[1, 1], xlabel=x_label, ylabel="CL [-]", 
               title="Lift Coefficient")
    
    # Plot CD vs angle
    ax2 = Axis(fig[1, 2], xlabel=x_label, ylabel="CD [-]", 
               title="Drag Coefficient")
    
    # Plot CS vs angle
    ax3 = Axis(fig[2, 1], xlabel=x_label, ylabel="CS [-]", 
               title="Side Force Coefficient")
    
    # Plot CL vs CD (polar)
    ax4 = Axis(fig[2, 2], xlabel="CD [-]", ylabel="CL [-]", 
               title="Lift-Drag Polar")
    
    # Color palette
    colors = Makie.wong_colors()
    
    # Plot data for each case
    for (i, (polar_data, label)) in enumerate(zip(polar_data_list, labels_with_re))
        color = colors[mod1(i, length(colors))]
        
        # CL vs angle
        lines!(ax1, polar_data[:alpha], polar_data[:CL], 
               label=label, color=color, linewidth=2)
        scatter!(ax1, polar_data[:alpha], polar_data[:CL], 
                color=color, markersize=8)
        
        # CD vs angle
        lines!(ax2, polar_data[:alpha], polar_data[:CD], 
               label=label, color=color, linewidth=2)
        scatter!(ax2, polar_data[:alpha], polar_data[:CD], 
                color=color, markersize=8)
        
        # CS vs angle
        lines!(ax3, polar_data[:alpha], polar_data[:CS], 
               label=label, color=color, linewidth=2)
        scatter!(ax3, polar_data[:alpha], polar_data[:CS], 
                color=color, markersize=8)
        
        # CL vs CD polar
        lines!(ax4, polar_data[:CD], polar_data[:CL], 
               label=label, color=color, linewidth=2)
        scatter!(ax4, polar_data[:CD], polar_data[:CL], 
                color=color, markersize=8)
    end
    
    # Add legends
    axislegend(ax1, position=:lt)
    
    # Add overall title
    Label(fig[0, :], title, fontsize=20, font=:bold)
    
    return fig
end

"""
    plot_polars_makie(solver, body_aero, label; kwargs...)

Convenience method for single solver/body_aero case.
"""
function plot_polars_makie(solver, body_aero, label; kwargs...)
    return plot_polars_makie([solver], [body_aero], [label]; kwargs...)
end

# Export the main function
export plot_polars_makie, generate_polar_data
