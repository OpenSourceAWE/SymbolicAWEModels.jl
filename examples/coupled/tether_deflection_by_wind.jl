"""
Validation example: vertical tether under pure horizontal wind drag.

This script sets up a simple system with multiple tether segments, hung vertically between two points
with gravity set to zero, under pure horizontal wind loading. This isolates the aerodynamic drag effects
without gravitational interference. The simulation compares the final deflected shape to
an analytical solution for pure drag loading.

This approach cleverly works within the existing horizontal wind system by:
1. Hanging the tether vertically (fixed points at different Z heights)
2. Setting gravity to zero (g_earth = 0)
3. Applying horizontal wind for pure drag loading
4. Comparing against analytical drag deflection theory

Requirements:
- SymbolicAWEModels
- ControlPlots
"""

using SymbolicAWEModels, VortexStepMethod, ControlPlots
using YAML
# include("load_settings.jl")

# --- Analytical wind deflection solution functions ---

"""
    wind_deflection_analytical_solution(vertical_span, total_length, wind_speed, line_diameter_mm, set; n_points=1000)

Analytical solution for a vertical tether under pure horizontal wind drag loading.
Uses small deflection theory for a cable under distributed lateral loading.
Returns analytical y, z coordinates and key parameters.

For a vertical cable of length L under uniform lateral loading q per unit length:
- Maximum deflection = q*L⁴/(8*E*I) for beam theory (small deflections)
- For cable theory with tension T: deflection ≈ q*L²/(8*T) where T is cable tension
"""
function wind_deflection_analytical_solution(vertical_span, total_length, wind_speed, line_diameter_mm, set; n_points=1000)
    # Calculate system parameters
    cross_section_area = π * ((1e-3)*line_diameter_mm/2)^2
    mass_line = set.rho_tether * cross_section_area * total_length
    
    # Wind drag parameters (horizontal wind on vertical cable)
    air_density = set.rho_0  # Use air density from settings (same as simulation)
    drag_coefficient = set.cd_tether  # Use same drag coefficient as simulation
    
    # Area calculation: length × diameter (not cross-sectional area)
    # This matches the simulation: area[segment.idx] ~ len[segment.idx] * get_diameter(psys, segment.idx)
    cable_diameter = (1e-3) * line_diameter_mm  # Convert mm to m
    area_per_length = cable_diameter  # area = length × diameter, so area_per_length = diameter
    wind_force_per_length = 0.5 * air_density * wind_speed^2 * area_per_length * drag_coefficient
    
    println("Pure Wind Drag System Parameters:")
    println("   Vertical span: ", vertical_span, " m")
    println("   Total cable length: ", total_length, " m")
    println("   Line diameter: ", line_diameter_mm, " mm")
    println("   Wind speed (horizontal): ", wind_speed, " m/s")
    println("   E modulus: ", set.e_tether, " Pa")
    println("   Density: ", set.rho_tether, " kg/m^3")
    println("   Drag coefficient: ", drag_coefficient, " (from set.cd_tether)")
    println("   Air density: ", air_density, " kg/m³ (from set.rho_0)")
    println("   Gravity: ", set.g_earth, " m/s^2 (should be 0, to isolate drag effects)")
    println("   Cross section area: ", round(cross_section_area, digits=6), " m^2")
    println("   Total cable mass: ", round(mass_line, digits=4), " kg")
    println("   Wind force per length: ", round(wind_force_per_length, digits=4), " N/m")
    
    # For a vertical cable under lateral loading, use cable theory
    # Tension in cable from self-weight (should be minimal with g=0) plus elastic tension
    # For pure drag with no gravity: assume initial tension from pre-tension or elastic effects
    
    # Estimate tension from elastic effects (cable stretched by small amount)
    stretch_ratio = (total_length - vertical_span) / vertical_span  # Small stretch
    elastic_tension = set.e_tether * cross_section_area * stretch_ratio
    
    # Ensure minimum tension for stability
    if elastic_tension < 1.0
        elastic_tension = 1.0  # Minimum 1N tension
        println("   Warning: Using minimum tension of 1 N for stability")
    end
    
    # For small deflections of a cable under lateral load:
    # Deflection y(z) = (q*z/(2*T)) * (L*z - z²) where q is lateral load per length
    # Maximum deflection at z = L/2: y_max = q*L³/(8*T)
    
    # Generate analytical curve (lateral deflection of vertical cable)
    z = range(0.0, vertical_span; length=n_points)  # Vertical position
    y = (wind_force_per_length .* z ./ (2 * elastic_tension)) .* (vertical_span .* z .- z.^2)
    
    max_deflection = maximum(abs.(y))
    
    println("   Estimated elastic tension: ", round(elastic_tension, digits=2), " N")
    println("   Cable stretch ratio: ", round(stretch_ratio * 100, digits=3), " %")
    println("   Maximum wind deflection: ", round(max_deflection, digits=4), " m")
    
    return y, z, elastic_tension, max_deflection
end

println("\n\nVertical Tether: Pure Wind Drag Example\n", "="^50)

# --- Settings
# set = load_settings(joinpath(@__DIR__, "..", "data", "base", "settings.yaml"))  # Loads as Dict
set = se("base")  # Loads from data/base using our helper function
# set = Setttings("base/system.yaml")  # Loads from data/base/settings.yaml
set.abs_tol = 1e-8
set.rel_tol = 1e-8
set.l_tether = 8        # Set tether length for plot size

# CRITICAL: Set gravity to zero to isolate wind drag effects
set.g_earth = 0.0       # Zero gravity - pure wind drag only

# Vertical tether parameters
vertical_span = 8.0            # Vertical distance between anchors [m]
n_segments = 10                # Number of tether segments (n points = n_segments+1)
total_length = 8.2             # Total unstretched length (slightly longer for tension)
compression_frac = 0.01        # Compression fraction for segments
world_frame_damping = 1
line_diameter_mm = 4.0         # Diameter for analytical calculations
wind_speed = 5.0               # Wind speed [m/s]

# Set horizontal wind
set.v_wind = wind_speed
# set.upwind_dir = 0.0         # Wind direction (commented out to use default)

# --- Points (nodes) - Set up HORIZONTALLY first, transform will make them vertical
# Note: First point at vertical_span will become top, last point at 0 will become bottom
points = Point[]
push!(points, Point(1, [vertical_span, 0.0, 0.0], STATIC))                     # Will become top anchor after transform

# Add dynamic points along horizontal line (will become vertical after transform)
for i in 1:n_segments-1
    x = vertical_span - i * vertical_span / n_segments  # Reverse horizontal spacing
    push!(points, Point(i+1, [x, 0.0, 0.0], DYNAMIC; world_frame_damping=world_frame_damping))    # Points 2,3,4,...
end
push!(points, Point(n_segments+1, [0.0, 0.0, 0.0], STATIC))                    # Will become bottom anchor after transform

# --- Segments (springs/dampers)
segments = Segment[]
l0_per_segment = total_length / n_segments  # Rest length per segment

for i in 1:n_segments
    # Connect consecutive points
    push!(segments, Segment(i, set, (i, i+1), POWER_LINE; l0=l0_per_segment, diameter_mm=line_diameter_mm,
        compression_frac=compression_frac))
end

# --- Transforms
# Transform: Rotate horizontal setup to vertical
# Key: -90° elevation rotates horizontal (X-axis) tether to vertical (Z-axis)
# base_pos sets where the base point (point 1, currently at x=vertical_span) should be positioned
# After transform: point 1 (at vertical_span) becomes top, point n_segments+1 (at 0) becomes bottom
transforms = [Transform(1, -deg2rad(90.0), 0.0, 0.0; base_pos=[0.0, 0.0, vertical_span], base_point_idx=1, rot_point_idx=n_segments+1)]

# --- System structure
sys_struct = SymbolicAWEModels.SystemStructure("wind_drag", set; points, segments, transforms)

# Analyze the system, to find optimal damping
include("damping_analysis.jl")
freq, zeta, recommended_damping = analyze_damping_response(
    sys_struct, set; verbose=true, perturbation_dir=[0.0, 1.0, 0.0]  # Perturbation in wind direction
    )
println("\n Setting recommended damping to: ", recommended_damping)
set.damping = recommended_damping  # Update settings with recommended damping

# Plot initial state
plot(sys_struct, 0.0; zoom=false)

# --- Construct symbolic model
sam = SymbolicAWEModel(set, sys_struct)
init!(sam; remake=false)

# --- Simulate until static equilibrium under wind drag
println("Simulating vertical tether under pure horizontal wind drag...")
n_steps = 800  # More steps for wind deflection convergence
for i in 1:n_steps
    next_step!(sam)
    # Plot every 10 steps to show evolution
    if i % 10 == 0
        plot(sam, i/set.sample_freq; zoom=false)
        # println("Step $i/$n_steps")
    end
end

# --- Final comparison with analytical solution ---
println("\n\nSimulation vs Analytical Comparison\n", "="^60)

### Calculate analytical solution
y_analytical, z_analytical, elastic_tension, max_analytical_deflection = wind_deflection_analytical_solution(
    vertical_span, total_length, wind_speed, line_diameter_mm, set)

# --- Final plot with analytical comparison
plot(sam, n_steps/set.sample_freq; zoom=false)

# Add analytical solution information
try
    # Extract simulation points for analysis
    sim_x = [sam.sys_struct.points[i].pos_w[1] for i in 1:length(points)]
    sim_y = [sam.sys_struct.points[i].pos_w[2] for i in 1:length(points)]  # Y is wind direction
    sim_z = [sam.sys_struct.points[i].pos_w[3] for i in 1:length(points)]
    
    println("\nSimulation analysis:")
    println("Simulation points: $(length(sim_x))")
    println("Analytical points: $(length(y_analytical))")
    
    # Calculate deflection statistics
    max_sim_deflection_y = maximum(abs.(sim_y))
    max_sim_drift_x = maximum(abs.(sim_x))  # Should be minimal
    
    println("Maximum Y deflection (simulation): $(round(max_sim_deflection_y, digits=4)) m")
    println("Maximum X drift: $(round(max_sim_drift_x, digits=6)) m (should be minimal)")
    
catch e
    println("Note: Could not complete simulation analysis: $e")
end

println("Simulation complete. The vertical tether should show horizontal deflection from wind drag.")

# --- Extract final positions for validation and comparison ---
println("\nFinal point positions - Simulation vs Analytical (wind drag):")
println("Point |Simulation (y,z)|Analytical (y,z)|Delta (y,z)")
println("------|----------------|----------------|-----------")

for i in 1:length(points)
    pos = sam.sys_struct.points[i].pos_w
    sim_y, sim_z = pos[2], pos[3]  # Focus on y,z (horizontal deflection, vertical position)
    
    # Find corresponding analytical point
    if i == 1
        anal_y, anal_z = y_analytical[end], vertical_span  # Top anchor (first point after transform)
    elseif i == length(points)
        anal_y, anal_z = 0.0, 0.0  # Bottom anchor (last point after transform)
    else
        # Interpolate analytical solution at simulation point z-coordinate
        z_pos = vertical_span - (i-1) * vertical_span / n_segments  # Reverse mapping due to transform
        anal_idx = argmin(abs.(z_analytical .- z_pos))
        anal_y = y_analytical[anal_idx]
        anal_z = z_analytical[anal_idx]
    end
    
    delta_y = sim_y - anal_y
    delta_z = sim_z - anal_z
    
    println("  $i   | ($(rpad(round(sim_y, digits=3), 5)), $(rpad(round(sim_z, digits=3), 5))) | ($(rpad(round(anal_y, digits=3), 5)), $(rpad(round(anal_z, digits=3), 5))) | ($(rpad(round(delta_y, digits=3), 5)), $(rpad(round(delta_z, digits=3), 5)))")
end

# Calculate RMS error for dynamic points
global rms_error_y = 0.0
global rms_error_z = 0.0
global n_dynamic_points = 0

for i in 2:(length(points)-1)  # Only dynamic points
    pos = sam.sys_struct.points[i].pos_w
    sim_y, sim_z = pos[2], pos[3]
    
    # Interpolate analytical solution
    z_pos = vertical_span - (i-1) * vertical_span / n_segments  # Reverse mapping due to transform
    anal_idx = argmin(abs.(z_analytical .- z_pos))
    anal_y = y_analytical[anal_idx]
    anal_z = z_analytical[anal_idx]
    
    global rms_error_y += (sim_y - anal_y)^2
    global rms_error_z += (sim_z - anal_z)^2
    global n_dynamic_points += 1
end

global rms_error_y = sqrt(rms_error_y / n_dynamic_points)
global rms_error_z = sqrt(rms_error_z / n_dynamic_points)

println("\nRMS Error (dynamic points only):")
println("  Y-direction: $(round(rms_error_y, digits=6)) m")
println("  Z-direction: $(round(rms_error_z, digits=6)) m")
println("  Total RMS:   $(round(sqrt(rms_error_y^2 + rms_error_z^2), digits=6)) m")
