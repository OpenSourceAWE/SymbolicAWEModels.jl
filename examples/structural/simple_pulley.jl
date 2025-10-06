using SymbolicAWEModels, VortexStepMethod, ControlPlots, KiteUtils
using YAML

# Initialize data path for loading settings
data_path = joinpath(@__DIR__, "..", "..", "data")
if isdir(data_path)
    KiteUtils.set_data_path(data_path)
else
    @warn "Data directory not found at $data_path"
end

include("load_settings.jl")

# --- Analytical solution functions for simple pulley ---
"""
    pulley_equilibrium_position(anchor1, anchor2, point_mass, rest_length_per_segment, line_diameter_mm, set)

Analytical equilibrium position for a mass hanging from a pulley system.
Uses proper pulley physics with straight segments and force balance.
Each side is a straight segment; line is linear-elastic.
"""
function pulley_equilibrium_position(anchor1, anchor2, point_mass, rest_length_per_segment, line_diameter_mm, set)
    # Calculate system parameters
    cross_section_area = π * ((1e-3)*line_diameter_mm/2)^2
    span = sqrt(sum((anchor2 .- anchor1).^2))
    k_spring = set.e_tether * cross_section_area / rest_length_per_segment  # CORRECTED: k = EA/L0
    
    # Handle rope mass (crude approximation: add full rope mass to point load)
    total_rope_length = 2 * rest_length_per_segment
    mass_line = set.rho_tether * cross_section_area * total_rope_length
    mass_system = point_mass + mass_line
    total_weight = mass_system * set.g_earth
    
    # Geometric parameters
    midpoint_x = (anchor1[1] + anchor2[1]) / 2
    midpoint_y = (anchor1[2] + anchor2[2]) / 2
    midpoint_z = (anchor1[3] + anchor2[3]) / 2
    
    println("Pulley System Parameters:")
    println("   Point mass: ", point_mass, " kg")
    println("   Rest length per segment: ", rest_length_per_segment, " m")
    println("   Total rope length: ", total_rope_length, " m")
    println("   Line diameter: ", line_diameter_mm, " mm")
    println("   E modulus: ", set.e_tether, " Pa")
    println("   Density: ", set.rho_tether, " kg/m^3")
    println("   Gravity: ", set.g_earth, " m/s^2")
    println("   Cross section area: ", round(cross_section_area, digits=6), " m^2")
    println("   Spring constant per segment: ", round(k_spring, digits=2), " N/m")
    println("   Mass system: ", round(mass_system, digits=4), " kg")
    println("   Anchor span: ", round(span, digits=2), " m")
    println("   Total weight: ", round(total_weight, digits=2), " N")
    
    # Solve for sag s >= 0 using proper pulley physics
    # Equilibrium equation: 2 * T(s) * sin(θ) = W
    # where T(s) = k * (ℓ(s) - L0), sin(θ) = s/ℓ(s), ℓ(s) = sqrt((L/2)^2 + s^2)
    function equilibrium_equation(s)
        if s < 0; return Inf; end
        segment_length = hypot(span/2, s)  # ℓ(s) = sqrt((L/2)^2 + s^2)
        extension = segment_length - rest_length_per_segment
        tension = k_spring * extension
        sin_theta = s / segment_length
        vertical_force = 2 * tension * sin_theta
        return vertical_force - total_weight
    end
    
    # Bracket the solution
    s_lo = 0.0
    s_hi = max(rest_length_per_segment, span) + total_weight / max(k_spring, 1e-9)
    
    # Ensure we have a sign change
    f_lo = equilibrium_equation(s_lo)
    f_hi = equilibrium_equation(s_hi)
    
    if f_lo > 0
        error("Rest length per segment appears too short for equilibrium (try larger rest_length_per_segment or smaller load).")
    end
    
    while f_hi < 0
        s_hi *= 2
        f_hi = equilibrium_equation(s_hi)
        if s_hi > 1e6 * max(span, rest_length_per_segment)
            error("Could not bracket a root - system may be unstable.")
        end
    end
    
    # Bisection method to solve for sag
    for _ in 1:200
        s_mid = 0.5 * (s_lo + s_hi)
        f_mid = equilibrium_equation(s_mid)
        
        if abs(f_mid) < 1e-10 || (s_hi - s_lo) < 1e-10
            equilibrium_z = midpoint_z - s_mid
            println("   Equilibrium sag: ", round(s_mid, digits=4), " m")
            return [midpoint_x, midpoint_y, equilibrium_z]
        end
        
        if sign(f_mid) == sign(f_lo)
            s_lo, f_lo = s_mid, f_mid
        else
            s_hi, f_hi = s_mid, f_mid
        end
    end
    
    error("Bisection did not converge - check system parameters.")
end

println("\n\nSimple Pulley Example\n", "="^40)
set = load_settings(joinpath(@__DIR__, "..", "data", "base", "settings.yaml"))  # Loads as Dict
set.v_wind = 0
set.l_tether = 5.0 # set l_tether as it affects the plot size
dynamics_type = DYNAMIC
set.sample_freq = 5  # Increase to 100 Hz for better visualization (dt = 0.01s)

# System parameters
point_mass = 2.0  # Mass of the hanging point in kg
rest_length_per_segment = 3.5  # Rest length per segment in meters
line_diameter_mm = 5.0  # Diameter of the segments in mm

# pulley point was placed in the middle, as the side-to-side movement converges very slowly
points = Point[]
push!(points, Point(1, [0.0, 0.0, 5.0], STATIC))
push!(points, Point(2, [5.0, 0.0, 5.0], STATIC))
push!(points, Point(3, [2.5, 0.0, 1], DYNAMIC; mass=point_mass))

segments = Segment[]
push!(segments, Segment(1, set, (3,1), BRIDLE,; l0=rest_length_per_segment, compression_frac=0.01, diameter_mm=line_diameter_mm))  
push!(segments, Segment(2, set, (3,2), BRIDLE,; l0=rest_length_per_segment, compression_frac=0.01, diameter_mm=line_diameter_mm))

pulleys = Pulley[]
push!(pulleys, Pulley(1, (1,2), DYNAMIC))

transforms = [Transform(1, -deg2rad(0.0), 0.0, 0.0; base_pos=[0.0, 0.0, 5.0], base_point_idx=1, rot_point_idx=2)]
sys_struct = SymbolicAWEModels.SystemStructure("pulley", set; points, segments, pulleys, transforms)
plot(sys_struct, 0.0; zoom=false, l_tether=set.l_tether)

# Analyze the system, to find optimal damping
include("damping_analysis.jl")
freq, zeta, recommended_damping = analyze_damping_response(sys_struct, set; verbose=true)
println("\n Setting recommended damping to: ", recommended_damping)
set.damping = recommended_damping  # Update settings with recommended damping

# Create the symbolic model
sam = SymbolicAWEModel(set, sys_struct)
init!(sam; remake=false)

for i in 1:100
    current_time = i/set.sample_freq
    plot(sam, current_time; zoom=false)
    next_step!(sam)
end

# --- Final comparison with analytical solution ---
println("\n\nSimulation vs Analytical Comparison\n", "="^45)

### Calculate analytical solution
anchor1 = [0.0, 0.0, 5.0]  # First anchor position
anchor2 = [5.0, 0.0, 5.0]  # Second anchor position
equilibrium_pos = pulley_equilibrium_position(anchor1, anchor2, point_mass, rest_length_per_segment, line_diameter_mm, set)

final_pos = sam.sys_struct.points[3].pos_w
final_x, final_y, final_z = final_pos[1], final_pos[2], final_pos[3]

println("Final Results:")
println("  Simulation final z-position: $(round(final_z, digits=4)) m")
println("  Analytical equilibrium z:    $(round(equilibrium_pos[3], digits=4)) m")
println("  Position error (Δz):         $(round(final_z - equilibrium_pos[3], digits=4)) m")

# Calculate relative error
z_error_percent = abs(final_z - equilibrium_pos[3]) / abs(equilibrium_pos[3]) * 100
println("  Relative error:              $(round(z_error_percent, digits=2))%")
