using SymbolicAWEModels, VortexStepMethod, ControlPlots, KiteUtils
using YAML

# # Initialize data path for loading settings
# data_path = joinpath(@__DIR__, "..", "..", "data")
# if isdir(data_path)
#     KiteUtils.set_data_path(data_path)
# else
#     @warn "Data directory not found at $data_path"
# end

# --- Analytical solution functions for hanging mass ---
"""
    hanging_mass_equilibrium(l0, mass, g, k)

Analytical equilibrium position for a hanging mass on a spring.
Returns the extension from rest length: Δl = mg/k
"""
function hanging_mass_equilibrium(point_mass, rest_length,line_diameter_mm, set)
    cross_section_area = π * ((1e-3)*line_diameter_mm/2)^2
    k_spring = set.e_tether * cross_section_area  # N/m
    mass_line = set.rho_tether * cross_section_area * rest_length  # kg/m (mass per unit length of the line)
    mass_system = point_mass + mass_line  # Total mass of the system
    println("Hanging Mass Parameters:")
    println("   Point mass: ", point_mass, " kg")
    println("   Rest length: ", rest_length, " m")
    println("   Line diameter: ", line_diameter_mm, " mm")
    println("   E modulus: ", set.e_tether, " Pa")
    println("   Density: ", set.rho_tether, " kg/m^3")
    println("   Gravity: ", set.g_earth, " m/s^2")
    println("   Cross section area: ", round(cross_section_area, digits=6), " m^2")
    println("   Spring constant: ", round(k_spring, digits=2), " N/m")
    println("   Mass system: ", round(mass_system, digits=4), " kg")
    # Calculate extension based on mass and spring constant
    extension = mass_system * set.g_earth / k_spring
    return rest_length + extension
end

println("\n\nHanging Mass Example\n", "="^40)
### Loading Settings
set = SymbolicAWEModels.load_settings("base")  # Loads from the data/ dir using our helper function

# Example usage: settings.v_wind = 10  # Set wind speed to 10 m/s
set.v_wind = 0  # No wind
set.sample_freq = 1  # Increase to 100 Hz for better visualization (dt = 0.01s)
set.abs_tol = 1e-6     # Higher precision for better dynamics resolution
set.rel_tol = 1e-6     # Higher precision for better dynamics resolution
point_mass = 1.0  # Mass of the hanging point in kg
rest_length = 4.0  # Rest length of the segment in meters
line_diameter_mm = 5.0  # Diameter of the segment in mm


# Create two points: anchor point (static) and hanging mass (dynamic)
points = Point[]
push!(points, Point(1, [2.0, 0.0, 5.0], STATIC))          # Anchor point at height 5m
push!(points, Point(2, [2.0, 0.0, 2], DYNAMIC; mass=point_mass)) # Hanging mass at height 2m, 1kg

### Create single segment connecting the points
# l0 is the rest length a bit shorter than the distances between the initial points
# compression_frac is set to 0.001, meaning the spring has 0.1% compresive stiffness compared to elongation stiffness
# diameter_mm is the diameter of the bridle segment in millimeters
# As the same E modulus (e_tether) is used, this determines the stiffness:
#     axial_stiffness = set.e_tether * (diameter_m/2)^2 * π
# and the damping:
#     axial_damping = (set.damping / set.c_spring) * axial_stiffness
# where the set. refers to defined values in settings.yaml
segments = Segment[]
push!(segments, Segment(1, set, (1, 2), BRIDLE; l0=rest_length, compression_frac=0.001, diameter_mm=line_diameter_mm))  # 5mm diameter, 4m rest length

### Transform to position the system
# The base position is set to [2.0, 0.0, 5.0], which is the anchor point position.
# The rot_point_idx is set to 2, which refers to the hanging mass point.
# The orientation from base to rot point is vertical, meaning the Z-axis points downwards.
# To transfer this back to an x-axis aligned with elev=0 and azimuth=0, we need to rotate the system.
# This is done by the Transform constructor: using -90 degrees around the Z-axis and translating it to the anchor point position.

transforms = [Transform(1, -deg2rad(90.0), 0.0, 0.0; base_pos=[2.0, 0.0, 5.0], base_point_idx=1, rot_point_idx=2)]

### Create system structure
# The system structure consists of:
# - name: "hanging_mass"
# - settings: `set`
# - points: `points`
# - segments: `segments`
# - transforms: `transforms`

sys_struct = SymbolicAWEModels.SystemStructure("hanging_mass", set; points, segments, transforms)

### Analyze damping response
include("damping_analysis.jl")
freq, zeta, recommended_damping = analyze_damping_response(sys_struct, set; verbose=true)
println("\n Setting recommended damping to: ", recommended_damping)
set.damping = recommended_damping  # Update settings with recommended damping

### Plot initial state
# even though the tether is not used here, it defines the size of the plot
# and therefore we must set it to a reasonable length
set.l_tether = 5.0
plot(sys_struct, 0.0; zoom=false)

# Create and initialize the symbolic model
sam = SymbolicAWEModel(set, sys_struct)
init!(sam; remake=false)

# Run simulation for longer time with smaller steps for better convergence visualization
for i in 1:30
    current_time = i/set.sample_freq
    plot(sam, current_time; zoom=false)
    next_step!(sam)
end

# --- Final comparison with analytical solution ---
println("\n\nSimulation vs Analytical Comparison\n", "="^45)

### Calculate analytical solution
equilibrium_length = hanging_mass_equilibrium(point_mass, rest_length, line_diameter_mm, set)
equilibrium_z = points[1].pos_w[3] - equilibrium_length  # anchor_z - total_length
final_z = sam.sys_struct.points[2].pos_w[3]
println("Final Results:")
println("  Simulation final z-position: $(round(final_z, digits=4)) m")
println("  Analytical equilibrium z:    $(round(equilibrium_z, digits=4)) m")
println("  Position error (Δz):         $(round(final_z - equilibrium_z, digits=4)) m")

# Calculate relative error
pos_error_percent = abs(final_z - equilibrium_z) / abs(equilibrium_z) * 100

println("  Relative error:              $(round(pos_error_percent, digits=2))%")
