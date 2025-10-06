"""
Validation example: tether fixed at both ends, under gravity,
resulting in a catenary line.

This script sets up a simple system with multiple tether segments, fixed at two points
at equal height, under gravity. The simulation compares the final shape to
the analytical catenary solution.

Requirements:
- SymbolicAWEModels
- ControlPlots
"""

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

# --- Analytical catenary solution functions ---

# ---- root finder (bisection) ----
function _bisect(f, lo, hi; tol=1e-12, maxiter=10_000)
    flo, fhi = f(lo), f(hi)
    if isnan(flo) || isnan(fhi)
        error("Function returned NaN in the bracket.")
    end
    if flo == 0.0; return lo; end
    if fhi == 0.0; return hi; end
    if sign(flo) == sign(fhi)
        error("Bisection requires f(lo) and f(hi) with opposite signs. Got f(lo)=$(flo), f(hi)=$(fhi).")
    end
    for _ in 1:maxiter
        mid = 0.5*(lo + hi)
        fmid = f(mid)
        if abs(fmid) < tol || 0.5*(hi - lo) < tol
            return mid
        end
        if sign(fmid) == sign(flo)
            lo, flo = mid, fmid
        else
            hi, fhi = mid, fmid
        end
    end
    error("Bisection did not converge in $maxiter iterations.")
end

# ---- solve a from total length: t_l = 2a*sinh(L/(2a)) ----
function _solve_a_from_length(L, t_l)
    if t_l <= L
        error("Total length must be greater than span: got t_l=$t_l, L=$L.")
    end
    # As a → 0+, length → ∞; as a → ∞, length → L. So f(a)=2a*sinh(L/(2a)) - t_l
    f(a) = 2a*sinh(L/(2a)) - t_l
    lo = eps(Float64)          # near zero
    hi = max(L, t_l)*1e6       # very large a -> length ~ L
    # ensure sign change
    while f(lo) < 0
        lo /= 10
    end
    while f(hi) > 0
        hi *= 10
    end
    _bisect(f, lo, hi)
end

"""
    catenary_analytical_solution(horizontal_span, total_length, line_diameter_mm, set; n_points=1000)

Exact catenary solution for a hanging cable between two fixed points.
Returns analytical x, y coordinates and catenary parameter.
"""
function catenary_analytical_solution(horizontal_span, total_length, line_diameter_mm, set; n_points=1000)
    # Calculate system parameters
    cross_section_area = π * ((1e-3)*line_diameter_mm/2)^2
    mass_line = set.rho_tether * cross_section_area * total_length
    
    println("Catenary System Parameters:")
    println("   Horizontal span: ", horizontal_span, " m")
    println("   Total cable length: ", total_length, " m")
    println("   Line diameter: ", line_diameter_mm, " mm")
    println("   E modulus: ", set.e_tether, " Pa")
    println("   Density: ", set.rho_tether, " kg/m^3")
    println("   Gravity: ", set.g_earth, " m/s^2")
    println("   Cross section area: ", round(cross_section_area, digits=6), " m^2")
    println("   Total cable mass: ", round(mass_line, digits=4), " kg")
    
    # Solve for catenary parameter 'a'
    a = _solve_a_from_length(horizontal_span, total_length)
    
    # Generate analytical catenary curve
    x = range(0.0, horizontal_span; length=n_points)
    y = a .* cosh.((x .- horizontal_span/2) ./ a)
    y .-= y[1]  # set supports to zero height
    
    println("   Catenary parameter a: ", round(a, digits=4), " m")
    println("   Maximum sag: ", round(abs(minimum(y)), digits=4), " m")
    
    return x, y, a
end

println("\n\nCatenary Line Example\n", "="^40)
# --- Settings
set = load_settings(joinpath(@__DIR__, "..", "data", "base", "settings.yaml"))  # Loads as Dict
# Example usage of set as Dict:
# set["v_wind"] = 0.0               # No wind for pure catenary
# set["abs_tol"] = 1e-8
# set["rel_tol"] = 1e-8
# set["l_tether"] = 8        # Set tether length for plot size

# Catenary parameters
horizontal_span = 8          # Horizontal distance between anchors [m]
n_segments = 10                # Number of tether segments (n points = n_segments+1)
total_length = 10.0            # Total unstretched length of tether [m]
compression_frac = 0.01      # Compression fraction for segments
world_frame_damping = 1
line_diameter_mm = 4.0       # Diameter for analytical calculations

# --- Points (nodes)
points = Point[]
push!(points, Point(1, [0.0, 0.0, 5.0], STATIC))                               # Left anchor

# Add dynamic points along initial straight line
for i in 1:n_segments-1
    x = i * horizontal_span / n_segments
    push!(points, Point(i+1, [x, 0.0, 5.0], DYNAMIC; world_frame_damping=world_frame_damping))       # Points 2,3,4,...
end
push!(points, Point(n_segments+1, [horizontal_span, 0.0, 5.0], STATIC))        # Right anchor

# --- Segments (springs/dampers)
segments = Segment[]
l0_per_segment = total_length / n_segments  # Rest length per segment

for i in 1:n_segments
    # Connect consecutive points
    push!(segments, Segment(i, set, (i, i+1), POWER_LINE; l0=l0_per_segment, diameter_mm=line_diameter_mm,
        compression_frac=compression_frac))
end

# --- Transforms
transforms = [Transform(1, 0.0, 0.0, 0.0; base_pos=[0.0, 0.0, 5.0], base_point_idx=1, rot_point_idx=2)]

# --- System structure
sys_struct = SymbolicAWEModels.SystemStructure("catenary", set; points, segments, transforms)


# Analyze the system, to find optimal damping
include("damping_analysis.jl")
freq, zeta, recommended_damping = analyze_damping_response(
    sys_struct, set; verbose=true, perturbation_dir=[0.0, 0.0, 1.0]
    )
println("\n Setting recommended damping to: ", recommended_damping)
set.damping = recommended_damping  # Update settings with recommended damping


# Plot initial state
plot(sys_struct, 0.0; zoom=false)

# --- Construct symbolic model
sam = SymbolicAWEModel(set, sys_struct)
init!(sam; remake=false)

# --- Simulate until static equilibrium
println("Simulating catenary formation...")
n_steps = 500
for i in 1:n_steps
    next_step!(sam)
    # plot(sam, i/set.sample_freq; zoom=false)
    ## Plot every 10 steps to show evolution
    if i % 5 == 0
        plot(sam, i/set.sample_freq; zoom=false)
        # println("Step $i/$n_steps")
    end
end

# --- Final comparison with analytical solution ---
println("\n\nSimulation vs Analytical Comparison\n", "="^50)

### Calculate analytical solution
x_analytical, y_analytical, catenary_a = catenary_analytical_solution(horizontal_span, total_length, line_diameter_mm, set)

# --- Final plot with analytical comparison
plot(sam, n_steps/set.sample_freq; zoom=false)

# Add analytical solution to the plot (if ControlPlots supports it)
try
    # Extract simulation points for plotting
    sim_x = [sam.sys_struct.points[i].pos_w[1] for i in 1:length(points)]
    sim_z = [sam.sys_struct.points[i].pos_w[3] for i in 1:length(points)]
    
    # Create comparison plot data
    println("\nPlotting analytical vs simulation comparison...")
    println("Simulation points: $(length(sim_x))")
    println("Analytical points: $(length(x_analytical))")
    
    # Note: Additional plotting with analytical overlay would require
    # extending the ControlPlots functionality or using a different plotting package
    
catch e
    println("Note: Could not overlay analytical solution on plot: $e")
end
println("Simulation complete. The tether should now show a catenary shape under gravity.")

# --- Extract final positions for validation and comparison ---
println("\nFinal point positions - Simulation vs Analytical:")
println("Point |Simulation (x,z)|Analytical (x,z)|Delta (x,z)")
println("------|----------------|----------------|-----------")

for i in 1:length(points)
    pos = sam.sys_struct.points[i].pos_w
    sim_x, sim_z = pos[1], pos[3]
    
    # Find closest analytical point for comparison
    # For end points, use exact positions; for intermediate points, interpolate
    if i == 1
        anal_x, anal_z = 0.0, 5.0 + y_analytical[1]  # Left anchor
    elseif i == length(points)
        anal_x, anal_z = horizontal_span, 5.0 + y_analytical[end]  # Right anchor
    else
        # Interpolate analytical solution at simulation point x-coordinate
        anal_idx = argmin(abs.(x_analytical .- sim_x))
        anal_x = x_analytical[anal_idx]
        anal_z = 5.0 + y_analytical[anal_idx]  # Add base height of 5.0
    end
    
    delta_x = sim_x - anal_x
    delta_z = sim_z - anal_z
    
    println("  $i   | ($(rpad(round(sim_x, digits=3), 5)), $(rpad(round(sim_z, digits=3), 5))) | ($(rpad(round(anal_x, digits=3), 5)), $(rpad(round(anal_z, digits=3), 5))) | ($(rpad(round(delta_x, digits=3), 5)), $(rpad(round(delta_z, digits=3), 5)))")
end

# Calculate RMS error
global rms_error_x = 0.0
global rms_error_z = 0.0
global n_dynamic_points = 0

for i in 2:(length(points)-1)  # Only dynamic points
    pos = sam.sys_struct.points[i].pos_w
    sim_x, sim_z = pos[1], pos[3]
    
    # Interpolate analytical solution
    anal_idx = argmin(abs.(x_analytical .- sim_x))
    anal_x = x_analytical[anal_idx]
    anal_z = 5.0 + y_analytical[anal_idx]
    
    global rms_error_x += (sim_x - anal_x)^2
    global rms_error_z += (sim_z - anal_z)^2
    global n_dynamic_points += 1
end

global rms_error_x = sqrt(rms_error_x / n_dynamic_points)
global rms_error_z = sqrt(rms_error_z / n_dynamic_points)

println("\nRMS Error (dynamic points only):")
println("  X-direction: $(round(rms_error_x, digits=6)) m")
println("  Z-direction: $(round(rms_error_z, digits=6)) m")
println("  Total RMS:   $(round(sqrt(rms_error_x^2 + rms_error_z^2), digits=6)) m")