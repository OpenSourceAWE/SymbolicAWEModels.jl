# Copyright (c) 2025 Bart van de Lint, Jelle Poland
# SPDX-License-Identifier: LGPL-3.0-only

data_dir = joinpath(@__DIR__, "..", "data", "v3", "2D_polars_CFD_extrapolated")

function read_polar_csv(file_path)
    lines = readlines(file_path)
    header = split(lines[1], ",")
    
    data = []
    for line in lines[2:end]
        vals = split(line, ",")
        push!(data, (
            alpha=parse(Float64, vals[1]),
            Cd=parse(Float64, vals[2]),
            Cs=parse(Float64, vals[3]),
            Cl=parse(Float64, vals[4]),
            Cm=parse(Float64, vals[5])
        ))
    end
    return data
end

function write_polar_csv(file_path, data)
    open(file_path, "w") do f
        write(f, "alpha,Cd,Cs,Cl,Cm\n")
        for row in data
            write(f, "$(row.alpha),$(row.Cd),$(row.Cs),$(row.Cl),$(row.Cm)\n")
        end
    end
end

# Process each file
for i in 1:19
    file_path = joinpath(data_dir, "$i.csv")
    
    # Read the CSV
    data = read_polar_csv(file_path)
    
    # Get the last angle
    last_alpha = data[end].alpha
    
    if last_alpha >= 40.0
        @info "File $i already covers up to 40°, skipping."
        continue
    end
    
    # Extract the last 2 rows for slope calculation
    row_n = data[end]
    row_n_minus_1 = data[end-1]
    
    # Calculate slopes for CD and CM (CL stays constant)
    cd_slope = (row_n.Cd - row_n_minus_1.Cd) / (row_n.alpha - row_n_minus_1.alpha)
    cm_slope = (row_n.Cm - row_n_minus_1.Cm) / (row_n.alpha - row_n_minus_1.alpha)
    cl_constant = row_n.Cl
    
    # Generate new rows from last_alpha + 0.5 to 40.0 in 0.5° increments
    new_alphas = (last_alpha + 0.5):0.5:40.0
    
    for alpha in new_alphas
        cd_new = row_n.Cd + cd_slope * (alpha - row_n.alpha)
        cm_new = row_n.Cm + cm_slope * (alpha - row_n.alpha)
        cl_new = cl_constant
        cs_new = row_n.Cs  # Keep Cs as is
        
        push!(data, (alpha=alpha, Cd=cd_new, Cs=cs_new, Cl=cl_new, Cm=cm_new))
    end
    
    # Write back
    write_polar_csv(file_path, data)
    
    @info "Extrapolated file $i.csv from $(last_alpha)° to 40.0°"
end

@info "All polar files extrapolated successfully!"
