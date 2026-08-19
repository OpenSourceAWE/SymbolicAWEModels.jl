# Copyright (c) 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# The synthetic surface-aero fixture the AeroPressure test and the cold-start
# inference benchmark both build their model from.

"""
Write a synthetic per-node surface aero fixture (`.dat` contour + Cp/cf CSVs) for
airfoil id 1 into `data_path/airfoils/`. `camber` bends the mean line, so both
surfaces can sit on one side of the chord; `uniform_cp` writes one Cp everywhere,
whose net traction a closed contour has to cancel. Patches the `info_dict` in the
2plate `aero_geometry.yaml` to reference them. Returns the contour node count.
"""
function write_pressure_fixture(data_path; camber=0.0, uniform_cp=nothing)
    afdir = joinpath(data_path, "airfoils")
    mkpath(afdir)
    n_half = 11
    x_top = [0.5 * (1 + cos(pi * (k - 1) / (n_half - 1))) for k in 1:n_half]
    thick(x) = 0.6 * (0.2969 * sqrt(x) - 0.1260 * x - 0.3516 * x^2 +
                      0.2843 * x^3 - 0.1015 * x^4)
    # A canopy's camber can exceed its half thickness, putting both surfaces on
    # one side of the chord line — the case a symmetric section cannot produce.
    mean_line(x) = 4 * camber * x * (1 - x)
    xs = Float64[]
    ys = Float64[]
    for x in x_top
        push!(xs, x); push!(ys, mean_line(x) + thick(x))
    end
    for x in reverse(x_top)[2:end]
        push!(xs, x); push!(ys, mean_line(x) - thick(x))
    end
    n_node = length(xs)
    open(joinpath(afdir, "1.dat"), "w") do io
        println(io, "synthetic")
        for k in 1:n_node
            println(io, "$(round(xs[k], digits=5)) $(round(ys[k], digits=5))")
        end
    end
    alphas_deg = [-30.0, -15.0, 0.0, 15.0, 30.0]
    open(joinpath(afdir, "1_cp.csv"), "w") do io
        println(io, "alpha,delta," * join(["n$(k-1)" for k in 1:n_node], ","))
        for a in alphas_deg
            cp = isnothing(uniform_cp) ?
                [(-3.0 * ys[k] - 0.05 * a * xs[k]) for k in 1:n_node] :
                fill(Float64(uniform_cp), n_node)
            println(io, "$a,0.0," * join(round.(cp, digits=4), ","))
        end
    end
    open(joinpath(afdir, "1_cf.csv"), "w") do io
        println(io, "alpha,delta," * join(["n$(k-1)" for k in 1:n_node], ","))
        for a in alphas_deg
            # Skin friction has no closed-contour identity, so a uniform-Cp
            # fixture drops it and leaves only the pressure term to cancel.
            cf = isnothing(uniform_cp) ? 0.006 : 0.0
            println(io, "$a,0.0," * join(fill(cf, n_node), ","))
        end
    end
    geom = joinpath(data_path, "aero_geometry.yaml")
    txt = read(geom, String)
    txt = replace(txt,
        "[1, polars, {csv_file_path: \"polars/1.csv\"}]" =>
        "[1, polars, {csv_file_path: \"polars/1.csv\", " *
        "dat_file: \"airfoils/1.dat\", cp_file: \"airfoils/1_cp.csv\", " *
        "cf_file: \"airfoils/1_cf.csv\"}]")
    write(geom, txt)
    return n_node
end

