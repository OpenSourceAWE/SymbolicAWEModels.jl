# Copyright (c) 2022-2025
# SPDX-License-Identifier: MIT

using Pkg
if ! ("GLMakie" ∈ keys(Pkg.project().dependencies))
    Pkg.activate(@__DIR__)
end

using GLMakie
using SymbolicAWEModels
using REPL.TerminalMenus

options = [
    "hanging_mass = include(\"hanging_mass.jl\")",
    "catenary_line = include(\"catenary_line.jl\")",
    "pulley = include(\"pulley.jl\")",
    "saddle_form = include(\"saddle_form.jl\")",
    "coupled_2plate_kite = include(\"coupled_2plate_kite.jl\")",
    "coupled_2plate_kite_linear_vsm = include(\"coupled_2plate_kite_linear_vsm.jl\")",
    "coupled_tether_deflection = include(\"coupled_tether_deflection.jl\")",
    "coupled_realtime_visualization = include(\"coupled_realtime_visualization.jl\")",
    "coupled_linearize = include(\"coupled_linearize.jl\")",
    "coupled_simple_lin_model = include(\"coupled_simple_lin_model.jl\")",
    "static_load_2plate_kite = include(\"static_load_2plate_kite.jl\")",
    "sam_tutorial = include(\"sam_tutorial.jl\")",
    "quit",
]

function example_menu()
    active = true
    while active
        menu = RadioMenu(options, pagesize=16)
        choice = request(
            "\nChoose example to run or `q` to quit: ",
            menu)
        if choice != -1 && choice != length(options)
            eval(Meta.parse(options[choice]))
        else
            println("Left menu. Press <ctrl><d> to quit Julia!")
            active = false
        end
    end
end

example_menu()
