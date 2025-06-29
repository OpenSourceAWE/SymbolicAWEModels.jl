# Copyright (c) 2022, 2024 Uwe Fechner
# SPDX-License-Identifier: MIT
using REPL.TerminalMenus

options = [
        "ram_air_kite = SIMPLE=false; include(\"ram_air_kite.jl\")",
        "simple_ram_air_kite = SIMPLE=true; include(\"ram_air_kite.jl\")",
        "lin_ram_model = include(\"lin_ram_model.jl\")",
        "quit"
]

function example_menu()
    active = true
    while active
        menu = RadioMenu(options, pagesize=8)
        choice = request("\nChoose function to execute or `q` to quit: ", menu)

        if choice != -1 && choice != length(options)
            eval(Meta.parse(options[choice]))
        else
            println("Left menu. Press <ctrl><d> to quit Julia!")
            active = false
        end
    end
end

example_menu()
