# # Copyright (c) 2022, 2024 Uwe Fechner
# # SPDX-License-Identifier: MIT
# using REPL.TerminalMenus

# options = [
#         "pyramid_model_aerodynamic = include(\"pyramid_model_aerodynamic.jl\")",
#         "TUDELFT_V3_KITE_aerodynamic = include(\"TUDELFT_V3_KITE_aerodynamic.jl\")",
#         "saddle_form = include(\"saddle_form.jl\")",
#         "tether_deflection_by_wind = include(\"tether_deflection_by_wind.jl\")",
#         "catenary_line = include(\"catenary_line.jl\")",
#         "hanging_mass = include(\"hanging_mass.jl\")",
#         "simple_pulley = include(\"simple_pulley.jl\")",
#         "lin_ram_model = include(\"lin_ram_model.jl\")",
#         "lin_simple_tuned_model = include(\"lin_simple_tuned_model.jl\")",
#         "pulley = include(\"pulley.jl\")",
#         "ram_air_kite = include(\"ram_air_kite.jl\")",
#         "sam_tutorial = include(\"sam_tutorial.jl\")",
#         "simple_lin_model = include(\"simple_lin_model.jl\")",
#         "simple_tuned_model = include(\"simple_tuned_model.jl\")",
#         "tether_props = include(\"tether_props.jl\")",
#         "quit"
# ]

# function example_menu()
#     active = true
#     while active
#         menu = RadioMenu(options, pagesize=8)
#         choice = request("\nChoose function to execute or `q` to quit: ", menu)

#         if choice != -1 && choice != length(options)
#             eval(Meta.parse(options[choice]))
#         else
#             println("Left menu. Press <ctrl><d> to quit Julia!")
#             active = false
#         end
#     end
# end

# example_menu()

# Copyright (c) 2022–2025
# SPDX-License-Identifier: MIT

using REPL.TerminalMenus
using OrderedCollections: OrderedDict
using KiteUtils

"""
    run_examples_menu()

Top-level menu to choose a category (by subdirectory), then a
specific example to run. Use `q` to exit any menu.
"""
function run_examples_menu()
    # # Initialize data path once for all examples
    # data_path = joinpath(@__DIR__, "..", "data")
    # if isdir(data_path)
    #     KiteUtils.set_data_path(data_path)
    # else
    #     @warn "Data directory not found at $data_path"
    # end
    categories = OrderedDict(
        "aerodynamic/" => OrderedDict(
            "pyramid_model" => """include(joinpath(@__DIR__, "examples", "aerodynamic", "pyramid_model.jl"))""",
            "TUDELFT_V3_KITE" => """include(joinpath(@__DIR__, "examples", "aerodynamic", "TUDELFT_V3_KITE.jl"))""",
        ),
        "structural/" => OrderedDict(
            "hanging_mass" => """include(joinpath(@__DIR__, "examples", "structural", "hanging_mass.jl"))""",
            "simple_pulley" => """include(joinpath(@__DIR__, "examples", "structural", "simple_pulley.jl"))""",
            "pulley" => """include(joinpath(@__DIR__, "examples", "structural", "pulley.jl"))""",
            "catenary_line" => """include(joinpath(@__DIR__, "examples", "structural", "catenary_line.jl"))""",
            "saddle_form" => """include(joinpath(@__DIR__, "examples", "structural", "saddle_form.jl"))""",
            "tether_props" => """include(joinpath(@__DIR__, "examples", "structural", "tether_props.jl"))""",
        ),
        "coupled/" => OrderedDict(
            "tether_deflection_by_wind" => """include(joinpath(@__DIR__, "examples", "coupled", "tether_deflection_by_wind.jl"))""",
            "ram_air_kite" => """include(joinpath(@__DIR__, "examples", "coupled", "ram_air_kite.jl"))""",
            "lin_ram_model" => """include(joinpath(@__DIR__, "examples", "coupled", "lin_ram_model.jl"))""",
            "lin_simple_tuned_model" => """include(joinpath(@__DIR__, "examples", "coupled", "lin_simple_tuned_model.jl"))""",
            "simple_lin_model" => """include(joinpath(@__DIR__, "examples", "coupled", "simple_lin_model.jl"))""",
            "simple_tuned_model" => """include(joinpath(@__DIR__, "examples", "coupled", "simple_tuned_model.jl"))""",
        ),
        "tutorials/" => OrderedDict(
            "sam_tutorial" => """include(joinpath(@__DIR__, "examples", "tutorials", "sam_tutorial.jl"))""",
        ),
    )

    main_loop(categories)
end

"""
    main_loop(categories)

Show the top-level category menu and dispatch to the selected category.
"""
function main_loop(categories::OrderedDict{String,<:Any})
    active = true
    while active
        cat_labels = collect(keys(categories))
        main_opts = vcat(cat_labels, "quit")
        choice = pick("Choose a category or `q` to quit:", main_opts)
        if choice == -1 || choice == length(main_opts)   # q or "quit"
            println("Left menu. Press <ctrl><d> to quit Julia!")
            active = false
        else
            cat = main_opts[choice]
            sub_loop(cat, categories[cat])
        end
    end
end

"""
    sub_loop(title, items)

Show the second-level menu (examples within a category).
`items` maps display names → code strings to `eval`.
"""
function sub_loop(title::AbstractString, items::OrderedDict{String,String})
    active = true
    while active
        labels = collect(keys(items))
        opts = vcat(labels, "back", "quit")
        prompt = "\n[$title] Choose example, `back` to return, or `q` to quit:"
        choice = pick(prompt, opts)
        if choice == -1                   # user pressed q
            println("Left menu. Press <ctrl><d> to quit Julia!")
            return
        elseif opts[choice] == "back"
            return
        elseif opts[choice] == "quit"
            println("Left menu. Press <ctrl><d> to quit Julia!")
            exit()  # or `return` if you prefer not to exit Julia here
        else
            # Run the selected example
            code = items[opts[choice]]
            eval(Meta.parse(code))
        end
    end
end

"""
    pick(prompt, options; pagesize=8) -> Int

Helper that shows a `RadioMenu` with `options` and returns the 1-based index
selected, or `-1` if the user pressed `q`.
"""
function pick(prompt::AbstractString, options::Vector{String}; pagesize::Int=8)
    menu = RadioMenu(options; pagesize)
    return request("\n$prompt ", menu)
end

# Entry point - run the examples menu once
run_examples_menu()
