# SPDX-FileCopyrightText: 2026 Uwe Fechner
# SPDX-License-Identifier: MPL-2.0

ENV["MPLBACKEND"] = "Agg"

using Pkg
if Base.active_project() != joinpath(@__DIR__, "Project.toml")
    Pkg.activate(@__DIR__)
end

using KiteUtils
using REPL.TerminalMenus
using SymbolicAWEModels

const EXCLUDE = Set(["test_for_precompile.jl", "test_menu.jl"])

function setup_test_data_path()
    pkg_root = dirname(@__DIR__)
    src_data_path = joinpath(pkg_root, "data", "2plate_kite")
    tmpdir = mktempdir()
    data_path = joinpath(tmpdir, "2plate_kite")
    cp(src_data_path, data_path; force=true)
    set_data_path(data_path)
end

function collect_test_files()
    files = filter(readdir(@__DIR__)) do f
        startswith(f, "test_") && endswith(f, ".jl") && f ∉ EXCLUDE
    end
    sort!(files)
    return files
end

function run_test_file(test_file)
    println("\n--> Running $test_file")
    setup_test_data_path()
    include(test_file)
end

function test_menu()
    test_files = collect_test_files()

    if isempty(test_files)
        println("No test files found.")
        return nothing
    end

    active = true
    while active
        options = vcat(test_files, ["quit"])
        menu = RadioMenu(options, pagesize=14)
        choice = request("\nChoose test to run or `q` to quit: ", menu)

        if choice != -1 && choice <= length(test_files)
            run_test_file(test_files[choice])
        else
            println("Left menu. Press <ctrl><d> to quit Julia!")
            active = false
        end
    end

    return nothing
end

if isinteractive() || abspath(PROGRAM_FILE) == abspath(@__FILE__)
    test_menu()
end
