# SPDX-FileCopyrightText: 2022, 2024, 2025 Uwe Fechner
# SPDX-License-Identifier: MIT

using SymbolicAWEModels, KiteUtils
using Test

cd("..")
KiteUtils.set_data_path("") 
@testset verbose = true "Testing SymbolicAWEModels..." begin
    println("--> 1")
    include("test_ram_air_kite.jl")
    println("--> 2")
    include("test_helpers.jl")
    println("--> 3")
    include("aqua.jl")
end
