# SPDX-FileCopyrightText: 2022, 2024, 2025 Uwe Fechner
# SPDX-License-Identifier: MIT

using KiteModels, KiteUtils
using Test

cd("..")
KiteUtils.set_data_path("") 
@testset verbose = true "Testing KiteModels..." begin
    include("test_orientation.jl")
    println("--> 1")
    include("test_ram_air_kite.jl")
    println("--> 2")
    include("test_helpers.jl")
    println("--> 3")
    include("aqua.jl")
end
