# SPDX-FileCopyrightText: 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_create_vsm_wing.jl - how create_vsm_wing cuts a mesh into aero sections.

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

using Test
using SymbolicAWEModels
using SymbolicAWEModels: VortexStepMethod, create_vsm_wing
using KiteUtils

data_path = joinpath(mktempdir(), "2plate_kite")
cp(joinpath(dirname(@__DIR__), "data", "2plate_kite"), data_path; force=true)
cp(joinpath(pkgdir(VortexStepMethod), "data", "ram_air_kite",
            "ram_air_kite_body.obj"), joinpath(data_path, "wing.obj"))
set_data_path(data_path)

set = Settings("system.yaml")
set.model = "wing.obj"

vsm_set = VortexStepMethod.VSMSettings(joinpath(data_path, "vsm_settings.yaml");
                                       data_prefix=false)
wing_set = vsm_set.wings[1]
wing_set.n_panels = 20
wing_set.spanwise_panel_distribution = VortexStepMethod.LINEAR
wing_set.mesh.n_sections = 4

@testset "mesh is sliced into the sections the wing settings ask for" begin
    vsm_wing = create_vsm_wing(set, vsm_set; prn=false)
    @test vsm_wing.n_unrefined_sections == 4
    @test length(vsm_wing.unrefined_sections) == 4
    @test vsm_wing.n_panels == 20
    @test length(vsm_wing.refined_sections) == 21
    @test isdir(joinpath(data_path, "obj_geometry", "4_sections"))
end
nothing
