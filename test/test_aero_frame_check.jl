# SPDX-FileCopyrightText: 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_aero_frame_check.jl - Test the load-time check that a wing's structural
# points and its aerodynamic sections share one CAD frame
#
# Loads the rigid 2plate kite with every structural point moved by a rigid
# motion, a scale or a mirror, and with a seeded COM, and checks which loads
# are refused.

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

using Test
using LinearAlgebra
using SymbolicAWEModels
using SymbolicAWEModels: VortexStepMethod
using KiteUtils

"""
    move_triple(text, motion)

Return the YAML vector `text`, `"[x, y, z]"`, moved by `motion`.
"""
function move_triple(text, motion)
    pos = motion(parse.(Float64, split(strip(text, ['[', ']']), ",")))
    return "[" * join(pos, ", ") * "]"
end

"""
    move_points(yaml, motion)

Return the structure YAML text with every `pos_cad` in its `points` block replaced
by `motion(pos_cad)`.
"""
function move_points(yaml, motion)
    head, rest = split(yaml, "\npoints:"; limit=2)
    points, tail = split(rest, "\nsegments:"; limit=2)
    number = raw"(-?[0-9.]+)"
    triple = Regex("\\[\\s*$number,\\s*$number,\\s*$number\\s*\\]")
    moved = replace(points, triple => match -> move_triple(match, motion))
    return head * "\npoints:" * moved * "\nsegments:" * tail
end

"""
    load_moved(data_path, set, vsm_set, motion; wing_columns="")

Load the rigid 2plate kite with its points moved by `motion` and `wing_columns`
appended to its wing row.
"""
function load_moved(data_path, set, vsm_set, motion; wing_columns="")
    yaml = read(joinpath(data_path, "rigid_structural_geometry.yaml"), String)
    yaml = replace(move_points(yaml, motion),
        "      aero_z_offset: 0.0\n" => "      aero_z_offset: 0.0\n" * wing_columns)
    path = joinpath(data_path, "moved_structural_geometry.yaml")
    write(path, yaml)
    return load_sys_struct_from_yaml(path; system_name="frame_check", set, vsm_set)
end

@testset "structural and aero frames agree at load" begin
    src_data_path = joinpath(dirname(@__DIR__), "data", "2plate_kite")
    data_path = joinpath(mktempdir(), "2plate_kite")
    cp(src_data_path, data_path; force=true)
    set_data_path(data_path)
    set = Settings("system.yaml")
    vsm_set = VortexStepMethod.VSMSettings(
        joinpath(data_path, "vsm_settings.yaml"); data_prefix=false)
    refused(motion; kw...) = @test_throws r"Wing main_wing: .*CAD frame" load_moved(
        data_path, set, vsm_set, motion; kw...)

    @testset "a matching frame loads, and so does a 0.06 m shift" begin
        @test load_moved(data_path, set, vsm_set, identity) isa SystemStructure
        @test load_moved(data_path, set, vsm_set,
            pos -> pos + [0.06, 0.0, -0.06]) isa SystemStructure
    end

    @testset "points outside the aero bounding box are refused" begin
        refused(pos -> pos + [1.0, 0.0, 0.0])
        refused(pos -> 1000 * pos)
        refused(pos -> [-pos[1], pos[2], -pos[3]])
    end

    @testset "a mirrored span is refused" begin
        refused(pos -> [pos[1], -pos[2], pos[3]])
    end

    @testset "a seeded COM outside the aero bounding box is refused" begin
        unit_inertia = "      unit_inertia: [0.1, 0.1, 0.1, 0.0, 0.0, 0.0]\n"
        @test load_moved(data_path, set, vsm_set, identity;
            wing_columns=unit_inertia * "      com: [0.0, 0.0, 2.4]\n") isa
            SystemStructure
        refused(identity; wing_columns=unit_inertia * "      com: [0.0, 0.0, 0.0]\n")
    end
end
