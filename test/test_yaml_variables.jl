# SPDX-FileCopyrightText: 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_yaml_variables.jl - Test the `variables` block of the YAML loader
#
# Verifies:
# 1. Scalar variables substituted in columns and inside nested lists
# 2. Variables defined in terms of other variables
# 3. Multi-variables filling several columns of a row at once, with a
#    material shared by segments of different diameter
# 4. Multi-variables merged into a dict row
# 5. Errors on cycles, on shadowing a component name, on columns that do not
#    match the fields of a multi-variable and on conflicting stiffness inputs

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

using Test
using SymbolicAWEModels
using KiteUtils

VARIABLES_YAML = """
variables:
  bridle_comp: 0.01
  same_comp: bridle_comp
  top_height: -10.0
  anchor_pos: [1.0, 2.0, 0.0]
  dyneema:
    youngs_modulus: 55.0e9
    damping_per_stiffness: 0.00077
    density: 724.0

points:
  headers: [name, pos_cad, type]
  data:
    - [anchor, anchor_pos, STATIC]
    - [top, [1.0, 2.0, top_height], DYNAMIC]

segments:
  headers: [name, point_i, point_j, l0, diameter_mm,
            youngs_modulus, damping_per_stiffness, density, compression_frac]
  data:
    - [thin_line, anchor, top, 10.0, 1.0, dyneema, bridle_comp]
    - [thick_line, anchor, top, 10.0, 2.0, dyneema, same_comp]
    - [plain_line, anchor, top, 10.0, 2.0, 5.0e9, 0.002, 900.0, same_comp]
"""

DICT_ROW_YAML = """
variables:
  bridle_comp: 0.01
  dyneema:
    youngs_modulus: 55.0e9
    damping_per_stiffness: 0.00077
    density: 724.0

points:
  headers: [name, pos_cad, type]
  data:
    - [anchor, [0.0, 0.0, 0.0], STATIC]
    - [top, [0.0, 0.0, -10.0], DYNAMIC]

segments:
  data:
    - {name: dict_line, point_i: anchor, point_j: top, l0: 10.0,
       diameter_mm: 1.0, material: dyneema, damping_per_stiffness: 0.001,
       compression_frac: bridle_comp}
"""

CYCLIC_YAML = """
variables:
  a: b
  b: a

points:
  headers: [name, pos_cad, type]
  data:
    - [anchor, [0.0, 0.0, 0.0], STATIC]
"""

SHADOW_YAML = """
variables:
  anchor: 1.0

points:
  headers: [name, pos_cad, type]
  data:
    - [anchor, [0.0, 0.0, 0.0], STATIC]
"""

MISMATCH_YAML = """
variables:
  dyneema:
    youngs_modulus: 55.0e9
    damping_per_stiffness: 0.00077

points:
  headers: [name, pos_cad, type]
  data:
    - [anchor, [0.0, 0.0, 0.0], STATIC]
    - [top, [0.0, 0.0, -10.0], DYNAMIC]

segments:
  headers: [name, point_i, point_j, youngs_modulus, compression_frac]
  data:
    - [line, anchor, top, dyneema, 0.01]
"""

CONFLICT_YAML = """
points:
  headers: [name, pos_cad, type]
  data:
    - [anchor, [0.0, 0.0, 0.0], STATIC]
    - [top, [0.0, 0.0, -10.0], DYNAMIC]

segments:
  headers: [name, point_i, point_j, diameter_mm, youngs_modulus,
            unit_stiffness]
  data:
    - [line, anchor, top, 1.0, 55.0e9, 43196.9]
"""

MATERIALS_YAML = """
materials:
  headers: [name, youngs_modulus, density, damping_per_stiffness]
  data:
    - [dyneema, 55e9, 724, 0.00077]

points:
  headers: [name, pos_cad, type]
  data:
    - [anchor, [0.0, 0.0, 0.0], STATIC]
"""

write_yaml = function (dir, name, text)
    path = joinpath(dir, name)
    write(path, text)
    return path
end

@testset "YAML variables" begin
    tmpdir = mktempdir()
    data_path = joinpath(tmpdir, "2plate_kite")
    cp(joinpath(dirname(@__DIR__), "data", "2plate_kite"), data_path;
       force=true)
    set_data_path(data_path)
    set = Settings("system.yaml")

    sys = load_sys_struct_from_yaml(
        write_yaml(tmpdir, "variables.yaml", VARIABLES_YAML);
        system_name="variables_test", set)

    @test sys.points[:anchor].pos_cad ≈ [1.0, 2.0, 0.0]
    @test sys.points[:top].pos_cad ≈ [1.0, 2.0, -10.0]

    # A multi-variable fills the three columns at its position
    thin_line = sys.segments[:thin_line]
    thick_line = sys.segments[:thick_line]
    @test thin_line.unit_stiffness ≈ 55.0e9 * π * 0.0005^2
    @test thin_line.unit_damping ≈ 0.00077 * thin_line.unit_stiffness
    @test thin_line.density ≈ 724.0
    @test thin_line.compression_frac ≈ 0.01

    # The same material on a thicker segment is four times as stiff
    @test thick_line.unit_stiffness ≈ 4 * thin_line.unit_stiffness
    @test thick_line.unit_damping ≈ 4 * thin_line.unit_damping

    # A variable may name another variable
    @test sys.segments[:plain_line].compression_frac ≈ 0.01
    @test sys.segments[:plain_line].unit_stiffness ≈ 5.0e9 * π * 0.001^2

    # In a dict row the fields are merged, without overriding the row
    dict_sys = load_sys_struct_from_yaml(
        write_yaml(tmpdir, "dict_row.yaml", DICT_ROW_YAML);
        system_name="variables_test", set)
    dict_line = dict_sys.segments[:dict_line]
    @test dict_line.unit_stiffness ≈ 55.0e9 * π * 0.0005^2
    @test dict_line.unit_damping ≈ 0.001 * dict_line.unit_stiffness
    @test dict_line.density ≈ 724.0

    @test_throws ErrorException load_sys_struct_from_yaml(
        write_yaml(tmpdir, "cyclic.yaml", CYCLIC_YAML);
        system_name="variables_test", set)
    @test_throws ErrorException load_sys_struct_from_yaml(
        write_yaml(tmpdir, "shadow.yaml", SHADOW_YAML);
        system_name="variables_test", set)
    @test_throws ErrorException load_sys_struct_from_yaml(
        write_yaml(tmpdir, "mismatch.yaml", MISMATCH_YAML);
        system_name="variables_test", set)
    @test_throws ErrorException load_sys_struct_from_yaml(
        write_yaml(tmpdir, "materials.yaml", MATERIALS_YAML);
        system_name="variables_test", set)
    @test_throws ErrorException load_sys_struct_from_yaml(
        write_yaml(tmpdir, "conflict.yaml", CONFLICT_YAML);
        system_name="variables_test", set)
end
nothing
