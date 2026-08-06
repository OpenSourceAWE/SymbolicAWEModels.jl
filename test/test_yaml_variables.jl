# SPDX-FileCopyrightText: 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_yaml_variables.jl - Test the `variables` block of the YAML loader
#
# Verifies:
# 1. Scalar variables substituted in columns and inside nested lists
# 2. Variables defined in terms of other variables
# 3. Multi-variables filling several columns of a row at once
# 4. Multi-variables merged into a dict row
# 5. Errors on cycles, on shadowing a component name and on columns that
#    do not match the fields of a multi-variable

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
    unit_stiffness: 43196.9
    unit_damping: 33.2
    density: 724.0

points:
  headers: [name, pos_cad, type]
  data:
    - [anchor, anchor_pos, STATIC]
    - [top, [1.0, 2.0, top_height], DYNAMIC]

segments:
  headers: [name, point_i, point_j, l0, diameter_mm,
            unit_stiffness, unit_damping, density, compression_frac]
  data:
    - [material_line, anchor, top, 10.0, 1.0, dyneema, bridle_comp]
    - [plain_line, anchor, top, 10.0, 2.0, 5000.0, 10.0, 900.0, same_comp]
"""

DICT_ROW_YAML = """
variables:
  bridle_comp: 0.01
  dyneema:
    unit_stiffness: 43196.9
    unit_damping: 33.2
    density: 724.0

points:
  headers: [name, pos_cad, type]
  data:
    - [anchor, [0.0, 0.0, 0.0], STATIC]
    - [top, [0.0, 0.0, -10.0], DYNAMIC]

segments:
  data:
    - {name: dict_line, point_i: anchor, point_j: top, l0: 10.0,
       diameter_mm: 1.0, material: dyneema, unit_damping: 7.0,
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
    unit_stiffness: 43196.9
    unit_damping: 33.2

points:
  headers: [name, pos_cad, type]
  data:
    - [anchor, [0.0, 0.0, 0.0], STATIC]
    - [top, [0.0, 0.0, -10.0], DYNAMIC]

segments:
  headers: [name, point_i, point_j, unit_stiffness, compression_frac]
  data:
    - [line, anchor, top, dyneema, 0.01]
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
    material_line = sys.segments[:material_line]
    @test material_line.unit_stiffness ≈ 43196.9
    @test material_line.unit_damping ≈ 33.2
    @test material_line.density ≈ 724.0
    @test material_line.compression_frac ≈ 0.01

    # A variable may name another variable
    @test sys.segments[:plain_line].compression_frac ≈ 0.01
    @test sys.segments[:plain_line].unit_stiffness ≈ 5000.0

    # In a dict row the fields are merged, without overriding the row
    dict_sys = load_sys_struct_from_yaml(
        write_yaml(tmpdir, "dict_row.yaml", DICT_ROW_YAML);
        system_name="variables_test", set)
    dict_line = dict_sys.segments[:dict_line]
    @test dict_line.unit_stiffness ≈ 43196.9
    @test dict_line.unit_damping ≈ 7.0
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
end
nothing
