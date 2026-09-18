# SPDX-FileCopyrightText: 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_yaml_fields.jl - Test that the loader rejects what it cannot read
#
# Verifies:
# 1. A field the loader reads is applied, and the same field misspelled errors
#    instead of leaving the component on its default; the retired `idx` column
#    says to rename it to `name`
# 2. An unknown top-level block errors
# 3. A row carrying more values than the table has headers errors

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

using Test
using SymbolicAWEModels
using KiteUtils

FIELDS_YAML = """
points:
  headers: [name, pos_cad, type, extra_mass]
  data:
    - [anchor, [0.0, 0.0, 0.0], STATIC, 0.0]
    - [top, [1.0, 2.0, 10.0], DYNAMIC, 2.5]
"""

TYPO_FIELD_YAML = replace(FIELDS_YAML, "extra_mass]" => "extra_masss]")

IDX_COLUMN_YAML = replace(FIELDS_YAML, "[name," => "[idx,")

UNKNOWN_BLOCK_YAML = FIELDS_YAML * """
pointz:
  headers: [name, pos_cad, type]
  data:
    - [stray, [0.0, 0.0, 1.0], STATIC]
"""

LONG_ROW_YAML = """
points:
  headers: [name, pos_cad, type]
  data:
    - [anchor, [0.0, 0.0, 0.0], STATIC, 0.0]
"""

"""
    load_yaml_text(dir, name, text, set) -> SystemStructure

Write `text` to `dir/name` and load it as a structural YAML.
"""
function load_yaml_text(dir, name, text, set)
    path = joinpath(dir, name)
    write(path, text)
    return load_sys_struct_from_yaml(path; system_name="yaml_fields_test", set)
end

@testset "YAML fields" begin
    tmpdir = mktempdir()
    data_path = joinpath(tmpdir, "2plate_kite")
    cp(joinpath(dirname(@__DIR__), "data", "2plate_kite"), data_path;
       force=true)
    set_data_path(data_path)
    set = Settings("system.yaml")

    sys = load_yaml_text(tmpdir, "fields.yaml", FIELDS_YAML, set)
    @test sys.points[:top].extra_mass ≈ 2.5

    @test_throws "extra_masss" load_yaml_text(
        tmpdir, "typo.yaml", TYPO_FIELD_YAML, set)
    @test_throws "Rename `idx` to `name`" load_yaml_text(
        tmpdir, "idx.yaml", IDX_COLUMN_YAML, set)
    @test_throws "pointz" load_yaml_text(
        tmpdir, "unknown_block.yaml", UNKNOWN_BLOCK_YAML, set)
    @test_throws "headers" load_yaml_text(
        tmpdir, "long_row.yaml", LONG_ROW_YAML, set)
end
nothing
