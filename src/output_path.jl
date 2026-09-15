# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

const OUTPUT_PATH = ["output"]

"""
    set_output_path(output_path="output")

Set the directory that simulation results are written to. A relative path is
taken from the working directory, where the default puts them in `output`.
"""
function set_output_path(output_path="output")
    OUTPUT_PATH[1] = output_path
    return nothing
end

"""
    get_output_path()

The directory that simulation results are written to, created if it is missing.
It holds nothing that re-running the simulation would not write again.
"""
function get_output_path()
    return mkpath(OUTPUT_PATH[1])
end
