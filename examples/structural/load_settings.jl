# SPDX-FileCopyrightText: 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

"""
Helper function to load settings from a YAML file path.

This is used by examples that need to load settings directly from a file path.
It ensures the data path is set correctly before loading.
"""
function load_settings(yaml_path::String)
    using KiteUtils
    
    # Ensure data path is set to the project's data directory
    # The yaml_path should be relative to examples/structural/
    # e.g., "../data/base/settings.yaml" or "../../data/base/settings.yaml"
    data_dir = joinpath(@__DIR__, "..", "..", "data")
    
    if isdir(data_dir)
        KiteUtils.set_data_path(data_dir)
    else
        error("Data directory not found at: $data_dir")
    end
    
    # Extract the relative path from the full path
    # The Settings constructor expects a path relative to data_path
    # e.g., if yaml_path is "../../data/base/settings.yaml" -> "base/system.yaml"
    abs_yaml = abspath(yaml_path)
    abs_data = abspath(data_dir)
    
    if !startswith(abs_yaml, abs_data)
        error("Settings file must be within the data directory. Got: $yaml_path")
    end
    
    # Get the relative path from data dir and find the corresponding system.yaml
    rel_path = relpath(abs_yaml, abs_data)
    dir_name = dirname(rel_path)
    system_yaml = joinpath(dir_name, "system.yaml")
    
    # Load using KiteUtils Settings constructor
    return Settings(system_yaml)
end
