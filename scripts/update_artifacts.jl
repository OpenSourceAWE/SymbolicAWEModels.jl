# SPDX-FileCopyrightText: 2025 Bart van de Lint <bart@vandelint.net>
#
# SPDX-License-Identifier: MPL-2.0

using SHA, Tar, TOML, CodecZlib

# --- Get version from Project.toml ---
project_toml_path = joinpath(@__DIR__, "..", "Project.toml")

if !isfile(project_toml_path)
    error("Project.toml not found at $project_toml_path. Please ensure this script is in a 'build' or 'scripts' directory next to your Project.toml.")
end

# Parse the Project.toml file to extract the version
project_data = TOML.parsefile(project_toml_path)
version = project_data["version"]
println("Found project version: $version")
# ---

data_dir = joinpath(@__DIR__, "..", "data")
artifacts_toml_path = joinpath(@__DIR__, "..", "Artifacts-v1.$(VERSION.minor).toml.default")
artifacts = Dict()

for entry in readdir(data_dir)
    if endswith(entry, ".tar.gz")
        fname = entry
        local_path = joinpath(data_dir, fname)
        println("Processing local asset: $fname")

        # Compute sha256 hash of the local file
        sha256hash = open(local_path) do io
            bytes2hex(sha256(io))
        end

        # Compute git-tree-sha1 hash string of the decompressed tarball contents
        gtree_sha1 = open(local_path) do file_io
            # Call the stream constructor directly on the open file stream
            gz_io = GzipDecompressorStream(file_io)
            try
                Tar.tree_hash(gz_io; algorithm="git-sha1", skip_empty=false)
            finally
                close(gz_io)
            end
        end

        # First remove the extension, then replace dots with underscores for valid TOML keys
        artifact_name = replace(replace(fname, r"\.tar\.gz" => ""), "." => "_")

        # Construct the correct URL using the version from Project.toml
        url = "https://github.com/OpenSourceAWE/SymbolicAWEModels.jl/releases/download/v$version/$fname"

        artifacts[artifact_name] = Dict(
            "git-tree-sha1" => gtree_sha1,
            "download" => [
                Dict(
                    "sha256" => sha256hash,
                    "url" => url
                )
            ]
        )
    end
end

# Write the Artifacts.toml.default file
open(artifacts_toml_path, "w") do io
    TOML.print(io, artifacts)
end

println("Artifacts.toml successfully written at $artifacts_toml_path")
