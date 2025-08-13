# Copyright (c) 2025 Uwe Fechner, Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

"""
    create_model_archive(source_dir, archive_path)

Finds all `model*.bin` files in the `source_dir`, copies them to a temporary
directory, and compresses that directory into a `.tar.gz` archive at the
specified `archive_path`.
"""
function create_model_archive(source_dir, archive_path; prn=true)
    # Find all files matching the pattern "model*.bin"
    version = VERSION.minor
    model_files = filter(
        x -> startswith(x, "model_1.$version") && endswith(x, ".bin"),
        readdir(source_dir)
    )
    if isempty(model_files)
        prn && @warn "No 'model*.bin' files found in '$source_dir'. Archive will be empty."
        return
    end
    mktempdir() do tmp_dir
        prn && @info "Staging files for compression in: $tmp_dir"
        # Copy the relevant .bin files to the temporary directory
        for file_name in model_files
            full_path = joinpath(source_dir, file_name)
            cp(full_path, joinpath(tmp_dir, file_name))
        end
        # Create the .tar.gz archive from the temporary directory
        prn && @info "Compressing files into archive: $archive_path"
        open(archive_path, "w") do io
            stream = GzipCompressorStream(io)
            Tar.create(tmp_dir, stream)
            close(stream)
        end
    end
end

"""
    extract_model_archive(archive_path, dest_dir)

Safely decompress a `.tar.gz` file by first extracting to a temporary
directory and then copying the contents to the final destination.

# Arguments
- `archive_path::String`: The path to the `.tar.gz` file to be extracted.
- `dest_dir::String`: The path to the target directory.
"""
function extract_model_archive(archive_path::String, dest_dir::String; prn=true)
    if !isfile(archive_path)
        error("Archive file not found: $archive_path")
    end
    prn && @info "Extracting '$archive_path' to '$dest_dir'..."
    mktempdir() do temp_dir
        # 1. Extract the archive to the temporary directory
        open(archive_path) do io
            stream = GzipDecompressorStream(io)
            Tar.extract(stream, temp_dir)
            close(stream)
        end
        # 2. Copy the extracted contents to the final destination
        for item in readdir(temp_dir)
            source_path = joinpath(temp_dir, item)
            dest_path = joinpath(dest_dir, item)
            cp(source_path, dest_path, force=true)
        end
    end
    prn && @info "Extraction complete."
end

"""
    filecmp(path1::AbstractString, path2::AbstractString) -> Bool

Compare two files byte-by-byte to check if they are identical.
"""
function filecmp(path1::AbstractString, path2::AbstractString)
    stat1, stat2 = stat(path1), stat(path2)
    if !(isfile(stat1) && isfile(stat2)) || filesize(stat1) != filesize(stat2)
        return false
    end
    open(path1, "r") do file1
        open(path2, "r") do file2
            buf1 = Vector{UInt8}(undef, 32768)
            buf2 = similar(buf1)
            while !eof(file1) && !eof(file2)
                n1 = readbytes!(file1, buf1)
                n2 = readbytes!(file2, buf2)
                if n1 != n2 || buf1[1:n1] != buf2[1:n2]
                    return false
                end
            end
            return eof(file1) && eof(file2)
        end
    end
end

"""
    create_default_models(; prn=true)

Create and initialize a set of default `SymbolicAWEModel` instances for precompilation.
"""
function create_default_models(; prn=true)
    function create_model(name; segments=3)
        set = Settings("system.yaml")
        set.segments = segments
        set.physical_model = name
        s = SymbolicAWEModel(set)
        init!(s; prn)
        return s
    end
    sam = create_model("ram")
    tether_sam = create_model("tether")
    simple_sam = create_model("simple_ram")
    one_seg_sam = create_model("ram"; segments=1)
    one_seg_tether_sam = create_model("tether"; segments=1)
    return sam, tether_sam, simple_sam, one_seg_sam, one_seg_tether_sam
end

@setup_workload begin
    using Pkg.Artifacts

    local will_precompile
    will_precompile = get(ENV, "SAM_PRECOMPILE", "true") != "false"
    if will_precompile
        try
            path = dirname(dirname(pathof(@__MODULE__)))
            data_path = joinpath(path, "data")
            set_data_path(data_path)

            # Copy .default files to expected files
            cp(joinpath(path, "Artifacts.toml.default"), joinpath(path, "Artifacts.toml"); force=true)

            version_minor = VERSION.minor
            manifest_default = joinpath(path, "Manifest-v1.$version_minor.toml.default")
            manifest_actual = joinpath(path, "Manifest-v1.$version_minor.toml")
            cp(manifest_default, manifest_actual; force=true)

            # Use explicit Artifacts API to get the artifact and extract it
            artifact_toml = joinpath(path, "Artifacts.toml")
            artifact_name = "models_v1_$version_minor"
            model_hash = artifact_hash(artifact_name, artifact_toml)

            if !isnothing(model_hash) && artifact_exists(model_hash)
                model_dir = artifact_path(model_hash)
                mkpath(data_path)
                for f in readdir(model_dir)
                    cp(joinpath(model_dir, f), joinpath(data_path, f); force=true)
                end
                @info "Downloaded and extracted $artifact_name to $data_path"
            else
                @warn "Artifact $artifact_name not found in Artifacts.toml or does not exist!"
                will_precompile = false
            end
        catch e
            will_precompile = false
            @info "Not precompiling because of error: $e"
        end
    end

    @compile_workload begin
        if will_precompile
            prn=true
            sam, tether_sam, simple_sam, _, _ = create_default_models(; prn)
            init!(sam; prn=false, reload=true)
            init!(sam; prn=false, reload=false)
            sim_oscillate!(sam; total_time=1.0)
            copy_to_simple!(sam, tether_sam, simple_sam; prn=false)
            find_steady_state!(sam)
            ss = SysState(sam)
            next_step!(sam)
            update_sys_state!(ss, sam)
        end
    end
end
