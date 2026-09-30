# Copyright (c) 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# The runner side of the job API in `api/openapi.yaml`: heartbeat, claim, run, upload.

const HEARTBEAT_INTERVAL = 20.0  # [s], stated in api/openapi.yaml
const CLAIM_INTERVAL = 5.0       # [s] between claims while nothing is queued
const PROGRESS_INTERVAL = 5.0    # [s] between progress reports of one job
const REQUEST_TIMEOUT = 60.0     # [s] for every call to the site but an upload
const REPORT_ATTEMPTS = 5
const REPORT_RETRY_PAUSE = 5.0   # [s]
const FAILURE_LOG_LINES = 60
const ARROW_TYPE = "application/vnd.apache.arrow.file"

"""A running job an admin cancelled."""
struct JobCancelled <: Exception end

"""
    JobRunner(kites; site, token, slots=1)

The state of [`serve_jobs`](@ref): the kites it serves by name, the site and token it
calls with, and a cancellation flag per running job. Setting `stopped[]` ends its
claim loop.
"""
struct JobRunner
    kites::Dict{String, JobKite}
    site::String
    token::String
    slots::Int
    id::String
    version::Dict{String, Any}
    cancelled::Dict{String, Threads.Atomic{Bool}}
    lock::ReentrantLock
    stopped::Threads.Atomic{Bool}
end

function JobRunner(kites; site, token, slots=1)
    names = sort([kite.name for kite in kites])
    allunique(names) || throw(ArgumentError("two kites share a name: $names"))
    return JobRunner(Dict(kite.name => kite for kite in kites), rstrip(site, '/'), token,
                     slots, join(names, "+"), runner_version(kites),
                     Dict{String, Threads.Atomic{Bool}}(), ReentrantLock(),
                     Threads.Atomic{Bool}(false))
end

"""
    serve_jobs(kites; site, token, slots=1)

Run the jobs `site` hands out for `kites`, a vector of [`JobKite`](@ref)s, at most
`slots` at a time, until interrupted. A heartbeat goes out every 20 s from the
interactive thread pool. Each job's parameters are checked against its kite's menu
before the kite runs, and each run is uploaded with a `run` key beside `topology` in
its metadata. A job that fails is reported to the site, and the runner carries on.
"""
serve_jobs(kites; site, token, slots=1) =
    serve_jobs(JobRunner(kites; site, token, slots))

function serve_jobs(runner::JobRunner)
    heartbeat!(runner)
    Threads.@spawn :interactive keep_beating(runner)
    free_slots = Base.Semaphore(runner.slots)
    try
        while !runner.stopped[]
            Base.acquire(free_slots)
            job = claim(runner)
            if isnothing(job)
                Base.release(free_slots)
                sleep(CLAIM_INTERVAL)
            else
                Threads.@spawn run_in_slot(runner, job, free_slots)
            end
        end
    finally
        runner.stopped[] = true
    end
    return nothing
end

"""The package versions `kites` compute with: Julia, SymbolicAWEModels and theirs."""
function runner_version(kites)
    julia_commit = Base.GIT_VERSION_INFO.commit
    version = Dict{String, Any}("julia" => Dict(
        "version" => string(VERSION),
        "commit" => isempty(julia_commit) ? nothing : julia_commit))
    kite_modules = [Base.moduleroot(parentmodule(run))
                    for kite in kites
                    for run in (kite.steady, kite.simulate, kite.sweep_point)]
    for mod in unique([SymbolicAWEModels; kite_modules])
        dir = pkgdir(mod)
        isnothing(dir) && continue
        version[string(nameof(mod))] = Dict("version" => string(pkgversion(mod)),
                                            "commit" => git_commit(dir))
    end
    return version
end

function git_commit(dir)
    (ispath(joinpath(dir, ".git")) && !isnothing(Sys.which("git"))) || return nothing
    return readchomp(`git -C $dir rev-parse HEAD`)
end

"""
Call the site with `body` as JSON or the contents of `file`; the parsed JSON reply,
or `nothing` for an empty one.
"""
function site_call(runner::JobRunner, method, path; body=nothing, file=nothing,
                   content_type="application/json", timeout=REQUEST_TIMEOUT)
    headers = ["Authorization" => "Bearer $(runner.token)"]
    input = if !isnothing(file)
        push!(headers, "Content-Type" => content_type,
              "Content-Length" => string(filesize(file)))
        open(file)
    elseif !isnothing(body)
        push!(headers, "Content-Type" => content_type)
        IOBuffer(JSON.json(body))
    end
    output = IOBuffer()
    response = try
        Downloads.request(runner.site * path; method, headers, input, output, timeout)
    finally
        isnothing(input) || close(input)
    end
    reply = String(take!(output))
    200 <= response.status < 300 ||
        error("$method $path: HTTP $(response.status) $reply")
    return isempty(reply) ? nothing : JSON.parse(reply)
end

"""
Post a report that ends a job, its result or its failure, trying again while the
connection fails.
"""
function send_report(runner::JobRunner, path; kwargs...)
    for attempt in 1:REPORT_ATTEMPTS
        try
            return site_call(runner, "POST", path; kwargs...)
        catch err
            (err isa Downloads.RequestError && attempt < REPORT_ATTEMPTS) || rethrow()
            @warn "Sending $path failed; trying again" exception=err
            sleep(REPORT_RETRY_PAUSE)
        end
    end
end

function heartbeat!(runner::JobRunner)
    running = @lock runner.lock collect(keys(runner.cancelled))
    kites = [Dict("name" => kite.name, "menu" => menu_document(kite))
             for kite in values(runner.kites)]
    body = Dict("version" => runner.version, "kites" => kites,
                "slots" => runner.slots, "running" => running)
    reply = try
        site_call(runner, "PUT", "/api/runners/$(runner.id)"; body)
    catch err
        @warn "Heartbeat failed" exception=err
        return nothing
    end
    @lock runner.lock for id in reply["cancel"]
        haskey(runner.cancelled, id) && (runner.cancelled[id][] = true)
    end
    return nothing
end

function keep_beating(runner::JobRunner)
    while true
        sleep(HEARTBEAT_INTERVAL)
        runner.stopped[] && return nothing
        heartbeat!(runner)
    end
end

function claim(runner::JobRunner)
    body = Dict("runner" => runner.id, "kites" => sort(collect(keys(runner.kites))))
    try
        return site_call(runner, "POST", "/api/jobs/claim"; body)
    catch err
        @warn "Claiming a job failed" exception=err
        return nothing
    end
end

function run_in_slot(runner::JobRunner, job, free_slots)
    try
        run_job(runner, job)
    finally
        Base.release(free_slots)
    end
end

"""
Progress reporter handed to a kite's `simulate`: posts at most every
`PROGRESS_INTERVAL`, and throws [`JobCancelled`](@ref) once the job is cancelled.
"""
mutable struct JobProgress
    const runner::JobRunner
    const id::String
    const cancelled::Threads.Atomic{Bool}
    const started::Float64
    reported::Float64
end

function (progress::JobProgress)(sim_time=nothing; phase="simulating")
    progress.cancelled[] && throw(JobCancelled())
    time() - progress.reported < PROGRESS_INTERVAL && return nothing
    progress.reported = time()
    body = Dict{String, Any}("wall_time" => time() - progress.started, "phase" => phase)
    isnothing(sim_time) || (body["sim_time"] = sim_time)
    reply = try
        site_call(progress.runner, "POST", "/api/jobs/$(progress.id)/progress"; body)
    catch err
        @warn "Reporting progress of job $(progress.id) failed" exception=err
        return nothing
    end
    reply["cancel"] && throw(JobCancelled())
    return nothing
end

function run_job(runner::JobRunner, job)
    id = job["id"]
    cancelled = Threads.Atomic{Bool}(false)
    @lock runner.lock runner.cancelled[id] = cancelled
    @info "Job $id: $(job["kind"]) of $(job["kite"])"
    try
        progress = JobProgress(runner, id, cancelled, time(), -Inf)
        mktempdir(dir -> run_and_upload(runner, job, progress, dir))
        @info "Job $id done"
    catch err
        report_failure(runner, id, err, catch_backtrace())
    finally
        @lock runner.lock delete!(runner.cancelled, id)
    end
    return nothing
end

function run_and_upload(runner::JobRunner, job, progress::JobProgress, dir)
    file, content_type = compute_job(runner, job, progress, dir)
    send_report(runner, "/api/jobs/$(job["id"])/result"; file, content_type, timeout=Inf)
    return nothing
end

"""Run `job` and write its result into `dir`; the file to upload and its type."""
function compute_job(runner::JobRunner, job, progress::JobProgress, dir)
    kite = get(runner.kites, job["kite"], nothing)
    isnothing(kite) && throw(JobRejected("this runner does not serve $(job["kite"])"))
    kind = job["kind"]
    kind == "sweep" && return run_sweep(runner, kite, job, progress, dir)
    kind in ("steady", "sim") || throw(JobRejected("a runner does not run $kind jobs"))
    menu = kind == "steady" ? kite.steady_menu : kite.sim_menu
    parameters = menu_values(menu, get(job, "parameters", Dict()))
    key = job_key(job["sandbox"], kite.name, kind, parameters, runner.version)
    started = time()
    log = kind == "steady" ? kite.steady(parameters) :
          kite.simulate(parameters; progress)
    run = run_metadata(runner, job, parameters, key, time() - started)
    return write_run(dir, "run", log, run), ARROW_TYPE
end

function run_sweep(runner::JobRunner, kite::JobKite, job, progress::JobProgress, dir)
    axes = sweep_axes(kite, get(job, "axes", Dict()))
    fixed = menu_values(kite.sweep_menu, get(job, "parameters", Dict());
                        swept=collect(keys(axes)))
    key = job_key(job["sandbox"], kite.name, "sweep", fixed, runner.version; axes)
    points = vec(collect(Iterators.product(values(axes)...)))
    manifest_points = Dict{String, Any}[]
    for (i, point) in enumerate(points)
        progress(; phase="point $i of $(length(points))")
        parameters = merge(fixed, Dict(zip(keys(axes), point)))
        started = time()
        log = kite.sweep_point(parameters)
        run = run_metadata(runner, job, parameters, key, time() - started)
        name = "point_" * lpad(i, ndigits(length(points)), '0')
        write_run(dir, name, log, run)
        push!(manifest_points, Dict("parameters" => parameters, "run" => "$name.arrow"))
    end
    manifest = Dict("axes" => axes, "parameters" => fixed, "points" => manifest_points)
    return write_multipart(dir, manifest)
end

"""The `run` metadata key of a run, as JSON."""
function run_metadata(runner::JobRunner, job, parameters, key, wall_time)
    return JSON.json(Dict("kind" => job["kind"], "kite" => job["kite"],
                          "parameters" => parameters, "job_key" => key,
                          "version" => runner.version, "wall_time" => wall_time))
end

"""Save `log` as `dir/name.arrow`, uncompressed, with `run` beside its topology."""
function write_run(dir, name, log::SysLog, run)
    haskey(log.metadata, "topology") ||
        error("the kite's log carries no structure document under `topology`; " *
              "build it with `sys_log(logger, sys)`")
    log.name = name
    save_log(log, false; path=dir, metadata=merge(log.metadata, Dict("run" => run)))
    return joinpath(dir, "$name.arrow")
end

"""A sweep's upload in `dir`: its manifest and the run each point names."""
function write_multipart(dir, manifest)
    boundary = "symawe-" * bytes2hex(rand(UInt8, 12))
    path = joinpath(dir, "result.multipart")
    open(io -> write_parts(io, boundary, dir, manifest), path, "w")
    return path, "multipart/form-data; boundary=$boundary"
end

function write_parts(io, boundary, dir, manifest)
    write_part(io, boundary, "manifest", "manifest.json", "application/json",
               JSON.json(manifest))
    for point in manifest["points"]
        write_part(io, boundary, "runs", point["run"], ARROW_TYPE,
                   read(joinpath(dir, point["run"])))
    end
    write(io, "--$boundary--\r\n")
    return nothing
end

function write_part(io, boundary, name, filename, content_type, content)
    write(io, "--$boundary\r\n",
          "Content-Disposition: form-data; name=\"$name\"; filename=\"$filename\"\r\n",
          "Content-Type: $content_type\r\n\r\n")
    write(io, content)
    write(io, "\r\n")
    return nothing
end

function report_failure(runner::JobRunner, id, err, backtrace)
    if err isa JobCancelled
        @info "Job $id cancelled"
        return nothing
    end
    if err isa JobRejected
        reason, log = "rejected", err.msg
        @warn "Job $id rejected: $log"
    else
        reason = "error"
        lines = split(sprint(showerror, err, backtrace), '\n')
        log = join(first(lines, FAILURE_LOG_LINES), '\n')
        @error "Job $id failed" exception=(err, backtrace)
    end
    try
        send_report(runner, "/api/jobs/$id/failure";
                    body=Dict("reason" => reason, "log" => log))
    catch report_err
        @error "Reporting the failure of job $id failed" exception=report_err
    end
    return nothing
end
