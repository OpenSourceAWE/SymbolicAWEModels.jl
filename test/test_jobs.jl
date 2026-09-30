# SPDX-FileCopyrightText: 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

using Test
using JSON
using JSONSchema
using KiteUtils
using Sockets
using SymbolicAWEModels
using YAML
using KiteUtils: init!, next_step!, update_sys_state!

"""A request the fake site received."""
struct SiteRequest
    method::String
    path::String
    headers::Dict{String, String}
    body::Vector{UInt8}
end

"""
The site's side of `api/openapi.yaml`, in memory, one request per connection: a
queue, what the runner sent, and the jobs an admin cancelled. Every JSON body is
checked against the spec's schema. Setting `bare_replies[]` answers heartbeats and
progress reports with an empty 200.
"""
struct FakeSite
    schemas::Dict{String, Schema}
    queue::Vector{Dict{String, Any}}
    heartbeats::Vector{Any}
    progress::Vector{Any}
    results::Dict{String, SiteRequest}
    failures::Dict{String, Any}
    cancelled::Set{String}
    bare_replies::Threads.Atomic{Bool}
    lock::ReentrantLock
end

function FakeSite()
    api_dir = joinpath(dirname(@__DIR__), "api")
    components = YAML.load_file(joinpath(api_dir, "openapi.yaml"))["components"]
    schemas = Dict(name => Schema(Dict("\$ref" => "#/components/schemas/$name",
                                       "components" => components);
                                  parent_dir=api_dir)
                   for name in ("Heartbeat", "Claim", "Progress", "Failure", "Job",
                                "JobStatus", "RunMetadata", "SweepManifest"))
    return FakeSite(schemas, [], [], [], Dict(), Dict(), Set(), Threads.Atomic{Bool}(false),
                    ReentrantLock())
end

function conforming_body(site::FakeSite, request::SiteRequest, schema_name)
    body = JSON.parse(String(copy(request.body)))
    failure = JSONSchema.validate(site.schemas[schema_name], body)
    isnothing(failure) || error("$(request.path) does not conform: $failure")
    return body
end

function serve_site(site::FakeSite, listener)
    while isopen(listener)
        socket = try
            accept(listener)
        catch
            return nothing
        end
        Threads.@spawn answer(site, socket)
    end
end

function answer(site::FakeSite, socket)
    method, path = split(readline(socket))
    headers = Dict{String, String}()
    for line in eachline(socket)
        isempty(line) && break
        name, value = split(line, ':'; limit=2)
        headers[lowercase(name)] = strip(value)
    end
    get(headers, "expect", "") == "100-continue" &&
        write(socket, "HTTP/1.1 100 Continue\r\n\r\n")
    content_length = parse(Int, get(headers, "content-length", "0"))
    body = read(socket, content_length)
    sizeof(body) == content_length || return close(socket)
    status, reply = handle(site, SiteRequest(method, path, headers, body))
    write(socket, "HTTP/1.1 $status -\r\nContent-Type: application/json\r\n",
          "Content-Length: $(sizeof(reply))\r\nConnection: close\r\n\r\n", reply)
    close(socket)
end

function handle(site::FakeSite, request::SiteRequest)
    get(request.headers, "authorization", "") == "Bearer runner-token" ||
        return 401, ""
    path = request.path
    @lock site.lock if startswith(path, "/api/runners/")
        push!(site.heartbeats, conforming_body(site, request, "Heartbeat"))
        site.bare_replies[] && return 200, ""
        return 200, JSON.json(Dict("cancel" => collect(site.cancelled)))
    elseif path == "/api/jobs/claim"
        kites = conforming_body(site, request, "Claim")["kites"]
        index = findfirst(job -> job["kite"] in kites, site.queue)
        isnothing(index) && return 204, ""
        return 200, JSON.json(popat!(site.queue, index))
    end
    id = split(path, '/')[4]
    @lock site.lock if endswith(path, "/progress")
        push!(site.progress, (id, conforming_body(site, request, "Progress")))
        site.bare_replies[] && return 200, ""
        return 200, JSON.json(Dict("cancel" => id in site.cancelled))
    elseif endswith(path, "/result")
        site.results[id] = request
        return 201, JSON.json(Dict("id" => id, "state" => "done"))
    elseif endswith(path, "/failure")
        site.failures[id] = conforming_body(site, request, "Failure")
        return 204, ""
    end
    return 404, ""
end

"""The parts of a multipart upload, by filename."""
function multipart_parts(request::SiteRequest)
    boundary = codeunits("\r\n--" * split(request.headers["content-type"], "=")[2])
    body = [codeunits("\r\n"); request.body]
    delimiters = findall(boundary, body)
    parts = Dict{String, Vector{UInt8}}()
    for (before, after) in zip(delimiters[1:end-1], delimiters[2:end])
        part = body[before.stop+1:after.start-1]
        head_end = findfirst(codeunits("\r\n\r\n"), part)
        filename = match(r"filename=\"([^\"]+)\"", String(part[1:head_end.start]))[1]
        parts[filename] = part[head_end.stop+1:end]
    end
    return parts
end

"""
Queue a job whose `parameters` hold its whole menu, keyed as the site keys it for
`runner`, or with `key` instead.
"""
function queue!(site::FakeSite, runner, id, kite, kind; parameters=Dict(), axes=nothing,
                key=job_key("open", kite, kind, parameters, runner.version; axes))
    job = Dict{String, Any}("id" => id, "sandbox" => "open", "kite" => kite,
                            "kind" => kind, "key" => key, "parameters" => parameters)
    isnothing(axes) || (job["axes"] = axes)
    @lock site.lock push!(site.queue, job)
    return job
end

"""Wait until `done()` holds, for at most `timeout` seconds."""
function wait_until(done; timeout=120)
    started = time()
    while !done()
        time() - started > timeout && error("timed out")
        sleep(0.2)
    end
    return true
end

"""A run the runner uploaded, read back with `load_log`."""
function read_run(dir, name, bytes)
    write(joinpath(dir, "$name.arrow"), bytes)
    return load_log(name; path=dir)
end

"""Log `steps` steps of the hanging mass at `params["mass"]`, as a job runs it."""
function fly!(sam, kite_runs, params, steps; progress=time -> nothing)
    kite_runs[] += 1
    sam.sys_struct.points[:mass].extra_mass = params["mass"]
    init!(sam; prn=false)
    logger = Logger(sam, steps)
    state = SysState(sam)
    for step in 1:steps
        next_step!(sam)
        update_sys_state!(state, sam)
        state.time = step / sam.set.sample_freq
        log!(logger, state)
        progress(state.time)
    end
    return sys_log(logger, sam.sys_struct, "hanging_mass")
end

"""The hanging mass of `test_for_precompile.jl` as a kite with `mass` on its menus;
`kite_runs` counts its runs."""
function hanging_mass_kite(set, kite_runs)
    points = [Point(:anchor, [2, 0, 5], STATIC),
              Point(:mass, [2, 0, 2], DYNAMIC; extra_mass=1.0)]
    segments = [Segment(:spring, :anchor, :mass, 500.0, 50.0, 0.005; l0=4.0)]
    transforms = [Transform(:tf, -deg2rad(90), 0, 0; base_pos=[2, 0, 5],
                            base_point=:anchor, rot_point=:mass)]
    sys = SystemStructure("hanging_mass", set; points, segments, transforms)
    sam = SymbolicAWEModel(set, sys)
    init!(sam; prn=false)
    sim_steps(params) = round(Int, params["duration"] * set.sample_freq)
    mass = JobParameter("mass", "kg", 0.5, 10.0, 1.0; step=0.5)
    duration = JobParameter("duration", "s", 0.05, 5.0, 0.5)
    return JobKite("hanging_mass"; steady_menu=[mass], sim_menu=[mass, duration],
                   sweep_menu=[mass],
                   steady=params -> fly!(sam, kite_runs, params, 1),
                   simulate=(params; progress) ->
                       fly!(sam, kite_runs, params, sim_steps(params); progress),
                   sweep_point=params -> fly!(sam, kite_runs, params, 1))
end

"""Report progress until the job is cancelled."""
function run_until_cancelled(params; progress)
    for time in 0:0.1:600
        progress(time)
        sleep(0.1)
    end
    error("never cancelled")
end

"""A kite whose sim runs until cancelled and whose steady job throws."""
stub_kite() = JobKite("stub"; steady_menu=JobParameter[], sim_menu=JobParameter[],
                      sweep_menu=JobParameter[],
                      steady=params -> error("no steady state"),
                      simulate=run_until_cancelled,
                      sweep_point=params -> error("no sweep"))

@testset verbose = true "Job runner" begin
    previous_data_path = get_data_path()
    set_data_path(joinpath(dirname(@__DIR__), "data"))
    set = Settings("base/system.yaml")
    set.v_wind = 0
    kite_runs = Ref(0)
    kites = [hanging_mass_kite(set, kite_runs), stub_kite()]
    runs_dir = mktempdir()

    @testset "the job key hashes the canonical JSON the spec's test vector gives" begin
        version = Dict("SymbolicAWEModels" => Dict("version" => "0.19.0",
                                                   "commit" => nothing))
        key = job_key("open", "toy", "steady", Dict("mass" => 2.0, "drop" => 0.5),
                      version)
        @test occursin("with key `$key`", read(joinpath(dirname(@__DIR__), "api",
                                                        "openapi.yaml"), String))
        @test SymbolicAWEModels.canonical_number.([1e21, 1e-7, 1e-6, 100.0, 1 / 3]) ==
              ["1e+21", "1e-7", "0.000001", "100", "0.3333333333333333"]
    end

    @testset "each job kind's schema accepts a request of that kind" begin
        jobs_dir = joinpath(dirname(@__DIR__), "api", "jobs")
        schema(kind) = Schema(JSON.parsefile(joinpath(jobs_dir, "$kind.json"));
                              parent_dir=jobs_dir)
        @test isnothing(JSONSchema.validate(schema("steady"),
                        Dict("kind" => "steady", "parameters" => Dict("mass" => 2))))
        @test isnothing(JSONSchema.validate(schema("sweep"),
                        Dict("kind" => "sweep", "axes" => Dict("mass" => [1, 2]))))
        @test !isnothing(JSONSchema.validate(schema("sweep"), Dict("kind" => "sweep")))
        @test isnothing(JSONSchema.validate(schema("upload"), Dict("kind" => "upload")))
    end

    site = FakeSite()

    @testset "an upload job has a status, and no job without a key is claimed" begin
        upload = Dict("id" => "7", "sandbox" => "open", "kite" => "toy", "kind" => "upload",
                      "key" => nothing, "state" => "running",
                      "requested_at" => "2026-09-30T08:00:00Z")
        @test isnothing(JSONSchema.validate(site.schemas["JobStatus"], upload))
        @test !isnothing(JSONSchema.validate(site.schemas["Job"], upload))
    end

    @testset "a checkout with uncommitted changes versions as its commit plus -dirty" begin
        repo = mktempdir()
        git(args...) = run(`git -C $repo -c user.name=test -c user.email=test@test $args`)
        git("init", "-q")
        write(joinpath(repo, "file"), "committed")
        git("add", "file")
        git("commit", "-qm", "first")
        commit = SymbolicAWEModels.git_commit(repo)
        @test occursin(r"^[0-9a-f]{40}$", commit)
        write(joinpath(repo, "file"), "changed")
        @test SymbolicAWEModels.git_commit(repo) == commit * "-dirty"
    end

    port, listener = listenany(ip"127.0.0.1", 49152)
    Threads.@spawn serve_site(site, listener)
    runner = SymbolicAWEModels.JobRunner(kites; site="http://127.0.0.1:$port",
                                         token="runner-token")
    queue!(site, runner, "steady", "hanging_mass", "steady";
           parameters=Dict("mass" => 2.0))
    queue!(site, runner, "sim", "hanging_mass", "sim";
           parameters=Dict("mass" => 3.0, "duration" => 0.25))
    queue!(site, runner, "sweep", "hanging_mass", "sweep";
           axes=Dict("mass" => [1.0, 4.0]))
    queue!(site, runner, "too_heavy", "hanging_mass", "steady";
           parameters=Dict("mass" => 50.0))
    queue!(site, runner, "stale_key", "hanging_mass", "steady";
           parameters=Dict("mass" => 2.0), key=repeat("0", 64))
    queue!(site, runner, "broken", "stub", "steady")
    serving = Threads.@spawn serve_jobs(runner)

    @testset "the runner heartbeats its kites' menus and version" begin
        wait_until(() -> !isempty(site.heartbeats))
        heartbeat = site.heartbeats[1]
        kite = only(filter(kite -> kite["name"] == "hanging_mass", heartbeat["kites"]))
        @test kite["menu"]["sim"]["parameters"][2]["name"] == "duration"
        @test kite["menu"]["sweep"]["max_axes"] == 1
        @test haskey(heartbeat["version"], "SymbolicAWEModels")
    end

    @testset "a steady job uploads one frame whose run metadata round-trips" begin
        wait_until(() -> haskey(site.results, "steady"))
        request = site.results["steady"]
        @test request.headers["content-type"] == "application/vnd.apache.arrow.file"
        log = read_run(runs_dir, "steady", request.body)
        run = JSON.parse(log.metadata["run"])
        @test isnothing(JSONSchema.validate(site.schemas["RunMetadata"], run))
        @test run["kind"] == "steady" && run["kite"] == "hanging_mass"
        @test run["parameters"] == Dict("mass" => 2.0)
        @test run["job_key"] == job_key("open", "hanging_mass", "steady",
                                        Dict("mass" => 2.0), runner.version)
        @test haskey(JSON.parse(log.metadata["topology"]), "metadata")
        @test length(log.syslog) == 1
    end

    @testset "a sim job reports progress and uploads its whole run" begin
        wait_until(() -> haskey(site.results, "sim"))
        log = read_run(runs_dir, "sim", site.results["sim"].body)
        @test length(log.syslog) == round(Int, 0.25 * set.sample_freq)
        @test JSON.parse(log.metadata["run"])["parameters"] ==
              Dict("mass" => 3.0, "duration" => 0.25)
        @test any(((id, report),) -> id == "sim" && haskey(report, "sim_time"),
                  site.progress)
    end

    @testset "a sweep job uploads a manifest and a run per point" begin
        wait_until(() -> haskey(site.results, "sweep"))
        parts = multipart_parts(site.results["sweep"])
        manifest = JSON.parse(String(parts["manifest.json"]))
        @test isnothing(JSONSchema.validate(site.schemas["SweepManifest"], manifest))
        @test [point["parameters"]["mass"] for point in manifest["points"]] == [1.0, 4.0]
        for point in manifest["points"]
            log = read_run(runs_dir, splitext(point["run"])[1], parts[point["run"]])
            @test JSON.parse(log.metadata["run"])["parameters"] == point["parameters"]
        end
    end

    @testset "an out-of-range parameter is refused before the kite runs" begin
        wait_until(() -> haskey(site.failures, "too_heavy"))
        @test site.failures["too_heavy"]["reason"] == "rejected"
        @test occursin("mass = 50", site.failures["too_heavy"]["log"])
        wait_until(() -> haskey(site.failures, "stale_key"))
        @test kite_runs[] == 4  # the steady job, the sim and two sweep points
        @test !haskey(site.results, "too_heavy")
    end

    @testset "a job whose key is not the runner's is refused before the kite runs" begin
        wait_until(() -> haskey(site.failures, "stale_key"))
        @test site.failures["stale_key"]["reason"] == "rejected"
        @test occursin("is not this runner's", site.failures["stale_key"]["log"])
        @test !haskey(site.results, "stale_key")
    end

    @testset "a job that throws is reported, and the runner carries on" begin
        wait_until(() -> haskey(site.failures, "broken"))
        @test site.failures["broken"]["reason"] == "error"
        @test occursin("no steady state", site.failures["broken"]["log"])
    end

    @testset "a cancelled job stops running" begin
        queue!(site, runner, "endless", "stub", "sim")
        wait_until(() -> any(((id, _),) -> id == "endless", site.progress))
        @lock site.lock push!(site.cancelled, "endless")
        wait_until(() -> @lock(runner.lock, isempty(runner.cancelled)); timeout=30)
        @test !haskey(site.results, "endless")
        @test !haskey(site.failures, "endless")
    end

    @testset "a 200 without `cancel` is a failed call, not a crash" begin
        site.bare_replies[] = true
        @test_logs (:warn, r"Heartbeat failed") SymbolicAWEModels.heartbeat!(runner)
        progress = SymbolicAWEModels.JobProgress(runner, "endless",
                                                 Threads.Atomic{Bool}(false), time(), -Inf)
        @test_logs (:warn, r"Reporting progress") progress(1.0)
        site.bare_replies[] = false
    end

    runner.stopped[] = true
    wait(serving)
    close(listener)
    set_data_path(previous_data_path)
end
