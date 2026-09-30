# Copyright (c) 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Fly the V3 kite under each aerodynamic coupling and log the runs the paper
# compares. Needs the NeuralFoil aero geometry: run make_aero_geometry.jl first.
# Writes output/logs/<case>.arrow and a row of output/cases.csv as each case ends;
# `RUNS` keeps each case's model in memory for the section figures of make_figures.jl.
# Case names as arguments fly only those.

using Pkg
Pkg.activate(@__DIR__)

using V3Kite
using SymbolicAWEModels
using Logging
using Printf

include(joinpath(@__DIR__, "paper_setup.jl"))

"""
    FailureCounter(parent)

Logger that forwards to `parent` and counts the VSM solves reported as failed.
"""
mutable struct FailureCounter <: AbstractLogger
    parent::AbstractLogger
    failed_solves::Int
end
FailureCounter(parent) = FailureCounter(parent, 0)
Logging.min_enabled_level(logger::FailureCounter) = Logging.min_enabled_level(logger.parent)
Logging.shouldlog(logger::FailureCounter, args...) = true
Logging.catch_exceptions(logger::FailureCounter) = Logging.catch_exceptions(logger.parent)
function Logging.handle_message(logger::FailureCounter, level, message, mod, group,
                                id, file, line; kwargs...)
    occursin("did not converge", string(message)) && (logger.failed_solves += 1)
    return Logging.handle_message(logger.parent, level, message, mod, group, id,
                                  file, line; kwargs...)
end

"""
    fly(case; sim_time) -> NamedTuple

Build `case.project` under `case.mode`, fly it for `sim_time` [s] with the steering
of `steering_offset`, and log it. Returns the run's summary, with integration and
VSM time [s] summed from the second step on, and its model.
"""
function fly(case; sim_time)
    kite_set = load_kite(case.project; data_path=DATA_PATH)
    kite_set.aero_mode = case.mode
    sam, sys = build_v3_model(case.project; data_path=DATA_PATH, kite_set)
    nominal = V3Kite.get_steering(sys, kite_set.geom)
    n_steps = round(Int, sim_time / case.dt)
    logger, sys_state = create_logger(sam, n_steps)
    counter = FailureCounter(current_logger())
    t_step = t_vsm = 0.0
    steps = 0
    started = time()
    with_logger(counter) do
        for step in 1:n_steps
            t = step * case.dt
            steering = nominal + steering_offset(t)
            set_steering!(sys, steering, kite_set.geom)
            sam.integrator.opts.maxiters = sam.integrator.iter + MAX_SOLVER_STEPS
            flew = try
                sim_step!(sam; dt=case.dt, vsm_interval=case.vsm_interval,
                          vsm_warn_on_fail=true)
            catch failure
                failure isa InterruptException && rethrow()
                @warn "Case stopped" case.name t exception=failure
                false
            end
            flew || break
            all(isfinite, sys.wings[1].pos_w) || break
            if step > 1
                t_step += sam.t_step
                t_vsm += sam.t_vsm
            end
            steps = step
            log_state!(logger, sys_state, sam, t; steering)
            counter.failed_solves > MAX_FAILED_SOLVES && break
            time() - started > WALL_LIMIT && break
            step % round(Int, 1 / case.dt) == 0 &&
                @info "Flown" case.name t wall=round(time() - started; digits=1)
        end
    end
    save_log(logger, case.name; path=LOG_PATH)
    flown = steps * case.dt
    return (; case.name, mode=mode_label(case.mode), case.dt, case.vsm_interval,
            flown, failed_solves=counter.failed_solves, t_step, t_vsm, sam)
end

"""
    record_case(run)

Replace case `run.name`'s row of output/cases.csv with `run`'s summary.
"""
function record_case(run)
    path = joinpath(OUTPUT_PATH, "cases.csv")
    header = join(CASE_COLUMNS, ",")
    rows = isfile(path) ? readlines(path)[2:end] : String[]
    filter!(row -> !startswith(row, run.name * ","), rows)
    push!(rows, @sprintf("%s,%s,%.3f,%d,%.2f,%d,%.2f,%.2f", run.name, run.mode, run.dt,
                         run.vsm_interval, run.flown, run.failed_solves, run.t_step,
                         run.t_vsm))
    write(path, join([header; rows], "\n") * "\n")
    return nothing
end

"""
    selected(cases)

The `cases` named in the script's arguments, or all of them when it has none.
"""
selected(cases) = filter(case -> isempty(ARGS) || case.name in ARGS, cases)

mkpath(LOG_PATH)
@isdefined(RUNS) || (RUNS = Dict{String, Any}())
for (cases, sim_time) in ((CASES, SIM_TIME), (SECTION_CASES, SECTION_TIME))
    for case in selected(cases)
        @info "Flying" case.name
        RUNS[case.name] = fly(case; sim_time)
        record_case(RUNS[case.name])
    end
end
nothing
