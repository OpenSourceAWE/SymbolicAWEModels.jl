# Copyright (c) 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Write the paper's figures to figures/ and the numbers it quotes to results.tex.
# Runs in the session run_simulations.jl flew in: the section figures read the
# models `RUNS` holds at the end of each case.

using Pkg
Pkg.activate(@__DIR__)

using CairoMakie
using MakieControlPlots
using SymbolicAWEModels
using KiteUtils
using LinearAlgebra
using DelimitedFiles
using Printf
using Statistics

include(joinpath(@__DIR__, "paper_setup.jl"))
@isdefined(RUNS) || error("Run run_simulations.jl in this session first")

CairoMakie.activate!(type="pdf")
set_theme!(theme_latexfonts(); fontsize=10, linewidth=1.2)
COLUMN_WIDTH = 252      # [pt] one column of the paper
TEXT_WIDTH = 453        # [pt] the paper's text width
MODE_COLORS = Dict("AeroDirect" => :black, "ContinuousAero" => :dodgerblue3,
                   "AeroPressure" => :darkorange, "AeroPressure+live" => :forestgreen)
COMPARED = ["direct_dt10", "continuous_dt10", "live_dt10", "pressure_dt10"]
SECTIONS = ["continuous_particle", "pressure_beam", "live_beam"]
STATION = 5             # a mid-span structural station of both V3 models
PANEL = 20              # a mid-span VSM panel
SETTLE_TIME = 5.0       # [s] start of the window the flight statistics are taken over

"""
    flight_log(name) -> SysLog

The log run_simulations.jl wrote for case `name`.
"""
flight_log(name) = load_log(name; path=LOG_PATH).syslog

"""
    settled(log) -> BitVector

Rows of `log` after `SETTLE_TIME`.
"""
settled(log) = log.time .>= SETTLE_TIME

"""
    high_pass(signal, dt; window=0.5) -> Vector

`signal` minus its centred moving mean over `window` [s]: the content above roughly
`1/window` Hz.
"""
function high_pass(signal, dt; window=0.5)
    half = max(1, round(Int, window / 2dt))
    n = length(signal)
    return [signal[i] - mean(signal[max(1, i - half):min(n, i + half)]) for i in 1:n]
end

"""
    unwrapped(angle) -> Vector

`angle` [rad] with its jumps of a full turn removed.
"""
function unwrapped(angle)
    turns = cumsum([0; round.(diff(angle) ./ 2pi)])
    return angle .- 2pi .* turns
end

"""
    tracking_figure() -> Figure

Heading, tether force and angle of attack of the compared cases against time, with
the steering offset they fly.
"""
function tracking_figure()
    fig = Figure(size=(TEXT_WIDTH, 380))
    labels = [L"u_s", L"\psi~[°]", L"F_t~[\mathrm{kN}]", L"\alpha~[°]"]
    axes = [Axis(fig[row, 1]; ylabel=labels[row], xticklabelsvisible=row == 4)
            for row in 1:4]
    axes[4].xlabel = L"t~[\mathrm{s}]"
    linkxaxes!(axes...)
    time = 0:0.05:SIM_TIME
    lines!(axes[1], time, steering_offset.(time); color=:gray)
    for name in COMPARED
        run = RUNS[name]
        log = flight_log(name)
        color = MODE_COLORS[run.mode]
        lines!(axes[2], log.time, rad2deg.(unwrapped(log.heading)); color, label=run.mode)
        lines!(axes[3], log.time, first.(log.winch_force) ./ 1e3; color)
        lines!(axes[4], log.time, rad2deg.(log.AoA); color)
    end
    Legend(fig[0, 1], axes[2]; orientation=:horizontal, framevisible=false,
           tellwidth=false, nbanks=1)
    return fig
end

"""
    trajectory_figure() -> Figure

Elevation against azimuth of the compared cases: the loop each coupling flies.
"""
function trajectory_figure()
    fig = Figure(size=(COLUMN_WIDTH, 220))
    ax = Axis(fig[1, 1]; xlabel=L"\phi~[°]", ylabel=L"\beta~[°]", aspect=DataAspect())
    for name in COMPARED
        run = RUNS[name]
        log = flight_log(name)
        lines!(ax, rad2deg.(log.azimuth), rad2deg.(log.elevation);
               color=MODE_COLORS[run.mode], label=run.mode)
    end
    return fig
end

"""
    station_loads_figure() -> Figure

The aero point loads of structural station `STATION` at the end of each of
`SECTIONS`.
"""
function station_loads_figure()
    fig = Figure(size=(TEXT_WIDTH, 260))
    for (col, name) in enumerate(SECTIONS)
        run = RUNS[name]
        ax = Axis(fig[1, col]; title=run.mode, aspect=DataAspect(),
                  xlabel=L"x_b~[\mathrm{m}]", ylabel=L"z_b~[\mathrm{m}]",
                  ylabelvisible=col == 1, yticklabelsvisible=col == 1)
        plot_station_loads!(ax, run.sam.sys_struct, STATION; force_scale=1.5e-3,
                            color=MODE_COLORS[run.mode])
        col > 1 && linkaxes!(ax, content(fig[1, 1]))
    end
    return fig
end

"""
    polar_figure(panel_idx) -> Figure

Lift and moment polars of VSM panel `panel_idx` of the beam wing at the end of the
tabulated and the live-polar `AeroPressure` section cases.
"""
function polar_figure(panel_idx)
    fig = Figure(size=(TEXT_WIDTH, 200))
    for (col, coefficient) in enumerate((:cl, :cm))
        ax = Axis(fig[1, col]; xlabel=L"\alpha~[°]",
                  ylabel=coefficient === :cl ? L"C_l" : L"C_m")
        for name in ("pressure_beam", "live_beam")
            run = RUNS[name]
            plot_panel_polar!(ax, run.sam.sys_struct, panel_idx; coefficient,
                              color=MODE_COLORS[run.mode], label=run.mode)
        end
        col == 1 && axislegend(ax; position=:rb, framevisible=false)
    end
    return fig
end

"""
    pressure_figure(panel_idx) -> Figure

The `-Cp` pattern of VSM panel `panel_idx` of the beam wing at the end of the
tabulated and the live-polar `AeroPressure` section cases.
"""
function pressure_figure(panel_idx)
    fig = Figure(size=(COLUMN_WIDTH, 200))
    ax = Axis(fig[1, 1]; xlabel=L"x/c", ylabel=L"-C_p")
    for name in ("pressure_beam", "live_beam")
        run = RUNS[name]
        plot_panel_pressure!(ax, run.sam.sys_struct, panel_idx;
                             color=MODE_COLORS[run.mode], label=run.mode)
    end
    axislegend(ax; position=:rt, framevisible=false)
    return fig
end

"""
    case_table() -> Vector{NamedTuple}

The rows output/cases.csv holds, in the order of `CASES`.
"""
function case_table()
    rows = readdlm(joinpath(OUTPUT_PATH, "cases.csv"), ','; header=true)[1]
    table = [(name=String(row[1]), mode=String(row[2]), dt=Float64(row[3]),
              vsm_interval=Int(row[4]), flown=Float64(row[5]), completed=row[6] == "true",
              failed_solves=Int(row[7]), wall=Float64(row[8]), t_step=Float64(row[9]),
              t_vsm=Float64(row[10])) for row in eachrow(rows)]
    order = [case.name for case in CASES]
    return sort!(table; by=row -> findfirst(==(row.name), order))
end

"""
    wing_points(sys) -> Vector{Int}

Indices of the points of `sys` that belong to its first wing.
"""
wing_points(sys) = [point.idx for point in sys.points if point.wing_idx == 1]

"""
    flight_statistics(row) -> NamedTuple

Mean tether force [kN], high-pass RMS angle of attack [deg] and high-pass RMS mean
wing-node speed [m/s] of case `row` after `SETTLE_TIME`, or `nothing` when it did not
fly that long.
"""
function flight_statistics(row)
    log = flight_log(row.name)
    rows = settled(log)
    any(rows) || return nothing
    wing = wing_points(RUNS[row.name].sam.sys_struct)
    aoa = rad2deg.(log.AoA)
    node_speed = [mean(norm((log.VX[k][i], log.VY[k][i], log.VZ[k][i])) for i in wing)
                  for k in eachindex(log.time)]
    return (force=mean(first.(log.winch_force)[rows]) / 1e3,
            aoa_ripple=sqrt(mean(abs2, high_pass(aoa, row.dt)[rows])),
            node_ripple=sqrt(mean(abs2, high_pass(node_speed, row.dt)[rows])))
end

"""
    model_label(name) -> String

The structural model case `name` flies: `particle` or `beam`.
"""
function model_label(name)
    project = only(case.project for case in CASES if case.name == name)
    return project == BEAM ? "beam" : "particle"
end

"""
    case_row(row, stats) -> String

One line of the paper's case table: `row` of `case_table` and its `flight_statistics`.
"""
function case_row(row, stats)
    timed = row.flown - row.dt
    cost = timed > 0 ? @sprintf("%.1f & %.0f", (row.t_step + row.t_vsm) / timed,
                                100 * row.t_vsm / (row.t_step + row.t_vsm)) : "-- & --"
    flight = isnothing(stats) ? "-- & -- & --" :
             @sprintf("%.2f & %.2f & %.3f", stats.force, stats.aoa_ripple,
                      stats.node_ripple)
    return @sprintf("\\code{%s} & %s & %.0f & %d & %.1f & %d & %s & %s \\\\",
                    row.mode, model_label(row.name), 1e3 * row.dt, row.vsm_interval,
                    row.flown, row.failed_solves, cost, flight)
end

"""
    aero_size_macros(io, suffix, sys)

Write `\\NumSections<suffix>` and `\\NumPanels<suffix>`, the sections and VSM
panels of `sys`'s wing.
"""
function aero_size_macros(io, suffix, sys)
    wing = sys.wings[1]
    println(io, "\\newcommand{\\NumSections", suffix, "}{",
            length(wing.vsm_wing.unrefined_sections), "}")
    println(io, "\\newcommand{\\NumPanels", suffix, "}{", length(wing.vsm_aero.panels), "}")
    return nothing
end

"""
    write_results(path, table)

Write the macros the paper reads: the manoeuvre, the size of both wings' aero models,
and `\\CaseRows`, one table row per case of `table`. Returns each row's
`flight_statistics`.
"""
function write_results(path, table)
    stats = flight_statistics.(table)
    open(path, "w") do io
        println(io, "% Written by make_figures.jl")
        println(io, "\\newcommand{\\Steering}{", STEERING, "}")
        println(io, "\\newcommand{\\RampStart}{", RAMP[1], "}")
        println(io, "\\newcommand{\\RampEnd}{", RAMP[2], "}")
        println(io, "\\newcommand{\\SimTime}{", SIM_TIME, "}")
        println(io, "\\newcommand{\\SectionTime}{", SECTION_TIME, "}")
        println(io, "\\newcommand{\\SettleTime}{", SETTLE_TIME, "}")
        aero_size_macros(io, "Particle", RUNS["continuous_particle"].sam.sys_struct)
        aero_size_macros(io, "Beam", RUNS["pressure_beam"].sam.sys_struct)
        println(io, "\\newcommand{\\CaseRows}{")
        foreach((row, stat) -> println(io, case_row(row, stat)), table, stats)
        println(io, "}")
    end
    return stats
end

mkpath(FIGURE_PATH)
TABLE = case_table()
STATS = write_results(joinpath(@__DIR__, "results.tex"), TABLE)
CairoMakie.save(joinpath(FIGURE_PATH, "tracking.pdf"), tracking_figure())
CairoMakie.save(joinpath(FIGURE_PATH, "trajectory.pdf"), trajectory_figure())
CairoMakie.save(joinpath(FIGURE_PATH, "station_loads.pdf"), station_loads_figure())
CairoMakie.save(joinpath(FIGURE_PATH, "polar.pdf"), polar_figure(PANEL))
CairoMakie.save(joinpath(FIGURE_PATH, "pressure.pdf"), pressure_figure(PANEL))
nothing
