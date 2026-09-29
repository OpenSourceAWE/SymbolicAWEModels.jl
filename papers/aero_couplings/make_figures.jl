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
COMPARED = ["direct_dt10", "continuous_dt10", "pressure_dt10", "live_dt10"]
SECTION_CASES = ["continuous_dt10", "pressure_dt10", "live_dt10"]
STATION = 5             # mid-span structural station of the beam V3
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
    tracking_figure() -> Figure

Heading, tether force and angle of attack of the compared cases against time.
"""
function tracking_figure()
    fig = Figure(size=(TEXT_WIDTH, 330))
    labels = [L"\psi~[°]", L"F_t~[\mathrm{kN}]", L"\alpha~[°]"]
    axes = [Axis(fig[row, 1]; ylabel=labels[row], xticklabelsvisible=row == 3)
            for row in 1:3]
    axes[3].xlabel = L"t~[\mathrm{s}]"
    linkxaxes!(axes...)
    for name in COMPARED
        run = RUNS[name]
        log = flight_log(name)
        color = MODE_COLORS[run.mode]
        lines!(axes[1], log.time, rad2deg.(log.heading); color, label=run.mode)
        lines!(axes[2], log.time, first.(log.winch_force) ./ 1e3; color)
        lines!(axes[3], log.time, rad2deg.(log.AoA); color)
    end
    time = 0:0.05:SIM_TIME
    lines!(axes[1], time, rad2deg.(heading_setpoint.(time)); color=:gray,
           linestyle=:dash, label="setpoint")
    Legend(fig[0, 1], axes[1]; orientation=:horizontal, framevisible=false,
           tellwidth=false, nbanks=1)
    return fig
end

"""
    station_loads_figure() -> Figure

The aero point loads of structural station `STATION` at the end of each of
`SECTION_CASES`.
"""
function station_loads_figure()
    fig = Figure(size=(TEXT_WIDTH, 260))
    for (col, name) in enumerate(SECTION_CASES)
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

Lift and moment polars of VSM panel `panel_idx` at the end of the tabulated and the
live-polar `AeroPressure` cases.
"""
function polar_figure(panel_idx)
    fig = Figure(size=(TEXT_WIDTH, 200))
    for (col, coefficient) in enumerate((:cl, :cm))
        ax = Axis(fig[1, col]; xlabel=L"\alpha~[°]",
                  ylabel=coefficient === :cl ? L"C_l" : L"C_m")
        for name in ("pressure_dt10", "live_dt10")
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

The `-Cp` pattern of VSM panel `panel_idx` at the end of the tabulated and the
live-polar `AeroPressure` cases.
"""
function pressure_figure(panel_idx)
    fig = Figure(size=(COLUMN_WIDTH, 200))
    ax = Axis(fig[1, 1]; xlabel=L"x/c", ylabel=L"-C_p")
    for name in ("pressure_dt10", "live_dt10")
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
    flight_statistics(row, wing) -> NamedTuple

Mean tether force [kN], high-pass RMS angle of attack [deg], RMS heading error [deg]
and high-pass RMS mean speed of the `wing` points [m/s] of case `row` after
`SETTLE_TIME`, or `nothing` when it did not fly that long.
"""
function flight_statistics(row, wing)
    log = flight_log(row.name)
    rows = settled(log)
    any(rows) || return nothing
    aoa = rad2deg.(log.AoA)
    heading_error = rad2deg.(log.heading .- log.bearing)
    node_speed = [mean(norm((log.VX[k][i], log.VY[k][i], log.VZ[k][i])) for i in wing)
                  for k in eachindex(log.time)]
    return (force=mean(first.(log.winch_force)[rows]) / 1e3,
            aoa_ripple=sqrt(mean(abs2, high_pass(aoa, row.dt)[rows])),
            heading_error=sqrt(mean(abs2, heading_error[rows])),
            node_ripple=sqrt(mean(abs2, high_pass(node_speed, row.dt)[rows])))
end

"""
    case_row(row, stats) -> String

One line of the paper's case table: `row` of `case_table` and its `flight_statistics`.
"""
function case_row(row, stats)
    timed = row.flown - row.dt
    cost = timed > 0 ? @sprintf("%.1f & %.0f", (row.t_step + row.t_vsm) / timed,
                                100 * row.t_vsm / (row.t_step + row.t_vsm)) : "-- & --"
    flight = isnothing(stats) ? "-- & -- & -- & --" :
             @sprintf("%.2f & %.2f & %.1f & %.3f", stats.force, stats.aoa_ripple,
                      stats.heading_error, stats.node_ripple)
    return @sprintf("\\code{%s} & %.0f & %d & %.1f & %d & %s & %s \\\\",
                    row.mode, 1e3 * row.dt, row.vsm_interval, row.flown,
                    row.failed_solves, cost, flight)
end

"""
    write_results(path, table, sys)

Write the macros the paper reads: the flight condition, the size of `sys`'s aero
model, and `\\CaseRows`, one table row per case of `table`.
"""
function write_results(path, table, sys)
    wing = sys.wings[1]
    stats = [flight_statistics(row, wing_points(sys)) for row in table]
    open(path, "w") do io
        println(io, "% Written by make_figures.jl")
        println(io, "\\newcommand{\\MaxHeading}{", MAX_HEADING, "}")
        println(io, "\\newcommand{\\Period}{", PERIOD, "}")
        println(io, "\\newcommand{\\SimTime}{", SIM_TIME, "}")
        println(io, "\\newcommand{\\SettleTime}{", SETTLE_TIME, "}")
        println(io, "\\newcommand{\\NumSections}{",
                length(wing.vsm_wing.unrefined_sections), "}")
        println(io, "\\newcommand{\\NumPanels}{", length(wing.vsm_aero.panels), "}")
        println(io, "\\newcommand{\\CaseRows}{")
        foreach((row, stat) -> println(io, case_row(row, stat)), table, stats)
        println(io, "}")
    end
    return stats
end

mkpath(FIGURE_PATH)
TABLE = case_table()
STATS = write_results(joinpath(@__DIR__, "results.tex"), TABLE,
                      RUNS[first(SECTION_CASES)].sam.sys_struct)
save(joinpath(FIGURE_PATH, "tracking.pdf"), tracking_figure())
save(joinpath(FIGURE_PATH, "station_loads.pdf"), station_loads_figure())
save(joinpath(FIGURE_PATH, "polar.pdf"), polar_figure(PANEL))
save(joinpath(FIGURE_PATH, "pressure.pdf"), pressure_figure(PANEL))
nothing
