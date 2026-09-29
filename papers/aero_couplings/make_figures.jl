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
STATION = 5             # mid-span structural station of the beam V3
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
    log = flight_log(first(COMPARED))
    lines!(axes[1], log.time, rad2deg.(log.bearing); color=:gray, linestyle=:dash,
           label="setpoint")
    Legend(fig[0, 1], axes[1]; orientation=:horizontal, framevisible=false,
           tellwidth=false, nbanks=1)
    return fig
end

"""
    station_loads_figure() -> Figure

The aero point loads of structural station `STATION` at the end of each compared case.
"""
function station_loads_figure()
    fig = Figure(size=(TEXT_WIDTH, 260))
    for (col, name) in enumerate(COMPARED)
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
    flight_statistics(name) -> NamedTuple

Mean tether force [N], mean and high-pass RMS angle of attack [deg], RMS heading
error [deg] and RMS high-pass wing-node speed [m/s] of case `name` after
`SETTLE_TIME`.
"""
function flight_statistics(name)
    log = flight_log(name)
    rows = settled(log)
    any(rows) || return nothing
    dt = RUNS[name].dt
    aoa = rad2deg.(log.AoA)
    heading_error = rad2deg.(log.heading .- log.bearing)
    wing = [i for i in eachindex(first(log.X)) if i <= length(RUNS[name].sam.sys_struct.points) &&
            RUNS[name].sam.sys_struct.points[i].wing_idx == 1]
    node_speed = [mean(norm((log.VX[k][i], log.VY[k][i], log.VZ[k][i])) for i in wing)
                  for k in eachindex(log.time)]
    return (force=mean(first.(log.winch_force)[rows]), aoa=mean(aoa[rows]),
            aoa_ripple=sqrt(mean(abs2, high_pass(aoa, dt)[rows])),
            heading_error=sqrt(mean(abs2, heading_error[rows])),
            node_ripple=sqrt(mean(abs2, high_pass(node_speed, dt)[rows])))
end

mkpath(FIGURE_PATH)
save(joinpath(FIGURE_PATH, "tracking.pdf"), tracking_figure())
save(joinpath(FIGURE_PATH, "station_loads.pdf"), station_loads_figure())
nothing
