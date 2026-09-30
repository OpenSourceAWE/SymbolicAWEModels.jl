# Copyright (c) 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Write the paper's figures to figures/ and the numbers it quotes to results.tex.
# The table reads the logs and output/cases.csv; the section figures and numbers read
# the models `RUNS` holds at the end of the section cases, so run
# `run_simulations.jl continuous_particle pressure_beam live_beam` in this session first.

using Pkg
Pkg.activate(@__DIR__)

using CairoMakie
import MakieControlPlots   # with GeometryBasics, loads the package's Makie extension
using V3Kite
using SymbolicAWEModels
using VortexStepMethod
using LinearAlgebra
using DelimitedFiles
using Printf
using Statistics

include(joinpath(@__DIR__, "paper_setup.jl"))
@isdefined(RUNS) || error("Fly the section cases with run_simulations.jl first")

CairoMakie.activate!(type="pdf")
set_theme!(theme_latexfonts(); fontsize=10, linewidth=1.2)
COLUMN_WIDTH = 252      # [pt] one column of the paper
TEXT_WIDTH = 453        # [pt] the paper's text width
MODE_COLORS = Dict("AeroDirect" => :black, "ContinuousAero" => :dodgerblue3,
                   "AeroPressure" => :darkorange, "AeroPressure+live" => :forestgreen)
COMPARED = ["direct_dt10", "continuous_dt10", "live_dt10", "pressure_dt10"]
SECTIONS = ["continuous_particle", "pressure_beam", "live_beam"]
STATION = 5             # a mid-span structural station of both V3 models
PANEL = 20              # a VSM panel a quarter of the span from the tip
SETTLE_TIME = 5.0       # [s] start of the window the flight statistics are taken over

"""
    find_case(name)

The case of `CASES` or `SECTION_CASES` named `name`.
"""
find_case(name) = only(case for case in [CASES; SECTION_CASES] if case.name == name)

"""
    flight_log(name) -> SysLog

The log run_simulations.jl wrote for case `name`.
"""
flight_log(name) = load_log(name; path=LOG_PATH).syslog

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

Steering offset, tether force and angle of attack of the compared cases against time.
"""
function tracking_figure()
    fig = Figure(size=(TEXT_WIDTH, 300))
    labels = [L"u_s", L"F_t~[\mathrm{kN}]", L"\alpha~[°]"]
    axes = [Axis(fig[row, 1]; ylabel=labels[row], xticklabelsvisible=row == 3)
            for row in 1:3]
    axes[3].xlabel = L"t~[\mathrm{s}]"
    linkxaxes!(axes...)
    time = 0:0.05:SIM_TIME
    lines!(axes[1], time, steering_offset.(time); color=:gray)
    for name in COMPARED
        mode = mode_label(find_case(name).mode)
        log = flight_log(name)
        color = MODE_COLORS[mode]
        lines!(axes[2], log.time, first.(log.winch_force) ./ 1e3; color, label=mode)
        lines!(axes[3], log.time, rad2deg.(log.AoA); color)
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
    axis = Axis(fig[1, 1]; xlabel=L"\phi~[°]", ylabel=L"\beta~[°]", aspect=DataAspect())
    for name in COMPARED
        log = flight_log(name)
        lines!(axis, rad2deg.(log.azimuth), rad2deg.(log.elevation);
               color=MODE_COLORS[mode_label(find_case(name).mode)])
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
        axis = Axis(fig[1, col]; title=run.mode, aspect=DataAspect(),
                  xlabel=L"x_b~[\mathrm{m}]", ylabel=L"z_b~[\mathrm{m}]",
                  ylabelvisible=col == 1, yticklabelsvisible=col == 1)
        plot_station_loads!(axis, run.sam.sys_struct, STATION; force_scale=1.5e-3,
                            color=MODE_COLORS[run.mode])
        col > 1 && linkaxes!(axis, content(fig[1, 1]))
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
        axis = Axis(fig[1, col]; xlabel=L"\alpha~[°]",
                  ylabel=coefficient === :cl ? L"C_l" : L"C_m")
        for name in ("pressure_beam", "live_beam")
            run = RUNS[name]
            plot_panel_polar!(axis, run.sam.sys_struct, panel_idx; coefficient,
                              color=MODE_COLORS[run.mode], label=run.mode)
        end
        col == 1 && axislegend(axis; position=:rb, framevisible=false)
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
    axis = Axis(fig[1, 1]; xlabel=L"x/c", ylabel=L"-C_p")
    for name in ("pressure_beam", "live_beam")
        run = RUNS[name]
        plot_panel_pressure!(axis, run.sam.sys_struct, panel_idx;
                             color=MODE_COLORS[run.mode], label=run.mode)
    end
    axislegend(axis; position=:rt, framevisible=false)
    return fig
end

"""
    case_rows() -> Dict{String, NamedTuple}

The rows output/cases.csv holds, by case name.
"""
function case_rows()
    rows = readdlm(joinpath(OUTPUT_PATH, "cases.csv"), ','; header=true)[1]
    table = [NamedTuple{CASE_COLUMNS}(Tuple(value isa AbstractString ? String(value) : value
                                            for value in row)) for row in eachrow(rows)]
    return Dict(row.name => row for row in table)
end

"""
    unbuilt_structure(case) -> SystemStructure

The structure `case` flies, loaded from its project without compiling a model.
"""
function unbuilt_structure(case)
    kite_set = load_kite(case.project; data_path=DATA_PATH)
    kite_set.aero_mode = case.mode
    return last(V3Kite.create_v3_model(case.project; data_path=DATA_PATH, kite_set))
end

"""
    flight_statistics(row, wing_nodes) -> NamedTuple

Mean tether force [kN], high-pass RMS angle of attack [deg] and high-pass RMS mean
speed [m/s] of the points `wing_nodes` of case `row` after `SETTLE_TIME`, or `nothing`
when it did not fly that long.
"""
function flight_statistics(row, wing_nodes)
    log = flight_log(row.name)
    rows = log.time .>= SETTLE_TIME
    any(rows) || return nothing
    aoa = rad2deg.(log.AoA)
    node_speed = [mean(norm((log.VX[k][i], log.VY[k][i], log.VZ[k][i])) for i in wing_nodes)
                  for k in eachindex(log.time)]
    return (force=mean(first.(log.winch_force)[rows]) / 1e3,
            aoa_ripple=sqrt(mean(abs2, high_pass(aoa, row.dt)[rows])),
            node_ripple=sqrt(mean(abs2, high_pass(node_speed, row.dt)[rows])))
end

"""
    case_row(row, model, stats) -> String

One line of the paper's case table: `row` of `case_rows`, flown on `model`, and its
`flight_statistics`.
"""
function case_row(row, model, stats)
    timed = row.flown - row.dt
    cost = timed > 0 ? @sprintf("%.1f & %.0f", (row.t_step + row.t_vsm) / timed,
                                100 * row.t_vsm / (row.t_step + row.t_vsm)) : "-- & --"
    flight = isnothing(stats) ? "-- & -- & --" :
             @sprintf("%.2f & %.2f & %.3f", stats.force, stats.aoa_ripple,
                      stats.node_ripple)
    return @sprintf("\\code{%s} & %s & %.0f & %d & %.1f & %d & %s & %s \\\\",
                    row.mode, model, 1e3 * row.dt, row.vsm_interval, row.flown,
                    row.failed_solves, cost, flight)
end

"""
    operating_coefficient(run, coefficient) -> Float64

Coefficient `coefficient` (`:cl` or `:cm`) of panel `PANEL` at its angle of attack at
the end of section case `run`.
"""
function operating_coefficient(run, coefficient)
    wing = run.sam.sys_struct.wings[1]
    calculate = getproperty(VortexStepMethod, Symbol(:calculate_, coefficient))
    return calculate(wing.vsm_aero.panels[PANEL], wing.vsm_solver.sol.alpha_dist[PANEL])
end

"""
    wing_node_idxs(sys) -> Vector{Int}

Indices of the points of the wing of `sys`.
"""
wing_node_idxs(sys) = [point.idx
                       for point in SymbolicAWEModels.wing_points(sys, sys.wings[1])]

"""
    ratio_text(numerator, denominator) -> String

`numerator / denominator` to one decimal.
"""
ratio_text(numerator, denominator) = @sprintf("%.1f", numerator / denominator)

"""
    tex_macro(io, name, value)

Write the LaTeX macro `\\<name>`, expanding to `value`.
"""
tex_macro(io, name, value) = println(io, "\\newcommand{\\", name, "}{", value, "}")

"""
    write_results(path, rows)

Write the macros the paper reads: the manoeuvre, the size of both wings' aero models,
`\\CaseRows` with one table row per case of `CASES` in `rows`, the ripple of a solve
every fifth step over every step, and what the live-polar section case found.
"""
function write_results(path, rows)
    table = [rows[case.name] for case in CASES if haskey(rows, case.name)]
    structures = Dict(case.project => unbuilt_structure(case)
                      for case in unique(case -> case.project, CASES))
    wing_nodes = Dict(project => wing_node_idxs(sys) for (project, sys) in structures)
    stats = Dict(row.name => flight_statistics(row, wing_nodes[find_case(row.name).project])
                 for row in table)
    stale, fresh = stats["continuous_dt10_vsm5"], stats["continuous_dt10"]
    live, tabulated = RUNS["live_beam"], RUNS["pressure_beam"]
    open(path, "w") do io
        println(io, "% Written by make_figures.jl")
        for (name, value) in (("Steering", STEERING), ("RampStart", RAMP[1]),
                              ("RampEnd", RAMP[2]), ("SimTime", SIM_TIME),
                              ("SectionTime", SECTION_TIME), ("SettleTime", SETTLE_TIME),
                              ("MaxFailedSolves", MAX_FAILED_SOLVES))
            tex_macro(io, name, value)
        end
        for (suffix, project) in (("Particle", PARTICLE), ("Beam", BEAM))
            wing = structures[project].wings[1]
            tex_macro(io, "NumSections" * suffix, length(wing.vsm_wing.unrefined_sections))
            tex_macro(io, "NumPanels" * suffix, length(wing.vsm_aero.panels))
        end
        println(io, "\\newcommand{\\CaseRows}{")
        for row in table
            model = find_case(row.name).project == BEAM ? "beam" : "particle"
            println(io, case_row(row, model, stats[row.name]))
        end
        println(io, "}")
        tex_macro(io, "StaleAoARatio", ratio_text(stale.aoa_ripple, fresh.aoa_ripple))
        tex_macro(io, "StaleSpeedRatio", ratio_text(stale.node_ripple, fresh.node_ripple))
        tex_macro(io, "LiftRatio", ratio_text(operating_coefficient(live, :cl),
                                              operating_coefficient(tabulated, :cl)))
        tex_macro(io, "MomentRatio", ratio_text(operating_coefficient(live, :cm),
                                                operating_coefficient(tabulated, :cm)))
        section = rows["live_beam"]
        tex_macro(io, "SectionSolves", round(Int, section.flown / section.dt))
        tex_macro(io, "SectionFailedSolves", section.failed_solves)
    end
    return nothing
end

mkpath(FIGURE_PATH)
write_results(joinpath(@__DIR__, "results.tex"), case_rows())
CairoMakie.save(joinpath(FIGURE_PATH, "tracking.pdf"), tracking_figure())
CairoMakie.save(joinpath(FIGURE_PATH, "trajectory.pdf"), trajectory_figure())
CairoMakie.save(joinpath(FIGURE_PATH, "station_loads.pdf"), station_loads_figure())
CairoMakie.save(joinpath(FIGURE_PATH, "polar.pdf"), polar_figure(PANEL))
CairoMakie.save(joinpath(FIGURE_PATH, "pressure.pdf"), pressure_figure(PANEL))
nothing
