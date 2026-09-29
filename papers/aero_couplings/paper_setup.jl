# Copyright (c) 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# The flight condition and the case matrix shared by the paper's scripts.

OUTPUT_PATH = joinpath(@__DIR__, "output")
DATA_PATH = joinpath(OUTPUT_PATH, "data")
LOG_PATH = joinpath(OUTPUT_PATH, "logs")
FIGURE_PATH = joinpath(@__DIR__, "figures")
PROJECT = "system_beam.yaml"

SIM_TIME = 20.0          # flown time per case [s]
MAX_HEADING = 40.0       # amplitude of the heading setpoint [deg]
PERIOD = 15.0            # period of the heading setpoint [s]
MAX_FAILED_SOLVES = 20   # a case stops once this many VSM solves have failed
WALL_LIMIT = 1200.0      # a case stops once it has run this long [s]

"""
    heading_setpoint(t)

Heading the controller tracks at time `t` [s], in radians.
"""
heading_setpoint(t) = deg2rad(MAX_HEADING) * sin(2pi * t / PERIOD)

"""
    mode_label(mode)

Short name of an aerodynamic coupling for tables and legends.
"""
mode_label(::AeroDirect) = "AeroDirect"
mode_label(::ContinuousAero) = "ContinuousAero"
mode_label(mode::AeroPressure) = mode.live_polars ? "AeroPressure+live" : "AeroPressure"

CASES = [
    (name="direct_dt50", mode=AeroDirect(), dt=0.05, vsm_interval=1),
    (name="direct_dt10", mode=AeroDirect(), dt=0.01, vsm_interval=1),
    (name="direct_dt2", mode=AeroDirect(), dt=0.002, vsm_interval=1),
    (name="direct_dt10_vsm5", mode=AeroDirect(), dt=0.01, vsm_interval=5),
    (name="continuous_dt50", mode=ContinuousAero(), dt=0.05, vsm_interval=1),
    (name="continuous_dt10", mode=ContinuousAero(), dt=0.01, vsm_interval=1),
    (name="continuous_dt10_vsm5", mode=ContinuousAero(), dt=0.01, vsm_interval=5),
    (name="pressure_dt50", mode=AeroPressure(), dt=0.05, vsm_interval=1),
    (name="pressure_dt10", mode=AeroPressure(), dt=0.01, vsm_interval=1),
    (name="live_dt50", mode=AeroPressure(; live_polars=true), dt=0.05, vsm_interval=1),
    (name="live_dt10", mode=AeroPressure(; live_polars=true), dt=0.01, vsm_interval=1),
]
