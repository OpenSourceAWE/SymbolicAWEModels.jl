# Copyright (c) 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# The flight condition and the case matrix shared by the paper's scripts.

OUTPUT_PATH = joinpath(@__DIR__, "output")
DATA_PATH = joinpath(OUTPUT_PATH, "data")
LOG_PATH = joinpath(OUTPUT_PATH, "logs")
FIGURE_PATH = joinpath(@__DIR__, "figures")
PARTICLE = "system_psm_paper.yaml"
PARTICLE_LIVE = "system_psm_live_paper.yaml"
BEAM = "system_beam.yaml"

SIM_TIME = 20.0          # flown time per case [s]
SECTION_TIME = 1.0       # flown time of the cases the section figures read [s]
STEERING = 0.1           # steering step on top of the nominal steering [-]
RAMP = (0.5, 2.5)        # start and end of the steering ramp [s]
MAX_FAILED_SOLVES = 200  # a case stops once this many VSM solves have failed
WALL_LIMIT = 1200.0      # a case stops once it has run this long [s]
MAX_SOLVER_STEPS = 1000  # a case stops once one output step takes more solver steps
CASE_COLUMNS = (:name, :mode, :dt, :vsm_interval, :flown, :completed, :failed_solves,
                :wall, :t_step, :t_vsm)   # the columns of output/cases.csv

"""
    steering_offset(t)

Steering added to the nominal steering at time `t` [s]: zero, ramped linearly to
`STEERING` over `RAMP`, then held.
"""
steering_offset(t) = STEERING * clamp((t - RAMP[1]) / (RAMP[2] - RAMP[1]), 0, 1)

"""
    mode_label(mode)

Short name of an aerodynamic coupling for tables and legends.
"""
mode_label(::AeroDirect) = "AeroDirect"
mode_label(::ContinuousAero) = "ContinuousAero"
mode_label(mode::AeroPressure) = mode.live_polars ? "AeroPressure+live" : "AeroPressure"

CASES = [
    (name="direct_dt50", project=PARTICLE, mode=AeroDirect(), dt=0.05, vsm_interval=1),
    (name="direct_dt10", project=PARTICLE, mode=AeroDirect(), dt=0.01, vsm_interval=1),
    (name="direct_dt2", project=PARTICLE, mode=AeroDirect(), dt=0.002, vsm_interval=1),
    (name="continuous_dt50", project=PARTICLE, mode=ContinuousAero(), dt=0.05,
     vsm_interval=1),
    (name="continuous_dt10", project=PARTICLE, mode=ContinuousAero(), dt=0.01,
     vsm_interval=1),
    (name="continuous_dt10_vsm5", project=PARTICLE, mode=ContinuousAero(), dt=0.01,
     vsm_interval=5),
    (name="live_dt50", project=PARTICLE_LIVE, mode=AeroPressure(; live_polars=true),
     dt=0.05, vsm_interval=1),
    (name="live_dt10", project=PARTICLE_LIVE, mode=AeroPressure(; live_polars=true),
     dt=0.01, vsm_interval=1),
    (name="pressure_dt50", project=BEAM, mode=AeroPressure(), dt=0.05, vsm_interval=1),
    (name="pressure_dt10", project=BEAM, mode=AeroPressure(), dt=0.01, vsm_interval=1),
]

"""
Cases flown for `SECTION_TIME`, whose end states the section figures compare.
"""
SECTION_CASES = [
    (name="continuous_particle", project=PARTICLE, mode=ContinuousAero(), dt=0.01,
     vsm_interval=1),
    (name="pressure_beam", project=BEAM, mode=AeroPressure(), dt=0.01, vsm_interval=1),
    (name="live_beam", project=BEAM, mode=AeroPressure(; live_polars=true), dt=0.01,
     vsm_interval=1),
]
