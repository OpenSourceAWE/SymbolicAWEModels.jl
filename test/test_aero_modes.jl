# Copyright (c) 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_aero_modes.jl
# Unified tests for the VSM-family aero modes (AeroNone, AeroDirect,
# ContinuousAero, AeroLinearized) on the 2plate kite, for every supported
# (aero mode x dynamics type) combination. Two drivers per case:
#
#   (A) strict pose sweep — init! at a grid of transform poses + wind speeds
#       (no ODE integration, so no solver residual noise), then compare the
#       model's total aero force and moment to a full VSM solve! on the
#       realised panel geometry.
#   (B) loose dynamic run — a short next_step! loop; assert the force/moment
#       stay finite and bounded every step (catches blow-ups and that each
#       combination actually steps).
#
# Reference point for the moment is the WING BODY ORIGIN (wing.pos_w), not the
# COM: VSM's sol.moment is about reference_point=(0,0,0)=body origin, and the
# rigid-body equations transport origin->COM via F x com_offset_b (wing_eqs.jl).
# A COM-referenced sum would differ by that term and hide couple bugs.
#
# AeroPlate is intentionally excluded: it is a flat-plate PlateWing with no VSM
# ground truth and needs its own surfaces YAML, so it cannot share this
# VSM-referenced contract. It keeps its own coverage elsewhere.

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

@isdefined(test_init!) || include(joinpath(@__DIR__, "util.jl"))

using Test
using SymbolicAWEModels
using SymbolicAWEModels: VortexStepMethod, PARTICLE_DYNAMICS,
    RIGID_DYNAMICS, AERO_SCALE_CHORD
using KiteUtils
using LinearAlgebra

"""
    particle_force_scale(wing)

Multiplier applied to distributed `PARTICLE_DYNAMICS` point forces relative to
the raw VSM solution (`1 + aero_scale_chord`, falling back to the package
default). `RIGID_DYNAMICS` uses the body force directly, so its scale is 1.
"""
function particle_force_scale(wing)
    wing.dynamics_type == PARTICLE_DYNAMICS || return 1.0
    return 1.0 + (isfinite(wing.aero_scale_chord) ?
        wing.aero_scale_chord : AERO_SCALE_CHORD)
end

"""
    model_force_moment(sam, wing)

Total aerodynamic force and moment about the wing body origin, in body frame, as
the model produces them. Both dynamics types expose the wing-level aggregates:
`RIGID_DYNAMICS` from the rigid body, `PARTICLE_DYNAMICS` summed from the
distributed point forces and their `r x F` (see `aero_eqs.jl`).
"""
model_force_moment(sam, wing) =
    (Vector(wing.aero_force_b), Vector(wing.aero_moment_b))

"""
    vsm_reference_force_moment(wing)

Full nonlinear VSM `solve!` on the wing's current panel geometry; returns the
total body-frame force and moment about the reference point (= body origin),
each pre-scaled by `particle_force_scale` so it is directly comparable to the
model output.
"""
function vsm_reference_force_moment(wing)
    VortexStepMethod.solve!(wing.vsm_solver, wing.vsm_aero)
    scale = particle_force_scale(wing)
    sol = wing.vsm_solver.sol
    return scale .* Vector(sol.force), scale .* Vector(sol.moment)
end

"""
    rel_error(value, reference)

Relative Euclidean error `|value - reference| / |reference|`.
"""
rel_error(value, reference) = norm(value .- reference) / norm(reference)

"""
    apply_pose!(sam, set, pose)

Set the main transform (elevation, azimuth, heading) and wind speed from a
`(elevation, azimuth, heading, v_wind)` tuple, then re-init the model. Used to
drive a controlled rigid pose without running the ODE.
"""
function apply_pose!(sam, set, pose)
    elevation, azimuth, heading, v_wind = pose
    transform = sam.sys_struct.transforms[:main_transform]
    transform.elevation = elevation
    transform.azimuth = azimuth
    transform.heading = heading
    transform.elevation_vel = 0.0
    transform.azimuth_vel = 0.0
    set.v_wind = v_wind
    init!(sam; prn=false)
    next_step!(sam; dt=1e-5, vsm_interval=1)
    return nothing
end

aero_poses = [
    (deg2rad(60), 0.0,           0.0,           15.0),
    (deg2rad(70), deg2rad(12),   0.0,           16.0),
    (deg2rad(52), deg2rad(-15),  deg2rad(8),    18.0),
    (deg2rad(64), 0.0,           deg2rad(-10),  20.0),
]

@testset "Aero modes" begin
    pkg_root = dirname(@__DIR__)
    src_data_path = joinpath(pkg_root, "data", "2plate_kite")
    tmpdir = mktempdir()
    data_path = joinpath(tmpdir, "2plate_kite")
    cp(src_data_path, data_path; force=true)
    set_data_path(data_path)
    set = Settings("system.yaml")

    vsm_set = VortexStepMethod.VSMSettings(
        joinpath(data_path, "vsm_settings.yaml"); data_prefix=false)
    # ContinuousAero requires the BILLOWING spanwise distribution; the other
    # modes keep the file default (SPLIT_PROVIDED).
    vsm_set_billow = VortexStepMethod.VSMSettings(
        joinpath(data_path, "vsm_settings.yaml"); data_prefix=false)
    for vsm_wing_settings in vsm_set_billow.wings
        vsm_wing_settings.spanwise_panel_distribution =
            VortexStepMethod.BILLOWING
        vsm_wing_settings.billowing_percentage = 8.0
    end
    particle_yaml = joinpath(data_path, "particle_structural_geometry.yaml")
    rigid_yaml = joinpath(data_path, "rigid_structural_geometry.yaml")

    # The moment tolerance is `moment_rtol * |M_ref| + moment_lever * |F_ref|`:
    # the relative term plus a small force-proportional floor. The floor
    # accounts for PARTICLE point-distribution representing a distributed load
    # (and its pitching couple) by a few structural-point forces, whose moment
    # residual scales with force times a fraction of the chord.
    cases = [
        (name="none particle", make=() -> AeroNone(), yaml=particle_yaml,
            vsm_set=vsm_set, dynamics=PARTICLE_DYNAMICS, reference=:zero,
            force_rtol=1e-6, moment_rtol=1e-6, moment_lever=0.0),
        (name="none rigid", make=() -> AeroNone(), yaml=rigid_yaml,
            vsm_set=vsm_set, dynamics=RIGID_DYNAMICS, reference=:zero,
            force_rtol=1e-6, moment_rtol=1e-6, moment_lever=0.0),
        (name="direct particle", make=() -> AeroDirect(), yaml=particle_yaml,
            vsm_set=vsm_set, dynamics=PARTICLE_DYNAMICS, reference=:vsm,
            force_rtol=0.006, moment_rtol=0.10, moment_lever=0.06),
        (name="direct rigid", make=() -> AeroDirect(), yaml=rigid_yaml,
            vsm_set=vsm_set, dynamics=RIGID_DYNAMICS, reference=:vsm,
            force_rtol=0.001, moment_rtol=0.001, moment_lever=0.0),
        (name="continuous particle", make=() -> ContinuousAero(),
            yaml=particle_yaml, vsm_set=vsm_set_billow,
            dynamics=PARTICLE_DYNAMICS, reference=:vsm,
            force_rtol=0.006, moment_rtol=0.06, moment_lever=0.04),
        (name="linearized rigid", make=() -> AeroLinearized(), yaml=rigid_yaml,
            vsm_set=vsm_set, dynamics=RIGID_DYNAMICS, reference=:vsm,
            force_rtol=0.001, moment_rtol=0.001, moment_lever=0.0),
    ]

    for (idx, case) in enumerate(cases)
        @testset "$(case.name)" begin
            sys = load_sys_struct_from_yaml(case.yaml;
                system_name="aero_modes_$(idx)", set, vsm_set=case.vsm_set,
                aero_mode=case.make())
            wing = sys.wings[1]
            @test wing.dynamics_type == case.dynamics

            sam = SymbolicAWEModel(set, sys)
            test_init!(sam)

            @testset "pose sweep" begin
                max_relF = 0.0; max_dir = 0.0
                max_relM = 0.0; max_mom_use = 0.0
                max_zeroF = 0.0; max_zeroM = 0.0
                for pose in aero_poses
                    apply_pose!(sam, set, pose)
                    force, moment = model_force_moment(sam, wing)
                    @test all(isfinite, force)
                    @test all(isfinite, moment)

                    if case.reference == :zero
                        max_zeroF = max(max_zeroF, norm(force))
                        max_zeroM = max(max_zeroM, norm(moment))
                        @test norm(force) < case.force_rtol
                        @test norm(moment) < case.moment_rtol
                        continue
                    end

                    force_ref, moment_ref = vsm_reference_force_moment(wing)
                    @test norm(force_ref) > 1.0
                    relF = rel_error(force, force_ref)
                    max_relF = max(max_relF, relF)
                    @test relF < case.force_rtol
                    cos_force = dot(force, force_ref) /
                        (norm(force) * norm(force_ref))
                    max_dir = max(max_dir,
                        rad2deg(acos(clamp(cos_force, -1.0, 1.0))))
                    @test cos_force > cos(deg2rad(1))

                    moment_tol = case.moment_rtol * norm(moment_ref) +
                        case.moment_lever * norm(force_ref)
                    max_relM = max(max_relM,
                        norm(moment .- moment_ref) / norm(moment_ref))
                    max_mom_use = max(max_mom_use,
                        norm(moment .- moment_ref) / moment_tol)
                    @test norm(moment .- moment_ref) <= moment_tol
                end
                pct(use) = "$(round(100 * use; digits=1))% of budget"
                if case.reference == :zero
                    println("  [$(case.name)] max|F|=",
                        "$(round(max_zeroF; sigdigits=2)) (tol $(case.force_rtol)), ",
                        "max|M|=$(round(max_zeroM; sigdigits=2)) ",
                        "(tol $(case.moment_rtol))")
                else
                    println("  [$(case.name)] ",
                        "rel_F=$(round(max_relF; sigdigits=3)) ",
                        "(tol $(case.force_rtol), $(pct(max_relF/case.force_rtol))); ",
                        "dir=$(round(max_dir; digits=3))° (tol 1°); ",
                        "rel_M=$(round(max_relM; sigdigits=3)), ",
                        "moment $(pct(max_mom_use))")
                end
            end

            @testset "dynamic run" begin
                apply_pose!(sam, set, aero_poses[1])
                dt = 0.05
                force0, _ = model_force_moment(sam, wing)
                bound = 50.0 * max(norm(force0), 1.0)
                for _ in 1:20
                    next_step!(sam; dt, vsm_interval=1)
                    force, moment = model_force_moment(sam, wing)
                    @test all(isfinite, force)
                    @test all(isfinite, moment)
                    @test norm(force) < bound
                end
            end
        end
    end

    rm(tmpdir; recursive=true)
end
nothing
