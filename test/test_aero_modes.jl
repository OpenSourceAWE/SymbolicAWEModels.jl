# Copyright (c) 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_aero_modes.jl
# Unified tests for the VSM-family aero modes (AeroNone, AeroDirect,
# ContinuousAero, AeroPressure, AeroLinearized) on the 2plate kite, for every
# supported (aero mode x dynamics type) combination. Two drivers per case:
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
@isdefined(write_pressure_fixture) ||
    include(joinpath(@__DIR__, "pressure_fixture.jl"))

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
    apply_pose!(sam, set, pose; twist=nothing)

Set the main transform (elevation, azimuth, heading) and wind speed from a
`(elevation, azimuth, heading, v_wind)` tuple, then re-init the model. Used to
drive a controlled rigid pose without running the ODE. `twist` prescribes the
station angles, so the pose is held in a deformed state. It is applied
*before* `init!`: `update_sys_struct!` writes `twist` back from the model every
step, so a value set afterwards is overwritten and never reaches the dynamics.
"""
function apply_pose!(sam, set, pose; twist=nothing)
    elevation, azimuth, heading, v_wind = pose
    transform = sam.sys_struct.transforms[:main_transform]
    transform.elevation = elevation
    transform.azimuth = azimuth
    transform.heading = heading
    transform.elevation_vel = 0.0
    transform.azimuth_vel = 0.0
    set.v_wind = v_wind
    isnothing(twist) || set_twist!(sam.sys_struct, twist)
    init!(sam; prn=false)
    next_step!(sam; dt=1e-5, vsm_interval=1)
    return nothing
end

"""
    component_moment_tol(case, moment_ref, force_ref)

Per-component moment budget: the case's relative term on `|M_ref|` plus the
force-proportional floor that pays for a distributed load being represented by a
few structural-point forces.
"""
component_moment_tol(case, moment_ref, force_ref) =
    case.moment_rtol * norm(moment_ref) + case.moment_lever * norm(force_ref)

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

    # Every mode reads the same patched geometry. AeroPressure is the only one that
    # needs the per-node surface aero the patch adds, and the others ignore it, but
    # sharing one fixture is what makes their numbers comparable: the yaw contract
    # below is a symmetry with no VSM reference, so a mode failing it while the
    # others pass says something only when all four fly the same wing.
    write_pressure_fixture(data_path)

    # The moment tolerance is `moment_rtol * |M_ref| + moment_lever * |F_ref|`:
    # the relative term plus a small force-proportional floor. The floor
    # accounts for PARTICLE point-distribution representing a distributed load
    # (and its pitching couple) by a few structural-point forces, whose moment
    # residual scales with force times a fraction of the chord.
    cases = [
        (name="none particle", make=() -> AeroNone(), yaml=particle_yaml,
            data=data_path, vsm_set=vsm_set, dynamics=PARTICLE_DYNAMICS,
            reference=:zero,
            force_rtol=1e-6, moment_rtol=1e-6, moment_lever=0.0,
            drag_rtol=1e-6, lift_rtol=1e-6, side_atol=1e-6,
            yaw_rtol=0.0, motion_rtol=0.0, motion_drag_rtol=0.0),
        (name="none rigid", make=() -> AeroNone(), yaml=rigid_yaml,
            data=data_path, vsm_set=vsm_set, dynamics=RIGID_DYNAMICS,
            reference=:zero,
            force_rtol=1e-6, moment_rtol=1e-6, moment_lever=0.0,
            drag_rtol=1e-6, lift_rtol=1e-6, side_atol=1e-6,
            yaw_rtol=0.0, motion_rtol=0.0, motion_drag_rtol=0.0),
        (name="direct particle", make=() -> AeroDirect(), yaml=particle_yaml,
            data=data_path, vsm_set=vsm_set, dynamics=PARTICLE_DYNAMICS,
            reference=:vsm,
            force_rtol=0.006, moment_rtol=0.10, moment_lever=0.06,
            drag_rtol=0.005, lift_rtol=0.05, side_atol=0.02,
            yaw_rtol=0.08, motion_rtol=0.001, motion_drag_rtol=0.001),
        (name="direct rigid", make=() -> AeroDirect(), yaml=rigid_yaml,
            data=data_path, vsm_set=vsm_set, dynamics=RIGID_DYNAMICS,
            reference=:vsm,
            force_rtol=0.001, moment_rtol=0.001, moment_lever=0.0,
            drag_rtol=0.001, lift_rtol=0.01, side_atol=0.01,
            yaw_rtol=0.08, motion_rtol=0.01, motion_drag_rtol=0.01),
        (name="continuous particle", make=() -> ContinuousAero(),
            yaml=particle_yaml, data=data_path, vsm_set=vsm_set_billow,
            dynamics=PARTICLE_DYNAMICS, reference=:vsm,
            force_rtol=0.006, moment_rtol=0.06, moment_lever=0.04,
            drag_rtol=0.005, lift_rtol=0.05, side_atol=0.02,
            yaw_rtol=0.08, motion_rtol=0.08, motion_drag_rtol=0.002),
        # Twice the particle moment budget because AeroPressure anchors the force
        # to VSM but no couple: the moment follows the Cp distribution, and this
        # fixture's synthetic Cp implies a load centre ~0.13 chord from the Cm in
        # polars/1.csv. That fixture disagreement, not the nearest-node lumping
        # (~0.01 chord, bounded by "moment placement" in test_pressure_aero.jl),
        # is what the budget pays for.
        (name="pressure particle", make=() -> AeroPressure(),
            yaml=particle_yaml, data=data_path, vsm_set=vsm_set,
            dynamics=PARTICLE_DYNAMICS, reference=:vsm,
            force_rtol=0.006, moment_rtol=0.20, moment_lever=0.12,
            drag_rtol=0.005, lift_rtol=0.05, side_atol=0.02,
            yaw_rtol=0.08, motion_rtol=0.005, motion_drag_rtol=0.005),
        (name="linearized rigid", make=() -> AeroLinearized(), yaml=rigid_yaml,
            data=data_path, vsm_set=vsm_set, dynamics=RIGID_DYNAMICS,
            reference=:vsm,
            force_rtol=0.001, moment_rtol=0.001, moment_lever=0.0,
            drag_rtol=0.001, lift_rtol=0.01, side_atol=0.01,
            yaw_rtol=0.08, motion_rtol=0.05, motion_drag_rtol=0.05),
    ]

    for (idx, case) in enumerate(cases)
        @testset "$(case.name)" begin
            set_data_path(case.data)
            sys = load_sys_struct_from_yaml(case.yaml;
                system_name="aero_modes_$(idx)", set, vsm_set=case.vsm_set,
                aero_mode=case.make())
            wing = sys.bodies[1]
            @test wing.dynamics_type == case.dynamics

            sam = SymbolicAWEModel(set, sys)
            test_init!(sam)

            @testset "pose sweep" begin
                max_relF = 0.0; max_dir = 0.0
                max_relM = 0.0; max_mom_use = 0.0
                max_zeroF = 0.0; max_zeroM = 0.0
                max_drag = 0.0; max_lift = 0.0; max_side = 0.0
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

                    # Drag, lift and side each within their own budget: |F| is
                    # lift-dominated, so its norm cannot see a drag error.
                    split = wind_frame_force(force, wing.va_b)
                    split_ref = wind_frame_force(force_ref, wing.va_b)
                    drag_err = abs(split.drag - split_ref.drag) /
                        abs(split_ref.drag)
                    lift_err = abs(split.lift - split_ref.lift) /
                        abs(split_ref.lift)
                    side_err = abs(split.side - split_ref.side) /
                        norm(force_ref)
                    max_drag = max(max_drag, drag_err)
                    max_lift = max(max_lift, lift_err)
                    max_side = max(max_side, side_err)
                    @test drag_err < case.drag_rtol
                    @test lift_err < case.lift_rtol
                    @test side_err < case.side_atol

                    # Positive drag and a physical L/D, reference or not.
                    @test split.drag > 0
                    @test 1.0 < abs(split.lift / split.drag) < 30.0

                    # Componentwise, so a wrong Mz cannot hide inside |M|.
                    moment_tol = component_moment_tol(case, moment_ref, force_ref)
                    for axis in 1:3
                        @test abs(moment[axis] - moment_ref[axis]) <= moment_tol
                    end
                    max_relM = max(max_relM,
                        norm(moment .- moment_ref) / norm(moment_ref))
                    max_mom_use = max(max_mom_use,
                        maximum(abs.(moment .- moment_ref)) / moment_tol)
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
                        "drag=$(round(max_drag; sigdigits=3)) ",
                        "(tol $(case.drag_rtol)); ",
                        "lift=$(round(max_lift; sigdigits=3)) ",
                        "(tol $(case.lift_rtol)); ",
                        "side=$(round(max_side; sigdigits=3)) ",
                        "(tol $(case.side_atol)); ",
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

            # A KINEMATIC wing's body frame is fitted from reference points, so
            # its angular velocity is the rate of that fit — not zero, which is
            # what it used to be hardcoded to. Diagnostic only, but
            # sys_state.turn_rates reads it.
            if case.dynamics == PARTICLE_DYNAMICS
                @testset "particle wing ω_b" begin
                    apply_pose!(sam, set, aero_poses[1])
                    settle_aero!(sam; steps=10)
                    rot_before = copy(wing.R_b_to_w)
                    next_step!(sam; dt=1e-4, vsm_interval=1)
                    rot_after = copy(wing.R_b_to_w)
                    spin = ((rot_after .- rot_before) ./ 1e-4) * rot_after'
                    omega_w = 0.5 .* [spin[3, 2] - spin[2, 3],
                                      spin[1, 3] - spin[3, 1],
                                      spin[2, 1] - spin[1, 2]]
                    omega_measured = rot_after' * omega_w
                    @test norm(omega_measured) > 1e-3
                    @test norm(wing.ω_b) > 1e-3
                    @test isapprox(wing.ω_b, omega_measured;
                        atol=0.05 * norm(omega_measured) + 1e-4)
                    println("  [$(case.name)] ω_b=",
                        "$(round(norm(wing.ω_b); sigdigits=3)) vs measured ",
                        "$(round(norm(omega_measured); sigdigits=3)) rad/s")
                end
            end

            case.reference == :zero && continue

            # The yaw contract, stated as a symmetry so it needs no VSM
            # reference: mirroring an antisymmetric deformation must mirror Mz, a
            # symmetric wing must barely yaw, and the deformation must excite yaw
            # at all. Deformation is a trailing-edge displacement with the nodes
            # pinned (`fix_static`), because prescribed twist is inert on a
            # PARTICLE_DYNAMICS wing — there the free points carry it.
            # Particle only: a RIGID_DYNAMICS wing takes its aero from the body
            # pose and has no shape DOF, so displacing its nodes cannot deform it.
            case.dynamics == PARTICLE_DYNAMICS &&
            @testset "yaw antisymmetry" begin
                apply_pose!(sam, set, aero_poses[1])
                pin_wing_nodes!(sys, wing)
                base_pos = wing_node_positions(sys, wing)
                @test !isempty(base_pos)

                function held(delta)
                    deform_trailing_edge!(sys, wing, base_pos, delta)
                    init!(sam; reinit_sys=false, prn=false)
                    next_step!(sam; dt=1e-5, vsm_interval=1)
                    return model_force_moment(sam, wing)
                end
                force_pos, moment_pos = held(0.05)
                _, moment_neg = held(-0.05)
                _, moment_flat = held(0.0)
                pin_wing_nodes!(sys, wing, false)

                @test abs(moment_pos[3]) > 1e-3 * norm(force_pos)
                @test abs(moment_flat[3]) < 0.25 * abs(moment_pos[3])
                @test isapprox(moment_pos[3], -moment_neg[3]; rtol=case.yaw_rtol)
                println("  [$(case.name)] yaw: Mz(+d)=",
                    "$(round(moment_pos[3]; sigdigits=3)), ",
                    "Mz(-d)=$(round(moment_neg[3]; sigdigits=3)), ",
                    "Mz(flat)=$(round(moment_flat[3]; sigdigits=3))")
            end

            # Parity in a moved, deformed state. The pose sweep calls init!
            # before every measurement, so it only ever samples a pristine
            # zero-velocity wing; this is the driver that exercises the frozen
            # state, which is where a 14% AeroPressure error once hid behind a
            # green 0.6% budget.
            @testset "under-motion parity" begin
                apply_pose!(sam, set, aero_poses[1])
                settle_aero!(sam)
                force, moment = model_force_moment(sam, wing)
                @test all(isfinite, force)
                # The refresh set the per-panel inflow, so this is matched.
                force_ref, moment_ref = vsm_reference_force_moment(wing)
                relF = rel_error(force, force_ref)
                split = wind_frame_force(force, wing.va_b)
                split_ref = wind_frame_force(force_ref, wing.va_b)
                drag_err = abs(split.drag - split_ref.drag) / abs(split_ref.drag)
                @test relF < case.motion_rtol
                @test drag_err < case.motion_drag_rtol
                println("  [$(case.name)] under motion: ",
                    "rel_F=$(round(relF; sigdigits=3)) ",
                    "(tol $(case.motion_rtol)); ",
                    "drag=$(round(drag_err; sigdigits=3)) ",
                    "(tol $(case.motion_drag_rtol))")
            end
        end
    end

    rm(tmpdir; recursive=true)
end
nothing
