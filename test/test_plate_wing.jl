# SPDX-FileCopyrightText: 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_plate_wing.jl - Flat-plate (polar group) aerodynamics.
#
# Builds the same kps4-style flat-plate kite as a RIGID_DYNAMICS and a
# PARTICLE_DYNAMICS wing and checks the two paths agree: under a uniform
# wind field and from an identical initial configuration the per-plate
# angle of attack and the total aerodynamic force (body frame) must match.

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

@isdefined(test_init!) || include(joinpath(@__DIR__, "util.jl"))

using Test
using SymbolicAWEModels
using SymbolicAWEModels: Point, Segment, Tether, Winch,
    BaseWing, plate_group, Transform, SystemStructure
using KiteUtils
using LinearAlgebra

"""Build a flat-plate kps4-style model with the given wing dynamics type."""
function build_plate_sam(mode, set)
    particles = KiteUtils.get_particles(
        set.height_k, set.h_bridle, set.width, set.m_k)
    pos_kcu, pos_nose, pos_top, pos_right, pos_left =
        particles[2], particles[3], particles[4],
        particles[5], particles[6]

    kite_mass = set.mass
    k_nose = set.rel_nose_mass * kite_mass
    k_top = set.rel_top_mass *
        (1.0 - set.rel_nose_mass) * kite_mass
    k_side = 0.5 * (1.0 - set.rel_top_mass) *
        (1.0 - set.rel_nose_mass) * kite_mass
    set.mass = 0.0

    pre_stress = 0.9975
    pos_map = Dict(:kcu => pos_kcu, :nose => pos_nose,
        :top => pos_top, :right => pos_right,
        :left => pos_left)
    bridle_l0(a, b) =
        norm(pos_map[b] - pos_map[a]) * pre_stress

    points = [
        Point(:ground, zeros(3), STATIC),
        Point(:kcu, pos_kcu, DYNAMIC;
            extra_mass=set.kcu_mass, transform=:main_tf),
        Point(:nose, pos_nose, DYNAMIC;
            extra_mass=k_nose, transform=:main_tf),
        Point(:top, pos_top, WING;
            extra_mass=k_top, wing=:plate_wing,
            transform=:main_tf),
        Point(:right, pos_right, WING;
            extra_mass=k_side, wing=:plate_wing,
            transform=:main_tf),
        Point(:left, pos_left, WING;
            extra_mass=k_side, wing=:plate_wing,
            transform=:main_tf),
    ]

    segments = [
        Segment(:kcu_nose, set, :kcu, :nose;
            l0=bridle_l0(:kcu, :nose), diameter_mm=set.d_line),
        Segment(:right_nose, set, :right, :nose;
            l0=bridle_l0(:right, :nose), diameter_mm=set.d_line),
        Segment(:right_left, set, :right, :left;
            l0=bridle_l0(:right, :left), diameter_mm=set.d_line),
        Segment(:top_right, set, :top, :right;
            l0=bridle_l0(:top, :right), diameter_mm=set.d_line),
        Segment(:left_kcu, set, :left, :kcu;
            l0=bridle_l0(:left, :kcu), diameter_mm=set.d_line),
        Segment(:right_kcu, set, :right, :kcu;
            l0=bridle_l0(:right, :kcu), diameter_mm=set.d_line),
        Segment(:top_left, set, :top, :left;
            l0=bridle_l0(:top, :left), diameter_mm=set.d_line),
        Segment(:left_nose, set, :left, :nose;
            l0=bridle_l0(:left, :nose), diameter_mm=set.d_line),
        Segment(:nose_top, set, :nose, :top;
            l0=bridle_l0(:nose, :top), diameter_mm=set.d_line),
    ]

    tethers = [Tether(:main_tether, set.l_tethers[1];
        start_point=:ground, end_point=:kcu,
        n_segments=set.segments)]
    winches = [Winch(:winch, set, [:main_tether];
        winch_point=:ground)]

    rel_side = set.rel_side_area / 100.0
    K = 1.0 - rel_side
    cl, cd = create_plate_interpolations(
        set.alpha_cl, set.cl_list, set.cd_list;
        alpha_cd=set.alpha_cd)

    groups = [
        plate_group(:main, :top; x_airf=[1,0,0], y_airf=[0,1,0],
            area=set.area, calc_cl=cl, calc_cd=cd,
            drag_corr=0.93 * K, twist=deg2rad(set.alpha_zero)),
        plate_group(:right_tip, :right; x_airf=[1,0,0], y_airf=[0,0,-1],
            area=set.area * rel_side, calc_cl=cl, calc_cd=cd,
            drag_corr=0.93 * K, twist=deg2rad(set.alpha_ztip)),
        plate_group(:left_tip, :left; x_airf=[1,0,0], y_airf=[0,0,1],
            area=set.area * rel_side, calc_cl=cl, calc_cd=cd,
            drag_corr=0.93 * K, twist=deg2rad(set.alpha_ztip)),
    ]

    wing = BaseWing(:plate_wing, [:main, :right_tip, :left_tip],
        [1.0 0 0; 0 1 0; 0 0 1], zeros(3), ones(3);
        dynamics_type=mode, aero_mode=AERO_PLATE,
        z_ref_points=([:right, :left], :top),
        y_ref_points=(:left, :right), origin=:kcu)

    transforms = [Transform(:main_tf,
        deg2rad(set.elevation), 0.0, 0.0;
        base_pos=zeros(3), base_point=:ground,
        wing=:plate_wing)]

    sys = SystemStructure("plate_test", set;
        points, groups, segments, tethers, winches,
        wings=[wing], transforms)
    return SymbolicAWEModel(set, sys)
end

@testset "Plate Wing Tests" begin
    pkg_root = dirname(@__DIR__)
    tmpdir = mktempdir()
    cp(joinpath(pkg_root, "data", "kps4"),
        joinpath(tmpdir, "kps4"); force=true)
    set_data_path(joinpath(tmpdir, "kps4"))

    # Compared at the identical initial configuration (same transform
    # placement, zero velocity) under a uniform wind field. `init!`
    # refreshes the struct from the ODE state, so `wing.aero_force_b`
    # holds the body-frame total aero force at that state — no stepping,
    # which would let the two (different) dynamics diverge.
    results = Dict{WingType, Any}()
    for mode in (RIGID_DYNAMICS, PARTICLE_DYNAMICS)
        set = Settings("system.yaml")
        set.upwind_dir = rad2deg(-pi/2)
        set.profile_law = 0   # constant (uniform) wind field

        sam = build_plate_sam(mode, set)
        init!(sam; prn=false)

        wing = sam.sys_struct.wings[1]
        results[mode] = (
            aoa = [g.aoa for g in sam.sys_struct.groups],
            force = Vector(wing.aero_force_b),
            wing_va = Vector(wing.va_b),
        )
        @testset "$mode builds & inits" begin
            @test all(isfinite, results[mode].force)
            @test norm(results[mode].force) > 0
        end
    end

    rig = results[RIGID_DYNAMICS]
    par = results[PARTICLE_DYNAMICS]

    @testset "Rigid vs particle agree" begin
        # Per-plate AoA uses the wing-level apparent wind and twist,
        # so it is identical between the two dynamics paths.
        @test rig.aoa ≈ par.aoa atol=1e-4
        # Same wing apparent wind at the (identical) initial state.
        @test rig.wing_va ≈ par.wing_va rtol=1e-4
        # Total body-frame aero force matches between the two paths.
        @test rig.force ≈ par.force rtol=1e-3
        println("  rigid    aoa=$(round.(rig.aoa, digits=3))")
        println("  particle aoa=$(round.(par.aoa, digits=3))")
        println("  rigid    F_b=$(round.(rig.force, digits=2))")
        println("  particle F_b=$(round.(par.force, digits=2))")
    end
end
