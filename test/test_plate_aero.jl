# Copyright (c) 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_plate_aero.jl
# AeroPlate's only coverage. It is a flat-plate CL/CD lookup with no VSM ground
# truth, so it cannot join the VSM-referenced contract in test_aero_modes.jl and
# instead owes the reference-free half of it:
#   - it compiles, steps, and stays finite and bounded
#   - a symmetric wing produces no yaw and no side force
#   - mirroring an antisymmetric tip twist mirrors the yaw moment
#
# Geometry is the programmatic kps4 flat-plate kite (1-point STATIC twist
# surfaces, which is what AeroPlate requires), mirroring the construction in
# test_transform.jl.

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

@isdefined(test_init!) || include(joinpath(@__DIR__, "util.jl"))

using Test
using SymbolicAWEModels
using SymbolicAWEModels: Point, Segment, Tether, Winch, PlateWing, TwistSurface,
    Transform, SystemStructure, create_plate_interpolations, PARTICLE_DYNAMICS
using KiteUtils
using LinearAlgebra

"""
    build_plate_kite(data_root)

The kps4 flat-plate kite built programmatically: three 1-point `STATIC` twist
surfaces (`main`, `right_tip`, `left_tip`) carrying one shared CL/CD polar set.
Returns `(set, sys)`.
"""
function build_plate_kite(data_root)
    data_path = joinpath(data_root, "kps4")
    cp(joinpath(dirname(@__DIR__), "data", "kps4"), data_path; force=true)
    set_data_path(data_path)
    set = Settings("system.yaml")
    set.upwind_dir = rad2deg(-pi / 2)

    particles = KiteUtils.get_particles(set.height_k, set.h_bridle,
        set.width, set.m_k)
    pos_kcu, pos_nose, pos_top = particles[2], particles[3], particles[4]
    pos_right, pos_left = particles[5], particles[6]

    kite_mass = set.mass
    k_nose = set.rel_nose_mass * kite_mass
    k_top = set.rel_top_mass * (1.0 - set.rel_nose_mass) * kite_mass
    k_side = 0.5 * (1.0 - set.rel_top_mass) * (1.0 - set.rel_nose_mass) * kite_mass
    set.mass = 0.0

    pos_map = Dict(:kcu => pos_kcu, :nose => pos_nose, :top => pos_top,
        :right => pos_right, :left => pos_left)
    bridle_l0(a, b) = norm(pos_map[b] - pos_map[a]) * 0.9975

    points = [
        Point(:ground, zeros(3), STATIC),
        Point(:kcu, pos_kcu, DYNAMIC; extra_mass=set.kcu_mass, transform=:main_tf),
        Point(:nose, pos_nose, DYNAMIC; extra_mass=k_nose, transform=:main_tf),
        Point(:top, pos_top, DYNAMIC; extra_mass=k_top, wing=:plate_wing,
            transform=:main_tf),
        Point(:right, pos_right, DYNAMIC; extra_mass=k_side, wing=:plate_wing,
            transform=:main_tf),
        Point(:left, pos_left, DYNAMIC; extra_mass=k_side, wing=:plate_wing,
            transform=:main_tf),
    ]
    pairs = [(:kcu, :nose), (:right, :nose), (:right, :left), (:top, :right),
             (:left, :kcu), (:right, :kcu), (:top, :left), (:left, :nose),
             (:nose, :top)]
    segments = [Segment(Symbol(a, :_, b), set, a, b; l0=bridle_l0(a, b),
                        diameter_mm=set.d_line) for (a, b) in pairs]
    tethers = [Tether(:main_tether, set.l_tethers[1]; start_point=:ground,
                      end_point=:kcu, n_segments=set.segments)]
    winches = [Winch(:winch, set, [:main_tether]; winch_point=:ground)]

    rel_side = set.rel_side_area / 100.0
    twist_surfaces = [
        TwistSurface(:main, [:top], STATIC, 0.0; x_airf=[1, 0, 0],
            y_airf=[0, 1, 0], area=set.area, twist=deg2rad(set.alpha_zero)),
        TwistSurface(:right_tip, [:right], STATIC, 0.0; x_airf=[1, 0, 0],
            y_airf=[0, 0, -1], area=set.area * rel_side,
            twist=deg2rad(set.alpha_ztip)),
        TwistSurface(:left_tip, [:left], STATIC, 0.0; x_airf=[1, 0, 0],
            y_airf=[0, 0, 1], area=set.area * rel_side,
            twist=deg2rad(set.alpha_ztip)),
    ]
    cl_interp, cd_interp = create_plate_interpolations(set.alpha_cl, set.cl_list,
        set.cd_list; alpha_cd=set.alpha_cd)
    wing = PlateWing(:plate_wing, [:main, :right_tip, :left_tip],
        cl_interp, cd_interp; dynamics_type=PARTICLE_DYNAMICS,
        z_ref_points=([:right, :left], :top), y_ref_points=(:left, :right),
        origin=:kcu, drag_corr=0.93 * (1.0 - rel_side))
    transforms = [Transform(:main_tf, deg2rad(set.elevation), 0.0, 0.0;
        base_pos=zeros(3), base_point=:ground, wing=:plate_wing)]

    sys = SystemStructure("plate_aero_test", set; points, twist_surfaces,
        segments, tethers, winches, wings=[wing], transforms)
    return set, sys
end

"""
    tip_twist!(sys, base, delta)

Antisymmetric tip twist about the wing's nominal `base` angle: `+delta` on the
left tip, `-delta` on the right. Passing `-delta` gives the mirrored wing.
"""
function tip_twist!(sys, base, delta)
    sys.twist_surfaces[:left_tip].twist = base + delta
    sys.twist_surfaces[:right_tip].twist = base - delta
    return nothing
end

"""
    plate_force_moment(sam, wing, sys, base, delta)

Re-init with the given tip twist and read the wing's total body-frame aero force
and moment.
"""
function plate_force_moment(sam, wing, sys, base, delta)
    tip_twist!(sys, base, delta)
    init!(sam; prn=false)
    next_step!(sam; dt=1e-5)
    return Vector(wing.aero_force_b), Vector(wing.aero_moment_b)
end

@testset "AeroPlate" begin
    root = mktempdir()
    set, sys = build_plate_kite(root)
    wing = sys.wings[1]
    base = deg2rad(set.alpha_ztip)

    sam = SymbolicAWEModel(set, sys)
    test_init!(sam)

    @testset "steps and stays bounded" begin
        force0 = Vector(wing.aero_force_b)
        @test all(isfinite, force0)
        @test norm(force0) > 1.0
        bound = 50.0 * norm(force0)
        for _ in 1:20
            next_step!(sam; dt=0.05)
            @test all(isfinite, wing.aero_force_b)
            @test all(isfinite, wing.aero_moment_b)
            @test norm(wing.aero_force_b) < bound
        end
    end

    # AeroPlate has no VSM reference, so its aero contract is a symmetry. The
    # tips are mirror images (y_airf ±z, equal area), so an antisymmetric tip
    # twist is the one input that must produce antisymmetric yaw; a wrong sign or
    # a swapped tip is invisible without it.
    @testset "yaw antisymmetry" begin
        delta = deg2rad(4.0)
        force_pos, moment_pos = plate_force_moment(sam, wing, sys, base, delta)
        _, moment_neg = plate_force_moment(sam, wing, sys, base, -delta)
        _, moment_flat = plate_force_moment(sam, wing, sys, base, 0.0)

        @test all(isfinite, moment_pos)
        # The twist must actually excite yaw, else the rest is vacuous.
        @test abs(moment_pos[3]) > 1e-3 * norm(force_pos)
        # A symmetric wing yaws far less than a twisted one.
        @test abs(moment_flat[3]) < 0.25 * abs(moment_pos[3])
        # Mirroring the twist mirrors the yaw moment.
        @test isapprox(moment_pos[3], -moment_neg[3]; rtol=0.20)
        println("  [AeroPlate] yaw: Mz(+d)=$(round(moment_pos[3]; sigdigits=3)), ",
            "Mz(-d)=$(round(moment_neg[3]; sigdigits=3)), ",
            "Mz(flat)=$(round(moment_flat[3]; sigdigits=3))")
    end
end
nothing
