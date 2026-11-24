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
using SymbolicAWEModels: Transform, SystemDefinition, create_plate_interpolations,
    PARTICLE_DYNAMICS
using KiteUtils
using LinearAlgebra

"""
    build_plate_kite(data_root)

The kps4 flat-plate kite built programmatically: three 1-point `STATIC` twist
surfaces (`main`, `right_tip`, `left_tip`) carrying one shared CL/CD polar set.
Returns `(set, sys, transforms)`.
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

    points = KiteGeometry.Point[
        Point(; name=:ground, pos_ENU=zeros(3), fixed=true),
        Point(; name=:kcu, pos_ENU=pos_kcu, extra_mass=set.kcu_mass),
        Point(; name=:nose, pos_ENU=pos_nose, extra_mass=k_nose),
        Point(; name=:top, pos_ENU=pos_top, extra_mass=k_top),
        Point(; name=:right, pos_ENU=pos_right, extra_mass=k_side),
        Point(; name=:left, pos_ENU=pos_left, extra_mass=k_side),
    ]
    pairs = [(:kcu, :nose), (:right, :nose), (:right, :left), (:top, :right),
             (:left, :kcu), (:right, :kcu), (:top, :left), (:left, :nose),
             (:nose, :top)]
    materials = KiteGeometry.MaterialTable(KiteGeometry.Material[], Dict{String, Any}())
    segments = KiteGeometry.Segment[
        KiteGeometry.line_segment(Symbol(a, :_, b), (a, b), set, materials;
                                  rest_length=bridle_l0(a, b), diameter=0.001 * set.d_line,
                                  model=SegmentSpring())
        for (a, b) in pairs]
    tethers = [KiteGeometry.expand_tether!(points, segments, fill(nothing, length(points)),
        :main_tether, :ground, :kcu, set.segments, set, materials;
        init=(; init_stretched_length=set.l_tethers[1]), model=KiteGeometry.NoModel(),
        segment_model=SegmentSpring())]
    winches = [Winch(; name=:winch, tethers=[:main_tether], winch_point=:ground)]

    rel_side = set.rel_side_area / 100.0
    stations = [
        SpanwiseStation(; name=:main, wing=:plate_wing, points=[:top],
            model=StationModel(; type="STATIC", moment_frac=0.0, chord=[1.0, 0.0, 0.0],
                               y_airf=[0.0, 1.0, 0.0], area=set.area,
                               twist=deg2rad(set.alpha_zero))),
        SpanwiseStation(; name=:right_tip, wing=:plate_wing, points=[:right],
            model=StationModel(; type="STATIC", moment_frac=0.0, chord=[1.0, 0.0, 0.0],
                               y_airf=[0.0, 0.0, -1.0], area=set.area * rel_side,
                               twist=deg2rad(set.alpha_ztip))),
        SpanwiseStation(; name=:left_tip, wing=:plate_wing, points=[:left],
            model=StationModel(; type="STATIC", moment_frac=0.0, chord=[1.0, 0.0, 0.0],
                               y_airf=[0.0, 0.0, 1.0], area=set.area * rel_side,
                               twist=deg2rad(set.alpha_ztip))),
    ]
    wing = KiteGeometry.Wing(; name=:plate_wing, pos_ENU=zeros(3), body=:plate_wing,
        model=WingModel(; aero_mode="AeroPlate", dynamics_type="PARTICLE_DYNAMICS",
            stations=[:main, :right_tip, :left_tip],
            z_ref_points=[[:right, :left], :top], y_ref_points=[:left, :right],
            origin_idx=:kcu, drag_corr=0.93 * (1.0 - rel_side)))
    body = Body(; name=:plate_wing, pos_ENU=zeros(3), extra_mass=0.0,
                extra_inertia_KA=Matrix(Diagonal(ones(3))),
                model=BodyModel(; angular_damping=0.0))
    sys = system_definition("plate_aero_test", set; points, spanwise_stations=stations,
        segments, materials, tethers, winches, wings=[wing], bodies=[body])
    transforms = [Transform(; name=:main_tf, elevation=deg2rad(set.elevation),
        azimuth=0.0, heading=0.0, base_pos=zeros(3), base_point=1, body=1,
        points=collect(2:length(sys.points)), bodies=[1])]
    place!(sys, transforms)
    return set, sys, transforms
end

"""
    tip_twist!(sys, base, delta)

Antisymmetric tip twist about the wing's nominal `base` angle: `+delta` on the
left tip, `-delta` on the right. Passing `-delta` gives the mirrored wing.
"""
function tip_twist!(sys, base, delta)
    sys.spanwise_stations[:left_tip].twist = base + delta
    sys.spanwise_stations[:right_tip].twist = base - delta
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
    set, sys, transforms = build_plate_kite(root)
    base = deg2rad(set.alpha_ztip)

    sam = SymbolicAWEModel(set, sys; transforms)
    wing = sam.state.wings[1]
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
        force_pos, moment_pos = plate_force_moment(sam, wing, sam.state, base, delta)
        _, moment_neg = plate_force_moment(sam, wing, sam.state, base, -delta)
        _, moment_flat = plate_force_moment(sam, wing, sam.state, base, 0.0)

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
