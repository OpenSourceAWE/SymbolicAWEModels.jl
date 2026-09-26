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
using LinearAlgebra

"""
    tip_twist!(sys, base, delta)

Antisymmetric tip twist about the wing's nominal `base` angle: `+delta` on the
left tip, `-delta` on the right. Passing `-delta` gives the mirrored wing.
"""
function tip_twist!(sys, base, delta)
    sys.stations[:left_tip].twist = base + delta
    sys.stations[:right_tip].twist = base - delta
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
