# SPDX-FileCopyrightText: 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_tube_laws.jl - Inflated-tube rigidity laws
#
# The reference numbers are produced by the Breukels (2011) correlations as
# implemented in awegroup/kite_fem (src/kite_fem/BeamElement.py), so a change in
# either the constants or the formula structure shows up here.

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

using Test
using SymbolicAWEModels
using SymbolicAWEModels: breukels_tip_force_coefficients, check_tube_geometry

# (radius [m], pressure [bar]) with the kite_fem tip force at two deflections,
# the collapse deflection, and the torsion factors c1 [N·m] and c2 [m].
TUBE_REFERENCE = (
    (radius=0.075, pressure=0.25, force_030=16.3895448301,
     force_050=23.8757680392, collapse=0.0849074834, moment_collapse=32.6120939374,
     c1=70.7542500000, c2=1.0193629278),
    (radius=0.05, pressure=0.35, force_030=5.7767141347,
     force_050=8.9061320262, collapse=0.0686288281, moment_collapse=11.3915501714,
     c1=77.8063000000, c2=0.2527125067),
    (radius=0.12, pressure=0.5, force_030=101.4918985206,
     force_050=140.7522897670, collapse=0.1690115040, moment_collapse=208.5155140805,
     c1=132.8640000000, c2=3.2549802817),
)

"""kite_fem's secant bending rigidity at tip `deflection`, for the limit check."""
function kite_fem_bending_rigidity(deflection, radius, pressure)
    area = π * radius^2
    shear_modulus = tube_linear_rigidities(radius, pressure)[4] / (π * radius^4 / 2)
    force = breukels_tip_force(deflection, radius, pressure)
    return force / (3 * (deflection -
        force / (TUBE_SHEAR_COEFF * area * shear_modulus)))
end

@testset "Tube rigidity laws" begin

@testset "Breukels correlations match kite_fem" begin
    for case in TUBE_REFERENCE
        radius, pressure = case.radius, case.pressure
        @test breukels_tip_force(0.03, radius, pressure) ≈ case.force_030 rtol=1e-9
        @test breukels_tip_force(0.05, radius, pressure) ≈ case.force_050 rtol=1e-9
        @test breukels_collapse_deflection(radius, pressure) ≈ case.collapse rtol=1e-9

        torsion = tube_torsion_law(radius, pressure)
        @test torsion.c1 ≈ case.c1 rtol=1e-9
        @test torsion.c2 ≈ case.c2 rtol=1e-9

        # GJ0 is the zero-twist limit of kite_fem's GJ = c1·atan(c2·θ)/θ.
        @test tube_linear_rigidities(radius, pressure)[4] ≈ case.c1 * case.c2 rtol=1e-9
    end
end

@testset "Tip force curve" begin
    for case in TUBE_REFERENCE
        radius, pressure = case.radius, case.pressure
        asymptote, slope = breukels_tip_force_coefficients(radius, pressure)
        @test breukels_tip_force(0.0, radius, pressure) == 0.0
        @test breukels_tip_force(1e6, radius, pressure) ≈ asymptote rtol=1e-9
        # `slope` is dP/dδ at zero deflection.
        @test breukels_tip_force(1e-9, radius, pressure) / 1e-9 ≈ slope rtol=1e-6
        @test breukels_tip_force(case.collapse, radius, pressure) ≈
            case.moment_collapse rtol=1e-9
    end
end

@testset "Linear rigidities are the zero-deflection limit" begin
    for case in TUBE_REFERENCE
        radius, pressure = case.radius, case.pressure
        EA, GA, EI0, GJ0 = tube_linear_rigidities(radius, pressure)
        @test all(>(0), (EA, GA, EI0, GJ0))

        # EI0 must be what kite_fem's secant rigidity approaches as δ → 0.
        @test EI0 ≈ kite_fem_bending_rigidity(1e-7, radius, pressure) rtol=1e-5
        # kite_fem clamps δ to 0.03, so its rigidity there is softer than EI0.
        @test kite_fem_bending_rigidity(0.03, radius, pressure) < EI0

        area = π * radius^2
        inertia = π * radius^4 / 4
        @test EA ≈ EI0 * area / inertia rtol=1e-12
        @test GA ≈ GJ0 / (2 * inertia) * area rtol=1e-12
    end
end

@testset "Torsion law" begin
    for case in TUBE_REFERENCE
        law = tube_torsion_law(case.radius, case.pressure)
        @test law.mode === :torsion
        @test law(0.0) ≈ law.c1 * law.c2 rtol=1e-12
        @test law(0.2) ≈ law.c1 * atan(law.c2 * 0.2) / 0.2 rtol=1e-12
        # atan saturates, so the tube softens as it twists.
        @test law(0.5) < law(0.2) < law(1e-9)
    end
end

@testset "Bending law" begin
    for case in TUBE_REFERENCE
        radius, pressure = case.radius, case.pressure
        EI0 = tube_linear_rigidities(radius, pressure)[3]
        law = tube_bending_law(radius, pressure)

        @test law.mode === :bending
        @test law.EI0 ≈ EI0 rtol=1e-12
        @test law.moment_collapse ≈ case.moment_collapse rtol=1e-9
        @test law.moment_knee ≈ 0.5 * case.moment_collapse rtol=1e-9
        @test law.curvature_knee ≈ law.moment_knee / EI0 rtol=1e-12

        # Constant below the knee, continuous across it.
        @test law(0.0) == EI0
        @test law(0.5 * law.curvature_knee) == EI0
        @test law(law.curvature_knee) ≈ EI0 rtol=1e-12
        @test law(1.001 * law.curvature_knee) ≈ EI0 rtol=1e-2

        # Softening above the knee, moment approaching collapse from below.
        curvatures = law.curvature_knee .* (1.5, 3.0, 10.0, 100.0)
        rigidities = [law(κ) for κ in curvatures]
        moments = [κ * law(κ) for κ in curvatures]
        @test issorted(rigidities; rev=true)
        @test issorted(moments)
        @test all(<(law.moment_collapse), moments)
        @test moments[end] ≈ law.moment_collapse rtol=1e-2
        @test law(-0.5) == law(0.5)
    end
end

# The two-branch law (constant EI0, then a power-law approach to the collapse
# moment) cannot follow the Breukels exponential saturation closely: it runs up
# to 30% stiff just above the knee, where the real tube has already softened but
# the law is still linear. The bound below pins that known gap, it is not a target.
@testset "Bending law tracks the Breukels curve" begin
    for case in TUBE_REFERENCE
        radius, pressure = case.radius, case.pressure
        _, GA, _, _ = tube_linear_rigidities(radius, pressure)
        law = tube_bending_law(radius, pressure)
        errors = Tuple{Float64, Float64}[]
        for frac in 0.35:0.025:0.975
            deflection = frac * case.collapse
            moment = breukels_tip_force(deflection, radius, pressure)
            moment > law.moment_knee || continue
            curvature = 3 * (deflection - moment / (TUBE_SHEAR_COEFF * GA))
            push!(errors, (frac, abs(curvature * law(curvature) - moment) / moment))
        end
        worst = maximum(last, errors)
        println("  r=$radius m p=$pressure bar: worst fit error " *
            "$(round(100 * worst; digits=1))% at $(100 * first(argmax(last, errors)))% " *
            "of the collapse deflection")
        @test worst < 0.35
    end
end

@testset "Out-of-range inputs throw" begin
    # Breukels' bending slope goes negative for a thin tube: EI0 would be < 0.
    @test_throws ErrorException tube_linear_rigidities(0.02, 0.1)
    @test_throws ErrorException breukels_tip_force(0.03, 0.02, 0.1)
    @test_throws ErrorException breukels_tip_force_coefficients(0.02, 0.1)
    @test_throws ErrorException tube_bending_law(0.02, 0.1)
    @test_throws ErrorException breukels_membrane_stiffness(0.02, 0.1)

    # c2 = C17·r⁴·log(p) + … flips sign for a fat tube above 1 bar: GJ would be < 0.
    @test_throws ErrorException tube_torsion_law(0.12, 2.0)
    @test_throws ErrorException tube_linear_rigidities(0.12, 2.0)

    for (radius, pressure) in ((0.0, 0.25), (-0.075, 0.25), (0.075, 0.0), (0.075, -0.25))
        @test_throws ErrorException check_tube_geometry(radius, pressure)
        @test_throws ErrorException tube_torsion_law(radius, pressure)
        @test_throws ErrorException breukels_tip_force_coefficients(radius, pressure)
    end

    # The message names the correlation that left its range.
    message = try
        tube_linear_rigidities(0.02, 0.1)
    catch exception
        sprint(showerror, exception)
    end
    @test occursin("bending correlation", message)
    @test occursin("0.02", message)
end

@testset "Comer-Levy section" begin
    radius, pressure_pa, membrane = 0.075, 25_000.0, 5.0e4
    moment_wrinkle = comer_levy_wrinkling_moment(radius, pressure_pa)
    moment_collapse = comer_levy_collapse_moment(radius, pressure_pa)

    @test moment_wrinkle ≈ pressure_pa * π * radius^3 / 2 rtol=1e-12
    @test moment_collapse ≈ 2 * moment_wrinkle rtol=1e-12
    @test comer_levy_bending_stiffness(radius, membrane) ≈ membrane * π * radius^3 rtol=1e-12

    law = comer_levy_bending_law(radius, pressure_pa, membrane)
    @test law.mode === :bending
    @test law.EI0 ≈ comer_levy_bending_stiffness(radius, membrane) rtol=1e-6
    @test law.moment_knee ≈ moment_wrinkle rtol=1e-12
    @test law.moment_collapse ≈ moment_collapse rtol=1e-12
    # Unlike Breukels, the law stays valid past collapse: moment saturates, EI → 0.
    @test 1e4 * law.curvature_knee * law(1e4 * law.curvature_knee) < moment_collapse
    @test law(1e4 * law.curvature_knee) < law(law.curvature_knee)

    # A Breukels-sourced E·t reproduces the Breukels linear rigidity.
    membrane_breukels = breukels_membrane_stiffness(0.075, 0.25)
    @test comer_levy_bending_stiffness(0.075, membrane_breukels) ≈
        tube_linear_rigidities(0.075, 0.25)[3] rtol=1e-12
end

@testset "Membrane rigidities and mass" begin
    radius, membrane = 0.075, 5.0e4
    EA, GA, EI0 = membrane_linear_rigidities(radius, membrane)
    @test EA ≈ membrane * 2π * radius rtol=1e-12
    @test GA ≈ membrane / (2 * (1 + TUBE_POISSON_RATIO)) * 2π * radius rtol=1e-12
    @test EI0 ≈ membrane * π * radius^3 rtol=1e-12

    @test tube_mass(0.075, 3.0; areal_density=0.15) ≈ 2π * 0.075 * 3.0 * 0.15 rtol=1e-12
    @test tube_mass(0.075, 0.0; areal_density=0.15) == 0.0
end

end
nothing
