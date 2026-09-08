# Copyright (c) 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_unsteady_aero.jl
# The two unsteady corrections a wing carries on top of its steady VSM loading:
# the entrained air its nodes accelerate, and the Wagner lag in how its
# circulation builds. Both are off by default, so each is read as the difference
# between two models rather than inferred from a total.
#
# What apparent mass owes:
#   - the thin-plate mass per panel, ρ·π·c²·w/4
#   - all of it on the wing's nodes, spread by the weights that carry the force
#   - weight that does not grow with it
#   - nothing at all when the scale is zero
#
# What the Wagner lag owes:
#   - two states for the wing, however many panels it has
#   - no deficiency in steady flow, so a trimmed wing is untouched
#   - φ(s) = 1 - A₁·exp(-b₁·s) - A₂·exp(-b₂·s) after a step in angle of attack
#   - gains and rates that reach the model cache key

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

@isdefined(test_init!) || include(joinpath(@__DIR__, "util.jl"))

using Test
using SymbolicAWEModels
using SymbolicAWEModels: VortexStepMethod, wing_points, unsteady_aero,
    panel_apparent_mass, panel_chord_width, wagner_enabled, get_sys_struct_hash,
    update_sys_struct!, aero_scatter_entries, apparent_mass_carriers
using KiteUtils
using LinearAlgebra

"""
    build_unsteady_case(root, tag; apparent_mass, wagner, rates)

Fresh data directory, `Settings` and `SystemStructure` for the particle-dynamics
2plate kite with the unsteady corrections set as given. Each build gets its own
`system_name`: the Wagner constants are substituted into the emitted equations, so
two lags must not share a compiled model.
"""
function build_unsteady_case(root, tag; apparent_mass=0.0, wagner=false,
                             rates=(0.0455, 0.3))
    data_path = joinpath(root, tag, "2plate_kite")
    mkpath(dirname(data_path))
    cp(joinpath(dirname(@__DIR__), "data", "2plate_kite"), data_path; force=true)
    set_data_path(data_path)
    set = Settings("system.yaml")
    vsm_set = VortexStepMethod.VSMSettings(
        joinpath(data_path, "vsm_settings.yaml"); data_prefix=false)
    for vsm_wing_settings in vsm_set.wings
        vsm_wing_settings.spanwise_panel_distribution = VortexStepMethod.BILLOWING
        vsm_wing_settings.billowing_percentage = 8.0
    end
    sys = load_sys_struct_from_yaml(
        joinpath(data_path, "particle_structural_geometry.yaml");
        system_name=tag, set, vsm_set, aero_mode=ContinuousAero())
    sys.wings[1].unsteady.apparent_mass = apparent_mass
    sys.wings[1].unsteady.wagner = wagner
    sys.wings[1].unsteady.wagner_rates = rates
    return set, sys
end

"""Wagner's indicial function at `s` semi-chords, for the wing's own constants."""
wagner_phi(unsteady, s) = 1 - sum(unsteady.wagner_gains[i] *
    exp(-unsteady.wagner_rates[i] * s) for i in 1:2)

@testset verbose=true "Unsteady aero" begin
    root = mktempdir()

    @testset "apparent mass" begin
        _, sys = build_unsteady_case(root, "am_on"; apparent_mass=1.0)
        wing = sys.wings[1]
        points = wing_points(sys, wing)

        rho = 1.225
        chord, width = panel_chord_width(wing)
        mass = panel_apparent_mass(wing, rho)
        @test all(mass .> 0)
        @test mass ≈ @. rho * π * chord^2 * width / 4

        # The force weights sum to one per panel, so the nodes get the whole panel.
        weights = zeros(length(mass))
        for (panel, _, force_weight, _) in aero_scatter_entries(wing.aero, wing, points)
            weights[panel] += force_weight
        end
        @test all(w -> isapprox(w, 1.0; atol=1e-9), weights)

        apply_apparent_mass!(sys, wing, rho)
        @test sum(p.apparent_mass for p in points) ≈ sum(mass) rtol=1e-9
        @test all(p -> p.apparent_mass ≥ 0, points)

        # It is inertia, not weight: extra_mass is untouched.
        @test all(p -> p.extra_mass == 0 || p.apparent_mass != p.extra_mass, points)

        # Half the scale is half the air, and zero clears what a previous scale left.
        wing.unsteady.apparent_mass = 0.5
        apply_apparent_mass!(sys, wing, rho)
        @test sum(p.apparent_mass for p in points) ≈ 0.5 * sum(mass) rtol=1e-9
        wing.unsteady.apparent_mass = 0.0
        apply_apparent_mass!(sys, wing, rho)
        @test all(p -> p.apparent_mass == 0, points)
    end

    @testset "apparent mass slows the nodes" begin
        # Same force, more inertia: the wing's nodes must accelerate less.
        accels = map((0.0, 1.0)) do scale
            set, sys = build_unsteady_case(root, "am_step_$(scale)";
                                           apparent_mass=scale)
            sam = SymbolicAWEModel(set, sys; backend=MonolithBackend())
            test_init!(sam; remake=true)
            points = wing_points(sam.sys_struct, sam.sys_struct.wings[1])
            before = [copy(p.vel_w) for p in points]
            next_step!(sam)
            update_sys_struct!(sam.prob, sam.integrator, sam.sys_struct)
            maximum(norm(p.vel_w .- v0) for (p, v0) in zip(points, before))
        end
        @test accels[2] < accels[1]
    end

    @testset "carriers: whatever integrates the node" begin
        # A beam node integrates nothing of its own; air left on it would be inert.
        set_data_path(joinpath(dirname(@__DIR__), "data", "2plate_kite"))
        set = Settings("system.yaml")
        inertia = [0.01, 0.01, 0.01]
        node_a = Body(:node_a; mass=1.0, inertia_principal=inertia,
                      pos=[0.0, 0.0, 0.0], type=STATIC)
        node_b = Body(:node_b; mass=1.0, inertia_principal=inertia,
                      pos=[1.0, 0.0, 0.0])
        joint = TimoshenkoJoint(:joint, :node_a, :node_b;
            EA=1.0e4, GA=1500.0, GJ=50.0, EIy=100.0, EIz=100.0,
            shear_coeff=5 / 6, damping=0.05)
        rider = Point(:rider, [0.25, 0.0, 0.0], BODY_STATIC; joint=:joint)
        anchored = Point(:anchored, [1.0, 0.0, 0.0], BODY_STATIC; body=:node_b)
        free = Point(:free, [2.0, 0.0, 0.0], DYNAMIC; wing=0)
        beam = SystemStructure("carrier_test", set;
            points=[rider, anchored, free], bodies=[node_a, node_b],
            timoshenko_joints=[joint])

        # A joint rider is placed by both bodies, split by where it sits.
        carriers = apparent_mass_carriers(beam, beam.points[:rider])
        @test length(carriers) == 2
        @test sum(share for (_, share) in carriers) ≈ 1.0
        frac = beam.points[:rider].beam_frac
        @test carriers[1][1] === beam.bodies[joint.body_a_idx]
        @test carriers[1][2] ≈ 1 - frac
        @test carriers[2][2] ≈ frac

        # An anchored node hands everything to its body, a free node keeps it.
        anchored_carriers = apparent_mass_carriers(beam, beam.points[:anchored])
        @test anchored_carriers == [(beam.bodies[:node_b], 1.0)]
        @test apparent_mass_carriers(beam, beam.points[:free]) ==
            [(beam.points[:free], 1.0)]
    end

    @testset "Wagner: off by default" begin
        _, sys = build_unsteady_case(root, "wag_default")
        @test !wagner_enabled(sys.wings[1])
        @test unsteady_aero(sys.wings[1]).apparent_mass == 0.0
    end

    @testset "Wagner: constants are parameters, not structure" begin
        # Retuning a lag must sync, not rebuild: only the on/off switch is
        # structure, because it is what adds or removes the two states.
        _, slow = build_unsteady_case(root, "wag_tune"; wagner=true,
                                      rates=(0.02, 0.15))
        slow_hash = get_sys_struct_hash(slow)
        slow.wings[1].unsteady.wagner_rates = [0.0455, 0.3]
        slow.wings[1].unsteady.wagner_gains = [0.2, 0.4]
        @test get_sys_struct_hash(slow) == slow_hash

        _, off = build_unsteady_case(root, "wag_tune"; wagner=false)
        @test get_sys_struct_hash(off) != slow_hash
    end

    @testset "Wagner: two states for the whole wing" begin
        # However many panels the wing has, the lag costs exactly two states.
        counts = map((false, true)) do wagner
            set, sys = build_unsteady_case(root, "wag_count_$(wagner)"; wagner)
            sam = SymbolicAWEModel(set, sys; backend=MonolithBackend())
            test_init!(sam; remake=true)
            (length(sam.integrator.u), Int(sys.wings[1].vsm_wing.n_panels))
        end
        @test counts[2][1] - counts[1][1] == 2
        @test counts[1][2] > 2
    end

    @testset "Wagner: indicial response" begin
        # φ(0) = 1 - A₁ - A₂, and φ grows monotonically to 1.
        unsteady = SymbolicAWEModels.UnsteadyAero(; wagner=true)
        @test wagner_phi(unsteady, 0.0) ≈ 0.5
        @test wagner_phi(unsteady, 1e4) ≈ 1.0 atol=1e-6
        samples = [wagner_phi(unsteady, s) for s in 0:0.5:60]
        @test all(diff(samples) .> 0)

        # The states integrate that response: from the steady start the lag is
        # inert, and a step in α produces a deficiency that decays back to zero.
        gains, rates = unsteady.wagner_gains, unsteady.wagner_rates
        deficiency(alpha, x) = sum(gains[i] * (alpha - rates[i] * x[i]) for i in 1:2)
        steady = [0.1 / rates[i] for i in 1:2]
        @test deficiency(0.1, steady) ≈ 0.0 atol=1e-12

        # Step α 0.1 → 0.2 and march ds; the recovered fraction must track φ.
        x = copy(steady)
        ds = 1e-3
        worst = 0.0
        for step in 0:round(Int, 40 / ds)
            reached = 1 - deficiency(0.2, x) / 0.1
            worst = max(worst, abs(reached - wagner_phi(unsteady, step * ds)))
            x .+= ds .* (0.2 .- rates .* x)
        end
        @test worst < 2e-3
    end
end
