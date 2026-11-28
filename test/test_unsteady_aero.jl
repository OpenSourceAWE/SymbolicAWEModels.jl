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
#   - a reference frame that is the one VSM gives the panels the lag shifts
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
    panel_apparent_mass, panel_chord_width, wagner_enabled, structure_hash,
    update_run_state!, aero_scatter_entries, apparent_mass_carriers,
    wagner_reference_frame, air_density
using KiteUtils
using LinearAlgebra

"""
    build_unsteady_case(root, tag; apparent_mass, wagner, rates)

Fresh data directory, `Settings`, `SystemDefinition`, its placing transforms, VSM settings
and the model's aero choice for the particle-dynamics 2plate kite, flown with the unsteady
corrections set as given. Each
build gets its own `system_name`: the Wagner constants are substituted into the
emitted equations, so two lags must not share a compiled model.
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
    sys, transforms = load_authoring(
        joinpath(data_path, "particle_structural_geometry.yaml");
        name=tag, set)
    unsteady = UnsteadyAero(; apparent_mass, wagner, wagner_rates=rates)
    aero = Dict{Symbol, SymbolicAWEModels.AbstractAeroModel}(
        :main_wing => ContinuousAero(; unsteady))
    return set, sys, transforms, vsm_set, aero
end

"""Wagner's indicial function at `s` semi-chords, for the wing's own constants."""
wagner_phi(unsteady, s) = 1 - sum(unsteady.wagner_gains[i] *
    exp(-unsteady.wagner_rates[i] * s) for i in 1:2)

@testset verbose=true "Unsteady aero" begin
    root = mktempdir()

    @testset "apparent mass" begin
        set, sys, transforms, vsm_set, aero = build_unsteady_case(root, "am_on";
                                                                  apparent_mass=1.0)
        sam = SymbolicAWEModel(set, sys; transforms, backend=MonolithBackend(), vsm_set,
                               aero)
        # The air goes on the parameters of a running model's integrator.
        test_init!(sam; remake=true)
        mode = sam.aero[1]
        points = wing_points(sam, 1)
        point(i, quantity) = run_value(sam, :points, i, quantity)

        rho = air_density(sam.am, run_value(sam, :wings, 1, :pos_w)[3])
        chord, width = panel_chord_width(mode)
        mass = panel_apparent_mass(mode, rho)
        @test all(mass .> 0)
        @test mass ≈ @. rho * π * chord^2 * width / 4

        # The force weights sum to one per panel, so the nodes get the whole panel.
        weights = zeros(length(mass))
        for (panel, _, force_weight, _) in aero_scatter_entries(mode, sam, 1)
            weights[panel] += force_weight
        end
        @test all(w -> isapprox(w, 1.0; atol=1e-9), weights)

        apply_apparent_mass!(sam)
        @test sum(point(i, :apparent_mass) for i in points) ≈ sum(mass) rtol=1e-9
        @test all(i -> point(i, :apparent_mass) ≥ 0, points)

        # It is inertia, not weight: extra_mass is untouched.
        @test all(i -> point(i, :extra_mass) == 0 ||
                       point(i, :apparent_mass) != point(i, :extra_mass), points)

        # Half the scale is half the air, and zero clears what a previous scale left.
        unsteady_aero(mode).apparent_mass = 0.5
        apply_apparent_mass!(sam)
        @test sum(point(i, :apparent_mass) for i in points) ≈ 0.5 * sum(mass) rtol=1e-9
        unsteady_aero(mode).apparent_mass = 0.0
        apply_apparent_mass!(sam)
        @test all(i -> point(i, :apparent_mass) == 0, points)
    end

    @testset "apparent mass slows the nodes" begin
        # Same force, more inertia: the wing's nodes must accelerate less.
        accels = map((0.0, 1.0)) do scale
            set, sys, transforms, vsm_set, aero = build_unsteady_case(root,
                "am_step_$(scale)"; apparent_mass=scale)
            sam = SymbolicAWEModel(set, sys; transforms, backend=MonolithBackend(), vsm_set,
                                   aero)
            test_init!(sam; remake=true)
            points = wing_points(sam, 1)
            before = [run_value(sam, :points, i, :vel_w) for i in points]
            next_step!(sam)
            update_run_state!(sam.prob, sam.integrator, sam)
            maximum(norm(run_value(sam, :points, i, :vel_w) .- v0)
                    for (i, v0) in zip(points, before))
        end
        @test accels[2] < accels[1]
    end

    @testset "carriers: whatever integrates the node" begin
        # A beam node integrates nothing of its own; air left on it would be inert.
        set_data_path(joinpath(dirname(@__DIR__), "data", "2plate_kite"))
        set = Settings("system.yaml")
        inertia = [0.01, 0.01, 0.01]
        node_a = Body(; name=:node_a, pos_ENU=[0.0, 0.0, 0.0], extra_mass=1.0,
                      extra_inertia_KA=Matrix(Diagonal(inertia)), fixed=true)
        node_b = Body(; name=:node_b, pos_ENU=[1.0, 0.0, 0.0], extra_mass=1.0,
                      extra_inertia_KA=Matrix(Diagonal(inertia)))
        materials = [Material(; name=:breukels2011, law="breukels2011",
            parameters=KiteGeometry.tube_law_parameters(:breukels2011))]
        tube = Tube(; name=:tube, bodies=(:node_a, :node_b), diameter=0.1, pressure=3e4,
            material=:breukels2011,
            model=TimoshenkoTubeModel(; EA=1.0e4, GA=1500.0, GJ=50.0, EIy=100.0,
                                      EIz=100.0, shear_coeff=5 / 6, damping=0.05))
        rider = Point(; name=:rider, pos_ENU=[0.25, 0.0, 0.0], tube=:tube)
        anchored = Point(; name=:anchored, pos_ENU=[1.0, 0.0, 0.0], body=:node_b)
        free = Point(; name=:free, pos_ENU=[2.0, 0.0, 0.0])
        beam = system_definition("carrier_test", set;
            points=[rider, anchored, free], bodies=[node_a, node_b],
            tubes=[tube], materials)

        beam_sam = SymbolicAWEModel(nothing, beam)
        point_idx(name) = KiteGeometry.row_index(beam.points, name)
        body_idx(name) = KiteGeometry.row_index(beam.bodies, name)
        # A tube rider is placed by both bodies, split by where it sits.
        carriers = apparent_mass_carriers(beam_sam, point_idx(:rider))
        @test length(carriers) == 2
        @test sum(share for (_, share) in carriers) ≈ 1.0
        frac = run_value(beam_sam, :points, :rider, :beam_frac)
        @test carriers[1][1] == (:bodies, body_idx(tube.bodies[1]))
        @test carriers[1][2] ≈ 1 - frac
        @test carriers[2][2] ≈ frac

        # An anchored node hands everything to its body, a free node keeps it.
        anchored_carriers = apparent_mass_carriers(beam_sam, point_idx(:anchored))
        @test anchored_carriers[1][1] == (:bodies, body_idx(:node_b))
        @test anchored_carriers[1][2] == 1.0
        free_carriers = apparent_mass_carriers(beam_sam, point_idx(:free))
        @test free_carriers == [((:points, point_idx(:free)), 1.0)]
    end

    @testset "Wagner: off by default" begin
        set, sys, _, vsm_set, aero = build_unsteady_case(root, "wag_default")
        mode = SymbolicAWEModel(set, sys; vsm_set, aero).aero[1]
        @test !wagner_enabled(mode)
        @test unsteady_aero(mode).apparent_mass == 0.0
    end

    @testset "Wagner: the reference frame is the panels' own" begin
        # The lag shifts every panel's α, so its own α is read in their frame.
        set, sys, _, vsm_set, aero = build_unsteady_case(root, "wag_frame")
        sam = SymbolicAWEModel(set, sys; vsm_set, aero)
        panels = sam.aero[1].vsm_aero.panels
        x_ref, z_ref, chord_ref = wagner_reference_frame(sam, 1)
        @test x_ref ≈ normalize(sum(Vector(panel.x_airf) for panel in panels))
        @test z_ref ≈ normalize(sum(Vector(panel.z_airf) for panel in panels))
        @test chord_ref ≈ sum(panel.chord for panel in panels) / length(panels)
    end

    @testset "Wagner: constants are parameters, not structure" begin
        # Retuning a lag must sync, not rebuild: only the on/off switch is
        # structure, because it is what adds or removes the two states.
        set, slow, _, vsm_set, aero = build_unsteady_case(root, "wag_tune"; wagner=true,
                                           rates=(0.02, 0.15))
        slow_sam = SymbolicAWEModel(set, slow; vsm_set, aero)
        slow_hash = structure_hash(slow_sam)
        unsteady_aero(slow_sam.aero[1]).wagner_rates = [0.0455, 0.3]
        unsteady_aero(slow_sam.aero[1]).wagner_gains = [0.2, 0.4]
        @test structure_hash(slow_sam) == slow_hash

        off_set, off, _, vsm_set, off_aero = build_unsteady_case(root, "wag_tune";
                                                                 wagner=false)
        @test structure_hash(SymbolicAWEModel(off_set, off; vsm_set,
                                              aero=off_aero)) != slow_hash
    end

    @testset "Wagner: two states for the whole wing" begin
        # However many panels the wing has, the lag costs exactly two states.
        counts = map((false, true)) do wagner
            set, sys, transforms, vsm_set, aero = build_unsteady_case(root,
                "wag_count_$(wagner)"; wagner)
            sam = SymbolicAWEModel(set, sys; transforms, backend=MonolithBackend(), vsm_set,
                                   aero)
            test_init!(sam; remake=true)
            (length(sam.integrator.u), Int(sam.aero[1].vsm_wing.n_panels))
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
