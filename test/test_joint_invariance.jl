# SPDX-FileCopyrightText: 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_joint_invariance.jl - A joint may only resist *deformation*. Rigid motion
# of the connected pair must produce no wrench at all. Every joint type runs
# through the same invariants, so a new type inherits the check by being added to
# `joint_cases`.
#
# The corotational argument (Yaw, "2D Corotational Beam Formulation", Walla Walla
# University, 2009): with respect to the co-rotating frame "the rigid body
# rotations and translations are zero and only local strain producing
# deformations remain" (§2). Eq. (16) gives the axial deformation as the e₁
# projection of the relative nodal movement, δuₗ = e₁ᵀδd₂₁, while eq. (18) gives
# the *frame rotation* as the transverse projection, δα = (1/L)e₂ᵀδd₂₁. So the
# transverse part of the relative nodal velocity is rigid frame spin carrying no
# strain: damping it brakes the assembly. Eq. (21), δθₗ = δθ − δα, is the same
# subtraction for the nodal rotations.
#
# Rayleigh damping (Alipour & Zareian, "Study Rayleigh Damping in Structures;
# Uncertainties and Treatments", 14th WCEE 2008): C = αM + βK, eq. (1.4), with β
# in seconds; the modal ratio is ζₙ = (α/2)(1/ωₙ) + (β/2)ωₙ, eq. (2.1). The βK
# half inherits K's rigid-body null space, so invariance is structural. The mass
# term αM is deliberately unavailable — it damps absolute velocity, which for an
# internal joint is exactly the defect these tests exist to catch.
#
# 1. Rigid motion (common translation + spin) conserves linear and angular
#    momentum, for every joint type.
# 2. Logarithmic decrement: one full damped period of a free oscillation decays
#    the velocity by exp(−2πζ/√(1−ζ²)) with ζ = βω/2 from eq. (2.1). Run on the
#    axial and torsional modes, this is what pins β to K rather than to anything
#    else with the right units.

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

@isdefined(test_init!) || include(joinpath(@__DIR__, "util.jl"))

using Test
using SymbolicAWEModels
using KiteUtils
using LinearAlgebra

INVARIANCE_SETTINGS_YAML = """
system:
    log_file: "data/joint_invariance_test"
    g_earth: 0.0
solver:
    solver: "FBDF"
    abs_tol: 1.0e-10
    rel_tol: 1.0e-10
kite:
    model: ""
    foil_file: "ram_air_kite/ram_air_kite_foil.dat"
    physical_model: "joint_invariance_test"
    mass: 0.0
tether:
    cd_tether: 0.958
    unit_damping: 0.0
    unit_stiffness: 0.0
    rho_tether: 724.0
    e_tether: 5.5e10
winch:
    winch_model: "TorqueControlledMachine"
    drum_radius: 0.110
    gear_ratio: 1.0
    inertia_total: 0.024
    f_coulomb: 122.0
    c_vf: 30.6
environment:
    rho_0: 1.225
    v_wind: 0.0
    upwind_dir: -90.0
    upwind_elevation: 0.0
    wind_vec: [0.0, 0.0, 0.0]
    profile_law: 0
"""

"""
    set_rigid_motion!(sys, omega, vel, ref)

Put every body of `sys` on one rigid-body motion: velocity `vel` at `ref`, spinning
at `omega`. Body frames start axis-aligned, so `ω_b == omega`.
"""
function set_rigid_motion!(sys, omega, vel, ref)
    for body in sys.bodies
        body.vel_w .= vel .+ omega × (Vector(body.com_w) .- ref)
        body.ω_b .= omega
    end
    return nothing
end

"""
    pair_momenta(sys, ref) -> (; linear, angular)

Linear momentum and angular momentum about `ref` of every body in `sys`.
"""
function pair_momenta(sys, ref)
    linear = zeros(3)
    angular = zeros(3)
    for body in sys.bodies
        vel = Vector(body.vel_w)
        R_b_to_w = SymbolicAWEModels.quaternion_to_rotation_matrix(body.Q_b_to_w)
        linear .+= body.mass .* vel
        angular .+= body.mass .* ((Vector(body.com_w) .- ref) × vel) .+
                    R_b_to_w * (Vector(body.inertia_principal) .*
                                Vector(body.ω_b))
    end
    return (; linear, angular)
end

"""
    log_decrement(zeta)

Velocity ratio across one full damped period of a free oscillation started from
rest at zero deflection: `exp(−2πζ/√(1−ζ²))`.
"""
log_decrement(zeta) = exp(-2π * zeta / sqrt(1 - zeta^2))

@testset "Joint rigid-motion invariance" begin
    pkg_root = dirname(@__DIR__)
    tmpdir = mktempdir()
    data_path = joinpath(tmpdir, "2plate_kite")
    cp(joinpath(pkg_root, "data", "2plate_kite"), data_path; force=true)
    write(joinpath(data_path, "settings.yaml"), INVARIANCE_SETTINGS_YAML)
    write(joinpath(data_path, "system.yaml"),
        "system:\n  sim_settings: settings.yaml\n")
    set_data_path(data_path)
    set = Settings("system.yaml")

    beam_length = 1.0
    body_mass = 1.5
    inertia = [0.02, 0.02, 0.02]
    EA = 1.0e4
    GJ = 50.0
    axial_stiffness = 500.0

    """Two bodies `beam_length` apart, node A optionally clamped."""
    function pair_bodies(type_a=DYNAMIC)
        return (Body(:nodeA; mass=body_mass, inertia_principal=inertia,
                     pos=[0.0, 0.0, 0.0], type=type_a),
                Body(:nodeB; mass=body_mass, inertia_principal=inertia,
                     pos=[beam_length, 0.0, 0.0]))
    end

    """Two bodies joined by one Timoshenko element."""
    function timoshenko_pair(; damping=0.0, type_a=DYNAMIC)
        node_a, node_b = pair_bodies(type_a)
        joint = TimoshenkoJoint(:joint, :nodeA, :nodeB;
            EA, GA=1500.0, GJ, EIy=100.0, EIz=100.0, shear_coeff=5 / 6, damping)
        return SystemStructure("joint_invariance_test", set;
            bodies=[node_a, node_b], timoshenko_joints=[joint])
    end

    """Two bodies joined by one 6-DOF elastic joint with separated anchors."""
    function elastic_pair(; damping=0.0, type_a=DYNAMIC)
        node_a, node_b = pair_bodies(type_a)
        # Anchors deliberately separated; coincident ones make Δv vanish anyway,
        # and only separated ones exercise the transport couple on body A.
        joint = ElasticJoint(:joint, :nodeA, :nodeB;
            anchor_a=[0.25, 0.0, 0.0], anchor_b=[-0.25, 0.0, 0.0],
            stiffness_axial=axial_stiffness, stiffness_shear=axial_stiffness,
            stiffness_torsion=20.0, stiffness_bending=20.0, damping)
        return SystemStructure("joint_invariance_test", set;
            bodies=[node_a, node_b], elastic_joints=[joint])
    end

    # Add a joint type here and it inherits the rigid-motion invariants.
    joint_cases = [
        ("TimoshenkoJoint, undamped", () -> timoshenko_pair()),
        ("TimoshenkoJoint, Rayleigh", () -> timoshenko_pair(; damping=0.02)),
        ("ElasticJoint, undamped", () -> elastic_pair()),
        ("ElasticJoint, Rayleigh", () -> elastic_pair(; damping=0.02)),
    ]

    for (label, build) in joint_cases
        @testset "$label: rigid motion produces no wrench" begin
            sam = SymbolicAWEModel(set, build())
            sys = sam.sys_struct
            ref = [beam_length / 2, 0.0, 0.0]
            # Off-axis spin so every damping direction is exercised at once.
            omega = 0.7 .* normalize([1.0, 2.0, -3.0])
            vel = [0.3, -0.2, 0.5]
            set_rigid_motion!(sys, omega, vel, ref)
            test_init!(sam; prn=false, reset_vel=false)
            before = pair_momenta(sys, ref)
            spin_before = norm(Vector(sys.bodies[:nodeB].ω_b))
            elapsed = 1.0
            for _ in 1:200
                next_step!(sam; dt=elapsed / 200, vsm_interval=0)
            end
            after = pair_momenta(sys, ref .+ vel .* elapsed)
            spin_after = norm(Vector(sys.bodies[:nodeB].ω_b))
            @info "$label: rigid spin must not be damped." spin_before spin_after
            @test after.linear ≈ before.linear rtol=1e-4
            @test after.angular ≈ before.angular rtol=1e-4
            # |ω| itself is not conserved: the pair needs centripetal force, so it
            # stretches slightly and slows. A dashpot on the transverse relative
            # velocity gives τ = I/(c_t·|d|²) ≈ 4 ms, which wipes the spin out.
            @test spin_after > 0.97 * spin_before
        end
    end

    @testset "Timoshenko: axial logarithmic decrement" begin
        omega_n = sqrt(EA / (beam_length * body_mass))
        zeta = 0.05
        sam = SymbolicAWEModel(set,
            timoshenko_pair(; damping=2 * zeta / omega_n, type_a=STATIC))
        node_b = sam.sys_struct.bodies[:nodeB]
        speed = 0.2
        node_b.vel_w .= [speed, 0.0, 0.0]
        node_b.ω_b .= 0.0
        test_init!(sam; prn=false, reset_vel=false)
        period = 2π / (omega_n * sqrt(1 - zeta^2))
        for _ in 1:2000
            next_step!(sam; dt=period / 2000, vsm_interval=0)
        end
        expected = speed * log_decrement(zeta)
        @info "Axial: ζ = βω/2, δ = 2πζ/√(1−ζ²)." measured=node_b.vel_w[1] expected
        @test node_b.vel_w[1] ≈ expected rtol=5e-3
    end

    @testset "Timoshenko: torsional logarithmic decrement" begin
        omega_n = sqrt(GJ / (beam_length * inertia[1]))
        zeta = 0.05
        sam = SymbolicAWEModel(set,
            timoshenko_pair(; damping=2 * zeta / omega_n, type_a=STATIC))
        node_b = sam.sys_struct.bodies[:nodeB]
        spin = 0.2
        node_b.vel_w .= 0.0
        node_b.ω_b .= [spin, 0.0, 0.0]
        test_init!(sam; prn=false, reset_vel=false)
        period = 2π / (omega_n * sqrt(1 - zeta^2))
        for _ in 1:2000
            next_step!(sam; dt=period / 2000, vsm_interval=0)
        end
        expected = spin * log_decrement(zeta)
        @info "Torsion: ζ = βω/2 on GJ." measured=node_b.ω_b[1] expected
        @test node_b.ω_b[1] ≈ expected rtol=5e-3
    end

    @testset "ElasticJoint: axial logarithmic decrement" begin
        omega_n = sqrt(axial_stiffness / body_mass)
        zeta = 0.05
        sam = SymbolicAWEModel(set,
            elastic_pair(; damping=2 * zeta / omega_n, type_a=STATIC))
        node_b = sam.sys_struct.bodies[:nodeB]
        speed = 0.2
        node_b.vel_w .= [speed, 0.0, 0.0]
        node_b.ω_b .= 0.0
        test_init!(sam; prn=false, reset_vel=false)
        period = 2π / (omega_n * sqrt(1 - zeta^2))
        for _ in 1:2000
            next_step!(sam; dt=period / 2000, vsm_interval=0)
        end
        expected = speed * log_decrement(zeta)
        @info "ElasticJoint axial: ζ = βω/2." measured=node_b.vel_w[1] expected
        @test node_b.vel_w[1] ≈ expected rtol=1e-2
    end

    rm(tmpdir; recursive=true)
end
nothing
