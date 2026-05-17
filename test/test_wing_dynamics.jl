# SPDX-FileCopyrightText: 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_wing_dynamics.jl - Wing rigid body dynamics tests
#
# Tests quaternion dynamics with AERO_NONE (no aero forces):
# 1. Free fall: translational acc = g
# 2. Constant spin: ω about principal axis stays constant,
#    Q(t) matches analytical solution
# 3. Torque-free precession: transverse ω oscillates at
#    predicted frequency from linearized Euler equations

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

@isdefined(test_init!) || include(joinpath(@__DIR__, "util.jl"))

using Test
using SymbolicAWEModels
using SymbolicAWEModels: KVec3, VortexStepMethod
using KiteUtils
using LinearAlgebra
using Statistics: mean

# ==================== YAML DEFINITIONS ==================== #

WING_FREEFALL_YAML = """
points:
  headers: [name, pos_cad, type, wing_idx, transform_idx,
            extra_mass, body_frame_damping,
            world_frame_damping, area, drag_coeff]
  data:
    - [le_left,   [-0.5, 1.0, 2.0], WING, main_wing,
       ~, 0.5, 0.0, 0.0, 0.0, 0.0]
    - [te_left,   [0.5,  1.0, 2.2], WING, main_wing,
       ~, 0.5, 0.0, 0.0, 0.0, 0.0]
    - [le_center, [-0.5, 0.0, 2.5], WING, main_wing,
       ~, 0.5, 0.0, 0.0, 0.0, 0.0]
    - [te_center, [0.5,  0.0, 2.7], WING, main_wing,
       ~, 0.5, 0.0, 0.0, 0.0, 0.0]
    - [le_right,  [-0.5,-1.0, 2.0], WING, main_wing,
       ~, 0.5, 0.0, 0.0, 0.0, 0.0]
    - [te_right,  [0.5, -1.0, 2.2], WING, main_wing,
       ~, 0.5, 0.0, 0.0, 0.0, 0.0]
    - [ground,    [0.0, 0.0, 0.0],  STATIC, ~,
       ~, 0.0, 0.0, 0.0, 0.0, 0.0]

wings:
  data:
    - name: main_wing
      type: QUATERNION
      aero_mode: AERO_NONE
      transform_idx: 0
      y_damping: 0.0
      aero_z_offset: 0.0
"""

SETTINGS_YAML = """
system:
    log_file: "data/wing_test"
    g_earth: 9.81

solver:
    solver: "FBDF"
    abs_tol: 0.0001
    rel_tol: 0.0001

kite:
    model: ""
    foil_file: "ram_air_kite/ram_air_kite_foil.dat"
    physical_model: "wing_test"
    mass: 0.0
    quasi_static: false

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

# ==================== HELPERS ==================== #

"""Hamilton product of quaternions [w, x, y, z]."""
function quat_multiply(p, q)
    return [
        p[1]*q[1] - p[2]*q[2] - p[3]*q[3] - p[4]*q[4],
        p[1]*q[2] + p[2]*q[1] + p[3]*q[4] - p[4]*q[3],
        p[1]*q[3] - p[2]*q[4] + p[3]*q[1] + p[4]*q[2],
        p[1]*q[4] + p[2]*q[3] - p[3]*q[2] + p[4]*q[1],
    ]
end

"""
Analytical quaternion for constant body-frame spin.

For dQ/dt = 0.5 * Omega(omega) * Q with constant omega,
the solution is Q(t) = Q0 * q_rot(|omega|*t, omega/|omega|).
"""
function analytical_quat_spin(Q0, omega, t)
    omega_mag = norm(omega)
    theta = omega_mag * t / 2
    n = omega / omega_mag
    q_rot = [cos(theta); sin(theta) .* n]
    return quat_multiply(Q0, q_rot)
end

"""
Linearized Euler equation coefficients for torque-free
precession about spin axis `k`.

Returns (p, q, C_pq, C_qp, Omega) where p, q are
transverse axis indices and the linearized dynamics are:
    d omega_p/dt = C_pq * omega_q
    d omega_q/dt = C_qp * omega_p
with precession frequency Omega = sqrt(-C_pq * C_qp).
"""
function precession_coeffs(I_b, k, omega0)
    p, q = sort(collect(setdiff(1:3, k)))
    if k == 1
        # d omega_2/dt = (I3-I1)*omega0/I2 * omega_3
        # d omega_3/dt = (I1-I2)*omega0/I3 * omega_2
        C_pq = (I_b[3] - I_b[1]) * omega0 / I_b[2]
        C_qp = (I_b[1] - I_b[2]) * omega0 / I_b[3]
    elseif k == 2
        # d omega_1/dt = (I2-I3)*omega0/I1 * omega_3
        # d omega_3/dt = (I1-I2)*omega0/I3 * omega_1
        C_pq = (I_b[2] - I_b[3]) * omega0 / I_b[1]
        C_qp = (I_b[1] - I_b[2]) * omega0 / I_b[3]
    else  # k == 3
        # d omega_1/dt = (I2-I3)*omega0/I1 * omega_2
        # d omega_2/dt = (I3-I1)*omega0/I2 * omega_1
        C_pq = (I_b[2] - I_b[3]) * omega0 / I_b[1]
        C_qp = (I_b[3] - I_b[1]) * omega0 / I_b[2]
    end
    @assert C_pq * C_qp < 0 "Rotation about axis $k " *
        "is unstable (intermediate inertia axis)"
    Omega = sqrt(-C_pq * C_qp)
    return (; p, q, C_pq, C_qp, Omega)
end

# ==================== TESTS ==================== #

@testset "Wing Dynamics" begin
    pkg_root = dirname(@__DIR__)
    src_data_path = joinpath(
        pkg_root, "data", "2plate_kite")

    tmpdir = mktempdir()
    data_path = joinpath(tmpdir, "2plate_kite")
    cp(src_data_path, data_path; force=true)

    settings_path = joinpath(data_path, "settings.yaml")
    write(settings_path, SETTINGS_YAML)
    system_path = joinpath(data_path, "system.yaml")
    write(system_path,
        "system:\n  sim_settings: settings.yaml\n")

    set_data_path(data_path)
    set = Settings("system.yaml")

    vsm_set = VortexStepMethod.VSMSettings(
        joinpath(data_path, "vsm_settings.yaml");
        data_prefix=false
    )

    yaml_path = joinpath(data_path, "wing_freefall.yaml")
    write(yaml_path, WING_FREEFALL_YAML)

    sys = load_sys_struct_from_yaml(
        yaml_path;
        system_name="wing_freefall",
        set, vsm_set,
        aero_mode=AERO_NONE
    )

    @testset "Model setup" begin
        @test length(sys.wings) == 1
        wing = sys.wings[:main_wing]
        @test wing.wing_type == SymbolicAWEModels.QUATERNION
        @test wing.aero_mode == AERO_NONE
        @test wing.mass ≈ 3.0  # 6 points * 0.5 kg
        @test length(sys.segments) == 0
    end

    sam = SymbolicAWEModel(set, sys)
    test_init!(sam)

    # ================ Free fall ================ #

    @testset "Free fall acceleration = g" begin
        wing = sam.sys_struct.wings[:main_wing]
        dt = 0.01
        for _ in 1:5
            next_step!(sam; dt, vsm_interval=0)
        end
        vel_before = copy(wing.vel_w)
        t_before = sam.integrator.t
        for _ in 1:10
            next_step!(sam; dt, vsm_interval=0)
        end
        vel_after = copy(wing.vel_w)
        elapsed = sam.integrator.t - t_before
        acc = (vel_after - vel_before) / elapsed

        @test acc[3] ≈ -9.81 atol=0.1
        @test abs(acc[1]) < 0.1
        @test abs(acc[2]) < 0.1
        @test norm(acc) ≈ 9.81 atol=0.15
    end

    # ============= Constant spin ============= #
    # Free-floating wing with initial omega about a
    # principal axis. No torque (AERO_NONE, no tethers)
    # so omega stays constant and Q(t) evolves
    # analytically as Q0 * q_rot(omega*t).

    @testset "Constant spin about principal axis" begin
        wing = sam.sys_struct.wings[:main_wing]
        I_b = collect(wing.inertia_principal)
        k = argmax(I_b)  # max I axis: always stable

        omega0 = 3.0  # rad/s
        omega_init = zeros(3)
        omega_init[k] = omega0
        wing.ω_b .= omega_init
        wing.vel_w .= 0.0
        test_init!(sam; prn=false, reset_vel=false)

        Q0 = copy(wing.Q_b_to_w)
        println("  I_b = $(round.(I_b; digits=4))")
        println("  Spin axis = $k, omega0 = $omega0")
        println("  Q0 = $(round.(Q0; digits=4))")
        println("  ω_b after init = $(wing.ω_b)")

        dt = 0.005
        T_rot = 2pi / omega0
        n_steps = round(Int, 3 * T_rot / dt)
        max_omega_err = 0.0
        max_q_err = 0.0
        max_norm_err = 0.0

        for _ in 1:n_steps
            next_step!(sam; dt, vsm_interval=0)
            t = sam.integrator.t

            # omega should stay constant
            omega_err = norm(wing.ω_b - omega_init)
            max_omega_err = max(max_omega_err, omega_err)

            # Q should match analytical solution
            Q_exp = analytical_quat_spin(Q0, omega_init, t)
            # Q and -Q represent the same rotation
            q_err = min(
                norm(wing.Q_b_to_w - Q_exp),
                norm(wing.Q_b_to_w + Q_exp))
            max_q_err = max(max_q_err, q_err)

            # Quaternion norm should be preserved
            norm_err = abs(norm(wing.Q_b_to_w) - 1.0)
            max_norm_err = max(max_norm_err, norm_err)
        end

        println("  Max |omega - omega0|: $max_omega_err")
        println("  Max |Q - Q_analytical|: $max_q_err")
        println("  Max ||Q| - 1|: $max_norm_err")

        @test max_omega_err < 1e-4
        @test max_q_err < 0.01
        @test max_norm_err < 1e-4
    end

    # ========== Torque-free precession ========== #
    # Spin about max-I axis with small perturbation on
    # a transverse axis. Linearized Euler equations
    # predict oscillation at frequency:
    #   Omega = omega0 * sqrt(|(I_b-I_k)(I_a-I_k)|
    #                         / (I_a * I_b))

    @testset "Torque-free precession" begin
        wing = sam.sys_struct.wings[:main_wing]
        I_b = collect(wing.inertia_principal)
        k = argmax(I_b)  # stable spin axis

        omega0 = 5.0
        eps = 0.05  # small perturbation (eps/omega0 = 1%)
        omega_init = zeros(3)
        omega_init[k] = omega0

        pc = precession_coeffs(I_b, k, omega0)
        omega_init[pc.p] = eps

        wing.ω_b .= omega_init
        wing.vel_w .= 0.0
        test_init!(sam; prn=false, reset_vel=false)

        T_prec = 2pi / pc.Omega
        amp_q = abs(eps * pc.Omega / pc.C_pq)
        println("  I_b = $(round.(I_b; digits=4))")
        println("  Spin axis: $k, transverse: " *
            "($(pc.p), $(pc.q))")
        println("  Omega_predicted = " *
            "$(round(pc.Omega; digits=3)) rad/s")
        println("  T_prec = " *
            "$(round(T_prec; digits=3)) s")
        println("  Expected amp[$(pc.q)] = " *
            "$(round(amp_q; digits=5))")

        dt = 0.002
        t_total = 5 * T_prec
        n_steps = round(Int, t_total / dt)
        times = Float64[]
        omega_p_vals = Float64[]
        omega_q_vals = Float64[]
        omega_k_vals = Float64[]

        for _ in 1:n_steps
            next_step!(sam; dt, vsm_interval=0)
            push!(times, sam.integrator.t)
            push!(omega_p_vals, wing.ω_b[pc.p])
            push!(omega_q_vals, wing.ω_b[pc.q])
            push!(omega_k_vals, wing.ω_b[k])
        end

        # Spin-axis omega stays constant
        omega_k_drift = maximum(abs.(omega_k_vals .- omega0))
        println("  Max spin-axis drift: $omega_k_drift")
        @test omega_k_drift < 0.01

        # Measure precession period from zero crossings
        # of omega_p (which oscillates around 0)
        crossings = Int[]
        for i in 2:lastindex(omega_p_vals)
            if omega_p_vals[i-1] * omega_p_vals[i] < 0
                push!(crossings, i)
            end
        end
        @test length(crossings) >= 4  # at least 2 full periods
        half_periods = [
            times[crossings[i+1]] - times[crossings[i]]
            for i in 1:(length(crossings)-1)]
        T_measured = 2 * mean(half_periods)
        println("  T_measured = " *
            "$(round(T_measured; digits=3)) s " *
            "(predicted: $(round(T_prec; digits=3)))")
        @test T_measured ≈ T_prec rtol=0.05

        # Rotational kinetic energy conservation:
        # E = 0.5 * sum(I_i * omega_i^2)
        E_initial = 0.5 * sum(I_b .* omega_init .^ 2)
        omega_final = collect(wing.ω_b)
        E_final = 0.5 * sum(I_b .* omega_final .^ 2)
        println("  E_initial = $(round(E_initial; digits=4))" *
            ", E_final = $(round(E_final; digits=4))")
        @test E_final ≈ E_initial rtol=0.01

        # Quaternion norm preserved
        @test norm(wing.Q_b_to_w) ≈ 1.0 atol=1e-4
    end

    rm(tmpdir; recursive=true)
end
nothing
