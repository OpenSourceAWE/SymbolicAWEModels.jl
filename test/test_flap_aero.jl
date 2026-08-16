# Copyright (c) 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_flap_aero.jl — deflecting trailing edge coupled to the continuous pressure
# aero (`AeroPressure`), driven by a real force on the flap.
#
# A STATIC LE body and a DYNAMIC TE (flap) body are hinged at mid-chord; the
# AeroPressure wing's LE/TE points ride the two bodies; a KINEMATIC flap
# twist_surface carries δ = the LE→TE hinge angle into the (α, δ) polars. Putting a
# spanwise moment on the TE body deflects the flap — δ changes and the live RHS
# per-panel force responds (δ enters cl/cd/cm every step). Runs both the lumped
# `ElasticJoint` hinge and the distributed `TimoshenkoJoint` beam, with a δ-swept
# POLAR_MATRICES polar so the flap actually bites.

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

@isdefined(test_init!) || include(joinpath(@__DIR__, "util.jl"))
using Test
using SymbolicAWEModels
const S = SymbolicAWEModels
using SymbolicAWEModels: VortexStepMethod, SimFloat, has_flap,
    twist_surface_deltas
using KiteUtils, LinearAlgebra

function write_delta_swept_fixture(data_path)
    afdir = joinpath(data_path, "airfoils"); mkpath(afdir)
    n_half = 11
    x_top = [0.5 * (1 + cos(pi * (k - 1) / (n_half - 1))) for k in 1:n_half]
    thick(x) = 0.6 * (0.2969*sqrt(x) - 0.1260*x - 0.3516*x^2 + 0.2843*x^3 - 0.1015*x^4)
    xs = Float64[]; ys = Float64[]
    for x in x_top; push!(xs, x); push!(ys, thick(x)); end
    for x in reverse(x_top)[2:end]; push!(xs, x); push!(ys, -thick(x)); end
    n_node = length(xs)
    open(joinpath(afdir, "1.dat"), "w") do io
        println(io, "synthetic")
        for k in 1:n_node; println(io, "$(round(xs[k],digits=5)) $(round(ys[k],digits=5))"); end
    end
    alphas = [-30.0, -15.0, 0.0, 15.0, 30.0]
    open(joinpath(afdir, "1_cp.csv"), "w") do io
        println(io, "alpha,delta," * join(["n$(k-1)" for k in 1:n_node], ","))
        for a in alphas
            cp = [(-3.0*ys[k] - 0.05*a*xs[k]) for k in 1:n_node]
            println(io, "$a,0.0," * join(round.(cp, digits=4), ","))
        end
    end
    open(joinpath(afdir, "1_cf.csv"), "w") do io
        println(io, "alpha,delta," * join(["n$(k-1)" for k in 1:n_node], ","))
        for a in alphas; println(io, "$a,0.0," * join(fill(0.006, n_node), ",")); end
    end
    deltas = [-20.0, -10.0, 0.0, 10.0, 20.0]
    open(joinpath(data_path, "polars", "1.csv"), "w") do io
        println(io, "alpha,delta,cl,cd,cm")
        for a in alphas, d in deltas
            cl = 0.09*a + 0.05*d          # strong δ slope so the flap clearly bites
            cd = 0.02 + 0.0004*a^2
            cm = -0.004*d + 0.001*a
            println(io, "$a,$d,$(round(cl,digits=5)),$(round(cd,digits=5)),$(round(cm,digits=5))")
        end
    end
    geom = joinpath(data_path, "aero_geometry.yaml")
    txt = read(geom, String)
    txt = replace(txt,
        "[1, polars, {csv_file_path: \"polars/1.csv\"}]" =>
        "[1, polars, {csv_file_path: \"polars/1.csv\", dat_file: \"airfoils/1.dat\", " *
        "cp_file: \"airfoils/1_cp.csv\", cf_file: \"airfoils/1_cf.csv\"}]")
    write(geom, txt)
end

function struct_yaml(joint_block)
    """
variables:
  dyneema:
    youngs_modulus: 55.0e9
    damping_per_stiffness: 0.00077
    density: 724.0
points:
  headers: [name, pos_cad, type, wing_idx, transform_idx, body, extra_mass, area, drag_coeff]
  data:
    - [le_left,   [-0.5, 1.0, 2.0], BODY_STATIC, main_wing, main_transform, le_body, 0.1, 0.0, 0.0]
    - [te_left,   [0.5,  1.0, 2.3], BODY_STATIC, main_wing, main_transform, te_body, 0.1, 0.0, 0.0]
    - [le_center, [-0.5, 0.0, 2.5], BODY_STATIC, main_wing, main_transform, le_body, 0.1, 0.0, 0.0]
    - [te_center, [0.5,  0.0, 2.8], BODY_STATIC, main_wing, main_transform, te_body, 0.1, 0.0, 0.0]
    - [le_right,  [-0.5,-1.0, 2.0], BODY_STATIC, main_wing, main_transform, le_body, 0.1, 0.0, 0.0]
    - [te_right,  [0.5, -1.0, 2.3], BODY_STATIC, main_wing, main_transform, te_body, 0.1, 0.0, 0.0]
    - [kcu,       [0.0,  0.0, 0.0],   DYNAMIC, main_wing, main_transform, nothing, 1.0, 0.1, 1.0]
    - [ground,    [0.0,  0.0,-20.0],  STATIC,  main_wing, main_transform, nothing, 0.0, 0.0, 0.0]
segments:
  headers: [name, point_i, point_j, l0, diameter_mm, youngs_modulus, damping_per_stiffness, density, compression_frac]
  data:
    - [bridle_le, le_center, kcu, nothing, 1.0, dyneema, 0.010]
    - [bridle_te, te_center, kcu, nothing, 1.0, dyneema, 0.010]
tethers:
  headers: [name, start_point, end_point, n_segments, youngs_modulus,
            damping_per_stiffness, density, diameter_mm,
            init_stretched_length]
  data:
    - [main_tether, kcu, ground, 4, dyneema, 1.0, 20.0]
winches:
  headers: [name, tether_idxs, winch_point]
  data:
    - [main_winch, [main_tether], ground]
bodies:
  headers: [name, pos, type, transform_idx, mass, inertia_principal]
  data:
    - [le_body, [-0.5, 0.0, 2.17], STATIC,  main_transform, 0.5, [0.05, 0.05, 0.05]]
    - [te_body, [0.5,  0.0, 2.47], DYNAMIC, main_transform, 0.5, [0.05, 0.05, 0.05]]
$(joint_block)twist_surfaces:
  headers: [name, wing, type, points, flap_bodies, flap_axis]
  data:
    - [flap, main_wing, KINEMATIC, [le_left, te_left, le_center, te_center, le_right, te_right], [le_body, te_body], [0.0, 1.0, 0.0]]
wings:
  data:
    - name: main_wing
      dynamics_type: PARTICLE_DYNAMICS
      origin_idx: kcu
      z_ref_points: [kcu, le_center]
      y_ref_points: [le_right, le_left]
transforms:
  data:
    - name: main_transform
      elevation: 50
      azimuth: 0.0
      heading: 0.0
      wing_idx: main_wing
      base_pos: [0.0, 0.0, 0.0]
      base_point_idx: ground
"""
end

ELASTIC_JOINT = """elastic_joints:
  headers: [name, body_a, body_b, anchor_a, anchor_b, stiffness_axial, stiffness_shear, stiffness_torsion, stiffness_bending, damping]
  data:
    - [flap_hinge, le_body, te_body, [0.5, 0.0, 0.15], [-0.5, 0.0, -0.15], 100000.0, 100000.0, 20.0, 6.0, 0.05]
"""
TIMO_JOINT = """timoshenko_joints:
  headers: [name, body_a, body_b, anchor_a, anchor_b, EA, GA, GJ, EIy, EIz, damping]
  data:
    - [flap_beam, le_body, te_body, [0.45, 0.0, 0.135], [-0.45, 0.0, -0.135], 100000.0, 50000.0, 20.0, 6.0, 6.0, 0.05]
"""

function build_flap_model(data_path, joint_block)
    write(joinpath(data_path, "particle_flap_geometry.yaml"), struct_yaml(joint_block))
    set_data_path(data_path)
    set = Settings("system.yaml")
    vsm_set = VortexStepMethod.VSMSettings(
        joinpath(data_path, "vsm_settings.yaml"); data_prefix=false)
    sys = load_sys_struct_from_yaml(
        joinpath(data_path, "particle_flap_geometry.yaml");
        system_name="flap_press_test", set, vsm_set, aero_mode=AeroPressure())
    return SymbolicAWEModel(set, sys), sys
end

flap_angle(sam) = twist_surface_deltas(sam.sys_struct)[sam.sys_struct.twist_surfaces[:flap].idx]

@testset "flap + continuous pressure aero" begin
    pkg_root = dirname(@__DIR__)
    for (label, joint_block) in (("elastic", ELASTIC_JOINT), ("timoshenko", TIMO_JOINT))
        @testset "$label joint" begin
            tmpdir = mktempdir()
            data_path = joinpath(tmpdir, "2plate_kite")
            cp(joinpath(pkg_root, "data", "2plate_kite"), data_path; force=true)
            write_delta_swept_fixture(data_path)
            sam, sys = build_flap_model(data_path, joint_block)
            wing = sys.wings[1]; mode = wing.aero
            te = sys.bodies[:te_body]
            @test mode isa AeroPressure
            @test has_flap(sys.twist_surfaces[:flap])
            @test any(!=(0), mode.panel_twist_surface)

            test_init!(sam; prn=false)
            # A spanwise moment on the TE body deflects the flap about the mid-chord
            # hinge. Compare a +moment settle vs a −moment settle: δ swings, and the
            # live RHS aero force responds because δ enters cl/cd/cm every step.
            function settle_with_moment(moment)
                te.ext_moment_b .= [0.0, moment, 0.0]
                for _ in 1:700; next_step!(sam; dt=0.005, vsm_interval=1); end
                return flap_angle(sam), copy(wing.aero_force_b)
            end
            δ0, F0 = settle_with_moment(0.0)     # aero-only flap equilibrium
            δ1, F1 = settle_with_moment(25.0)    # a strong spanwise moment on the TE
            @test all(isfinite, F0) && all(isfinite, F1)

            println("$label: δ0=$(round(rad2deg(δ0),digits=2))° δ1=$(round(rad2deg(δ1),digits=2))° ",
                    "|F0|=$(round(norm(F0),digits=2)) |F1-F0|=$(round(norm(F1-F0),digits=2))")
            @test abs(δ1 - δ0) > deg2rad(3)               # the flap actually deflected
            @test norm(F1 - F0) > 0.02 * norm(F0)         # aero responded to δ
            rm(tmpdir; recursive=true)
        end
    end
    println("flap + continuous pressure aero: OK")
end
nothing
