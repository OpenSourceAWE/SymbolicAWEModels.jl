# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Under twist_surface twist (steering), the structural strut trailing-edge
# points should stay aligned with the deformed VSM panel trailing
# edges for a RIGID_DYNAMICS wing.

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

using Test
using SymbolicAWEModels
using SymbolicAWEModels: VortexStepMethod
using KiteUtils: init!, next_step!
using LinearAlgebra

pkg_root = dirname(@__DIR__)
set_data_path(joinpath(pkg_root, "data", "2plate_kite"))
struc_yaml = joinpath(
    get_data_path(), "rigid_structural_geometry.yaml")

set = Settings("system.yaml")
set.g_earth = 0.0
vsm_set = VortexStepMethod.VSMSettings(
    joinpath(get_data_path(), "vsm_settings.yaml");
    data_prefix=false)

sys = load_sys_struct_from_yaml(struc_yaml;
    system_name="2plate_kite", set, vsm_set)
sys.winches[:main_winch].brake = true
sam = SymbolicAWEModel(set, sys)
l0_left = sam.sys_struct.segments[:kcu_steering_left].l0
l0_right = sam.sys_struct.segments[:kcu_steering_right].l0

dt = 0.05

init!(sam; prn=false, remake=false, remake_vsm=false)
SymbolicAWEModels.find_steady_state!(sam)

"""
Compare each twist_surface's structural strut chord (body frame, under the
current twist) with the matching deformed VSM panel chord.

Reports the angle between the two chord vectors and their relative length
difference rather than the TE-point distance: both endpoints share the wing
frame, so subtracting the LE cancels the rigid-body constraint drift that
otherwise dominates the raw distance during fast transients.
"""
function twist_chord_diffs(sam)
    wing = sam.sys_struct.wings[1]
    points = sam.sys_struct.points
    twist_surfaces = sam.sys_struct.twist_surfaces
    vsm = wing.vsm_wing
    R = wing.R_b_to_w
    origin = wing.pos_w

    n_unref = vsm.n_unrefined_sections
    theta = zeros(Float64, n_unref)
    for g in twist_surfaces, u in g.unrefined_section_idxs
        theta[u] = g.twist
    end
    VortexStepMethod.unrefined_deform!(vsm, theta)
    refined = vsm.refined_sections

    rows = NamedTuple[]
    for g in twist_surfaces
        i1, i2 = g.point_idxs[1], g.point_idxs[end]
        le_idx, te_idx =
            points[i1].pos_cad[1] < points[i2].pos_cad[1] ?
            (i1, i2) : (i2, i1)
        le_b = R' * (points[le_idx].pos_w - origin)
        te_b = R' * (points[te_idx].pos_w - origin)
        k = argmin([norm(Vector(s.LE_point) - le_b)
                    for s in refined])
        sec_le = Vector(refined[k].LE_point)
        sec_te = Vector(refined[k].TE_point)
        chord_struct = te_b - le_b
        chord_vsm = sec_te - sec_le
        cos_angle = dot(chord_struct, chord_vsm) /
                    (norm(chord_struct) * norm(chord_vsm))
        push!(rows, (name=g.name, twist=g.twist,
            angle=rad2deg(acos(clamp(cos_angle, -1, 1))),
            length_err=abs(norm(chord_vsm) - norm(chord_struct)) /
                       norm(chord_struct)))
    end
    VortexStepMethod.unrefined_deform!(vsm, zeros(Float64, n_unref))
    return rows
end

function report(label, rows)
    for r in rows
        println("[$label] twist_surface $(r.name): twist=",
            round(r.twist; digits=5),
            "  angle=", round(r.angle; digits=7),
            " deg  length_err=", round(r.length_err; digits=9))
    end
end

# At the settled state the coupling is exact to solver precision, so this is the
# sharp check; the steered case below only has to survive constraint drift.
static_rows = twist_chord_diffs(sam)
report("static", static_rows)

@testset "twist chord alignment at rest" begin
    for r in static_rows
        @test r.angle < 1e-3
        @test r.length_err < 1e-5
    end
end

# Ramp steering over a fixed trajectory and check the alignment there. Holding
# this steering runs the twist_surface twist away until the VSM solve diverges,
# so a settled operating point does not exist for this 2-plate config.
steer_mag = 0.03
for step in 1:60
    steer = steer_mag * clamp(step * dt / 2.0, 0.0, 1.0)
    sam.sys_struct.segments[:kcu_steering_left].l0 =
        l0_left - steer
    sam.sys_struct.segments[:kcu_steering_right].l0 =
        l0_right + steer
    next_step!(sam; dt, vsm_interval=1)
end
dyn_rows = twist_chord_diffs(sam)
report("dynamic", dyn_rows)

# Where the trajectory lands is not reproducible across platforms, so the
# thresholds are set by what a broken coupling looks like (degrees of chord
# misalignment at ~40 deg of twist), not by the observed spread.
@testset "twist chord alignment under steering" begin
    @test any(r -> abs(r.twist) > 0.1, dyn_rows)
    for r in dyn_rows
        @test r.angle < 0.5
        @test r.length_err < 0.01
    end
end
