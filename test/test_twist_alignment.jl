# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Under station twist (steering), the structural strut trailing-edge
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

sys, transforms = load_authoring(struc_yaml;
    name="2plate_kite", set)
sys.winches[:main_winch].brake = true
sam = SymbolicAWEModel(set, sys; transforms, vsm_set)
l0_left = run_value(sam, :segments, :kcu_steering_left, :rest_length)
l0_right = run_value(sam, :segments, :kcu_steering_right, :rest_length)

dt = 0.05

init!(sam; prn=false, remake=false, remake_vsm=false)
SymbolicAWEModels.find_steady_state!(sam)

"""
Compare each station's structural strut chord (body frame, under the
current twist) with the matching deformed VSM panel chord.

Reports the angle between the two chord vectors and their relative length
difference rather than the TE-point distance: both endpoints share the wing
frame, so subtracting the LE cancels the rigid-body constraint drift that
otherwise dominates the raw distance during fast transients.

The angle comes from `atan(|a × b|, a ⋅ b)` rather than `acos`, which near
zero resolves no finer than `sqrt(eps)` rad (8.5e-7 deg per ULP of the
cosine) and would quantize an exact match to a spurious ~1e-6 deg.
"""
function twist_chord_diffs(sam)
    wing_idx = 1
    body_idx = SymbolicAWEModels.wing_body(sam, wing_idx)
    points = sam.definition.points
    stations = SymbolicAWEModels.wing_station_idxs(sam, wing_idx)
    station(i, quantity) = run_value(sam, :spanwise_stations, i, quantity)
    vsm = sam.aero[wing_idx].vsm_wing
    R = run_value(sam, :wings, wing_idx, :R_b_to_w)
    origin = run_value(sam, :wings, wing_idx, :pos_w)
    pos_b(i) = R' * (run_value(sam, :points, i, :pos_w) - origin)
    chord_position(i) =
        SymbolicAWEModels.initial_body_position(sam, body_idx, points[i].pos_ENU)[1]

    n_unref = vsm.n_unrefined_sections
    theta = zeros(Float64, n_unref)
    for g in stations, u in station(g, :unrefined_section_idxs)
        theta[u] = station(g, :twist)
    end
    VortexStepMethod.unrefined_deform!(vsm, theta)
    refined = vsm.refined_sections

    rows = NamedTuple[]
    for g in stations
        station_points = sam.definition.spanwise_stations[g].points
        i1, i2 = station_points[1], station_points[end]
        le_idx, te_idx = chord_position(i1) < chord_position(i2) ? (i1, i2) : (i2, i1)
        le_b = pos_b(le_idx)
        te_b = pos_b(te_idx)
        k = argmin([norm(Vector(s.LE_point) - le_b)
                    for s in refined])
        sec_le = Vector(refined[k].LE_point)
        sec_te = Vector(refined[k].TE_point)
        chord_struct = te_b - le_b
        chord_vsm = sec_te - sec_le
        push!(rows, (name=sam.definition.spanwise_stations[g].name,
            twist=station(g, :twist),
            angle=rad2deg(atan(norm(chord_struct × chord_vsm),
                               dot(chord_struct, chord_vsm))),
            length_err=abs(norm(chord_vsm) - norm(chord_struct)) /
                       norm(chord_struct)))
    end
    VortexStepMethod.unrefined_deform!(vsm, zeros(Float64, n_unref))
    return rows
end

function report(label, rows)
    for r in rows
        println("[$label] station $(r.name): twist=",
            round(r.twist; digits=5),
            "  angle=", r.angle, " deg  length_err=", r.length_err)
    end
end

static_rows = twist_chord_diffs(sam)
report("static", static_rows)

@testset "twist chord alignment at rest" begin
    for r in static_rows
        @test r.angle < 1e-9
        @test r.length_err < 1e-12
    end
end

# Ramp steering until a station twists past 0.1 rad and check the alignment
# there. Holding this steering runs the station twist away until the VSM solve
# diverges, so a settled operating point does not exist for this 2-plate config.
steer_mag = 0.03
for step in 1:100
    steer = steer_mag * clamp(step * dt / 2.0, 0.0, 1.0)
    set_run_value!(sam, :segments, :kcu_steering_left, :rest_length, l0_left - steer)
    set_run_value!(sam, :segments, :kcu_steering_right, :rest_length, l0_right + steer)
    next_step!(sam; dt, vsm_interval=1)
    any(i -> abs(run_value(sam, :spanwise_stations, i, :twist)) > 0.1,
        eachindex(sys.spanwise_stations)) && break
end
dyn_rows = twist_chord_diffs(sam)
report("dynamic", dyn_rows)

# Both coupled paths agree to machine precision here too (~1e-13 deg), so the
# bounds only need enough slack for the trajectory landing elsewhere on another
# platform; a broken coupling misaligns the chord by degrees.
@testset "twist chord alignment under steering" begin
    @test any(r -> abs(r.twist) > 0.1, dyn_rows)
    for r in dyn_rows
        @test r.angle < 1e-5
        @test r.length_err < 1e-10
    end
end
