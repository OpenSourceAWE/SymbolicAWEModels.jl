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
init!(sam; prn=false, remake=false, remake_vsm=false)
SymbolicAWEModels.find_steady_state!(sam)

# Distance between each twist_surface's structural strut TE (body frame, under
# the current twist) and the matching deformed VSM panel TE.
function twist_te_diffs(sam)
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
        sec_te = Vector(refined[k].TE_point)
        push!(rows, (name=g.name, twist=g.twist,
            struct_te=te_b, vsm_te=sec_te,
            diff=norm(sec_te - te_b)))
    end
    VortexStepMethod.unrefined_deform!(vsm, zeros(Float64, n_unref))
    return rows
end

function report(label, rows)
    for r in rows
        println("[$label] twist_surface $(r.name): twist=",
            round(r.twist; digits=5),
            "  diff=", round(r.diff; digits=7))
    end
end

# Ramp steering to a moderate, stable level. Larger amplitudes drive
# the twist_surface twist of this 2-plate config unstable (the twist runs away
# until the VSM solve diverges), so it is not a valid operating point.
dt = 0.05
steer_mag = 0.03
for step in 1:60
    steer = steer_mag * clamp(step * dt / 2.0, 0.0, 1.0)
    sam.sys_struct.segments[:kcu_steering_left].l0 =
        l0_left - steer
    sam.sys_struct.segments[:kcu_steering_right].l0 =
        l0_right + steer
    next_step!(sam; dt, vsm_interval=1)
end
dyn_rows = twist_te_diffs(sam)
report("dynamic", dyn_rows)

sam.sys_struct.segments[:kcu_steering_left].l0 = l0_left - steer_mag
sam.sys_struct.segments[:kcu_steering_right].l0 = l0_right + steer_mag
SymbolicAWEModels.find_steady_state!(sam)
eq_rows = twist_te_diffs(sam)
report("steady", eq_rows)

@testset "twist TE alignment" begin
    @test any(r -> abs(r.twist) > 0.1, eq_rows)
    for r in dyn_rows
        @test r.diff < 1e-3
    end
    for r in eq_rows
        @test r.diff < 1e-3
    end
end
