# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# AeroPressure: distribute VSM's per-section force onto arbitrary structural points
# using the airfoil surface traction (−Cp·n̂ + cf·ŝ) as the placement pattern. The
# per-panel total force is re-derived symbolically each RHS step (live), like
# ContinuousAero; only the scatter differs (surface pattern vs strut couple).

"""
    AeroPressure(; frame_tol_frac=2.0)

VSM aerodynamics whose per-section force is distributed onto arbitrary structural
points (a chord-line skeleton, a strut, or a double-skin membrane of point masses)
by the airfoil **surface traction** pattern (`PARTICLE_DYNAMICS` only). VSM owns the
totals; the traction pattern owns only placement and direction.

Like [`ContinuousAero`](@ref) this is a **continuous** mode: the VSM solve runs every
`vsm_interval` (the full `solve!`), freezing the circulation-derived induced velocity
`v_ind` and the per-node surface traction pattern; but the per-panel **total** force is
re-derived symbolically every RHS step from `v_ind` and the live apparent wind (shared
[`build_panel_force_eqs`](@ref)), so the force responds to apparent-wind changes
between solves. The frozen traction only sets the distribution *shape*; each wing
node's force is its frozen traction plus an equal share of `(live panel force − frozen
pattern net)`, so per panel the point forces sum to the **live** total exactly.

Like [`ContinuousAero`](@ref) the mesh **deforms with the structure**: its unrefined
sections are rebuilt onto the structural LE/TE stations and refined (polars and Cp
contours spanwise-interpolated onto the refined panels), and each refined section's
LE/TE is a live function of the station points' `pos_w`. `frame_tol_frac` is the
frame-alignment guard: construction errors if any surface node maps to a point farther
than `frame_tol_frac ×` the local chord. Carries a [`VSMEngine`](@ref); the no-arg form
is the engine-less marker filled in during wing construction.
"""
mutable struct AeroPressure{E} <: AbstractVSMAero
    engine::Union{Nothing, E}
    "Nearest wing-node index for every (panel, contour node): `station_point[panel][node]`."
    station_point::Vector{Vector{Int64}}
    "Frame-guard tolerance: max node→point distance as a multiple of local chord."
    frame_tol_frac::SimFloat
    "Frozen body-frame induced velocity per refined panel (3 × n_panels)."
    v_ind::Matrix{SimFloat}
    "Frozen section-1 share of each refined panel's chord blend (n_panels)."
    chord_weight::Vector{SimFloat}
    "Frozen body-frame traction per contour node, panel-major (3 × Σ nodes)."
    traction::Matrix{SimFloat}
    "Each wing point's constant force: its nodes' `traction` less its share of `traction_net`."
    point_offset::Dict{Int64, Vector{SimFloat}}
    "Frozen per-panel net traction `Σ_nodes` (3 × n_panels)."
    traction_net::Matrix{SimFloat}
    "Polar callables `(panel_idx, α[, δ])` for cl/cd/cm, read as callable flat params."
    cl::Any
    cd::Any
    cm::Any
    "Per-panel flap twist_surface index (global), or 0 for no flap. Structural (enters the cache hash)."
    panel_twist_surface::Vector{Int64}
    "Left unrefined strut of each refined section (n_panels + 1)."
    section_left_strut::Vector{Int64}
    "Left-strut weight: section = w·strut[left] + (1−w)·strut[left+1]."
    section_left_weight::Vector{SimFloat}
    "Frozen body-frame LE billow offset off the strut line (3 × n_sections)."
    section_le_offset::Matrix{SimFloat}
    "Frozen body-frame TE billow offset off the strut line (3 × n_sections)."
    section_te_offset::Matrix{SimFloat}
    AeroPressure{E}(engine, station_point, frame_tol_frac, v_ind, chord_weight,
        traction, point_offset, traction_net, cl, cd, cm, panel_twist_surface,
        section_left_strut, section_left_weight, section_le_offset,
        section_te_offset) where {E} =
        new{E}(engine, station_point, frame_tol_frac, v_ind, chord_weight,
               traction, point_offset, traction_net, cl, cd, cm,
               panel_twist_surface,
               section_left_strut, section_left_weight, section_le_offset,
               section_te_offset)
end

AeroPressure(; frame_tol_frac=2.0) =
    AeroPressure{VSMEngine}(nothing, Vector{Vector{Int64}}(),
        SimFloat(frame_tol_frac), zeros(SimFloat, 3, 0), SimFloat[],
        zeros(SimFloat, 3, 0),
        Dict{Int64, Vector{SimFloat}}(), zeros(SimFloat, 3, 0),
        nothing, nothing, nothing, Int64[],
        Int64[], SimFloat[], zeros(SimFloat, 3, 0), zeros(SimFloat, 3, 0))
attach_engine!(mode::AeroPressure, engine::VSMEngine) =
    AeroPressure{typeof(engine)}(engine, mode.station_point, mode.frame_tol_frac,
        mode.v_ind, mode.chord_weight,
        mode.traction, mode.point_offset, mode.traction_net,
        mode.cl, mode.cd, mode.cm,
        mode.panel_twist_surface, mode.section_left_strut,
        mode.section_left_weight, mode.section_le_offset, mode.section_te_offset)

is_builtin_aero(::AeroPressure) = true
aero_mode_tag(::AeroPressure) = "press"

"""
    aero_hash_id(mode::AeroPressure)

The surface-node→point map is baked into the generated scatter equations, so it is
structural and enters the model-cache hash (distinguishing it from any stale
frozen-force build that used the default empty id). The panel→flap-twist_surface
map is likewise baked into the δ wiring.
"""
aero_hash_id(mode::AeroPressure) = (mode.station_point, mode.panel_twist_surface,
    mode.section_left_strut, round.(mode.section_left_weight; digits=8),
    round.(mode.section_le_offset; digits=8),
    round.(mode.section_te_offset; digits=8))

# ==================== equation builder ==================== #

"""
    aero_component(mode::AeroPressure, wing::ParticleWing, sys_struct; name, params)

Live per-refined-panel VSM force (shared [`build_panel_force_eqs`](@ref) on the
fixed loaded mesh + live apparent wind) scattered over the frozen surface traction
pattern. Each wing node's force is `Σ over mapped (panel, node) of
traction[node] + (panel_force[panel] − traction_net[panel]) / n_nodes`, so per panel
the point forces sum to the live `panel_force` exactly. Apparent wind and density are
per refined section ([`reconstruct_inflow_sym`](@ref)); `v_ind`, `traction`,
`traction_net` and the `cl/cd/cm` polars are flat params refreshed every
`vsm_interval`.
"""
function aero_component(mode::AeroPressure, wing::ParticleWing, sys_struct;
                        name, params=nothing)
    wing_idx = wing.idx
    vind_p = params.wings[wing_idx].aero.v_ind
    chord_w = params.wings[wing_idx].aero.chord_weight
    cl = params.wings[wing_idx].aero.cl
    cd = params.wings[wing_idx].aero.cd
    cm = params.wings[wing_idx].aero.cm
    traction_p = params.wings[wing_idx].aero.traction
    traction_net_p = params.wings[wing_idx].aero.traction_net

    points = wing_points(sys_struct, wing)
    num_points = length(points)
    panels = wing.vsm_aero.panels
    n_panels = length(panels)
    length(mode.station_point) == n_panels || error(
        "AeroPressure wing $(wing.name): surface→point map not built " *
        "($(length(mode.station_point)) entries for $n_panels panels).")

    # Add a live per-panel flap deflection connector only when this wing has a flap.
    has_flap_coupling = length(mode.panel_twist_surface) == n_panels &&
        any(!=(0), mode.panel_twist_surface)
    connectors = particle_aero_connectors(num_points;
        n_delta=has_flap_coupling ? n_panels : 0)

    # Live geometry: refined-section LE/TE from the structural stations' pos_w.
    column = aero_section_columns(wing, points)
    sec_le, sec_te =
        reconstruct_sections_sym(mode, wing, points, connectors, column)

    sec_va, sec_rho = reconstruct_inflow_sym(mode, wing, connectors, column)

    spanwise = collect(SimFloat, wing.vsm_wing.spanwise_direction)
    scale = 1.0 + (isfinite(wing.aero_scale_chord) ?
        wing.aero_scale_chord : AERO_SCALE_CHORD)

    orient = panel_span_signs(wing, spanwise)
    delta = has_flap_coupling ? collect(connectors.delta) : nothing
    eqs, panel_vars, panel_force, _ = build_panel_force_eqs(
        sec_le, sec_te, sec_va, sec_rho, vind_p, chord_w, cl, cd, cm, spanwise,
        scale, orient; delta)
    vars = particle_unknowns(connectors)
    append!(vars, panel_vars)

    point_num = Dict(point.idx => k for (k, point) in enumerate(points))
    point_force = [zeros(Num, 3) for _ in 1:num_points]
    column = 0
    for i in 1:n_panels
        assigned = mode.station_point[i]
        node_forces = surface_node_forces(traction_p, column + 1, length(assigned),
            panel_force[:, i], [traction_net_p[c, i] for c in 1:3])
        for node in eachindex(assigned)
            column += 1
            k = point_num[assigned[node]]
            point_force[k] = point_force[k] .+ node_forces[node]
        end
    end

    for k in 1:num_points
        eqs = [eqs; connectors.point_force[:, k] ~ point_force[k]]
    end
    return System(eqs, t, vars,
        [vind_p, chord_w, cl, cd, cm, traction_p, traction_net_p]; name)
end

# ==================== build-time panel→flap map + live δ ==================== #

"""
    build_panel_twist_surface_map!(mode, wing, sys_struct)

Map each refined VSM panel to the flap `KINEMATIC` twist_surface of `wing` whose
spanwise station (midpoint of its two flap bodies) is nearest the panel's — the
nearest-station philosophy of [`build_station_point_map!`](@ref). Stored in
`mode.panel_twist_surface` (global twist_surface idx, `0` = no flap); structural,
so it enters [`aero_hash_id`](@ref). All-zeros (no coupling) when the wing has no
flap surface. Generic modes are a no-op.
"""
build_panel_twist_surface_map!(::AbstractAeroModel, wing, sys_struct) = nothing
function build_panel_twist_surface_map!(mode::AeroPressure, wing, sys_struct)
    panels = wing.vsm_aero.panels
    n_panels = length(panels)
    flaps = [ts for ts in sys_struct.twist_surfaces
             if ts.type == KINEMATIC && has_flap(ts) && ts.wing_idx == wing.idx]
    if isempty(flaps)
        mode.panel_twist_surface = zeros(Int64, n_panels)
        return nothing
    end
    rot_cad_to_body = wing.R_b_to_c'
    origin_cad = wing.pos_cad
    spanwise = collect(SimFloat, wing.vsm_wing.spanwise_direction)
    flap_station = map(flaps) do ts
        mid_cad = 0.5 .* (sys_struct.bodies[ts.flap_body_idxs[1]].pos_cad .+
                          sys_struct.bodies[ts.flap_body_idxs[2]].pos_cad)
        dot(rot_cad_to_body * (mid_cad .- origin_cad), spanwise)
    end
    assignment = zeros(Int64, n_panels)
    for (panel_idx, panel) in enumerate(panels)
        le_mid = 0.5 .* (Vector(panel.LE_point_1) .+ Vector(panel.LE_point_2))
        te_mid = 0.5 .* (Vector(panel.TE_point_1) .+ Vector(panel.TE_point_2))
        station = dot(0.5 .* (le_mid .+ te_mid), spanwise)
        best = argmin(abs.(flap_station .- station))
        assignment[panel_idx] = flaps[best].idx
    end
    mode.panel_twist_surface = assignment
    return nothing
end

"""
    twist_surface_deltas(sys_struct) -> Vector{SimFloat}

Live flap deflection δ [rad] per twist_surface (indexed by `twist_surface.idx`;
`0` for non-flap surfaces), from the current flap-body orientations via
[`flap_delta`](@ref). The Julia ground truth mirroring the symbolic
`twist_surface_delta_eqs!`, used to drive `panel.delta` at refresh and by tests.
"""
function twist_surface_deltas(sys_struct)
    twist_surfaces = sys_struct.twist_surfaces
    bodies = sys_struct.bodies
    deltas = zeros(SimFloat, length(twist_surfaces))
    for twist_surface in twist_surfaces
        has_flap(twist_surface) || continue
        R_main = quaternion_to_rotation_matrix(
            bodies[twist_surface.flap_body_idxs[1]].Q_b_to_w)
        R_flap = quaternion_to_rotation_matrix(
            bodies[twist_surface.flap_body_idxs[2]].Q_b_to_w)
        deltas[twist_surface.idx] = flap_delta(twist_surface, R_main, R_flap)
    end
    return deltas
end

"""
    apply_flap_delta!(mode, wing, sys_struct)

Set each VSM panel's `delta` from the live flap deflection of its mapped
twist_surface ([`twist_surface_deltas`](@ref)), before the VSM solve, so the
converged forces and the frozen traction contour track the flap at refresh
cadence. No-op for modes/wings without flap coupling.
"""
apply_flap_delta!(::AbstractAeroModel, wing, sys_struct) = nothing
function apply_flap_delta!(mode::AeroPressure, wing, sys_struct)
    panels = wing.vsm_aero.panels
    (length(mode.panel_twist_surface) == length(panels) &&
     any(!=(0), mode.panel_twist_surface)) || return nothing
    deltas = twist_surface_deltas(sys_struct)
    for (panel_idx, panel) in enumerate(panels)
        ts_idx = mode.panel_twist_surface[panel_idx]
        ts_idx == 0 && continue
        panel.delta = deltas[ts_idx]
    end
    return nothing
end

function restore_flap_delta!(mode::AeroPressure, wing, sys_state)
    isempty(sys_state.flap_angle) && return nothing
    wing.dynamics_type == PARTICLE_DYNAMICS || return nothing
    panels = mode.vsm_aero.panels
    length(mode.panel_twist_surface) == length(panels) || return nothing
    for (panel_idx, panel) in enumerate(panels)
        ts_idx = mode.panel_twist_surface[panel_idx]
        (ts_idx == 0 || ts_idx > length(sys_state.flap_angle)) && continue
        panel.delta = sys_state.flap_angle[ts_idx]
    end
    return nothing
end

# ==================== build-time panel→point map ==================== #

"""
    loft_contour_node(panel, xc, yc, k) -> Vector{SimFloat}

Body-frame position of contour node `k` of a panel's airfoil: the mid-leading-edge
plus `xc[k]` chords along the chord axis (`x_airf`) and `yc[k]` chords along the
surface normal (`z_airf`). The `(xc, yc)` contour is chord-normalised and depends
only on the panel's `delta`, so the loft geometry is static.
"""
function loft_contour_node(panel, xc, yc, k)
    le_mid = 0.5 .* (Vector(panel.LE_point_1) .+ Vector(panel.LE_point_2))
    return le_mid .+ Vector(panel.x_airf) .* (xc[k] * panel.chord) .+
           Vector(panel.z_airf) .* (yc[k] * panel.chord)
end

"""
    build_station_point_map!(mode::AeroPressure, wing, points; prn=false)

Build the static map from each refined panel's surface contour nodes to their
nearest wing-node structural point (body frame), and apply the frame-alignment guard.
Stored in `mode.station_point`; errors on a missing `section_aero` or a mesh that
sits farther than `frame_tol_frac` chords from the points.
"""
function build_station_point_map!(mode::AeroPressure, wing, points; prn=false)
    wing_pts = [p for p in points if p.is_wing_node && p.wing_idx == wing.idx]
    isempty(wing_pts) && error(
        "AeroPressure wing $(wing.name): no wing nodes to receive forces.")
    rot_cad_to_body = wing.R_b_to_c'
    origin_cad = wing.pos_cad
    point_idx = [p.idx for p in wing_pts]
    point_pos_b = [rot_cad_to_body * (p.pos_cad - origin_cad) for p in wing_pts]

    panels = wing.vsm_aero.panels
    station_point = Vector{Vector{Int64}}(undef, length(panels))
    max_chord_ratio = 0.0
    for (panel_idx_local, panel) in enumerate(panels)
        panel.section_aero === nothing && error(
            "AeroPressure wing $(wing.name): panel $panel_idx_local has no " *
            "section_aero (surface contour). Load geometry with .dat/Cp/cf tables.")
        xc, yc, _, _ = VortexStepMethod.section_surface(
            panel.section_aero, 0.0, panel.delta)
        assigned = Vector{Int64}(undef, length(xc))
        for k in eachindex(xc)
            node = loft_contour_node(panel, xc, yc, k)
            min_dist = Inf
            best = point_idx[1]
            for (m, pos_b) in enumerate(point_pos_b)
                dist = norm(node - pos_b)
                if dist < min_dist
                    min_dist = dist
                    best = point_idx[m]
                end
            end
            assigned[k] = best
            max_chord_ratio = max(max_chord_ratio,
                                  min_dist / max(panel.chord, eps()))
        end
        station_point[panel_idx_local] = assigned
    end
    max_chord_ratio > mode.frame_tol_frac && error(
        "AeroPressure wing $(wing.name): max surface-node→point distance is " *
        "$(round(max_chord_ratio, digits=2))× the local chord, exceeding " *
        "frame_tol_frac=$(mode.frame_tol_frac). The VSM mesh and structural " *
        "points are likely misaligned (frame mismatch).")
    mode.station_point = station_point
    prn && println("✓ AeroPressure wing $(wing.name): mapped " *
        "$(sum(length, station_point)) surface nodes across $(length(panels)) " *
        "panels to $(length(wing_pts)) points " *
        "(max $(round(max_chord_ratio, digits=2))×chord).")
    return nothing
end

"""
    init_pressure_buffers!(mode::AeroPressure, wing)

Size the frozen state buffers (`v_ind`, `chord_weight`, `traction`,
`traction_net`) to the current
mesh and set the `cl/cd/cm` polar callables. Run after
[`build_station_point_map!`](@ref) (needs the per-panel node counts).
"""
function init_pressure_buffers!(mode::AeroPressure, wing)
    body_aero = wing.vsm_aero
    n_panels = length(body_aero.panels)
    total_nodes = sum(length, mode.station_point)
    size_frozen_panel_buffers!(mode, n_panels)
    size(mode.traction) == (3, total_nodes) ||
        (mode.traction = zeros(SimFloat, 3, total_nodes))
    size(mode.traction_net) == (3, n_panels) ||
        (mode.traction_net = zeros(SimFloat, 3, n_panels))
    mode.cl = ContinuousPolar(body_aero, VortexStepMethod.calculate_cl)
    mode.cd = ContinuousPolar(body_aero, VortexStepMethod.calculate_cd)
    mode.cm = ContinuousPolar(body_aero, VortexStepMethod.calculate_cm)
    return nothing
end

"""
    setup_aero!(mode::AeroPressure, wing, points, twist_surfaces; prn=false)

Run the generic particle VSM setup (rebuild the unrefined sections onto the
structural LE/TE stations, refine, build `point_to_vsm_point`), freeze the
strut-interpolation caches ([`build_section_interp`](@ref)), then build the
surface→point traction map ([`build_station_point_map!`](@ref)) and size the frozen
buffers/polars ([`init_pressure_buffers!`](@ref)). `PARTICLE_DYNAMICS` only.
"""
function setup_aero!(mode::AeroPressure, wing, points, twist_surfaces; prn=false)
    wing.dynamics_type == PARTICLE_DYNAMICS || error(
        "AeroPressure supports PARTICLE_DYNAMICS wings only; wing " *
        "$(wing.name) is $(wing.dynamics_type).")
    invoke(setup_aero!, Tuple{AbstractVSMAero, Any, Any, Any},
           mode, wing, points, twist_surfaces; prn)
    mode.section_left_strut, mode.section_left_weight,
        mode.section_le_offset, mode.section_te_offset =
        build_section_interp(wing.vsm_wing)
    build_station_point_map!(mode, wing, points; prn)
    init_pressure_buffers!(mode, wing)
    return nothing
end

"""
    validate_aero_structure(mode::AeroPressure, wing, points; prn=false)

Check the generic section↔point mapping plus the surface→point traction map.
"""
function validate_aero_structure(mode::AeroPressure, wing, points; prn=false)
    invoke(validate_aero_structure, Tuple{AbstractVSMAero, Any, Any},
           mode, wing, points; prn)
    isempty(mode.station_point) && error(
        "AeroPressure wing $(wing.idx): surface→point map not built.")
    prn && println("✓ AeroPressure wing $(wing.idx): " *
        "$(length(mode.station_point)) panels mapped.")
    return nothing
end

"""
    remake_aero!(mode::AeroPressure, wing, set, vsm_set, points, twist_surfaces)

Rebuild the VSM objects from edited settings via the generic particle remake
(coarsen onto the structural stations, refine, rebuild `point_to_vsm_point`), then
rebuild the strut-interpolation caches, surface→point map and frozen buffers/polars.
"""
function remake_aero!(mode::AeroPressure, wing, set, vsm_set, points,
                      twist_surfaces)
    invoke(remake_aero!, Tuple{AbstractVSMAero, Any, Any, Any, Any, Any},
           mode, wing, set, vsm_set, points, twist_surfaces)
    mode.section_left_strut, mode.section_left_weight,
        mode.section_le_offset, mode.section_te_offset =
        build_section_interp(wing.vsm_wing)
    build_station_point_map!(mode, wing, points)
    init_pressure_buffers!(mode, wing)
    return nothing
end

# ==================== refresh: solve + freeze pattern ==================== #

"""
    refresh_particle_aero!(mode::AeroPressure, wing, points, va_point_b_vals;
                           vsm_min_wind=0.5, cold_start=false)

Solve at the per-panel apparent wind the symbolic RHS uses
([`set_refined_panel_va!`](@ref)), freezing the induced velocity and the per-node
surface traction pattern ([`freeze_traction_pattern!`](@ref)). Point forces are
re-derived symbolically each RHS step and read back into each node's
`aero_force_b`. Below `vsm_min_wind` all frozen buffers are zeroed, the per-point
offsets included, so no stale constant force survives.
"""
function refresh_particle_aero!(mode::AeroPressure, wing, points,
                                va_point_b_vals; vsm_min_wind=0.5,
                                cold_start=false)
    if norm(wing.va_b) < vsm_min_wind
        fill!(mode.v_ind, 0.0)
        fill!(mode.traction, 0.0)
        fill!(mode.traction_net, 0.0)
        accumulate_point_offset!(mode)
        return nothing
    end
    update_vsm_wing_from_structure!(wing, points)
    set_refined_panel_va!(mode, wing, points, va_point_b_vals)
    solve_and_freeze_circulation!(mode, wing; cold_start)
    freeze_traction_pattern!(mode, wing)
    any(!isfinite, mode.traction) && throw(AssertionError(
        "AeroPressure: non-finite traction pattern on wing $(wing.idx)"))
    return nothing
end

"""
    contour_winding(section) -> Float64

`+1` when the closed contour `section` of `(chordwise, normal)` pairs runs
counter-clockwise, `-1` when it runs clockwise, by the shoelace area. Orienting a
normal off the loop itself works on a cambered canopy, where both surfaces can lie
on one side of the chord line and "away from the chord" points inward.
"""
function contour_winding(section)
    n = length(section)
    area = sum(section[k][1] * section[mod1(k + 1, n)][2] -
               section[mod1(k + 1, n)][1] * section[k][2] for k in 1:n)
    return area >= 0 ? 1.0 : -1.0
end

"""
    freeze_traction_pattern!(mode::AeroPressure, wing)

Fill the frozen traction params from the converged VSM solve: for each refined panel
loft the airfoil surface contour at the effective sectional α (`sol.alpha_dist`), build
the per-segment traction `(−Cp·n̂ + cf·ŝ)·q·dA` (`n̂` outward, `ŝ` along the chord) into
`mode.traction` (panel-major node order, matching [`aero_component`](@ref)), and store
the per-panel net in `mode.traction_net`. The net anchors the symbolic scatter so each
panel's point forces sum to the live VSM total.
"""
function freeze_traction_pattern!(mode::AeroPressure, wing)
    sol = wing.vsm_solver.sol
    panels = wing.vsm_aero.panels
    dynamic_pressure = 0.5 * wing.vsm_solver.density * dot(wing.va_b, wing.va_b)
    column = 0
    for (panel_idx, panel) in enumerate(panels)
        xc, yc, cp, cf = VortexStepMethod.section_surface(
            panel.section_aero, sol.alpha_dist[panel_idx], panel.delta)
        n_nodes = length(xc)
        chord_axis = Vector(panel.x_airf)
        normal_axis = Vector(panel.z_airf)
        pos = [loft_contour_node(panel, xc, yc, k) for k in 1:n_nodes]
        le_mid = 0.5 .* (Vector(panel.LE_point_1) .+ Vector(panel.LE_point_2))
        section = [(dot(p .- le_mid, chord_axis), dot(p .- le_mid, normal_axis))
                   for p in pos]
        winding = contour_winding(section)
        net = zeros(SimFloat, 3)
        for k in 1:n_nodes
            column += 1
            edge = pos[mod1(k + 1, n_nodes)] - pos[mod1(k - 1, n_nodes)]
            edge_len = norm(edge)
            if edge_len < eps()
                @views mode.traction[:, column] .= 0.0
                continue
            end
            tangent = edge ./ edge_len
            surface_dir = dot(tangent, chord_axis) >= 0 ? tangent : -tangent
            normal = normalize(winding .*
                (dot(tangent, normal_axis) .* chord_axis .-
                 dot(tangent, chord_axis) .* normal_axis))
            segment_area = 0.5 * edge_len * panel.width
            traction = (-cp[k] .* normal .+ cf[k] .* surface_dir) .*
                       (dynamic_pressure * segment_area)
            @views mode.traction[:, column] .= traction
            net .+= traction
        end
        @views mode.traction_net[:, panel_idx] .= net
    end
    accumulate_point_offset!(mode)
    return nothing
end

"""
    surface_node_forces(traction, first_column, n_nodes, panel_force, panel_net)

Force on each contour node of one panel: its frozen traction plus an equal share of
`panel_force - panel_net`. Sums to `panel_force` over the panel, the frozen part
cancelling. `traction` is column-indexed panel-major as
[`freeze_traction_pattern!`](@ref) fills it, `first_column` being this panel's
first. Works on the symbolic params and on the numeric buffers alike, so the
scatter has one definition.
"""
function surface_node_forces(traction, first_column, n_nodes, panel_force, panel_net)
    residual = (collect(panel_force) .- collect(panel_net)) ./ n_nodes
    return [[traction[c, first_column + node - 1] for c in 1:3] .+ residual
            for node in 1:n_nodes]
end

"""
    accumulate_point_offset!(mode)

Reduce the frozen pattern to one constant force per wing point, keyed by global point
index: the tractions of the nodes it owns, less the share of each panel's frozen net
that [`aero_scatter_entries`](@ref) re-adds through the live panel force. The scatter
is then a pure weighted sum of panel forces, which is what the wiring layer carries.
Derived wholly from `traction`, `traction_net` and `station_point`, and refreshed with
them.
"""
function accumulate_point_offset!(mode::AeroPressure)
    for total in values(mode.point_offset)
        fill!(total, 0.0)
    end
    column = 0
    for (panel, assigned) in enumerate(mode.station_point)
        share = 1 / length(assigned)
        for point_idx in assigned
            column += 1
            total = get!(zero_aero_force, mode.point_offset, point_idx)
            @views total .+= mode.traction[:, column] .-
                             share .* mode.traction_net[:, panel]
        end
    end
    return nothing
end

"""A fresh zero body-frame force, the accumulator [`accumulate_point_offset!`](@ref)
starts each point from."""
zero_aero_force() = zeros(SimFloat, 3)

# ==================== per-panel decomposition ==================== #

supports_panel_decomposition(::AeroPressure) = true

"""
    aero_scatter_entries(mode::AeroPressure, wing, points)

The surface pattern: a point takes an equal share of its panel's live force for each
contour node of that panel assigned to it. The frozen traction itself is constant and
lives in [`aero_point_offset`](@ref), so only the share remains here.
"""
function aero_scatter_entries(mode::AeroPressure, wing, points)
    point_num = Dict(point.idx => k for (k, point) in enumerate(points))
    totals = Dict{Tuple{Int, Int}, Vector{SimFloat}}()
    for (panel, assigned) in enumerate(mode.station_point)
        share = SimFloat(1 / length(assigned))
        for point_idx in assigned
            scatter_totals!(totals, panel, point_num[point_idx], share, 0.0)
        end
    end
    return scatter_entry_list(totals)
end

aero_point_offset(::AeroPressure, params, wing_idx, point_idx) =
    params.wings[wing_idx].aero.points[point_idx].offset

write_aero_log_points!(mode::AeroPressure, wing, sys_struct, sys_state,
                       point_idx, zoom) =
    write_live_aero_log_points!(mode, wing, sys_struct, sys_state,
                                point_idx, zoom)
