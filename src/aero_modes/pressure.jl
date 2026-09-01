# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# AeroPressure: distribute VSM's per-section force onto arbitrary structural points
# using the airfoil surface traction (−Cp·n̂ + cf·ŝ) as the placement pattern. The
# per-panel total force is re-derived symbolically each RHS step (live), like
# ContinuousAero; only the scatter differs (surface pattern vs strut couple).

"""
    AeroPressure(; frame_tol_frac=2.0, live_polars=false)

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

`live_polars` replaces the tabulated `(α, δ)` polars with polars regenerated from the
deformed shape every solve (see [`solve_with_live_polars!`](@ref)): each panel's
chordwise deformation deforms its Kulfan fit, NeuralFoil is evaluated on a grid of
angles about the panel's own angle of attack, and those values become its rewritten
polar. The flap angle δ then carries no information and is dropped from the equations
entirely, so the RHS has no live deflection left in it — only α stays live. Use it
where the chord bends into a shape a single hinge angle cannot stand for.

A live polar spans only the angles it was sampled over and is held flat past them, so
bring-up still wants damping enough to keep the solve inside that range — releasing the
SK100 from its placed geometry at `start_world_damping = 20` accelerates the wing to
20 m/s in 10 ms and slews α some 20° inside one 0.05 s step, where `300` settles. Past
the range the answer is bounded rather than lost, which is the difference from a
tabulated polar spanning the whole range: stale, not wrong-signed.
"""
mutable struct AeroPressure{E} <: AbstractVSMAero
    engine::Union{Nothing, E}
    "Near-station wing-node index for every (panel, contour node): `station_point[panel][node]`."
    station_point::Vector{Vector{Int64}}
    "Far-station wing-node index for the same (panel, contour node)."
    blend_point::Vector{Vector{Int64}}
    "The two stations every panel lies between and its share of the second."
    panel_blend::Vector{Tuple{Int64, Int64, SimFloat}}
    "Pure-moment chordwise share of every (panel, contour node), see [`couple_shape`](@ref)."
    node_couple_shape::Vector{Vector{SimFloat}}
    "Residual share of every (panel, contour node), see [`node_residual_shares`](@ref)."
    node_residual_share::Vector{Vector{SimFloat}}
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
    "Frozen moment of the traction pattern about each panel's leading-edge midpoint (3 × n_panels)."
    traction_moment::Matrix{SimFloat}
    "Frozen share-weighted node offset from each panel's leading-edge midpoint (3 × n_panels)."
    residual_arm::Matrix{SimFloat}
    "Polar callables `(panel_idx, α[, δ])` for cl/cd/cm, read as callable flat params."
    cl::Any
    cd::Any
    cm::Any
    "Per-panel flap station index (global), or 0 for no flap. Structural (enters the cache hash)."
    panel_station::Vector{Int64}
    "Left unrefined strut of each refined section (n_panels + 1)."
    section_left_strut::Vector{Int64}
    "Left-strut weight: section = w·strut[left] + (1−w)·strut[left+1]."
    section_left_weight::Vector{SimFloat}
    "Frozen body-frame LE billow offset off the strut line (3 × n_sections)."
    section_le_offset::Matrix{SimFloat}
    "Frozen body-frame TE billow offset off the strut line (3 × n_sections)."
    section_te_offset::Matrix{SimFloat}
    "Regenerate the polars from the deformed shape each solve instead of reading δ tables."
    live_polars::Bool
    "Live polar source and control-point map, `nothing` unless `live_polars`."
    live::Any
    AeroPressure{E}(engine, station_point, blend_point, panel_blend,
        node_couple_shape, node_residual_share,
        frame_tol_frac, v_ind, chord_weight, traction, point_offset, traction_net,
        traction_moment, residual_arm,
        cl, cd, cm, panel_station, section_left_strut, section_left_weight,
        section_le_offset, section_te_offset, live_polars, live) where {E} =
        new{E}(engine, station_point, blend_point, panel_blend,
               node_couple_shape, node_residual_share,
               frame_tol_frac, v_ind, chord_weight, traction, point_offset,
               traction_net, traction_moment, residual_arm,
               cl, cd, cm, panel_station, section_left_strut,
               section_left_weight, section_le_offset, section_te_offset,
               live_polars, live)
end

AeroPressure(; frame_tol_frac=2.0, live_polars=false) =
    AeroPressure{VSMEngine}(nothing, Vector{Vector{Int64}}(),
        Vector{Vector{Int64}}(), Tuple{Int64, Int64, SimFloat}[],
        Vector{Vector{SimFloat}}(), Vector{Vector{SimFloat}}(),
        SimFloat(frame_tol_frac),
        zeros(SimFloat, 3, 0), SimFloat[], zeros(SimFloat, 3, 0),
        Dict{Int64, Vector{SimFloat}}(), zeros(SimFloat, 3, 0),
        zeros(SimFloat, 3, 0), zeros(SimFloat, 3, 0),
        nothing, nothing, nothing, Int64[],
        Int64[], SimFloat[], zeros(SimFloat, 3, 0), zeros(SimFloat, 3, 0),
        live_polars, nothing)
attach_engine!(mode::AeroPressure, engine::VSMEngine) =
    AeroPressure{typeof(engine)}(engine, mode.station_point, mode.blend_point,
        mode.panel_blend,
        mode.node_couple_shape, mode.node_residual_share, mode.frame_tol_frac,
        mode.v_ind,
        mode.chord_weight, mode.traction, mode.point_offset, mode.traction_net,
        mode.traction_moment, mode.residual_arm, mode.cl, mode.cd, mode.cm,
        mode.panel_station, mode.section_left_strut,
        mode.section_left_weight, mode.section_le_offset, mode.section_te_offset,
        mode.live_polars, mode.live)

is_builtin_aero(::AeroPressure) = true
aero_mode_tag(::AeroPressure) = "press"

"""
    aero_hash_id(mode::AeroPressure)

The surface-node→point map is baked into the generated scatter equations, so it is
structural and enters the model-cache hash (distinguishing it from any stale
frozen-force build that used the default empty id). The panel→flap-station
map is likewise baked into the δ wiring.
"""
aero_hash_id(mode::AeroPressure) = (mode.station_point, mode.blend_point,
    [round(w; digits=8) for (_, _, w) in mode.panel_blend], mode.panel_station,
    mode.live_polars,
    mode.section_left_strut, round.(mode.section_left_weight; digits=8),
    round.(mode.section_le_offset; digits=8),
    round.(mode.section_te_offset; digits=8))

# ==================== equation builder ==================== #

"""
    aero_component(mode::AeroPressure, wing::ParticleWing, sys_struct; name, params)

Live per-refined-panel VSM force (shared [`build_panel_force_eqs`](@ref) on the
fixed loaded mesh + live apparent wind) scattered over the frozen surface traction
pattern. Each wing node's force is `Σ over mapped (panel, node) of
traction[node] + share[node]·(panel_force[panel] − traction_net[panel])`, its share
from [`node_residual_shares`](@ref), so per panel the point forces sum to the live
`panel_force` exactly. Apparent wind and density are
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
    traction_moment_p = params.wings[wing_idx].aero.traction_moment
    residual_arm_p = params.wings[wing_idx].aero.residual_arm

    points = wing_points(sys_struct, wing)
    num_points = length(points)
    panels = wing.vsm_aero.panels
    n_panels = length(panels)
    length(mode.station_point) == n_panels || error(
        "AeroPressure wing $(wing.name): surface→point map not built " *
        "($(length(mode.station_point)) entries for $n_panels panels).")

    # Add a live per-panel flap deflection connector only when this wing has a flap.
    has_flap_coupling = length(mode.panel_station) == n_panels &&
        any(!=(0), mode.panel_station)
    connectors = particle_aero_connectors(num_points;
        n_delta=has_flap_coupling ? n_panels : 0)

    # Live geometry: refined-section LE/TE from the structural stations' pos_w.
    column = aero_section_columns(wing, points)
    sec_le, sec_te =
        reconstruct_sections_sym(mode, wing, points, connectors, column)

    sec_va, sec_rho, sec_dva =
        reconstruct_inflow_sym(mode, wing, connectors, column)

    spanwise = collect(SimFloat, wing.vsm_wing.spanwise_direction)
    scale = 1.0 + (isfinite(wing.aero_scale_chord) ?
        wing.aero_scale_chord : AERO_SCALE_CHORD)

    orient = panel_span_signs(wing, spanwise)
    delta = has_flap_coupling ? collect(connectors.delta) : nothing
    wagner_eqs, wagner_vars, deficiency, wagner_defaults =
        wagner_wing_eqs(wing, sec_va, params)
    eqs, panel_vars, panel_force, panel_couple, curvature_couple, slots =
        build_panel_force_eqs(sec_le, sec_te, sec_va, sec_rho, vind_p, chord_w,
            cl, cd, cm, spanwise, scale, orient; delta, sec_dva, deficiency)
    append!(eqs, wagner_eqs)
    column_of(p, i) = [p[c, i] for c in 1:3]
    couple = [pressure_couple(slots.panel_couple[:, i], slots.panel_force[:, i],
                  column_of(traction_net_p, i), column_of(traction_moment_p, i),
                  column_of(residual_arm_p, i), slots.x_airf[:, i],
                  slots.y_airf[:, i], slots.z_airf[:, i], slots.chord[i])
              for i in 1:n_panels]
    vars = particle_unknowns(connectors)
    append!(vars, panel_vars)
    append!(vars, wagner_vars)

    point_num = Dict(point.idx => k for (k, point) in enumerate(points))
    point_force = [zeros(Num, 3) for _ in 1:num_points]
    column = 0
    for i in 1:n_panels
        assigned = mode.station_point[i]
        node_forces = surface_node_forces(traction_p, column + 1, length(assigned),
            panel_force[:, i], column_of(traction_net_p, i),
            couple[i], mode.node_couple_shape[i], mode.node_residual_share[i])
        share_far = mode.panel_blend[i][3]
        for node in eachindex(assigned)
            column += 1
            k = point_num[assigned[node]]
            point_force[k] = point_force[k] .+ (1 - share_far) .* node_forces[node]
            if share_far > 0
                k_far = point_num[mode.blend_point[i][node]]
                point_force[k_far] = point_force[k_far] .+ share_far .* node_forces[node]
            end
        end
    end

    for k in 1:num_points
        eqs = [eqs; connectors.point_force[:, k] ~ point_force[k]]
    end
    return System(eqs, t, vars,
        [Any[vind_p, chord_w, cl, cd, cm, traction_p, traction_net_p,
             traction_moment_p, residual_arm_p];
         wagner_params(wing, params)];
        name, initial_conditions=wagner_defaults)
end

# ==================== build-time panel→flap map + live δ ==================== #

"""
    build_panel_station_map!(mode, wing, sys_struct)

Map each refined VSM panel to the flap `KINEMATIC` station of `wing` whose
spanwise station (midpoint of its two flap bodies) is nearest the panel's — the
nearest-station philosophy of [`build_station_point_map!`](@ref). Stored in
`mode.panel_station` (global station idx, `0` = no flap); structural,
so it enters [`aero_hash_id`](@ref). All-zeros (no coupling) when the wing has no
flap surface. Generic modes are a no-op.
"""
build_panel_station_map!(::AbstractAeroModel, wing, sys_struct) = nothing
function build_panel_station_map!(mode::AeroPressure, wing, sys_struct)
    panels = wing.vsm_aero.panels
    n_panels = length(panels)
    if mode.live_polars
        mode.panel_station = zeros(Int64, n_panels)
        return nothing
    end
    flaps = [ts for ts in sys_struct.stations
             if ts.type == KINEMATIC && has_flap(ts) && ts.wing_idx == wing.idx]
    if isempty(flaps)
        mode.panel_station = zeros(Int64, n_panels)
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
    mode.panel_station = assignment
    return nothing
end

"""
    station_deltas(sys_struct) -> Vector{SimFloat}

Live flap deflection δ [rad] per station (indexed by `station.idx`;
`0` for non-flap surfaces), from the current flap-body orientations via
[`flap_delta`](@ref). The Julia ground truth mirroring the symbolic
`station_delta_eqs!`, used to drive `panel.delta` at refresh and by tests.
"""
function station_deltas(sys_struct)
    stations = sys_struct.stations
    bodies = sys_struct.bodies
    deltas = zeros(SimFloat, length(stations))
    for station in stations
        has_flap(station) || continue
        R_main = quaternion_to_rotation_matrix(
            bodies[station.flap_body_idxs[1]].Q_b_to_w)
        R_flap = quaternion_to_rotation_matrix(
            bodies[station.flap_body_idxs[2]].Q_b_to_w)
        deltas[station.idx] = flap_delta(station, R_main, R_flap)
    end
    return deltas
end

"""
    set_panel_deltas!(mode, wing, deltas) -> Bool

Put one flap deflection [rad] per panel onto `wing`'s VSM mesh, `deltas` indexed by
station. `false` when the wing carries no panel-to-surface map to write
through.

Both the panels and the wing's `delta_dist` are written, and the wing's is the one
that matters: a panel's `delta` is derived, and `VortexStepMethod.reinit!` re-seeds
every panel from `delta_dist` whenever the mesh is rebuilt — which
[`refresh_particle_aero!`](@ref) does from the deformed structure before each
solve. Writing the panels alone therefore lasts until the next refresh and no
longer, leaving the solve and the traction contour at δ = 0 on a wing whose flaps
are deflected.
"""
function set_panel_deltas!(mode, wing, deltas)
    panels = wing.vsm_aero.panels
    surfaces = mode.panel_station
    (length(surfaces) == length(panels) && any(!=(0), surfaces)) || return false
    distribution = wing.vsm_wing.delta_dist
    for (panel_idx, panel) in enumerate(panels)
        ts_idx = surfaces[panel_idx]
        (ts_idx == 0 || ts_idx > length(deltas)) && continue
        panel.delta = deltas[ts_idx]
        panel_idx <= length(distribution) &&
            (distribution[panel_idx] = deltas[ts_idx])
    end
    return true
end

"""
    apply_flap_delta!(mode, wing, sys_struct)

Set the VSM mesh's flap deflection from the live deflection of each mapped
station ([`station_deltas`](@ref)), before the VSM solve, so the
converged forces and the frozen traction contour track the flap at refresh
cadence. No-op for modes/wings without flap coupling.
"""
apply_flap_delta!(::AbstractAeroModel, wing, sys_struct) = nothing
function apply_flap_delta!(mode::AeroPressure, wing, sys_struct)
    set_panel_deltas!(mode, wing, station_deltas(sys_struct))
    return nothing
end

function restore_flap_delta!(mode::AeroPressure, wing, sys_state)
    isempty(sys_state.flap_angle) && return nothing
    wing.dynamics_type == PARTICLE_DYNAMICS || return nothing
    set_panel_deltas!(mode, wing, sys_state.flap_angle)
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
    couple_shape(chord_fracs) -> Vector{SimFloat}

Chordwise share of a pure pitching couple over a panel's surface contour nodes, from
their chord fractions. Weights sum to zero and their first moment is `-1` in chords,
so `Σ shape[k] · panel_couple` places no net force and exactly the moment
[`panel_force_eqs`](@ref) put in `panel_couple` — the same normalisation
[`ContinuousAero`](@ref) gets from its `±couple` pair at the LE and TE.

The shape is thin-airfoil theory's `sin 2θ` loading (`θ` from `x/c = (1-cos θ)/2`),
the one mode of the chordwise distribution that carries a moment and no lift. That is
what the increment needs: flow curvature is a parabolic-camber `A₁` term whose lift is
already in the panel force through the three-quarter-chord inflow
([`COLLOCATION_CHORD_FRAC`](@ref)), leaving only its moment to place. Fitting the
constraints within `α·sin 2θ + β` keeps the placement in that family; `β` is nonzero
only because the contour's nodes are not evenly spaced.
"""
function couple_shape(chord_fracs)
    frac = clamp.(SimFloat.(chord_fracs), 0.0, 1.0)
    basis = [sin(2 * acos(1 - 2 * xi)) for xi in frac]
    n = length(frac)
    # α·Σbasis + β·n = 0 and α·Σ(frac·basis) + β·Σfrac = -1.
    det = sum(basis) * sum(frac) - n * sum(frac .* basis)
    abs(det) > 1e-9 * max(n, 1) || error(
        "AeroPressure: the surface contour gives a degenerate pure-moment shape " *
        "(det=$det over $n nodes); its chord fractions carry no first moment.")
    alpha = n / det
    beta = -alpha * sum(basis) / n
    return alpha .* basis .+ beta
end

"""
    node_residual_shares(xc, yc) -> Vector{SimFloat}

Share of its panel's residual every contour node carries, summing to one over the
panel: the node's area times the chordwise shape an added section force takes.

The residual is the panel force the VSM converged on less what the frozen `Cp`
pattern integrates to, and it has to go somewhere. An equal split per node is the
one choice that cannot be right, because it is a property of the mesh rather than
of the flow: the contour is clustered where the airfoil turns, so the trailing edge
owns two fifths of a section's nodes while the Kutta condition leaves it almost
none of the load, and an equal split hands it two fifths of the correction — a
download on the trailing edge of a lifting wing, which folds the flap the way the
flap is already going.

Area alone is mesh-independent but still spreads the correction evenly along the
chord. What an *added* section force actually looks like is thin-airfoil theory's
additional load, `sqrt((1-ξ)/ξ)`: largest at the leading edge, zero at the trailing
edge. That puts the correction where a two-dimensional pattern goes wrong — the
nose — and leaves the trailing edge at the zero load its sharp edge demands.
"""
function node_residual_shares(xc, yc)
    n = length(xc)
    n > 0 || return SimFloat[]
    frac = clamp.(SimFloat.(xc), 0.0, 1.0)
    span(k) = 0.5 * hypot(xc[mod1(k + 1, n)] - xc[mod1(k - 1, n)],
                          yc[mod1(k + 1, n)] - yc[mod1(k - 1, n)])
    weights = [span(k) * sqrt((1 - frac[k]) / max(frac[k], 1.0e-3)) for k in 1:n]
    total = sum(weights)
    return total > 0 ? SimFloat.(weights ./ total) : fill(SimFloat(1 / n), n)
end

"""
    panel_station_candidates(mode, panels, stations, points, wing,
                             point_idx, point_pos_b)
        -> Vector{Vector{Tuple{Int64, KVec3}}}

The two spanwise stations each panel's load is shared between, and the far one's
share: the station it sits on, the next one out, and how far between them it lies.

Which strut a panel sits on comes from the refined-section interpolation the loft
carries, not from measuring the panel's position, so it holds for a swept or dihedral
wing as well as a flat one. A panel's load is split between its two neighbouring
stations rather than rounded onto the nearer, because rounding is a discontinuity: two
panels that mirror each other land on centres that agree to the last bit but not
exactly, and the nearer station is then a different one for each, moving a whole
panel's load a station across the span. Sharing it moves that disagreement back into
the weights, where it stays the size it actually is. A panel is a slice of one airfoil, and it has to stay on one
strut. Left to search the
whole wing, the nodes near a station boundary defect to the neighbouring strut, so half
an airfoil hangs off one station and half off the next — the panel then spans two
spanwise planes, and a chordwise profile read off it has a step in it that no airfoil
basis can represent. Choosing the station first and mapping within it keeps each panel
whole. Returns `nothing` when the wing declares no stations, which leaves
the caller on the plain nearest-point search: matching by chord fraction alone is only
safe once the candidates are confined to one station, since the same fraction occurs at
every station across the span.
"""
function panel_station_candidates(mode, panels, stations, points, wing,
                                  point_idx, point_pos_b)
    control = station_control_points(stations, points, wing)
    length(control) >= 2 || return nothing
    lookup = Dict(zip(point_idx, point_pos_b))
    grouped = [[(idx, lookup[idx]) for idx in group if haskey(lookup, idx)]
               for group in control]
    any(isempty, grouped) && return nothing
    return [(grouped[lo], grouped[hi], lo == hi ? 0.0 : weight)
            for (lo, hi, weight) in mode.panel_blend]
end

"""
    build_station_point_map!(mode::AeroPressure, wing, points; prn=false)

Build the static map from each refined panel's surface contour nodes to their
nearest wing-node structural point (body frame), and apply the frame-alignment guard.
Stored in `mode.station_point`; errors on a missing `section_aero` or a mesh that
sits farther than `frame_tol_frac` chords from the points.
"""
function build_station_point_map!(mode::AeroPressure, wing, points, stations;
                                  prn=false)
    wing_pts = [p for p in points if p.is_wing_node && p.wing_idx == wing.idx]
    isempty(wing_pts) && error(
        "AeroPressure wing $(wing.name): no wing nodes to receive forces.")
    rot_cad_to_body = wing.R_b_to_c'
    origin_cad = wing.pos_cad
    point_idx = [p.idx for p in wing_pts]
    point_pos_b = [rot_cad_to_body * (p.pos_cad - origin_cad) for p in wing_pts]

    panels = wing.vsm_aero.panels
    control = station_control_points(stations, points, wing)
    mode.panel_blend = length(control) >= 2 ?
        panel_strut_blend(mode, length(control), length(panels)) :
        [(1, 1, 0.0) for _ in eachindex(panels)]
    candidates = panel_station_candidates(mode, panels, stations, points,
                                          wing, point_idx, point_pos_b)
    station_point = Vector{Vector{Int64}}(undef, length(panels))
    blend_point = Vector{Vector{Int64}}(undef, length(panels))
    node_couple_shape = Vector{Vector{SimFloat}}(undef, length(panels))
    node_residual_share = Vector{Vector{SimFloat}}(undef, length(panels))
    max_chord_ratio = 0.0
    for (panel_idx_local, panel) in enumerate(panels)
        panel.section_aero === nothing && error(
            "AeroPressure wing $(wing.name): panel $panel_idx_local has no " *
            "section_aero (surface contour). Load geometry with .dat/Cp/cf tables.")
        xc, yc, _, _ = VortexStepMethod.section_surface(
            panel.section_aero, 0.0, panel.delta)
        assigned = Vector{Int64}(undef, length(xc))
        shared = Vector{Int64}(undef, length(xc))
        near, far, weight = isnothing(candidates) ? (nothing, nothing, 0.0) :
            candidates[panel_idx_local]
        fractions(group) = isnothing(group) ? nothing :
            [chord_frame_coordinates(panel, pos_b)[1] for (_, pos_b) in group]
        near_fraction, far_fraction = fractions(near), fractions(far)
        for k in eachindex(xc)
            node = loft_contour_node(panel, xc, yc, k)
            if isnothing(near)
                min_dist = Inf
                best = point_idx[1]
                for (m, pos_b) in enumerate(point_pos_b)
                    dist = norm(node - pos_b)
                    dist < min_dist && (min_dist = dist; best = point_idx[m])
                end
                assigned[k] = best
                shared[k] = best
                max_chord_ratio = max(max_chord_ratio,
                                      min_dist / max(panel.chord, eps()))
            else
                pick = argmin(abs.(near_fraction .- xc[k]))
                far_pick = argmin(abs.(far_fraction .- xc[k]))
                assigned[k] = near[pick][1]
                shared[k] = far[far_pick][1]
                # The guard asks whether the mesh and the structure line up, so it
                # measures the station actually carrying the panel, not the one that
                # happens to come first.
                held = weight < 0.5 ? near[pick][2] : far[far_pick][2]
                max_chord_ratio = max(max_chord_ratio,
                    norm(node - held) / max(panel.chord, eps()))
            end
        end
        station_point[panel_idx_local] = assigned
        blend_point[panel_idx_local] = shared
        node_couple_shape[panel_idx_local] = couple_shape(xc)
        node_residual_share[panel_idx_local] = node_residual_shares(xc, yc)
    end
    max_chord_ratio > mode.frame_tol_frac && error(
        "AeroPressure wing $(wing.name): max surface-node→point distance is " *
        "$(round(max_chord_ratio, digits=2))× the local chord, exceeding " *
        "frame_tol_frac=$(mode.frame_tol_frac). The VSM mesh and structural " *
        "points are likely misaligned (frame mismatch).")
    mode.station_point = station_point
    mode.blend_point = blend_point
    mode.node_couple_shape = node_couple_shape
    mode.node_residual_share = node_residual_share
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
    size(mode.traction_moment) == (3, n_panels) ||
        (mode.traction_moment = zeros(SimFloat, 3, n_panels))
    size(mode.residual_arm) == (3, n_panels) ||
        (mode.residual_arm = zeros(SimFloat, 3, n_panels))
    mode.cl = ContinuousPolar(body_aero, VortexStepMethod.calculate_cl)
    mode.cd = ContinuousPolar(body_aero, VortexStepMethod.calculate_cd)
    mode.cm = ContinuousPolar(body_aero, VortexStepMethod.calculate_cm)
    return nothing
end

"""
    setup_aero!(mode::AeroPressure, wing, points, stations; prn=false)

Run the generic particle VSM setup (rebuild the unrefined sections onto the
structural LE/TE stations, refine, build `point_to_vsm_point`), freeze the
strut-interpolation caches ([`build_section_interp`](@ref)), then build the
surface→point traction map ([`build_station_point_map!`](@ref)) and size the frozen
buffers/polars ([`init_pressure_buffers!`](@ref)). `PARTICLE_DYNAMICS` only.
"""
function setup_aero!(mode::AeroPressure, wing, points, stations; prn=false)
    wing.dynamics_type == PARTICLE_DYNAMICS || error(
        "AeroPressure supports PARTICLE_DYNAMICS wings only; wing " *
        "$(wing.name) is $(wing.dynamics_type).")
    invoke(setup_aero!, Tuple{AbstractVSMAero, Any, Any, Any},
           mode, wing, points, stations; prn)
    mode.section_left_strut, mode.section_left_weight,
        mode.section_le_offset, mode.section_te_offset =
        build_section_interp(wing.vsm_wing)
    build_station_point_map!(mode, wing, points, stations; prn)
    init_pressure_buffers!(mode, wing)
    mode.live_polars && build_live_polars!(mode, wing, points, stations)
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
    remake_aero!(mode::AeroPressure, wing, set, vsm_set, points, stations)

Rebuild the VSM objects from edited settings via the generic particle remake
(coarsen onto the structural stations, refine, rebuild `point_to_vsm_point`), then
rebuild the strut-interpolation caches, surface→point map and frozen buffers/polars.
"""
function remake_aero!(mode::AeroPressure, wing, set, vsm_set, points,
                      stations)
    invoke(remake_aero!, Tuple{AbstractVSMAero, Any, Any, Any, Any, Any},
           mode, wing, set, vsm_set, points, stations)
    mode.section_left_strut, mode.section_left_weight,
        mode.section_le_offset, mode.section_te_offset =
        build_section_interp(wing.vsm_wing)
    build_station_point_map!(mode, wing, points, stations)
    init_pressure_buffers!(mode, wing)
    mode.live_polars && build_live_polars!(mode, wing, points, stations)
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
        fill!(mode.traction_moment, 0.0)
        fill!(mode.residual_arm, 0.0)
        accumulate_point_offset!(mode)
        return nothing
    end
    update_vsm_wing_from_structure!(wing, points)
    set_refined_panel_va!(mode, wing, points, va_point_b_vals)
    if mode.live_polars
        solve_with_live_polars!(mode, wing, points; cold_start)
    else
        solve_and_freeze_circulation!(mode, wing; cold_start)
    end
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
loft the airfoil surface contour and take its traction pattern from
[`surface_pattern`](@ref) — the section's own tables at the effective sectional α
(`sol.alpha_dist`), or, on live polars, the deformed shape's — build
the per-segment traction `(−Cp·n̂ + cf·ŝ)·q·dA` (`n̂` outward, `ŝ` along the chord) into
`mode.traction` (panel-major node order, matching [`aero_component`](@ref)), and store
the per-panel net in `mode.traction_net`. The net anchors the symbolic scatter so each
panel's point forces sum to the live VSM total.

`mode.traction_moment` and `mode.residual_arm` are the moment the frozen pattern
already carries about the panel's leading-edge midpoint and the share-weighted lever
arm the residual will be spread on; [`pressure_couple`](@ref) subtracts both from the
polar's moment so the couple it places is only what is missing.
"""
function freeze_traction_pattern!(mode::AeroPressure, wing)
    sol = wing.vsm_solver.sol
    panels = wing.vsm_aero.panels
    dynamic_pressure = 0.5 * wing.vsm_solver.density * dot(wing.va_b, wing.va_b)
    column = 0
    for (panel_idx, panel) in enumerate(panels)
        xc, yc, cp, cf = surface_pattern(mode, panel, panel_idx,
                                         sol.alpha_dist[panel_idx])
        n_nodes = length(xc)
        chord_axis = Vector(panel.x_airf)
        normal_axis = Vector(panel.z_airf)
        pos = [loft_contour_node(panel, xc, yc, k) for k in 1:n_nodes]
        le_mid = 0.5 .* (Vector(panel.LE_point_1) .+ Vector(panel.LE_point_2))
        section = [(dot(p .- le_mid, chord_axis), dot(p .- le_mid, normal_axis))
                   for p in pos]
        winding = contour_winding(section)
        share = mode.node_residual_share[panel_idx]
        net = zeros(SimFloat, 3)
        moment = zeros(SimFloat, 3)
        arm = zeros(SimFloat, 3)
        for k in 1:n_nodes
            column += 1
            edge = pos[mod1(k + 1, n_nodes)] - pos[mod1(k - 1, n_nodes)]
            edge_len = norm(edge)
            arm .+= share[k] .* (pos[k] .- le_mid)
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
            moment .+= cross(pos[k] .- le_mid, traction)
        end
        @views mode.traction_net[:, panel_idx] .= net
        @views mode.traction_moment[:, panel_idx] .= moment
        @views mode.residual_arm[:, panel_idx] .= arm
    end
    accumulate_point_offset!(mode)
    return nothing
end

"""
    surface_node_forces(traction, first_column, n_nodes, panel_force, panel_net,
                        couple, shape)

Force on each contour node of one panel: its frozen traction, its
[`node_residual_shares`](@ref) share of `panel_force - panel_net`, and
`shape[node] · couple`. Sums to `panel_force` over the panel — the frozen part
cancels, `share` sums to one and [`couple_shape`](@ref) sums to zero — while the
couple adds the pitching moment the shape is normalised to. `traction` is
column-indexed panel-major as [`freeze_traction_pattern!`](@ref) fills it,
`first_column` being this panel's first. Works on the symbolic params and on the
numeric buffers alike, so the scatter has one definition.
"""
function surface_node_forces(traction, first_column, n_nodes, panel_force, panel_net,
                             couple, shape, share)
    residual = collect(panel_force) .- collect(panel_net)
    moment = collect(couple)
    return [[traction[c, first_column + node - 1] for c in 1:3] .+
            share[node] .* residual .+ shape[node] .* moment
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
        share = mode.node_residual_share[panel]
        far = mode.panel_blend[panel][3]
        for (node, point_idx) in enumerate(assigned)
            column += 1
            for (idx, w) in ((point_idx, 1 - far),
                             (mode.blend_point[panel][node], far))
                w > 0 || continue
                total = get!(zero_aero_force, mode.point_offset, idx)
                @views total .+= w .* (mode.traction[:, column] .-
                                       share[node] .* mode.traction_net[:, panel])
            end
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
    pressure_couple(panel_couple, panel_force, traction_net, traction_moment,
                    residual_arm, x_airf, y_airf, z_airf, chord) -> Vector

The couple [`couple_shape`](@ref) has to place for the panel's pitching moment to be
the polar's. The target about the panel's leading-edge midpoint is the one
[`ContinuousAero`](@ref) places from its `±couple` pair: `panel_force` acting at the
quarter chord plus `panel_couple` as a pure couple. Against it stand the moment the
frozen traction already carries (`traction_moment`) and the one the residual
`panel_force - traction_net` will carry on its share-weighted arm `residual_arm`. What
is left is the deficit, and dividing by the chord returns it to the
force-along-`z_airf` currency `panel_couple` is written in.

The airfoil axes are not exactly orthogonal on a billowed panel, so the divisor carries
the frame's own triple product rather than assuming one; on the SK100 that alone is
most of a percent.

Without this the panel's pitching moment is whatever its `Cp` pattern integrates to and
the polar's `cm` never reaches the structure — the two agree closely on a live polar,
where both come from one solve of one shape, but nothing was holding them together.
"""
function pressure_couple(panel_couple, panel_force, traction_net, traction_moment,
                         residual_arm, x_airf, y_airf, z_airf, chord)
    force = collect(panel_force)
    x_axis, y_axis, z_axis = collect(x_airf), collect(y_airf), collect(z_airf)
    residual = force .- collect(traction_net)
    present = collect(traction_moment) .+ cross(collect(residual_arm), residual)
    target = (0.25 * chord) .* cross(x_axis, force) .-
             chord .* cross(x_axis, collect(panel_couple))
    gain = chord * dot(cross(z_axis, x_axis), y_axis)
    return ((dot(target .- present, y_axis)) / gain) .* z_axis
end

scatter_couple(::AeroPressure, slots, i, panel) =
    pressure_couple(slots.panel_couple[:, i], slots.panel_force[:, i],
                    panel.traction_net, panel.traction_moment, panel.residual_arm,
                    slots.x_airf[:, i], slots.y_airf[:, i], slots.z_airf[:, i],
                    slots.chord[i])

"""
    aero_scatter_entries(mode::AeroPressure, wing, points)

The surface pattern: a point takes each of its contour nodes'
[`node_residual_shares`](@ref) share of the panel's live force, and that node's
[`couple_shape`](@ref) share of the panel's couple. The frozen traction itself is
constant and lives in [`aero_point_offset`](@ref), so only the shares remain here.
"""
function aero_scatter_entries(mode::AeroPressure, wing, points)
    point_num = Dict(point.idx => k for (k, point) in enumerate(points))
    totals = Dict{Tuple{Int, Int}, Vector{SimFloat}}()
    for (panel, assigned) in enumerate(mode.station_point)
        share = mode.node_residual_share[panel]
        shape = mode.node_couple_shape[panel]
        far = mode.panel_blend[panel][3]
        for (node, point_idx) in enumerate(assigned)
            scatter_totals!(totals, panel, point_num[point_idx],
                            (1 - far) * share[node], (1 - far) * shape[node])
            far > 0 && scatter_totals!(totals, panel,
                point_num[mode.blend_point[panel][node]],
                far * share[node], far * shape[node])
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
