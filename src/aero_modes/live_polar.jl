# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Live polars for AeroPressure: re-generate every panel's 2D polar from the deformed
# structure each VSM solve, instead of reading a table on a flap angle. The chordwise
# deformation reaches the aerodynamics as the shape it is, so nothing has to be
# projected onto a single hinge angle first.

using VortexStepMethod: AirfoilAero

"""
    LivePolarState(source, control_point, control_fraction, control_offset)

Per-wing state of the live polar path: the
[`LivePolars`](@ref VortexStepMethod.AirfoilAero.LivePolars) source that owns the base
airfoils and the fit, plus the control points each panel reads its chordwise
deformation from and where those points sat in the reference geometry.

The control points are the distinct structural points a panel's surface nodes map to
(`AeroPressure.station_point`), so they are exactly the wing nodes that already carry
that panel's load — a chord-line beam's nodes today, a membrane's nodes later. Their
current offset off the deformed chord line, minus the reference offset stored here, is
the camber increment handed to
[`deform_kulfan`](@ref VortexStepMethod.AirfoilAero.deform_kulfan).
"""
mutable struct LivePolarState
    "The polar source: base airfoils, CST basis, reference angles and network scratch."
    source::AirfoilAero.LivePolars
    "Structural point indices of each spanwise station."
    control_point::Vector{Vector{Int64}}
    "Reference chord fraction of every control point, station by station."
    control_fraction::Vector{Vector{SimFloat}}
    "Reference offset off the chord line, over chord, of every control point."
    control_offset::Vector{Vector{SimFloat}}
    "Panel each station reads its chord frame from."
    station_panel::Vector{Int64}
    "Stations bracketing each panel and how far between them it lies."
    panel_blend::Vector{Tuple{Int64, Int64, SimFloat}}
    "Camber increment per station on the CST basis stations, rewritten every refresh."
    station_deflection::Vector{Vector{SimFloat}}
    "Camber increment per panel, interpolated from the stations either side."
    deflection::Vector{Vector{SimFloat}}
    "Chord fraction of every contour node, panel by panel, at δ = 0."
    contour_x::Vector{Vector{SimFloat}}
    "Undeformed normal offset of every contour node, panel by panel, at δ = 0."
    contour_y_ref::Vector{Vector{SimFloat}}
    "Normal offset of every contour node with the deformation added, per solve."
    contour_y::Vector{Vector{SimFloat}}
    "CST shape matrix at each panel's contour fractions, built once."
    contour_shape::Vector{Matrix{SimFloat}}
    "Contour node index of each panel's leading edge, splitting upper from lower."
    leading_edge::Vector{Int64}
    "`Cp` per contour node, regenerated from the deformed shape every solve."
    contour_cp::Vector{Vector{SimFloat}}
    "Skin friction per contour node, regenerated every solve."
    contour_cf::Vector{Vector{SimFloat}}
end

"""
    chord_frame_coordinates(panel, pos_b) -> (fraction, offset)

Where a body-frame position sits in a panel's chord frame: along-chord fraction off the
mid leading edge and offset off the chord line, both over the panel chord. The frame is
built from the panel's four `corner_points`, which is the one part of a panel's geometry
a replayed log frame restores, so a replay measures the deformation of the frame it
draws rather than of whatever the last solve left behind. Its normal is sign-aligned to
`z_airf`, so a control point and the traction pattern cannot disagree about which way is
up.
"""
function chord_frame_coordinates(panel, pos_b)
    corners = panel.corner_points
    le_1 = SVector{3}(@view corners[:, 1])
    te_1 = SVector{3}(@view corners[:, 2])
    te_2 = SVector{3}(@view corners[:, 3])
    le_2 = SVector{3}(@view corners[:, 4])
    le_mid = 0.5 .* (le_1 .+ le_2)
    chord_vec = 0.5 .* (te_1 .+ te_2) .- le_mid
    chord = max(norm(chord_vec), eps())
    x_axis = chord_vec ./ chord
    normal = cross(x_axis, le_1 .- le_2)
    normal_length = norm(normal)
    z_axis = normal_length < eps() ? SVector{3}(panel.z_airf) : normal ./ normal_length
    dot(z_axis, SVector{3}(panel.z_airf)) < 0 && (z_axis = -z_axis)
    offset = SVector{3}(pos_b) .- le_mid
    return (dot(offset, x_axis) / chord, dot(offset, z_axis) / chord)
end

"""
    station_control_points(stations, points, wing) -> Vector{Vector{Int64}}

The wing's control points grouped into the spanwise stations the structure declares.
A station already names the points of one station — its leading edge, its trailing
edge and the control points between them — so the grouping is read rather than inferred
from geometry or from the beam graph, and it holds for a wing with no beam joints at all.

Ownership is taken from the points rather than from the surface's `wing_idx`, which a
geometry that never names a wing leaves at zero. A surface names its leading and
trailing edge as well as the control points running between them, and the first and last
of those often sit on the very same nodes, so coincident points are kept once: they are
one place on the chord, and a fit handed the same station twice has no answer.
"""
function station_control_points(stations, points, wing)
    owned = [surface.point_idxs for surface in stations
             if length(surface.point_idxs) >= 2 &&
                all(points[i].wing_idx == wing.idx for i in surface.point_idxs)]
    return [unique(idx -> round.(points[idx].pos_cad; digits = 9), group)
            for group in owned]
end

"""
    panel_strut_blend(mode::AeroPressure, n_stations, n_panels)
        -> Vector{Tuple{Int64, Int64, SimFloat}}

The two stations each panel lies between and how far it lies towards the second, taken
from the refined-section interpolation the loft carries rather than measured off the
panel's position. Load and deformation both read this one answer: a panel is deformed by
the points it is loaded through, so the two cannot be allowed to disagree. A panel
spans two refined sections, so it sits at the mean of their places.

A section's weight is its share of `strut[left]`, which puts it at `left + (1 - w)` in
strut units. Averaging the weights themselves only holds while both sections name the
same left strut: across a strut the two are shares of different pairs, and their mean
says nothing. A panel straddling a strut then lands on the far one.
"""
function panel_strut_blend(mode::AeroPressure, n_stations, n_panels)
    left, weight = mode.section_left_strut, mode.section_left_weight
    length(left) == n_panels + 1 || error(
        "AeroPressure: the refined-section interpolation has $(length(left)) entries " *
        "for $(n_panels + 1) sections; the loft and the mesh disagree.")
    maximum(left) <= n_stations || error(
        "AeroPressure: the loft names strut $(maximum(left)) but the wing declares " *
        "only $n_stations stations; sections and stations do not correspond.")
    place(i) = left[i] + (1 - weight[i])
    centre = [0.5 * (place(i) + place(i + 1)) for i in 1:n_panels]
    return map(centre) do c
        lo = clamp(floor(Int, c), 1, n_stations)
        (lo, min(lo + 1, n_stations), clamp(c - lo, 0.0, 1.0))
    end
end

"""
    build_live_polars!(mode::AeroPressure, wing, points, stations; settings)

Build the wing's [`LivePolarState`](@ref) from the reference geometry: fit each panel's
undeformed Kulfan parameters off its surface contour, take the spanwise stations the
wing declares, and record where each of their control points sits in that station's
chord frame. Run after [`build_station_point_map!`](@ref), which binds each panel to
one of those same stations, and while the mesh still stands on the reference (CAD)
structure — the offsets stored here are the zero the live deformation is measured
against.

Deformation and load read the same declared stations, so the points a panel is deformed
by are the ones it is loaded through. Measuring per station and interpolating onto the
panels is what keeps the deformation smooth across the span: a panel between two stations
has no structural points of its own, and binding it to a single station would make the
deformation a staircase.
"""
function build_live_polars!(mode::AeroPressure, wing, points, stations;
                            settings=AirfoilAero.LivePolarSettings())
    panels = wing.vsm_aero.panels
    rot_cad_to_body = wing.R_b_to_c'
    origin_cad = wing.pos_cad
    source = AirfoilAero.LivePolars(AirfoilAero.panel_kulfan_parameters(panels);
                                    settings)
    n_panels = length(panels)
    body_of(idx) = rot_cad_to_body * (points[idx].pos_cad - origin_cad)

    control = station_control_points(stations, points, wing)
    n_stations = length(control)
    n_stations >= 2 || error(
        "AeroPressure wing $(wing.name): live polars need at least two spanwise " *
        "stations; this wing declares $n_stations.")
    panel_blend = mode.panel_blend
    length(panel_blend) == n_panels || error(
        "AeroPressure wing $(wing.name): the station blend has $(length(panel_blend)) " *
        "entries for $n_panels panels; run build_station_point_map! first.")
    centre = [lo + weight for (lo, _, weight) in panel_blend]
    station_panel = [argmin(abs.(centre .- s)) for s in 1:n_stations]

    control_fraction = Vector{Vector{SimFloat}}(undef, n_stations)
    control_offset = Vector{Vector{SimFloat}}(undef, n_stations)
    for (s, group) in enumerate(control)
        frames = [chord_frame_coordinates(panels[station_panel[s]], body_of(i))
                  for i in group]
        control_fraction[s] = SimFloat[f[1] for f in frames]
        control_offset[s] = SimFloat[f[2] for f in frames]
    end

    maximum(length, control_fraction) <= 2 && @warn(
        "Live polars see no chordwise deformation: every station names at most two " *
        "points, which are the chord ends the frame is built on. Add chordwise " *
        "control points to the wing, or the polars only track α and Reynolds.",
        wing=wing.name)
    n_basis = length(source.basis.x)
    contour = [VortexStepMethod.section_surface(panel.section_aero, 0.0, 0.0)
               for panel in panels]
    contour_x = [collect(SimFloat, c[1]) for c in contour]
    contour_y = [collect(SimFloat, c[2]) for c in contour]
    node_buffer() = [zeros(SimFloat, length(x)) for x in contour_x]
    mode.live = LivePolarState(source, control, control_fraction, control_offset,
        station_panel, panel_blend,
        [zeros(SimFloat, n_basis) for _ in 1:n_stations],
        [zeros(SimFloat, n_basis) for _ in 1:n_panels],
        contour_x, contour_y, deepcopy(contour_y),
        [AirfoilAero.contour_shape_matrix(x, source.basis.n_weights)
         for x in contour_x],
        [argmin(x) for x in contour_x], node_buffer(), node_buffer())
    return nothing
end

"""
    update_live_deflection!(mode::AeroPressure, wing, points)

Refresh every spanwise station's camber increment from the deformed structure, then
interpolate those increments onto the panels. Each control point contributes its current
offset off its station's deformed chord line minus its reference offset, placed at its
current chord fraction; the chord's own rotation and stretch are already absorbed by the
frame, which the mesh rebuild took from the deformed leading and trailing edges. What
travels spanwise is the increment, not the airfoil, so a panel between two stations keeps
its own base fit. Fills `mode.live.deflection`, ready for [`refit_live_polars!`](@ref).
"""
function update_live_deflection!(mode::AeroPressure, wing, points)
    state = mode.live::LivePolarState
    panels = wing.vsm_aero.panels
    rot_body_to_world = wing.R_b_to_w::Matrix{SimFloat}
    origin = wing.pos_w::KVec3
    for s in eachindex(state.control_point)
        assigned = state.control_point[s]
        panel = panels[state.station_panel[s]]
        fractions = Vector{SimFloat}(undef, length(assigned))
        deflections = Vector{SimFloat}(undef, length(assigned))
        for (k, idx) in enumerate(assigned)
            pos_b = rot_body_to_world' * (points[idx].pos_w - origin)
            fraction, offset = chord_frame_coordinates(panel, pos_b)
            fractions[k] = fraction
            deflections[k] = offset - state.control_offset[s][k]
        end
        state.station_deflection[s] .= AirfoilAero.control_point_deflection(
            state.source.basis, fractions, deflections)
    end
    for i in eachindex(panels)
        lo, hi, blend = state.panel_blend[i]
        @. state.deflection[i] = (1 - blend) * state.station_deflection[lo] +
            blend * state.station_deflection[hi]
    end
    return nothing
end

"""
    panel_reynolds(wing) -> Vector{SimFloat}

Reynolds number of every panel, from the solver's air properties and the panel's own
apparent wind and chord.
"""
panel_reynolds(wing) = [wing.vsm_solver.density * norm(panel.va) * panel.chord /
                        wing.vsm_solver.mu for panel in wing.vsm_aero.panels]

"""
    refit_live_polars!(mode::AeroPressure, wing, alpha) -> Float64

Regenerate every panel's polar at the stored deformation and write it into the
panel's own table about `alpha` [rad] per panel. Reynolds is per panel from the
solver's air properties and the panel's own apparent wind and chord. Returns the
lowest NeuralFoil analysis confidence over the wing; warns once below `0.5`, which
means a deformed section has left the region the network was trained on.
"""
function refit_live_polars!(mode::AeroPressure, wing, alpha)
    state = mode.live::LivePolarState
    confidence = AirfoilAero.refresh_live_polars!(state.source, wing.vsm_aero.panels,
        alpha, panel_reynolds(wing); deflection=state.deflection)
    confidence < 0.5 && @warn("Live polar off the trained shape range",
        wing=wing.idx, confidence, maxlog=1)
    return confidence
end

"""
    live_polar_alpha(wing) -> Vector{SimFloat}

The angles of attack [rad] to sample each panel's polar about: the previous solve's
converged `lr.alpha_dist`, or, before there is one, each panel's geometric angle from
its own apparent wind. Sampling the first solve of a run about zero would hand it
polars that say nothing at the real angle.
"""
function live_polar_alpha(wing)
    alpha = collect(SimFloat, wing.vsm_solver.lr.alpha_dist)
    all(iszero, alpha) || return alpha
    return [SimFloat(VortexStepMethod.calculate_relative_alpha_and_relative_velocity(
                panel, zeros(SimFloat, 3))[1])
            for panel in wing.vsm_aero.panels]
end

"""
    check_live_polar(mode::AeroPressure, wing; panel_idx=nothing) -> NamedTuple

XFoil against the live polar one panel is flying, at the state the model is in now: the
panel's own deformed shape, its converged angle of attack and its Reynolds. `panel_idx`
defaults to the panel NeuralFoil is least confident about, which is the one worth
asking about. See `AirfoilAero.compare_live_polar` for what comes back.

The `Live polar off the trained shape range` warning says the network is extrapolating.
This says whether the extrapolation is any good, which the confidence cannot: it scores
the input shape, not the answer. One viscous solve, so it is cheap enough to call from
the REPL mid-run and far too slow to call every step.
"""
function check_live_polar(mode::AeroPressure, wing; panel_idx=nothing)
    mode.live_polars || error(
        "check_live_polar: wing $(wing.idx) is not flying live polars.")
    state = mode.live::LivePolarState
    idx = isnothing(panel_idx) ? argmin(state.source.confidence) : panel_idx
    return AirfoilAero.compare_live_polar(state.source, wing.vsm_aero.panels[idx], idx,
        wing.vsm_solver.lr.alpha_dist[idx], panel_reynolds(wing)[idx])
end

"""
    refresh_live_pressure!(mode::AeroPressure, wing)

Regenerate the surface traction pattern from the deformed shape at the converged angle
of attack: the contour offset by the deformation's own camber increment, `Cp` from one
batched NeuralFoil pass over every panel, and skin friction from the flat-plate closure
at each panel's own Reynolds.

This is what keeps force and placement telling the same story. The polars already make a
panel's total force follow its deformed shape; without this the pattern spreading that
force over the structure would still be the undeformed section's. Both halves matter:
`Cp` sets how hard each node is pulled, and the contour sets which way — every traction
is `-Cp·n̂` with `n̂` a finite difference of neighbouring node positions, so an
undeformed contour points the whole load the wrong way wherever the section has moved.

Run it after the solve has converged, on the shapes [`refit_live_polars!`](@ref)
evaluated. The nodes move but their assignment to structural points does not: that map
is baked into the generated equations, which is what keeps the scatter continuous when
the section deforms.
"""
function refresh_live_pressure!(mode::AeroPressure, wing)
    state = mode.live::LivePolarState
    reynolds = panel_reynolds(wing)
    AirfoilAero.live_shape_offset!(state.contour_y, state.source, state.contour_shape)
    for i in eachindex(state.contour_y)
        state.contour_y[i] .+= state.contour_y_ref[i]
    end
    AirfoilAero.refresh_live_pressure!(state.contour_cp, state.source,
        state.contour_x, state.contour_y, state.leading_edge,
        collect(SimFloat, wing.vsm_solver.lr.alpha_dist), reynolds)
    AirfoilAero.live_surface_friction!(state.contour_cf, state.contour_x, reynolds)
    return nothing
end

"""
    surface_pattern(mode::AeroPressure, panel, panel_idx, alpha) -> (x, y, cp, cf)

The contour a panel's force is spread over and the traction pattern over it. Tabulated
polars read the section's own `(α, δ)` tables. Live polars read the same contour at
δ = 0 with the contour, `Cp` and skin friction all regenerated from the deformed shape
by [`refresh_live_pressure!`](@ref).
"""
surface_pattern(mode::AeroPressure, panel, panel_idx, alpha) =
    if mode.live_polars
        state = mode.live::LivePolarState
        (state.contour_x[panel_idx], state.contour_y[panel_idx],
         state.contour_cp[panel_idx], state.contour_cf[panel_idx])
    else
        VortexStepMethod.section_surface(panel.section_aero, alpha, panel.delta)
    end

"""
    solve_with_live_polars!(mode::AeroPressure, wing, points; cold_start=false)

Solve the wing with polars regenerated from its current shape, sampled about
[`live_polar_alpha`](@ref). One sampling and one solve: a sampled polar holds its
last value past either end rather than extrapolating, so a solve landing outside the
sampled range reads a bounded answer instead of a runaway one and needs no refit to
be safe. How far it landed is reported by `AirfoilAero.polar_drift`, which warns once
past the range — there the panel's polar is flat and no longer tracks the shape.

The traction pattern is regenerated from the same deformed shapes once the solve has
converged, so the forces and their placement both follow the deformation.
"""
function solve_with_live_polars!(mode::AeroPressure, wing, points; cold_start=false)
    state = mode.live::LivePolarState
    update_live_deflection!(mode, wing, points)
    refit_live_polars!(mode, wing, live_polar_alpha(wing))
    solve_and_freeze_circulation!(mode, wing; cold_start)
    refresh_live_pressure!(mode, wing)
    drift = AirfoilAero.polar_drift(state.source,
        collect(SimFloat, wing.vsm_solver.lr.alpha_dist))
    drift > 1.0 && @warn("Live polars solved past the sampled range; those panels " *
        "are flying the last sampled value", wing=wing.idx,
        drift=round(drift; digits=2), maxlog=1)
    return nothing
end

"""
    restore_live_shape!(mode::AeroPressure, wing, points)

Re-derive every panel's deformed airfoil from a replayed log frame, so the shape a
plot draws is the one that frame's structure gives rather than the one the last live
solve left on the panels. Only the shape is rebuilt: a replayed frame is never solved,
so its polars would say nothing and the network pass they cost would buy nothing.
No-op unless the wing is actually on live polars.
"""
function restore_live_shape!(mode::AeroPressure, wing, points)
    mode.live_polars || return nothing
    state = mode.live
    isnothing(state) && return nothing
    update_live_deflection!(mode, wing, points)
    AirfoilAero.apply_live_shapes!(state.source, wing.vsm_aero.panels;
                                   deflection=state.deflection)
    return nothing
end
