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
    "The polar source: base airfoils, CST basis, expansion points and network scratch."
    source::AirfoilAero.LivePolars
    "Structural point indices each panel reads its deformation from."
    control_point::Vector{Vector{Int64}}
    "Reference chord fraction of every control point, panel by panel."
    control_fraction::Vector{Vector{SimFloat}}
    "Reference offset off the chord line, over chord, of every control point."
    control_offset::Vector{Vector{SimFloat}}
    "Camber increment per panel on the CST basis stations, rewritten every refresh."
    deflection::Vector{Vector{SimFloat}}
end

"""
    chord_frame_coordinates(panel, pos_b) -> (fraction, offset)

Where a body-frame position sits in a panel's chord frame: along-chord fraction off the
mid leading edge and offset off the chord line, both over the panel chord. The frame is
the one [`freeze_traction_pattern!`](@ref) lofts the surface in, so a control point and
the traction pattern cannot disagree about which way is up.
"""
function chord_frame_coordinates(panel, pos_b)
    le_mid = 0.5 .* (Vector(panel.LE_point_1) .+ Vector(panel.LE_point_2))
    offset = pos_b .- le_mid
    chord = max(panel.chord, eps())
    return (dot(offset, Vector(panel.x_airf)) / chord,
            dot(offset, Vector(panel.z_airf)) / chord)
end

"""
    build_live_polars!(mode::AeroPressure, wing, points; settings)

Build the wing's [`LivePolarState`](@ref) from the reference geometry: fit each panel's
undeformed Kulfan parameters off its surface contour, and record where every control
point sits in that panel's chord frame. Run after
[`build_station_point_map!`](@ref), whose map picks the control points, and while the
mesh still stands on the reference (CAD) structure — the offsets stored here are the
zero the live deformation is measured against.
"""
function build_live_polars!(mode::AeroPressure, wing, points;
                            settings=AirfoilAero.LivePolarSettings())
    panels = wing.vsm_aero.panels
    rot_cad_to_body = wing.R_b_to_c'
    origin_cad = wing.pos_cad
    source = AirfoilAero.LivePolars(AirfoilAero.panel_kulfan_parameters(panels);
                                    settings)
    n_panels = length(panels)
    control_point = Vector{Vector{Int64}}(undef, n_panels)
    control_fraction = Vector{Vector{SimFloat}}(undef, n_panels)
    control_offset = Vector{Vector{SimFloat}}(undef, n_panels)
    for i in 1:n_panels
        assigned = unique(mode.station_point[i])
        length(assigned) >= 2 || error(
            "AeroPressure wing $(wing.name): panel $i reads only " *
            "$(length(assigned)) structural point(s), too few to measure a " *
            "chordwise deformation. Live polars need chordwise stations.")
        frames = [chord_frame_coordinates(panels[i],
                      rot_cad_to_body * (points[idx].pos_cad - origin_cad))
                  for idx in assigned]
        control_point[i] = assigned
        control_fraction[i] = SimFloat[f[1] for f in frames]
        control_offset[i] = SimFloat[f[2] for f in frames]
    end
    maximum(length, control_fraction) <= 2 && @warn(
        "Live polars see no chordwise deformation: every panel maps to at most two " *
        "structural points, which are the chord ends the frame is built on. Add " *
        "chordwise stations to the wing, or the polars only track α and Reynolds.",
        wing=wing.name)
    deflection = [zeros(SimFloat, length(source.basis.x)) for _ in 1:n_panels]
    mode.live = LivePolarState(source, control_point, control_fraction,
                               control_offset, deflection)
    return nothing
end

"""
    update_live_deflection!(mode::AeroPressure, wing, points)

Refresh every panel's camber increment from the deformed structure. Each control point
contributes its current offset off the deformed chord line minus its reference offset,
placed at its current chord fraction; the chord's own rotation and stretch are already
absorbed by the panel frame, which the mesh rebuild took from the deformed leading and
trailing edges. Fills `mode.live.deflection`, ready for
[`refit_live_polars!`](@ref).
"""
function update_live_deflection!(mode::AeroPressure, wing, points)
    state = mode.live::LivePolarState
    panels = wing.vsm_aero.panels
    rot_body_to_world = wing.R_b_to_w::Matrix{SimFloat}
    origin = wing.pos_w::KVec3
    for i in eachindex(panels)
        assigned = state.control_point[i]
        fractions = Vector{SimFloat}(undef, length(assigned))
        deflections = Vector{SimFloat}(undef, length(assigned))
        for (k, idx) in enumerate(assigned)
            pos_b = rot_body_to_world' * (points[idx].pos_w - origin)
            fraction, offset = chord_frame_coordinates(panels[i], pos_b)
            fractions[k] = fraction
            deflections[k] = offset - state.control_offset[i][k]
        end
        state.deflection[i] .= AirfoilAero.control_point_deflection(
            state.source.basis, fractions, deflections)
    end
    return nothing
end

"""
    refit_live_polars!(mode::AeroPressure, wing, alpha) -> Float64

Re-generate every panel's polar at the stored deformation and write it in as a
`TAYLOR` expansion about `alpha` [rad] per panel. Reynolds is per panel from the
solver's air properties and the panel's own apparent wind and chord. Returns the
lowest NeuralFoil analysis confidence over the wing; warns once below `0.5`, which
means a deformed section has left the region the network was trained on.
"""
function refit_live_polars!(mode::AeroPressure, wing, alpha)
    state = mode.live::LivePolarState
    solver = wing.vsm_solver
    panels = wing.vsm_aero.panels
    reynolds = [solver.density * norm(panel.va) * panel.chord / solver.mu
                for panel in panels]
    confidence = AirfoilAero.refresh_live_polars!(state.source, panels, alpha,
                                                 reynolds; deflection=state.deflection)
    confidence < 0.5 && @warn("Live polar off the trained shape range",
        wing=wing.idx, confidence, maxlog=1)
    return confidence
end

"""
    live_polar_alpha(wing) -> Vector{SimFloat}

The angles of attack [rad] to fit each panel's polar about: the previous solve's
converged `lr.alpha_dist`, or, before there is one, each panel's geometric angle from
its own apparent wind. Fitting the first solve of a run about zero would hand it
polars extrapolated the whole way to the real angle.
"""
function live_polar_alpha(wing)
    alpha = collect(SimFloat, wing.vsm_solver.lr.alpha_dist)
    all(iszero, alpha) || return alpha
    return [SimFloat(VortexStepMethod.calculate_relative_alpha_and_relative_velocity(
                panel, zeros(SimFloat, 3))[1])
            for panel in wing.vsm_aero.panels]
end

"""
    solve_with_live_polars!(mode::AeroPressure, wing, points; cold_start=false,
                            max_refits=3)

Solve the wing with polars regenerated from its current shape. The fit is built about
[`live_polar_alpha`](@ref), and the solve repeated while any panel converged outside
the fit window — a local expansion is only continued linearly there, so the answer is
an extrapolation until the fit is moved to where the solve went. Errors if the window
is still not held after `max_refits` fits.
"""
function solve_with_live_polars!(mode::AeroPressure, wing, points;
                                 cold_start=false, max_refits=3)
    state = mode.live::LivePolarState
    update_live_deflection!(mode, wing, points)
    alpha = live_polar_alpha(wing)
    for _ in 1:max_refits
        refit_live_polars!(mode, wing, alpha)
        solve_and_freeze_circulation!(mode, wing; cold_start)
        alpha = collect(SimFloat, wing.vsm_solver.lr.alpha_dist)
        AirfoilAero.polar_drift(state.source, alpha) <= 1.0 && return nothing
    end
    throw(AssertionError("AeroPressure wing $(wing.idx): live polars did not " *
        "settle inside the fit window in $max_refits fits; the solve is moving " *
        "more than $(round(rad2deg(state.source.settings.half_window); digits=1))° " *
        "of angle of attack per refresh."))
end
