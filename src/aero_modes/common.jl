# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Shared aero-mode code: dispatch interface, connector scaffolding, VSM refresh.

# ==================== interface: capability traits ==================== #

"""
    vsm_engine(mode::AbstractAeroModel) -> Union{Nothing, VSMEngine}

The mode's [`VSMEngine`](@ref) (VSM geometry + linearization state), or
`nothing` for modes without one. After construction every `AbstractVSMAero`
carries an engine ([`require_vsm_engine`](@ref) enforces this in
[`setup_aero!`](@ref)); `nothing` only occurs for non-VSM modes and bare
pre-construction markers. Used for the wing's VSM property forwarding and
the VSM-settings loading; per-mode behaviour goes through the dispatch
hooks instead.
"""
vsm_engine(::AbstractAeroModel) = nothing
vsm_engine(mode::AbstractVSMAero) = getfield(mode, :engine)

"""
    has_vsm_engine(mode::AbstractAeroModel) -> Bool

`true` if [`vsm_engine`](@ref)`(mode)` returns an engine. Derived, not a
dispatch point — implement `vsm_engine` instead.
"""
has_vsm_engine(mode::AbstractAeroModel) = vsm_engine(mode) !== nothing

"""
    couples_to_sections(mode::AbstractAeroModel) -> Bool

`true` if the mode needs per-section twist surfaces (auto-creation and
aero-section matching). VSM modes ([`AbstractVSMAero`](@ref)) do.
"""
couples_to_sections(::AbstractAeroModel) = false
couples_to_sections(::AbstractVSMAero) = true

"""
    has_vsm_wing(sys_struct::SystemStructure) -> Bool

`true` if any wing carries a VSM aero mode ([`AbstractVSMAero`](@ref)) whose
state must be refreshed each `vsm_interval`. When `false`, [`refresh_aero!`](@ref)
is a no-op for every wing, so callers skip it and the trailing parameter sync.
"""
has_vsm_wing(sys_struct::SystemStructure) =
    any(wing -> wing.aero isa AbstractVSMAero, sys_struct.wings)

"""
    provides_aero_override(mode::AbstractAeroModel) -> Bool

`true` if the mode supplies frozen body-frame force/moment overrides read by the
compiled RHS (the stored-force path of [`AeroDirect`](@ref)).
"""
provides_aero_override(::AbstractAeroModel) = false

# ==================== interface: required cache tag ==================== #

"""
    aero_mode_tag(mode::AbstractAeroModel) -> String

Short identifier for the mode in the compiled-model cache filename. Required: the
`AbstractAeroModel` fallback errors, so every aero mode must declare its own tag
(no silent default that could collide two distinct modes on one cache file).
Built-ins: `"lin"`, `"dir"`, `"cont"`, `"press"`, `"none"`, `"plate"`.
"""
aero_mode_tag(mode::AbstractAeroModel) = error(
    "aero_mode_tag is not defined for aero mode $(typeof(mode)); " *
    "every aero mode must provide its own cache tag.")

# ==================== interface: diagnostics ==================== #

"""
    calc_aoa(mode::AbstractAeroModel, wing) -> SimFloat

Angle of attack [rad] for `wing` under aero `mode`. Defaults to `NaN`
(undefined); VSM modes read the mid-span geometric AoA (wrapped to [-π, π]) and
[`AeroPlate`](@ref) derives it from the body-frame apparent wind.
"""
calc_aoa(::AbstractAeroModel, wing) = SimFloat(NaN)
function calc_aoa(mode::AbstractVSMAero, wing)
    dist = mode.vsm_solver.sol.alpha_geometric_dist
    n = length(dist)
    return mod(dist[n ÷ 2 + n % 2] + π, 2π) - π
end

"""
    calc_side_slip(wing) -> SimFloat

Side-slip angle [rad] from the body-frame apparent wind. Pure geometry —
the same formula for every aero mode, so it does not dispatch.
"""
calc_side_slip(wing) =
    atan(wing.va_b[2], hypot(wing.va_b[1], wing.va_b[3]))

"""
    normalized_inertia(mode::AbstractAeroModel, wing, points)
        -> (com_cad, inertia)

Normalized (per-unit-mass) inertia of the wing body about its COM in the CAD
frame, with `inertia` in [m²] — multiply by the wing's mass for the physical
tensor [kg·m²]. `inertia` is `nothing` when there is no mass to normalize
by. The default normalizes the wing nodes' point-mass inertia
([`normalized_point_inertia`](@ref)); VSM modes with an `ObjWing` mesh
return the per-unit-mass mesh tensor as-is (its COM is `-T_cad_body`) and
fall back to the point masses otherwise.
"""
normalized_inertia(::AbstractAeroModel, wing, points) =
    normalized_point_inertia(wing, points)

function normalized_inertia(mode::AbstractVSMAero, wing, points)
    tensor = mode.vsm_wing.inertia_tensor
    (isempty(tensor) || all(iszero, tensor)) &&
        return normalized_point_inertia(wing, points)
    return -mode.vsm_wing.T_cad_body, tensor
end

"""
    normalized_point_inertia(wing, points) -> (com_cad, inertia)

Per-unit-mass inertia of the wing's wing nodes treated as point masses
(`extra_mass`), normalized by their total mass. Exact under the construction
invariant `wing.mass == sum of wing-node masses` (the constructor
distributes `set.mass` onto the points). With zero total mass, `com_cad` is
the unweighted centroid and `inertia` is `nothing`.
"""
function normalized_point_inertia(wing, points)
    wing_points = [point for point in points
                   if wing_frame_member(point, wing.idx)]
    masses = [point.extra_mass for point in wing_points]
    total_mass = sum(masses)
    com_cad = total_mass > 0 ?
        sum(masses[j] .* wing_points[j].pos_cad
            for j in eachindex(wing_points)) / total_mass :
        mean([point.pos_cad for point in wing_points])
    total_mass > 0 || return com_cad, nothing
    inertia = zeros(3, 3)
    for (mass, point) in zip(masses, wing_points)
        r = point.pos_cad - com_cad
        inertia += mass * (dot(r, r) * I(3) - r * r')
    end
    return com_cad, inertia / total_mass
end

# ==================== connector scaffolding ==================== #

function rigid_aero_connectors(num_twist_surfaces::Int)
    @variables begin
        va(t)[1:3]
        rho(t)
        R_b_w(t)[1:3, 1:3]
        omega(t)[1:3]
        force(t)[1:3]
        moment(t)[1:3]
    end
    if num_twist_surfaces > 0
        @variables twist(t)[1:num_twist_surfaces] twist_vel(t)[1:num_twist_surfaces] twist_moment(t)[1:num_twist_surfaces]
    else
        twist = nothing
        twist_vel = nothing
        twist_moment = nothing
    end
    return (; va, rho, R_b_w, omega, force, moment,
            twist, twist_vel, twist_moment)
end

function rigid_unknowns(connectors)
    vars = Any[connectors.va, connectors.rho, connectors.R_b_w,
               connectors.omega, connectors.force, connectors.moment]
    connectors.twist === nothing ||
        append!(vars, Any[connectors.twist, connectors.twist_vel,
                          connectors.twist_moment])
    return vars
end

function particle_aero_connectors(num_points::Int; n_delta::Int=0)
    @variables begin
        point_pos(t)[1:3, 1:num_points]
        point_vel(t)[1:3, 1:num_points]
        va(t)[1:3, 1:num_points]
        rho(t)[1:num_points]
        point_force(t)[1:3, 1:num_points]
    end
    if n_delta > 0
        @variables delta(t)[1:n_delta]
    else
        delta = nothing
    end
    return (; point_pos, point_vel, va, rho, point_force, delta)
end

function particle_unknowns(connectors)
    vars = Any[connectors.point_pos, connectors.point_vel, connectors.va,
               connectors.rho, connectors.point_force]
    hasproperty(connectors, :delta) && connectors.delta !== nothing &&
        push!(vars, connectors.delta)
    return vars
end

"""
    frozen_point_force_component(wing::AbstractWing, sys_struct; name, params) -> System

Particle aero component binding each wing node's connector force to its frozen
`point.aero_force_b` flat parameter (synced every refresh). Used by
[`AeroDirect`](@ref), whose refresh precomputes that force. Touching the parameter
here is what registers it.
"""
function frozen_point_force_component(wing, sys_struct; name, params=nothing)
    points = wing_points(sys_struct, wing)
    connectors = particle_aero_connectors(length(points))
    flat_ps = Any[]
    eqs = Equation[]
    for (point_num, point) in enumerate(points)
        force_p = params.points[point.idx].aero_force_b
        push!(flat_ps, force_p)
        eqs = [eqs
               collect(connectors.point_force[:, point_num]) .~ collect(force_p)]
    end
    return System(eqs, t, particle_unknowns(connectors), flat_ps; name)
end

function wing_points(sys_struct, wing)
    return [point for point in sys_struct.points
            if point.is_wing_node && point.wing_idx == wing.idx]
end

"""
    panel_span_signs(wing, spanwise)

Per-panel sign (`±1`) that orients the local span/normal so each panel's `y_airf`
points along `+spanwise` and `z_airf` to the upper surface, independent of section
ordering. Baked at build time because the section order is fixed for a wing.
"""
function panel_span_signs(wing, spanwise)
    refined = wing.vsm_wing.refined_sections
    n = Int(wing.vsm_wing.n_panels)
    map(1:n) do i
        qc1 = 0.75 .* refined[i].LE_point .+ 0.25 .* refined[i].TE_point
        qc2 = 0.75 .* refined[i + 1].LE_point .+ 0.25 .* refined[i + 1].TE_point
        dot(qc1 .- qc2, spanwise) < 0 ? -1.0 : 1.0
    end
end

"""
    build_panel_force_eqs(sec_le, sec_te, sec_va, sec_rho, vind_p, chord_w,
                          cl, cd, cm, spanwise, scale,
                          orient) -> (eqs, vars, panel_force, panel_couple)

Shared per-refined-panel VSM force assembly for the live particle aero modes
([`ContinuousAero`](@ref), [`AeroPressure`](@ref)). Re-expresses
`VortexStepMethod.calc_forces!` symbolically from the frozen circulation: per panel
`i` (between section boundaries `i` and `i+1`) it builds the airfoil axes, chord,
width, effective angle of attack (live apparent wind `sec_va` + frozen induced
velocity `vind_p`), polar coefficients (`cl/cd/cm` callable params) and the lift/drag
directions, and emits the panel force and pitching-moment couple. Returns the panel
equations, the intermediate variables to register, and the
`panel_force`/`panel_couple` arrays for the caller's scatter.

`sec_le`/`sec_te`/`sec_va` are length-`n_panels+1` vectors of body-frame
3-vectors (positions and apparent wind at the section boundaries), `sec_rho`
the matching air densities; they may be live (interpolated from structure) or
constant (a fixed mesh). `spanwise` is the wing spanwise direction, `scale` the
chord-scale factor. `orient` is the per-panel `±1` span/normal sign from
[`panel_span_signs`](@ref) and `chord_w` the per-panel chord blend weight from
[`store_chord_weights!`](@ref).

`delta` is an optional length-`n_panels` vector of symbolic per-panel flap
deflections δ. When given the polars are evaluated as `cl(i, alpha[i], delta[i])`
(the `(α, δ)` tables); when `nothing` the 2-arg `cl(i, alpha[i])` is used, so a
mode without a flap ([`ContinuousAero`](@ref)) is untouched.
"""
function build_panel_force_eqs(sec_le, sec_te, sec_va, sec_rho,
                               vind_p, chord_w, cl, cd, cm, spanwise, scale,
                               orient; delta=nothing)
    n_panels = length(sec_le) - 1
    slots = panel_force_slots(n_panels)
    eqs = Equation[]
    for i in 1:n_panels
        append!(eqs, panel_force_eqs(slots, i,
            (sec_le[i], sec_te[i], sec_le[i + 1], sec_te[i + 1]),
            (sec_va[i], sec_va[i + 1], sec_rho[i], sec_rho[i + 1],
             [vind_p[c, i] for c in 1:3]),
            (PanelPolar(cl, i), PanelPolar(cd, i), PanelPolar(cm, i)),
            spanwise, scale, orient[i], chord_w[i],
            delta === nothing ? nothing : delta[i]))
    end

    return eqs, panel_force_vars(slots), slots.panel_force, slots.panel_couple
end

"""
    panel_force_slots(n_panels)

The symbolic intermediates [`panel_force_eqs`](@ref) writes into, as `n_panels`
columns. A whole wing's system builds every column and an [`AeroPanel`](@ref)
component builds one, so the two emit identical expressions per panel.
"""
function panel_force_slots(n_panels::Int)
    @variables begin
        x_airf(t)[1:3, 1:n_panels]
        y_airf(t)[1:3, 1:n_panels]
        z_airf(t)[1:3, 1:n_panels]
        v_eff(t)[1:3, 1:n_panels]
        chord(t)[1:n_panels]
        width(t)[1:n_panels]
        alpha(t)[1:n_panels]
        q_dyn(t)[1:n_panels]
        dir_lift(t)[1:3, 1:n_panels]
        dir_drag(t)[1:3, 1:n_panels]
        panel_force(t)[1:3, 1:n_panels]
        panel_couple(t)[1:3, 1:n_panels]
    end
    return (; x_airf, y_airf, z_airf, v_eff, chord, width, alpha, q_dyn,
            dir_lift, dir_drag, panel_force, panel_couple)
end

"""The slots of [`panel_force_slots`](@ref) as the variable list a `System` declares."""
panel_force_vars(slots) = Any[values(slots)...]

"""
    panel_force_eqs(slots, i, sections, flow, polars, spanwise, scale, orient,
                    chord_weight, delta)

One panel's aerodynamic equations, writing into column `i` of the symbolic arrays
in `slots`. `sections` is `(le_1, te_1, le_2, te_2)` in body frame, `flow` is
`(va_1, va_2, rho_1, rho_2, v_ind)` and `polars` the `(cl, cd, cm)` callables,
indexed by the panel number the polar tables were built for. `chord_weight` is the
section-1 share of the chord-direction blend ([`store_chord_weights!`](@ref)),
entering as an offset from the midpoint: the equivalent
`chord_weight * te_1 + (1 - chord_weight) * te_2` form leaves no constant term to
fold, which makes a continuous-mode wing 3.7x slower to build symbolically.
`delta` is the flap deflection or `nothing` for the 2-argument polars.

Every quantity it reads belongs to this panel alone, so the same equations serve a
whole-wing system (looped by [`build_panel_force_eqs`](@ref)) and a per-panel
component compiled once and instantiated for each panel.
"""
function panel_force_eqs(slots, i, sections, flow, polars, spanwise, scale,
                         orient, chord_weight, delta)
    (; x_airf, y_airf, z_airf, v_eff, chord, width, alpha, q_dyn,
       dir_lift, dir_drag, panel_force, panel_couple) = slots
    le_1, te_1, le_2, te_2 = sections
    va_1, va_2, rho_1, rho_2, vind = flow
    cl, cd, cm = polars

    lean = chord_weight - 0.5
    chord_vec = (0.5 * (te_1 + te_2) - 0.5 * (le_1 + le_2)) +
        lean * ((te_1 - te_2) - (le_1 - le_2))
    x_unit = chord_vec ./ smooth_norm(chord_vec)
    span_vec = (0.75 * le_1 + 0.25 * te_1) - (0.75 * le_2 + 0.25 * te_2)
    y_unit = orient .* (span_vec ./ smooth_norm(span_vec))
    z_cross = x_unit × (le_1 - le_2)
    z_unit = orient .* (z_cross ./ smooth_norm(z_cross))

    va_panel = 0.5 * (va_1 + va_2)
    v_eff_panel = va_panel + vind
    rho_panel = 0.5 * (rho_1 + rho_2)
    v_eff_crossy = v_eff_panel × y_unit

    lift = evaluate_polar(cl, alpha[i], delta) * q_dyn[i] * chord[i]
    drag = evaluate_polar(cd, alpha[i], delta) * q_dyn[i] * chord[i]
    panel_moment = evaluate_polar(cm, alpha[i], delta) * q_dyn[i] * chord[i]^2

    dir_iva = cos(alpha[i]) .* x_unit .+ sin(alpha[i]) .* z_unit
    lift_cross = dir_iva × y_unit
    drag_cross = spanwise × (lift_cross ./ smooth_norm(lift_cross))

    return [
        chord[i] ~ 0.5 * (smooth_norm(te_1 - le_1) +
                          smooth_norm(te_2 - le_2));
        width[i] ~ smooth_norm(span_vec);
        x_airf[:, i] ~ x_unit;
        y_airf[:, i] ~ y_unit;
        z_airf[:, i] ~ z_unit;
        v_eff[:, i] ~ v_eff_panel;
        alpha[i] ~ atan(v_eff_panel ⋅ z_unit, v_eff_panel ⋅ x_unit);
        q_dyn[i] ~ 0.5 * rho_panel * (v_eff_crossy ⋅ v_eff_crossy);
        dir_lift[:, i] ~ lift_cross ./ smooth_norm(lift_cross);
        dir_drag[:, i] ~ drag_cross ./ smooth_norm(drag_cross);
        panel_force[:, i] ~ (scale * width[i]) .*
            (lift .* collect(dir_lift[:, i]) .+
             drag .* collect(dir_drag[:, i]));
        panel_couple[:, i] ~
            (scale * width[i] * panel_moment / chord[i]) .* z_unit]
end

"""
    evaluate_polar(polar, alpha, delta)

The polar's coefficient at `alpha`, with the flap deflection when the tables carry
one. Keeps [`panel_force_eqs`](@ref) blind to whether its caller bound the panel
index into the polar ([`PanelPolar`](@ref)) or handed it a per-panel callable.
"""
evaluate_polar(polar, alpha, delta) =
    delta === nothing ? polar(alpha) : polar(alpha, delta)

"""
    PanelPolar(polar, panel)

`polar` with its panel index bound, so it is called as `p(alpha)` or
`p(alpha, delta)`. A whole-wing system holds one polar addressed by panel number;
a per-panel component holds a callable that is already only about its own panel.
Applied while the equations are built, so both produce the same expression.
"""
struct PanelPolar{P}
    polar::P
    panel::Int
end

(p::PanelPolar)(alpha) = p.polar(p.panel, alpha)
(p::PanelPolar)(alpha, delta) = p.polar(p.panel, alpha, delta)

# ============ per-panel / per-point parameter addressing ============ #

"""
    PanelAeroList(mode)

`mode`'s panels addressed one at a time, reached as `wings[w].aero.panels[i]`. It
holds the mode rather than its matrices, so a VSM refresh that reallocates them is
still seen, and nothing is copied.
"""
struct PanelAeroList{M}
    mode::M
end

"""
    PanelAero(mode, idx)

One refined panel's frozen aerodynamic data. [`AeroPanel`](@ref) reads its parameters
through this, so the [`remap_path`](@ref) index swap lands on the panel rather than on
the wing, and every read goes to the parent mode's live matrices.
"""
struct PanelAero{M}
    mode::M
    idx::Int
end

"""
    PointAeroList(mode)

`mode`'s wing points addressed one at a time, reached as `wings[w].aero.points[i]`,
where `i` is the global point index. The panel counterpart of [`PanelAeroList`](@ref).
"""
struct PointAeroList{M}
    mode::M
end

"""
    PointAero(mode, idx)

One wing point's frozen aerodynamic data, addressed by its global point index so
[`remap_path`](@ref) can swap it for each [`AeroPointForce`](@ref) instance.
"""
struct PointAero{M}
    mode::M
    idx::Int
end

Base.getindex(list::PanelAeroList, idx::Integer) =
    PanelAero(getfield(list, :mode), Int(idx))
Base.length(list::PanelAeroList) = size(getfield(list, :mode).v_ind, 2)
Base.getindex(list::PointAeroList, idx::Integer) =
    PointAero(getfield(list, :mode), Int(idx))

param_descend(::PanelAeroList) = true
param_descend(::PointAeroList) = true

# `panels` and `points` are addresses into the frozen data, not fields of the mode.
path_step(mode::AbstractVSMAero, key::Symbol) =
    key === :panels ? PanelAeroList(mode) :
    key === :points ? PointAeroList(mode) : getproperty(mode, key)

function Base.getproperty(panel::PanelAero, sym::Symbol)
    mode = getfield(panel, :mode)
    i = getfield(panel, :idx)
    sym === :v_ind && return mode.v_ind[:, i]
    sym === :chord_weight && return mode.chord_weight[i]
    sym === :le_offset_a && return mode.section_le_offset[:, i]
    sym === :te_offset_a && return mode.section_te_offset[:, i]
    sym === :le_offset_b && return mode.section_le_offset[:, i + 1]
    sym === :te_offset_b && return mode.section_te_offset[:, i + 1]
    sym === :cl && return PanelPolar(mode.cl, i)
    sym === :cd && return PanelPolar(mode.cd, i)
    sym === :cm && return PanelPolar(mode.cm, i)
    return error("panel has no aerodynamic field $sym")
end

"""The force a point receives from a mode that stores none — read as a zero of the
right shape while the frozen pattern is still empty."""
const NO_AERO_FORCE = zeros(SimFloat, 3)

function Base.getproperty(point::PointAero, sym::Symbol)
    mode = getfield(point, :mode)
    sym === :offset &&
        return get(mode.point_offset, getfield(point, :idx), NO_AERO_FORCE)
    return error("wing point has no aerodynamic field $sym")
end

# ==================== per-panel decomposition ==================== #

"""
    supports_panel_decomposition(mode) -> Bool

Whether the mode's `PARTICLE_DYNAMICS` aerodynamics can be built as one component per
panel and per point instead of one per wing. A mode opts in by defining
[`aero_inflow_groups`](@ref) and [`aero_scatter_entries`](@ref); everything else about
it — the panel equations, the strut geometry, the wrench — is already shared.
"""
supports_panel_decomposition(::AbstractAeroModel) = false

"""
    aero_inflow_groups(mode, wing, points) -> (groups, section_group)

How each refined section's apparent wind and density are averaged over the wing's
structural points. `groups[g]` is a list of `(point column, weight)` whose weights sum
to one, and `section_group[s]` names the group refined section `s` reads. Sections
sharing a group share the one [`AeroInflow`](@ref) that averages it.

The default is one group per refined section, gathered from the same bounding struts
[`aero_geometry_entries`](@ref) reconstructs that section's corners from, so a
section's inflow and its geometry read the same points. Overriding this with a
wing-wide mean drops the rotational part of the velocity field — a rigid rotation
about the point centroid averages to zero — leaving the panels without rate damping.
"""
function aero_inflow_groups(mode, wing, points)
    column = aero_section_columns(wing, points)
    sections = eachindex(section_interp_caches(mode)[1])
    groups = [strut_inflow_weights(mode, section, column) for section in sections]
    return groups, collect(sections)
end

"""
    aero_scatter_entries(mode, wing, points) -> Vector

The linear map from panel loads to point forces as `(panel, point column, force
weight, couple weight)`: point `k`'s body-frame force gains `force weight ·
panel_force + couple weight · panel_couple`. The weights are constants of the mesh, so
the wiring layer carries the whole scatter and no component holds it.
"""
function aero_scatter_entries end

"""
    aero_point_offset(mode, params, wing_idx, point_idx)

The constant body-frame force a wing point receives on top of the scattered panel
loads, or `nothing`. [`AeroPressure`](@ref) puts its frozen surface traction here, net
of the share of each panel's frozen total that the scatter already re-adds.
"""
aero_point_offset(::AbstractAeroModel, params, wing_idx, point_idx) = nothing

"""
    aero_geometry_entries(mode, wing, points) -> Vector

Where each panel's four section corners come from, as `(panel, corner, point column,
weight)` with `corner` one of `:le_a`, `:te_a`, `:le_b`, `:te_b`. This is
[`interp_sections`](@ref) split in two: the strut weights the wiring gathers, and the
frozen billow offset the panel adds as a parameter. Mode-independent — every
continuous mode reconstructs its sections the same way.
"""
function aero_geometry_entries(mode, wing, points)
    column = aero_section_columns(wing, points)
    left, weight, _, _ = section_interp_caches(mode)
    entries = Tuple{Int, Symbol, Int, SimFloat}[]
    for panel in 1:(length(left) - 1), (side, section) in ((:a, panel), (:b, panel + 1))
        for (strut, share) in ((left[section], weight[section]),
                               (left[section] + 1, 1.0 - weight[section]))
            share == 0.0 && continue
            push!(entries, (panel, Symbol(:le_, side), column[(strut, :LE)], share))
            push!(entries, (panel, Symbol(:te_, side), column[(strut, :TE)], share))
        end
    end
    return entries
end

"""
    const COLLOCATION_CHORD_FRAC = 0.75

Chordwise station, as a fraction from LE to TE, that a section's inflow is sampled
at. Thin-airfoil theory enforces flow tangency at the 3/4-chord point, so sampling
there lets a chord twist rate reach `alpha` as downwash and be damped. Sampling at
mid-chord instead leaves the twist rate invisible to the aero, or, for a chord
pivoting between mid- and 3/4-chord, damped with the wrong sign.
"""
const COLLOCATION_CHORD_FRAC = 0.75

"""
    strut_inflow_weights(mode, section, column) -> Vector{Tuple{Int, SimFloat}}

The `(point column, weight)` pairs averaging refined section `section`'s inflow over
its two bounding struts, each strut blending its LE and TE station at
[`COLLOCATION_CHORD_FRAC`](@ref). With no chord twist rate the two stations share a
velocity, so the blend is inert in steady flight whatever the fraction.
"""
function strut_inflow_weights(mode, section::Int, column)
    left, weight, _, _ = section_interp_caches(mode)
    pairs = Tuple{Int, SimFloat}[]
    for (strut, share) in ((left[section], weight[section]),
                           (left[section] + 1, 1.0 - weight[section]))
        share == 0.0 && continue
        push!(pairs, (column[(strut, :LE)],
                      (1 - COLLOCATION_CHORD_FRAC) * share))
        push!(pairs, (column[(strut, :TE)], COLLOCATION_CHORD_FRAC * share))
    end
    return pairs
end

"""
    scatter_totals!(totals, panel, point, force_weight, couple_weight)

Accumulate one contribution into the `(panel, point) → [force weight, couple weight]`
map [`aero_scatter_entries`](@ref) builds, so a panel that reaches the same point
twice becomes one wiring edge instead of two.
"""
function scatter_totals!(totals, panel::Int, point::Int, force_weight, couple_weight)
    total = get!(zero_weight_pair, totals, (panel, point))
    total[1] += force_weight
    total[2] += couple_weight
    return nothing
end

"""A fresh `[force weight, couple weight]` accumulator."""
zero_weight_pair() = zeros(SimFloat, 2)

"""The `(panel, point, force weight, couple weight)` list of a
[`scatter_totals!`](@ref) map, ordered so the wiring is reproducible."""
scatter_entry_list(totals) =
    sort!([(panel, point, total[1], total[2])
           for ((panel, point), total) in totals])

"""
    store_chord_weights!(chord_weight, body_aero)

Freeze each refined panel's chord blend weight into `chord_weight` (n_panels): the
share of section 1 in the leading- and trailing-edge blend whose difference gives
the panel's chord direction. `VortexStepMethod` weights that blend by panel
spacing rather than taking the midpoint, so a panel next to a narrower neighbour
leans towards it, and the weight follows the mesh as the wing deforms. Written at
the same refresh as [`store_induced_velocity!`](@ref).
"""
function store_chord_weights!(chord_weight, body_aero)
    panels = body_aero.panels
    n_panels = length(panels)
    length(chord_weight) == n_panels || error(
        "chord-weight buffer is stale ($(length(chord_weight)) for $n_panels " *
        "panels); reinitialize the model.")
    if n_panels < 2
        fill!(chord_weight, 0.5)
        return nothing
    end
    for i in 1:n_panels
        own = panels[i].width
        spacing_share = if i == 1
            own / (own + panels[2].width)
        elseif i == n_panels
            panels[n_panels - 1].width / (panels[n_panels - 1].width + own)
        else
            0.25 * (panels[i - 1].width / (panels[i - 1].width + own) +
                    own / (own + panels[i + 1].width) + 1)
        end
        chord_weight[i] = 1 - spacing_share
    end
    return nothing
end

"""
    size_frozen_panel_buffers!(mode, n_panels)

Size the per-panel buffers every live VSM mode freezes at a refresh — the induced
velocity and the chord blend weight — reallocating only when the refined mesh
changed. The chord weights start at the midpoint, which is what they are on a
uniformly panelled wing.
"""
function size_frozen_panel_buffers!(mode, n_panels::Int)
    size(mode.v_ind) == (3, n_panels) ||
        (mode.v_ind = zeros(SimFloat, 3, n_panels))
    length(mode.chord_weight) == n_panels ||
        (mode.chord_weight = fill(SimFloat(0.5), n_panels))
    return nothing
end

"""
    store_induced_velocity!(v_ind, body_aero, gamma)

Freeze the converged circulation into the buffer `v_ind` (3 × n_panels): each
refined panel's induced velocity is `AIC · gamma`, the same product the VSM
gamma loop converged on. Shared by the live particle modes.
"""
function store_induced_velocity!(v_ind, body_aero, gamma)
    n_panels = length(body_aero.panels)
    size(v_ind) == (3, n_panels) || error(
        "induced-velocity buffer is stale ($(size(v_ind)) for $n_panels " *
        "panels); reinitialize the model.")
    aic = body_aero.AIC
    for component in 1:3
        for i in 1:n_panels
            acc = 0.0
            for j in 1:n_panels
                acc += aic[i, j, component] * gamma[j]
            end
            v_ind[component, i] = acc
        end
    end
    return nothing
end

"""
    aero_component(mode::AbstractAeroModel, wing::AbstractWing, sys_struct; name, params) -> System

Build the aero subsystem for `wing`, selected by dispatch on both the wing's
`aero` model and its dynamics type ([`RigidWing`](@ref)/[`ParticleWing`](@ref)).
Returns a `System` exposing the connectors fixed by the dynamics type, all in
the wing body frame (the wiring layer `aero_eqs!` drives inputs and reads
outputs; connectors a mode ignores still exist for binding):
- `RIGID_DYNAMICS` (`num = length(wing.twist_surface_idxs)`): in `va[1:3]`,
  `rho`, `R_b_w[1:3,1:3]`, `omega[1:3]`, `twist[1:num]`, `twist_vel[1:num]`;
  out `force[1:3]`, `moment[1:3]`, `twist_moment[1:num]`.
- `PARTICLE_DYNAMICS` (`np = number of wing nodes`): in `point_pos[1:3,1:np]`,
  `point_vel[1:3,1:np]`, `va[1:3,1:np]`, `rho[1:np]`; out `point_force[1:3,1:np]`.

A mode supports a dynamics type by defining the matching method (rigid,
particle, or both). Add a method on a custom `AbstractAeroModel` subtype to
plug in your own aerodynamics.
"""
function aero_component end

"""
    validate_aero_component(subsys, wing)

Check the built aero `subsys` exposes the connectors the wiring layer needs for
the wing's `dynamics_type`; error naming the missing connector otherwise.
"""
function validate_aero_component(subsys, wing)
    if wing.dynamics_type == RIGID_DYNAMICS
        required = Symbol[:va, :rho, :R_b_w, :omega, :force, :moment]
        length(wing.twist_surface_idxs) > 0 &&
            append!(required, [:twist, :twist_vel, :twist_moment])
    else
        required = Symbol[:point_pos, :point_vel, :va, :rho, :point_force]
    end
    required_str = join(required, ", ")
    for con in required
        hasproperty(subsys, con) || error(
            "Wing $(wing.name): aero component is missing required " *
            "connector `$con`. Required: $required_str.")
    end
    return nothing
end

# ==================== refresh orchestrator ==================== #

"""
    sync_aero_density!(wing, am)

Set the wing's VSM solver air density to `air_density(am, wing.pos_w[3])`, the same
altitude-dependent density the symbolic RHS uses to dimensionalize aero forces
(see `aero_eqs.jl`). Keeps the VSM solve and the model consistent on dynamic
pressure. No-op for non-VSM aero modes.
"""
function sync_aero_density!(wing, am)
    wing.aero isa AbstractVSMAero || return nothing
    wing.vsm_solver.density = air_density(am, wing.pos_w[3])
    return nothing
end

"""
    point_acceleration_w(point, wing_frame, wing_vel) -> KVec3

A free particle's world-frame acceleration, rebuilt from what the struct already
carries: its net force per unit mass less the world-frame and body-frame damping.
The same expression `point_eqs!` binds to `acc`, which is where the monolith's
fitted wing reads its own `acc_w` from.
"""
function point_acceleration_w(point, wing_frame, wing_vel)
    damping = collect(point.world_frame_damping) .* collect(point.vel_w) .+
        body_frame_damp_accel(point.vel_w, point.body_frame_damping, wing_frame,
                              collect(wing_vel))
    return collect(point.force) ./ point.total_mass .- damping
end

"""
    wing_kinematics_from_points!(wing, points, set, am, wind_mode;
                                 zp1, zp2, yp1, yp2, origin, aero_points)

Recompute a KINEMATIC/PARTICLE wing's kinematic state directly from the current point
positions/velocities, for a backend that does not carry these quantities as state.
`zp1`, `zp2`, `yp1`, `yp2` and `origin` are [`WeightedRefPoints`](@ref), so each
reference is the weighted blend of its points. Writes the body frame `R_b_to_w`
([`wing_frame_columns`](@ref)), the origin pose `pos_w`/`vel_w`, the frame's own
`ω_b` ([`body_frame_omega`](@ref)), the origin's acceleration `acc_w`
([`point_acceleration_w`](@ref)), the reported scalars
([`write_wing_scalars!`](@ref)), the wing apparent wind `va_b`, and each aero point's
`va_b = R'·(wind − vel)`, the wind coming from the height profile or, under
[`PerPointWind`](@ref), from `point.wind_vec` and `wing.wind_vec`, which the caller
owns — the same quantities the monolith's `get_all_state` copies out of the integrator.
"""
function wing_kinematics_from_points!(wing, points, set, am, wind_mode::WindMode;
        zp1, zp2, yp1, yp2, origin, aero_points,
        base_point = 0, twist_surfaces = nothing)
    pos_z1 = get_ref_position_from_points(points, zp1)
    pos_z2 = get_ref_position_from_points(points, zp2)
    pos_y1 = get_ref_position_from_points(points, yp1)
    pos_y2 = get_ref_position_from_points(points, yp2)
    vel_z1 = get_ref_position_from_points(points, zp1; field = :vel_w)
    vel_z2 = get_ref_position_from_points(points, zp2; field = :vel_w)
    vel_y1 = get_ref_position_from_points(points, yp1; field = :vel_w)
    vel_y2 = get_ref_position_from_points(points, yp2; field = :vel_w)
    x, y, z = wing_frame_columns(pos_z1, pos_z2, pos_y1, pos_y2)
    R = hcat(x, y, z)
    wing.R_b_to_w .= R
    wing.pos_w .= get_ref_position_from_points(points, origin)
    wing.vel_w .= get_ref_position_from_points(points, origin; field = :vel_w)
    axes = (collect(x), collect(y), collect(z))
    wing.ω_b .= body_frame_omega(axes,
        wing_frame_rates(pos_z1, pos_z2, pos_y1, pos_y2,
            (vel_z1, vel_z2, vel_y1, vel_y2), axes))
    fill!(wing.acc_w, 0.0)
    for (idx, weight) in zip(origin.ids, origin.weights)
        wing.acc_w .+= weight .* point_acceleration_w(points[idx], R, wing.vel_w)
    end
    profile = wind_mode isa PerPointWind ? nothing : WindFactor(am, set.profile_law)
    isnothing(profile) || (wing.wind_vec .= profile(wing.pos_w[3]) .* set.wind_vec)
    wing.va_b .= R' * (wing.wind_vec .- wing.vel_w .+ wing.wind_disturb)
    for point_idx in aero_points
        point = points[point_idx]
        wind = isnothing(profile) ? point.wind_vec :
            profile(point.pos_w[3]) .* set.wind_vec
        point.va_b .= R' * (wind .- point.vel_w)
    end
    write_wing_scalars!(wing, points; base_point, twist_surfaces)
    return nothing
end

"""
    write_wing_scalars!(wing, points; base_point, alpha_b, twist_surfaces) -> nothing

Fill a fitted wing's reported scalars — heading, course, elevation, azimuth,
distance and angle of attack with their rates — from its freshly rebuilt pose,
via the one definition in [`wing_scalar_kinematics`](@ref).

`turn_rate` follows the `ω_b` the caller fitted from the frame's own ref points, the
`_acc` scalars the `wing.acc_w` it wrote, and `turn_acc` the `alpha_b` it passed.
`alpha_b` defaults to zero, which is what the monolith binds a fitted wing's to.
"""
function write_wing_scalars!(wing, points; base_point = 0, alpha_b = zeros(3),
                             twist_surfaces = nothing)
    rel_pos = base_point == 0 ? collect(wing.pos_w) :
        collect(wing.pos_w) .- collect(points[base_point].pos_w)
    idxs = wing.twist_surface_idxs
    twist_offset = if isnothing(twist_surfaces) || isempty(idxs)
        0.0
    else
        half = idxs[1] + length(idxs) ÷ 2 - 1
        0.5 * twist_surfaces[half].twist + 0.5 * twist_surfaces[half + 1].twist
    end
    R_b_to_w = collect(wing.R_b_to_w)
    e_x = R_b_to_w[:, 1]
    scalars = wing_scalar_kinematics(; rel_pos, e_x,
        R_t_to_w = sym_calc_R_t_to_w(rel_pos),
        R_v_to_w = calc_R_v_to_w(rel_pos, e_x),
        R_b_to_w, vel = collect(wing.vel_w), acc = collect(wing.acc_w),
        omega_b = collect(wing.ω_b), alpha_b = collect(alpha_b),
        va_b = collect(wing.va_b), twist_offset)
    wing.heading = scalars.heading
    wing.turn_rate .= scalars.turn_rate
    wing.turn_acc .= scalars.turn_acc
    wing.elevation = scalars.elevation
    wing.elevation_vel = scalars.elevation_vel
    wing.elevation_acc = scalars.elevation_acc
    wing.azimuth = scalars.azimuth
    wing.azimuth_vel = scalars.azimuth_vel
    wing.azimuth_acc = scalars.azimuth_acc
    wing.course = scalars.course
    wing.aoa = scalars.angle_of_attack
    return nothing
end

"""
    refresh_aero!(sam::SymbolicAWEModel; vsm_min_wind=0.5, cold_start=false,
                  vsm_warn_on_fail=false)

Refresh each wing's aerodynamic state, dispatching on the wing's aero mode
([`refresh_rigid_aero!`](@ref) / [`refresh_particle_aero!`](@ref)). Runs on the
low-frequency VSM-update schedule (`vsm_interval`), not the compiled RHS. Reads
per-point apparent wind from the `points` (populated by `update_sys_struct!`,
which always runs first), so it does not re-extract integrator state.

**RIGID_DYNAMICS VSM modes:** compute wind-axis coefficients (CL, CD, CS, CM, cm)
at the operating point, plus the `ForwardDiff` Jacobian over `[α, β, ω₁, ω₂, ω₃,
θ_twist…]` (`AeroLinearized`) or the frozen forces (`AeroDirect`).

**PARTICLE_DYNAMICS VSM modes:** full nonlinear VSM solve with per-point force
distribution. Non-VSM modes (`AeroNone`/`AeroPlate`) are no-ops, so callers
should gate this on [`has_vsm_wing`](@ref).

`cold_start` discards the previous solve's circulation as the warm start
([`safe_vsm_solve!`](@ref)), making the result a function of the current state
alone. The per-step refresh wants the warm start; the first solve after a
[`reinit!`](@ref) must not have it, or it inherits whatever ran on this model
before.

`vsm_warn_on_fail` downgrades a [`VSMSolveFailure`](@ref) to a warning, leaving
the failing wing at the aero state of its last converged solve.
"""
function refresh_aero!(sam::SymbolicAWEModel; vsm_min_wind=0.5,
                       cold_start=false, vsm_warn_on_fail=false)
    wings = sam.sys_struct.wings
    twist_surfaces = sam.sys_struct.twist_surfaces
    points = sam.sys_struct.points

    length(wings) == 0 && return nothing

    for wing in wings
        sync_aero_density!(wing, sam.am)
    end

    for wing in wings
        wing.dynamics_type == RIGID_DYNAMICS || continue
        try
            refresh_rigid_aero!(wing.aero, wing, sam.am, twist_surfaces;
                                vsm_min_wind)
        catch failure
            warn_or_rethrow(failure, vsm_warn_on_fail)
        end
    end

    any(w.dynamics_type === PARTICLE_DYNAMICS for w in wings) ||
        return nothing
    va_point_b_vals = Matrix{SimFloat}(undef, 3, length(points))
    for point in points
        @views va_point_b_vals[:, point.idx] .= point.va_b
    end
    for wing in wings
        wing.dynamics_type == PARTICLE_DYNAMICS || continue
        apply_flap_delta!(wing.aero, wing, sam.sys_struct)
        try
            refresh_particle_aero!(wing.aero, wing, points, va_point_b_vals;
                                   vsm_min_wind, cold_start)
        catch failure
            warn_or_rethrow(failure, vsm_warn_on_fail)
        end
    end

    nothing
end

"""
    warn_or_rethrow(failure, vsm_warn_on_fail)

Rethrow `failure`, unless it is a [`VSMSolveFailure`](@ref) and `vsm_warn_on_fail`
is set: then warn instead, which leaves that wing's frozen aero state and
circulation at its last converged solve.
"""
function warn_or_rethrow(failure, vsm_warn_on_fail)
    (vsm_warn_on_fail && failure isa VSMSolveFailure) || rethrow(failure)
    @warn "$(failure.msg) Reusing the last converged aero state."
    return nothing
end

"""
    refresh_rigid_aero!(mode, wing, am, twist_surfaces; vsm_min_wind=0.5)

Refresh a `RIGID_DYNAMICS` wing's aero state, dispatched on its aero `mode`:
- `AeroNone` / any non-VSM mode → no-op (fallback).
- `AeroLinearized` → compute the baseline coefficients ([`rigid_aero_baseline!`](@ref))
  and the `ForwardDiff` Jacobian `d(coeffs)/d(inputs)` into `wing.aero_jac`.
- `AeroDirect` → compute the baseline coefficients and apply the frozen body-frame
  force/moment; below `vsm_min_wind` everything is zeroed.
"""
refresh_rigid_aero!(::AbstractAeroModel, wing, am, twist_surfaces;
                    vsm_min_wind=0.5) = nothing

"""
    refresh_particle_aero!(mode, wing, points, va_point_b_vals;
                           vsm_min_wind=0.5, cold_start=false)

Refresh a `PARTICLE_DYNAMICS` wing's aero state, dispatched on its aero `mode`:
- `AeroNone` / any non-VSM mode → no-op (fallback).
- `AeroDirect` → full nonlinear VSM solve with per-section apparent wind, then
  distribute panel forces onto the wing's structural points
  ([`distribute_panel_forces_to_points!`](@ref)); below `vsm_min_wind` the point
  forces are zeroed.
- `AeroLinearized` → unsupported (errors).

`cold_start` forwards to [`safe_vsm_solve!`](@ref).
"""
refresh_particle_aero!(::AbstractAeroModel, wing, points, va_point_b_vals;
                       vsm_min_wind=0.5, cold_start=false) = nothing

# ==================== per-wing lifecycle ==================== #

"""
    remake_aero!(mode, wing, set, vsm_set, points, twist_surfaces)

Rebuild the mode's aero engine from `set`/`vsm_set` (the `remake_vsm` path in
`reinit!`, used after editing settings). Default no-op; VSM modes recreate the
VSM wing/aero/solver, re-transform sections to the body frame, re-match aero
sections to structure, and rebuild the twist-surface / point mappings.
"""
remake_aero!(::AbstractAeroModel, wing, set, vsm_set, points, twist_surfaces) =
    nothing

"""
    attach_engine!(mode, engine) -> mode

Attach a freshly built [`VSMEngine`](@ref) to a VSM aero `mode` during wing
construction. Built-in modes reconstruct (so the concrete engine type lands in
the wing's type parameter, removing the abstract-`engine` dispatch from the RHS);
the default mutates a custom mode in place.
"""
attach_engine!(mode::AbstractVSMAero, engine::VSMEngine) =
    (setfield!(mode, :engine, engine); mode)

function remake_aero!(mode::AbstractVSMAero, wing, set, vsm_set, points,
                      twist_surfaces)
    vsm_set isa VortexStepMethod.VSMSettings || error(
        "remake_aero!: VSM wing $(wing.idx) needs a VSMSettings, " *
        "got $(typeof(vsm_set)).")
    wing.vsm_wing = create_vsm_wing(set, vsm_set;
        prn=false, sort_sections=false)
    wing.vsm_aero = VortexStepMethod.BodyAerodynamics([wing.vsm_wing])
    wing.vsm_solver = VortexStepMethod.Solver(wing.vsm_aero, vsm_set)

    # Transform sections CAD → body frame (matches the SystemStructure constructor)
    transform_vsm_sections_to_body!(wing;
        aero_z_offset=(wing.dynamics_type == PARTICLE_DYNAMICS ? nothing :
                       wing.aero_z_offset))

    match_aero_sections_to_structure!(wing, points; twist_surfaces)

    if wing.dynamics_type == RIGID_DYNAMICS && !isempty(wing.twist_surface_idxs)
        compute_spatial_twist_surface_mapping!(wing, twist_surfaces, points)
    end
    if wing.dynamics_type == PARTICLE_DYNAMICS &&
       !isnothing(wing.wing_segments)
        wing.point_to_vsm_point =
            build_point_to_vsm_point_mapping(wing.wing_segments)
    end
    return nothing
end

"""
    validate_aero_structure(mode, wing, points; prn=false)

Check structural invariants the mode's compiled equations rely on (run at build).
Default no-op; VSM `PARTICLE_DYNAMICS` wings verify each unrefined section has its
LE and TE structural point mapped (interior chord control points are allowed and
need no mapping).
"""
validate_aero_structure(::AbstractAeroModel, wing, points; prn=false) = nothing

function validate_aero_structure(::AbstractVSMAero, wing, points; prn=false)
    wing.dynamics_type == PARTICLE_DYNAMICS || return nothing
    @assert !isnothing(wing.point_to_vsm_point) "PARTICLE_DYNAMICS wing $(wing.idx) missing point_to_vsm_point mapping"

    n_sections = length(wing.vsm_wing.unrefined_sections)
    mapped = values(wing.point_to_vsm_point)
    for section_idx in 1:n_sections, le_or_te in (:LE, :TE)
        @assert (section_idx, le_or_te) in mapped "PARTICLE_DYNAMICS wing $(wing.idx): section $(section_idx) missing its $(le_or_te) point"
    end

    n_wing_points = count(p -> p.is_wing_node && p.wing_idx == wing.idx, points)
    prn && println("✓ PARTICLE_DYNAMICS wing $(wing.idx) validated: $(n_wing_points) points, $(n_sections) sections, $(length(wing.vsm_aero.panels)) panels")
    return nothing
end

"""
    setup_aero!(mode, wing, points, twist_surfaces; prn=false)

Construction-time aero setup for `wing`, dispatched on its aero `mode` (default
no-op). VSM modes transform the VSM panels into the body frame and, for
section-coupled wings, auto-create twist surfaces, match aero sections to
structure, and build the twist-surface / structural↔panel mappings. A custom mode
adds a method to participate in construction without editing the SystemStructure
constructor. Runs after [`setup_wing_frame!`](@ref) (which sets the body frame).
"""
setup_aero!(::AbstractAeroModel, wing, points, twist_surfaces; prn=false) =
    nothing

"""
    require_vsm_engine(mode, wing) -> VSMEngine

Return the mode's [`VSMEngine`](@ref), erroring with construction advice when it
is missing (a bare `AeroDirect()`/`AeroLinearized()` marker attached to a wing
that was not built via [`VSMWing`](@ref)). Called once, in [`setup_aero!`](@ref);
after construction every `AbstractVSMAero` is guaranteed to carry an engine.
"""
function require_vsm_engine(mode, wing)
    engine = vsm_engine(mode)
    engine === nothing && error(
        "Wing $(wing.name): aero mode $(typeof(mode)) has no VSM engine. " *
        "Construct the wing via VSMWing (or attach a VSMEngine to the mode) " *
        "to use VSM aerodynamics.")
    return engine
end

"""
    transform_vsm_sections_to_body!(wing; aero_z_offset=nothing)

Move the wing's VSM sections/panels from the CAD frame into the body frame
(translate to `wing.pos_cad`, rotate by `wing.R_b_to_c'`) and reinit the panels.
With `aero_z_offset` set, also apply the chordwise aero z-offset (RIGID wings);
PARTICLE wings pass `nothing`. Shared by [`setup_aero!`](@ref) and
[`remake_aero!`](@ref).
"""
function transform_vsm_sections_to_body!(wing; aero_z_offset=nothing)
    vsm_wing = wing.vsm_wing
    vsm_wing.T_cad_body .= wing.pos_cad
    adjust_vsm_panels_to_origin!(vsm_wing, wing.pos_cad)
    rotate_vsm_sections!(vsm_wing, wing.R_b_to_c')
    vsm_wing.R_cad_body .= wing.R_b_to_c
    isnothing(aero_z_offset) || apply_aero_z_offset!(vsm_wing, aero_z_offset)
    VortexStepMethod.reinit!(wing.vsm_aero)
    return nothing
end

function setup_aero!(mode::AbstractVSMAero, wing, points, twist_surfaces;
                     prn=false)
    require_vsm_engine(mode, wing)
    if wing.dynamics_type == RIGID_DYNAMICS
        transform_vsm_sections_to_body!(wing; aero_z_offset=wing.aero_z_offset)

        if couples_to_sections(mode) && isempty(wing.twist_surface_idxs)
            error("Section-coupled aero on RIGID wing $(wing.idx) requires " *
                  "explicit twist_surfaces covering its LE/TE structural " *
                  "sections; none were declared. Add them to the wing.")
        end
        couples_to_sections(mode) &&
            match_aero_sections_to_structure!(wing, points; twist_surfaces)
        isempty(wing.twist_surface_idxs) ||
            compute_spatial_twist_surface_mapping!(wing, twist_surfaces, points)
        compute_twist_surface_geometry!(wing, twist_surfaces, points)
        for twist_surface_idx in wing.twist_surface_idxs
            twist_surfaces[twist_surface_idx].le_pos .-= wing.com_offset_b
        end
    else  # PARTICLE_DYNAMICS
        isnothing(wing.origin) || transform_vsm_sections_to_body!(wing)
        couples_to_sections(mode) &&
            match_aero_sections_to_structure!(wing, points; twist_surfaces)
        setup_particle_point_mapping!(wing, points, twist_surfaces)
    end
    return nothing
end

"""
    resize_aero_state!(mode, wing)

Resize the mode's per-wing aero state after `wing.twist_surface_idxs` is
resolved (name resolution can change the twist-surface count the initial
sizing estimated from `n_unrefined`). Default no-op; VSM modes resize
`aero_y`/`aero_x`/`aero_jac` for `RIGID_DYNAMICS` wings.
"""
resize_aero_state!(::AbstractAeroModel, wing) = nothing

function resize_aero_state!(mode::AbstractVSMAero, wing)
    wing.dynamics_type == RIGID_DYNAMICS || return nothing
    n_twist_surfaces = length(wing.twist_surface_idxs)
    num_aero_outputs = 6 + n_twist_surfaces
    num_aero_inputs = 5 + n_twist_surfaces
    if length(mode.aero_x) != num_aero_outputs ||
            length(mode.aero_y) != num_aero_inputs
        mode.aero_y = zeros(SimFloat, num_aero_inputs)
        mode.aero_x = zeros(SimFloat, num_aero_outputs)
        mode.aero_jac = zeros(
            SimFloat, num_aero_outputs, num_aero_inputs)
    end
    return nothing
end

"""
    init_aero_state!(mode, wing, va_b_init)

Initialize the mode's aero state from the initial body-frame apparent wind
`va_b_init` (runs in `update_sys_struct!`, before the first refresh). Default
no-op; VSM modes write the operating-point angles α, β into `aero_y`.
"""
init_aero_state!(::AbstractAeroModel, wing, va_b_init) = nothing

function init_aero_state!(mode::AbstractVSMAero, wing, va_b_init)
    aero_y = mode.aero_y
    length(aero_y) < 2 && return nothing
    aero_y .= 0.0
    aero_y[1] = atan(va_b_init[3], va_b_init[1])
    aero_y[2] = atan(va_b_init[2],
        hypot(va_b_init[1], va_b_init[3]))
    return nothing
end

# ============ live refined-section geometry (shared by continuous modes) ========

"""
    build_section_interp(vsm_wing) -> (left, weight, le_offset, te_offset)

Freeze the refined-section → unrefined-strut interpolation of `vsm_wing`: refined
section `s` sits at `weight[s]·strut[left[s]] + (1−weight[s])·strut[left[s]+1]`.
`le_offset`/`te_offset` (3 × n_sections) are each refined section's body-frame
displacement off that straight strut line (nonzero only for `BILLOWING`). All are
constants of the mesh, baked into the generated equations, so they enter the
model-cache hash via each mode's `aero_hash_id`.
"""
function build_section_interp(vsm_wing)
    n_panels = Int(vsm_wing.n_panels)
    n_sections = n_panels + 1
    n_struts = Int(vsm_wing.n_unrefined_sections)
    n_struts >= 2 || error(
        "continuous aero: need at least 2 unrefined sections, got $n_struts.")
    left_raw = vsm_wing.refined_section_left_idx
    weight_raw = vsm_wing.refined_section_weight
    if length(left_raw) == n_sections && length(weight_raw) == n_sections
        left = Int64.(left_raw)
        weight = SimFloat.(weight_raw)
    elseif n_struts == n_sections
        left = [min(Int64(s), Int64(n_struts - 1)) for s in 1:n_sections]
        weight = [s < n_sections ? 1.0 : 0.0 for s in 1:n_sections]
    else
        error("continuous aero: VSM wing has no refined-section interpolation " *
              "cache ($(length(left_raw)) entries for $n_sections sections).")
    end
    le_offset = zeros(SimFloat, 3, n_sections)
    te_offset = zeros(SimFloat, 3, n_sections)
    refined = vsm_wing.refined_sections
    unrefined = vsm_wing.unrefined_sections
    if length(refined) == n_sections
        for s in 1:n_sections
            strut = left[s]
            w = weight[s]
            line_le = w .* unrefined[strut].LE_point .+
                (1.0 - w) .* unrefined[strut + 1].LE_point
            line_te = w .* unrefined[strut].TE_point .+
                (1.0 - w) .* unrefined[strut + 1].TE_point
            le_offset[:, s] .= refined[s].LE_point .- line_le
            te_offset[:, s] .= refined[s].TE_point .- line_te
        end
    end
    return left, weight, le_offset, te_offset
end

"""
    section_interp_caches(mode) -> (left, weight, le_offset, te_offset)

The frozen strut-interpolation caches a continuous mode stored via
[`build_section_interp`](@ref).
"""
section_interp_caches(mode::AbstractVSMAero) =
    (mode.section_left_strut, mode.section_left_weight,
     mode.section_le_offset, mode.section_te_offset)

"""
    aero_section_columns(wing, points) -> Dict{Tuple{Int64,Symbol},Int}

Map `(unrefined_section, :LE/:TE)` to the connector column (position in `points`)
of the structural station point there, via `wing.point_to_vsm_point`. Interior
chord control points (absent from the map) are skipped.
"""
function aero_section_columns(wing, points)
    point_to_vsm = wing.point_to_vsm_point
    isnothing(point_to_vsm) && error(
        "continuous aero: wing $(wing.name) missing point_to_vsm_point mapping.")
    column = Dict{Tuple{Int64, Symbol}, Int}()
    for (k, point) in enumerate(points)
        entry = get(point_to_vsm, point.idx, nothing)
        entry === nothing && continue
        column[entry] = k
    end
    return column
end

"""
    interp_strut(values, left, weight, s)

Refined section `s` from its two bounding struts:
`weight[s]·values[left[s]] + (1−weight[s])·values[left[s]+1]`.
"""
interp_strut(values, left, weight, s) =
    weight[s] .* values[left[s]] .+ (1.0 - weight[s]) .* values[left[s] + 1]

"""
    interp_sections(strut_le, strut_te, left, weight, le_offset, te_offset)
        -> (sec_le, sec_te)

Interpolate per-strut LE/TE positions to every refined section and add the frozen
billow offset. Works on symbolic (from live connectors) or numeric (from `pos_w`)
per-strut vectors alike.
"""
function interp_sections(strut_le, strut_te, left, weight, le_offset, te_offset)
    n_sections = length(left)
    sec_le = [interp_strut(strut_le, left, weight, s) .+ le_offset[:, s]
              for s in 1:n_sections]
    sec_te = [interp_strut(strut_te, left, weight, s) .+ te_offset[:, s]
              for s in 1:n_sections]
    return sec_le, sec_te
end

"""
    reconstruct_sections_sym(mode, wing, points, connectors, column)
        -> (sec_le, sec_te)

Live symbolic body-frame LE/TE of every refined section: the strut LE/TE connector
positions (`connectors.point_pos`, derived from world `pos`) interpolated by the
frozen mesh weights plus the frozen billow offset. Shared by all continuous VSM
modes so the force model and the plot use identical geometry.
"""
function reconstruct_sections_sym(mode, wing, points, connectors, column)
    n_struts = Int(wing.vsm_wing.n_unrefined_sections)
    strut_le = [collect(connectors.point_pos[:, column[(s, :LE)]])
                for s in 1:n_struts]
    strut_te = [collect(connectors.point_pos[:, column[(s, :TE)]])
                for s in 1:n_struts]
    left, weight, le_offset, te_offset = section_interp_caches(mode)
    return interp_sections(strut_le, strut_te, left, weight, le_offset, te_offset)
end

"""
    reconstruct_inflow_sym(mode, wing, connectors, column) -> (sec_va, sec_rho)

Live symbolic body-frame apparent wind and density of every refined section: each
strut's LE/TE connector values blended at [`COLLOCATION_CHORD_FRAC`](@ref) into one
strut value, then interpolated by the frozen mesh weights — the same struts and
weights [`reconstruct_sections_sym`](@ref) builds that section's corners from, so a
section's inflow and its geometry read the same points. Symbolic twin of
[`strut_inflow_weights`](@ref). Shared by all continuous VSM modes.
"""
function reconstruct_inflow_sym(mode, wing, connectors, column)
    n_struts = Int(wing.vsm_wing.n_unrefined_sections)
    nose = 1 - COLLOCATION_CHORD_FRAC
    strut_va = [nose * collect(connectors.va[:, column[(s, :LE)]]) +
                COLLOCATION_CHORD_FRAC *
                    collect(connectors.va[:, column[(s, :TE)]])
                for s in 1:n_struts]
    strut_rho = [nose * connectors.rho[column[(s, :LE)]] +
                 COLLOCATION_CHORD_FRAC * connectors.rho[column[(s, :TE)]]
                 for s in 1:n_struts]
    left, weight, _, _ = section_interp_caches(mode)
    sec_va = [interp_strut(strut_va, left, weight, s) for s in eachindex(left)]
    sec_rho = [interp_strut(strut_rho, left, weight, s) for s in eachindex(left)]
    return sec_va, sec_rho
end

"""
    reconstruct_sections_b(mode, wing, points) -> (sec_le, sec_te)

Numeric body-frame refined-section LE/TE from the live structural `pos_w`, exactly
as [`reconstruct_sections_sym`](@ref) builds them symbolically. For plotting and
logging so the drawn panels are the ones the dynamics use.
"""
function reconstruct_sections_b(mode, wing, points)
    column = aero_section_columns(wing, points)
    n_struts = Int(wing.vsm_wing.n_unrefined_sections)
    rot_w_to_b = (wing.R_b_to_w::Matrix{SimFloat})'
    strut_le = [rot_w_to_b * (points[column[(s, :LE)]].pos_w - wing.pos_w)
                for s in 1:n_struts]
    strut_te = [rot_w_to_b * (points[column[(s, :TE)]].pos_w - wing.pos_w)
                for s in 1:n_struts]
    left, weight, le_offset, te_offset = section_interp_caches(mode)
    return interp_sections(strut_le, strut_te, left, weight, le_offset, te_offset)
end

"""
    write_live_aero_log_points!(mode, wing, sys_struct, sys_state, point_idx, zoom)

Log the panel corners the force model reconstructs (strut interpolation + frozen
billow offset), not the raw VSM mesh, so the plot shows the deforming geometry the
dynamics use. Shared by the continuous VSM modes.
"""
function write_live_aero_log_points!(mode, wing, sys_struct, sys_state,
                                     point_idx, zoom)
    points = wing_points(sys_struct, wing)
    sec_le, sec_te = reconstruct_sections_b(mode, wing, points)
    rot_b_to_w = wing.R_b_to_w::Matrix{SimFloat}
    n_panels = length(sec_le) - 1
    for i in 1:n_panels
        for corner_b in (sec_le[i], sec_te[i], sec_te[i + 1], sec_le[i + 1])
            point_idx += 1
            corner_w = wing.pos_w + rot_b_to_w * corner_b
            sys_state.X[point_idx] = corner_w[1] * zoom
            sys_state.Y[point_idx] = corner_w[2] * zoom
            sys_state.Z[point_idx] = corner_w[3] * zoom
        end
    end
    return point_idx
end

# ==================== logging / visualization hooks ==================== #

"""
    n_aero_log_points(mode, wing) -> Int

Number of extra `SysState` log slots the mode contributes for `wing`
(visualization geometry such as panel corners). Default 0; VSM modes log
4 corners per panel. Must match what [`write_aero_log_points!`](@ref) writes.
"""
n_aero_log_points(::AbstractAeroModel, wing) = 0
n_aero_log_points(mode::AbstractVSMAero, wing) =
    4 * length(mode.vsm_aero.panels)

"""
    write_aero_log_points!(mode, wing, sys_struct, sys_state, point_idx,
                           zoom) -> Int

Write the mode's log points (world frame, scaled by `zoom`) into
`sys_state.X/Y/Z` starting after `point_idx`; return the last index written.
Default writes nothing; VSM modes write the panel corners, [`AeroPlate`](@ref)
writes each section's display quad. Flap deflections are logged separately, per
aero segment, into `sys_state.flap_angle` (see [`write_flap_deflections!`](@ref)).
"""
write_aero_log_points!(::AbstractAeroModel, wing, sys_struct, sys_state,
                       point_idx, zoom) = point_idx

function write_aero_log_points!(mode::AbstractVSMAero, wing, sys_struct,
                                sys_state, point_idx, zoom)
    R_b_to_w = wing.R_b_to_w::Matrix{SimFloat}
    for panel in mode.vsm_aero.panels
        for j in 1:4
            point_idx += 1
            corner_w = wing.pos_w + R_b_to_w * panel.corner_points[:, j]
            sys_state.X[point_idx] = corner_w[1] * zoom
            sys_state.Y[point_idx] = corner_w[2] * zoom
            sys_state.Z[point_idx] = corner_w[3] * zoom
        end
    end
    return point_idx
end

"""
    read_aero_log_points!(mode, wing, sys_struct, sys_state, point_idx) -> Int

Inverse of [`write_aero_log_points!`](@ref): restore the mode's state from the
logged points starting after `point_idx`; return the last index consumed (the
slots must be skipped even when unused). Default consumes nothing; VSM
`PARTICLE_DYNAMICS` modes read the panel corners back (rigid wings recompute
panels from twist instead and only skip their slots). Flap deflections are
restored separately from `sys_state.flap_angle` (see [`restore_flap_delta!`](@ref)).
"""
read_aero_log_points!(::AbstractAeroModel, wing, sys_struct, sys_state,
                      point_idx) = point_idx

function read_aero_log_points!(mode::AbstractVSMAero, wing, sys_struct,
                               sys_state, point_idx)
    n_corners = 4 * length(mode.vsm_aero.panels)
    wing.dynamics_type == RIGID_DYNAMICS && return point_idx + n_corners
    R_w_to_b = (wing.R_b_to_w::Matrix{SimFloat})'
    for panel in mode.vsm_aero.panels
        for j in 1:4
            point_idx += 1
            corner_w = [sys_state.X[point_idx], sys_state.Y[point_idx],
                        sys_state.Z[point_idx]]
            panel.corner_points[:, j] .= R_w_to_b * (corner_w - wing.pos_w)
        end
    end
    return point_idx
end

"""
    restore_aero_twist!(mode, wing, twist_surfaces)

Re-apply the (already restored) twist-surface angles to the mode's geometry
when loading a `SysState` log frame. Default no-op; VSM `RIGID_DYNAMICS`
modes deform the unrefined sections and reinit the panels.
"""
restore_aero_twist!(::AbstractAeroModel, wing, twist_surfaces) = nothing

function restore_aero_twist!(mode::AbstractVSMAero, wing, twist_surfaces)
    wing.dynamics_type == RIGID_DYNAMICS || return nothing
    isempty(wing.twist_surface_idxs) && return nothing
    vsm = mode.vsm_wing
    isempty(vsm.non_deformed_sections) && return nothing
    theta = zeros(Float64, vsm.n_unrefined_sections)
    for twist_surface_idx in wing.twist_surface_idxs
        for section_idx in
                twist_surfaces[twist_surface_idx].unrefined_section_idxs
            theta[section_idx] = twist_surfaces[twist_surface_idx].twist
        end
    end
    VortexStepMethod.unrefined_deform!(vsm, theta)
    VortexStepMethod.reinit!(mode.vsm_aero; init_aero=false)
    return nothing
end

"""
    n_flap_deflections(sys_struct) -> Int

Number of aero segments logged in `SysState.flap_angle`: one flap deflection δ
per twist_surface. Sets the `D` type parameter of the model's `SysState`.
"""
n_flap_deflections(sys_struct) = length(sys_struct.twist_surfaces)

"""
    write_flap_deflections!(sys_state, sys_struct)

Write each twist_surface's flap deflection δ [rad] into `sys_state.flap_angle`
(indexed by twist_surface `idx`, via [`twist_surface_deltas`](@ref)). No-op when
the state carries no flap slots (`D == 0`).
"""
function write_flap_deflections!(sys_state, sys_struct)
    isempty(sys_state.flap_angle) && return nothing
    sys_state.flap_angle .= twist_surface_deltas(sys_struct)
    return nothing
end

"""
    restore_flap_delta!(mode, wing, sys_state)

Restore each VSM panel's flap deflection `δ` from `sys_state.flap_angle` (mapped
through the panel→twist_surface index) when loading a `SysState` frame, so a
replayed frame shows the logged flap state. Default no-op; [`AeroPressure`](@ref)
restores its `PARTICLE_DYNAMICS` panels.
"""
restore_flap_delta!(::AbstractAeroModel, wing, sys_state) = nothing

# ==================== shared VSM numerics ==================== #

"""
    finite_full(x) -> Bool

`true` if `x` is finite. For a `ForwardDiff.Dual` it also checks every partial, so
a NaN/Inf *derivative* is caught — used by [`safe_vsm_solve!`](@ref) to reject a
bad VSM solve during the Jacobian pass, not just the value pass.
"""
finite_full(x::Real) = isfinite(x)
finite_full(x::ForwardDiff.Dual) =
    isfinite(ForwardDiff.value(x)) &&
    all(isfinite, ForwardDiff.partials(x))

"""
    VSMSolveFailure(msg)

A VSM solve that did not converge or returned a non-finite result, thrown by
every VSM mode's refresh. `vsm_warn_on_fail` ([`refresh_aero!`](@ref),
[`next_step!`](@ref)) downgrades it to a warning; the assertions on a corrupted
frozen state are `AssertionError` and stay fatal.
"""
struct VSMSolveFailure <: Exception
    msg::String
end

Base.showerror(io::IO, failure::VSMSolveFailure) = print(io, failure.msg)

"""
NaN/Inf-guarded `solve!`. Checks both Dual value and partials. On a non-finite or
non-converged result, restore the aero state of the last converged solve — its
circulation and its two angle-of-attack distributions, which `solve!` overwrites
with the diverged ones — and return `false`.

Without `gamma_init`, `solve!` warm-starts from the circulation it left in
`solver.sol` last time. `cold_start` starts from the solver's configured initial
distribution instead: past stall the iteration has more than one fixed point, so
a warm start makes the answer a function of what ran before rather than of the
current state.
"""
function safe_vsm_solve!(solver, body_aero,
                          gamma_init=nothing; moment_frac=0.1, cold_start=false)
    gamma_converged = isnothing(solver.sol.gamma_distribution) ? nothing :
        copy(solver.sol.gamma_distribution)
    alpha_converged = copy(solver.lr.alpha_dist)
    alpha_corrected_converged = copy(solver.sol.alpha_dist)
    if cold_start
        VortexStepMethod.solve!(solver, body_aero, nothing;
            moment_frac, log=false)
    elseif isnothing(gamma_init)
        VortexStepMethod.solve!(solver, body_aero;
            moment_frac, log=false)
    else
        VortexStepMethod.solve!(solver, body_aero, gamma_init;
            moment_frac, log=false)
    end
    force_coeffs = solver.sol.force_coeffs
    moment_coeffs = solver.sol.moment_coeffs
    if !solver.lr.converged ||
            any(!finite_full, force_coeffs) ||
            any(!finite_full, moment_coeffs)
        isnothing(gamma_converged) ||
            copyto!(solver.sol.gamma_distribution, gamma_converged)
        copyto!(solver.lr.alpha_dist, alpha_converged)
        copyto!(solver.sol.alpha_dist, alpha_corrected_converged)
        return false
    end
    return true
end

"""
    solve_and_freeze_circulation!(mode, wing; cold_start=false)

Solve and freeze the per-refined-panel induced velocity into `mode.v_ind` and the
chord blend weights into `mode.chord_weight`, shared by the continuous VSM modes.
Warm starting is `VortexStepMethod.solve!`'s own, gated by `use_gamma_prev` and by
`cold_start` ([`safe_vsm_solve!`](@ref)). Throws a [`VSMSolveFailure`](@ref) on a
non-converged or non-finite solve, before anything frozen is written.
"""
function solve_and_freeze_circulation!(mode, wing; cold_start=false)
    solver = wing.vsm_solver
    body_aero = wing.vsm_aero
    if !safe_vsm_solve!(solver, body_aero; cold_start)
        throw(VSMSolveFailure("$(nameof(typeof(mode))) VSM solve failed " *
            "(non-converged or non-finite) on wing $(wing.idx)."))
    end
    store_induced_velocity!(mode.v_ind, body_aero, solver.lr.gamma_new)
    store_chord_weights!(mode.chord_weight, body_aero)
    return nothing
end

"""
    set_particle_panel_va!(wing, va_point_b_vals)

Set the per-panel apparent wind on the wing's VSM `BodyAerodynamics` from the
per-point body-frame apparent wind: each section's `va` is the mean of its
LE/TE point values, refined panels take their parent section's value. Falls
back to the wing-level `va_b` when the structural↔panel mapping is missing.
"""
function set_particle_panel_va!(wing, va_point_b_vals)
    if isnothing(wing.point_to_vsm_point)
        set_va!(wing.vsm_aero, wing.va_b)
        return nothing
    end
    n_sections = length(wing.vsm_wing.unrefined_sections)
    section_va = Vector{Vector{Float64}}(undef, n_sections)

    vsm_point_to_struct = Dict{Tuple{Int64, Symbol}, Int64}()
    for (point_idx, (section_idx, le_or_te)) in wing.point_to_vsm_point
        vsm_point_to_struct[(section_idx, le_or_te)] = point_idx
    end

    for section_idx in 1:n_sections
        le_pi = get(vsm_point_to_struct, (Int64(section_idx), :LE), nothing)
        te_pi = get(vsm_point_to_struct, (Int64(section_idx), :TE), nothing)
        if !isnothing(le_pi) && !isnothing(te_pi)
            va_le = va_point_b_vals[:, le_pi]
            va_te = va_point_b_vals[:, te_pi]
            section_va[section_idx] = 0.5 * (va_le + va_te)
        else
            section_va[section_idx] = wing.va_b
        end
    end

    n_panels = length(wing.vsm_aero.panels)
    va_dist = zeros(n_panels, 3)
    mapping = wing.vsm_wing.refined_panel_mapping
    for rpi in 1:n_panels
        va_dist[rpi, :] .= section_va[mapping[rpi]]
    end
    set_va!(wing.vsm_aero, va_dist)
    return nothing
end

"""
    set_refined_panel_va!(mode, wing, points, va_point_b_vals)

Per-panel apparent wind for the continuous VSM modes: refined sections from
[`strut_inflow_weights`](@ref) (numeric twin of [`reconstruct_inflow_sym`](@ref)),
each panel the mean of its two bounding sections. Matches the inflow the symbolic
RHS uses, so the frozen circulation belongs to it. Falls back to `va_b` without a
structural↔panel mapping.
"""
function set_refined_panel_va!(mode, wing, points, va_point_b_vals)
    if isnothing(wing.point_to_vsm_point)
        set_va!(wing.vsm_aero, wing.va_b)
        return nothing
    end
    column = aero_section_columns(wing, points)
    section_va = [
        sum(weight .* view(va_point_b_vals, :, points[col].idx)
            for (col, weight) in strut_inflow_weights(mode, section, column))
        for section in eachindex(section_interp_caches(mode)[1])
    ]
    n_panels = length(wing.vsm_aero.panels)
    va_dist = zeros(SimFloat, n_panels, 3)
    for panel in 1:n_panels
        va_dist[panel, :] .= 0.5 .* (section_va[panel] .+ section_va[panel + 1])
    end
    set_va!(wing.vsm_aero, va_dist)
    return nothing
end

"""
    vsm_solve_objects(wing, ::Type{T}, shadow_ref) -> (body_aero, solver, wing)

The VSM solve objects [`vsm_aero_coeffs`](@ref) runs on, selected by the input
eltype `T`. A value pass (`Float64`) uses the wing's real objects. A `ForwardDiff`
pass feeds `Dual` numbers through the solve to get the Jacobian, but the real
objects' buffers are `Float64` and can't hold Duals — so for a `Dual` eltype we
solve on a `Dual`-typed "shadow" of the solver/aero. The shadow is expensive, so
it is built lazily and cached in `shadow_ref`, keyed by the `Dual` eltype (rebuilt
if it changes); `use_gamma_prev` warm-starts each perturbed solve from the
previous circulation.
"""
vsm_solve_objects(wing, ::Type{Float64}, shadow_ref) =
    (wing.vsm_aero, wing.vsm_solver, wing.vsm_wing)

function vsm_solve_objects(wing, ::Type{T}, shadow_ref) where {T}
    shadow = shadow_ref[]
    if shadow === nothing || eltype(shadow[1]._va) !== T
        shadow = VortexStepMethod.make_dual_shadow(
            wing.vsm_solver, wing.vsm_aero, T)
        shadow[2].use_gamma_prev = true
        shadow_ref[] = shadow
    end
    body_aero, solver = shadow
    return body_aero, solver, body_aero.wings[1]
end

"""
    vsm_aero_coeffs(wing, y, va_mag, n_unrefined, n_twist_surfaces,
                     twist_surface_idxs, twist_surfaces, moment_frac, shadow_ref;
                     gamma_init=nothing) -> Vector

Run one VSM solve at operating-point input `y = [α, β, ω₁, ω₂, ω₃, θ_twist…]` and
return the wind-axis coefficient vector `[CL, CD, CS, CM₁, CM₂, CM₃, cm_twist…]`.
`ForwardDiff.Dual`-aware via `vsm_solve_objects`: for a Dual eltype it solves on a
cached dual shadow of the VSM solver, so the same routine yields the Jacobian
under AD.
"""
function vsm_aero_coeffs(wing, y::AbstractVector{T},
        va_mag, n_unrefined, n_twist_surfaces,
        twist_surface_idxs, twist_surfaces, moment_frac,
        shadow_ref::Ref;
        gamma_init=nothing) where {T}

    body_aero_c, solver_c, wing_c = vsm_solve_objects(wing, T, shadow_ref)

    α = y[1]
    β = y[2]
    ω = MVector{3, T}(y[3], y[4], y[5])

    # Body-frame apparent wind from (α, β, va_mag)
    cα, sα = cos(α), sin(α)
    cβ, sβ = cos(β), sin(β)
    va_b_local = MVector{3, T}(va_mag * cα * cβ,
                               va_mag * sβ,
                               va_mag * sα * cβ)

    # Per-twist_surface → per-section twist
    theta = zeros(T, n_unrefined)
    for (twist_surface_index, gidx) in enumerate(twist_surface_idxs)
        for unrefined_index in twist_surfaces[gidx].unrefined_section_idxs
            theta[unrefined_index] = y[5 + twist_surface_index]
        end
    end

    if n_unrefined > 0
        VortexStepMethod.unrefined_deform!(
            wing_c, theta; smooth=false)
        VortexStepMethod.reinit!(
            body_aero_c; init_aero=false)
    end
    set_va!(body_aero_c, va_b_local, ω)
    if !safe_vsm_solve!(solver_c, body_aero_c, gamma_init;
                         moment_frac)
        throw(VSMSolveFailure("VSM solve failed (non-converged or " *
            "non-finite) on wing $(wing.idx) [eltype=$T]."))
    end

    sol = solver_c.sol
    force_coeffs = sol.force_coeffs
    cm_body = sol.moment_coeffs
    moment_coeff_unrefined = sol.moment_coeff_unrefined_dist

    # Wind-axis basis (VSM): drag∥va, lift=norm(drag×span), side=lift×drag.
    span = SVector(zero(T), one(T), zero(T))
    drag_dir = va_b_local ./ va_mag
    lift_dir = smooth_normalize(cross(drag_dir, span))
    side_dir = cross(lift_dir, drag_dir)

    x = zeros(T, 6 + n_twist_surfaces)
    x[1] = dot(force_coeffs, lift_dir)
    x[2] = dot(force_coeffs, drag_dir)
    x[3] = dot(force_coeffs, side_dir)
    x[4] = cm_body[1]
    x[5] = cm_body[2]
    x[6] = cm_body[3]
    for (twist_surface_index, gidx) in enumerate(twist_surface_idxs)
        x[6 + twist_surface_index] = sum(
            moment_coeff_unrefined[unrefined_index]
            for unrefined_index in
                twist_surfaces[gidx].unrefined_section_idxs;
            init = zero(T))
    end
    return x
end

"""
    rigid_aero_baseline!(wing, twist_surfaces; vsm_min_wind=0.5)

Compute the operating point and baseline wind-axis coefficients for one wing:
writes `wing.aero_y` / `wing.aero_x` and updates `twist_surfaces[gidx].aero_moment`.
Returns the context (`va_mag`, section counts, `moment_frac`, `shadow_ref`, `y0`)
the mode-specific reduction (`refresh_rigid_aero!`) needs for the Jacobian.
"""
function rigid_aero_baseline!(wing, twist_surfaces;
                              vsm_min_wind=0.5)
    va_b = wing.va_b
    va_mag_actual = norm(va_b)
    omega_b = wing.ω_b

    twist_surface_idxs = wing.twist_surface_idxs
    n_twist_surfaces = length(twist_surface_idxs)
    n_unrefined = wing.vsm_wing.n_unrefined_sections

    moment_frac = isempty(twist_surface_idxs) ? 0.25 :
        twist_surfaces[first(twist_surface_idxs)].moment_frac

    va_mag = max(va_mag_actual, vsm_min_wind)
    alpha_0 = atan(va_b[3], va_b[1])
    beta_0 = atan(va_b[2], hypot(va_b[1], va_b[3]))
    if !isfinite(alpha_0)
        alpha_0 = 0.0
    end
    if !isfinite(beta_0)
        beta_0 = 0.0
    end

    # Operating-point input vector y₀ = [α, β, ω, θ_twist_surface]
    y0 = wing.aero_y
    y0[1] = alpha_0
    y0[2] = beta_0
    y0[3] = omega_b[1]
    y0[4] = omega_b[2]
    y0[5] = omega_b[3]
    for (twist_surface_index, gidx) in enumerate(twist_surface_idxs)
        y0[5 + twist_surface_index] = twist_surfaces[gidx].twist
    end

    shadow_ref = Ref{Any}(nothing)
    f_baseline = y -> vsm_aero_coeffs(wing, y, va_mag,
        n_unrefined, n_twist_surfaces, twist_surface_idxs, twist_surfaces,
        moment_frac, shadow_ref)

    wing.aero_x .= f_baseline(y0)
    for (twist_surface_index, gidx) in enumerate(twist_surface_idxs)
        twist_surfaces[gidx].aero_moment = wing.aero_x[6 + twist_surface_index]
    end

    return (; va_mag, n_unrefined, n_twist_surfaces,
            twist_surface_idxs, moment_frac, shadow_ref, y0)
end
