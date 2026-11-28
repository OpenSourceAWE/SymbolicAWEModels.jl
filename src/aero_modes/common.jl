# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Shared aero-mode code: dispatch interface, connector scaffolding, VSM refresh.

# ==================== interface: capability traits ==================== #

"""
    vsm_engine(mode::AbstractAeroModel) -> Union{Nothing, VSMEngine}

The mode's [`VSMEngine`](@ref) (VSM geometry + linearization state), or
`nothing` for modes without one. Used for the wing's VSM property forwarding and
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

`true` if the mode needs per-section stations (auto-creation and
aero-section matching). VSM modes ([`AbstractVSMAero`](@ref)) do.
"""
couples_to_sections(::AbstractAeroModel) = false
couples_to_sections(::AbstractVSMAero) = true

"""
    has_vsm_wing(sam) -> Bool

`true` if any wing of the run `sam` carries a VSM aero mode ([`AbstractVSMAero`](@ref))
whose state must be refreshed each `vsm_interval`. When `false`,
[`refresh_aero!`](@ref) is a no-op for every wing, so callers skip it and the trailing
parameter sync.
"""
has_vsm_wing(sam) = any(has_vsm_engine, sam.aero)

# ==================== interface: required cache tag ==================== #

"""
    aero_mode_tag(mode::AbstractAeroModel) -> String

Short identifier for the mode in the compiled-model cache filename, the same for a mode
and for the [`AeroSpec`](@ref) a structure holds of it. Required: an aero mode declares
its tag on its instance or its type, else this errors, so no silent default collides two
distinct modes on one cache file. Built-ins: `"lin"`, `"dir"`, `"cont"`, `"press"`,
`"none"`, `"plate"`.
"""
aero_mode_tag(mode::AbstractAeroModel) = aero_mode_tag(typeof(mode))
aero_mode_tag(::AeroSpec{M}) where {M} = aero_mode_tag(M)
aero_mode_tag(M::Type) = error(
    "aero_mode_tag is not defined for aero mode $M; " *
    "every aero mode must provide its own cache tag.")

# ==================== interface: diagnostics ==================== #

"""
    calc_aoa(mode::AbstractAeroModel, sam, wing_idx) -> SimFloat

Angle of attack [rad] of wing `wing_idx` of `sam` under its aero `mode`. Defaults to `NaN`
(undefined); VSM modes read the mid-span geometric AoA (wrapped to [-π, π]) and
[`AeroPlate`](@ref) derives it from the body-frame apparent wind.
"""
calc_aoa(::AbstractAeroModel, sam, wing_idx) = SimFloat(NaN)
function calc_aoa(mode::AbstractVSMAero, sam, wing_idx)
    dist = mode.vsm_solver.sol.alpha_geometric_dist
    n = length(dist)
    return mod(dist[n ÷ 2 + n % 2] + π, 2π) - π
end

"""The apparent wind on the body of wing `wing_idx` of `sam`, body frame [m/s]."""
wing_apparent_wind(sam, wing_idx) =
    run_value(sam, Val(:bodies), wing_body(sam, wing_idx), Val(:va_b))

"""
    calc_side_slip(sam, wing_idx) -> SimFloat

Side-slip angle [rad] of wing `wing_idx` of `sam` from its body-frame apparent wind, the
same for every aero mode.
"""
function calc_side_slip(sam, wing_idx)
    va_b = wing_apparent_wind(sam, wing_idx)
    return atan(va_b[2], hypot(va_b[1], va_b[3]))
end

# ==================== connector scaffolding ==================== #

function rigid_aero_connectors(num_stations::Int)
    @variables begin
        va(t)[1:3]
        rho(t)
        R_b_w(t)[1:3, 1:3]
        omega(t)[1:3]
        force(t)[1:3]
        moment(t)[1:3]
    end
    if num_stations > 0
        @variables begin
            twist(t)[1:num_stations]
            twist_vel(t)[1:num_stations]
            twist_moment(t)[1:num_stations]
        end
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
    frozen_point_force_component(sam, wing_idx; name, params) -> System

Particle aero component binding each node of wing `wing_idx` of `sam` to its frozen
`aero_force_b` flat parameter, set at each refresh ([`sync_aero_params!`](@ref)). Used by
[`AeroDirect`](@ref), whose refresh precomputes that force. Touching the parameter
here is what registers it.
"""
function frozen_point_force_component(sam, wing_idx; name, params=nothing)
    point_idxs = wing_points(sam, wing_idx)
    connectors = particle_aero_connectors(length(point_idxs))
    flat_ps = Any[]
    eqs = Equation[]
    for (point_num, i) in enumerate(point_idxs)
        force_p = params.points[i].aero_force_b
        push!(flat_ps, force_p)
        append!(eqs, collect(connectors.point_force[:, point_num]) .~ collect(force_p))
    end
    return System(eqs, t, particle_unknowns(connectors), flat_ps; name)
end

"""
    wing_points(sam, wing_idx) -> Vector{Int}

The indices of the nodes of wing `wing_idx` of `sam`, in point order: the connector columns
of its particle aero component.
"""
wing_points(sam, wing_idx) = [i for i in eachindex(sam.definition.points)
                       if wing_node(sam.definition, i) && point_wing(sam, i) == wing_idx]

"""
    build_panel_force_eqs(sec_le, sec_te, sec_va, sec_rho, vind_p, chord_w,
                          cl, cd, cm, spanwise, scale)
        -> (eqs, vars, panel_force, panel_couple, curvature_couple, slots)

Shared per-refined-panel VSM force assembly for the live particle aero modes
([`ContinuousAero`](@ref), [`AeroPressure`](@ref)). Traces the `VortexStepMethod`
panel aerodynamics symbolically from the frozen circulation: per panel
`i` (between section boundaries `i` and `i+1`) it builds the airfoil axes, chord,
width, effective angle of attack (live apparent wind `sec_va` + frozen induced
velocity `vind_p`), polar coefficients (`cl/cd/cm` callable params) and the lift/drag
directions, and emits the panel force and pitching-moment couple. Returns the panel
equations, the intermediate variables to register, the `panel_force`/`panel_couple`
arrays for the caller's scatter, and the whole [`panel_force_slots`](@ref) named
tuple for a scatter that needs the panel axes too. `curvature_couple` is the
`VortexStepMethod.flow_curvature_cm` part of `panel_couple` alone.

`sec_le`/`sec_te`/`sec_va` are length-`n_panels+1` vectors of body-frame
3-vectors (positions and apparent wind at the section boundaries), `sec_rho`
the matching air densities; they may be live (interpolated from structure) or
constant (a fixed mesh). `spanwise` is the wing spanwise direction, `scale` the
chord-scale factor and `chord_w` the per-panel chord blend weight from
[`store_chord_weights!`](@ref).

`delta` is an optional length-`n_panels` vector of symbolic per-panel flap
deflections δ. When given the polars are evaluated as `cl(i, alpha[i], delta[i])`
(the `(α, δ)` tables); when `nothing` the 2-arg `cl(i, alpha[i])` is used, so a
mode without a flap ([`ContinuousAero`](@ref)) is untouched.

`sec_dva` is the matching per-section trailing minus leading edge apparent wind that
carries the `VortexStepMethod.flow_curvature_cm` increment, or `nothing` to leave it out.

`deficiency` is the wing's Wagner lag ([`wagner_lag_eqs`](@ref)), subtracted from
every panel's angle of attack before the polars are read; `0.0` leaves them steady.
"""
function build_panel_force_eqs(sec_le, sec_te, sec_va, sec_rho,
                               vind_p, chord_w, cl, cd, cm, spanwise, scale;
                               delta=nothing, sec_dva=nothing, deficiency=0.0)
    n_panels = length(sec_le) - 1
    slots = panel_force_slots(n_panels)
    eqs = Equation[]
    for i in 1:n_panels
        append!(eqs, panel_force_eqs(slots, i,
            (sec_le[i], sec_te[i], sec_le[i + 1], sec_te[i + 1]),
            (sec_va[i], sec_va[i + 1], sec_rho[i], sec_rho[i + 1],
             [vind_p[c, i] for c in 1:3],
             sec_dva === nothing ? nothing : sec_dva[i],
             sec_dva === nothing ? nothing : sec_dva[i + 1]),
            (PanelPolar(cl, i), PanelPolar(cd, i), PanelPolar(cm, i)),
            spanwise, scale, chord_w[i],
            delta === nothing ? nothing : delta[i], deficiency))
    end

    return eqs, panel_force_vars(slots), slots.panel_force, slots.panel_couple,
           slots.curvature_couple, slots
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
        alpha_eff(t)[1:n_panels]
        q_dyn(t)[1:n_panels]
        pitch_rate(t)[1:n_panels]
        dir_lift(t)[1:3, 1:n_panels]
        dir_drag(t)[1:3, 1:n_panels]
        panel_force(t)[1:3, 1:n_panels]
        panel_couple(t)[1:3, 1:n_panels]
        curvature_couple(t)[1:3, 1:n_panels]
    end
    return (; x_airf, y_airf, z_airf, v_eff, chord, width, alpha, alpha_eff, q_dyn,
            pitch_rate, dir_lift, dir_drag, panel_force, panel_couple,
            curvature_couple)
end

"""The slots of [`panel_force_slots`](@ref) as the variable list a `System` declares."""
panel_force_vars(slots) = Any[values(slots)...]

"""
    panel_force_eqs(slots, i, sections, flow, polars, spanwise, scale,
                    chord_weight, delta)

One panel's aerodynamic equations, writing into column `i` of the symbolic arrays
in `slots`. `sections` is `(le_1, te_1, le_2, te_2)` in body frame, `flow` is
`(va_1, va_2, rho_1, rho_2, v_ind, dva_1, dva_2)` and `polars` the `(cl, cd, cm)`
callables, indexed by the panel number the polar tables were built for.

The physics is not written here: every expression comes from tracing the
`VortexStepMethod` panel aerodynamics (`panel_axes`, `panel_inflow`,
`panel_force_directions`, `panel_loads` and friends) with symbolic arguments, so
the equations are the ones the numeric solver evaluates. All this function
decides is where to tear the expression graph, binding each stage to a slot
before feeding it to the next.

`chord_weight` is the section-1 share of the chord-direction blend
([`store_chord_weights!`](@ref)). `delta` is the flap deflection or `nothing` for
the 2-argument polars. `deficiency` is the wing-wide Wagner lag subtracted from
`alpha` to give the `alpha_eff` the polars read; the geometric `alpha` still sets
the force directions.

`dva_1`/`dva_2` are the sections' trailing minus leading edge apparent wind, giving
the panel the pitch rate its `VortexStepMethod.flow_curvature_cm` moment increment
is read at. They are `nothing` when the wing's solver has the term disabled
([`flow_curvature_enabled`](@ref)), and the rate is then bound to zero rather than
dropped, so the slot means the same thing either way.

Every quantity it reads belongs to this panel alone, so the same equations serve a
whole-wing system (looped by [`build_panel_force_eqs`](@ref)) and a per-panel
component compiled once and instantiated for each panel.
"""
function panel_force_eqs(slots, i, sections, flow, polars, spanwise, scale,
                         chord_weight, delta, deficiency=0.0)
    (; x_airf, y_airf, z_airf, v_eff, chord, width, alpha, alpha_eff, q_dyn,
       pitch_rate, dir_lift, dir_drag, panel_force, panel_couple,
       curvature_couple) = slots
    le_1, te_1, le_2, te_2 = sections
    va_1, va_2, rho_1, rho_2, vind, dva_1, dva_2 = flow
    cl, cd, cm = polars

    axes = VortexStepMethod.panel_axes(le_1, te_1, le_2, te_2, chord_weight)
    inflow = VortexStepMethod.panel_inflow(axes, va_1, va_2, vind)
    bound_axes = (x_airf=collect(x_airf[:, i]), y_airf=collect(y_airf[:, i]),
                  z_airf=collect(z_airf[:, i]), chord=chord[i], width=width[i])
    dirs = VortexStepMethod.panel_force_directions(bound_axes, alpha[i], spanwise)
    bound_dirs = (dir_lift=collect(dir_lift[:, i]),
                  dir_drag=collect(dir_drag[:, i]))

    pitch_rate_rhs = dva_1 === nothing ? 0.0 :
        VortexStepMethod.section_pitch_rate(0.5 .* (dva_1 .+ dva_2),
                                            bound_axes.z_airf, chord[i])
    curvature = dva_1 === nothing ? 0.0 :
        VortexStepMethod.flow_curvature_cm(pitch_rate[i], chord[i],
                                           smooth_norm(inflow.v_span))
    moment_coeff = evaluate_polar(cm, alpha_eff[i], delta) + curvature
    loads = VortexStepMethod.panel_loads(bound_axes, bound_dirs, q_dyn[i],
        evaluate_polar(cl, alpha_eff[i], delta),
        evaluate_polar(cd, alpha_eff[i], delta), moment_coeff, scale)

    return [
        chord[i] ~ axes.chord;
        width[i] ~ axes.width;
        x_airf[:, i] ~ axes.x_airf;
        y_airf[:, i] ~ axes.y_airf;
        z_airf[:, i] ~ axes.z_airf;
        v_eff[:, i] ~ inflow.v_eff;
        alpha[i] ~ inflow.alpha;
        alpha_eff[i] ~ VortexStepMethod.effective_alpha(alpha[i], deficiency);
        q_dyn[i] ~ VortexStepMethod.dynamic_pressure(rho_1, rho_2, inflow.v_span);
        pitch_rate[i] ~ pitch_rate_rhs;
        dir_lift[:, i] ~ dirs.dir_lift;
        dir_drag[:, i] ~ dirs.dir_drag;
        panel_force[:, i] ~ loads.force;
        panel_couple[:, i] ~ VortexStepMethod.panel_couple_force(
            moment_coeff, q_dyn[i], chord[i], width[i], scale) .*
            bound_axes.z_airf;
        curvature_couple[:, i] ~ VortexStepMethod.panel_couple_force(
            curvature, q_dyn[i], chord[i], width[i], scale) .* bound_axes.z_airf]
end

"""
    scatter_couple(mode, slots, i, panel) -> Vector

The couple panel `i`'s scatter places, from the [`panel_force_slots`](@ref) the panel
equations wrote. The default is the whole `panel_couple`, which is what a scatter that
carries no moment of its own needs. [`AeroPressure`](@ref) overrides it: its frozen
traction already carries most of the moment, so it places only the deficit, see
[`pressure_couple`](@ref). `panel` addresses that mode's frozen per-panel params and is
`nothing` for a mode that has none.
"""
scatter_couple(::AbstractAeroModel, slots, i, panel) =
    collect(slots.panel_couple[:, i])

"""
    flow_curvature_enabled(mode) -> Bool

Whether the VSM solver of aero `mode` adds the `VortexStepMethod.flow_curvature_cm`
increment. The aero modes all read the one flag, so a model cannot carry the term through
`AeroDirect` and lose it in [`ContinuousAero`](@ref). A mode without a VSM engine
([`AeroPlate`](@ref)) has no flag and no increment.
"""
function flow_curvature_enabled(mode)
    engine = vsm_engine(mode)
    return !isnothing(engine) && engine.vsm_solver.flow_curvature
end

# ==================== unsteady corrections ==================== #

"""
    unsteady_aero(mode) -> UnsteadyAero

The [`UnsteadyAero`](@ref) corrections of aero `mode`, or an all-off default for a mode
that carries no VSM engine.
"""
function unsteady_aero(mode)
    engine = vsm_engine(mode)
    return isnothing(engine) ? UnsteadyAero() : engine.unsteady
end

"""
    wagner_enabled(mode) -> Bool

Whether a wing flown with aero `mode` carries the two-state Wagner lift lag
([`wagner_lag_eqs`](@ref)).
"""
wagner_enabled(mode) = unsteady_aero(mode).wagner

"""
    panel_chord_width(mode) -> (chord, width)

Each refined panel's mean chord and quarter-chord width [m] read off the frozen VSM
mesh of aero `mode`, from the `VortexStepMethod` helpers the `chord` and `width`
equations in [`panel_force_eqs`](@ref) are traced through.
"""
function panel_chord_width(mode)
    refined = mode.vsm_wing.refined_sections
    n = Int(mode.vsm_wing.n_panels)
    chord = zeros(SimFloat, n)
    width = zeros(SimFloat, n)
    for i in 1:n
        section = (refined[i].LE_point, refined[i].TE_point,
                   refined[i + 1].LE_point, refined[i + 1].TE_point)
        chord[i] = VortexStepMethod.panel_chord(section...)
        width[i] = smooth_norm(VortexStepMethod.panel_span_vector(section...))
    end
    return chord, width
end

"""
    chord_blend_weights!(chord_weight, width)

Every panel's section-1 share of the chord-direction blend, written into
`chord_weight` (n_panels): `VortexStepMethod.panel_chord_weight` over the panel
widths `width`, with `nothing` standing in for the neighbour a tip panel lacks.
"""
function chord_blend_weights!(chord_weight, width)
    n_panels = length(width)
    for i in 1:n_panels
        chord_weight[i] = VortexStepMethod.panel_chord_weight(
            i == 1 ? nothing : width[i - 1], width[i],
            i == n_panels ? nothing : width[i + 1])
    end
    return chord_weight
end

"""
    panel_apparent_mass(mode, rho) -> Vector{SimFloat}

The entrained air [kg] of each refined panel of aero `mode` at air density `rho`: the
thin-plate added mass per unit span `ρ·π·(c/2)²` taken over the panel's width, which is
the fluid a flat plate carries with it when it accelerates normal to itself.
"""
function panel_apparent_mass(mode, rho)
    chord, width = panel_chord_width(mode)
    return @. rho * π * chord^2 * width / 4
end

"""
    apply_apparent_mass!(sam)

Give the nodes of each wing of the running `sam` the air they entrain, scaled by its
[`UnsteadyAero`](@ref) `apparent_mass`, as the `apparent_mass` of whatever integrates
them ([`apparent_mass_carriers`](@ref)). Each panel's [`panel_apparent_mass`](@ref), at
the air density at the wing, is spread over the nodes by the weights that carry its
force ([`aero_scatter_entries`](@ref)). A carrier of no entrained air gets zero. Writes
only values that change, and nothing before `init!` builds the integrator.
"""
function apply_apparent_mass!(sam)
    isnothing(sam.integrator) && return nothing
    totals = Dict{Tuple{Symbol, Int}, SimFloat}()
    for wing_idx in eachindex(sam.aero)
        add_apparent_mass!(totals, sam, wing_idx)
    end
    for ((collection, idx), mass) in totals
        run_value(sam, collection, idx, :apparent_mass) == mass && continue
        set_run_value!(sam, collection, idx, :apparent_mass, mass)
    end
    return nothing
end

"""
    add_apparent_mass!(totals, sam, wing_idx)

Add the air the nodes of wing `wing_idx` of `sam` entrain to `totals`, by carrier
`(collection, index)`, every carrier of its nodes starting at zero.
"""
function add_apparent_mass!(totals, sam, wing_idx)
    point_idxs = wing_points(sam, wing_idx)
    for i in point_idxs, (carrier, _) in apparent_mass_carriers(sam, i)
        get!(totals, carrier, 0.0)
    end
    mode = sam.aero[wing_idx]
    scale = unsteady_aero(mode).apparent_mass
    scale == 0 && return nothing
    name = sam.definition.wings[wing_idx].name
    body_idx = wing_body(sam, wing_idx)
    body_dynamics_type(sam, body_idx) == PARTICLE_DYNAMICS || error(
        "Wing $name: apparent_mass acts on the nodes' own acceleration, " *
        "which a RIGID_DYNAMICS wing does not integrate.")
    supports_panel_decomposition(mode) || error(
        "Wing $name: apparent_mass needs an aero mode that scatters per " *
        "panel; $(typeof(mode)) does not.")
    height = run_value(sam, Val(:bodies), body_idx, Val(:pos_w))[3]
    mass = panel_apparent_mass(mode, air_density(sam.am, height))
    for (panel, node, force_weight, _) in aero_scatter_entries(mode, sam, wing_idx)
        node_mass = (scale * force_weight) * mass[panel]
        for (carrier, share) in apparent_mass_carriers(sam, point_idxs[node])
            totals[carrier] += share * node_mass
        end
    end
    return nothing
end

"""
    apparent_mass_carriers(sam, i) -> Vector{Tuple{Tuple{Symbol, Int}, SimFloat}}

Whatever integrates the translation of point `i` of `sam`, as `(collection, index)`, with
each one's share — the only things the air it entrains can slow down:

- a free `DYNAMIC` node integrates itself and takes all of it;
- a node anchored to a rigid body is placed by that body, which takes all of it;
- a node riding a Timoshenko tube's deformed centerline is placed by the two
  bodies the tube spans, which split it by where along the element it sits
  (`beam_frac`).
"""
function apparent_mass_carriers(sam, i)
    point = sam.definition.points[i]
    if !isnothing(point.tube)
        tube_bodies = sam.definition.tubes[point.tube].bodies
        frac = clamp(run_value(sam, Val(:points), i, Val(:beam_frac)), 0.0, 1.0)
        return [((:bodies, Int(tube_bodies[1])), 1.0 - frac),
                ((:bodies, Int(tube_bodies[2])), frac)]
    end
    !isnothing(point.body) && return [((:bodies, Int(point.body)), 1.0)]
    return [((:points, i), 1.0)]
end

"""
    wagner_reference_frame(sam, wing_idx) -> (x_ref, z_ref, chord_ref)

The mean chordwise and normal directions in body frame of wing `wing_idx` of `sam` and its
mean chord [m], averaged over the frozen VSM mesh. [`wagner_lag_eqs`](@ref) reads its
one angle of attack against these, so the lag follows the whole wing rather than any
one panel. Each panel contributes the axes `VortexStepMethod.panel_axes` builds at its
[`chord_blend_weights!`](@ref) blend weight, which are the axes the panel itself
carries.
"""
function wagner_reference_frame(sam, wing_idx)
    mode = sam.aero[wing_idx]
    refined = mode.vsm_wing.refined_sections
    n = Int(mode.vsm_wing.n_panels)
    chord, width = panel_chord_width(mode)
    weight = chord_blend_weights!(similar(width), width)
    x_ref = zeros(SimFloat, 3)
    z_ref = zeros(SimFloat, 3)
    for i in 1:n
        left, right = refined[i], refined[i + 1]
        axes = VortexStepMethod.panel_axes(left.LE_point, left.TE_point,
                                           right.LE_point, right.TE_point,
                                           weight[i])
        x_ref .+= axes.x_airf
        z_ref .+= axes.z_airf
    end
    return x_ref ./ norm(x_ref), z_ref ./ norm(z_ref), sum(chord) / n
end

"""
    wagner_wing_eqs(sam, wing_idx, sec_va, params) -> (eqs, vars, deficiency, initial)

The Wagner lag of wing `wing_idx` of `sam` from its live section inflow, or an empty set
when the wing carries no lag. Wraps [`wagner_lag_eqs`](@ref) with the mean of `sec_va` as
the wing's apparent wind and the frozen [`wagner_reference_frame`](@ref) as the angle
it is measured against, so both particle aero modes build the lag the same way.
"""
function wagner_wing_eqs(sam, wing_idx, sec_va, params)
    wagner_enabled(sam.aero[wing_idx]) || return (Equation[], Any[], 0.0, Dict())
    x_ref, z_ref, chord_ref = wagner_reference_frame(sam, wing_idx)
    va_ref = sum(sec_va) ./ length(sec_va)
    return wagner_lag_eqs(wagner_gain_params(params, wing_idx),
                          wagner_rate_params(params, wing_idx),
                          va_ref, x_ref, z_ref, chord_ref)
end

"""
    wagner_gain_params(params, wing_idx)

The registered Wagner gains `(A₁, A₂)` of wing `wing_idx`. An ordinary parameter like every
other tunable scalar, which is what keeps a change of lag out of the model cache key.
"""
wagner_gain_params(params, wing_idx) = params.wings[wing_idx].aero.unsteady.wagner_gains

"""
    wagner_rate_params(params, wing_idx)

The registered Wagner rates `(b₁, b₂)` of wing `wing_idx`, tunable in the same way as
[`wagner_gain_params`](@ref).
"""
wagner_rate_params(params, wing_idx) = params.wings[wing_idx].aero.unsteady.wagner_rates

"""
    wagner_params(sam, wing_idx, params) -> Vector{Any}

The lag's parameters of wing `wing_idx` of `sam` for a component's declared parameter list,
empty for a wing without one. A whole-wing aero component names its parameters
explicitly, so reading these in [`wagner_wing_eqs`](@ref) is not enough to declare them.
"""
wagner_params(sam, wing_idx, params) =
    wagner_enabled(sam.aero[wing_idx]) ?
        Any[wagner_gain_params(params, wing_idx),
            wagner_rate_params(params, wing_idx)] : Any[]

"""
    wagner_lag_eqs(gains, rates, va_ref, x_ref, z_ref, chord_ref)
        -> (eqs, vars, deficiency, initial)

The wing's two-state Wagner lift lag. A step in angle of attack does not build its
full circulatory lift at once; Wagner's indicial function
`φ(s) = 1 - A₁·exp(-b₁·s) - A₂·exp(-b₂·s)` gives the fraction reached after `s`
semi-chords of travel, and starts at `φ(0) = 1 - A₁ - A₂ = 0.5`.

Driving two states with `dxᵢ/ds = α - bᵢ·xᵢ` and reading the deficiency
`d = Σ Aᵢ·(α - bᵢ·xᵢ)` reproduces `α_E = α - d = α·φ(s)` after a step while leaving
`d = 0` in steady flow, and needs no `dα/dt` on the right-hand side. Converting to
time with `s = ∫2·v/c·dt` gives the emitted `D(xᵢ)`.

Steady flow puts `xᵢ` at `α/bᵢ`, which is what the returned defaults set, so a model
starts trimmed instead of building its lift up from `φ(0) = 0.5`.

`gains` and `rates` are the wing's registered `(A₁, A₂)` and `(b₁, b₂)`
parameters, so retuning the lag syncs rather than rebuilds. `va_ref` is the wing's
mean body-frame apparent wind, measured against `x_ref`/`z_ref`/`chord_ref` from
[`wagner_reference_frame`](@ref). It carries no induced velocity, unlike a panel's
`v_eff`, so this angle sits a downwash off the panels'; only changes in it drive
the lag, which leaves a steady offset at zero deficiency. The one `deficiency`
shifts every panel's angle of attack in [`panel_force_eqs`](@ref), so the whole
wing lags together and the spanwise shape of the loading is untouched.
"""
function wagner_lag_eqs(gains, rates, va_ref, x_ref, z_ref, chord_ref)
    vars = @variables(wagner_state(t)[1:2], wagner_alpha(t), wagner_speed(t),
                      wagner_deficiency(t))
    state, alpha_ref, speed, deficiency = vars
    alpha_rhs = atan(dot(va_ref, z_ref), dot(va_ref, x_ref))
    eqs = [
        alpha_ref ~ alpha_rhs;
        speed ~ smooth_norm(collect(va_ref));
        [D(state[i]) ~ (2 * speed / chord_ref) *
            (alpha_ref - rates[i] * state[i]) for i in 1:2];
        deficiency ~ sum(gains[i] * (alpha_ref - rates[i] * state[i])
                         for i in 1:2)]
    defaults = Dict(state => [alpha_rhs / rates[i] for i in 1:2])
    return eqs, Any[state, alpha_ref, speed, deficiency], deficiency, defaults
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

`mode`'s panels addressed one at a time, reached as `wings[wing_idx].aero.panels[i]`. It
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

`mode`'s wing points addressed one at a time, reached as `wings[wing_idx].aero.points[i]`,
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
    sym === :traction_net && return mode.traction_net[:, i]
    sym === :traction_moment && return mode.traction_moment[:, i]
    sym === :residual_arm && return mode.residual_arm[:, i]
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
    aero_inflow_groups(mode, sam, wing_idx) -> (groups, section_group)

How each refined section's apparent wind and density are averaged over the structural
points of wing `wing_idx` of `sam`. `groups[g]` is a list of `(point column, weight)` whose
weights sum to one, and `section_group[s]` names the group refined section `s` reads.
Sections sharing a group share the one [`AeroInflow`](@ref) that averages it.

The default is one group per refined section, gathered from the same bounding struts
[`aero_geometry_entries`](@ref) reconstructs that section's corners from, so a
section's inflow and its geometry read the same points. Overriding this with a
wing-wide mean drops the rotational part of the velocity field — a rigid rotation
about the point centroid averages to zero — leaving the panels without rate damping.
"""
aero_inflow_groups(mode, sam, wing_idx) =
    aero_station_groups(mode, sam, wing_idx, strut_inflow_weights)

"""
    aero_pitch_groups(mode, sam, wing_idx) -> (groups, section_group)

[`aero_inflow_groups`](@ref) with the trailing minus leading edge weights instead of
the chordwise blend, so the same gather yields each section's pitch-rate input
([`strut_pitch_weights`](@ref)). Only built for a wing with
[`flow_curvature_enabled`](@ref).
"""
aero_pitch_groups(mode, sam, wing_idx) =
    aero_station_groups(mode, sam, wing_idx, strut_pitch_weights)

"""
    aero_station_groups(mode, sam, wing_idx, station_weights) -> (groups, section_group)

One group per refined section of wing `wing_idx` of `sam`, each gathered by
`station_weights`. Shared by [`aero_inflow_groups`](@ref) and [`aero_pitch_groups`](@ref),
which differ only in the chordwise weights they ask for.
"""
function aero_station_groups(mode, sam, wing_idx, station_weights)
    column = aero_section_columns(sam, wing_idx)
    sections = eachindex(section_interp_caches(mode)[1])
    groups = [station_weights(mode, section, column) for section in sections]
    return groups, collect(sections)
end

"""
    aero_scatter_entries(mode, sam, wing_idx) -> Vector

The linear map from the panel loads of wing `wing_idx` of `sam` to its point forces as
`(panel, point column, force weight, couple weight)`: point `k`'s body-frame force
gains `force weight · panel_force + couple weight · panel_couple`. The weights are
constants of the mesh, so the wiring layer carries the whole scatter and no component
holds it.
"""
function aero_scatter_entries end

"""
    aero_point_offset(mode, params, wing_idx, i)

The constant body-frame force point `i` of wing `wing_idx` receives on top of the scattered
panel loads, or `nothing`. [`AeroPressure`](@ref) puts its frozen surface traction here,
net of the share of each panel's frozen total that the scatter already re-adds.
"""
aero_point_offset(::AbstractAeroModel, params, wing_idx, i) = nothing

"""
    aero_geometry_entries(mode, sam, wing_idx) -> Vector

Where each panel's four section corners of wing `wing_idx` of `sam` come from, as `(panel,
corner, point column, weight)` with `corner` one of `:le_a`, `:te_a`, `:le_b`, `:te_b`.
This is [`interp_sections`](@ref) split in two: the strut weights the wiring gathers,
and the frozen billow offset the panel adds as a parameter. Mode-independent — every
continuous mode reconstructs its sections the same way.
"""
function aero_geometry_entries(mode, sam, wing_idx)
    column = aero_section_columns(sam, wing_idx)
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
    strut_station_weights(mode, section, column, le_weight, te_weight)
        -> Vector{Tuple{Int, SimFloat}}

The `(point column, weight)` pairs that combine refined section `section`'s two
bounding struts, each strut weighting its LE station by `le_weight` and its TE
station by `te_weight`. The strut shares are the section's own interpolation
weights, so whatever chordwise combination the caller asks for is taken at the
same struts [`aero_geometry_entries`](@ref) builds that section's corners from.
"""
function strut_station_weights(mode, section::Int, column, le_weight, te_weight)
    left, weight, _, _ = section_interp_caches(mode)
    pairs = Tuple{Int, SimFloat}[]
    for (strut, share) in ((left[section], weight[section]),
                           (left[section] + 1, 1.0 - weight[section]))
        share == 0.0 && continue
        push!(pairs, (column[(strut, :LE)], le_weight * share))
        push!(pairs, (column[(strut, :TE)], te_weight * share))
    end
    return pairs
end

"""
    strut_inflow_weights(mode, section, column) -> Vector{Tuple{Int, SimFloat}}

The `(point column, weight)` pairs averaging refined section `section`'s inflow over
its two bounding struts, each strut blending its LE and TE station at
[`COLLOCATION_CHORD_FRAC`](@ref). With no chord twist rate the two stations share a
velocity, so the blend is inert in steady flight whatever the fraction.
"""
strut_inflow_weights(mode, section::Int, column) =
    strut_station_weights(mode, section, column,
                          1 - COLLOCATION_CHORD_FRAC, COLLOCATION_CHORD_FRAC)

"""
    strut_pitch_weights(mode, section, column) -> Vector{Tuple{Int, SimFloat}}

The `(point column, weight)` pairs taking refined section `section`'s trailing minus
leading edge apparent wind, which `VortexStepMethod.section_pitch_rate` turns into the
section's rotation rate about its own spanwise axis. Weights sum to zero, so a wing
in uniform translation contributes nothing.
"""
strut_pitch_weights(mode, section::Int, column) =
    strut_station_weights(mode, section, column, -1.0, 1.0)

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
    panel_corners(panel) -> (le_1, te_1, le_2, te_2)

A panel's four `corner_points` in the argument order
`VortexStepMethod.panel_axes` takes them; they are stored leading edge 1,
trailing edge 1, trailing edge 2, leading edge 2.
"""
panel_corners(panel) = (SVector{3}(@view panel.corner_points[:, 1]),
                        SVector{3}(@view panel.corner_points[:, 2]),
                        SVector{3}(@view panel.corner_points[:, 4]),
                        SVector{3}(@view panel.corner_points[:, 3]))

"""
    panel_span_width(panel) -> SimFloat

The panel's span width: the length of the quarter-chord step between its two
sections, from its `corner_points`.
"""
panel_span_width(panel) = smooth_norm(
    VortexStepMethod.panel_span_vector(panel_corners(panel)...))

"""
    store_chord_weights!(chord_weight, body_aero)

Freeze the live panels' [`chord_blend_weights!`](@ref) into `chord_weight`
(n_panels), so the weight follows the mesh as the wing deforms. Written at the
same refresh as [`store_induced_velocity!`](@ref).
"""
function store_chord_weights!(chord_weight, body_aero)
    panels = body_aero.panels
    length(chord_weight) == length(panels) || error(
        "chord-weight buffer is stale ($(length(chord_weight)) for " *
        "$(length(panels)) panels); reinitialize the model.")
    chord_blend_weights!(chord_weight, [panel.width for panel in panels])
    return nothing
end

"""
    corner_chord_weights(mode) -> Vector{SimFloat}

Every panel's [`chord_blend_weights!`](@ref) of aero `mode` over the span widths of its
`corner_points` ([`panel_span_width`](@ref)).
"""
function corner_chord_weights(mode)
    width = map(panel_span_width, mode.vsm_aero.panels)
    return chord_blend_weights!(similar(width), width)
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
    aero_component(mode, dynamics::WingDynamics, sam, wing_idx; name, params) -> System

Build the aero subsystem of wing `wing_idx` of `sam`, selected by dispatch on both its aero
`mode` and its [`WingDynamics`](@ref). Returns a `System` exposing the connectors fixed
by the dynamics type, all in the wing body frame (the wiring layer `aero_eqs!` drives
inputs and reads outputs; connectors a mode ignores still exist for binding):
- `RIGID_DYNAMICS` (`num` stations, [`wing_station_idxs`](@ref)): in `va[1:3]`,
  `rho`, `R_b_w[1:3,1:3]`, `omega[1:3]`, `twist[1:num]`, `twist_vel[1:num]`;
  out `force[1:3]`, `moment[1:3]`, `twist_moment[1:num]`.
- `PARTICLE_DYNAMICS` (`np` wing nodes, [`wing_points`](@ref)): in
  `point_pos[1:3,1:np]`, `point_vel[1:3,1:np]`, `va[1:3,1:np]`, `rho[1:np]`; out
  `point_force[1:3,1:np]`.

A mode supports a dynamics type by defining the matching method (rigid,
particle, or both). Add a method on a custom `AbstractAeroModel` subtype to
plug in your own aerodynamics.
"""
function aero_component end

"""
    validate_aero_component(subsys, sam, wing_idx)

Check the built aero `subsys` of wing `wing_idx` of `sam` exposes the connectors the wiring
layer needs for its dynamics type; error naming the missing connector otherwise.
"""
function validate_aero_component(subsys, sam, wing_idx)
    if body_dynamics_type(sam, wing_body(sam, wing_idx)) == RIGID_DYNAMICS
        required = Symbol[:va, :rho, :R_b_w, :omega, :force, :moment]
        isempty(wing_station_idxs(sam, wing_idx)) ||
            append!(required, [:twist, :twist_vel, :twist_moment])
    else
        required = Symbol[:point_pos, :point_vel, :va, :rho, :point_force]
    end
    required_str = join(required, ", ")
    for con in required
        hasproperty(subsys, con) || error(
            "Wing $(sam.definition.wings[wing_idx].name): aero component is missing " *
            "required connector `$con`. Required: $required_str.")
    end
    return nothing
end

# ==================== refresh orchestrator ==================== #

"""
    sync_aero_density!(mode, sam, wing_idx)

Set the VSM solver air density of wing `wing_idx` of `sam` to the density at the height of
its body, the same altitude-dependent density the symbolic RHS uses to dimensionalize aero
forces (see `aero_eqs.jl`). Keeps the VSM solve and the model consistent on dynamic
pressure. No-op for non-VSM aero modes.
"""
sync_aero_density!(::AbstractAeroModel, sam, wing_idx) = nothing
function sync_aero_density!(mode::AbstractVSMAero, sam, wing_idx)
    height = run_value(sam, Val(:bodies), wing_body(sam, wing_idx), Val(:pos_w))[3]
    mode.vsm_solver.density = air_density(sam.am, height)
    return nothing
end

"""
    point_acceleration_w(sam, i, R_b_to_w, frame) -> SVector

The world-frame acceleration of free point `i` of `sam`, rebuilt from its run: its net
force per unit mass less the world-frame damping and the body-frame damping on its wing's
axes `R_b_to_w` against the wing's [`small_earth_frame`](@ref) `frame`. The same
expression `point_eqs!` binds to `acc`, which is where the monolith's fitted wing reads
its own `acc_w` from.
"""
function point_acceleration_w(sam, i, R_b_to_w, frame)
    pos_w = run_value(sam, Val(:points), i, Val(:pos_w))
    vel_w = run_value(sam, Val(:points), i, Val(:vel_w))
    damping = run_value(sam, Val(:points), i, Val(:world_frame_damping)) .* vel_w .+
        body_frame_damp_accel(vel_w,
                              run_value(sam, Val(:points), i, Val(:body_frame_damping)),
                              R_b_to_w, frame, pos_w)
    return run_value(sam, Val(:points), i, Val(:force)) ./
        run_value(sam, Val(:points), i, Val(:total_mass)) .- damping
end

"""
    wing_kinematics_from_points!(sam, wing_idx, frame, aero_points; base_point=0,
                                 frame_base=nothing)

Recompute the kinematic state of `PARTICLE_DYNAMICS` wing `wing_idx` of `sam` from the
current point positions and velocities, for a backend that does not carry these quantities
as state. `frame` is the wing's KiteGeometry `WingFrame`, each of whose references is the
weighted blend of its points. Writes its body's orientation
([`wing_frame_columns`](@ref)), origin `pos_w`/`vel_w`, the frame's own `ω_b`
([`fitted_frame_omega`](@ref)), the origin's acceleration `acc_w`
([`point_acceleration_w`](@ref), damped against the [`small_earth_frame`](@ref) about
the [`transform_base`](@ref) `frame_base`), the reported scalars
([`write_wing_scalars!`](@ref)) and, from the height profile, the wind of the body and
of each of `aero_points`; under [`PerPointWind`](@ref) those winds are the caller's. The
same quantities the monolith's `get_all_state` copies out of the integrator.
"""
function wing_kinematics_from_points!(sam, wing_idx, frame, aero_points; base_point=0,
                                      frame_base=nothing)
    body_idx = wing_body(sam, wing_idx)
    positions, velocities = run_positions(sam), run_velocities(sam)
    (zp1, zp2), (yp1, yp2), origin = frame.z_ref_points, frame.y_ref_points, frame.origin
    pos_z1 = KiteGeometry.ref_position(positions, zp1)
    pos_z2 = KiteGeometry.ref_position(positions, zp2)
    pos_y1 = KiteGeometry.ref_position(positions, yp1)
    pos_y2 = KiteGeometry.ref_position(positions, yp2)
    rates = (KiteGeometry.ref_position(velocities, zp1),
             KiteGeometry.ref_position(velocities, zp2),
             KiteGeometry.ref_position(velocities, yp1),
             KiteGeometry.ref_position(velocities, yp2))
    x, y, z = wing_frame_columns(pos_z1, pos_z2, pos_y1, pos_y2)
    R_b_to_w = hcat(x, y, z)
    pos_w = SVec3(KiteGeometry.ref_position(positions, origin))
    vel_w = SVec3(KiteGeometry.ref_position(velocities, origin))
    axes = (collect(x), collect(y), collect(z))
    ω_b = fitted_frame_omega(pos_z1, pos_z2, pos_y1, pos_y2, rates, axes)
    base_w = damping_base_position(frame_base, point -> positions[point])
    flight = small_earth_frame(base_w, pos_w, vel_w, R_b_to_w[:, 2], R_b_to_w * ω_b)
    acc_w = zero(SVec3)
    for (i, weight) in origin
        acc_w += weight .* point_acceleration_w(sam, i, R_b_to_w, flight)
    end
    set_run_value!(sam, Val(:bodies), body_idx, Val(:Q_b_to_w),
                   rotation_matrix_to_quaternion(R_b_to_w))
    set_run_value!(sam, Val(:bodies), body_idx, Val(:pos_w), pos_w)
    set_run_value!(sam, Val(:bodies), body_idx, Val(:vel_w), vel_w)
    set_run_value!(sam, Val(:bodies), body_idx, Val(:ω_b), ω_b)
    set_run_value!(sam, Val(:bodies), body_idx, Val(:acc_w), acc_w)
    if !per_point_wind(sam)
        profile = WindFactor(sam.am, sam.set.profile_law)
        wind_vec = sam.set.wind_vec
        set_run_value!(sam, Val(:bodies), body_idx, Val(:wind_vec),
                       profile(pos_w[3]) .* wind_vec)
        for i in aero_points
            height = run_value(sam, Val(:points), i, Val(:pos_w))[3]
            set_run_value!(sam, Val(:points), i, Val(:wind_vec),
                           profile(height) .* wind_vec)
        end
    end
    write_wing_scalars!(sam, wing_idx; base_point)
    return nothing
end

"""
    write_wing_scalars!(sam, wing_idx; base_point=0, alpha_b=zeros(3)) -> nothing

Fill the reported scalars of wing `wing_idx` of `sam` — heading, course, elevation, azimuth,
angle of attack and `turn_acc` — from the pose and motion of its body, via the one
definition in [`wing_scalar_kinematics`](@ref). `base_point` is the point its elevation
and azimuth are measured from, 0 for the origin; `turn_acc` follows the body-frame
angular acceleration `alpha_b`, zero for a fitted wing as the monolith binds it.
"""
function write_wing_scalars!(sam, wing_idx; base_point=0, alpha_b=zeros(3))
    body_idx = wing_body(sam, wing_idx)
    pos_w = run_value(sam, Val(:bodies), body_idx, Val(:pos_w))
    rel_pos = collect(base_point == 0 ? pos_w :
        pos_w - run_value(sam, Val(:points), base_point, Val(:pos_w)))
    twist_offset = mid_span_twist(
        station_idx -> run_value(sam, Val(:spanwise_stations), station_idx, Val(:twist)),
        wing_station_idxs(sam, wing_idx))
    R_b_to_w = collect(body_rotation(sam, body_idx))
    e_x = R_b_to_w[:, 1]
    scalars = wing_scalar_kinematics(; rel_pos, e_x,
        R_t_to_w = sym_calc_R_t_to_w(rel_pos),
        R_v_to_w = calc_R_v_to_w(rel_pos, e_x), R_b_to_w,
        vel = collect(run_value(sam, Val(:bodies), body_idx, Val(:vel_w))),
        acc = collect(run_value(sam, Val(:bodies), body_idx, Val(:acc_w))),
        omega_b = collect(run_value(sam, Val(:bodies), body_idx, Val(:ω_b))),
        alpha_b = collect(alpha_b),
        va_b = collect(run_value(sam, Val(:bodies), body_idx, Val(:va_b))), twist_offset)
    set_run_value!(sam, Val(:wings), wing_idx, Val(:turn_acc), scalars.turn_acc)
    set_run_value!(sam, Val(:wings), wing_idx, Val(:heading), scalars.heading)
    set_run_value!(sam, Val(:wings), wing_idx, Val(:elevation), scalars.elevation)
    set_run_value!(sam, Val(:wings), wing_idx, Val(:azimuth), scalars.azimuth)
    set_run_value!(sam, Val(:wings), wing_idx, Val(:course), scalars.course)
    set_run_value!(sam, Val(:wings), wing_idx, Val(:aoa), scalars.angle_of_attack)
    return nothing
end

"""
    refresh_aero!(sam::SymbolicAWEModel; vsm_min_wind=0.5, cold_start=false,
                  vsm_warn_on_fail=false)

Refresh each wing's aerodynamic state, dispatching on the wing's aero mode
([`refresh_rigid_aero!`](@ref) / [`refresh_particle_aero!`](@ref)), and re-spread
the entrained air ([`apply_apparent_mass!`](@ref)) at the density it just read. Runs
on the low-frequency VSM-update schedule (`vsm_interval`), not the compiled RHS. Reads
the apparent wind of each point from the run (written by `update_run_state!`, which
always runs first), so it does not re-extract integrator state.

**RIGID_DYNAMICS VSM modes:** compute wind-axis coefficients (CL, CD, CS, CM, cm)
at the operating point, plus the `ForwardDiff` Jacobian over `[α, β, ω₁, ω₂, ω₃,
θ_twist…]` (`AeroLinearized`) or the frozen forces (`AeroDirect`).

**PARTICLE_DYNAMICS VSM modes:** full nonlinear VSM solve with per-point force
distribution. Non-VSM modes (`AeroNone`/`AeroPlate`) are no-ops, so callers
should gate this on [`has_vsm_wing`](@ref).

Each VSM solve starts from the wing's [`circulation`](@ref) in the run's state and
leaves its converged circulation there. `cold_start` clears it first, so the solves start
from the solver's initial distribution: past stall the iteration has more than one fixed
point, and the first solve after a [`reinit!`](@ref) must not depend on what ran before.

`vsm_warn_on_fail` downgrades a `VortexStepMethod.SolveFailure` to a warning
([`warn_or_rethrow`](@ref)).
"""
function refresh_aero!(sam::SymbolicAWEModel; vsm_min_wind=0.5,
                       cold_start=false, vsm_warn_on_fail=false)
    isempty(sam.aero) && return nothing
    if cold_start
        foreach(wing_idx -> fill!(circulation(sam, wing_idx), 0.0), eachindex(sam.aero))
    end
    for (wing_idx, mode) in enumerate(sam.aero)
        sync_aero_density!(mode, sam, wing_idx)
    end
    apply_apparent_mass!(sam)

    particle_wings = Int[]
    for (wing_idx, mode) in enumerate(sam.aero)
        if body_dynamics_type(sam, wing_body(sam, wing_idx)) == PARTICLE_DYNAMICS
            push!(particle_wings, wing_idx)
            continue
        end
        try
            refresh_rigid_aero!(mode, sam, wing_idx; vsm_min_wind)
        catch failure
            warn_or_rethrow(failure, vsm_warn_on_fail, sam, wing_idx)
        end
    end

    isempty(particle_wings) && return nothing
    va_point_b = point_apparent_winds(sam)
    for wing_idx in particle_wings
        mode = sam.aero[wing_idx]
        apply_flap_delta!(mode, sam, wing_idx)
        try
            refresh_particle_aero!(mode, sam, wing_idx, va_point_b; vsm_min_wind)
        catch failure
            warn_or_rethrow(failure, vsm_warn_on_fail, sam, wing_idx)
        end
    end
    return nothing
end

"""The apparent wind of every point of `sam` in the body frame of its wing, one column
per point [m/s], zero for a point of no wing."""
function point_apparent_winds(sam)
    va_point_b = Matrix{SimFloat}(undef, 3, length(sam.definition.points))
    for i in axes(va_point_b, 2)
        va_point_b[:, i] .= run_value(sam, Val(:points), i, Val(:va_b))
    end
    return va_point_b
end

"""
    warn_or_rethrow(failure, vsm_warn_on_fail, sam, wing_idx)

Rethrow `failure`, unless it is a `VortexStepMethod.SolveFailure` and
`vsm_warn_on_fail` is set: then warn instead, and wing `wing_idx` of `sam` flies on with the
circulation, the angles of attack and the frozen forces of its last converged solve.
"""
function warn_or_rethrow(failure, vsm_warn_on_fail, sam, wing_idx)
    (vsm_warn_on_fail && failure isa VortexStepMethod.SolveFailure) ||
        rethrow(failure)
    @warn "$(nameof(typeof(sam.aero[wing_idx]))) on wing $wing_idx: $(failure.msg) " *
          "Reusing the last converged aero state."
    return nothing
end

"""
    refresh_rigid_aero!(mode, sam, wing_idx; vsm_min_wind=0.5)

Refresh the aero state of `RIGID_DYNAMICS` wing `wing_idx` of `sam`, dispatched on its aero
`mode`:
- `AeroNone` / any non-VSM mode → no-op (fallback).
- `AeroLinearized` → compute the baseline coefficients ([`rigid_aero_baseline!`](@ref))
  and the `ForwardDiff` Jacobian `d(coeffs)/d(inputs)` into the mode's `aero_jac`.
- `AeroDirect` → compute the baseline coefficients and apply the frozen body-frame
  force/moment; below `vsm_min_wind` everything is zeroed.
"""
refresh_rigid_aero!(::AbstractAeroModel, sam, wing_idx; vsm_min_wind=0.5) = nothing

"""
    refresh_particle_aero!(mode, sam, wing_idx, va_point_b; vsm_min_wind=0.5)

Refresh the aero state of `PARTICLE_DYNAMICS` wing `wing_idx` of `sam`, dispatched on its
aero `mode`, `va_point_b` holding each point's body-frame apparent wind as a column:
- `AeroNone` / any non-VSM mode → no-op (fallback).
- `AeroDirect` → full nonlinear VSM solve with per-section apparent wind, then
  distribute panel forces onto the wing's structural points
  ([`distribute_panel_forces_to_points!`](@ref)); below `vsm_min_wind` the point
  forces are zeroed.
- `AeroLinearized` → unsupported (errors).
"""
refresh_particle_aero!(::AbstractAeroModel, sam, wing_idx, va_point_b; vsm_min_wind=0.5) =
    nothing

# ==================== per-wing lifecycle ==================== #

"""
    remake_aero!(mode, sam, wing_idx)

Rebuild the VSM objects of the aero `mode` of wing `wing_idx` of `sam` from its `set` and
`vsm_set` (the `remake_vsm` path of `init!`, after editing settings) and set them up
again ([`setup_model_aero!`](@ref)). A no-op for a mode without an engine.
"""
remake_aero!(::AbstractAeroModel, sam, wing_idx) = nothing

function remake_aero!(mode::AbstractVSMAero, sam, wing_idx)
    mode.vsm_wing = create_vsm_wing(sam.set, sam.vsm_set; prn=false, sort_sections=false)
    mode.vsm_aero = VortexStepMethod.BodyAerodynamics([mode.vsm_wing])
    mode.vsm_solver = build_vsm_solver(mode.vsm_aero, sam.vsm_set,
                                       body_dynamics_type(sam, wing_body(sam, wing_idx)),
                                       length(wing_station_idxs(sam, wing_idx)))
    return setup_model_aero!(sam, wing_idx)
end

"""
    setup_model_aero!(sam, wing_idx)

Set the aero mode of wing `wing_idx` of `sam` up on its structure, once per engine: size its
aero state ([`resize_aero_state!`](@ref)), put its panels in the body frame and its
sections and stations on the structure's ([`setup_aero!`](@ref)), and map its flap
stations ([`build_panel_station_map!`](@ref)). A no-op for a mode without VSM geometry.
"""
function setup_model_aero!(sam, wing_idx)
    mode = sam.aero[wing_idx]
    resize_aero_state!(mode, sam, wing_idx)
    setup_aero!(mode, sam, wing_idx)
    build_panel_station_map!(mode, sam, wing_idx)
    return nothing
end

"""
    validate_aero_structure(mode, sam, wing_idx; prn=false)

Check structural invariants the compiled equations of wing `wing_idx` of `sam` rely on under
its aero `mode` (run at build). Default no-op; VSM `PARTICLE_DYNAMICS` wings verify each
unrefined section has its LE and TE structural point mapped (interior chord control
points are allowed and need no mapping).
"""
validate_aero_structure(::AbstractAeroModel, sam, wing_idx; prn=false) = nothing

function validate_aero_structure(mode::AbstractVSMAero, sam, wing_idx; prn=false)
    body_dynamics_type(sam, wing_body(sam, wing_idx)) == PARTICLE_DYNAMICS || return nothing
    @assert !isnothing(mode.point_to_vsm_point) "PARTICLE_DYNAMICS wing $wing_idx " *
        "missing point_to_vsm_point mapping"
    n_sections = length(mode.vsm_wing.unrefined_sections)
    mapped = values(mode.point_to_vsm_point)
    for section_idx in 1:n_sections, le_or_te in (:LE, :TE)
        @assert (section_idx, le_or_te) in mapped "PARTICLE_DYNAMICS wing $wing_idx: " *
            "section $section_idx missing its $le_or_te point"
    end
    prn && println("✓ PARTICLE_DYNAMICS wing $wing_idx validated: " *
        "$(length(wing_points(sam, wing_idx))) points, $n_sections sections, " *
        "$(length(mode.vsm_aero.panels)) panels")
    return nothing
end

"""
    setup_aero!(mode, sam, wing_idx; prn=false)

Construction-time aero setup of wing `wing_idx` of `sam`, dispatched on its aero `mode`
(default no-op). VSM modes transform the VSM panels into the body frame, match aero
sections to the structure and map the stations onto them in the engine
([`compute_spatial_station_mapping!`](@ref), [`compute_station_geometry!`](@ref)), or
for a `PARTICLE_DYNAMICS` wing build the structural↔panel mapping. A custom mode adds a
method to participate in construction. A mode that reads how its sections were
generated finds the model's `VSMSettings` in `sam.vsm_set`.
"""
setup_aero!(::AbstractAeroModel, sam, wing_idx; prn=false) = nothing

"""
    transform_vsm_sections_to_body!(mode, sam, wing_idx; aero_z_offset=nothing)

Move the VSM sections/panels of the aero `mode` of wing `wing_idx` of `sam` from the frame
its aero geometry is authored in into the body frame, given the body origin and orientation
in that frame (`mesh_origin_b`, `R_b_to_mesh`), record both in `T_cad_body` /
`R_cad_body`, and reinit the panels. With `aero_z_offset` set, also apply the chordwise
aero z-offset (RIGID wings); PARTICLE wings pass `nothing`. Not idempotent: an engine is
transformed once.
"""
function transform_vsm_sections_to_body!(mode, sam, wing_idx; aero_z_offset=nothing)
    (; origin_b, R_b_to_mesh) = wing_mesh_pose(sam, wing_idx)
    vsm_wing = mode.vsm_wing
    vsm_wing.T_cad_body .= origin_b
    adjust_vsm_panels_to_origin!(vsm_wing, origin_b)
    rotate_vsm_sections!(vsm_wing, R_b_to_mesh')
    vsm_wing.R_cad_body .= R_b_to_mesh
    isnothing(aero_z_offset) || apply_aero_z_offset!(vsm_wing, aero_z_offset)
    check_span_order(mode, sam, wing_idx)
    VortexStepMethod.reinit!(mode.vsm_aero)
    return nothing
end

function setup_aero!(mode::AbstractVSMAero, sam, wing_idx; prn=false)
    if body_dynamics_type(sam, wing_body(sam, wing_idx)) == RIGID_DYNAMICS
        transform_vsm_sections_to_body!(mode, sam, wing_idx;
                                        aero_z_offset=mode.aero_z_offset)
        station_idxs = wing_station_idxs(sam, wing_idx)
        couples_to_sections(mode) && isempty(station_idxs) &&
            error("Section-coupled aero on RIGID wing $wing_idx requires explicit " *
                  "stations covering its LE/TE structural sections; none were " *
                  "declared. Add them to the wing.")
        couples_to_sections(mode) && match_aero_sections_to_structure!(mode, sam, wing_idx)
        isempty(station_idxs) || compute_spatial_station_mapping!(sam, wing_idx)
        compute_station_geometry!(sam, wing_idx)
        return nothing
    end
    isnothing(wing_origin(sam, wing_idx)) ||
        transform_vsm_sections_to_body!(mode, sam, wing_idx)
    couples_to_sections(mode) && match_aero_sections_to_structure!(mode, sam, wing_idx)
    setup_particle_point_mapping!(mode, sam, wing_idx)
    return nothing
end

"""
    resize_aero_state!(mode, sam, wing_idx)

Resize the per-wing aero state of the aero `mode` of wing `wing_idx` of `sam` to its
stations (the engine was sized from `n_unrefined`). Default no-op; VSM modes resize
`aero_y`/`aero_x`/`aero_jac` for `RIGID_DYNAMICS` wings.
"""
resize_aero_state!(::AbstractAeroModel, sam, wing_idx) = nothing

function resize_aero_state!(mode::AbstractVSMAero, sam, wing_idx)
    body_dynamics_type(sam, wing_body(sam, wing_idx)) == RIGID_DYNAMICS || return nothing
    n_stations = length(wing_station_idxs(sam, wing_idx))
    num_aero_outputs = 6 + n_stations
    num_aero_inputs = 5 + n_stations
    if length(mode.aero_x) != num_aero_outputs ||
            length(mode.aero_y) != num_aero_inputs
        mode.aero_y = zeros(SimFloat, num_aero_inputs)
        mode.aero_x = zeros(SimFloat, num_aero_outputs)
        mode.aero_jac = zeros(SimFloat, num_aero_outputs, num_aero_inputs)
    end
    return nothing
end

"""
    init_aero_state!(mode, sam, wing_idx, va_b)

Initialize the aero state of the aero `mode` of wing `wing_idx` of `sam` from its initial
body-frame apparent wind `va_b` (before the first refresh). Default no-op; VSM modes
write the operating-point angles α, β into `aero_y`.
"""
init_aero_state!(::AbstractAeroModel, sam, wing_idx, va_b) = nothing

function init_aero_state!(mode::AbstractVSMAero, sam, wing_idx, va_b)
    aero_y = mode.aero_y
    length(aero_y) < 2 && return nothing
    aero_y .= 0.0
    aero_y[1] = atan(va_b[3], va_b[1])
    aero_y[2] = atan(va_b[2], hypot(va_b[1], va_b[3]))
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
            share = weight[s]
            line_le = share .* unrefined[strut].LE_point .+
                (1.0 - share) .* unrefined[strut + 1].LE_point
            line_te = share .* unrefined[strut].TE_point .+
                (1.0 - share) .* unrefined[strut + 1].TE_point
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
    aero_section_columns(sam, wing_idx) -> Dict{Tuple{Int64,Symbol},Int}

Map `(unrefined_section, :LE/:TE)` of wing `wing_idx` of `sam` to the connector column
(position in [`wing_points`](@ref)) of the structural station point there, via its engine's
`point_to_vsm_point`. Interior chord control points (absent from the map) are skipped.
"""
function aero_section_columns(sam, wing_idx)
    point_to_vsm = sam.aero[wing_idx].point_to_vsm_point
    isnothing(point_to_vsm) && error("continuous aero: wing " *
        "$(sam.definition.wings[wing_idx].name) missing point_to_vsm_point mapping.")
    column = Dict{Tuple{Int64, Symbol}, Int}()
    for (k, i) in enumerate(wing_points(sam, wing_idx))
        entry = get(point_to_vsm, i, nothing)
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
    reconstruct_sections_sym(mode, connectors, column) -> (sec_le, sec_te)

Live symbolic body-frame LE/TE of every refined section of a wing flown with the
continuous VSM aero `mode`: the strut LE/TE connector positions (`connectors.point_pos`,
derived from world `pos`) interpolated by the frozen mesh weights plus the frozen billow
offset. Shared by all continuous VSM modes so the force model and the plot use identical
geometry.
"""
function reconstruct_sections_sym(mode, connectors, column)
    n_struts = Int(mode.vsm_wing.n_unrefined_sections)
    strut_le = [collect(connectors.point_pos[:, column[(s, :LE)]])
                for s in 1:n_struts]
    strut_te = [collect(connectors.point_pos[:, column[(s, :TE)]])
                for s in 1:n_struts]
    left, weight, le_offset, te_offset = section_interp_caches(mode)
    return interp_sections(strut_le, strut_te, left, weight, le_offset, te_offset)
end

"""
    reconstruct_inflow_sym(mode, connectors, column) -> (sec_va, sec_rho, sec_dva)

Live symbolic body-frame apparent wind and density of every refined section of a wing
flown with the continuous VSM aero `mode`: each strut's LE/TE connector values blended
at [`COLLOCATION_CHORD_FRAC`](@ref) into one strut value, then interpolated by the
frozen mesh weights — the same struts and weights [`reconstruct_sections_sym`](@ref)
builds that section's corners from, so a section's inflow and its geometry read the
same points. `sec_dva` is the same interpolation of each strut's trailing minus leading
edge apparent wind, the pitch-rate input of [`panel_force_eqs`](@ref); it is `nothing`
unless the solver asks for the flow curvature term. Symbolic twin of
[`strut_inflow_weights`](@ref) and [`strut_pitch_weights`](@ref).
"""
function reconstruct_inflow_sym(mode, connectors, column)
    n_struts = Int(mode.vsm_wing.n_unrefined_sections)
    nose = 1 - COLLOCATION_CHORD_FRAC
    va_at(s, edge) = collect(connectors.va[:, column[(s, edge)]])
    strut_va = [nose * va_at(s, :LE) + COLLOCATION_CHORD_FRAC * va_at(s, :TE)
                for s in 1:n_struts]
    strut_rho = [nose * connectors.rho[column[(s, :LE)]] +
                 COLLOCATION_CHORD_FRAC * connectors.rho[column[(s, :TE)]]
                 for s in 1:n_struts]
    left, weight, _, _ = section_interp_caches(mode)
    sec_va = [interp_strut(strut_va, left, weight, s) for s in eachindex(left)]
    sec_rho = [interp_strut(strut_rho, left, weight, s) for s in eachindex(left)]
    flow_curvature_enabled(mode) || return sec_va, sec_rho, nothing
    strut_dva = [va_at(s, :TE) - va_at(s, :LE) for s in 1:n_struts]
    sec_dva = [interp_strut(strut_dva, left, weight, s) for s in eachindex(left)]
    return sec_va, sec_rho, sec_dva
end

"""The position of point `i` of `sam` from the origin of the body `body_idx`, in its
frame of the run [m]."""
point_position_b(sam, body_idx, i) = body_rotation(sam, body_idx)' *
    (run_value(sam, Val(:points), i, Val(:pos_w)) -
     run_value(sam, Val(:bodies), body_idx, Val(:pos_w)))

"""
    reconstruct_sections_b(mode, sam, wing_idx) -> (sec_le, sec_te)

Numeric body-frame refined-section LE/TE of wing `wing_idx` of `sam` from the live
structural positions, exactly as [`reconstruct_sections_sym`](@ref) builds them
symbolically. For plotting and logging so the drawn panels are the ones the dynamics use.
"""
function reconstruct_sections_b(mode, sam, wing_idx)
    column = aero_section_columns(sam, wing_idx)
    point_idxs = wing_points(sam, wing_idx)
    body_idx = wing_body(sam, wing_idx)
    n_struts = Int(mode.vsm_wing.n_unrefined_sections)
    strut_le = [point_position_b(sam, body_idx, point_idxs[column[(s, :LE)]])
                for s in 1:n_struts]
    strut_te = [point_position_b(sam, body_idx, point_idxs[column[(s, :TE)]])
                for s in 1:n_struts]
    left, weight, le_offset, te_offset = section_interp_caches(mode)
    return interp_sections(strut_le, strut_te, left, weight, le_offset, te_offset)
end

"""
    write_log_corners!(sys_state, point_idx, corners_b, sam, wing_idx, zoom) -> Int

Write the body-frame `corners_b` of wing `wing_idx` of `sam`, turned into the world and
scaled by `zoom`, into `sys_state.X/Y/Z` after `point_idx`; return the last index written.
"""
function write_log_corners!(sys_state, point_idx, corners_b, sam, wing_idx, zoom)
    body_idx = wing_body(sam, wing_idx)
    R_b_to_w = body_rotation(sam, body_idx)
    pos_w = run_value(sam, Val(:bodies), body_idx, Val(:pos_w))
    for corner_b in corners_b
        point_idx += 1
        corner_w = pos_w + R_b_to_w * corner_b
        sys_state.X[point_idx] = corner_w[1] * zoom
        sys_state.Y[point_idx] = corner_w[2] * zoom
        sys_state.Z[point_idx] = corner_w[3] * zoom
    end
    return point_idx
end

"""
    write_live_aero_log_points!(mode, sam, wing_idx, sys_state, point_idx, zoom)

Log the panel corners of wing `wing_idx` of `sam` the force model reconstructs (strut
interpolation + frozen billow offset), not the raw VSM mesh, so the plot shows the
deforming geometry the dynamics use. Shared by the continuous VSM modes.
"""
function write_live_aero_log_points!(mode, sam, wing_idx, sys_state, point_idx, zoom)
    sec_le, sec_te = reconstruct_sections_b(mode, sam, wing_idx)
    corners = (corner for i in 1:(length(sec_le) - 1)
               for corner in (sec_le[i], sec_te[i], sec_te[i + 1], sec_le[i + 1]))
    return write_log_corners!(sys_state, point_idx, corners, sam, wing_idx, zoom)
end

# ==================== logging / visualization hooks ==================== #

"""
    n_aero_log_points(mode, sys, wing_idx) -> Int

Number of extra `SysState` log slots the aero `mode` of wing `wing_idx` of definition `sys`
contributes (visualization geometry such as panel corners). Default 0; VSM modes log
4 corners per panel. Must match what [`write_aero_log_points!`](@ref) writes.
"""
n_aero_log_points(::AbstractAeroModel, sys, wing_idx) = 0
n_aero_log_points(mode::AbstractVSMAero, sys, wing_idx) = 4 * length(mode.vsm_aero.panels)

"""
    write_aero_log_points!(mode, sam, wing_idx, sys_state, point_idx, zoom) -> Int

Write the log points of the aero `mode` of wing `wing_idx` of `sam` (world frame, scaled by
`zoom`) into `sys_state.X/Y/Z` starting after `point_idx`; return the last index
written. Default writes nothing; VSM modes write the panel corners,
[`AeroPlate`](@ref) writes each section's display quad. Flap deflections are logged
separately, per station, into `sys_state.flap_angle` (see
[`write_flap_deflections!`](@ref)).
"""
write_aero_log_points!(::AbstractAeroModel, sam, wing_idx, sys_state, point_idx, zoom) =
    point_idx

function write_aero_log_points!(mode::AbstractVSMAero, sam, wing_idx, sys_state, point_idx,
                                zoom)
    corners = (SVec3(panel.corner_points[:, j])
               for panel in mode.vsm_aero.panels for j in 1:4)
    return write_log_corners!(sys_state, point_idx, corners, sam, wing_idx, zoom)
end

"""
    read_aero_log_points!(mode, sam, wing_idx, sys_state, point_idx) -> Int

Inverse of [`write_aero_log_points!`](@ref): restore the state of the aero `mode` of
wing `wing_idx` of `sam` from the logged points starting after `point_idx`; return the last
index consumed (the slots must be skipped even when unused). Default consumes nothing;
VSM `PARTICLE_DYNAMICS` modes read the panel corners back (rigid wings recompute panels
from twist instead and only skip their slots). Flap deflections are restored separately
from `sys_state.flap_angle` (see [`restore_flap_delta!`](@ref)).
"""
read_aero_log_points!(::AbstractAeroModel, sam, wing_idx, sys_state, point_idx) = point_idx

function read_aero_log_points!(mode::AbstractVSMAero, sam, wing_idx, sys_state, point_idx)
    n_corners = 4 * length(mode.vsm_aero.panels)
    body_idx = wing_body(sam, wing_idx)
    body_dynamics_type(sam, body_idx) == RIGID_DYNAMICS && return point_idx + n_corners
    R_w_to_b = body_rotation(sam, body_idx)'
    pos_w = run_value(sam, Val(:bodies), body_idx, Val(:pos_w))
    for panel in mode.vsm_aero.panels, j in 1:4
        point_idx += 1
        corner_w = SVec3(sys_state.X[point_idx], sys_state.Y[point_idx],
                         sys_state.Z[point_idx])
        panel.corner_points[:, j] .= R_w_to_b * (corner_w - pos_w)
    end
    return point_idx
end

"""
    restore_aero_twist!(mode, sam, wing_idx)

Re-apply the (already restored) station twists of wing `wing_idx` of `sam` to the geometry
of its aero `mode` when loading a `SysState` log frame. Default no-op; VSM
`RIGID_DYNAMICS` modes deform the unrefined sections and reinit the panels.
"""
restore_aero_twist!(::AbstractAeroModel, sam, wing_idx) = nothing

function restore_aero_twist!(mode::AbstractVSMAero, sam, wing_idx)
    body_dynamics_type(sam, wing_body(sam, wing_idx)) == RIGID_DYNAMICS || return nothing
    station_idxs = wing_station_idxs(sam, wing_idx)
    isempty(station_idxs) && return nothing
    vsm = mode.vsm_wing
    isempty(vsm.non_deformed_sections) && return nothing
    theta = zeros(Float64, vsm.n_unrefined_sections)
    for (k, station_idx) in enumerate(station_idxs)
        twist = run_value(sam, Val(:spanwise_stations), station_idx, Val(:twist))
        for section_idx in mode.engine.station_sections[k]
            theta[section_idx] = twist
        end
    end
    VortexStepMethod.unrefined_deform!(vsm, theta)
    VortexStepMethod.reinit!(mode.vsm_aero; init_aero=false)
    return nothing
end

"""
    n_flap_deflections(sys) -> Int

Number of stations logged in `SysState.flap_angle` for a run of definition `sys`: one
flap deflection δ per station. Sets the `D` type parameter of the run's `SysState`.
"""
n_flap_deflections(sys) = length(sys.spanwise_stations)

"""
    write_flap_deflections!(sys_state, sam)

Write each station's flap deflection δ [rad] of the run `sam` into
`sys_state.flap_angle`, by station index, via [`station_deltas`](@ref). No-op when the
state carries no flap slots (`D == 0`).
"""
function write_flap_deflections!(sys_state, sam)
    isempty(sys_state.flap_angle) && return nothing
    sys_state.flap_angle .= station_deltas(sam)
    return nothing
end

"""
    restore_flap_delta!(mode, sam, wing_idx, sys_state)

Restore each VSM panel's flap deflection `δ` of wing `wing_idx` of `sam` from
`sys_state.flap_angle` (mapped through the panel→station index) when loading a
`SysState` frame, so a replayed frame shows the logged flap state. Default no-op;
[`AeroPressure`](@ref) restores its `PARTICLE_DYNAMICS` panels.
"""
restore_flap_delta!(::AbstractAeroModel, sam, wing_idx, sys_state) = nothing

"""
    restore_live_shape!(mode, sam, wing_idx)

Re-derive the live airfoil shapes of the aero `mode` of wing `wing_idx` of `sam` from a
replayed log frame's structure. Default no-op; [`AeroPressure`](@ref) rebuilds them
when it is on live polars.
"""
restore_live_shape!(::AbstractAeroModel, sam, wing_idx) = nothing

# ==================== shared VSM numerics ==================== #

"""
    circulation(sam, wing_idx) -> AbstractVector

The circulation of wing `wing_idx` of `sam` per VSM panel [m²/s], its slice of the run's
`SysState` `gamma_distribution`: the aerodynamic state its solves start from and leave,
all zero before the first.
"""
function circulation(sam, wing_idx)
    start = sum(n_vsm_panels(sam.aero[other]) for other in 1:(wing_idx - 1); init=0)
    panels = start .+ (1:n_vsm_panels(sam.aero[wing_idx]))
    return view(sam.sysstate.gamma_distribution, panels)
end

"""
    circulation_alpha(sam, wing_idx) -> Vector{SimFloat}

The angle of attack [rad] each VSM panel of wing `wing_idx` of `sam` flies at its
[`circulation`](@ref), in the panels' current apparent wind and the influence its engine
last worked out; each panel's geometric angle before the engine has solved once.
"""
function circulation_alpha(sam, wing_idx)
    body_aero = sam.aero[wing_idx].vsm_aero
    v_ind = zeros(SimFloat, 3, length(body_aero.panels))
    all(iszero, body_aero.AIC) ||
        store_induced_velocity!(v_ind, body_aero, circulation(sam, wing_idx))
    return [SimFloat(first(VortexStepMethod.calculate_relative_alpha_and_relative_velocity(
                panel, v_ind[:, i]))) for (i, panel) in enumerate(body_aero.panels)]
end

"""The number of VSM panels of aero `mode`, none without an engine."""
n_vsm_panels(mode::AbstractVSMAero) = length(mode.vsm_aero.panels)
n_vsm_panels(::AbstractAeroModel) = 0

"""The start of a solve from circulation `gamma`: a copy of it, `nothing` for the solver's
initial distribution where it is all zero."""
circulation_start(gamma) = all(iszero, gamma) ? nothing : copy(gamma)

"""
    solve_circulation!(sam, wing_idx, solver, body_aero; moment_frac=0.1)

Solve the VSM `solver` of wing `wing_idx` of `sam` on `body_aero` from its
[`circulation`](@ref) ([`safe_vsm_solve!`](@ref)) and keep the converged circulation
there.
"""
function solve_circulation!(sam, wing_idx, solver, body_aero; moment_frac=0.1)
    gamma = circulation(sam, wing_idx)
    safe_vsm_solve!(solver, body_aero, circulation_start(gamma); moment_frac)
    gamma .= solver.sol.gamma_distribution
    return nothing
end

"""
    safe_vsm_solve!(solver, body_aero, gamma_init; moment_frac=0.1)

Solve under `throw_on_fail`, so a solve that missed the solver's tolerances or
came back non-finite throws `VortexStepMethod.SolveFailure`. The solve starts from the
circulation `gamma_init`, else, where it is `nothing`, from the solver's configured initial
distribution.
"""
safe_vsm_solve!(solver, body_aero, gamma_init; moment_frac=0.1) =
    VortexStepMethod.solve!(solver, body_aero, gamma_init;
        moment_frac, log=false, throw_on_fail=true)

"""
    solve_and_freeze_circulation!(mode, sam, wing_idx)

Solve wing `wing_idx` of `sam` and freeze the per-refined-panel induced velocity into
`mode.v_ind` and the chord blend weights into `mode.chord_weight`, shared by the
continuous VSM modes. The solve starts from the wing's [`circulation`](@ref) and throws
on a failure before anything frozen is written ([`solve_circulation!`](@ref)).
"""
function solve_and_freeze_circulation!(mode, sam, wing_idx)
    solver = mode.vsm_solver
    body_aero = mode.vsm_aero
    solve_circulation!(sam, wing_idx, solver, body_aero)
    store_induced_velocity!(mode.v_ind, body_aero, solver.lr.gamma_new)
    store_chord_weights!(mode.chord_weight, body_aero)
    return nothing
end

"""
    set_particle_panel_va!(mode, sam, wing_idx, va_point_b)

Set the per-panel apparent wind on the VSM `BodyAerodynamics` of wing `wing_idx` of `sam`
from the body-frame apparent wind of each point, a column of `va_point_b`: each section's
`va` is the mean of its LE/TE point values, refined panels take their parent section's
value. Falls back to the wing's own apparent wind for a section, or the whole wing, the
structural↔panel mapping does not reach.
"""
function set_particle_panel_va!(mode, sam, wing_idx, va_point_b)
    va_b = collect(wing_apparent_wind(sam, wing_idx))
    if isnothing(mode.point_to_vsm_point)
        set_va!(mode.vsm_aero, va_b)
        return nothing
    end
    n_sections = length(mode.vsm_wing.unrefined_sections)
    vsm_point_to_struct = Dict{Tuple{Int64, Symbol}, Int64}()
    for (point_idx, (section_idx, le_or_te)) in mode.point_to_vsm_point
        vsm_point_to_struct[(section_idx, le_or_te)] = point_idx
    end
    section_va = Vector{Vector{Float64}}(undef, n_sections)
    for section_idx in 1:n_sections
        le_idx = get(vsm_point_to_struct, (Int64(section_idx), :LE), nothing)
        te_idx = get(vsm_point_to_struct, (Int64(section_idx), :TE), nothing)
        section_va[section_idx] = isnothing(le_idx) || isnothing(te_idx) ? va_b :
            0.5 * (va_point_b[:, le_idx] + va_point_b[:, te_idx])
    end
    n_panels = length(mode.vsm_aero.panels)
    va_dist = zeros(n_panels, 3)
    mapping = mode.vsm_wing.refined_panel_mapping
    for panel_idx in 1:n_panels
        va_dist[panel_idx, :] .= section_va[mapping[panel_idx]]
    end
    set_va!(mode.vsm_aero, va_dist)
    return nothing
end

"""
    set_refined_panel_va!(mode, sam, wing_idx, va_point_b)

Per-panel apparent wind of wing `wing_idx` of `sam` for the continuous VSM modes: refined
sections from [`strut_inflow_weights`](@ref) (numeric twin of
[`reconstruct_inflow_sym`](@ref)) over the columns of `va_point_b`, each panel the mean
of its two bounding sections. Matches the inflow the symbolic RHS uses, so the frozen
circulation belongs to it. Falls back to the wing's own apparent wind without a
structural↔panel mapping.
"""
function set_refined_panel_va!(mode, sam, wing_idx, va_point_b)
    if isnothing(mode.point_to_vsm_point)
        set_va!(mode.vsm_aero, collect(wing_apparent_wind(sam, wing_idx)))
        return nothing
    end
    column = aero_section_columns(sam, wing_idx)
    point_idxs = wing_points(sam, wing_idx)
    section_va = [
        sum(weight .* view(va_point_b, :, point_idxs[col])
            for (col, weight) in strut_inflow_weights(mode, section, column))
        for section in eachindex(section_interp_caches(mode)[1])
    ]
    n_panels = length(mode.vsm_aero.panels)
    va_dist = zeros(SimFloat, n_panels, 3)
    for panel in 1:n_panels
        va_dist[panel, :] .= 0.5 .* (section_va[panel] .+ section_va[panel + 1])
    end
    set_va!(mode.vsm_aero, va_dist)
    return nothing
end

"""
    vsm_solve_objects(mode, ::Type{T}, shadow_ref) -> (body_aero, solver, vsm_wing)

The VSM solve objects [`vsm_aero_coeffs`](@ref) runs on, selected by the input
eltype `T`. A value pass (`Float64`) uses the real objects of aero `mode`. A
`ForwardDiff` pass feeds `Dual` numbers through the solve to get the Jacobian, but the
real objects' buffers are `Float64` and can't hold Duals — so for a `Dual` eltype we
solve on a `Dual`-typed "shadow" of the solver/aero. The shadow is expensive, so
it is built lazily and cached in `shadow_ref`, keyed by the `Dual` eltype (rebuilt
if it changes); `use_gamma_prev` warm-starts each perturbed solve from the
previous circulation.
"""
vsm_solve_objects(mode, ::Type{Float64}, shadow_ref) =
    (mode.vsm_aero, mode.vsm_solver, mode.vsm_wing)

function vsm_solve_objects(mode, ::Type{T}, shadow_ref) where {T}
    shadow = shadow_ref[]
    if shadow === nothing || eltype(shadow[1]._va) !== T
        shadow = VortexStepMethod.make_dual_shadow(mode.vsm_solver, mode.vsm_aero, T)
        shadow[2].use_gamma_prev = true
        shadow_ref[] = shadow
    end
    body_aero, solver = shadow
    return body_aero, solver, body_aero.wings[1]
end

"""
    vsm_aero_coeffs(mode, y, va_mag, moment_frac, shadow_ref; gamma_init=nothing)
        -> Vector

Run one VSM solve of a wing flown with the VSM aero `mode` at operating-point input
`y = [α, β, ω₁, ω₂, ω₃, θ_twist…]` and return the wind-axis coefficient vector
`[CL, CD, CS, CM₁, CM₂, CM₃, cm_twist…]`, normalised by the dynamic pressure of `va_mag`,
one twist and one `cm` per station, mapped onto the sections and panels of the engine's
`station_sections` and `station_panels`.
`ForwardDiff.Dual`-aware via `vsm_solve_objects`: for a Dual eltype it solves on a
cached dual shadow of the VSM solver, so the same routine yields the Jacobian
under AD.
"""
function vsm_aero_coeffs(mode, y::AbstractVector{T}, va_mag, moment_frac,
                         shadow_ref::Ref; gamma_init=nothing) where {T}
    body_aero_c, solver_c, wing_c = vsm_solve_objects(mode, T, shadow_ref)
    station_sections = mode.engine.station_sections
    station_panels = mode.engine.station_panels
    n_stations = length(station_sections)
    n_unrefined = mode.vsm_wing.n_unrefined_sections

    α = y[1]
    β = y[2]
    ω = MVector{3, T}(y[3], y[4], y[5])

    # Body-frame apparent wind from (α, β, va_mag)
    cα, sα = cos(α), sin(α)
    cβ, sβ = cos(β), sin(β)
    va_b_local = MVector{3, T}(va_mag * cα * cβ,
                               va_mag * sβ,
                               va_mag * sα * cβ)

    theta = zeros(T, n_unrefined)
    for (k, sections) in enumerate(station_sections), section_idx in sections
        theta[section_idx] = y[5 + k]
    end

    if n_unrefined > 0
        VortexStepMethod.unrefined_deform!(wing_c, theta; smooth=false)
        VortexStepMethod.reinit!(body_aero_c; init_aero=false)
    end
    set_va!(body_aero_c, va_b_local, ω)
    safe_vsm_solve!(solver_c, body_aero_c, gamma_init; moment_frac)

    # VSM's own coefficients are scaled by the panel inflow, which carries ω × r.
    sol = solver_c.sol
    q_area = 0.5 * solver_c.density * va_mag^2 * body_aero_c.projected_area
    force_coeffs = sol.force ./ q_area
    cm_body = sol.moment ./ (q_area * body_aero_c.c_ref)

    # Wind-axis basis (VSM): drag∥va, lift=norm(drag×span), side=lift×drag.
    span = SVector(zero(T), one(T), zero(T))
    drag_dir = va_b_local ./ va_mag
    lift_dir = smooth_normalize(cross(drag_dir, span))
    side_dir = cross(lift_dir, drag_dir)

    x = zeros(T, 6 + n_stations)
    x[1] = dot(force_coeffs, lift_dir)
    x[2] = dot(force_coeffs, drag_dir)
    x[3] = dot(force_coeffs, side_dir)
    x[4] = cm_body[1]
    x[5] = cm_body[2]
    x[6] = cm_body[3]
    for (k, panels) in enumerate(station_panels)
        x[6 + k] = sum(sol.moment_dist[panel_idx] for panel_idx in panels;
                       init=zero(T)) / (q_area * body_aero_c.c_ref)
    end
    return x
end

"""
    rigid_aero_baseline!(mode, sam, wing_idx; vsm_min_wind=0.5)

Compute the operating point and baseline wind-axis coefficients of wing `wing_idx` of `sam`
under its VSM aero `mode`: writes the mode's `aero_y` / `aero_x`. Returns the context
(`va_mag`, `moment_frac`, `shadow_ref`, `y0`) the mode-specific reduction
(`refresh_rigid_aero!`) needs for the Jacobian.
"""
function rigid_aero_baseline!(mode, sam, wing_idx; vsm_min_wind=0.5)
    body_idx = wing_body(sam, wing_idx)
    va_b = run_value(sam, Val(:bodies), body_idx, Val(:va_b))
    omega_b = run_value(sam, Val(:bodies), body_idx, Val(:ω_b))
    station_idxs = wing_station_idxs(sam, wing_idx)
    moment_frac = isempty(station_idxs) ? 0.25 :
        run_value(sam, Val(:spanwise_stations), first(station_idxs), Val(:moment_frac))

    va_mag = max(norm(va_b), vsm_min_wind)
    alpha_0 = atan(va_b[3], va_b[1])
    beta_0 = atan(va_b[2], hypot(va_b[1], va_b[3]))
    isfinite(alpha_0) || (alpha_0 = 0.0)
    isfinite(beta_0) || (beta_0 = 0.0)

    # Operating-point input vector y₀ = [α, β, ω, θ_station]
    y0 = mode.aero_y
    y0[1] = alpha_0
    y0[2] = beta_0
    y0[3] = omega_b[1]
    y0[4] = omega_b[2]
    y0[5] = omega_b[3]
    for (k, station_idx) in enumerate(station_idxs)
        y0[5 + k] = run_value(sam, Val(:spanwise_stations), station_idx, Val(:twist))
    end

    shadow_ref = Ref{Any}(nothing)
    gamma = circulation(sam, wing_idx)
    mode.aero_x .= vsm_aero_coeffs(mode, y0, va_mag, moment_frac, shadow_ref;
                                   gamma_init=circulation_start(gamma))
    gamma .= mode.vsm_solver.sol.gamma_distribution
    return (; va_mag, moment_frac, shadow_ref, y0)
end
