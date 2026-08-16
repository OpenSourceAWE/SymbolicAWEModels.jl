# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

"""
    ContinuousAero()

Frozen-circulation VSM aerodynamics with live symbolic force assembly
(`PARTICLE_DYNAMICS` only). The VSM solver runs every `vsm_interval` steps and
solves only the circulation distribution (`VortexStepMethod.solve_base!`); the
resulting per-refined-panel induced velocity is frozen. Each RHS step then
evaluates the `calc_forces!` chain symbolically per refined panel — geometry
interpolated from the live strut points with the frozen mesh weights,
effective angle of attack, polar lookups, and lift/drag directions all live —
capturing aerodynamic damping between refreshes. Carries a
[`VSMEngine`](@ref); the no-arg form is the engine-less marker filled in
during wing construction.
"""
mutable struct ContinuousAero{E} <: AbstractVSMAero
    engine::Union{Nothing, E}
    "Frozen body-frame induced velocity per refined panel (3 × n_panels)."
    v_ind::Matrix{SimFloat}
    "Left strut (unrefined section) of each refined section (n_panels + 1)."
    section_left_strut::Vector{Int64}
    "Left-strut weight: section = w·strut[left] + (1−w)·strut[left+1]."
    section_left_weight::Vector{SimFloat}
    "Frozen body-frame billow offset of each refined-section LE off the strut line (3 × n_sections)."
    section_le_offset::Matrix{SimFloat}
    "Frozen body-frame billow offset of each refined-section TE off the strut line (3 × n_sections)."
    section_te_offset::Matrix{SimFloat}
    "Polar callables `(panel_idx, α)` for cl/cd/cm, set in `build_mesh_maps!`; read as callable flat params."
    cl::Any
    cd::Any
    cm::Any
    ContinuousAero{E}(engine, v_ind, section_left_strut, section_left_weight,
        section_le_offset, section_te_offset, cl, cd, cm) where {E} =
        new{E}(engine, v_ind, section_left_strut, section_left_weight,
            section_le_offset, section_te_offset, cl, cd, cm)
end

ContinuousAero() =
    ContinuousAero{VSMEngine}(nothing, zeros(SimFloat, 3, 0), Int64[], SimFloat[],
                   zeros(SimFloat, 3, 0), zeros(SimFloat, 3, 0),
                   nothing, nothing, nothing)
attach_engine!(mode::ContinuousAero, engine::VSMEngine) =
    ContinuousAero{typeof(engine)}(engine, mode.v_ind, mode.section_left_strut,
        mode.section_left_weight, mode.section_le_offset, mode.section_te_offset,
        mode.cl, mode.cd, mode.cm)

is_builtin_aero(::ContinuousAero) = true
aero_mode_tag(::ContinuousAero) = "cont"

"""
    aero_hash_id(mode::ContinuousAero)

The frozen mesh-interpolation weights and billow offsets are baked into the
generated equations, so they are structural and enter the model-cache hash.
"""
aero_hash_id(mode::ContinuousAero) =
    (mode.section_left_strut, round.(mode.section_left_weight; digits=8),
     round.(mode.section_le_offset; digits=8),
     round.(mode.section_te_offset; digits=8))

# ==================== polar callable ==================== #
"""
    AeroHandle(body_aero)

Mutable holder for the live `BodyAerodynamics` a polar reads, and the reason a
serialized model does not carry one. A polar is a callable parameter, so its
concrete type is baked into the parameter store and into what the solver
specialized on and has to survive a round trip — but its contents are dead weight,
because `sync_params!` rebinds every polar from the `SystemStructure` before the
problem is evaluated. Serializing the holder empty keeps the type and drops the
wing's whole Cp/cf surface tables, which were 2.4 GB of bin on a large kite. It is
undefined between `deserialize` and that sync, and reading it there is a bug.
"""
mutable struct AeroHandle{BA}
    body_aero::BA
    AeroHandle{BA}() where {BA} = new{BA}()
    AeroHandle(body_aero::BA) where {BA} = new{BA}(body_aero)
end

function Serialization.serialize(s::Serialization.AbstractSerializer,
                                 handle::AeroHandle)
    Serialization.serialize_type(s, typeof(handle))
    return nothing
end

Serialization.deserialize(s::Serialization.AbstractSerializer,
                          ::Type{AeroHandle{BA}}) where {BA} = AeroHandle{BA}()

"""
    ContinuousPolar(body_aero, coef)

Callable polar for [`ContinuousAero`](@ref), used as a callable flat parameter
`p(panel_idx, α)`: looks up refined panel `panel_idx` and evaluates the VSM
coefficient function `coef` (`calculate_cl`/`calculate_cd`/`calculate_cm`) at
angle of attack `α`. The panel is typeasserted concrete so the polar dispatches
statically with no boxing in the compiled RHS; `ForwardDiff.Dual`-safe in `α`. The
aerodynamics is held through an [`AeroHandle`](@ref) so that serializing the polar
does not serialize the wing.
"""
struct ContinuousPolar{BA, F}
    handle::AeroHandle{BA}
    coef::F
end
ContinuousPolar(body_aero, coef) = ContinuousPolar(AeroHandle(body_aero), coef)

(p::ContinuousPolar)(panel_idx, alpha) = p.coef(
    p.handle.body_aero.panels[round(Int, panel_idx)]::VortexStepMethod.Panel{SimFloat},
    alpha)
(p::ContinuousPolar)(panel_idx, alpha, delta) = p.coef(
    p.handle.body_aero.panels[round(Int, panel_idx)]::VortexStepMethod.Panel{SimFloat},
    alpha, delta)

# ==================== mesh maps ==================== #

"""
    build_mesh_maps!(mode::ContinuousAero)

Size the frozen induced-velocity buffer, set the polar callables, and freeze the
refined-section → strut interpolation caches via [`build_section_interp`](@ref).
"""
function build_mesh_maps!(mode::ContinuousAero)
    vsm_wing = mode.vsm_wing
    n_panels = Int(vsm_wing.n_panels)
    mode.section_left_strut, mode.section_left_weight,
        mode.section_le_offset, mode.section_te_offset =
        build_section_interp(vsm_wing)
    size(mode.v_ind) == (3, n_panels) ||
        (mode.v_ind = zeros(SimFloat, 3, n_panels))
    body_aero = mode.vsm_aero
    mode.cl = ContinuousPolar(body_aero, VortexStepMethod.calculate_cl)
    mode.cd = ContinuousPolar(body_aero, VortexStepMethod.calculate_cd)
    mode.cm = ContinuousPolar(body_aero, VortexStepMethod.calculate_cm)
    return nothing
end

"""
    setup_aero!(mode::ContinuousAero, wing, points, twist_surfaces; prn=false)

The generic VSM particle setup plus the [`ContinuousAero`](@ref) mesh maps
([`build_mesh_maps!`](@ref)).
"""
function setup_aero!(mode::ContinuousAero, wing, points, twist_surfaces;
                     prn=false)
    wing.dynamics_type == PARTICLE_DYNAMICS || error(
        "ContinuousAero supports PARTICLE_DYNAMICS wings only; wing " *
        "$(wing.name) is $(wing.dynamics_type).")
    invoke(setup_aero!, Tuple{AbstractVSMAero, Any, Any, Any},
           mode, wing, points, twist_surfaces; prn)
    mode.vsm_wing.spanwise_distribution == VortexStepMethod.BILLOWING || error(
        "ContinuousAero requires the BILLOWING spanwise distribution so the " *
        "refined panels carry the canopy billow shape; wing $(wing.name) uses " *
        "$(mode.vsm_wing.spanwise_distribution). Set " *
        "spanwise_panel_distribution: BILLOWING in the VSM settings.")
    build_mesh_maps!(mode)
    return nothing
end

"""
    remake_aero!(mode::ContinuousAero, wing, set, vsm_set, points,
                 twist_surfaces)

The generic VSM remake plus a rebuild of the mesh maps (the VSM wing geometry
objects are replaced, invalidating the panel indexing).
"""
function remake_aero!(mode::ContinuousAero, wing, set, vsm_set, points,
                      twist_surfaces)
    invoke(remake_aero!, Tuple{AbstractVSMAero, Any, Any, Any, Any, Any},
           mode, wing, set, vsm_set, points, twist_surfaces)
    build_mesh_maps!(mode)
    return nothing
end

# ==================== equation builder ==================== #

"""
    aero_component(mode::ContinuousAero, wing::ParticleWing, sys_struct; name)

Symbolic per-refined-panel re-expression of `VortexStepMethod.calc_forces!` on
the `PARTICLE_DYNAMICS` connector contract. Refined-section positions,
apparent wind, and density are interpolated from the live strut points with
the frozen mesh weights; per panel, the axes, chord, width, effective angle of
attack (live apparent wind + frozen induced velocity), polar coefficients, and
lift/drag directions are symbolic variables of the component (observable
through the integrator, e.g. `aero_1.alpha`). Each panel force acts on the
quarter-chord line (75 % LE / 25 % TE) with the pitching moment as an LE/TE
force couple, distributed to the bounding struts by the mesh weights.
"""
function aero_component(mode::ContinuousAero, wing::ParticleWing, sys_struct;
                        name, params=nothing)
    wing_idx = wing.idx
    vind_p = params.wings[wing_idx].aero.v_ind
    cl = params.wings[wing_idx].aero.cl   # callable flat params: `cl(panel_idx, α)`
    cd = params.wings[wing_idx].aero.cd
    cm = params.wings[wing_idx].aero.cm

    points = wing_points(sys_struct, wing)
    num_points = length(points)
    connectors = particle_aero_connectors(num_points)

    column = aero_section_columns(wing, points)
    n_panels = Int(wing.vsm_wing.n_panels)
    left = mode.section_left_strut
    lweight = mode.section_left_weight
    length(left) == n_panels + 1 || error(
        "ContinuousAero: mesh maps not built for wing $(wing.name).")
    spanwise = collect(SimFloat, wing.vsm_wing.spanwise_direction)
    scale = 1.0 + (isfinite(wing.aero_scale_chord) ?
        wing.aero_scale_chord : AERO_SCALE_CHORD)

    sec_le, sec_te =
        reconstruct_sections_sym(mode, wing, points, connectors, column)
    sec_va, sec_rho = reconstruct_inflow_sym(mode, wing, connectors, column)

    orient = panel_span_signs(wing, spanwise)
    eqs, panel_vars, panel_force, panel_couple = build_panel_force_eqs(
        sec_le, sec_te, sec_va, sec_rho, vind_p, cl, cd, cm, spanwise, scale, orient)
    vars = particle_unknowns(connectors)
    append!(vars, panel_vars)

    point_force = [zeros(Num, 3) for _ in 1:num_points]
    for i in 1:n_panels
        force = collect(panel_force[:, i])
        couple = collect(panel_couple[:, i])
        force_le = 0.75 * force + couple
        force_te = 0.25 * force - couple
        for s in (i, i + 1), (strut, w) in
                ((left[s], lweight[s]), (left[s] + 1, 1.0 - lweight[s]))
            w == 0.0 && continue
            kle = column[(strut, :LE)]
            kte = column[(strut, :TE)]
            point_force[kle] = point_force[kle] + (0.5 * w) * force_le
            point_force[kte] = point_force[kte] + (0.5 * w) * force_te
        end
    end

    for k in 1:num_points
        eqs = [eqs; connectors.point_force[:, k] ~ point_force[k]]
    end
    return System(eqs, t, vars, [vind_p, cl, cd, cm]; name)
end

# ==================== per-panel decomposition ==================== #

supports_panel_decomposition(::ContinuousAero) = true

"""
    aero_scatter_entries(mode::ContinuousAero, wing, points)

The strut couple: each panel's force acts on the quarter-chord line, so its bounding
sections' struts take `0.75·force + couple` at the LE station and `0.25·force −
couple` at the TE station, halved between the two sections and split by the mesh
weights.
"""
function aero_scatter_entries(mode::ContinuousAero, wing, points)
    column = aero_section_columns(wing, points)
    left, weight, _, _ = section_interp_caches(mode)
    totals = Dict{Tuple{Int, Int}, Vector{SimFloat}}()
    for panel in 1:(length(left) - 1), section in (panel, panel + 1)
        for (strut, share) in ((left[section], weight[section]),
                               (left[section] + 1, 1.0 - weight[section]))
            share == 0.0 && continue
            scatter_totals!(totals, panel, column[(strut, :LE)],
                            0.375 * share, 0.5 * share)
            scatter_totals!(totals, panel, column[(strut, :TE)],
                            0.125 * share, -0.5 * share)
        end
    end
    return scatter_entry_list(totals)
end

# ==================== refresh ==================== #

"""
    refresh_particle_aero!(::ContinuousAero, wing, points, va_point_b_vals;
                           vsm_min_wind=0.5)

Refresh: update the VSM geometry from the structure, set the per-panel apparent wind
([`set_refined_panel_va!`](@ref)), solve and freeze the induced velocity. The forces
are re-derived symbolically each RHS step; `calc_forces!` runs only so `sol`
(`alpha_dist`, `f_body_3D`) is not left stale. Below `vsm_min_wind` the induced
velocity is zeroed.
"""
function refresh_particle_aero!(mode::ContinuousAero, wing, points,
                                va_point_b_vals; vsm_min_wind=0.5)
    if norm(wing.va_b) < vsm_min_wind
        fill!(mode.v_ind, 0.0)
        return nothing
    end

    update_vsm_wing_from_structure!(wing, points)
    set_refined_panel_va!(mode, wing, points, va_point_b_vals)
    solve_and_freeze_circulation!(mode, wing)
    return nothing
end

# ==================== visualization ==================== #

write_aero_log_points!(mode::ContinuousAero, wing, sys_struct, sys_state,
                       point_idx, zoom) =
    write_live_aero_log_points!(mode, wing, sys_struct, sys_state,
                                point_idx, zoom)
