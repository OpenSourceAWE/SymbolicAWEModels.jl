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
    ContinuousPolar(body_aero, coef)

Callable polar for [`ContinuousAero`](@ref), used as a callable flat parameter
`p(panel_idx, α)`: looks up refined panel `panel_idx` and evaluates the VSM
coefficient function `coef` (`calculate_cl`/`calculate_cd`/`calculate_cm`) at
angle of attack `α`. The panel is typeasserted concrete so the polar dispatches
statically with no boxing in the compiled RHS; `ForwardDiff.Dual`-safe in `α`.
"""
struct ContinuousPolar{BA, F}
    body_aero::BA
    coef::F
end
(p::ContinuousPolar)(panel_idx, alpha) = p.coef(
    p.body_aero.panels[round(Int, panel_idx)]::VortexStepMethod.Panel{SimFloat}, alpha)

# ==================== mesh maps ==================== #

"""
    build_mesh_maps!(mode::ContinuousAero)

Size the frozen induced-velocity buffer and copy the refined-section →
strut interpolation (`refined_section_left_idx` / `refined_section_weight`
from the VSM wing): refined section `s` lies at
`w·strut[left] + (1−w)·strut[left+1]`. Also freeze the per-section billow
offset — the refined-section position minus that straight-line interpolation,
in the body frame — which is nonzero only for the `BILLOWING` distribution
(the trailing edge bulges off the strut line). All are constants of the mesh
and are baked into the symbolic equations ([`store_billow_offsets!`](@ref)).
"""
function build_mesh_maps!(mode::ContinuousAero)
    vsm_wing = mode.vsm_wing
    n_panels = Int(vsm_wing.n_panels)
    n_sections = n_panels + 1
    n_struts = Int(vsm_wing.n_unrefined_sections)
    n_struts >= 2 || error(
        "ContinuousAero: need at least 2 unrefined sections, got $n_struts.")
    left = vsm_wing.refined_section_left_idx
    weight = vsm_wing.refined_section_weight
    if length(left) == n_sections && length(weight) == n_sections
        mode.section_left_strut = Int64.(left)
        mode.section_left_weight = SimFloat.(weight)
    elseif n_struts == n_sections
        # No refinement: refined section s is strut s
        mode.section_left_strut = [min(Int64(s), Int64(n_struts - 1))
                                   for s in 1:n_sections]
        mode.section_left_weight = [s < n_sections ? 1.0 : 0.0
                                    for s in 1:n_sections]
    else
        error("ContinuousAero: VSM wing has no refined-section " *
              "interpolation cache ($(length(left)) entries for " *
              "$n_sections refined sections).")
    end
    store_billow_offsets!(mode)
    size(mode.v_ind) == (3, n_panels) ||
        (mode.v_ind = zeros(SimFloat, 3, n_panels))
    body_aero = mode.vsm_aero
    mode.cl = ContinuousPolar(body_aero, VortexStepMethod.calculate_cl)
    mode.cd = ContinuousPolar(body_aero, VortexStepMethod.calculate_cd)
    mode.cm = ContinuousPolar(body_aero, VortexStepMethod.calculate_cm)
    return nothing
end

"""
    store_billow_offsets!(mode::ContinuousAero)

Freeze each refined section's body-frame displacement off the straight strut
line: `refined_pos − (w·strut[left] + (1−w)·strut[left+1])` for both LE and TE.
Zero for in-line distributions (`SPLIT_PROVIDED`, `COSINE`, …); nonzero only
where `BILLOWING` bulges the trailing edge between ribs.
"""
function store_billow_offsets!(mode::ContinuousAero)
    vsm_wing = mode.vsm_wing
    left = mode.section_left_strut
    weight = mode.section_left_weight
    n_sections = length(left)
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
    mode.section_le_offset = le_offset
    mode.section_te_offset = te_offset
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

    point_to_vsm = wing.point_to_vsm_point
    isnothing(point_to_vsm) && error(
        "ContinuousAero: wing $(wing.name) is missing the structural↔panel " *
        "point mapping.")
    column = Dict{Tuple{Int64, Symbol}, Int}()
    for (k, point) in enumerate(points)
        strut_idx, le_or_te = point_to_vsm[point.idx]
        column[(strut_idx, le_or_te)] = k
    end

    n_panels = Int(wing.vsm_wing.n_panels)
    n_struts = Int(wing.vsm_wing.n_unrefined_sections)
    left = mode.section_left_strut
    lweight = mode.section_left_weight
    length(left) == n_panels + 1 || error(
        "ContinuousAero: mesh maps not built for wing $(wing.name).")
    spanwise = collect(SimFloat, wing.vsm_wing.spanwise_direction)
    scale = 1.0 + (isfinite(wing.aero_scale_chord) ?
        wing.aero_scale_chord : AERO_SCALE_CHORD)

    strut_le = [collect(connectors.point_pos[:, column[(s, :LE)]])
                for s in 1:n_struts]
    strut_te = [collect(connectors.point_pos[:, column[(s, :TE)]])
                for s in 1:n_struts]
    strut_va = [0.5 * (collect(connectors.va[:, column[(s, :LE)]]) +
                       collect(connectors.va[:, column[(s, :TE)]]))
                for s in 1:n_struts]
    strut_rho = [0.5 * (connectors.rho[column[(s, :LE)]] +
                        connectors.rho[column[(s, :TE)]])
                 for s in 1:n_struts]

    le_offset = mode.section_le_offset
    te_offset = mode.section_te_offset
    interp(values, s) = lweight[s] * values[left[s]] +
                        (1.0 - lweight[s]) * values[left[s] + 1]
    # Add the frozen billow offset (body frame) to the straight strut line.
    sec_le = [interp(strut_le, s) .+ le_offset[:, s] for s in 1:(n_panels + 1)]
    sec_te = [interp(strut_te, s) .+ te_offset[:, s] for s in 1:(n_panels + 1)]
    sec_va = [interp(strut_va, s) for s in 1:(n_panels + 1)]
    sec_rho = [interp(strut_rho, s) for s in 1:(n_panels + 1)]

    eqs, panel_vars, panel_force, panel_couple = build_panel_force_eqs(
        sec_le, sec_te, sec_va, sec_rho, vind_p, cl, cd, cm, spanwise, scale)
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

# ==================== refresh ==================== #

"""
    refresh_particle_aero!(::ContinuousAero, wing, points, va_point_b_vals;
                           vsm_min_wind=0.5)

Refresh: update the VSM geometry from the structure, set the per-panel apparent
wind, run the full `VortexStepMethod.solve!` (via [`safe_vsm_solve!`](@ref)), and
freeze the per-refined-panel induced velocity ([`store_induced_velocity!`](@ref)).
The forces themselves are re-derived symbolically each RHS step, so `calc_forces!`
is only run so the shared solve path also populates `sol` (`alpha_dist`,
`f_body_3D`) rather than leaving it stale. Below `vsm_min_wind` the induced
velocity is zeroed; the symbolic forces remain live (and vanish with the dynamic
pressure).
"""
function refresh_particle_aero!(mode::ContinuousAero, wing, points,
                                va_point_b_vals; vsm_min_wind=0.5)
    if norm(wing.va_b) < vsm_min_wind
        fill!(mode.v_ind, 0.0)
        return nothing
    end

    update_vsm_wing_from_structure!(wing, points)
    set_particle_panel_va!(wing, va_point_b_vals)

    solver = wing.vsm_solver
    body_aero = wing.vsm_aero
    if !safe_vsm_solve!(solver, body_aero, solver.sol.gamma_distribution)
        throw(AssertionError(
            "ContinuousAero VSM solve failed (non-converged or non-finite) " *
            "on wing $(wing.idx)"))
    end
    gamma = solver.lr.gamma_new
    if isnothing(solver.sol.gamma_distribution)
        solver.sol.gamma_distribution = copy(gamma)
    else
        solver.sol.gamma_distribution .= gamma
    end
    store_induced_velocity!(mode.v_ind, body_aero, gamma)
    return nothing
end

# ==================== visualization ==================== #

"""
    reconstruct_sections_b(mode::ContinuousAero, wing, points)

Body-frame refined-section LE/TE positions exactly as the force model builds
them: the live strut points interpolated by the frozen mesh weights plus the
frozen billow offset. Mirrors the symbolic geometry in [`aero_component`](@ref)
so the plotted panels are the ones the dynamics actually use (and a wrong
billow offset shows up visually).
"""
function reconstruct_sections_b(mode::ContinuousAero, wing, points)
    point_to_vsm = wing.point_to_vsm_point
    column = Dict{Tuple{Int64, Symbol}, Int}()
    for (k, point) in enumerate(points)
        strut_idx, le_or_te = point_to_vsm[point.idx]
        column[(strut_idx, le_or_te)] = k
    end
    rot_w_to_b = (wing.R_b_to_w::Matrix{SimFloat})'
    n_struts = Int(wing.vsm_wing.n_unrefined_sections)
    strut_le = [rot_w_to_b * (points[column[(s, :LE)]].pos_w - wing.pos_w)
                for s in 1:n_struts]
    strut_te = [rot_w_to_b * (points[column[(s, :TE)]].pos_w - wing.pos_w)
                for s in 1:n_struts]
    left = mode.section_left_strut
    weight = mode.section_left_weight
    interp(values, s) = weight[s] .* values[left[s]] .+
                        (1.0 - weight[s]) .* values[left[s] + 1]
    sec_le = [interp(strut_le, s) .+ mode.section_le_offset[:, s]
              for s in eachindex(left)]
    sec_te = [interp(strut_te, s) .+ mode.section_te_offset[:, s]
              for s in eachindex(left)]
    return sec_le, sec_te
end

"""
    write_aero_log_points!(mode::ContinuousAero, wing, sys_struct, sys_state,
                           point_idx, zoom)

Log the panel corners the force model reconstructs (strut interpolation +
frozen billow offset), not the raw VSM mesh, so the plot shows the geometry
the dynamics use.
"""
function write_aero_log_points!(mode::ContinuousAero, wing, sys_struct,
                                sys_state, point_idx, zoom)
    points = wing_points(sys_struct, wing)
    sec_le, sec_te = reconstruct_sections_b(mode, wing, points)
    rot_b_to_w = wing.R_b_to_w::Matrix{SimFloat}
    n_panels = Int(wing.vsm_wing.n_panels)
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
