# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Code shared by all aero modes: the dispatch interface (generic functions +
# abstract-type defaults), the MTK connector scaffolding the `aero_component`
# builders share, and the low-frequency VSM refresh orchestrator + numerics.
# The concrete aero modes live one-per-file alongside this (none/direct/
# linearized/plate.jl) and add their own methods. The abstract aero types
# themselves (AbstractAeroModel, AbstractVSMAero, VSMEngine) live in
# system_structure/types.jl because the `Wing.aero` field references them.

# ==================== interface: capability traits ==================== #

"""
    has_vsm_engine(mode::AbstractAeroModel) -> Bool

`true` if the mode carries a [`VSMEngine`](@ref) (VSM geometry + linearization
state). Drives the wing-level [`is_vsm`](@ref) helper and every VSM-geometry loop.
"""
has_vsm_engine(::AbstractAeroModel) = false
has_vsm_engine(::AbstractVSMAero) = true

"""
    couples_to_sections(mode::AbstractAeroModel) -> Bool

`true` if the mode needs per-section twist surfaces (auto-creation and
aero-section matching). VSM modes do, except [`AeroNone`](@ref).
"""
couples_to_sections(::AbstractAeroModel) = false
couples_to_sections(::AbstractVSMAero) = true

"""
    provides_aero_override(mode::AbstractAeroModel) -> Bool

`true` if the mode supplies frozen body-frame force/moment overrides read by the
compiled RHS (the stored-force path of [`AeroDirect`](@ref)).
"""
provides_aero_override(::AbstractAeroModel) = false

"""
    stores_point_force(mode::AbstractAeroModel) -> Bool

`true` if a WING point's `aero_force_b` is meaningful for this mode.
[`AeroNone`](@ref) returns `false` (it produces no force).
"""
stores_point_force(::AbstractAeroModel) = true

"""
    exposes_aero_input(mode::AbstractAeroModel) -> Bool

`true` if the mode's subsystem exposes an `aero_input` connector (the linearized
operating-point vector of [`AeroLinearized`](@ref)), which is logged as wing state.
"""
exposes_aero_input(::AbstractAeroModel) = false

# ==================== interface: required cache tag ==================== #

"""
    aero_mode_tag(mode::AbstractAeroModel) -> String

Short identifier for the mode in the compiled-model cache filename. Required: the
`AbstractAeroModel` fallback errors, so every aero mode must declare its own tag
(no silent default that could collide two distinct modes on one cache file).
Built-ins: `"lin"`, `"dir"`, `"none"`, `"plate"`.
"""
aero_mode_tag(mode::AbstractAeroModel) = error(
    "aero_mode_tag is not defined for aero mode $(typeof(mode)); " *
    "every aero mode must provide its own cache tag.")

# ==================== interface: diagnostics ==================== #

"""
    calc_aoa(mode::AbstractAeroModel, wing) -> SimFloat

Angle of attack [rad] for `wing` under aero `mode`. Defaults to `NaN`
(undefined); VSM modes read the mid-span geometric AoA (wrapped to [-π, π]) and
[`AeroPlate`](@ref) derives it from the body-frame apparent wind. [`AeroNone`](@ref)
produces no solve, so it stays `NaN`.
"""
calc_aoa(::AbstractAeroModel, wing) = SimFloat(NaN)
function calc_aoa(mode::AbstractVSMAero, wing)
    dist = mode.vsm_solver.sol.alpha_geometric_dist
    n = length(dist)
    return mod(dist[n ÷ 2 + n % 2] + π, 2π) - π
end

"""
    calc_side_slip(mode::AbstractAeroModel, wing) -> SimFloat

Side-slip angle [rad] from the body-frame apparent wind. Mode-independent
(pure geometry), so the default serves every mode.
"""
calc_side_slip(::AbstractAeroModel, wing) =
    atan(wing.va_b[2], hypot(wing.va_b[1], wing.va_b[3]))

# ==================== connector scaffolding ==================== #
#
# A wing carries an `aero::AbstractAeroModel`; `aero_component(mode, …)` is
# dispatched on its type and returns a `System` whose connectors are fixed by
# the wing's `dynamics_type`:
#
#   RIGID_DYNAMICS (num_twist_surfaces = length(wing.twist_surface_idxs)):
#     inputs:  va[1:3], rho, R_b_w[1:3,1:3], omega[1:3],
#              twist[1:num_twist_surfaces], twist_vel[1:num_twist_surfaces]
#     outputs: force[1:3], moment[1:3], twist_moment[1:num_twist_surfaces]
#
#   PARTICLE_DYNAMICS (num_points = number of WING points):
#     inputs:  point_pos[1:3,1:np], point_vel[1:3,1:np], va[1:3,1:np], rho[1:np]
#     outputs: point_force[1:3,1:np]
#
# Everything is in the wing body frame. The wiring layer (aero_eqs!) drives the
# inputs and reads the outputs. Connectors are declared as array variables and
# passed to `System` unflattened, so input connectors a mode ignores still exist
# for the wiring layer to bind.

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

function particle_aero_connectors(num_points::Int)
    @variables begin
        point_pos(t)[1:3, 1:num_points]
        point_vel(t)[1:3, 1:num_points]
        va(t)[1:3, 1:num_points]
        rho(t)[1:num_points]
        point_force(t)[1:3, 1:num_points]
    end
    return (; point_pos, point_vel, va, rho, point_force)
end

particle_unknowns(connectors) =
    Any[connectors.point_pos, connectors.point_vel, connectors.va,
        connectors.rho, connectors.point_force]

function wing_points(sys_struct, wing)
    return [point for point in sys_struct.points
            if point.type == WING && point.wing_idx == wing.idx]
end

"""
    aero_component(mode::AbstractAeroModel, sys_struct, wing_idx; name) -> System

Build the aero subsystem for `sys_struct.wings[wing_idx]`, selected by dispatch
on the wing's `aero` model. Returns a `System` exposing the connectors fixed by
the wing's `dynamics_type` (see above). Add a method on a custom
`AbstractAeroModel` subtype to plug in your own aerodynamics.
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
#
# The low-frequency VSM-update path (every `vsm_interval` steps). `refresh_aero!`
# orchestrates; per-mode work is dispatched on the wing's aero mode via
# `refresh_rigid_aero!` / `refresh_particle_aero!` (in the per-mode files). This
# is NOT the compiled RHS, so dynamic dispatch on the abstract `wing.aero` field
# is free.

"""
    refresh_aero!(sam::SymbolicAWEModel, prob::ProbWithAttributes,
                  integ=sam.integrator; vsm_min_wind=0.5)

Refresh each wing's aerodynamic state, dispatching on the wing's aero mode
([`refresh_rigid_aero!`](@ref) / [`refresh_particle_aero!`](@ref)). Runs on the
low-frequency VSM-update schedule (`vsm_interval`), not the compiled RHS.

**RIGID_DYNAMICS VSM modes:** compute wind-axis coefficients (CL, CD, CS, CM, cm)
at the operating point, plus the `ForwardDiff` Jacobian over `[α, β, ω₁, ω₂, ω₃,
θ_twist…]` (`AeroLinearized`) or the frozen forces (`AeroDirect`).

**PARTICLE_DYNAMICS VSM modes:** full nonlinear VSM solve with per-point force
distribution. Non-VSM modes (`AeroNone`/`AeroPlate`) are no-ops.
"""
function refresh_aero!(sam::SymbolicAWEModel,
                       prob::ProbWithAttributes,
                       integ=sam.integrator;
                       vsm_min_wind=0.5)
    wings = sam.sys_struct.wings
    twist_surfaces = sam.sys_struct.twist_surfaces
    points = sam.sys_struct.points

    length(wings) == 0 && return nothing

    for wing in wings
        wing.dynamics_type == RIGID_DYNAMICS || continue
        refresh_rigid_aero!(wing.aero, wing, sam.am, twist_surfaces;
                            vsm_min_wind)
    end

    any(w.dynamics_type === PARTICLE_DYNAMICS for w in wings) ||
        return nothing
    point_state = prob.get_point_state(integ)
    va_point_b_vals = point_state[4]
    for wing in wings
        wing.dynamics_type == PARTICLE_DYNAMICS || continue
        refresh_particle_aero!(wing.aero, wing, points, va_point_b_vals;
                               vsm_min_wind)
    end

    nothing
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
    refresh_particle_aero!(mode, wing, points, va_point_b_vals; vsm_min_wind=0.5)

Refresh a `PARTICLE_DYNAMICS` wing's aero state, dispatched on its aero `mode`:
- `AeroNone` / any non-VSM mode → no-op (fallback).
- `AeroDirect` → full nonlinear VSM solve with per-section apparent wind, then
  distribute panel forces onto the wing's structural points
  ([`distribute_panel_forces_to_points!`](@ref)); below `vsm_min_wind` the point
  forces are zeroed.
- `AeroLinearized` → unsupported (errors).
"""
refresh_particle_aero!(::AbstractAeroModel, wing, points, va_point_b_vals;
                       vsm_min_wind=0.5) = nothing

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

function remake_aero!(::AbstractVSMAero, wing, set, vsm_set, points,
                      twist_surfaces)
    vsm_set isa VortexStepMethod.VSMSettings || error(
        "remake_aero!: VSM wing $(wing.idx) needs a VSMSettings, " *
        "got $(typeof(vsm_set)).")
    wing.vsm_wing = create_vsm_wing(set, vsm_set;
        prn=false, sort_sections=false)
    wing.vsm_aero = VortexStepMethod.BodyAerodynamics([wing.vsm_wing])
    wing.vsm_solver = VortexStepMethod.Solver(wing.vsm_aero, vsm_set)

    # Transform sections CAD → body frame (matches the SystemStructure constructor)
    vsm_wing = wing.vsm_wing
    vsm_wing.T_cad_body .= wing.pos_cad
    adjust_vsm_panels_to_origin!(vsm_wing, wing.pos_cad)
    rotate_vsm_sections!(vsm_wing, wing.R_b_to_c')
    vsm_wing.R_cad_body .= wing.R_b_to_c
    if wing.dynamics_type != PARTICLE_DYNAMICS
        apply_aero_z_offset!(vsm_wing, wing.aero_z_offset)
    end
    VortexStepMethod.reinit!(wing.vsm_aero)

    match_aero_sections_to_structure!(wing, points; twist_surfaces)

    if wing.dynamics_type == RIGID_DYNAMICS && !isempty(wing.twist_surface_idxs)
        compute_spatial_twist_surface_mapping!(wing, twist_surfaces, points)
    end
    if wing.dynamics_type == PARTICLE_DYNAMICS &&
       !isnothing(wing.point_to_vsm_point)
        wing_point_idxs = collect(keys(something(wing.point_to_vsm_point)))
        wing_pts = [points[idx] for idx in wing_point_idxs]
        wing.point_to_vsm_point =
            build_point_to_vsm_point_mapping(wing_pts, wing)
    end
    return nothing
end

"""
    validate_aero_structure(mode, wing, points; prn=false)

Check structural invariants the mode's compiled equations rely on (run at build).
Default no-op; VSM `PARTICLE_DYNAMICS` wings verify the structural↔panel point
mapping exists, covers every WING point, and matches `2 × n_sections` points.
"""
validate_aero_structure(::AbstractAeroModel, wing, points; prn=false) = nothing

function validate_aero_structure(::AbstractVSMAero, wing, points; prn=false)
    wing.dynamics_type == PARTICLE_DYNAMICS || return nothing
    @assert !isnothing(wing.point_to_vsm_point) "PARTICLE_DYNAMICS wing $(wing.idx) missing point_to_vsm_point mapping"

    wing_point_idxs = [p.idx for p in points if p.type == WING && p.wing_idx == wing.idx]
    for point_idx in wing_point_idxs
        @assert haskey(wing.point_to_vsm_point, point_idx) "PARTICLE_DYNAMICS wing $(wing.idx) missing mapping for point $(point_idx)"
    end

    n_sections = length(wing.vsm_wing.unrefined_sections)
    @assert length(wing_point_idxs) == 2 * n_sections "PARTICLE_DYNAMICS wing $(wing.idx): expected $(2*n_sections) points for $(n_sections) sections, got $(length(wing_point_idxs))"

    prn && println("✓ PARTICLE_DYNAMICS wing $(wing.idx) validated: $(length(wing_point_idxs)) points, $(n_sections) sections, $(length(wing.vsm_aero.panels)) panels")
    return nothing
end

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
NaN/Inf-guarded `solve!`. Checks both Dual value and partials. On
non-finite or non-converged result, zero gamma and return `false`.
"""
function safe_vsm_solve!(solver, body_aero,
                          gamma_init=nothing; moment_frac=0.1)
    if isnothing(gamma_init)
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
        if !isnothing(solver.sol.gamma_distribution)
            fill!(solver.sol.gamma_distribution, 0)
        end
        return false
    end
    return true
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
        throw(AssertionError("VSM solve failed (non-converged or non-finite) on wing $(wing.idx) [eltype=$T]"))
    end

    sol = solver_c.sol
    force_coeffs = sol.force_coeffs
    cm_body = sol.moment_coeffs
    moment_coeff_unrefined = sol.moment_coeff_unrefined_dist

    # Wind-axis basis (matches VSM): drag along va,
    # lift = normalize(drag × span), side = lift × drag.
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
