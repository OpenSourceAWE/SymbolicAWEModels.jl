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
by. The default normalizes the WING-point point-mass inertia
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

Per-unit-mass inertia of the wing's WING points treated as point masses
(`extra_mass`), normalized by their total mass. Exact under the construction
invariant `wing.mass == sum of WING-point masses` (the constructor
distributes `set.mass` onto the points). With zero total mass, `com_cad` is
the unweighted centroid and `inertia` is `nothing`.
"""
function normalized_point_inertia(wing, points)
    wing_points = [point for point in points
                   if point.type == WING && point.wing_idx == wing.idx]
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

function setup_aero!(mode::AbstractVSMAero, wing, points, twist_surfaces;
                     prn=false)
    require_vsm_engine(mode, wing)
    vsm_wing = wing.vsm_wing
    if wing.dynamics_type == RIGID_DYNAMICS
        # Transform VSM sections CAD → body (with aero z-offset)
        vsm_wing.T_cad_body .= wing.pos_cad
        adjust_vsm_panels_to_origin!(vsm_wing, wing.pos_cad)
        rotate_vsm_sections!(vsm_wing, wing.R_b_to_c')
        vsm_wing.R_cad_body .= wing.R_b_to_c
        apply_aero_z_offset!(vsm_wing, wing.aero_z_offset)
        VortexStepMethod.reinit!(wing.vsm_aero)

        if couples_to_sections(mode) && isempty(wing.twist_surface_idxs)
            auto_create_twist_surfaces!(wing, points, twist_surfaces; prn)
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
        if !isnothing(wing.origin)
            # Transform VSM sections CAD → body (no z-offset for particle)
            vsm_wing.T_cad_body .= wing.pos_cad
            adjust_vsm_panels_to_origin!(vsm_wing, wing.pos_cad)
            rotate_vsm_sections!(vsm_wing, wing.R_b_to_c')
            vsm_wing.R_cad_body .= wing.R_b_to_c
            VortexStepMethod.reinit!(wing.vsm_aero)
        end
        couples_to_sections(mode) &&
            match_aero_sections_to_structure!(wing, points; twist_surfaces)
        isempty(wing.twist_surface_idxs) || empty!(wing.twist_surface_idxs)
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

# ==================== logging / visualization hooks ==================== #
#
# A mode can contribute extra `SysState` log slots (after the structural
# points, before the per-wing position slots) for visualization — VSM modes
# log 4 corners per panel; a custom mode adds methods for its own geometry.

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
writes each section's display quad.
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
panels from twist instead and only skip their slots).
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
