# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

"""
Plot the linearised VSM response (ForwardDiff tangents) around
the current operating point with GLMakie.

For every input direction the panel shows the line
  x[i] = x0[i] + J[i, j] · d
on `d ∈ [-0.1, +0.1]`, where `J = ForwardDiff.jacobian` of the
wind-axis coefficient vector at the operating point.

Inputs:
    α, β   [rad]   apparent-flow direction perturbations
    |va|   [—]     fractional wind-speed perturbation
    ω₁..₃  [rad/s] body angular rate perturbations
    θ_g    [rad]   per-group twist perturbations

Outputs (wind-axis coefficients):
    CL, CD, CS, CM₁..₃, cm_g per group
"""

using Pkg
if Base.active_project() != joinpath(@__DIR__, "Project.toml")
    Pkg.activate(@__DIR__)
end

using GLMakie
using LinearAlgebra
using Printf
using ForwardDiff
using StaticArrays
using KiteUtils: init!
using SymbolicAWEModels, VortexStepMethod

# ─── Model setup ─────────────────────────────────────
MODEL_NAME = "2plate_kite"

pkg_root = dirname(@__DIR__)
set_data_path(joinpath(pkg_root, "data", MODEL_NAME))

struc_yaml = joinpath(
    get_data_path(), "rigid_structural_geometry.yaml")

set = Settings("system.yaml")
vsm_set = VortexStepMethod.VSMSettings(
    joinpath(get_data_path(), "vsm_settings.yaml");
    data_prefix=false)

sys = load_sys_struct_from_yaml(struc_yaml;
    system_name=MODEL_NAME, set, vsm_set)
sam = SymbolicAWEModel(set, sys)
init!(sam; remake=false)

# ─── Extract VSM objects ─────────────────────────────
wing = sam.sys_struct.wings[1]
solver = wing.vsm_solver
body_aero = wing.vsm_aero
vsm_wing = wing.vsm_wing
groups = sam.sys_struct.groups

@info "Solver" solver_type=solver.solver_type use_gamma_prev=solver.use_gamma_prev

va_b_0 = copy(wing.va_b)
omega_b_0 = copy(wing.ω_b)
n_unrefined = vsm_wing.n_unrefined_sections

moment_frac =
    groups[first(wing.group_idxs)].moment_frac

# Current twist per unrefined section
theta_0 = zeros(n_unrefined)
for gidx in wing.group_idxs
    g = groups[gidx]
    for ui in g.unrefined_section_idxs
        theta_0[ui] = g.twist
    end
end

# Operating point in (alpha, beta, va_mag) space
va_mag_0 = norm(va_b_0)
alpha_0 = atan(va_b_0[3], va_b_0[1])
beta_0 = asin(clamp(va_b_0[2] / va_mag_0, -1, 1))

@info "Operating point" alpha_deg=rad2deg(alpha_0) beta_deg=rad2deg(
    beta_0) va_mag_0 omega_b_0 theta_0

# ─── ForwardDiff coefficient evaluator ───────────────
#
# Pure function `d → x` on a perturbation vector `d` (one entry
# per input). For `Float64` it drives the live wing's VSM
# objects; for `Dual` it uses a `make_dual_shadow` so
# derivatives propagate correctly.

shadow_ref = Ref{Any}(nothing)

function coeffs_at_perturbation(d::AbstractVector{T},
        wing_grp_idxs) where {T}
    if T === Float64
        body_aero_c = body_aero
        solver_c = solver
        wing_c = vsm_wing
    else
        sh = shadow_ref[]
        if sh === nothing || eltype(sh[1]._va) !== T
            shadow_ref[] = VortexStepMethod.make_dual_shadow(
                solver, body_aero, T)
            sh = shadow_ref[]
        end
        body_aero_c, solver_c = sh
        wing_c = body_aero_c.wings[1]
    end

    α = alpha_0 + d[1]
    β = beta_0  + d[2]
    va_mag = va_mag_0 * (1 + d[3])
    cα, sα = cos(α), sin(α)
    cβ, sβ = cos(β), sin(β)
    va_b = MVector{3, T}(va_mag * cα * cβ,
                          va_mag * sβ,
                          va_mag * sα * cβ)
    ω = MVector{3, T}(omega_b_0[1] + d[4],
                      omega_b_0[2] + d[5],
                      omega_b_0[3] + d[6])

    theta = Vector{T}(undef, n_unrefined)
    @inbounds for i in 1:n_unrefined
        theta[i] = theta_0[i]
    end
    @inbounds for (gi, gidx) in enumerate(wing_grp_idxs)
        for ui in groups[gidx].unrefined_section_idxs
            theta[ui] += d[6 + gi]
        end
    end

    if n_unrefined > 0
        VortexStepMethod.unrefined_deform!(
            wing_c, theta; smooth=false)
        VortexStepMethod.reinit!(
            body_aero_c; init_aero=false)
    end
    set_va!(body_aero_c, va_b, ω)
    VortexStepMethod.solve!(
        solver_c, body_aero_c; moment_frac, log=false)

    cf = solver_c.sol.force_coeffs
    cm_body = solver_c.sol.moment_coeffs
    cm_unr = solver_c.sol.cm_unrefined_dist

    drag_dir = va_b ./ va_mag
    side_dir = SVector(zero(T), one(T), zero(T))
    lift_dir = normalize(cross(drag_dir, side_dir))

    n_g = length(wing_grp_idxs)
    x = Vector{T}(undef, 6 + n_g)
    x[1] = dot(cf, lift_dir)
    x[2] = dot(cf, drag_dir)
    x[3] = dot(cf, side_dir)
    x[4] = cm_body[1]
    x[5] = cm_body[2]
    x[6] = cm_body[3]
    @inbounds for (gi, gidx) in enumerate(wing_grp_idxs)
        s = zero(T)
        for ui in groups[gidx].unrefined_section_idxs
            s += cm_unr[ui]
        end
        x[6 + gi] = s
    end
    return x
end

# ─── Baseline + Jacobian via ForwardDiff ─────────────

n_g = length(wing.group_idxs)
n_inputs = 6 + n_g

x0 = coeffs_at_perturbation(zeros(n_inputs), wing.group_idxs)

@info "Baseline" CL=x0[1] CD=x0[2] CS=x0[3] CM=x0[4:6]

@info "Computing ForwardDiff Jacobian at operating point …"
J = ForwardDiff.jacobian(
    d -> coeffs_at_perturbation(d, wing.group_idxs),
    zeros(n_inputs))

# ─── Plot configuration ──────────────────────────────

group_names = [
    groups[gidx].name for gidx in wing.group_idxs]

input_labels = [
    "α", "β", "|va|",
    "ω₁", "ω₂", "ω₃",
    ["θ_$n" for n in group_names]...
]
output_labels = [
    "CL", "CD", "CS",
    "CM₁", "CM₂", "CM₃",
    ["cm_$n" for n in group_names]...
]

n_outputs = length(output_labels)

sweep_range = range(-0.1, 0.1; length=21)

# ─── Plot grid: rows = outputs, cols = inputs ────────

fig = Figure(size=(180 * n_inputs + 80,
                   90 * n_outputs + 80))

Label(fig[0, 1:n_inputs],
    "VSM linearisation (ForwardDiff tangents) " *
    "around operating point";
    fontsize=18, font=:bold,
    tellwidth=false)

for ri in 1:n_outputs
    for ci in 1:n_inputs
        ax = Axis(fig[ri, ci];
            xlabel=ri == n_outputs ? input_labels[ci] : "",
            ylabel=ci == 1 ? output_labels[ri] : "",
            xticklabelsvisible=ri == n_outputs,
            yticklabelsvisible=ci == 1,
            xticksvisible=ri == n_outputs,
            yticksvisible=ci == 1)

        # ForwardDiff tangent: y = x0[ri] + J[ri,ci] * d
        tangent = x0[ri] .+ J[ri, ci] .* sweep_range
        lines!(ax, sweep_range, tangent;
            color=:darkorange, linewidth=1.5)
        scatter!(ax, [0.0], [x0[ri]];
            color=:crimson, markersize=6)
    end
end

# Tighten layout
colgap!(fig.layout, 6)
rowgap!(fig.layout, 4)

display(fig)
