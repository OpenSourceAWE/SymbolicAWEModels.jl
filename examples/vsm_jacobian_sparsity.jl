# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

"""
Sweep every VSM input from -0.1 to +0.1 around the operating
point and plot all output coefficients with GLMakie.

Inputs:
    α, β   [rad]   apparent-flow direction perturbations
    |va|   [—]     fractional wind-speed perturbation
    ω₁..₃  [rad/s] body angular rate perturbations
    θ_g    [rad]   per-group twist perturbations

Outputs (wind-axis coefficients):
    CL, CD, CS, CM₁..₃, cm_g per group

Solver type and `use_gamma_prev` come from `vsm_settings.yaml`.
"""

using Pkg
if Base.active_project() != joinpath(@__DIR__, "Project.toml")
    Pkg.activate(@__DIR__)
end

using GLMakie
using LinearAlgebra
using Printf
using KiteUtils: init!
using SymbolicAWEModels, VortexStepMethod

# ─── Model setup ─────────────────────────────────────
MODEL_NAME = "2plate_kite"

pkg_root = dirname(@__DIR__)
set_data_path(joinpath(pkg_root, "data", MODEL_NAME))

struc_yaml = joinpath(
    get_data_path(), "quat_struc_geometry.yaml")

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

# ─── Helpers ─────────────────────────────────────────

function va_from_angles(alpha, beta, va_mag)
    return [va_mag * cos(alpha) * cos(beta),
            va_mag * sin(beta),
            va_mag * sin(alpha) * cos(beta)]
end

last_theta = fill(NaN, 0)

function solve_coeffs!(
    solver, body_aero, vsm_wing,
    va_b, omega_b, theta;
    moment_frac, groups, group_idxs
)
    if length(last_theta) != length(theta) ||
            !all(theta .== last_theta)
        VortexStepMethod.unrefined_deform!(
            vsm_wing, theta; smooth=false)
        reinit!(body_aero; init_aero=false)
        global last_theta = copy(theta)
    end
    set_va!(body_aero, va_b, omega_b)
    VortexStepMethod.solve!(
        solver, body_aero;
        moment_frac, log=false)

    cf = solver.sol.force_coeffs
    cm_body = solver.sol.moment_coeffs
    cm_unr = solver.sol.cm_unrefined_dist

    drag_dir = va_b / norm(va_b)
    side_dir = [0.0, 1.0, 0.0]
    lift_dir = normalize(cross(drag_dir, side_dir))

    CL = dot(cf, lift_dir)
    CD = dot(cf, drag_dir)
    CS = dot(cf, side_dir)

    cm_groups = zeros(length(group_idxs))
    for (i, gidx) in enumerate(group_idxs)
        g = groups[gidx]
        cm_groups[i] = sum(
            cm_unr[ui]
            for ui in g.unrefined_section_idxs)
    end

    return [CL, CD, CS,
            cm_body[1], cm_body[2], cm_body[3],
            cm_groups...]
end

# ─── Baseline solve ──────────────────────────────────

x0 = solve_coeffs!(
    solver, body_aero, vsm_wing,
    va_b_0, omega_b_0, theta_0;
    moment_frac, groups,
    group_idxs=wing.group_idxs)

@info "Baseline" CL=x0[1] CD=x0[2] CS=x0[3] CM=x0[4:6]

# ─── Sweep configuration ─────────────────────────────

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

n_inputs = length(input_labels)
n_outputs = length(output_labels)

sweep_range = range(-0.1, 0.1; length=21)
n_sweep = length(sweep_range)

# results[input_idx][sweep_idx, output_idx]
results = [zeros(n_sweep, n_outputs) for _ in 1:n_inputs]

# ─── Run sweeps ──────────────────────────────────────

for (col, label) in enumerate(input_labels)
    @info "Sweeping $label …"
    for (si, d) in enumerate(sweep_range)
        va_b = copy(va_b_0)
        omega = copy(omega_b_0)
        theta = copy(theta_0)

        if col == 1       # alpha
            va_b .= va_from_angles(
                alpha_0 + d, beta_0, va_mag_0)
        elseif col == 2   # beta
            va_b .= va_from_angles(
                alpha_0, beta_0 + d, va_mag_0)
        elseif col == 3   # wind scale (fractional)
            va_b .= va_b_0 .* (1 + d)
        elseif col <= 6   # omega_b component
            omega[col - 3] += d
        else              # group twist
            gi = col - 6
            gidx = wing.group_idxs[gi]
            for ui in groups[gidx].unrefined_section_idxs
                theta[ui] += d
            end
        end

        x = solve_coeffs!(
            solver, body_aero, vsm_wing,
            va_b, omega, theta;
            moment_frac, groups,
            group_idxs=wing.group_idxs)
        results[col][si, :] .= x
    end
end

# Restore baseline state
VortexStepMethod.unrefined_deform!(
    vsm_wing, theta_0; smooth=false)
reinit!(body_aero; init_aero=false)
set_va!(body_aero, va_b_0, omega_b_0)

# ─── Plot grid: rows = outputs, cols = inputs ────────

fig = Figure(size=(180 * n_inputs + 80,
                   90 * n_outputs + 80))

Label(fig[0, 1:n_inputs],
    "VSM input sweep [-0.1, 0.1] around operating point";
    fontsize=18, font=:bold,
    tellwidth=false)

axes_grid = Matrix{Axis}(undef, n_outputs, n_inputs)

for ri in 1:n_outputs
    for ci in 1:n_inputs
        ax = Axis(fig[ri, ci];
            xlabel=ri == n_outputs ? input_labels[ci] : "",
            ylabel=ci == 1 ? output_labels[ri] : "",
            xticklabelsvisible=ri == n_outputs,
            yticklabelsvisible=ci == 1,
            xticksvisible=ri == n_outputs,
            yticksvisible=ci == 1)
        axes_grid[ri, ci] = ax

        ys = results[ci][:, ri]
        lines!(ax, sweep_range, ys; color=:steelblue,
            linewidth=1.5)
        scatter!(ax, [0.0], [x0[ri]];
            color=:crimson, markersize=6)
    end
end

# Tighten layout
colgap!(fig.layout, 6)
rowgap!(fig.layout, 4)

display(fig)
