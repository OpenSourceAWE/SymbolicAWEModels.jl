# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

"""
Rigidity laws of pressurised fabric tubes, supplying the `EA`, `GA`, `GJ`, `EIy`
and `EIz` of a [`TimoshenkoJoint`](@ref) that models an inflated leading edge or
strut.

Two independent sources are available. The Breukels (2011) empirical
correlations fit a 1 m cantilever tip-force curve and cover the linear regime up
to collapse. The Comer-Levy (1963) wrinkled-section theory is analytical, needs
the fabric membrane stiffness `E·t` on top of radius and pressure, and stays
valid past collapse. Both are exposed as a callable [`TubeRigidityLaw`](@ref)
that a joint evaluates at its current curvature.
"""

# Breukels (2011) empirical correlation constants for inflated fabric tubes,
# as used by ASKITE/kite_fem. Radius in m, pressure in bar.
"""Breukels bending correlation constants `C1`-`C8`."""
const BREUKELS_BENDING = (
    C1=6582.82, C2=-272.43, C3=40852.38, C4=14.31,
    C5=271865251.42, C6=215.93, C7=14021.79, C8=-589.05,
)

"""Breukels collapse-deflection correlation constants `C9`-`C12`."""
const BREUKELS_COLLAPSE = (C9=322.55, C10=0.0239, C11=5.3833, C12=0.0461)

"""Breukels torsion correlation constants `C13`-`C19`."""
const BREUKELS_TORSION = (
    C13=1467.0, C14=40.908, C15=-191.8, C16=47.406,
    C17=-17703.0, C18=358.05, C19=0.0918,
)

"""Shear correction factor of a thin-walled inflated tube (kite_fem)."""
const TUBE_SHEAR_COEFF = 8 / 9

"""Poisson ratio of the isotropic tube fabric, used for the membrane shear modulus."""
const TUBE_POISSON_RATIO = 0.3

"""
    TubeRigidityLaw

Callable effective rigidity of an inflated tube for one [`TimoshenkoJoint`](@ref)
mode, evaluated at that mode's curvature [1/m]. `mode == :torsion` returns
`GJ(κ) = c1·atan(c2·κ)/κ` (Breukels). `mode == :bending` returns `EI(κ)`: the
linear rigidity below the knee curvature and the smooth power-law approach to the
collapse moment above it. One struct covers both modes because all callable
rigidities of a joint must share a type.

$(TYPEDFIELDS)
"""
struct TubeRigidityLaw
    "Which rigidity this law returns: `:bending` or `:torsion`."
    mode::Symbol
    "Linear bending rigidity below the knee [N·m²]."
    EI0::Float64
    "Moment at the knee, where softening starts [N·m]."
    moment_knee::Float64
    "Asymptotic collapse moment [N·m]."
    moment_collapse::Float64
    "Curvature at the knee [1/m]."
    curvature_knee::Float64
    "Exponent of the post-knee power-law approach to `moment_collapse`."
    exponent::Float64
    "Torsion coefficient `c1` [N·m]."
    c1::Float64
    "Torsion coefficient `c2` [m]."
    c2::Float64
end

function (law::TubeRigidityLaw)(curvature)
    κ = abs(curvature)
    if law.mode === :torsion
        κ < 1e-6 && return law.c1 * law.c2
        return law.c1 * atan(law.c2 * κ) / κ
    end
    κ <= law.curvature_knee && return law.EI0
    moment = law.moment_collapse -
        (law.moment_collapse - law.moment_knee) *
        (law.curvature_knee / κ)^law.exponent
    return moment / κ
end

"""
    breukels_tip_force(deflection, radius, pressure) -> Float64

Tip force [N] of a 1 m inflated-tube cantilever at normalized tip `deflection`,
per the Breukels correlation (`radius` [m], `pressure` [bar]).
"""
function breukels_tip_force(deflection, radius, pressure)
    (; C1, C2, C3, C4, C5, C6, C7, C8) = BREUKELS_BENDING
    denom = (C1 * radius + C2) * pressure^2 + (C3 * radius^3 + C4)
    numer = (C5 * radius^5 + C6) * pressure + (C7 * radius + C8)
    return denom * (1 - exp(-(numer / denom) * deflection))
end

"""
    breukels_collapse_deflection(radius, pressure) -> Float64

Normalized tip deflection at which the Breukels 1 m cantilever of `radius` [m] at
`pressure` [bar] collapses.
"""
function breukels_collapse_deflection(radius, pressure)
    (; C9, C10, C11, C12) = BREUKELS_COLLAPSE
    return (C9 * radius^4 + C10) * pressure + C11 * radius^2 + C12
end

"""
    tube_torsion_law(radius, pressure) -> TubeRigidityLaw

Breukels torsional rigidity law `GJ(κ)` of an inflated tube (`radius` [m],
`pressure` [bar]).
"""
function tube_torsion_law(radius, pressure)
    (; C13, C14, C15, C16, C17, C18, C19) = BREUKELS_TORSION
    c1 = (C13 * radius + C14) * pressure + (C15 * radius + C16)
    c2 = C17 * radius^4 * log(pressure) + C18 * radius^3 + C19
    return TubeRigidityLaw(:torsion, 0.0, 0.0, 0.0, 0.0, 0.0, c1, c2)
end

"""
    tube_linear_rigidities(radius, pressure) -> (EA, GA, EI0, GJ0)

Small-deformation rigidities of an inflated tube from the Breukels laws:
`GJ0 = c1·c2`, shear modulus from `GJ0` and the section (`J = πr⁴/2`, `A = πr²`),
`EI0` from the initial slope of the tip-force curve with the shear share of the
deflection removed, and `EA = EI0·A/I`.
"""
function tube_linear_rigidities(radius, pressure)
    area = π * radius^2
    inertia = π * radius^4 / 4
    torsion_constant = 2 * inertia
    torsion = tube_torsion_law(radius, pressure)
    GJ0 = torsion.c1 * torsion.c2
    GA = GJ0 / torsion_constant * area
    (; C5, C6, C7, C8) = BREUKELS_BENDING
    numer = (C5 * radius^5 + C6) * pressure + (C7 * radius + C8)
    shear_frac = numer / (TUBE_SHEAR_COEFF * GA)
    shear_frac < 1 || error(
        "Breukels shear share ≥ 1 for radius $radius, pressure $pressure")
    EI0 = numer / (3 * (1 - shear_frac))
    EA = EI0 * area / inertia
    return EA, GA, EI0, GJ0
end

"""
    tube_bending_law(radius, pressure; n_samples=120) -> TubeRigidityLaw

Breukels bending rigidity law `EI(κ)` of an inflated tube (`radius` [m],
`pressure` [bar]): moment-curvature samples from the tip-force curve (curvature
`κ = 3·(δ - δ_shear)` of the 1 m cantilever), knee seeded at half the collapse
moment, and the post-knee exponent from a least-squares fit of
`M_c - M = (M_c - M_knee)·(κ_knee/κ)^a` on the sampled tail.
"""
function tube_bending_law(radius, pressure; n_samples=120)
    _, GA, EI0, _ = tube_linear_rigidities(radius, pressure)
    collapse = breukels_collapse_deflection(radius, pressure)
    moment_collapse = breukels_tip_force(collapse, radius, pressure)
    moment_knee = 0.5 * moment_collapse
    curvature_knee = moment_knee / EI0
    deflections = range(0.0, collapse; length=n_samples)[2:end]
    x = Float64[]
    y = Float64[]
    for δ in deflections
        moment = breukels_tip_force(δ, radius, pressure)
        moment_knee < moment < 0.999 * moment_collapse || continue
        δ_shear = moment / (TUBE_SHEAR_COEFF * GA)
        κ = 3 * (δ - δ_shear)
        push!(x, log(curvature_knee / κ))
        push!(y, log((moment_collapse - moment) /
                     (moment_collapse - moment_knee)))
    end
    exponent = isempty(x) ? 1.0 : x \ y
    return TubeRigidityLaw(:bending, EI0, moment_knee, moment_collapse,
        curvature_knee, exponent, 0.0, 0.0)
end

# ---- Comer-Levy analytical inflated-cylinder bending (post-collapse capable) ----
# Isotropic wrinkled-section theory (Comer & Levy 1963): closed-form moment-
# curvature of a pressurised fabric tube from the linear regime, through wrinkling
# onset, to the collapse plateau — and valid past collapse, which the empirical
# Breukels tip-force correlation is not. Needs the fabric membrane stiffness E·t
# [N/m] on top of (radius, pressure); `breukels_membrane_stiffness` derives that
# E·t from Breukels so no external coupon value is required.

"""
    comer_levy_bending_stiffness(radius, membrane_stiffness) -> Float64

Pre-wrinkling bending stiffness `EI = E·t·π·r³` [N·m²] (Fichter) of a thin-walled
inflated tube of `radius` [m] and membrane stiffness `membrane_stiffness` (`E·t`,
[N/m]). Pressure enters the wrinkling/collapse limits, not this small-deflection
slope.
"""
comer_levy_bending_stiffness(radius, membrane_stiffness) =
    membrane_stiffness * π * radius^3

"""
    comer_levy_wrinkling_moment(radius, pressure) -> Float64

Wrinkling-onset moment `M_w = p·π·r³/2` [N·m] (`radius` [m], `pressure` [Pa]): the
compressed fibre's pressure pre-stress is cancelled and the fabric first goes
slack (Comer-Levy slack angle θ₀ = 0).
"""
comer_levy_wrinkling_moment(radius, pressure) = pressure * π * radius^3 / 2

"""
    comer_levy_collapse_moment(radius, pressure) -> Float64

Collapse moment `M_c = p·π·r³` [N·m] (`radius` [m], `pressure` [Pa]), the θ₀ → π
asymptote; equals `2·M_w`.
"""
comer_levy_collapse_moment(radius, pressure) = pressure * π * radius^3

"""Slack-section geometry factor `h(θ) = sin θ + (π − θ) cos θ`, monotone π → 0."""
comer_levy_h(θ) = sin(θ) + (π - θ) * cos(θ)

"""Slack-section geometry factor `n(θ) = 2π − 2θ + sin 2θ`."""
comer_levy_n(θ) = 2π - 2θ + sin(2θ)

"""
    comer_levy_point(θ, radius, pressure, membrane_stiffness) -> (κ, M)

Curvature [1/m] and moment [N·m] of the Comer-Levy wrinkled section at slack angle
`θ ∈ (0, π)`: `κ = π·p / (2·E·t·h(θ))` (radius-independent) and
`M = (p·r³/4)·π·n(θ)/h(θ)`.
"""
function comer_levy_point(θ, radius, pressure, membrane_stiffness)
    h = comer_levy_h(θ)
    return π * pressure / (2 * membrane_stiffness * h),
        pressure * radius^3 / 4 * π * comer_levy_n(θ) / h
end

"""
    comer_levy_sample_curve(radius, pressure, membrane_stiffness; n=60, frac_max=0.97)
        -> (κ, M)

Synthetic moment-curvature table [1/m], [N·m] for a tube of `radius`: linear
samples below wrinkling plus exact Comer-Levy wrinkled-section samples up to
`frac_max` of the collapse moment.
"""
function comer_levy_sample_curve(radius, pressure, membrane_stiffness;
        n=60, frac_max=0.97)
    EI = comer_levy_bending_stiffness(radius, membrane_stiffness)
    M_w = comer_levy_wrinkling_moment(radius, pressure)
    M_c = comer_levy_collapse_moment(radius, pressure)
    κ_lin = collect(range(0.0, M_w / EI; length=max(2, n ÷ 3)))
    M_lin = EI .* κ_lin
    M_cap = frac_max * M_c
    points = [comer_levy_point(θ, radius, pressure, membrane_stiffness)
              for θ in range(1.0e-3, π - 1.0e-3; length=4n)]
    keep = [M <= M_cap for (_, M) in points]
    κ_post = [κ for (κ, _) in points[keep]]
    M_post = [M for (_, M) in points[keep]]
    return vcat(κ_lin, κ_post), vcat(M_lin, M_post)
end

"""
    fit_bending_law(κ, M, moment_knee, moment_collapse)
        -> (EI0, curvature_knee, exponent)

Least-squares fit of the smooth [`TubeRigidityLaw`](@ref) bending law to sampled
`(κ, M)` given the wrinkling/collapse anchors: `EI0` from the linear region
(`M < moment_knee`), `curvature_knee = moment_knee/EI0`, and the post-wrinkling
`exponent` from the log-linearised tail
`M_c − M = (M_c − M_w)(κ_w/κ)^exponent`. Backslash, no optimisation dependency.
"""
function fit_bending_law(κ, M, moment_knee, moment_collapse)
    linear = M .< moment_knee
    EI0 = κ[linear] \ M[linear]
    curvature_knee = moment_knee / EI0
    post = (M .> moment_knee) .& (M .< 0.999 * moment_collapse)
    x = log.(curvature_knee ./ κ[post])
    y = log.((moment_collapse .- M[post]) ./ (moment_collapse - moment_knee))
    return EI0, curvature_knee, (x \ y)
end

"""
    comer_levy_bending_law(radius, pressure, membrane_stiffness; n_samples=60)
        -> TubeRigidityLaw

Curvature-softening bending rigidity `EI(κ)` of an inflated tube from the
analytical Comer-Levy wrinkled section (`radius` [m], `pressure` [Pa],
`membrane_stiffness` `E·t` [N/m]): linear `EI0 = E·t·π·r³` below the wrinkling
knee, then the smooth power-law approach to the collapse moment `M_c = p·π·r³`
fitted from the closed-form section curve. The returned
[`TubeRigidityLaw`](@ref) `:bending` callable is the drop-in replacement for
[`tube_bending_law`](@ref) and, unlike it, stays valid past collapse (moment
saturates at `M_c`, rigidity → 0).
"""
function comer_levy_bending_law(radius, pressure, membrane_stiffness; n_samples=60)
    κ, M = comer_levy_sample_curve(radius, pressure, membrane_stiffness; n=n_samples)
    moment_knee = comer_levy_wrinkling_moment(radius, pressure)
    moment_collapse = comer_levy_collapse_moment(radius, pressure)
    EI0, curvature_knee, exponent = fit_bending_law(κ, M, moment_knee, moment_collapse)
    return TubeRigidityLaw(:bending, EI0, moment_knee, moment_collapse,
        curvature_knee, exponent, 0.0, 0.0)
end

"""
    breukels_membrane_stiffness(radius, pressure_bar) -> Float64

Membrane stiffness `E·t` [N/m] implied by the Breukels linear bending rigidity,
`E·t = EI0_breukels / (π·r³)` (inverting `EI0 = E·t·π·r³`). This sources the
Comer-Levy fabric input from the same empirical law the model already uses, so the
linear regime is unchanged and only the pressure-driven collapse branch is added.
`radius` [m], `pressure_bar` [bar]. Trustworthy only in Breukels' fitted range
(`radius ≲ 0.15 m`).
"""
function breukels_membrane_stiffness(radius, pressure_bar)
    _, _, EI0, _ = tube_linear_rigidities(radius, pressure_bar)
    return EI0 / (π * radius^3)
end

"""
    membrane_linear_rigidities(radius, membrane_stiffness) -> (EA, GA, EI0)

Thin-walled inflated-tube axial, shear and bending rigidities from the membrane
stiffness `membrane_stiffness` (`E·t` [N/m]) and `radius` [m], all with one
consistent provenance: `EA = E·t·2πr`, `GA = G·t·2πr` with
`G·t = E·t / (2(1 + ν))` ([`TUBE_POISSON_RATIO`](@ref)), and `EI0 = E·t·π·r³`. The
Comer-Levy counterpart of [`tube_linear_rigidities`](@ref) (Breukels); torsion
`GJ` is not covered here and stays on [`tube_torsion_law`](@ref).
"""
function membrane_linear_rigidities(radius, membrane_stiffness)
    circumference = 2π * radius
    shear_stiffness = membrane_stiffness / (2 * (1 + TUBE_POISSON_RATIO))
    EA = membrane_stiffness * circumference
    GA = shear_stiffness * circumference
    EI0 = comer_levy_bending_stiffness(radius, membrane_stiffness)
    return EA, GA, EI0
end

"""
    tube_mass(radius, len; areal_density) -> Float64

Fabric mass [kg] of an inflated tube of `radius` and length `len` [m].
"""
tube_mass(radius, len; areal_density) = 2π * radius * len * areal_density
