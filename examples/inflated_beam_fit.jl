# Copyright (c) 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Fit the bending stiffness of an inflatable tube (kite leading edge) WITHOUT
# real bend tests. The isotropic Comer-Levy wrinkled-section theory acts as a
# virtual experiment: it gives the moment-curvature curve M(κ) of the pressurised
# tube in closed form (parametrised by the slack-arc angle θ₀), from the linear
# regime through wrinkling onset to collapse. We sample that curve as synthetic
# "measurements", fit a smooth per-joint bending law to it, and build a
# Body + ElasticJoint chain. It is validated as a cantilever under a ramped
# downward tip force, reproducing the P-vs-tip-deflection curve up to collapse.
#
# Isotropic by choice: Comer-Levy assumes one modulus E. A woven (orthotropic)
# tube has no closed-form post-wrinkling law — we approximate it by using the
# axial membrane stiffness E·t here, which is the standard pragmatic move.
#
# Caveat: this does not remove measurement, it RELOCATES it — from "bend a whole
# tube" to "fabric coupon / datasheet" (E·t, t, internal pressure, radius
# profile), which remain required inputs.
#
# This script also overlays the Breukels (2011) empirical fits used by the TU
# Delft V3 kite_fem/ASKITE pipeline, to compare against the Comer-Levy curve at
# real V3 leading-edge/strut radii. Breukels needs only (r, p) — no E·t or wall
# thickness — but its constants are opaque and only valid within its fitted range.
#
# Reference:
#   Comer, R.L. & Levy, S. (1963). Deflections of an inflated circular-cylindrical
#     cantilever beam. AIAA Journal 1(7), 1652-1655.
#     https://doi.org/10.2514/3.1873
#     Gives everything used here: linear EI = E·t·π·r³, wrinkling onset
#     M_w = pπr³/2, collapse M_c = pπD³/8 = pπr³, and the post-wrinkling
#     M(θ₀), κ(θ₀) parametrised by the slack-arc angle θ₀.

using Pkg
Pkg.activate(@__DIR__)
using SymbolicAWEModels
using SymbolicAWEModels: Point
using KiteUtils
using GLMakie

# ----- tube geometry and fabric (the only required measured inputs) -----
pressure_bar = 0.3               # internal gauge pressure [bar]
pressure = pressure_bar * 1.0e5  # [Pa] (1 bar = 1e5 Pa), used by all formulas
# Radius distribution along the span as (fractional position s∈[0,1], radius [m])
# control points (s=0 root, s=1 tip), piecewise-linearly interpolated. Add points
# for a taper or a local pinch; here the tube is uniform (d = 130 mm → r = 65 mm).
radius_profile = [(0.0, 0.065), (1.0, 0.065)]
beam_length = 1.0                # [m] total tube length
n_hf = 30                        # links in the high-fidelity (radius-resolved) beam
membrane_stiffness = 2.0e5   # E·t [N/m]

"""
    station_radius(s) -> Float64

Tube radius [m] at fractional span `s ∈ [0,1]` (0 = root, 1 = tip), piecewise-
linearly interpolated from `radius_profile`'s control points (assumed sorted by
`s`). Clamps to the end values outside the control range.
"""
function station_radius(s)
    s <= radius_profile[1][1] && return radius_profile[1][2]
    s >= radius_profile[end][1] && return radius_profile[end][2]
    for k in 2:length(radius_profile)
        s_lo, r_lo = radius_profile[k - 1]
        s_hi, r_hi = radius_profile[k]
        s <= s_hi && return r_lo + (r_hi - r_lo) * (s - s_lo) / (s_hi - s_lo)
    end
    return radius_profile[end][2]
end

# ----- analytical inflated-beam model (isotropic Comer-Levy wrinkled section) -----

"""
    bending_stiffness(radius) -> Float64

Pre-wrinkling bending stiffness `EI = E·t·π·r³` [N·m²] (Fichter). Pressure enters
the wrinkling/collapse limits, not this small-deflection slope.
"""
bending_stiffness(radius) = membrane_stiffness * π * radius^3

"""
    wrinkling_moment(radius) -> Float64

Wrinkling-onset moment `M_w = p·π·r³/2` [N·m]: the compressed fibre's pressure
pre-stress is cancelled and the fabric first goes slack (Comer-Levy, θ₀ = 0).
"""
wrinkling_moment(radius) = pressure * π * radius^3 / 2

"""
    collapse_moment(radius) -> Float64

Collapse moment `M_c = p·π·D³/8 = p·π·r³` [N·m] (Comer-Levy, θ₀ → π). Equals 2·M_w.
"""
collapse_moment(radius) = pressure * π * radius^3

# Comer-Levy wrinkled cross-section, parametrised by the slack-arc angle θ₀∈[0,π].
# `comer_levy_h` is monotone decreasing π → 0; moment and curvature are M(θ₀), κ(θ₀).
comer_levy_h(θ) = sin(θ) + (π - θ) * cos(θ)
comer_levy_n(θ) = 2π - 2θ + sin(2θ)

"""
    comer_levy_point(θ, radius) -> (κ, M)

Curvature [1/m] and moment [N·m] of the wrinkled section at slack angle `θ`∈(0,π):
`κ = π·p / (2·E·t·h)` (radius-independent), `M = (p·r³/4)·π·n/h`.
"""
function comer_levy_point(θ, radius)
    h = comer_levy_h(θ)
    return π * pressure / (2 * membrane_stiffness * h),
        pressure * radius^3 / 4 * π * comer_levy_n(θ) / h
end

"""
    curvature_of_moment(moment, radius) -> Float64

Curvature [1/m] sustaining `moment` (exact Comer-Levy; used for the analytical
reference only, NOT the ODE RHS). Linear below wrinkling; bisection on θ₀ above,
where the ratio `n/h` rises monotonically 2 → 4 over θ₀ ∈ [0, π]. Diverges at M_c.
"""
function curvature_of_moment(moment, radius)
    moment <= wrinkling_moment(radius) && return moment / bending_stiffness(radius)
    moment >= collapse_moment(radius) && return Inf
    target = 4 * moment / (π * pressure * radius^3)        # = n/h ∈ (2, 4)
    lo, hi = 0.0, Float64(π)
    for _ in 1:80
        mid = 0.5 * (lo + hi)
        comer_levy_n(mid) / comer_levy_h(mid) < target ? (lo = mid) : (hi = mid)
    end
    return comer_levy_point(0.5 * (lo + hi), radius)[1]
end

# ----- Breukels (2011) empirical inflatable-beam fits (TU Delft kite_fem) -----
# Breukels gives the tip load P [N] vs tip deflection δ [m] of a 1 m cantilever
# tube as an empirical function of radius and pressure IN BAR (not Pa). This is
# the bending model the V3 kite_fem/ASKITE pipeline uses. It needs no fabric E·t
# or wall thickness, only (r, p) — but the constants are opaque and extrapolate
# poorly outside Breukels' fitted diameter/pressure range.

"""
    breukels_tip_force(deflection, radius) -> Float64

Breukels' empirical tip load P [N] for a 1 m inflatable cantilever of `radius`
deflected `deflection` [m] at `pressure_bar`. Saturating exponential: P rises to
the `denom` asymptote as deflection grows (collapse is a separate criterion).
"""
function breukels_tip_force(deflection, radius)
    p = pressure_bar
    denom = (6582.82 * radius - 272.43) * p^2 + (40852.38 * radius^3 + 14.31)
    numer = (2.7186525142e8 * radius^5 + 215.93) * p + (14021.79 * radius - 589.05)
    return denom * (1 - exp(-(numer / denom) * deflection))
end

"""
    breukels_collapse_deflection(radius) -> Float64

Breukels' tip deflection [m] at which a 1 m tube of `radius` wrinkles through
(collapses) at `pressure_bar`.
"""
breukels_collapse_deflection(radius) =
    (322.55 * radius^4 + 0.0239) * pressure_bar + 5.3833 * radius^2 + 0.0461

"""
    breukels_section_curve(radius; n=120) -> (κ, M)

Effective section moment-curvature [1/m], [N·m] implied by Breukels' 1 m
cantilever P(δ), via constant-EI cantilever kinematics: root curvature κ = 3δ/L²
(L = 1 m), root moment M = P(δ)·L. Sampled from δ = 0 to the Breukels collapse
deflection. Approximate (folds the shear share of δ into bending); for comparison
against the Comer-Levy section law only.
"""
function breukels_section_curve(radius; n = 120)
    δ = range(0.0, breukels_collapse_deflection(radius); length = n)
    return 3.0 .* δ, breukels_tip_force.(δ, radius)
end

# ----- generate synthetic measurements and fit the per-joint law -----

"""
    sample_curve(radius; n=60, frac_max=0.97) -> (κ, M)

Synthetic measurement table for a tube of `radius`: linear samples below wrinkling
plus exact Comer-Levy wrinkled-section samples up to `frac_max` of collapse.
"""
function sample_curve(radius; n = 60, frac_max = 0.97)
    κ_w = wrinkling_moment(radius) / bending_stiffness(radius)
    κ_lin = collect(range(0.0, κ_w; length = max(2, n ÷ 3)))
    M_lin = bending_stiffness(radius) .* κ_lin
    M_cap = frac_max * collapse_moment(radius)
    points = [comer_levy_point(θ, radius) for θ in range(1.0e-3, π - 1.0e-3; length = 4n)]
    keep = [M <= M_cap for (_, M) in points]
    κ_post = [κ for (κ, _) in points[keep]]
    M_post = [M for (_, M) in points[keep]]
    return vcat(κ_lin, κ_post), vcat(M_lin, M_post)
end

"""
    fit_law(κ, M, M_w, M_c) -> (EI, κ_w, exponent)

Least-squares fit of the smooth bending law to sampled `(κ, M)`, given the
wrinkling/collapse moment anchors. `EI` from the linear region (`M < M_w`),
`κ_w = M_w/EI`, and the post-wrinkling exponent `a` from the log-linearised tail
`M_c − M = (M_c − M_w)(κ_w/κ)^a` (Comer-Levy decays as κ^(−2/3)). Backslash, no
optimisation dependency. Works for one station or a many-segment aggregate.
"""
function fit_law(κ, M, M_w, M_c)
    linear = M .< M_w
    EI = κ[linear] \ M[linear]
    κ_w = M_w / EI
    post = (M .> M_w) .& (M .< 0.999 * M_c)
    x = log.(κ_w ./ κ[post])
    y = log.((M_c .- M[post]) ./ (M_c - M_w))
    return EI, κ_w, (x \ y)
end

"""
    WrinklingLaw(EI, M_w, M_c, κ_w, exponent, length)

Callable per-joint bending law: maps joint bend angle `Δθ` to curvature
`κ = Δθ/length` and returns the section moment [N·m]. Linear `EI·κ` below the
wrinkling curvature `κ_w`, then the smooth power-law approach to collapse
`M_c − (M_c − M_w)(κ_w/κ)^exponent`. Smooth with no internal solve, so it is
ForwardDiff-safe inside the ODE RHS Jacobian. Passed as `stiffness_bending`.
Mutable so the 2-seg ODE-sim fit can retune `EI`/`κ_w`/`exponent` in place — the
joint reads the law live each evaluation, so no model rebuild is needed.
"""
mutable struct WrinklingLaw
    EI::Float64
    M_w::Float64
    M_c::Float64
    κ_w::Float64
    exponent::Float64
    length::Float64
end

function (law::WrinklingLaw)(joint_angle)
    κ = abs(joint_angle) / law.length
    moment = κ < law.κ_w ? law.EI * κ :
        law.M_c - (law.M_c - law.M_w) * (law.κ_w / κ)^law.exponent
    return sign(joint_angle) * moment
end

"""
    linear_bending_stiffness(law) -> Float64

Small-angle joint bending stiffness `dM/dΔθ` at zero [N·m/rad], per law type.
Used to size near-critical joint damping regardless of the law representation.
"""
linear_bending_stiffness(law::WrinklingLaw) = law.EI / law.length

# Demonstrate the fit round-trip on the root station.
let radius = station_radius(0.0)
    κ, M = sample_curve(radius)
    M_w = wrinkling_moment(radius)
    M_c = collapse_moment(radius)
    EI_fit, κ_w, a_fit = fit_law(κ, M, M_w, M_c)
    law = WrinklingLaw(EI_fit, M_w, M_c, κ_w, a_fit, 1.0)
    rel_rms = sqrt(sum(((law.(κ) .- M) ./ M_c) .^ 2) / length(M))
    println("Comer-Levy station fit:  EI $(round(EI_fit; digits=2)) " *
            "(true $(round(bending_stiffness(radius); digits=2)))  exponent a=" *
            "$(round(a_fit; digits=3))  rel-RMS $(round(rel_rms; digits=4))")
end

# ----- compare Comer-Levy vs Breukels M(κ) at representative V3 radii -----
# V3 leading-edge tube radii span ~0.056 m (tip) to ~0.101 m (centre); struts are
# smaller (~0.04-0.07 m). The centre LE is the radius most likely OUTSIDE
# Breukels' fitted range, so watch the divergence there.
v3_stations = [("strut  r50 ", 0.050), ("LE tip r56 ", 0.056),
               ("LE mid r75 ", 0.075), ("LE ctr r101", 0.101)]

println("\n station       EI_CL    EI_Brk   ratio    Mc_CL   Mc_Brk   ratio")
fig_cmp = Figure(size = (900, 650))
for (k, (label, radius)) in enumerate(v3_stations)
    row, col = fldmod1(k, 2)
    cmp_ax = Axis(fig_cmp[row, col]; title = label,
        xlabel = "curvature κ [1/m]", ylabel = "moment M [N·m]")
    κ_cl, M_cl = sample_curve(radius)
    κ_br, M_br = breukels_section_curve(radius)
    lines!(cmp_ax, κ_cl, M_cl; label = "Comer-Levy")
    lines!(cmp_ax, κ_br, M_br; linestyle = :dash, label = "Breukels")
    axislegend(cmp_ax; position = :rb)

    EI_cl = bending_stiffness(radius)
    EI_br = (M_br[2] - M_br[1]) / (κ_br[2] - κ_br[1])
    Mc_cl, Mc_br = collapse_moment(radius), M_br[end]
    println(rpad(label, 12),
            lpad(round(EI_cl; digits = 2), 8), "  ",
            lpad(round(EI_br; digits = 2), 7), "  ",
            lpad(round(EI_br / EI_cl; digits = 2), 5), "  ",
            lpad(round(Mc_cl; digits = 2), 7), "  ",
            lpad(round(Mc_br; digits = 2), 7), "  ",
            lpad(round(Mc_br / Mc_cl; digits = 2), 6))
end
display(fig_cmp)

# ----- settings: reuse the committed beam data dir, no gravity -----
set_data_path(joinpath(dirname(@__DIR__), "data", "beam"))
set = Settings("system.yaml")
set.g_earth = 0.0                       # pure static bend test, no sag
set.physical_model = "inflated_beam_fit"

# Non-bending DOF stay stiff/linear so the tip pull bends only the joints.
stiffness_axial = 1.0e6
stiffness_shear = 1.0e6
stiffness_torsion = 1.0e5
seg_mass = 0.5                            # transient-only (g=0); not fitted

"""
    critical_damping(stiffness, inertia; ratio=1.0) -> Float64

Critical (ζ=`ratio`) damping `2·ratio·√(stiffness·inertia)` of a 2nd-order mode:
translational [N·s/m] from [N/m]·[kg], or rotational [N·m·s/rad] from
[N·m/rad]·[kg·m²]. Lets each joint damp itself near-critically (fast settle, no
creep) from its own stiffness and inertia, instead of a guessed constant.
"""
critical_damping(stiffness, inertia; ratio = 1.0) =
    2.0 * ratio * sqrt(stiffness * inertia)

# Tip-load winch+tether. A CascadedLengthWinch reels a vertical tether at a low
# v_max: tether length is the independent variable, so the test is displacement-
# controlled — it traces the whole force-deflection curve through the collapse
# limit point with no force-control runaway (see the collapse discussion).
winch_depth      = 1.0                    # ground point this far below the tip [m]
sweep_v_max      = 0.01                   # reel-in speed [m/s] (quasi-static)
tether_stiffness = 5.0e5                  # unit stiffness [N]
tether_damping   = 2.0e3                  # unit damping  [N·s]
tether_diameter  = 0.01                   # [m]

"""
    station_law(radius, length, backend) -> WrinklingLaw

Per-station bending law fitted at `radius` for a joint of `length`, from either
virtual experiment: `backend = :comer` (Comer-Levy wrinkled section) or
`:breukels` (Breukels empirical fit). Breukels has no sharp wrinkling knee, so
its linear/post split is seeded at half the collapse moment.
"""
function station_law(radius, length, backend)
    if backend === :comer
        κ, M = sample_curve(radius)
        M_w, M_c = wrinkling_moment(radius), collapse_moment(radius)
    else
        κ, M = breukels_section_curve(radius)
        M_c = M[end]
        M_w = 0.5 * M_c
    end
    EI, κ_w, exponent = fit_law(κ, M, M_w, M_c)
    return WrinklingLaw(EI, M_w, M_c, κ_w, exponent, length)
end

"""
    section_collapse_moment(radius, backend) -> Float64

Section collapse moment [N·m] for the chosen backend (sets the sweep force cap).
"""
section_collapse_moment(radius, backend) = backend === :comer ?
    collapse_moment(radius) : breukels_section_curve(radius)[2][end]

"""
    build_beam(name, n_seg, joint_law) -> SymbolicAWEModel

Build and init an `n_seg`-link Body chain along `beam_length` (root fixed,
per-station radius for inertia), joint `i` carrying `joint_law(i)`. A point is
anchored to the tip body and pulled straight down by a CascadedLengthWinch tether
whose ground anchor sits `winch_depth` below the tip — the displacement actuator.
"""
function build_beam(name, n_seg, joint_law)
    seg_len = beam_length / n_seg
    bodies = Body[]
    for i in 1:n_seg
        radius = station_radius((i - 0.5) / n_seg)
        inertia = [0.5 * seg_mass * radius^2,
                   seg_mass * seg_len^2 / 12, seg_mass * seg_len^2 / 12]
        push!(bodies, Body(Symbol("seg_$i"); mass = seg_mass,
            inertia_principal = inertia, pos = [(i - 0.5) * seg_len, 0.0, 0.0],
            type = i == 1 ? STATIC : DYNAMIC))
    end
    # Near-critical joint damping from each joint's stiffness and the inertia of
    # everything outboard of it about the hinge — that is the bending mode the
    # joint must damp (the root joint swings the whole beam, not one segment).
    damping_trans = critical_damping(stiffness_axial, seg_mass)
    function make_joint(i)
        law = joint_law(i)
        joint_x = i * seg_len
        inertia_out = sum(seg_mass * seg_len^2 / 12 +
                          seg_mass * ((j - 0.5) * seg_len - joint_x)^2
                          for j in (i + 1):n_seg)
        damping_rot = critical_damping(linear_bending_stiffness(law), inertia_out)
        ElasticJoint(Symbol("joint_$i"), Symbol("seg_$i"), Symbol("seg_$(i + 1)");
            anchor_a = [seg_len / 2, 0.0, 0.0], anchor_b = [-seg_len / 2, 0.0, 0.0],
            stiffness_axial, stiffness_shear, stiffness_torsion,
            stiffness_bending = law, damping_trans, damping_rot)
    end
    joints = [make_joint(i) for i in 1:(n_seg - 1)]
    tip_body = Symbol("seg_$n_seg")
    points = [
        Point(:ground, [beam_length, 0.0, -winch_depth], STATIC),
        Point(:tip_anchor, [beam_length, 0.0, 0.0], BODY_STATIC;
              body = tip_body, anchor_b = [seg_len / 2, 0.0, 0.0]),
    ]
    segments = [Segment(:tether_seg, :ground, :tip_anchor,
        tether_stiffness, tether_damping, tether_diameter; l0 = winch_depth)]
    tethers = [Tether(:tether, [:tether_seg], winch_depth)]
    winches = [Winch(:winch, set, [:tether]; winch_point = :ground,
        model = CascadedLengthWinch(v_max = sweep_v_max,
            position_gain = 20.0, velocity_gain = 20.0))]
    winches[1].inertia_total = 1.0e-4     # tiny rotor → near-instant tracking
    sys = SystemStructure(name, set; points, segments, tethers, winches,
        bodies = bodies, elastic_joints = joints)
    sam = SymbolicAWEModel(set, sys)
    init!(sam; remake = false, prn = true)
    return sam
end

"""
    winch_sweep(sam; force_cap, name, log_to=nothing) -> (δ_mm, force)

Reel the tip tether in continuously (displacement control) until the tether
tension reaches `force_cap`, recording tip deflection [mm] (at `:tip_anchor`, the
loaded point — same span location for any segment count) and tether tension [N]
at every step. Optionally logs for replay.
"""
function winch_sweep(sam; force_cap, name = "sweep", log_to = nothing)
    dt = 0.02
    max_steps = 6000
    tip = sam.sys_struct.points[:tip_anchor]
    tether_seg = sam.sys_struct.segments[:tether_seg]
    z0 = tip.pos_w[3]
    target = sam.sys_struct.tethers[:tether].len - 0.6  # command a long reel-in
    logger = isnothing(log_to) ? nothing : Logger(sam, max_steps + 1)
    sys_state = isnothing(log_to) ? nothing : SysState(sam)
    deflection_mm, tension, sim_time = Float64[], Float64[], 0.0
    for _ in 1:max_steps
        next_step!(sam; set_values = [target], dt, vsm_interval = 0)
        sim_time += dt
        push!(deflection_mm, 1.0e3 * (z0 - tip.pos_w[3]))
        push!(tension, abs(tether_seg.force))
        if !isnothing(logger)
            update_sys_state!(sys_state, sam)
            sys_state.time = sim_time
            log!(logger, sys_state)
        end
        tension[end] > force_cap && break
    end
    !isnothing(logger) && save_log(logger, log_to)
    return deflection_mm, tension
end

"""
    InterpJointLaw(angle, moment)

Callable joint bending law built by spline-interpolating sampled `(Δθ, M)` points
(odd-symmetric in `Δθ`, clamped to the sampled range). Captures an arbitrary
measured/derived `M(Δθ)` with no parametric knee. Use as `stiffness_bending`.
"""
struct InterpJointLaw{S}
    spline::S
    angle_max::Float64
    stiffness0::Float64
    moment_peak::Float64
end
function InterpJointLaw(angle, moment)
    spline = SymbolicAWEModels.CubicSpline(moment, angle)
    stiffness0 = (moment[2] - moment[1]) / (angle[2] - angle[1])
    return InterpJointLaw(spline, angle[end], stiffness0, maximum(moment))
end
function (law::InterpJointLaw)(joint_angle)
    φ = clamp(abs(joint_angle), 0.0, law.angle_max)
    return sign(joint_angle) * law.spline(φ)
end
linear_bending_stiffness(law::InterpJointLaw) = law.stiffness0

"""
    build_joint_law_static(δ_hf_mm, F_hf) -> InterpJointLaw

Build the 2-seg single-joint bending law directly by static inversion of the
high-fidelity sweep — no ODE-in-the-loop, no parametric fit. The 2-seg beam is
two rigid half-segments hinged at mid-span (root half fixed), so each high-fi
point fixes the geometry and the joint moment by static equilibrium:

- joint angle from the loaded-tip deflection (`:tip_anchor` drops `seg_len·sin φ`),
- joint moment `M = F · (lever arm of the tether force about the joint)`.

The `(φ, M)` samples are spline-interpolated into a smooth [`InterpJointLaw`](@ref),
so the roll-off comes straight from the data with no knee artefact. Points beyond
the 2-seg geometric reach (`δ > seg_len`) are dropped.
"""
function build_joint_law_static(δ_hf_mm, F_hf)
    seg_len = beam_length / 2
    joint = [seg_len, 0.0, 0.0]
    ground = [beam_length, 0.0, -winch_depth]
    angle, moment = [0.0], [0.0]
    for (δ_mm, force) in zip(δ_hf_mm, F_hf)
        arg = 1.0e-3 * δ_mm / seg_len            # tip drop = seg_len·sin φ
        arg >= 1.0 && break                      # beyond the 2-seg reach
        φ = asin(arg)
        φ <= angle[end] && continue              # keep spline knots increasing
        anchor = joint .+ seg_len .* [cos(φ), 0.0, -sin(φ)]
        to_ground = ground .- anchor
        dir = to_ground ./ sqrt(sum(abs2, to_ground))
        arm = anchor .- joint
        push!(angle, φ)
        push!(moment, force * abs(arm[3] * dir[1] - arm[1] * dir[3]))
    end
    return InterpJointLaw(angle, moment)
end

# ----- high-fidelity geometry -----
n_joints_hf = n_hf - 1
hf_joint_len = beam_length / n_joints_hf
hf_radii = [station_radius(i / n_hf) for i in 1:n_joints_hf]

# ----- for each backend: high-fi sweep, then static-invert a 2-seg law from it -----
backends = [(:comer, "Comer-Levy"), (:breukels, "Breukels")]
results = Dict{Symbol, NamedTuple}()
for (backend, label) in backends
    min_collapse = minimum(section_collapse_moment.(hf_radii, backend))
    force_cap = 0.9 * min_collapse / beam_length
    println("\n[$label] root M_c=$(round(min_collapse; digits=2)) N·m  " *
            "force cap≈$(round(force_cap; digits=2)) N")

    sam_hf = build_beam("inflated_beam_hf_$backend", n_hf,
        i -> station_law(hf_radii[i], hf_joint_len, backend))
    δ_hf, F_hf = winch_sweep(sam_hf; force_cap,
        log_to = "inflated_beam_hf_$backend")
    println("  high-fi sweep: $(length(δ_hf)) steps, " *
            "δ→$(round(δ_hf[end]; digits=1)) mm, F→$(round(F_hf[end]; digits=2)) N")

    # 2-seg: single mid-span joint, its law built by static inversion of the high-fi
    # sweep (no ODE-in-the-loop fit), then one sweep of the resulting reduced model.
    joint_law2 = build_joint_law_static(δ_hf, F_hf)
    sam2 = build_beam("inflated_beam_2seg_$backend", 2, _ -> joint_law2)
    δ_2, F_2 = winch_sweep(sam2; force_cap = 1.05 * force_cap)
    println("  2-seg law: peak M=$(round(joint_law2.moment_peak; digits=2)) N·m  " *
            "k0=$(round(joint_law2.stiffness0; digits=2)) N·m/rad")

    results[backend] = (; label, δ_hf, F_hf, δ_2, F_2, sam_hf)
end

# ----- plot force-deflection: 30-seg high-fi vs fitted 2-seg, both backends,
#       all in one axis (solid = high-fi, dashed = 2-seg fit, colour = backend) -----
fig = Figure(size = (720, 500))
ax = Axis(fig[1, 1]; title = "Inflated-beam tip load: $n_hf-seg high-fi vs fitted 2-seg",
    xlabel = "tip deflection δ [mm]", ylabel = "tether tension F [N]")
backend_color = Dict(:comer => :steelblue, :breukels => :darkorange)
for (backend, _) in backends
    r = results[backend]
    color = backend_color[backend]
    lines!(ax, r.δ_hf, r.F_hf; color, label = "$(r.label) $n_hf-seg high-fi")
    lines!(ax, r.δ_2, r.F_2; color, linestyle = :dash,
        label = "$(r.label) 2-seg fit")
end
axislegend(ax; position = :rb)
display(fig)

# ----- replay the Comer-Levy high-fi sweep -----
scene = replay(load_log("inflated_beam_hf_comer"),
    results[:comer].sam_hf.sys_struct; vector_scale = 0.03)
display(scene)

nothing
