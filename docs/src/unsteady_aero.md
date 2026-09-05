# Unsteady aerodynamics

The VSM solve is steady: it answers what the loading would be if the wing held its
current shape and inflow forever. Two corrections restore what a wing actually feels
while that state is changing. Both are off by default and live on the wing's
`UnsteadyAero`, held by its `VSMEngine`:

```julia
sys.wings[1].unsteady.apparent_mass = 1.0   # thin-plate value
sys.wings[1].unsteady.wagner = true
```

They are independent — either can be used alone.

## Apparent mass

### What it models

A plate cannot accelerate normal to itself without accelerating the air around it.
That entrained air resists the motion exactly as extra inertia would. For a 2D flat
plate of chord `c` the added mass per unit span is the mass of a cylinder of air of
diameter `c`:

```
m_a' = ρ·π·(c/2)² = ρ·π·c²/4
```

Over a panel of width `w` that is `m_a = ρ·π·c²·w/4`, which is
`panel_apparent_mass`. The chord and width are the build-time twins of the `chord`
and `width` equations in `panel_force_eqs`, read off the same frozen VSM mesh, so a
panel's added mass and its lift describe the same panel.

This matters more for kites than the term's obscurity suggests. A ram-air canopy is
light and broad, so the air it drags along is a large fraction of — often a multiple
of — its own mass, and leaving it out overstates every acceleration the aerodynamic
loads produce.

### How it is implemented

`apply_apparent_mass!` spreads each panel's `m_a` over the wing's structural nodes
using `aero_scatter_entries` — the *same* weights that already carry that panel's
force. The air therefore arrives exactly where the lift does, and no second
distribution scheme has to be kept in step with the first. The weights sum to one
per panel, so the wing's nodes receive the full panel mass and nothing more.

It lands in a field of its own, `apparent_mass`, rather than in `extra_mass`. This is
the whole point of the term and the reason it could not be bolted onto the existing
mass: **entrained air resists acceleration but has no weight**. Gravity is built from
`mass`, so inflating `extra_mass` would have hung a canopy's worth of nonexistent
weight off the wing. Instead the acceleration divides by the sum and gravity is left
alone:

```julia
accel = net_force ./ (mass .+ apparent_mass) .- world_damping .* vel
```

### Where the air actually goes

The air has to slow down whatever *integrates* the node's translation, and that is
often not the node. `apparent_mass_carriers` picks it, and there are three cases:

- a free `DYNAMIC` node carries its own air, through `point_acceleration`;
- a node anchored to a rigid body is *placed* by that body and integrates nothing of
  its own, so its air goes to the body, through `rigid_body_pose_expressions`;
- a node riding a `TimoshenkoJoint`'s deformed centerline is placed by the two bodies
  the joint spans, which split its air by `beam_frac`, where along the element it
  sits.

The last two are not edge cases. **No node of a beam wing integrates itself**: a beam
wing keeps its mass and its motion in its beam bodies, its leading- and trailing-edge
stations are anchored to those bodies, and the great majority of its nodes ride the
strut joints between them. Mass left on such a node would be silently inert — it
would look applied and change nothing.

Both backends reach every path through shared kernels: points via `point_eqs.jl`
(monolith) and `dynamic_point_dynamics`/`WingNodePoint` (kernel), bodies via
`rigid_body_pose_expressions` in either. All of them divide the net force while
gravity keeps reading `mass`.

Only the *translational* added mass is modelled. A plate also has an added moment of
inertia about its pitch axis (`ρπb⁴/8` per unit span); that is not applied, so a
body's rotational inertia is untouched.

`apply_apparent_mass!` runs from `refresh_aero!`, next to `sync_aero_density!` and
at the density that call just read, so the entrained air tracks altitude. It clears
every node before it distributes, so switching the correction off and rebuilding
cannot leave stale inertia behind. It is a parameter, not structure, and so does
not enter the model hash.

### Approximation

The added mass is applied **isotropically**. A flat plate's true added-mass tensor
is `ρπb²` normal to the plate and near zero in the chordwise and spanwise
directions; the nodes here get the normal-direction value on all three axes. The
normal direction is the one that carries the lift and the one the correction exists
for, so the term is right where it matters and too stiff in-plane.

The alternative — an acceleration-dependent force with the correct anisotropy —
makes `D(vel)` appear on both sides of its own equation. That is a linearly implicit
DAE, and rather than a mass-matrix change it becomes an algebraic loop the solver
must iterate at every RHS evaluation, over every node of the wing. The isotropic
lumped form is unconditionally stable, costs one extra scalar per point, and keeps
the RHS explicit. If the in-plane over-stiffness ever matters, the honest fix is a
per-node anisotropic mass matrix in the wing frame, not a scale factor.

The `apparent_mass` field is a *scale*, not a switch: `1.0` is the thin-plate value,
`0.0` disables the term, and values in between trade physical fidelity for a weaker
correction. It errors on a `RIGID_DYNAMICS` wing, where the whole wing is one body
already carrying the load and the per-node scatter has no meaning.

## Wagner lag

### What it models

Circulation does not appear the instant the angle of attack does. The vorticity a
change in incidence sheds has to convect downstream before the bound circulation
reaches its new steady value. Wagner's indicial function gives the fraction of the
final circulatory lift reached after `s` semi-chords of travel:

```
φ(s) = 1 - A₁·exp(-b₁·s) - A₂·exp(-b₂·s)
```

with R. T. Jones' fit `A₁ = 0.165`, `A₂ = 0.335`, `b₁ = 0.0455`, `b₂ = 0.3` as the
defaults. Note `φ(0) = 1 - A₁ - A₂ = 0.5`: a step in incidence produces **half** its
steady lift immediately, and the rest builds over the following chord lengths.

### How it is implemented

The obvious realisation of the Duhamel integral needs `dα/dt` on the right-hand
side, which is exactly the kind of state derivative that should not appear in an
ODE right-hand side. A change of variable removes it. Drive two states with

```
dxᵢ/ds = α - bᵢ·xᵢ
```

and read the deficiency off them as

```
d = Σ Aᵢ·(α - bᵢ·xᵢ),      α_E = α - d
```

After a step in `α` this reproduces `α_E = α·φ(s)`, and in steady flow every
`xᵢ → α/bᵢ`, so `d → 0` and the lag disappears — a trimmed wing is untouched. In
time, with `s = ∫2·v/c·dt`:

```
D(xᵢ) = (2·v/c)·(α - bᵢ·xᵢ)
```

which is what `wagner_lag_eqs` emits.

### One lag for the whole wing

There are **two states for the entire wing**, not two per panel. The reference angle
is the wing's own: `wagner_reference_frame` averages the frozen VSM mesh into a mean
chordwise direction, a mean normal and a mean chord, and the lag measures its `α`
and its speed against the wing's mean body-frame apparent wind.

The single resulting `deficiency` is then subtracted from *every* panel's angle of
attack before the polars are read:

```
alpha_eff[i] = alpha[i] - deficiency
```

Because the shift is common to all panels, the spanwise shape of the loading is
untouched — only its unsteady overall response lags. Two states buy the whole
correction and the cost does not grow with panel count.

This also keeps the dependency acyclic, which is what makes it work in the kernel
backend. The reference `α` is built from the wing's inflow, not from panel outputs,
so nothing feeds a panel's own result back into its input.

The lag is applied to the **polar lookup** — `cl`, `cd` and `cm` all read
`alpha_eff` — while the geometric `alpha` still sets the lift and drag *directions*.
Lagging the coefficients is the physical claim (circulation builds over time);
rotating the force vectors by a lagged angle is not, and the two are separable here
because `dir_lift`/`dir_drag` come from the geometry rather than from the polars.

### Both backends

The physics is written once, in `wagner_lag_eqs`. The monolith path calls it through
`wagner_wing_eqs` from both particle aero modes' `aero_component`. The kernel
backend wraps the same function in a `WagnerLag` component, one instance per wing,
whose `va_in` is gathered at equal weight over the wing's nodes and whose
`deficiency` output feeds every `AeroPanel` as an input. `ParticleWingAero` calls
`aero_component` and so inherits the monolith form.

The four lag constants are **registered parameters**, read from
`wing.unsteady` like any other tunable scalar in the model. Only the `wagner` switch
itself is structure — it is what adds or removes the two states — so it alone reaches
the model cache key. Retuning a lag is a `sync_params!`, not a rebuild, which is what
makes sweeping the constants cheap.

## Tuning the lag

The rates are per semi-chord travelled, so the physical time constants scale with
`c/(2·v)`. Slower rates mean a longer-lived lag. None of these force a rebuild, so a
sweep can set them between runs on one compiled model:

```julia
sys.wings[1].unsteady.wagner_rates = [0.0455, 0.3] .* 0.5   # twice as slow
sys.wings[1].unsteady.wagner_gains = [0.165, 0.335]         # φ(0) = 0.5
```

The gains set how much lift is withheld initially: `φ(0) = 1 - A₁ - A₂`. Raising
their sum toward 1 makes the wing start from nothing after a step; lowering it
toward 0 removes the lag's initial bite while keeping its tail.
