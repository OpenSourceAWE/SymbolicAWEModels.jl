# Copyright (c) 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

"""
    EquilibriumLayout

The coordinates [`find_equilibrium!`](@ref) moves, as slots of the ODE state, and the
inertia that turns each one's acceleration into the force on it.
"""
struct EquilibriumLayout
    "state slots of each translational coordinate: free points, then body centres of mass"
    position::Vector{Int}
    "state slots of the matching velocities"
    velocity::Vector{Int}
    "mass of each translational coordinate, apparent mass included [kg]"
    mass::Vector{Float64}
    "whether each translational coordinate belongs to the kite rather than a tether"
    on_kite::Vector{Bool}
    "state slots of each dynamic body's principal-frame quaternion, one column per body"
    quaternion::Matrix{Int}
    "state slots of each dynamic body's principal-frame angular velocity"
    spin::Matrix{Int}
    "principal moments of inertia of each dynamic body, one column per body [kg m²]"
    inertia::Matrix{Float64}
    "state slots of each scalar coordinate: free pulley rope splits, then station twists"
    scalar::Vector{Int}
    "state slots of the matching rates"
    scalar_rate::Vector{Int}
    "inertia of each scalar coordinate [kg, kg m²]"
    scalar_inertia::Vector{Float64}
end

"""
    state_slots(model, instance, name) -> Vector{Int}

The ODE state slots of the variable `name` of the kernel `instance` of `model`.
"""
state_slots(model, instance, name) =
    collect(buffer_slots(model.system, instance, :states, name))

"""
    tether_points(state) -> Set{Int}

The indices of the points the tether segments of `state` join.
"""
tether_points(state) = Set(Iterators.flatten(state.segments[idx].point_idxs
                                             for tether in state.tethers
                                             for idx in tether.segment_idxs))

"""The mass `component` accelerates with: its own plus the air it drags along [kg]."""
inertial_mass(component) = component.total_mass + component.apparent_mass

"""
    twist_inertia(sam, model, idx) -> Float64

The inertia of station `idx`'s twist about its hinge [kg m²], as the model's
[`station_dynamics`](@ref) has it at the current state.
"""
function twist_inertia(sam::SymbolicAWEModel, model, idx)
    (; f, u, p, t) = sam.integrator
    inputs = refresh_outputs!(f.f, u, p, t).input
    node_mass = inputs[only(buffer_slots(model.system, model.twist_instances[idx],
                                         :inputs, :node_mass_in))]
    station = sam.state.stations[idx]
    return station_dynamics(; free_angle = 0.0, twist_vel = 0.0, aero_moment = 0.0,
                            node_moment = 0.0, mass = node_mass + station.body_mass,
                            chord = collect(station.chord), damping = 0.0,
                            stiffness = 0.0).inertia
end

"""
    EquilibriumLayout(sam)

The coordinates of `sam`, read off its kernel model and weighed with the masses its
state holds.
"""
function EquilibriumLayout(sam::SymbolicAWEModel)
    sam.backend isa KernelBackend ||
        error("find_equilibrium! needs a model built with the KernelBackend.")
    (; points, bodies, pulleys, stations) = sam.state
    model = sam.prob.initial_sync.model
    free_points = [idx for (idx, role) in enumerate(model.point_roles)
                   if role.kind in (:particle, :wing_node, :pulley)]
    free_bodies = [idx for idx in eachindex(bodies) if bodies[idx].type == DYNAMIC]
    free_pulleys = [(idx, role) for (idx, role) in enumerate(model.point_roles)
                    if role.kind === :pulley && pulleys[role.pulley_idx].brake <= 0.5]
    twists = [idx for idx in eachindex(stations) if stations[idx].type == DYNAMIC]
    parts = [[(model.point_instances[idx], :pos, :vel) for idx in free_points];
             [(model.body_instances[idx], :com_w, :com_vel) for idx in free_bodies]]
    scalars = [[(model.point_instances[idx], :pulley_len, :pulley_vel)
                for (idx, _) in free_pulleys];
               [(model.twist_instances[idx], :free_twist_angle, :twist_omega)
                for idx in twists]]
    tether = tether_points(sam.state)
    return EquilibriumLayout(
        reduce(vcat, [state_slots(model, instance, pos) for (instance, pos, _) in parts];
               init = Int[]),
        reduce(vcat, [state_slots(model, instance, vel) for (instance, _, vel) in parts];
               init = Int[]),
        repeat([[inertial_mass(points[idx]) for idx in free_points];
                [inertial_mass(bodies[idx]) for idx in free_bodies]]; inner = 3),
        repeat([[!(idx in tether) for idx in free_points]; trues(length(free_bodies))];
               inner = 3),
        reduce(hcat, [state_slots(model, model.body_instances[idx], :Q)
                      for idx in free_bodies]; init = zeros(Int, 4, 0)),
        reduce(hcat, [state_slots(model, model.body_instances[idx], :omega_p)
                      for idx in free_bodies]; init = zeros(Int, 3, 0)),
        reduce(hcat, [collect(bodies[idx].inertia_principal) for idx in free_bodies];
               init = zeros(3, 0)),
        [only(state_slots(model, instance, len)) for (instance, len, _) in scalars],
        [only(state_slots(model, instance, rate)) for (instance, _, rate) in scalars],
        [[pulley_rope_mass(sam.state, role.pulley_idx, role.segment_idx)
          for (_, role) in free_pulleys];
         [twist_inertia(sam, model, idx) for idx in twists]])
end

"""Hamilton product of two scalar-first quaternions."""
quaternion_product(q, p) = [q[1] * p[1] - dot(q[2:4], p[2:4]);
                            q[1] .* p[2:4] .+ p[1] .* q[2:4] .+ cross(q[2:4], p[2:4])]

"""
    quaternion_turn(q) -> Matrix

How the scalar-first quaternion `q` changes per small rotation about each of its own
axes [rad]: the 4×3 matrix `dq/dθ`.
"""
quaternion_turn(q) =
    reduce(hcat, [quaternion_product(q, [0.0; 0.5 .* axis]) for axis in eachcol(I(3))])

"""
    skew(r) -> Matrix

The matrix that takes `x` to `r × x`.
"""
skew(r) = [0 -r[3] r[2]; r[3] 0 -r[1]; -r[2] r[1] 0]

"""
    coordinate_ranges(layout) -> (; translation, rotation, scalar)

Where the translational, rotational and scalar coordinates of `layout` sit in a
coordinate vector.
"""
function coordinate_ranges(layout::EquilibriumLayout)
    n_translation = length(layout.position)
    n_rotation = 3size(layout.quaternion, 2)
    return (; translation = 1:n_translation,
            rotation = n_translation + 1:n_translation + n_rotation,
            scalar = n_translation + n_rotation + 1:
                     n_translation + n_rotation + length(layout.scalar))
end

"""The consecutive triples of `range`, one per vector coordinate."""
triples(range) = (range[i:i+2] for i in 1:3:length(range))

"""The largest norm among the vectors of `vector` at the triples of `range`."""
largest_norm(vector, range) =
    maximum((norm(vector[slots]) for slots in triples(range)); init = 0.0)

"""
    EquilibriumProblem(sam)

The equilibrium of `sam` at its current pose, held about the base of its first transform.
"""
struct EquilibriumProblem{S}
    "the model whose right-hand side and Jacobian give the forces"
    sam::S
    "where the coordinates sit in the ODE state"
    layout::EquilibriumLayout
    "inertia of each coordinate, in the order of the generalized force"
    inertia::Vector{Float64}
    "rows of the ODE state whose rates are the accelerations"
    rows::Vector{Int}
    "ground-station position the pose turns about [m]"
    base::Vector{Float64}
end

function EquilibriumProblem(sam::SymbolicAWEModel)
    layout = EquilibriumLayout(sam)
    (; transforms) = sam.state
    base = collect(KiteGeometry.base_position(KiteGeometry.Pose(sam.state), transforms,
                                              first(transforms)))
    rows = [layout.velocity; vec(layout.spin); layout.scalar_rate]
    inertia = [layout.mass; vec(layout.inertia); layout.scalar_inertia]
    return EquilibriumProblem(sam, layout, inertia, rows, base)
end

"""
    pose_rows(problem, u) -> Matrix

The 3×n map from a coordinate step at `u` to the mass-weighted rotation it gives the
kite about the base, the angular momentum a unit-rate step would carry [kg m²]. Tethers
are left out.
"""
function pose_rows(problem::EquilibriumProblem, u)
    (; layout, base) = problem
    translational = [layout.mass[first(slots)] * layout.on_kite[first(slots)] .*
                     skew(u[layout.position[slots]] .- base)
                     for slots in triples(coordinate_ranges(layout).translation)]
    rotational = [quaternion_to_rotation_matrix(u[quat]) * Diagonal(inertia)
                  for (quat, inertia) in zip(eachcol(layout.quaternion),
                                             eachcol(layout.inertia))]
    return [reduce(hcat, translational; init = zeros(3, 0)) reduce(hcat, rotational;
            init = zeros(3, 0)) zeros(3, length(layout.scalar))]
end

"""
    at_rest(problem, u) -> Vector

`u` with every coordinate's rate set to zero.
"""
function at_rest(problem::EquilibriumProblem, u)
    rest = copy(u)
    rest[problem.rows] .= 0
    return rest
end

"""
    generalized_force(problem, u) -> Vector

The force on each coordinate at rest at `u` [N, N m]: the model's own acceleration times
the coordinate's inertia. Rest kills every damping term, so damping cannot shift it.
"""
function generalized_force(problem::EquilibriumProblem, u)
    integrator = problem.sam.integrator
    du = similar(u)
    integrator.f(du, at_rest(problem, u), integrator.p, integrator.t)
    return problem.inertia .* du[problem.rows]
end

"""
    stiffness(problem, u) -> Matrix

Minus the derivative of [`generalized_force`](@ref) with respect to the coordinates:
the tangent stiffness, from the model's own Jacobian.
"""
function stiffness(problem::EquilibriumProblem, u)
    (; layout, rows, inertia) = problem
    integrator = problem.sam.integrator
    rest = at_rest(problem, u)
    jacobian = similar(integrator.f.jac_prototype, Float64)
    integrator.f.jac(jacobian, rest, integrator.p, integrator.t)
    rates = jacobian[rows, :]
    turns = [Matrix(rates[:, quat]) * quaternion_turn(rest[quat])
             for quat in eachcol(layout.quaternion)]
    return -inertia .* [Matrix(rates[:, layout.position]) reduce(hcat, turns;
                        init = zeros(length(rows), 0)) Matrix(rates[:, layout.scalar])]
end

"""
    move(problem, u, step) -> Vector

`u` with the coordinates moved by `step`: translations and scalar coordinates added,
bodies turned about their own axes.
"""
function move(problem::EquilibriumProblem, u, step)
    (; layout) = problem
    ranges = coordinate_ranges(layout)
    moved = copy(u)
    moved[layout.position] .+= step[ranges.translation]
    for (quat, turn) in zip(eachcol(layout.quaternion), triples(ranges.rotation))
        turned = quaternion_product(u[quat], [1.0; 0.5 .* step[turn]])
        moved[quat] .= turned ./ norm(turned)
    end
    moved[layout.scalar] .+= step[ranges.scalar]
    return moved
end

"""
    pose_moment(problem, u, force) -> Vector

The moment about the base of the generalized `force` at `u` [N m]: what has to hold the
pose against the loads.
"""
function pose_moment(problem::EquilibriumProblem, u, force)
    (; layout, base) = problem
    ranges = coordinate_ranges(layout)
    moment = zeros(3)
    for slots in triples(ranges.translation)
        moment .+= cross(u[layout.position[slots]] .- base, force[slots])
    end
    for (quat, turn) in zip(eachcol(layout.quaternion), triples(ranges.rotation))
        moment .+= quaternion_to_rotation_matrix(u[quat]) * force[turn]
    end
    return moment
end

"""
    unbalanced(problem, u, force) -> Vector

The part of the generalized `force` at `u` that holding the pose does not take, in the
least squares sense: what still deforms the structure.
"""
function unbalanced(problem::EquilibriumProblem, u, force)
    pose = pose_rows(problem, u)
    return force .- pose' * ((pose * pose') \ (pose * force))
end

"""
    largest_load(problem, force) -> Float64

The largest force [N] or moment [N m] among the coordinates of `force`.
"""
function largest_load(problem::EquilibriumProblem, force)
    ranges = coordinate_ranges(problem.layout)
    return max(largest_norm(force, ranges.translation),
               largest_norm(force, ranges.rotation),
               maximum(abs, force[ranges.scalar]; init = 0.0))
end

"""
    newton_step(problem, u, tangent, force, damping) -> Vector

The step from `u` that zeroes `force` on a linear structure of stiffness `tangent` with
the pose held, `damping` [1/s²] times each coordinate's inertia added to the stiffness.
"""
function newton_step(problem::EquilibriumProblem, u, tangent, force, damping)
    (; inertia) = problem
    pose = pose_rows(problem, u)
    system = [tangent + damping * Diagonal(inertia) pose'; pose zeros(3, 3)]
    return (system \ [force; zeros(3)])[1:length(force)]
end

"""
    energy_slope(problem, u, step, α) -> Float64

The rate of change of the potential energy along `step` at `u` moved by `α` of it:
minus the work the generalized force does on it.
"""
energy_slope(problem::EquilibriumProblem, u, step, α) =
    -dot(generalized_force(problem, move(problem, u, α .* step)), step)

"""
    line_search(problem, u, step, slope) -> Float64

The fraction of `step` to take: all of it while the energy still falls at its end,
otherwise one where the energy's `slope` has shrunk to half its size at `u`.
"""
function line_search(problem::EquilibriumProblem, u, step, slope; max_evaluations = 8)
    low, high = 0.0, 1.0
    slope_low, slope_high = slope, energy_slope(problem, u, step, 1.0)
    slope_high <= 0 && return 1.0
    α = high
    for _ in 1:max_evaluations
        α = low + (high - low) * slope_low / (slope_low - slope_high)
        slope_α = energy_slope(problem, u, step, α)
        abs(slope_α) <= 0.5abs(slope) && return α
        slope_α < 0 ? (low, slope_low = α, slope_α) : (high, slope_high = α, slope_α)
    end
    return α
end

"""
    limit_step(problem, step, move_limit, turn_limit) -> Vector

`step` scaled down so no point moves further than `move_limit` [m] and no body turns
further than `turn_limit` [rad].
"""
function limit_step(problem::EquilibriumProblem, step, move_limit, turn_limit)
    ranges = coordinate_ranges(problem.layout)
    return step .* min(1.0, move_limit / largest_norm(step, ranges.translation),
                       turn_limit / largest_norm(step, ranges.rotation))
end

"""
    find_equilibrium!(sam; rtol=1e-4, max_rounds=30, max_iterations=100,
                      move_limit=0.25, turn_limit=0.2) -> NamedTuple

Move `sam` to its quasi-steady state at rest, with the pose of its first transform held
about the base: the structure at the minimum of its potential energy under the
aerodynamics solved on its own shape. Each round settles the structure under the
aerodynamics of the current shape as dead loads, in at most `max_iterations` Newton
steps of no more than `move_limit` [m] per point and `turn_limit` [rad] per body, and
solves the aerodynamics on the settled shape; the next round starts from the shape
[`next_input`](@ref) predicts from the last five. Stops once no coordinate carries more
than `rtol` times the kite's weight plus its aerodynamic force, in [N] or [N m], or
after `max_rounds` rounds. Needs the `KernelBackend`.

Returns `(; converged, rounds, iterations, residual, moment)`: the rounds and Newton
steps taken, the largest unbalanced load left [N] and the moment about the base [N m,
world frame] that holding the pose takes, zero where the pose is itself an equilibrium.
"""
function find_equilibrium!(sam::SymbolicAWEModel; rtol = 1e-4, max_rounds = 30,
                           max_iterations = 100, move_limit = 0.25, turn_limit = 0.2)
    problem = EquilibriumProblem(sam)
    reference = at_rest(problem, copy(sam.integrator.u))
    u, loads = reference, loads_at!(problem, reference, rtol)
    inputs, outputs = Vector{Float64}[], Vector{Float64}[]
    rounds, iterations = 0, 0
    while true
        rounds += 1
        push!(inputs, coordinate_change(problem, reference, u))
        settled = settle_structure(problem, u, loads.force, loads.tolerance;
                                   max_iterations, move_limit, turn_limit)
        iterations += settled.iterations
        u, loads = settled.u, loads_at!(problem, settled.u, rtol)
        (loads.residual <= loads.tolerance || rounds == max_rounds) && break
        push!(outputs, coordinate_change(problem, reference, u))
        length(inputs) > 5 && (popfirst!(inputs); popfirst!(outputs))
        u = move(problem, reference, next_input(inputs, outputs))
        loads = loads_at!(problem, u, rtol)
    end
    return (; converged = loads.residual <= loads.tolerance, rounds, iterations,
            residual = loads.residual, moment = pose_moment(problem, u, loads.force))
end

"""
    loads_at!(problem, u, rtol) -> (; force, residual, tolerance)

Put the model at `u`, solve its aerodynamics there and return the generalized force,
its largest unbalanced load and the residual it must fall below [N, N m].
"""
function loads_at!(problem::EquilibriumProblem, u, rtol)
    (; sam) = problem
    settle!(sam, u)
    if has_vsm_wing(sam.state)
        refresh_aero!(sam; cold_start = true, vsm_warn_on_fail = true)
        sync_params!(sam.prob.param_sync, sam.integrator, sam.state)
    end
    force = generalized_force(problem, u)
    return (; force, residual = largest_load(problem, unbalanced(problem, u, force)),
            tolerance = rtol * kite_load(sam))
end

"""
    settle_structure(problem, u, force, tolerance; max_iterations, move_limit,
                     turn_limit) -> (; u, iterations)

The minimum of the structure's potential energy from `u`, where it carries `force`,
under the loads the model holds now: damped Newton steps until no coordinate carries
more than `tolerance`.
"""
function settle_structure(problem::EquilibriumProblem, u, force, tolerance;
                          max_iterations, move_limit, turn_limit)
    residual = largest_load(problem, unbalanced(problem, u, force))
    damping, iterations = 0.0, 0
    while residual > tolerance && iterations < max_iterations
        iterations += 1
        (; step, damping, smallest) = descent_step(problem, u, force, damping)
        step = limit_step(problem, step, move_limit, turn_limit)
        α = line_search(problem, u, step, -dot(force, step))
        u = move(problem, u, α .* step)
        force = generalized_force(problem, u)
        residual = largest_load(problem, unbalanced(problem, u, force))
        damping = α < 0.25 ? max(4damping, smallest) : α == 1 ? damping / 2 : damping
    end
    return (; u, iterations)
end

"""
    coordinate_change(problem, u, target) -> Vector

The step that [`move`](@ref) takes from `u` to `target`.
"""
function coordinate_change(problem::EquilibriumProblem, u, target)
    (; layout) = problem
    ranges = coordinate_ranges(layout)
    change = zeros(length(problem.inertia))
    change[ranges.translation] .= target[layout.position] .- u[layout.position]
    for (quat, turn) in zip(eachcol(layout.quaternion), triples(ranges.rotation))
        relative = quaternion_product([u[quat[1]]; -u[quat[2:4]]], target[quat])
        change[turn] .= 2 .* relative[2:4] ./ relative[1]
    end
    change[ranges.scalar] .= target[layout.scalar] .- u[layout.scalar]
    return change
end

"""
    next_input(inputs, outputs; relaxation=0.1) -> Vector

The next input of the fixed point that maps each of `inputs` to the matching one of
`outputs`, by interface quasi-Newton with least squares (IQN-ILS): the latest output
corrected by how outputs and their residuals changed between rounds, or the first input
moved `relaxation` of the way to its output where there is one round.
"""
function next_input(inputs, outputs; relaxation = 0.1)
    residuals = outputs .- inputs
    n = length(inputs)
    n == 1 && return inputs[1] .+ relaxation .* residuals[1]
    residual_changes = reduce(hcat, [residuals[n] .- residuals[i] for i in 1:(n - 1)])
    output_changes = reduce(hcat, [outputs[n] .- outputs[i] for i in 1:(n - 1)])
    return outputs[n] .+ output_changes * (residual_changes \ -residuals[n])
end

"""
    descent_step(problem, u, force, damping) -> (; step, damping, smallest)

A [`newton_step`](@ref) along which the potential energy falls, raising `damping` from
its given value until it does, and the `smallest` damping worth trying: a millionth of
the stiffest coordinate's stiffness per inertia.
"""
function descent_step(problem::EquilibriumProblem, u, force, damping)
    tangent = stiffness(problem, u)
    smallest = 1e-6 * maximum(abs.(diag(tangent)) ./ problem.inertia)
    for _ in 1:30
        step = newton_step(problem, u, tangent, force, damping)
        all(isfinite, step) && dot(force, step) > 0 && return (; step, damping, smallest)
        damping = max(10damping, smallest)
    end
    error("find_equilibrium!: no step lowers the potential energy.")
end

"""
    kite_load(sam) -> Float64

The kite's weight plus the size of its wings' aerodynamic force [N], the scale a force
residual is judged against.
"""
kite_load(sam::SymbolicAWEModel) =
    sam.set.g_earth * sum(point -> point.total_mass, sam.state.points) +
    sum(wing -> norm(wing.aero_force_b), sam.state.wings; init = 0.0)

"""
    settle!(sam, u)

Put `sam`'s integrator and state at `u`.
"""
function settle!(sam::SymbolicAWEModel, u)
    SciMLBase.set_u!(sam.integrator, copy(u))
    update_sys_struct!(sam.prob, sam.integrator, sam.state)
    return sam
end
