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

"""
    twist_inertia(sam, model, idx) -> Float64

The inertia of station `idx`'s twist about its hinge [kg m²], as the model's
[`station_dynamics`](@ref) has it at the current state.
"""
function twist_inertia(sam::SymbolicAWEModel, model, idx)
    integrator = sam.integrator
    inputs = refresh_outputs!(integrator.f.f, integrator.u, integrator.p, integrator.t).input
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
        repeat([[points[idx].total_mass + points[idx].apparent_mass for idx in free_points];
                [bodies[idx].total_mass + bodies[idx].apparent_mass for idx in free_bodies]];
               inner = 3),
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
    "the three rows of linearised rotation about the base that hold the pose"
    pose::Matrix{Float64}
end

function EquilibriumProblem(sam::SymbolicAWEModel)
    layout = EquilibriumLayout(sam)
    (; transforms) = sam.state
    base = collect(KiteGeometry.base_position(KiteGeometry.Pose(sam.state), transforms,
                                              first(transforms)))
    rows = [layout.velocity; vec(layout.spin); layout.scalar_rate]
    inertia = [layout.mass; vec(layout.inertia); layout.scalar_inertia]
    return EquilibriumProblem(sam, layout, inertia, rows, base,
                              pose_rows(layout, sam.integrator.u, base))
end

"""
    pose_rows(layout, u, base) -> Matrix

The 3×n map from a coordinate step to the mass-weighted rotation it gives the kite about
`base`, the angular momentum a unit-rate step would carry [kg m²]. Tethers are left out.
"""
function pose_rows(layout::EquilibriumLayout, u, base)
    translational = [layout.mass[3i] * layout.on_kite[3i] .* skew(u[layout.position[3i-2:3i]] .- base)
                     for i in 1:length(layout.position)÷3]
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

`u` with the coordinates moved by `step`: translations and splits added, bodies turned
about their own axes.
"""
function move(problem::EquilibriumProblem, u, step)
    (; layout) = problem
    moved = copy(u)
    n_position, n_body = length(layout.position), size(layout.quaternion, 2)
    moved[layout.position] .+= step[1:n_position]
    for (k, quat) in enumerate(eachcol(layout.quaternion))
        turn = step[n_position + 3k - 2:n_position + 3k]
        turned = quaternion_product(u[quat], [1.0; 0.5 .* turn])
        moved[quat] .= turned ./ norm(turned)
    end
    moved[layout.scalar] .+= step[n_position + 3n_body + 1:end]
    return moved
end

"""
    pose_moment(problem, u, force) -> Vector

The moment about the base of the generalized `force` at `u` [N m]: what has to hold the
pose against the loads.
"""
function pose_moment(problem::EquilibriumProblem, u, force)
    (; layout, base) = problem
    n_position = length(layout.position)
    moment = zeros(3)
    for i in 1:n_position÷3
        slots = 3i-2:3i
        moment .+= cross(u[layout.position[slots]] .- base, force[slots])
    end
    for (k, quat) in enumerate(eachcol(layout.quaternion))
        moment .+= quaternion_to_rotation_matrix(u[quat]) *
                   force[n_position + 3k - 2:n_position + 3k]
    end
    return moment
end

"""
    unbalanced(problem, force) -> Vector

The part of the generalized `force` that holding the pose does not take, in the least
squares sense: what still deforms the structure.
"""
function unbalanced(problem::EquilibriumProblem, force)
    (; pose) = problem
    return force .- pose' * ((pose * pose') \ (pose * force))
end

"""
    largest_load(problem, force) -> Float64

The largest force [N] or moment [N m] among the coordinates of `force`.
"""
function largest_load(problem::EquilibriumProblem, force)
    (; layout) = problem
    n_position = length(layout.position)
    blocks = [norm(force[i:i+2]) for i in 1:3:n_position + 3size(layout.quaternion, 2)]
    return maximum([blocks; abs.(force[n_position + 3size(layout.quaternion, 2) + 1:end])];
                   init = 0.0)
end

"""
    newton_step(problem, tangent, force, damping) -> Vector

The step that zeroes `force` on a linear structure of stiffness `tangent` with the pose
held, `damping` [1/s²] times each coordinate's inertia added to the stiffness.
"""
function newton_step(problem::EquilibriumProblem, tangent, force, damping)
    (; pose, inertia) = problem
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
    (; layout) = problem
    n_position, n_body = length(layout.position), size(layout.quaternion, 2)
    largest_move = maximum([norm(step[i:i+2]) for i in 1:3:n_position]; init = 0.0)
    largest_turn = maximum([norm(step[n_position+i:n_position+i+2]) for i in 1:3:3n_body];
                           init = 0.0)
    return step .* min(1.0, move_limit / largest_move, turn_limit / largest_turn)
end

"""
    find_equilibrium!(sam; rtol=1e-4, max_iterations=100, move_limit=0.25,
                      turn_limit=0.2) -> NamedTuple

Move `sam` to the equilibrium of its structure under the loads it carries now, at rest
and with the pose of its first transform held about the base: the minimum of the
structure's potential energy, the aerodynamic loads acting as dead loads at their last
solve. Steps no point further than `move_limit` [m] and no body further than
`turn_limit` [rad]. Stops once no coordinate carries more than `rtol` times the kite's
weight plus its aerodynamic force, in [N] or [N m]. Needs the `KernelBackend`.

Returns `(; converged, iterations, residual, moment)`: the largest unbalanced load left
[N] and the moment about the base [N m, world frame] that holding the pose takes, zero
where the pose is itself an equilibrium.
"""
function find_equilibrium!(sam::SymbolicAWEModel; rtol = 1e-4, max_iterations = 100,
                           move_limit = 0.25, turn_limit = 0.2)
    problem = EquilibriumProblem(sam)
    tolerance = rtol * kite_load(sam)
    u = at_rest(problem, copy(sam.integrator.u))
    force = generalized_force(problem, u)
    residual = largest_load(problem, unbalanced(problem, force))
    damping, iterations = 0.0, 0
    while residual > tolerance && iterations < max_iterations
        iterations += 1
        (; step, damping, smallest) = descent_step(problem, u, force, damping)
        step = limit_step(problem, step, move_limit, turn_limit)
        α = line_search(problem, u, step, -dot(force, step))
        u = move(problem, u, α .* step)
        force = generalized_force(problem, u)
        residual = largest_load(problem, unbalanced(problem, force))
        damping = α < 0.25 ? max(4damping, smallest) : α == 1 ? damping / 2 : damping
    end
    settle!(sam, u)
    return (; converged = residual <= tolerance, iterations, residual,
            moment = pose_moment(problem, u, force))
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
        step = newton_step(problem, tangent, force, damping)
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
