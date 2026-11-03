# Copyright (c) 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

"""
    RigidLayout

Where [`find_quasi_steady!`](@ref) reads and writes the ODE state: the position and
velocity slots of every free point, and the pose and rate slots of every dynamic body.
"""
struct RigidLayout
    "indices in `sam.state.points` of the points carrying state"
    points::Vector{Int}
    "state slots of each point's position, one column per point"
    point_pos::Matrix{Int}
    "state slots of each point's velocity"
    point_vel::Matrix{Int}
    "indices in `sam.state.bodies` of the dynamic bodies"
    bodies::Vector{Int}
    "state slots of each body's centre of mass"
    body_com::Matrix{Int}
    "state slots of each body's centre-of-mass velocity"
    body_vel::Matrix{Int}
    "state slots of each body's principal-frame quaternion, scalar first"
    body_quat::Matrix{Int}
    "state slots of each body's principal-frame angular velocity"
    body_spin::Matrix{Int}
    "position slots of every state point, then of every body's centre of mass"
    part_pos::Matrix{Int}
    "velocity slots in the order of `part_pos`"
    part_vel::Matrix{Int}
end

"""
    state_slots(model, instances, name) -> Matrix{Int}

The state slots of the variable `name` of each of `instances`, one column each.
"""
state_slots(model, instances, name) = reduce(hcat,
    [collect(buffer_slots(model.system, instance, :states, name))
     for instance in instances]; init = Matrix{Int}(undef, 3 + (name === :Q), 0))

"""
    RigidLayout(sam)

The state slots of `sam`'s points and dynamic bodies, looked up on its kernel model.
"""
function RigidLayout(sam::SymbolicAWEModel)
    sam.backend isa KernelBackend ||
        error("find_quasi_steady! needs a model built with the KernelBackend.")
    model = sam.prob.initial_sync.model
    points = [idx for (idx, role) in enumerate(model.point_roles)
              if role.kind in (:particle, :wing_node, :pulley)]
    bodies = [idx for idx in eachindex(model.body_instances)
              if sam.state.bodies[idx].type == DYNAMIC]
    point_instances = model.point_instances[points]
    body_instances = model.body_instances[bodies]
    point_pos = state_slots(model, point_instances, :pos)
    point_vel = state_slots(model, point_instances, :vel)
    body_com = state_slots(model, body_instances, :com_w)
    body_vel = state_slots(model, body_instances, :com_vel)
    return RigidLayout(points, point_pos, point_vel, bodies, body_com, body_vel,
        state_slots(model, body_instances, :Q),
        state_slots(model, body_instances, :omega_p),
        [point_pos body_com], [point_vel body_vel])
end

"""The three-vector in the state slots `slots` of `u`."""
state_vector(u, slots) = SVector{3}(u[slots[1]], u[slots[2]], u[slots[3]])

"""Hamilton product of two scalar-first quaternions."""
quat_product(q, p) = [q[1] * p[1] - dot(q[2:4], p[2:4]);
                      q[1] .* p[2:4] .+ p[1] .* q[2:4] .+ cross(q[2:4], p[2:4])]

"""
    set_rigid_velocities!(u, layout, spin, base) -> u

Put every point and body of `u` on the rigid rotation `spin` [rad/s] about `base`:
velocity `spin × (x − base)`, and the body spin `spin` resolved on its principal axes.
"""
function set_rigid_velocities!(u, layout::RigidLayout, spin, base)
    for (pos, vel) in zip(eachcol(layout.part_pos), eachcol(layout.part_vel))
        u[vel] .= cross(spin, u[pos] .- base)
    end
    for (quat, body_spin) in zip(eachcol(layout.body_quat), eachcol(layout.body_spin))
        u[body_spin] .= quaternion_to_rotation_matrix(u[quat])' * spin
    end
    return u
end

"""
    rotate_state!(u, layout, rotation, base) -> u

Turn the whole of `u` rigidly by the rotation matrix `rotation` about `base`:
positions, velocities and body orientations.
"""
function rotate_state!(u, layout::RigidLayout, rotation, base)
    for (pos, vel) in zip(eachcol(layout.part_pos), eachcol(layout.part_vel))
        u[pos] .= base .+ rotation * (u[pos] .- base)
        u[vel] .= rotation * u[vel]
    end
    turn = rotation_matrix_to_quaternion(rotation)
    for quat in eachcol(layout.body_quat)
        u[quat] .= quat_product(turn, u[quat])
    end
    return u
end

"""
    best_rotation(from, to) -> Matrix{Float64}

The rotation matrix that best carries the columns of `from` onto those of `to`, in
the least-squares sense.
"""
function best_rotation(from, to)
    factors = svd(to * from')
    return factors.U * Diagonal([1.0, 1.0, sign(det(factors.U * factors.Vt))]) *
           factors.Vt
end

"""
    SteadyProblem(sam)

One quasi-steady solve of `sam`, about the base of its first transform.
"""
struct SteadyProblem{S}
    "the model whose integrator carries the state"
    sam::S
    "where the state of each point and body sits in `u`"
    layout::RigidLayout
    "ground-station position the rigid rotation turns about [m]"
    base::SVector{3, Float64}
    "mass of each part in the order of `layout.part_pos` [kg]"
    mass::Vector{Float64}
end

function SteadyProblem(sam::SymbolicAWEModel)
    (; points, bodies, transforms) = sam.state
    layout = RigidLayout(sam)
    base = SVector{3, Float64}(KiteGeometry.base_position(KiteGeometry.Pose(sam.state),
                                                          transforms, first(transforms)))
    mass = Float64[[points[idx].total_mass for idx in layout.points];
                   [bodies[idx].total_mass for idx in layout.bodies]]
    return SteadyProblem(sam, layout, base, mass)
end

"""
    fitted_spin(problem, u) -> Vector{Float64}

The rigid rotation about the base [rad/s] that best fits the velocities of `u`,
each point and body weighed by its mass.
"""
function fitted_spin(problem::SteadyProblem, u)
    (; layout, base, mass) = problem
    inertia, momentum = zeros(3, 3), zeros(3)
    for (part_mass, pos, vel) in zip(mass, eachcol(layout.part_pos),
                                     eachcol(layout.part_vel))
        arm = u[pos] .- base
        inertia .+= part_mass .* (dot(arm, arm) .* I(3) .- arm * arm')
        momentum .+= part_mass .* cross(arm, u[vel])
    end
    return inertia \ momentum
end

"""Write the ODE state `u` into `sam.state`, leaving the integrator alone."""
function read_state!(problem::SteadyProblem, u)
    (; sam) = problem
    sam.prob.get_all_state((; u, p = sam.integrator.p, t = sam.integrator.t), sam.state)
    return nothing
end

"""
    hold_pose!(problem, u) -> u

Turn the whole of `u` rigidly about the base back onto the pose its transforms hold,
as `reposition!` places it, velocities and orientations turned with it.
"""
function hold_pose!(problem::SteadyProblem, u)
    (; sam, layout, base) = problem
    points = sam.state.points
    read_state!(problem, u)
    before = reduce(hcat, [points[idx].pos_w .- base for idx in layout.points])
    reposition!(sam.state)
    after = reduce(hcat, [points[idx].pos_w .- base for idx in layout.points])
    return rotate_state!(u, layout, best_rotation(before, after), base)
end

"""
    refresh_aero!(problem::SteadyProblem, u; cold_start=false) -> Bool

Write `u` into the model's integrator and state, re-solve the aerodynamics there, from
scratch with `cold_start`, and push the result into the parameters. An aero solve that
gives up returns `false`.
"""
function refresh_aero!(problem::SteadyProblem, u; cold_start = false)
    (; sam) = problem
    sam.integrator.u .= u
    read_state!(problem, u)
    try
        refresh_aero!(sam; cold_start)
    catch failure
        failure isa VortexStepMethod.SolveFailure || rethrow()
        return false
    end
    sync_params!(sam.prob.param_sync, sam.integrator, sam.state)
    return true
end

"""
    RigidDampedRHS

The model's right-hand side with every point and body damped against the rigid
rotation `spin` about `base`: `damping` [1/s] on the velocity left beyond
`spin × (x − base)`, and on the body spin left beyond `spin`.
"""
struct RigidDampedRHS{F}
    rhs::F
    layout::RigidLayout
    base::SVector{3, Float64}
    spin::Vector{Float64}
    damping::Float64
end

function (model::RigidDampedRHS)(du, u, p, t)
    (; layout, base, damping) = model
    model.rhs(du, u, p, t)
    spin = SVector{3}(model.spin)
    for (pos, vel) in zip(eachcol(layout.part_pos), eachcol(layout.part_vel))
        drift = state_vector(u, vel) - cross(spin, state_vector(u, pos) - base)
        @views du[vel] .-= damping .* drift
    end
    for (quat, body_spin) in zip(eachcol(layout.body_quat), eachcol(layout.body_spin))
        orientation = quaternion_to_rotation_matrix(view(u, quat))
        drift = state_vector(u, body_spin) - SMatrix{3, 3}(orientation)' * spin
        @views du[body_spin] .-= damping .* drift
    end
    return nothing
end

"""
    RigidDampedJacobian

The model's analytic state Jacobian with the diagonal of [`RigidDampedRHS`](@ref)'s
damping added, on a pattern that holds both.
"""
struct RigidDampedJacobian{J}
    jac::J
    "scratch the model's own Jacobian is written into"
    model_jacobian::SparseMatrixCSC{Float64, Int}
    "position in the joint pattern of each stored entry of `model_jacobian`"
    model_entries::Vector{Int}
    "positions in the joint pattern of the damped diagonal"
    damped_entries::Vector{Int}
    damping::Float64
end

"""
    RigidDampedJacobian(model::RigidDampedRHS, jac, prototype) -> (jacobian, pattern)

The Jacobian of `model` from the model's own `jac` on its `prototype` pattern, and the
joint pattern it writes into.
"""
function RigidDampedJacobian(model::RigidDampedRHS, jac, prototype)
    damped = [vec(model.layout.part_vel); vec(model.layout.body_spin)]
    rows, columns, _ = findnz(prototype)
    pattern = sparse([rows; damped], [columns; damped],
                     zeros(length(rows) + length(damped)), size(prototype)...)
    jacobian = RigidDampedJacobian(jac, SparseMatrixCSC{Float64, Int}(prototype),
        [pattern_entry(pattern, row, column) for (row, column) in zip(rows, columns)],
        [pattern_entry(pattern, slot, slot) for slot in damped], model.damping)
    return jacobian, pattern
end

"""Position in `pattern.nzval` of the stored entry at `row`, `column`."""
function pattern_entry(pattern::SparseMatrixCSC, row, column)
    range = nzrange(pattern, column)
    return range[searchsortedfirst(view(pattern.rowval, range), row)]
end

function (jacobian::RigidDampedJacobian)(J, u, p, t)
    jacobian.jac(jacobian.model_jacobian, u, p, t)
    J.nzval .= 0.0
    @views J.nzval[jacobian.model_entries] .+= jacobian.model_jacobian.nzval
    @views J.nzval[jacobian.damped_entries] .-= jacobian.damping
    return nothing
end

"""
    set_tube_damping!(state, values) -> state

Set the Rayleigh β of each of `state`'s tubes to its entry in `values`.
"""
function set_tube_damping!(state, values)
    foreach((tube, value) -> tube.model.damping = value, state.tubes, values)
    return state
end

"""
    set_damping!(state, point_values, body_values, tube_values) -> state

Write per point its `(body_frame, world_frame)` damping, per body its `(body_frame,
world_frame, angular)` damping and per tube its Rayleigh β into `state`.
"""
function set_damping!(state, point_values, body_values, tube_values)
    for (idx, (body_frame, world_frame)) in enumerate(point_values)
        set_body_frame_damping(state.points, body_frame, idx:idx)
        set_world_frame_damping(state.points, world_frame, idx:idx)
    end
    for (idx, (body_frame, world_frame, angular)) in enumerate(body_values)
        set_body_frame_damping(state.bodies, body_frame, idx:idx)
        set_world_frame_damping(state.bodies, world_frame, idx:idx)
        set_angular_damping(state.bodies, angular, idx:idx)
    end
    return set_tube_damping!(state, tube_values)
end

"""
    without_damping(f, sam)

Call `f()` with every artificial damping term zero — body-frame, world-frame and
angular damping on points and bodies, and the tubes' Rayleigh β — and put them back
afterwards.
"""
function without_damping(f, sam)
    state = sam.state
    points = [(copy(point.body_frame_damping), copy(point.world_frame_damping))
              for point in state.points]
    bodies = [(copy(body.body_frame_damping), copy(body.world_frame_damping),
               copy(body.angular_damping)) for body in state.bodies]
    tubes = [tube.model.damping for tube in state.tubes]
    set_damping!(state, [zero.(values) for values in points],
                 [zero.(values) for values in bodies], zero(tubes))
    sync_params!(sam.prob.param_sync, sam.integrator, state)
    try
        return f()
    finally
        set_damping!(state, points, bodies, tubes)
        sync_params!(sam.prob.param_sync, sam.integrator, state)
    end
end

"""Indices of the points on a tether of `state`, the winch points among them."""
function tether_point_idxs(state)
    segment_ends = (state.segments[idx].point_idxs
                    for tether in state.tethers for idx in tether.segment_idxs)
    return Set{Int}(Iterators.flatten(segment_ends))
end

"""
    mid_span_chord(state, wing) -> Float64

The chord of `wing` at its middle station [m]: between that station's leading- and
trailing-edge points where it has both, else the chord it stores.
"""
function mid_span_chord(state, wing)
    station = state.stations[wing.station_idxs[(end + 1) ÷ 2]]
    length(station.point_idxs) >= 2 || return norm(station.chord)
    leading_edge = state.points[first(station.point_idxs)]
    trailing_edge = state.points[last(station.point_idxs)]
    return norm(trailing_edge.pos_w .- leading_edge.pos_w)
end

"""
    residual_at(problem, u) -> NamedTuple

The undamped acceleration `u` leaves beyond the rigid rotation fitted to it, see
[`quasi_steady_residual`](@ref).
"""
function residual_at(problem::SteadyProblem, u)
    (; sam, layout, base) = problem
    state = sam.state
    spin = fitted_spin(problem, u)
    rigid = set_rigid_velocities!(copy(u), layout, spin, base)
    du = similar(u)
    without_damping(() -> sam.integrator.f(du, rigid, sam.integrator.p, 0.0), sam)
    read_state!(problem, u)
    tether = tether_point_idxs(state)
    kite = [[!(idx in tether) for idx in layout.points]; trues(length(layout.bodies))]
    wing = state.wings[1]
    force, moment = zeros(3), zeros(3)
    for (on_kite, part_mass, pos, vel) in zip(kite, problem.mass,
            eachcol(layout.part_pos), eachcol(layout.part_vel))
        on_kite || continue
        load = part_mass .* (du[vel] .- cross(spin, cross(spin, rigid[pos] .- base)))
        force .+= load
        moment .+= cross(rigid[pos] .- wing.pos_w, load)
    end
    for (idx, quat, body_spin) in zip(layout.bodies, eachcol(layout.body_quat),
                                      eachcol(layout.body_spin))
        moment .+= quaternion_to_rotation_matrix(rigid[quat]) *
                   (state.bodies[idx].inertia_principal .* du[body_spin])
    end
    aero = norm(wing.aero_force_b)
    chord = mid_span_chord(state, wing)
    return (; spin, force = wing.R_b_to_w' * force, moment = wing.R_b_to_w' * moment,
            force_ratio = norm(force) / aero,
            moment_ratio = norm(moment) / (aero * chord))
end

"""
    quasi_steady_residual(sam) -> NamedTuple

The undamped acceleration the integrator's state leaves beyond the rigid rotation about
the first transform's base fitted to it, on the last aero solve's forces and with every
artificial damping term off; writes that state into `sam.state`. Returns `(; spin, force,
moment, force_ratio, moment_ratio)`: `force` [N] and `moment` [N·m, about the wing origin,
each body's spin torque included] summed over every point and body not on a tether, on the
first wing's axes, and their sizes over the aero force and that force times the mid-span
chord. Winch, pulley and twist motion is not judged.
"""
quasi_steady_residual(sam::SymbolicAWEModel) =
    residual_at(SteadyProblem(sam), copy(sam.integrator.u))

"""The larger of a residual's two ratios."""
worst_ratio(residual) = max(residual.force_ratio, residual.moment_ratio)

"""
    find_quasi_steady!(sam; dt=0.01, max_steps=1000, check_every=5, tol=0.01, damping=5.0,
                       tube_damping_scale=30.0, max_time=600, prn=true)

Find the kite's quasi-steady state at the pose its transforms hold: every point and body
moving with one rigid rotation about the first transform's base. Flies it in steps of `dt`
[s] under `damping` [1/s] against that rotation and the tubes' Rayleigh β at
`tube_damping_scale` times its own, turning it back onto the pose after each, until
[`quasi_steady_residual`](@ref), read every `check_every` steps, has both ratios within
`tol`. Call after `init!`, on the `KernelBackend`, with the winches held: winch, pulley and
twist motion is neither damped nor judged.

Leaves the model ready to step on from the best state judged, on the damping it had before.
Returns `(; converged, status, steps, force, moment, force_ratio, moment_ratio, spin, speed,
turn_rate, wall_time)`; `status` is `:converged`, `:max_steps`, `:max_time`, `:aero_failed`
or `:integration_failed`, the last two with the residual of the best state before them.
"""
function find_quasi_steady!(sam::SymbolicAWEModel; dt = 0.01, max_steps = 1000,
                            check_every = 5, tol = 0.01, damping = 5.0,
                            tube_damping_scale = 30.0, max_time = 600, prn::Bool = true)
    start = time()
    problem = SteadyProblem(sam)
    tube_values = [tube.model.damping for tube in sam.state.tubes]
    result = without_damping(sam) do
        set_tube_damping!(sam.state, tube_values .* tube_damping_scale)
        sync_params!(sam.prob.param_sync, sam.integrator, sam.state)
        relax!(problem; dt, max_steps, check_every, tol, damping,
               deadline = start + max_time, prn)
    end
    refresh_aero!(problem, result.u; cold_start = true) ||
        @warn "find_quasi_steady!: the aero solve fails at the best state; the model " *
              "keeps the aerodynamics of the last solve" result.status
    run_time = sam.integrator.t
    OrdinaryDiffEqCore.reinit!(sam.integrator, result.u; t0 = run_time, tf = run_time + dt,
                               reinit_dae = false)
    (; spin, force, moment, force_ratio, moment_ratio) = result.residual
    radial = sam.state.wings[1].pos_w .- problem.base
    return (; converged = result.status === :converged, result.status, result.steps,
            force, moment, force_ratio, moment_ratio, spin,
            speed = norm(cross(spin, radial)),
            turn_rate = dot(spin, radial) / norm(radial), wall_time = time() - start)
end

"""
    relax!(problem; dt, max_steps, check_every, tol, damping, deadline, prn)

The flight loop of [`find_quasi_steady!`](@ref) on whatever damping the model carries:
returns `(; status, steps, u, residual)` with the best state it judged.
"""
function relax!(problem::SteadyProblem; dt, max_steps, check_every, tol, damping,
                deadline, prn)
    (; sam, layout) = problem
    (; f, p) = sam.integrator
    u = hold_pose!(problem, copy(sam.integrator.u))
    refresh_aero!(problem, u) || return (; status = :aero_failed, steps = 0, u,
                                         residual = residual_at(problem, u))
    model = RigidDampedRHS(f.f, layout, problem.base, fitted_spin(problem, u), damping)
    jacobian, pattern = RigidDampedJacobian(model, f.jac, f.jac_prototype)
    equations = ODEFunction(model; mass_matrix = f.mass_matrix, jac = jacobian,
                            jac_prototype = pattern)
    flight = init(ODEProblem(equations, u, (0.0, dt), p), sam.integrator.alg;
                  abstol = sam.set.abs_tol, reltol = sam.set.rel_tol, save_on = false)
    best = (; u = copy(u), residual = residual_at(problem, u))
    status, steps = :max_steps, 0
    for step in 1:max_steps
        time() > deadline && (status = :max_time; break)
        OrdinaryDiffEqCore.step!(flight, dt, true)
        successful_retcode(flight.sol) || (status = :integration_failed; break)
        u = hold_pose!(problem, copy(flight.u))
        model.spin .= fitted_spin(problem, u)
        refresh_aero!(problem, u) || (status = :aero_failed; break)
        OrdinaryDiffEqCore.reinit!(flight, u; t0 = 0.0, tf = dt)
        steps = step
        step % check_every == 0 || continue
        residual = residual_at(problem, u)
        worst_ratio(residual) < worst_ratio(best.residual) &&
            (best = (; u = copy(u), residual))
        prn && @info "find_quasi_steady!" step residual.force_ratio residual.moment_ratio
        worst_ratio(residual) <= tol && (status = :converged; break)
    end
    return (; status, steps, best.u, best.residual)
end

"""
    quasi_steady_log(sam; name="quasi_steady") -> SysLog

A `SysLog` of one frame, the state `sam` stands in as [`find_quasi_steady!`](@ref) leaves
it, holding the structure document of `sam` under the metadata key `topology`.
"""
function quasi_steady_log(sam::SymbolicAWEModel; name = "quasi_steady")
    logger = Logger(sam, 1)
    log!(logger, SysState(sam))
    metadata = Dict("topology" => KiteGeometry.to_json(sam.sys_struct.definition))
    return KiteUtils.sys_log(logger, name; metadata)
end
