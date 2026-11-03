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
    return RigidLayout(points, state_slots(model, point_instances, :pos),
        state_slots(model, point_instances, :vel), bodies,
        state_slots(model, body_instances, :com_w),
        state_slots(model, body_instances, :com_vel),
        state_slots(model, body_instances, :Q),
        state_slots(model, body_instances, :omega_p))
end

"""Position slots of every state point, then of every body's centre of mass."""
part_positions(layout::RigidLayout) = [layout.point_pos layout.body_com]

"""Velocity slots in the order of [`part_positions`](@ref)."""
part_velocities(layout::RigidLayout) = [layout.point_vel layout.body_vel]

"""Hamilton product of two scalar-first quaternions."""
quat_product(q, p) = [q[1] * p[1] - dot(q[2:4], p[2:4]);
                      q[1] .* p[2:4] .+ p[1] .* q[2:4] .+ cross(q[2:4], p[2:4])]

"""
    set_rigid_velocities!(u, layout, spin, base) -> u

Put every point and body of `u` on the rigid rotation `spin` [rad/s] about `base`:
velocity `spin × (x − base)`, and the body spin `spin` resolved on its principal axes.
"""
function set_rigid_velocities!(u, layout::RigidLayout, spin, base)
    for (pos, vel) in zip(eachcol(part_positions(layout)), eachcol(part_velocities(layout)))
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
    for (pos, vel) in zip(eachcol(part_positions(layout)), eachcol(part_velocities(layout)))
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
    PassProblem(sam)

One quasi-steady solve of `sam`, about the base of its first transform.
"""
struct PassProblem{S}
    "the model whose integrator carries the state"
    sam::S
    "where the state of each point and body sits in `u`"
    layout::RigidLayout
    "ground-station position the rigid rotation turns about [m]"
    base::Vector{Float64}
    "mass of each part in the order of [`part_positions`](@ref) [kg]"
    mass::Vector{Float64}
end

function PassProblem(sam::SymbolicAWEModel)
    (; points, bodies, transforms) = sam.state
    layout = RigidLayout(sam)
    base = Vector{Float64}(KiteGeometry.base_position(KiteGeometry.Pose(sam.state),
                                                      transforms, first(transforms)))
    mass = Float64[[points[idx].total_mass for idx in layout.points];
                   [bodies[idx].total_mass for idx in layout.bodies]]
    return PassProblem(sam, layout, base, mass)
end

"""
    pass_spin(problem, u) -> Vector{Float64}

The rigid rotation about the base [rad/s] that best fits the velocities of `u`,
each point and body weighed by its mass.
"""
function pass_spin(problem::PassProblem, u)
    (; layout, base, mass) = problem
    inertia, momentum = zeros(3, 3), zeros(3)
    for (part_mass, pos, vel) in zip(mass, eachcol(part_positions(layout)),
                                     eachcol(part_velocities(layout)))
        arm = u[pos] .- base
        inertia .+= part_mass .* (dot(arm, arm) .* I(3) .- arm * arm')
        momentum .+= part_mass .* cross(arm, u[vel])
    end
    return inertia \ momentum
end

"""Write the ODE state `u` into `sam.state`, leaving the integrator alone."""
function read_state!(problem::PassProblem, u)
    (; sam) = problem
    sam.prob.get_all_state((; u, p = sam.integrator.p, t = sam.integrator.t), sam.state)
    return nothing
end

"""
    hold_pose!(problem, u) -> u

Turn the whole of `u` rigidly about the base back onto the pose its transforms hold,
as `reposition!` places it, velocities and orientations turned with it.
"""
function hold_pose!(problem::PassProblem, u)
    (; sam, layout, base) = problem
    points = sam.state.points
    read_state!(problem, u)
    before = reduce(hcat, [points[idx].pos_w .- base for idx in layout.points])
    reposition!(sam.state)
    after = reduce(hcat, [points[idx].pos_w .- base for idx in layout.points])
    return rotate_state!(u, layout, best_rotation(before, after), base)
end

"""
    refresh_aero!(problem::PassProblem, u; cold_start=false) -> Bool

Write `u` into the model's integrator and state, re-solve the aerodynamics there, from
scratch with `cold_start`, and push the result into the parameters. An aero solve that
gives up returns `false`.
"""
function refresh_aero!(problem::PassProblem, u; cold_start = false)
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
    PassDampedRHS

The model's right-hand side with every point and body damped against the rigid
rotation `spin` about `base`: `damping` [1/s] on the velocity left beyond
`spin × (x − base)`, and on the body spin left beyond `spin`.
"""
struct PassDampedRHS{F}
    rhs::F
    layout::RigidLayout
    base::Vector{Float64}
    spin::Vector{Float64}
    damping::Float64
end

function (model::PassDampedRHS)(du, u, p, t)
    (; layout, base, spin, damping) = model
    model.rhs(du, u, p, t)
    for (pos, vel) in zip(eachcol(part_positions(layout)), eachcol(part_velocities(layout)))
        du[vel] .-= damping .* (u[vel] .- cross(spin, u[pos] .- base))
    end
    for (quat, body_spin) in zip(eachcol(layout.body_quat), eachcol(layout.body_spin))
        orientation = quaternion_to_rotation_matrix(u[quat])
        du[body_spin] .-= damping .* (u[body_spin] .- orientation' * spin)
    end
    return nothing
end

"""
    PassDampedJacobian

The model's analytic state Jacobian with the diagonal of [`PassDampedRHS`](@ref)'s
damping added, on a pattern that holds both.
"""
struct PassDampedJacobian{J}
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
    PassDampedJacobian(model::PassDampedRHS, jac, prototype) -> (jacobian, pattern)

The Jacobian of `model` from the model's own `jac` on its `prototype` pattern, and the
joint pattern it writes into.
"""
function PassDampedJacobian(model::PassDampedRHS, jac, prototype)
    damped = [vec(part_velocities(model.layout)); vec(model.layout.body_spin)]
    rows, columns, _ = findnz(prototype)
    pattern = sparse([rows; damped], [columns; damped],
                     zeros(length(rows) + length(damped)), size(prototype)...)
    jacobian = PassDampedJacobian(jac, SparseMatrixCSC{Float64, Int}(prototype),
        [pattern_entry(pattern, row, column) for (row, column) in zip(rows, columns)],
        [pattern_entry(pattern, slot, slot) for slot in damped], model.damping)
    return jacobian, pattern
end

"""Position in `pattern.nzval` of the stored entry at `row`, `column`."""
function pattern_entry(pattern::SparseMatrixCSC, row, column)
    range = nzrange(pattern, column)
    return range[searchsortedfirst(view(pattern.rowval, range), row)]
end

function (jacobian::PassDampedJacobian)(J, u, p, t)
    jacobian.jac(jacobian.model_jacobian, u, p, t)
    J.nzval .= 0.0
    J.nzval[jacobian.model_entries] .+= jacobian.model_jacobian.nzval
    J.nzval[jacobian.damped_entries] .-= jacobian.damping
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
    le, te = state.points[first(station.point_idxs)], state.points[last(station.point_idxs)]
    return norm(te.pos_w .- le.pos_w)
end

"""
    pass_residual(problem, u) -> NamedTuple

The undamped acceleration `u` leaves beyond the rigid rotation fitted to it, see
[`quasi_steady_residual`](@ref).
"""
function pass_residual(problem::PassProblem, u)
    (; sam, layout, base) = problem
    state = sam.state
    spin = pass_spin(problem, u)
    rigid = set_rigid_velocities!(copy(u), layout, spin, base)
    du = similar(u)
    without_damping(() -> sam.integrator.f(du, rigid, sam.integrator.p, 0.0), sam)
    read_state!(problem, u)
    tether = tether_point_idxs(state)
    kite = [[!(idx in tether) for idx in layout.points]; trues(length(layout.bodies))]
    wing = state.wings[1]
    force, moment = zeros(3), zeros(3)
    for (on_kite, part_mass, pos, vel) in zip(kite, problem.mass,
            eachcol(part_positions(layout)), eachcol(part_velocities(layout)))
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

The undamped acceleration the state of `sam` leaves beyond the rigid rotation about its
first transform's base fitted to it, with every artificial damping term off:
`(; spin, force, moment, force_ratio, moment_ratio)`. `force` [N] and `moment` [N·m,
about the wing origin, each body's spin torque included] are summed over every point and
body not on a tether, on the axes of the first wing; the ratios divide their sizes by the
aerodynamic force and by that force times the mid-span chord.
"""
quasi_steady_residual(sam::SymbolicAWEModel) =
    pass_residual(PassProblem(sam), copy(sam.integrator.u))

"""The larger of a residual's two ratios."""
worst_ratio(residual) = max(residual.force_ratio, residual.moment_ratio)

"""
    find_quasi_steady!(sam; dt=0.05, max_steps=1000, check_every=5, tol=0.01, damping=20.0,
                       tube_damping_scale=30.0, max_time=600, prn=true)

Find the kite's quasi-steady state at the pose its transforms hold: the state in which
every point and body moves with one rigid rotation about the base of the first transform
and nothing accelerates beyond it. The elevation, azimuth and heading stay the
transforms'; the shape, the attitude and the rotation are free. Call after `init!`, on a
model built with the `KernelBackend`.

The kite is flown in steps of `dt` [s] with only the damping of [`PassDampedRHS`](@ref),
`damping` [1/s], and the tubes' Rayleigh β at `tube_damping_scale` times its own. After
each step it is turned back onto the pose, the rotation is refitted and the aerodynamics
re-solved, so they hold still within a step and a light kite needs a shorter `dt`. Every
`check_every` steps [`quasi_steady_residual`](@ref) judges it; it has converged when both
ratios are within `tol`.

Leaves the model ready to step on from the best state judged, on the damping it had
before, with the aerodynamics re-solved there where that solve converges, and returns
`(; converged, status, steps, force, moment, force_ratio, moment_ratio, spin, speed,
turn_rate, wall_time)`. `status` is `:converged`, `:max_steps` or `:max_time`, or
`:aero_failed` / `:integration_failed` where the flight left the model's range; the
residual is then that of the best state before it.
"""
function find_quasi_steady!(sam::SymbolicAWEModel; dt = 0.05, max_steps = 1000,
                            check_every = 5, tol = 0.01, damping = 20.0,
                            tube_damping_scale = 30.0, max_time = 600, prn::Bool = true)
    start = time()
    problem = PassProblem(sam)
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
function relax!(problem::PassProblem; dt, max_steps, check_every, tol, damping,
                deadline, prn)
    (; sam, layout) = problem
    (; f, p) = sam.integrator
    u = hold_pose!(problem, copy(sam.integrator.u))
    refresh_aero!(problem, u)
    model = PassDampedRHS(f.f, layout, problem.base, pass_spin(problem, u), damping)
    jacobian, pattern = PassDampedJacobian(model, f.jac, f.jac_prototype)
    equations = ODEFunction(model; mass_matrix = f.mass_matrix, jac = jacobian,
                            jac_prototype = pattern)
    flight = init(ODEProblem(equations, u, (0.0, dt), p), sam.integrator.alg;
                  abstol = sam.set.abs_tol, reltol = sam.set.rel_tol, save_on = false)
    best = (; u = copy(u), residual = pass_residual(problem, u))
    status, steps = :max_steps, 0
    for step in 1:max_steps
        time() > deadline && (status = :max_time; break)
        OrdinaryDiffEqCore.step!(flight, dt, true)
        successful_retcode(flight.sol) || (status = :integration_failed; break)
        u = hold_pose!(problem, copy(flight.u))
        model.spin .= pass_spin(problem, u)
        refresh_aero!(problem, u) || (status = :aero_failed; break)
        OrdinaryDiffEqCore.reinit!(flight, u; t0 = 0.0, tf = dt)
        steps = step
        step % check_every == 0 || continue
        residual = pass_residual(problem, u)
        worst_ratio(residual) < worst_ratio(best.residual) &&
            (best = (; u = copy(u), residual))
        prn && @info "find_quasi_steady!" step residual.force_ratio residual.moment_ratio
        worst_ratio(residual) <= tol && (status = :converged; break)
    end
    return (; status, steps, best.u, best.residual)
end

"""
    quasi_steady_log(sam; name="quasi_steady") -> SysLog

A `SysLog` of one frame: the state `sam` stands in, as [`find_quasi_steady!`](@ref)
leaves it.
"""
function quasi_steady_log(sam::SymbolicAWEModel; name = "quasi_steady")
    logger = Logger(sam, 1)
    log!(logger, SysState(sam))
    return KiteUtils.sys_log(logger, name)
end
