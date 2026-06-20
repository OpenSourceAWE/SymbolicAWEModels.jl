# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Main system creation entry point

"""
    create_sys!(s::SymbolicAWEModel, system::SystemStructure; prn=true)

Create the full `ModelingToolkit.ODESystem` for the AWE model.

This is the main top-level function that orchestrates the generation of the entire
set of differential-algebraic equations (DAEs). It calls specialized sub-functions
to build the equations for each part of the system (forces, wing dynamics, scalar
kinematics, linearized aerodynamics) and assembles them into a single `System`.

# Arguments
- `s::SymbolicAWEModel`: The main model object to be populated with the system.
- `system::SystemStructure`: The physical structure definition.
- `prn::Bool=true`: If true, print progress information during system creation.

# Returns
- `set_values`: The symbolic variable representing the control inputs (winch torques).
"""
function create_sys!(s::SymbolicAWEModel, system::SystemStructure;
                     prn=true, tunable_params::Bool=false)
    eqs = Equation[]
    defaults = Pair{Num, Any}[]

    (; points, twist_surfaces, segments, pulleys, tethers, winches, wings,
       rigid_bodies, elastic_joints) = system

    validate_twist_surface_modes(twist_surfaces, wings)

    # Per-mode structural validation (e.g. VSM particle point↔panel mapping)
    for wing in wings
        validate_aero_structure(wing.aero, wing, points; prn)
    end

    # Flattened-parameter registry + build-time `params` view (see flat_params.jl).
    param_registry = ParamRegistry(system)
    params = ParamView(param_registry)

    # Initial-condition registry + build-time `initial` view (initial_conditions.jl).
    initial_registry = InitialRegistry(system)
    initial = InitialView(initial_registry)

    if tunable_params
        @parameters fix_wing = false
    else
        @parameters (fix_wing = false), [tunable = false]
    end
    @variables begin
        # Control inputs
        set_values(t)[eachindex(winches)] = zeros(length(winches))
        # Wing body frame output
        wing_pos(t)[1:3, eachindex(wings)]
        wing_vel(t)[1:3, eachindex(wings)]
        wing_acc(t)[1:3, eachindex(wings)]
        ω_b(t)[1:3, eachindex(wings)]
        α_b(t)[1:3, eachindex(wings)]
        # Wing principal frame ODE state
        com_w(t)[1:3, eachindex(wings)]
        com_vel(t)[1:3, eachindex(wings)]
        com_acc(t)[1:3, eachindex(wings)]
        Q_p_to_w(t)[1:4, eachindex(wings)]
        ω_p(t)[1:3, eachindex(wings)]
        α_p(t)[1:3, eachindex(wings)]
        # Rotation matrices
        R_b_to_w(t)[1:3, 1:3, eachindex(wings)]
        R_p_to_w(t)[1:3, 1:3, eachindex(wings)]
        R_v_to_w(t)[1:3, 1:3, eachindex(wings)]
        # Aerodynamic forces and moments
        aero_force_b(t)[1:3, eachindex(wings)]
        aero_moment_b(t)[1:3, eachindex(wings)]
        twist_surface_aero_moment(t)[eachindex(twist_surfaces)]
        # Wing deformation states
        twist_angle(t)[eachindex(twist_surfaces)]
        twist_ω(t)[eachindex(twist_surfaces)]
        # Wind and apparent velocity
        wind_vec_gnd(t)[1:3]
        va_wing_b(t)[1:3, eachindex(wings)]
        # Standalone rigid body state/output
        body_pos_w(t)[1:3, eachindex(rigid_bodies)]
        body_vel_w(t)[1:3, eachindex(rigid_bodies)]
        body_acc_w(t)[1:3, eachindex(rigid_bodies)]
        body_ω_b(t)[1:3, eachindex(rigid_bodies)]
        body_α_b(t)[1:3, eachindex(rigid_bodies)]
        body_com_w(t)[1:3, eachindex(rigid_bodies)]
        body_com_vel(t)[1:3, eachindex(rigid_bodies)]
        body_com_acc(t)[1:3, eachindex(rigid_bodies)]
        body_Q_p_to_w(t)[1:4, eachindex(rigid_bodies)]
        body_Q_b_to_w(t)[1:4, eachindex(rigid_bodies)]
        body_ω_p(t)[1:3, eachindex(rigid_bodies)]
        body_α_p(t)[1:3, eachindex(rigid_bodies)]
        body_moment_p(t)[1:3, eachindex(rigid_bodies)]
        body_Q_p_vel(t)[1:4, eachindex(rigid_bodies)]
        body_R_b_to_w(t)[1:3, 1:3, eachindex(rigid_bodies)]
        body_R_p_to_w(t)[1:3, 1:3, eachindex(rigid_bodies)]
    end
    R_b_to_w = collect(R_b_to_w)
    R_p_to_w = collect(R_p_to_w)
    R_v_to_w = collect(R_v_to_w)
    body_R_b_to_w = collect(body_R_b_to_w)
    body_R_p_to_w = collect(body_R_p_to_w)

    # Rigid body load accumulators (filled by joint_eqs!, read by body_eqs!).
    body_force = zeros(Num, 3, length(rigid_bodies))
    body_moment = zeros(Num, 3, length(rigid_bodies))

    # ==================== INLINED FORCE_EQS! CONTENT ==================== #
    # The following variables and component calls were previously in force_eqs!

    # Declare twist_surface geometry symbolic variables
    if length(twist_surfaces) > 0
        @variables begin
            twist_surface_y_airf(t)[1:3, eachindex(twist_surfaces)]
            twist_surface_chord(t)[1:3, eachindex(twist_surfaces)]
            twist_surface_le_pos(t)[1:3, eachindex(twist_surfaces)]
        end
    else
        twist_surface_y_airf = nothing
        twist_surface_chord = nothing
        twist_surface_le_pos = nothing
    end

    # Aggregate forces and moments from tethers onto the wing's center of mass
    tether_wing_force = zeros(Num, 3, length(wings))
    tether_wing_moment = zeros(Num, 3, length(wings))

    # Check if we have any PARTICLE_DYNAMICS wings (need aero force per point)
    has_particle_dynamics_wings = any(wing.dynamics_type === PARTICLE_DYNAMICS for wing in wings)

    @variables begin
        # Point states
        pos(t)[1:3, eachindex(points)]
        vel(t)[1:3, eachindex(points)]
        acc(t)[1:3, eachindex(points)]
        # Point forces and geometry
        point_force(t)[1:3, eachindex(points)]
        point_drag_force(t)[1:3, eachindex(points)]
        total_drag(t)[1:3, eachindex(points)]
        spring_sum_force(t)[1:3, eachindex(points)]  # Accumulated spring/drag forces
        disturb_force(t)[1:3, eachindex(points)]
        tether_r(t)[1:3, eachindex(points)]
        point_mass(t)[eachindex(points)]
        chord_b(t)[1:3, eachindex(points)]
        fixed_pos(t)[1:3, eachindex(points)]
        normal(t)[1:3, eachindex(points)]
        pos_b(t)[1:3, eachindex(points)]
        fix_point_sphere(t)[eachindex(points)]
        fix_static(t)[eachindex(points)]
        # Point damping (per-axis)
        body_frame_damping(t)[1:3, eachindex(points)]
        world_frame_damping(t)[1:3, eachindex(points)]
        # Segment forces and rest length
        spring_force_vec(t)[1:3, eachindex(segments)]
        drag_force(t)[1:3, eachindex(segments)]
        l0(t)[eachindex(segments)]
    end

    # Per-point variables for all points
    @variables begin
        va_point_b(t)[1:3, eachindex(points)]
        va_point_w(t)[1:3, eachindex(points)]
        wind_at_point(t)[1:3, eachindex(points)]
        height(t)[eachindex(points)]
    end

    # PARTICLE_DYNAMICS-specific variables: per-point aero forces in body frame
    if has_particle_dynamics_wings
        @variables begin
            aero_force_point_b(t)[1:3, eachindex(points)]
        end
    else
        aero_force_point_b = nothing
    end

    # Pulley and tether length state variables
    @variables begin
        pulley_len(t)[eachindex(pulleys)]
        pulley_vel(t)[eachindex(pulleys)]
        tether_len(t)[eachindex(tethers)]
        winch_vel(t)[eachindex(winches)]
        winch_acc(t)[eachindex(winches)]
        winch_force_vec(t)[1:3, eachindex(winches)]
        winch_force(t)[eachindex(winches)]
        winch_friction(t)[eachindex(winches)]
    end

    # ==================== CALL COMPONENT FUNCTIONS ==================== #

    # 1. Point equations (generates point dynamics, modifies tether_wing_force/moment in-place)
    eqs, defaults = point_eqs!(
        s, eqs, defaults, points, segments, twist_surfaces, wings, params, initial;
        R_b_to_w, com_w,
        wing_vel, wind_vec_gnd, twist_angle,
        pos, vel, acc, point_force, point_mass, spring_force_vec, drag_force, l0,
        spring_sum_force, point_drag_force, total_drag,
        disturb_force, tether_r, chord_b, fixed_pos, normal, pos_b,
        fix_point_sphere, fix_static, body_frame_damping, world_frame_damping,
        va_point_b, va_point_w, wind_at_point, height,
        aero_force_point_b,
        twist_surface_y_airf, tether_wing_force, tether_wing_moment
    )

    # 2. TwistSurface equations (deformable wing sections with twist dynamics)
    eqs, defaults = twist_surface_eqs!(
        eqs, defaults, twist_surfaces, wings, params, initial;
        R_b_to_w, fix_wing, twist_angle, twist_ω, twist_surface_aero_moment,
        point_force, tether_wing_moment, twist_surface_y_airf, twist_surface_chord, twist_surface_le_pos
    )

    # 3. Segment equations (spring-damper forces, returns len and spring_force)
    eqs, len, spring_force = segment_eqs!(
        s, eqs, points, segments, pulleys, tethers, wings, params;
        pos, vel, wind_vec_gnd, spring_force_vec, drag_force, l0,
        pulley_len, tether_len
    )

    # 4. Pulley equations (rope distribution)
    eqs, defaults = pulley_eqs!(
        eqs, defaults, pulleys, segments, params, initial;
        spring_force, pulley_len, pulley_vel
    )

    # 5. Winch equations (motor dynamics, tether reeling)
    eqs, defaults, winch_subsystems = winch_eqs!(
        eqs, defaults, winches, tethers, segments, points,
        system, params, initial;
        spring_force_vec, set_values, tether_len,
        winch_vel, winch_acc, winch_force_vec, winch_force,
        winch_friction
    )

    # 6. Tether equations (stretched length, average force)
    eqs = tether_eqs!(eqs, tethers; len, spring_force)

    # ==================== END INLINED FORCE_EQS! CONTENT ==================== #

    # Build aerodynamic equations: each wing's aero component (including
    # flat-plate) is wired in winch-style and returned as a subsystem.
    eqs, aero_subsystems = aero_eqs!(
        s, eqs, params;
        aero_force_b, aero_moment_b, twist_surface_aero_moment,
        twist_angle, twist_ω, va_wing_b, wing_pos, ω_b, R_b_to_w,
        pos, vel, va_point_b, height, aero_force_point_b
    )

    # Build wing rigid body dynamics equations
    eqs, defaults = wing_eqs!(
        s, eqs, defaults, params, initial;
        tether_wing_force, tether_wing_moment,
        aero_force_b, aero_moment_b,
        ω_b, α_b, R_b_to_w, R_p_to_w,
        wing_pos, wing_vel, wing_acc,
        com_w, com_vel, com_acc, Q_p_to_w, ω_p, α_p,
        fix_wing, pos, vel, acc
    )

    # Elastic joints: accumulate restoring wrenches into body_force/body_moment
    # (must precede body_eqs!, which reads them).
    eqs = joint_eqs!(
        eqs, elastic_joints, params;
        body_force, body_moment,
        body_com_w, body_pos_w, body_com_vel, body_ω_b, body_R_b_to_w,
    )

    # Build standalone rigid body dynamics equations
    eqs, defaults = body_eqs!(
        eqs, defaults, rigid_bodies, params, initial;
        body_force, body_moment,
        body_com_w, body_com_vel, body_com_acc, body_Q_p_to_w, body_ω_p, body_α_p,
        body_pos_w, body_vel_w, body_acc_w, body_ω_b, body_α_b, body_Q_b_to_w,
        body_R_b_to_w, body_R_p_to_w, body_moment_p, body_Q_p_vel,
    )

    # Build scalar kinematic and apparent wind equations
    eqs = scalar_eqs!(
        s, eqs, params;
        R_b_to_w, wind_vec_gnd, va_wing_b, wing_pos, wing_vel,
        wing_acc, twist_angle, ω_b, α_b, R_v_to_w, pos
    )

    # Debug: Find which equation fails to scalarize
    for (i, eq) in enumerate(eqs)
        try
            Symbolics.scalarize(eq)
        catch e
            println("Failed to scalarize equation index: $i")
            println("Eq: ", eqs[i])
            rethrow(e)
        end
    end
    
    eqs = Symbolics.scalarize.(reduce(vcat, Symbolics.scalarize.(eqs)))

    # Debug: Look for any remaining slice references after scalarization
    for (i, eq) in enumerate(eqs)
        eq_str = string(eq)
        if occursin("Colon()", eq_str)
            @warn "Equation $i contains Colon() after scalarization: $eq"
        end
    end

    all_subsystems = [winch_subsystems; aero_subsystems]
    time = @elapsed begin
        if isempty(all_subsystems)
            @named sys = System(eqs, t)
        else
            @named sys = System(eqs, t; systems = all_subsystems)
        end
    end
    prn && println("\tCreated System in $time seconds.")

    # set_values is seeded from the struct at init and set every step; no default.

    # Debug: Check defaults for slice references
    for (i, d) in enumerate(defaults)
        d_str = string(d)
        if occursin("Colon()", d_str)
            @warn "Default $i contains Colon(): $d"
        end
    end

    s.defaults = defaults
    s.full_sys = sys
    s.param_registry = param_registry
    s.initial_registry = initial_registry
    return set_values
end
