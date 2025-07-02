
function linearize_vsm!(s::SymbolicAWEModel, integ=s.integrator)
    @unpack wings, y, x, jac = s.sys_struct
    if length(wings) > 0
        y .= s.get_vsm_y(integ)
        for wing in wings
            res = VortexStepMethod.linearize(
                s.vsm_solvers[wing.idx], 
                s.vsm_aeros[wing.idx], 
                y[wing.idx, :];
                va_idxs=1:3, 
                omega_idxs=4:6,
                theta_idxs=7:6+length(s.sys_struct.groups),
                moment_frac=s.sys_struct.groups[1].moment_frac
            )
            jac[wing.idx, :, :] .= res[1]
            x[wing.idx, :] .= res[2]
        end
        s.set_vsm(integ, [x, y, jac])
    end
    nothing
end

function linearize(s::SymbolicAWEModel; set_values=s.get_set_values(s.integrator))
    isnothing(s.lin_prob) && error("Run init_sim! with remake=true and lin_outputs=...")
    s.set_lin_vsm(s.lin_prob, s.get_vsm(s.integrator))
    s.set_lin_set_values(s.lin_prob, set_values)
    s.set_lin_unknowns(s.lin_prob, s.get_unknowns(s.integrator))
    return solve(s.lin_prob)
end

function simple_linearize(s::SymbolicAWEModel)
    init_sim!(s; remake=false, adaptive=true)
    find_steady_state!(s)
    simple_linearize!(s)
    return s.A, s.B, s.C, s.D
end

function set_measured!(s::SymbolicAWEModel, 
    heading, turn_rate,
    tether_length, tether_vel
)
    @unpack wings, winches = s.sys_struct
    wing = wings[1]

    # get variables from integrator
    distance = norm(wing.pos_w)
    R_t_w = calc_R_t_w(wing.elevation, wing.azimuth) # rotation of tether to world, similar to view rotation, but always pointing up
    
    # get wing_pos, rotate it by elevation and azimuth around the x and z axis
    wing.pos_w .= R_t_w * [0, 0, distance]
    # wing_vel from elevation_vel and azimuth_vel
    wing.vel_w .= R_t_w * [-wing.elevation_vel, wing.azimuth_vel, tether_vel[1]]
    # find quaternion orientation from heading, R_cad_body and R_t_w
    x = [cos(-heading), -sin(-heading), 0]
    y = [sin(-heading),  cos(-heading), 0]
    z = [0, 0, 1]
    wing.R_b_w .= R_t_w * s.vsm_wings[1].R_cad_body' * [x y z]
    # adjust the turn rates for observed turn rate
    wing.ω_b .= wing.R_b_w' * R_t_w * [0, 0, turn_rate]
    # directly set tether length
    for winch in winches
        winch.tether_length = tether_length[winch.idx]
        winch.tether_vel = tether_vel[winch.idx]
    end
    return nothing
end

function jacobian(f::Function, x::AbstractVector, ϵ::AbstractVector)
    n = length(x)
    fx = f(x)
    m = length(fx)
    J = zeros(m, n)
    for i in 1:n
        x_perturbed = copy(x)
        x_perturbed[i] += ϵ[i]
        J[:, i] = (f(x_perturbed) - fx) / ϵ[i]
    end
    return J
end

function simple_linearize!(s::SymbolicAWEModel; tstab=0.1)
    old_stab = s.get_stabilize(s.integrator)
    s.set_stabilize(s.integrator, true)
    integ = s.integrator
    init_sim!(s; reload=false)
    lin_x0 = s.get_lin_x(s.integrator)
    u0 = [winch.set_value for winch in s.sys_struct.winches]
    s.A .= 0.0
    s.B .= 0.0
    s.C .= 0.0
    s.D .= 0.0
    
    function f(x, u)
        heading = x[1]
        turn_rate = x[2]
        @show turn_rate
        tether_length = x[3:5]
        tether_vel = x[6:8]
        set_measured!(s, heading, turn_rate,
                      tether_length, tether_vel)
        s.set_psys(integ, s.sys_struct)
        @show s.sys_struct.wings[1].angular_vel
        s.set_set_values(s.integrator, u)
        OrdinaryDiffEqCore.reinit!(integ)
        OrdinaryDiffEqCore.step!(integ, tstab)
        @show s.get_lin_dx(integ)[1]
        return s.get_lin_dx(integ)
    end

    # yes it looks weird to step in an output function, but this is a steady state finder rather than output
    function h(x)
        heading = x[1]
        turn_rate = x[2]
        tether_length = x[3:5]
        tether_vel = x[6:8]
        set_measured!(s, heading, turn_rate,
                      tether_length, tether_vel)
        s.set_psys(integ, s.sys_struct)
        OrdinaryDiffEqCore.reinit!(integ)
        OrdinaryDiffEqCore.step!(integ, tstab)
        return s.get_lin_y(integ)
    end

    f_x(x) = f(x, u0)
    f_u(u) = f(lin_x0, u)

    # calculate jacobian
    ϵ_x = [0.01, 0.01, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
    ϵ_u = [1.0, 0.1, 0.1]
    s.A .= jacobian(f_x, lin_x0, ϵ_x)
    s.B .= jacobian(f_u, u0, ϵ_u)
    s.C .= jacobian(h,   lin_x0, ϵ_x)
    
    s.set_set_values(integ, u0)
    s.set_stabilize(integ, old_stab)
    nothing
end
