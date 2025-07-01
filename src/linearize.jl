

function create_callback(s::SymbolicAWEModel; atol=0.1)
    function affect!(integ::ODEIntegrator)
        linearize_vsm!(s, integ)
        nothing
    end
    function condition(u, t, integ)
        norm(s.get_dy(integ)) > atol
    end
    cb = DiscreteCallback(condition, affect!; initializealg=OrdinaryDiffEqCore.NoInit())
    return cb
end

function linearize_vsm!(s::SymbolicAWEModel, integ=s.integrator)
    @unpack wings, y, x, jac = s.sys_struct
    if length(wings) > 0
        y .= s.get_y(integ)
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

