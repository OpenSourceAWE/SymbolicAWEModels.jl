# SPDX-FileCopyrightText: 2025 Bart van de Lint
#
# SPDX-License-Identifier: MPL-2.0

function vsm_affect!(integ::SciMLBase.Integrator)
end

function SciMLBase.__solve(prob::SciMLBase.SteadyStateProblem, alg::VSMDynamicSS,
        args...; abstol = 1e-8, reltol = 1e-6, odesolve_kwargs = (;),
        save_idxs = nothing, termination_condition = NonlinearSolveBase.NormTerminationMode(infnorm),
        kwargs...)
    tspan = __get_tspan(prob.u0, alg)

    f = if prob isa SteadyStateProblem
        prob.f
    elseif prob isa NonlinearProblem
        if isinplace(prob)
            (du, u, p, t) -> prob.f(du, u, p)
        else
            (u, p, t) -> prob.f(u, p)
        end
    end

    if isinplace(prob)
        du = similar(prob.u0)
        f(du, prob.u0, prob.p, first(tspan))
    else
        du = f(prob.u0, prob.p, first(tspan))
    end

    tc_cache = init(prob, termination_condition, du, prob.u0; abstol, reltol)
    abstol = NonlinearSolveBase.get_abstol(tc_cache)
    reltol = NonlinearSolveBase.get_reltol(tc_cache)

    function terminate_function(u, t, integrator)
        return tc_cache(get_du(integrator), integrator.u, integrator.uprev, t)
    end

    callback = TerminateSteadyState(abstol, reltol, terminate_function;
        wrap_test = Val(false))

    haskey(kwargs, :callback) && (callback = CallbackSet(callback, kwargs[:callback]))
    haskey(odesolve_kwargs, :callback) &&
        (callback = CallbackSet(callback, odesolve_kwargs[:callback]))

    # Construct and solve the ODEProblem
    odeprob = ODEProblem{isinplace(prob), true}(f, prob.u0, tspan, prob.p)
    odesol = solve(odeprob, alg.alg, args...; abstol, reltol, kwargs...,
        odesolve_kwargs..., callback, save_end = true)

    resid, u, retcode = __get_result_from_sol(termination_condition, tc_cache, odesol)

    if save_idxs !== nothing
        u = u[save_idxs]
        resid = resid[save_idxs]
    end

    return SciMLBase.build_solution(prob, DynamicSS(odesol.alg, alg.tspan), u, resid;
        retcode, original = odesol)
end
