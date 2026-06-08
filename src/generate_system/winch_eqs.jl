# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

"""
    default_winch_component(sys_struct, winch_idx; name) -> ODESystem

Build the default winch motor component for
`sys_struct.winches[winch_idx]`.

The component is pure-algebraic at its connector boundary: it exposes
inputs (`vel`, `len`, `force`, `set_value`, `brake`) and outputs
(`acc`, `friction`) and contains no internal differential states. The
outer SymbolicAWEModels integrator owns `winch_vel` and `tether_len`;
`len` is wired to the mean of the connected tethers' lengths.

Custom winch components must expose the same seven connector variables
(see [`validate_winch_component`](@ref)) but may declare arbitrary
internal `D(x) ~ …` states. Researchers can plug in their own model by
passing a builder of the same signature to `Winch(...; model=...)`.
The default component declares `len` and ignores it.

# Arguments
- `sys_struct::SystemStructure`: Live system structure. Used as the
  default value for the component's `psys` parameter so that
  registered `get_winch_*` accessors read live struct fields.
- `winch_idx::Int`: Index of the winch being built.

# Keyword arguments
- `name::Symbol`: Subsystem name (required by `@named`-style usage).

# Default equations
```
ω_motor   = gear_ratio / drum_radius * vel
friction  = smooth_sign(ω_motor, friction_eps) * f_coulomb *
            drum_radius / gear_ratio
            + c_vf * ω_motor * drum_radius^2 / gear_ratio^2
tau_total = set_value + drum_radius / gear_ratio * force - friction
α_motor   = tau_total / inertia_total
acc       = ifelse(brake > 0.5, 0, drum_radius / gear_ratio * α_motor)
```

`set_value` is interpreted as motor torque [N·m].
"""
function default_winch_component(sys_struct::SystemStructure,
                                 winch_idx::Int; name)
    SST = typeof(sys_struct)
    @parameters (psys::SST = sys_struct), [tunable = false]
    @variables begin
        vel(t)
        len(t)
        force(t)
        set_value(t)
        brake(t)
        acc(t)
        friction(t)
        ω_motor(t)
        tau_total(t)
        α_motor(t)
    end

    gear_ratio    = get_winch_gear_ratio(psys, winch_idx)
    drum_radius   = get_winch_drum_radius(psys, winch_idx)
    f_coulomb     = get_winch_f_coulomb(psys, winch_idx)
    c_vf          = get_winch_c_vf(psys, winch_idx)
    inertia_total = get_winch_inertia_total(psys, winch_idx)
    friction_eps  = get_winch_friction_epsilon(psys, winch_idx)
    smooth_sign(x, eps) = x / sqrt(x * x + eps * eps)
    ratio = drum_radius / gear_ratio

    eqs = [
        ω_motor   ~ vel / ratio
        friction  ~ smooth_sign(ω_motor, friction_eps) * f_coulomb * ratio +
                    c_vf * ω_motor * ratio^2
        tau_total ~ set_value + ratio * force - friction
        α_motor   ~ tau_total / inertia_total
        acc       ~ ifelse(brake > 0.5, 0.0, ratio * α_motor)
    ]
    return System(eqs, t,
                  [vel, len, force, set_value, brake, acc, friction,
                   ω_motor, tau_total, α_motor],
                  [psys]; name)
end

"""
    validate_winch_component(subsys, winch)

Check that `subsys` (built by `winch.model(...)`) satisfies the
connector contract.

Required connector variables:
- `vel` (input, drum-perimeter velocity [m/s])
- `len` (input, mean of connected tether lengths [m])
- `force` (input, summed tether tension magnitude [N])
- `set_value` (input, abstract setpoint; component fixes meaning)
- `brake` (input, brake in [0, 1])
- `acc` (output, drum-perimeter acceleration [m/s²])
- `friction` (output, friction torque [N·m])

Forbidden:
- Equations whose LHS is `D(vel)` or `D(len)` (those derivatives
  belong to the outer SymbolicAWEModels system).

Internal `D(x) ~ …` equations for any other variable are allowed.
"""
function validate_winch_component(subsys, winch)
    required = (:vel, :len, :force, :set_value, :brake, :acc, :friction)
    required_str = join(required, ", ")
    for c in required
        hasproperty(subsys, c) || error(
            "Winch $(winch.name): component returned by `winch.model` " *
            "is missing required connector `$c`. Required connectors: " *
            "$required_str.")
    end
    for eq in ModelingToolkit.equations(subsys)
        lhs = ModelingToolkit.Symbolics.unwrap(eq.lhs)
        var_name = _differential_inner_name(lhs)
        if var_name === :vel || var_name === :len
            error("Winch $(winch.name): component must not define " *
                  "`D($var_name) ~ …`; that derivative is owned by the outer " *
                  "system.")
        end
    end
    return nothing
end

function _differential_inner_name(expr)
    try
        ModelingToolkit.iscall(expr) || return nothing
        mtk_operation = ModelingToolkit.operation(expr)
        mtk_operation isa ModelingToolkit.Differential || return nothing
        arg = ModelingToolkit.arguments(expr)[1]
        ModelingToolkit.iscall(arg) || return nothing
        inner = ModelingToolkit.operation(arg)
        return nameof(inner)
    catch
        return nothing
    end
end

"""
    winch_eqs!(eqs, defaults, winches, tethers, segments, points,
               sys_struct, psys;
               spring_force_vec, set_values, tether_len,
               winch_vel, winch_acc, winch_force_vec, winch_friction)

Generate equations for winch motor dynamics and per-tether length
state, and return the list of `ODESystem` subsystems to attach to the
parent system.

For each winch:
1. Sum spring force vectors of the segments meeting the winch point
   (sign-aware via segment `point_idxs`) → `winch_force_vec`.
2. Instantiate the user's winch component via `winch.model(...)`.
3. Validate the connector contract with
   [`validate_winch_component`](@ref).
4. Bind connectors to the parent variables (including
   `subsys.len ~ mean(tether_len[tether_idx] for tether_idx in winch.tether_idxs)`)
   and integrate `D(winch_vel) = ifelse(brake > 0.5, 0, winch_acc)`.
   When `winch.speed_controlled` is true, `winch_acc` is forced to 0
   (ignoring `subsys.acc`) so velocity is prescribed via `winch.vel`.

For each tether:
- With winch:    `D(tether_len) = ifelse(brake > 0.5, 0, winch_vel)`.
- Without winch: `D(tether_len) = 0`.
"""
function winch_eqs!(eqs, defaults, winches, tethers, segments, points,
                    sys_struct, psys;
                    spring_force_vec, set_values, tether_len,
                    winch_vel, winch_acc, winch_force_vec, winch_force,
                    winch_friction)
    tether_winch = Dict{Int, Int}()
    for winch in winches
        for tether_idx in winch.tether_idxs
            haskey(tether_winch, tether_idx) && error(
                "Tether $tether_idx is connected to winch " *
                "$(tether_winch[tether_idx]) and winch $(winch.idx). " *
                "Each tether can have at most one winch.")
            tether_winch[tether_idx] = winch.idx
        end
    end

    for tether in tethers
        if haskey(tether_winch, tether.idx)
            winch_idx = tether_winch[tether.idx]
            eqs = [eqs
                   D(tether_len[tether.idx]) ~
                       ifelse(get_brake(psys, winch_idx) > 0.5,
                              0, winch_vel[winch_idx])]
        else
            eqs = [eqs; D(tether_len[tether.idx]) ~ 0]
        end
        defaults = [defaults
                    tether_len[tether.idx] =>
                        get_tether_len(psys, tether.idx)]
    end

    winch_subsystems = Any[]
    for winch in winches
        isempty(winch.tether_idxs) &&
            error("Winch $(winch.name): no connected tethers; " *
                  "at least one is required.")
        winch_point_idx = winch.winch_point_idx
        (winch_point_idx > length(points)) &&
            error("Winch $(winch.name): point $winch_point_idx does not exist.")

        winch_seg_idxs = Set{Int}()
        for tether_idx in winch.tether_idxs
            union!(winch_seg_idxs, tethers[tether_idx].segment_idxs)
        end
        force_vec = zeros(Num, 3)
        for segment in segments
            segment.idx in winch_seg_idxs || continue
            if segment.point_idxs[1] == winch_point_idx
                force_vec .+= spring_force_vec[:, segment.idx]
            elseif segment.point_idxs[2] == winch_point_idx
                force_vec .-= spring_force_vec[:, segment.idx]
            end
        end

        subsys = winch.model(sys_struct, winch.idx;
                             name=Symbol("winch_$(winch.idx)"))
        validate_winch_component(subsys, winch)
        push!(winch_subsystems, subsys)

        brake_p = get_brake(psys, winch.idx)
        eqs = [eqs
               winch_force_vec[:, winch.idx] ~ force_vec
               winch_force[winch.idx] ~
                   smooth_norm(winch_force_vec[:, winch.idx])
               subsys.vel       ~ winch_vel[winch.idx]
               subsys.len       ~
                   sum(tether_len[tether_idx] for tether_idx in winch.tether_idxs) /
                   length(winch.tether_idxs)
               subsys.force     ~ winch_force[winch.idx]
               subsys.set_value ~ set_values[winch.idx]
               subsys.brake     ~ brake_p
               winch_acc[winch.idx]      ~
                   ifelse(get_speed_controlled(psys, winch.idx) == true,
                          0.0, subsys.acc)
               winch_friction[winch.idx] ~ subsys.friction
               D(winch_vel[winch.idx]) ~
                   ifelse(brake_p > 0.5, 0, winch_acc[winch.idx])]
        defaults = [defaults
                    winch_vel[winch.idx] =>
                        get_winch_vel(psys, winch.idx)]
    end
    return eqs, defaults, winch_subsystems
end
