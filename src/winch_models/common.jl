# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Code shared by all winch models: dispatch interface + connector validation.

"""
    winch_component(model::AbstractWinchModel, sys_struct, winch_idx; name, params)

Build the winch motor subsystem for `sys_struct.winches[winch_idx]`, selected by
dispatch on the winch's `model`. The component is pure-algebraic at its connector
boundary (see [`validate_winch_component`](@ref)) but may declare arbitrary
internal `D(x) ~ …` states. Common drum parameters are read live via the flat
`params` view (`params.winches[winch_idx].gear_ratio` etc.); a model's own extra
fields go through `params.winches[winch_idx].model.your_field`. Add a method on a
custom [`AbstractWinchModel`](@ref) subtype to plug in your own dynamics.
"""
function winch_component end

"""
    is_builtin_winch(model::AbstractWinchModel) -> Bool

`true` for models shipped with this package. A `false` (custom) model forces a
cache rebuild ([`has_custom_component`](@ref)). Custom subtypes inherit the
`false` fallback.
"""
is_builtin_winch(::AbstractWinchModel) = false

"""
    validate_winch_component(subsys, winch)

Check that `subsys` (built by [`winch_component`](@ref)) satisfies the connector
contract.

Required connector variables:
- `vel` (input, drum-perimeter velocity [m/s])
- `len` (input, mean of connected tether lengths [m])
- `force` (input, summed tether tension magnitude [N])
- `set_value` (input, abstract setpoint; component fixes meaning)
- `brake` (input, brake in [0, 1])
- `acc` (output, drum-perimeter acceleration [m/s²])
- `friction` (output, friction torque [N·m])

Forbidden:
- Equations whose LHS is `D(vel)` or `D(len)` (those derivatives belong to the
  outer SymbolicAWEModels system). Internal `D(x) ~ …` for any other variable is
  allowed.
"""
function validate_winch_component(subsys, winch)
    required = (:vel, :len, :force, :set_value, :brake, :acc, :friction)
    required_str = join(required, ", ")
    for c in required
        hasproperty(subsys, c) || error(
            "Winch $(winch.name): component returned by `winch_component` " *
            "is missing required connector `$c`. Required connectors: " *
            "$required_str.")
    end
    for eq in ModelingToolkit.equations(subsys)
        lhs = ModelingToolkit.Symbolics.unwrap(eq.lhs)
        var_name = differential_inner_name(lhs)
        if var_name === :vel || var_name === :len
            error("Winch $(winch.name): component must not define " *
                  "`D($var_name) ~ …`; that derivative is owned by the outer " *
                  "system.")
        end
    end
    return nothing
end

function differential_inner_name(expr)
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
