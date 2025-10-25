```@meta
CurrentModule = SymbolicAWEModels
```

# State Management

SymbolicAWEModels provides functions for managing and converting system state between different representations. The primary type for state representation is `SysState`, which provides a high-level view of the system's current condition.

## Overview

State management in SymbolicAWEModels involves:
- Converting between low-level integrator state and high-level `SysState`
- Accessing state information during simulation
- Setting initial conditions
- Interfacing with other packages in the Julia Kite Power Tools ecosystem

## SysState Type

The `SysState` type is defined in KiteUtils.jl and provides a convenient representation of the system state with fields for:
- Kite position and orientation
- Velocity information
- Tether lengths and forces
- Wing-specific parameters
- Time stamps

For complete documentation of the `SysState` type, see [KiteUtils.jl Documentation](https://opensourceawe.github.io/KiteUtils.jl/stable/).

## State Conversion Functions

### Creating a SysState from SymbolicAWEModel

```julia
# Create empty SysState compatible with the model
sys_state = SysState(sam)

# Update it with current integrator state
update_sys_state!(sys_state, sam)
```

### Setting Model State from SysState

```julia
# Set the integrator state from a SysState object
update_from_sysstate!(sam, sys_state)
```

## Usage Examples

### Example 1: Monitoring State During Simulation

```julia
using SymbolicAWEModels

set = Settings("system.yaml")
sam = SymbolicAWEModel(set, "ram")
init!(sam)

# Create SysState for state access
sys_state = SysState(sam)

# Simulate and monitor state
for i in 1:100
    next_step!(sam)
    update_sys_state!(sys_state, sam)

    # Access state information
    println("Time: $(sys_state.time)")
    println("Position: $(sys_state.X)")
    println("Heading: $(sys_state.heading)")
end
```

### Example 2: Setting Initial Conditions

```julia
# Load or create a SysState with desired initial conditions
sys_state = SysState(sam)
sys_state.X = [100.0, 0.0, 150.0]  # Set initial position
sys_state.heading = deg2rad(45)     # Set initial heading

# Apply to model
update_from_sysstate!(sam, sys_state)

# Now simulate from these initial conditions
(log, _) = sim_oscillate!(sam)
```

### Example 3: State-Based Control

```julia
# Create SysState for feedback
sys_state = SysState(sam)

for i in 1:1000
    # Read current state
    update_sys_state!(sys_state, sam)

    # Compute control based on state
    heading_error = target_heading - sys_state.heading
    control_input = heading_error * gain

    # Apply control
    next_step!(sam; set_values=[control_input])
end
```

## API Reference

```@docs
SysState
update_sys_state!
update_from_sysstate!
```

## Accessing State via Getter Functions

In addition to `SysState`, you can access specific state information using getter functions:

```julia
# Get wing state directly from integrator
wing_state = sam.prob.get_wing_state(sam.integrator)

# Get winch state
winch_state = sam.prob.get_winch_state(sam.integrator)
```

See the [Simulation Functions](exported_functions.md) page for more getter functions like `winch_force`, `unstretched_length`, and `tether_length`.

## Performance Considerations

- Creating a `SysState` allocates memory; reuse the same object when possible
- `update_sys_state!` modifies the existing `SysState` in-place for efficiency
- Direct getter functions may be faster for accessing specific values repeatedly

## See Also

- [`SymbolicAWEModel`](@ref) - Main simulation type
- [Simulation Functions](exported_functions.md) - High-level simulation interface and state accessor functions
- [KiteUtils.jl](https://opensourceawe.github.io/KiteUtils.jl/dev/) - SysState type documentation
