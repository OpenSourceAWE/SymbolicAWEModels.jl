# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

"""
Steering and depower configuration for segment-based control.

Maps normalized steering [-1, 1] and depower [0, 1] inputs to
segment rest length offsets.
"""

# ==================== STEERING CONFIG ==================== #

"""
    mutable struct SteeringConfig

Configuration for mapping normalized steering and depower inputs
to segment rest length offsets.

Steering input `s ∈ [-1, 1]` differentially adjusts left/right
segment lengths. Depower input `d ∈ [0, 1]` lengthens either a
dedicated segment or the steering segments.

$(TYPEDFIELDS)
"""
mutable struct SteeringConfig
    "Index of the left steering segment."
    steer_left_idx::Int64
    "Index of the right steering segment."
    steer_right_idx::Int64
    "Max length change per side for steering [m]."
    steer_gain::SimFloat

    "Index of the depower segment (0 = use steering segments)."
    depower_idx::Int64
    "Max length change for depower [m]."
    depower_gain::SimFloat

    "Base rest length of left steering segment [m]."
    l0_steer_left::SimFloat
    "Base rest length of right steering segment [m]."
    l0_steer_right::SimFloat
    "Base rest length of depower segment [m] (only when depower_idx > 0)."
    l0_depower::SimFloat

    "Current steering input [-1, 1]."
    steering::SimFloat
    "Current depower input [0, 1]."
    depower::SimFloat

    "Raw left steering segment reference (for YAML resolution)."
    steer_left_ref::Union{NameRef, Nothing}
    "Raw right steering segment reference (for YAML resolution)."
    steer_right_ref::Union{NameRef, Nothing}
    "Raw depower segment reference (for YAML resolution)."
    depower_ref::Union{NameRef, Nothing}
end

"""
    SteeringConfig(; steer_left, steer_right, steer_gain,
                     depower_segment=nothing, depower_gain=0.0)

Construct a `SteeringConfig` from segment references.

# Keyword Arguments
- `steer_left`: Reference to the left steering segment
  (name or index).
- `steer_right`: Reference to the right steering segment
  (name or index).
- `steer_gain::Float64`: Maximum length offset per side [m].
- `depower_segment=nothing`: Reference to a dedicated depower
  segment. If `nothing`, depower is applied to the steering
  segments.
- `depower_gain::Float64=0.0`: Maximum depower length
  offset [m].
"""
function SteeringConfig(; steer_left, steer_right,
                          steer_gain,
                          depower_segment=nothing,
                          depower_gain=0.0)
    sl = steer_left isa Integer ? Int(steer_left) :
         Symbol(steer_left)
    sr = steer_right isa Integer ? Int(steer_right) :
         Symbol(steer_right)
    dp = if isnothing(depower_segment)
        nothing
    elseif depower_segment isa Integer
        Int(depower_segment)
    else
        Symbol(depower_segment)
    end
    return SteeringConfig(
        0, 0, Float64(steer_gain),
        0, Float64(depower_gain),
        0.0, 0.0, 0.0,
        0.0, 0.0,
        sl, sr, dp)
end

"""
    apply_steering_config!(sys_struct)

Apply the current steering and depower values to segment rest
lengths. Called by `set_steering!` and `set_depower!`.

Steering offsets: left = -s * gain, right = +s * gain.
Depower lengthens by d * depower_gain.

If `depower_idx > 0` (dedicated segment):
- `seg_left.l0  = l0_left + left_offset`
- `seg_right.l0 = l0_right + right_offset`
- `seg_depower.l0 = l0_depower + d * depower_gain`

If `depower_idx == 0` (applied to steering lines):
- `seg_left.l0  = l0_left + left_offset + d * depower_gain`
- `seg_right.l0 = l0_right + right_offset + d * depower_gain`
"""
function apply_steering_config!(sys_struct)
    cfg = sys_struct.steering_config
    isnothing(cfg) && error(
        "No steering_config on this SystemStructure")
    @unpack segments = sys_struct
    s = cfg.steering
    d = cfg.depower

    left_offset  = -s * cfg.steer_gain
    right_offset =  s * cfg.steer_gain

    if cfg.depower_idx > 0
        # Dedicated depower segment
        segments[cfg.steer_left_idx].l0 =
            cfg.l0_steer_left + left_offset
        segments[cfg.steer_right_idx].l0 =
            cfg.l0_steer_right + right_offset
        segments[cfg.depower_idx].l0 =
            cfg.l0_depower + d * cfg.depower_gain
    else
        # Depower applied to steering lines
        depower_offset = d * cfg.depower_gain
        segments[cfg.steer_left_idx].l0 =
            cfg.l0_steer_left + left_offset +
            depower_offset
        segments[cfg.steer_right_idx].l0 =
            cfg.l0_steer_right + right_offset +
            depower_offset
    end
    return nothing
end

"""
    set_steering!(sys_struct, value)

Set the normalized steering input and update segment lengths.

# Arguments
- `sys_struct::SystemStructure`: The system structure.
- `value::Real`: Steering input in [-1, 1].
  Negative = left turn, positive = right turn.
"""
function set_steering!(sys_struct, value::Real)
    cfg = sys_struct.steering_config
    isnothing(cfg) && error(
        "No steering_config on this SystemStructure")
    cfg.steering = clamp(Float64(value), -1.0, 1.0)
    apply_steering_config!(sys_struct)
    return nothing
end

"""
    set_depower!(sys_struct, value)

Set the normalized depower input and update segment lengths.

# Arguments
- `sys_struct::SystemStructure`: The system structure.
- `value::Real`: Depower input in [0, 1].
  0 = fully powered, 1 = fully depowered.
"""
function set_depower!(sys_struct, value::Real)
    cfg = sys_struct.steering_config
    isnothing(cfg) && error(
        "No steering_config on this SystemStructure")
    cfg.depower = clamp(Float64(value), 0.0, 1.0)
    apply_steering_config!(sys_struct)
    return nothing
end

"""
    capture_steering_base_lengths!(sys_struct)

Capture the current segment l0 values as the base lengths
for the steering config. Called during SystemStructure
construction after segment l0 values are finalized.
"""
function capture_steering_base_lengths!(sys_struct)
    cfg = sys_struct.steering_config
    isnothing(cfg) && return nothing
    @unpack segments = sys_struct
    cfg.l0_steer_left = segments[cfg.steer_left_idx].l0
    cfg.l0_steer_right = segments[cfg.steer_right_idx].l0
    if cfg.depower_idx > 0
        cfg.l0_depower = segments[cfg.depower_idx].l0
    end
    return nothing
end

"""
    resolve_steering_config!(cfg::SteeringConfig, segment_names)

Resolve symbolic segment references to indices.
"""
function resolve_steering_config!(
    cfg::SteeringConfig,
    segment_names::Dict{Symbol, Int64}
)
    if !isnothing(cfg.steer_left_ref)
        cfg.steer_left_idx = resolve_ref(
            cfg.steer_left_ref, segment_names,
            "segment")
    end
    if !isnothing(cfg.steer_right_ref)
        cfg.steer_right_idx = resolve_ref(
            cfg.steer_right_ref, segment_names,
            "segment")
    end
    if !isnothing(cfg.depower_ref)
        cfg.depower_idx = resolve_ref(
            cfg.depower_ref, segment_names,
            "segment")
    end
    return nothing
end
