# SPDX-FileCopyrightText: 2025 Bart van de Lint
#
# SPDX-License-Identifier: MPL-2.0

"""
    @recipe function f(sys::PointMassSystem; zoom=false, front=false, kite_pos=nothing, reltime=nothing)

A plot recipe for a `PointMassSystem` type (assumed to be defined, e.g., in `KiteModels.jl`)
using `RecipesBase.jl`.

This recipe allows plotting a `PointMassSystem` in 2D, showing either a
side view (X-Z plane) or a front view (Y-Z plane). It supports zooming
in on the last point of the system (assumed to be the kite or primary object of interest).

## Type Assumptions for `PointMassSystem`

This recipe assumes `PointMassSystem` has fields:
- `points::Vector{Point}`: Where each `Point` object has a `pos_w` field
  (e.g., `p.pos_w`) that yields a 3D position vector (like `SVector{3, Float64}`).
- `segments::Vector{Segment}`: Where each `Segment` object has a `points` field
  (e.g., `seg.points`) that yields a tuple of two integers, representing
  1-based indices into the `sys.points` vector.

Adjust field access within the recipe if your struct definitions differ.

## Attributes

These can be passed as keyword arguments to the `plot` call:

* `zoom::Bool` (default: `get(plotattributes, :zoom, false)`):
  If `true`, the plot view zooms in on the last point in `sys.points`.
* `front::Bool` (default: `get(plotattributes, :front, false)`):
  If `true`, shows the front view (Y-Z). Otherwise, shows the side view (X-Z).
* `kite_pos::Union{AbstractVector{<:Real}, Nothing}` (default: `get(plotattributes, :kite_pos, nothing)`):
  Optionally, an additional 3D position vector (e.g., `SVector{3, Float64}(x,y,z)` or `[x,y,z]`)
  to be plotted as a distinct point. This point is included in zoom calculations if `zoom` is true
  and `kite_pos` is the last effective point.
* `reltime::Union{Real, Nothing}` (default: `get(plotattributes, :reltime, nothing)`):
  If provided, sets the plot title to "Time: [reltime] s".

## Example

```julia
# Make sure PointMassSystem, Point, Segment are defined and you have an instance `my_system`.
# using Plots # Or any other RecipesBase-compatible plotting backend

# Basic side view
# plot(my_system)

# Zoomed side view
# plot(my_system, zoom=true)

# Front view
# plot(my_system, front=true)

# Zoomed front view with a specific kite position and time display
# extra_kite_marker = SVector(10.0, 2.0, 30.0) # Or KVec3(10.0, 2.0, 30.0)
# current_time = 1.23
# plot(my_system, zoom=true, front=true, kite_pos=extra_kite_marker, reltime=current_time)
```
"""
@recipe function f(sys::PointMassSystem) # Now using the actual PointMassSystem type
    # --- Retrieve Plot Attributes ---
    user_zoom = get(plotattributes, :zoom, false)
    user_front = get(plotattributes, :front, false)
    user_kite_pos_input = get(plotattributes, :kite_pos, nothing)
    user_reltime = get(plotattributes, :reltime, nothing)

    # --- Prepare Data for Plotting ---
    # Extract 3D positions of all points in the system.
    # Assumes each Point object `p` in `sys.points` has a `pos_w` field.
    points_positions = [p.pos_w for p in sys.points] # Accesses p.pos_w

    # If an additional `kite_pos` is provided, convert it to SVector and add it to the list.
    if !isnothing(user_kite_pos_input)
        local_kite_svector::SVector{3, Float64} # Assuming KVec3 is compatible with SVector{3,Float64}
        if user_kite_pos_input isa SVector{3, Float64}
            local_kite_svector = user_kite_pos_input
        elseif typeof(user_kite_pos_input).name.name == :KVec3 # Check if it's KVec3 (adapt if KVec3 is an alias)
             local_kite_svector = SVector{3, Float64}(user_kite_pos_input[1], user_kite_pos_input[2], user_kite_pos_input[3])
        elseif user_kite_pos_input isa AbstractVector{<:Real} && length(user_kite_pos_input) == 3
            local_kite_svector = SVector{3, Float64}(user_kite_pos_input[1], user_kite_pos_input[2], user_kite_pos_input[3])
        else
            error("`kite_pos` must be a 3-element AbstractVector, SVector{3, Float64}, KVec3, or nothing.")
        end
        points_positions = [points_positions..., local_kite_svector]
    end

    # Determine which coordinate indices to use for 2D projection.
    x_plot_idx = user_front ? 2 : 1
    y_plot_idx = 3

    # --- Set Plot Limits (xlims, ylims) ---
    if !isempty(points_positions)
        if user_zoom
            last_point_for_zoom = points_positions[end]
            xlims --> (last_point_for_zoom[x_plot_idx] - 5, last_point_for_zoom[x_plot_idx] + 5)
            ylims --> (last_point_for_zoom[y_plot_idx] - 8, last_point_for_zoom[y_plot_idx] + 2)
        else
            if user_front
                xlims --> (-30, 30)
                ylims --> (0, 60)
            else
                xlims --> (0, 60)
                ylims --> (0, 60)
            end
        end
    else
        xlims --> (-1, 1)
        ylims --> (-1, 1)
    end


    # --- Set General Plot Attributes ---
    aspect_ratio := :equal
    legend := false
    xlabel --> (user_front ? "Y [m]" : "X [m]")
    ylabel --> "Z [m]"

    if !isnothing(user_reltime)
        title --> string("Time: ", round(user_reltime, digits=1), " s")
    end

    # --- Define Series for Plotting ---
    # Plot Segments
    if !isempty(sys.segments) && !isempty(sys.points)
        # Use original system points for segments, before potentially adding kite_pos
        original_system_points_pos = [p.pos_w for p in sys.points]

        for seg in sys.segments
            # `seg.points` is assumed to be a tuple of 1-based indices (e.g., (1, 2))
            p1_idx = seg.points[1] # Accesses seg.points
            p2_idx = seg.points[2]

            if (1 <= p1_idx <= length(original_system_points_pos)) && (1 <= p2_idx <= length(original_system_points_pos))
                @series begin
                    seriestype := :path
                    linecolor --> :black
                    linewidth --> 1.5
                    label := ""
                    (
                        [original_system_points_pos[p1_idx][x_plot_idx], original_system_points_pos[p2_idx][x_plot_idx]],
                        [original_system_points_pos[p1_idx][y_plot_idx], original_system_points_pos[p2_idx][y_plot_idx]]
                    )
                end
            else
                # @warn "Segment indices ($p1_idx, $p2_idx) out of bounds for sys.points."
            end
        end
    end

    # Plot Points (including the optional `kite_pos`)
    if !isempty(points_positions)
        @series begin
            seriestype := :scatter
            markercolor --> :blue
            markersize --> 3
            markerstrokecolor --> :match
            label := ""
            (
                [p[x_plot_idx] for p in points_positions],
                [p[y_plot_idx] for p in points_positions]
            )
        end
    end

    # --- Finalize Recipe ---
    # seriestype := :none
    return ()
end
