using SymbolicAWEModels
using VortexStepMethod
using LinearAlgebra
using Statistics
using GLMakie
using KiteUtils
using DiscretePIDs
using Dates

function resolve_log_path(log_name::String, base_dir::String)
    candidates = String[]
    push!(candidates, log_name)
    if !endswith(log_name, ".arrow")
        push!(candidates, log_name * ".arrow")
    end
    push!(candidates, joinpath(base_dir, log_name))
    if !endswith(log_name, ".arrow")
        push!(candidates, joinpath(base_dir, log_name * ".arrow"))
    end
    for candidate in candidates
        if isfile(candidate)
            return candidate
        end
    end
    return candidates[end]
end

function adjust_tether_length!(sam::SymbolicAWEModel, tether_length_raw; tether_point_idxs=39:44)
    tether_length = float(tether_length_raw)
    sys = sam.sys_struct
    set = sam.set

    if !isempty(set.l_tethers)
        set.l_tethers[1] = tether_length
    end

    n_points = length(tether_point_idxs)
    for (n, p_idx) in enumerate(tether_point_idxs)
        pos = (0.0, 0.0, -n * tether_length / n_points)
        sys.points[p_idx].pos_cad .= pos
        sys.points[p_idx].pos_b .= pos
    end

    if !isempty(sys.transforms)
        transform = sys.transforms[1]
        if !isempty(sys.wings) && norm(sys.wings[1].pos_w) > 0
            target_pos = normalize(sys.wings[1].pos_w) * tether_length
            transform.elevation = KiteUtils.calc_elevation(target_pos)
            transform.azimuth = KiteUtils.azimuth_east(target_pos)
        end
        SymbolicAWEModels.reinit!([transform], sys)
    end

    if !isempty(sys.winches)
        winch = sys.winches[1]
        winch.tether_len = tether_length
        winch.tether_vel = 0.0
        winch.brake = true
    end
    return nothing
end

function load_log_and_system(; log_name::String)
    # Extract up, us, v_wind, and lt from log_name. Support multi-us tags like ..._us_10_25_50_...
    m = match(r"_up_([0-9]+)_us_([0-9._-]+)_vw_([0-9]+)_lt_([0-9]+)", log_name)
    m === nothing && error("Could not parse up/us/vw/lt from log name: $log_name")
    up = parse(Float64, m.captures[1])
    us_tokens = split(m.captures[2], "_")
    us_vals = parse.(Float64, us_tokens)
    v_wind = parse(Int, m.captures[3])
    lt = parse(Int, m.captures[4])
    @info "up=$(up/100) [-]" us=us_vals ./ 100 v_wind=v_wind tether_length=lt
    initial_damping = 100.0

    # Load settings
    wing_type = SymbolicAWEModels.REFINE
    wing_type_str = "REFINE"
    @info "Running v3 kite simulation with REFINE wing type..."

    set_data_path("data/v3")
    set = Settings("system.yaml")
    set.v_wind = v_wind
    set.upwind_dir = -90.0
    if !isempty(set.l_tethers)
        set.l_tethers[1] = lt
    end

    # Load YAML structure path
    model_name = "v3_refine"
    struc_yaml_path = joinpath("data", "v3", "CORRECT_struc_geometry.yaml")

    # Load VSMSettings
    vsm_set_path = joinpath(get_data_path(), "CORRECT_vsm_settings.yaml")
    vsm_set = VortexStepMethod.VSMSettings(vsm_set_path; data_prefix=false)

    # Use 36 panels for both wing types (matches vsm_settings.yaml default)
    vsm_set.wings[1].n_panels = 36
    # Note: n_unrefined_sections is automatically inferred from YAML geometry

    # Load system structure with wing_type and vsm_set parameters
    sys = load_sys_struct_from_yaml(struc_yaml_path;
        system_name=model_name, set, wing_type, vsm_set)

    # Initialize damping
    SymbolicAWEModels.set_world_frame_damping(sys, initial_damping)

    wing_points = [p for p in sys.points if p.type == WING]
    n_unrefined = sys.wings[1].vsm_wing.n_unrefined_sections
    @info "REFINE wing setup:" n_wing_points=length(wing_points) n_groups=length(sys.groups) n_unrefined=n_unrefined n_panels=length(sys.wings[1].vsm_aero.panels) n_segments=length(sys.segments)

    # Create symbolic model
    sam = SymbolicAWEModel(set, sys)
    adjust_tether_length!(sam, lt)
    data_dir = normpath(joinpath(@__DIR__, "..", "..", "processed_data", "v3_kite"))
    log_path = resolve_log_path(log_name, data_dir)
    lg = KiteUtils.load_log(log_path)
    return lg, sam, up, us_vals, v_wind, lt
end

function print_and_plot_wing(lg, sam; is_print::Bool=false)

   # Grab the last recorded SysState from the log
   lg_last = lg.syslog[end]

   wing = sam.sys_struct.wings[1]
   origin_idx = wing.origin_idx
   origin_w = [
      lg_last.X[origin_idx],
      lg_last.Y[origin_idx],
      lg_last.Z[origin_idx],
   ]
   R_b_w = SymbolicAWEModels.quaternion_to_rotation_matrix(lg_last.orient)
   wing_point_idxs = [p.idx for p in sam.sys_struct.points if p.type == SymbolicAWEModels.WING]

   # println("Wing node positions (world frame):")
   # for idx in wing_point_idxs
   #    println("$(lg_last.X[idx]), $(lg_last.Y[idx]), $(lg_last.Z[idx])")
   # end

   if is_print
      println("\n# Wing node positions (world frame):")
      for idx in wing_point_idxs
         pos_w = [
            lg_last.X[idx],
            lg_last.Y[idx],
            lg_last.Z[idx],
         ]
         println("- [$idx, [$(Float64(pos_w[1])), $(Float64(pos_w[2])), $(Float64(pos_w[3]))], WING, 1, 1, 0.0, 10.0, 0.0]")
      end

      println("\n# Wing node positions (body frame):")
      for idx in wing_point_idxs
         pos_w = [
            lg_last.X[idx],
            lg_last.Y[idx],
            lg_last.Z[idx],
         ]
         pos_b = R_b_w' * (pos_w .- origin_w)
         println("- [$idx, [$(Float64(pos_b[1])), $(Float64(pos_b[2])), $(Float64(pos_b[3]))], WING, 1, 1, 0.0, 10.0, 0.0]")
      end

      ## print bridle nodes in body frame
      # Bridle points are points 22:38 in the v3 geometry (DYNAMIC points tied by bridle segments)
      # Symmetric pairs: (22,25), (23,24), (26,27), (28,31), (29,30), (32,33), (34,36), (35), (37,38)
      # For each pair, we use the minimum y-coordinate as baseline
      bridle_pairs = [
         (22, 25), (23, 24), (26, 27), (28, 31), (29, 30), 
         (32, 33), (34, 36), (37, 38)
      ]
      bridle_center = [35]
      
      println("\n# Bridle node positions (body frame):")
      
      # Process symmetric pairs
      for (idx_pos, idx_neg) in bridle_pairs
         # Get positions for positive y node
         pos_w_pos = [lg_last.X[idx_pos], lg_last.Y[idx_pos], lg_last.Z[idx_pos]]
         pos_b_pos = R_b_w' * (pos_w_pos .- origin_w)
         
         # Get positions for negative y node
         pos_w_neg = [lg_last.X[idx_neg], lg_last.Y[idx_neg], lg_last.Z[idx_neg]]
         pos_b_neg = R_b_w' * (pos_w_neg .- origin_w)
         
         # Calculate center point between the two y-coordinates
         y_center = (pos_b_pos[2] + pos_b_neg[2]) / 2.0
         
         # Calculate offset from center (half the distance between them)
         y_offset = (pos_b_pos[2] - pos_b_neg[2]) / 2.0
         
         # Print positive y node (adjusted to be symmetric)
         println("- [$idx_pos, [$(Float64(pos_b_pos[1])), $(Float64(y_center + y_offset)), $(Float64(pos_b_pos[3]))], DYNAMIC, 1, 1, 0.0, 30.000, 0.0, 0.0, 0.0]")
         
         # Print negative y node (adjusted to be symmetric)
         println("- [$idx_neg, [$(Float64(pos_b_neg[1])), $(Float64(y_center - y_offset)), $(Float64(pos_b_neg[3]))], DYNAMIC, 1, 1, 0.0, 30.000, 0.0, 0.0, 0.0]")
      end
      
      # Process center node
      for idx in bridle_center
         pos_w = [lg_last.X[idx], lg_last.Y[idx], lg_last.Z[idx]]
         pos_b = R_b_w' * (pos_w .- origin_w)
         println("- [$idx, [$(Float64(pos_b[1])), $(Float64(pos_b[2])), $(Float64(pos_b[3]))], DYNAMIC, 1, 1, 0.1, 30.000, 0.0, 0.0, 0.0]")
      end
   end

   # 2D scatter plots of wing nodes in body frame
   xs_b = Float64[]
   ys_b = Float64[]
   zs_b = Float64[]
   for idx in wing_point_idxs
      pos_w = [
         lg_last.X[idx],
         lg_last.Y[idx],
         lg_last.Z[idx],
      ]
      pos_b = R_b_w' * (pos_w .- origin_w)
      push!(xs_b, pos_b[1]); push!(ys_b, pos_b[2]); push!(zs_b, pos_b[3])
   end

   fig2 = Figure(resolution=(900, 300))
   ax_xy = Axis(fig2[1, 1], title="Top (x,y)", xlabel="y_b", ylabel="x_b")
   ax_xz = Axis(fig2[1, 2], title="Side (x,z)", xlabel="x_b", ylabel="z_b")
   ax_yz = Axis(fig2[1, 3], title="Rear (y,z)", xlabel="y_b", ylabel="z_b")

   ax_xy.aspect = DataAspect()
   ax_xz.aspect = DataAspect()
   ax_yz.aspect = DataAspect()

   xs_b_rot = ys_b
   ys_b_rot = [-x for x in xs_b]

   scatter!(ax_xy, xs_b_rot, ys_b_rot)
   scatter!(ax_xz, xs_b, zs_b)
   scatter!(ax_yz, ys_b, zs_b)

   # enforce equal scale/step size across all plots (square span)
   spans = [
      maximum(xs_b_rot) - minimum(xs_b_rot),
      maximum(ys_b_rot) - minimum(ys_b_rot),
      maximum(xs_b) - minimum(xs_b),
      maximum(zs_b) - minimum(zs_b),
      maximum(ys_b) - minimum(ys_b),
   ]
   global_span = maximum(spans)
   global_span = global_span > 0 ? global_span : 1.0
   global_span *= 1.05  # slight padding

   function set_square_limits!(ax, xs, ys, span)
      cx = mean(xs)
      cy = mean(ys)
      half = span / 2
      xlims!(ax, cx - half, cx + half)
      ylims!(ax, cy - half*0.5, cy + half*0.5)
   end

   set_square_limits!(ax_xy, xs_b_rot, ys_b_rot, global_span)
   set_square_limits!(ax_xz, xs_b, zs_b, global_span)
   set_square_limits!(ax_yz, ys_b, zs_b, global_span)

   # Draw wing structural segments for the plotted wing (indices 1:19 in YAML)
   wing_seg_idxs = 1:19  # leading edge + struts
   lines_xy = Point2f[]
   lines_xz = Point2f[]
   lines_yz = Point2f[]
   for seg_idx in wing_seg_idxs
      seg = sam.sys_struct.segments[seg_idx]
      p1, p2 = seg.point_idxs
      pos1_w = [lg_last.X[p1], lg_last.Y[p1], lg_last.Z[p1]]
      pos2_w = [lg_last.X[p2], lg_last.Y[p2], lg_last.Z[p2]]
      pos1_b = R_b_w' * (pos1_w .- origin_w)
      pos2_b = R_b_w' * (pos2_w .- origin_w)
      x1, y1, z1 = pos1_b
      x2, y2, z2 = pos2_b
      # rotate only the xy projection by 90 deg clockwise: (x,y)->(y,-x)
      push!(lines_xy, Point2f(y1, -x1), Point2f(y2, -x2), Point2f(NaN, NaN)) # NaN breaks segments
      push!(lines_xz, Point2f(x1, z1), Point2f(x2, z2), Point2f(NaN, NaN))
      push!(lines_yz, Point2f(y1, z1), Point2f(y2, z2), Point2f(NaN, NaN))
   end

   lines!(ax_xy, lines_xy, color=:gray)
   lines!(ax_xz, lines_xz, color=:gray)
   lines!(ax_yz, lines_yz, color=:gray)

   return fig2
end

function plot_time_series(lg, sam)
    fig = plot(sam.sys_struct, 
               lg;
               plot_turn_rates=false, 
               plot_reelout=false,
               plot_twist=false,
               plot_yaw_rate_paper=false,
               yaw_rate_paper_ylims=(-90.0, 90.0),
               yaw_rate_paper_compare=false, 
               plot_v_app=true,
               plot_kite_vel=false,
               plot_gk=false,
               gk_ylims=(0.0, 15.0),
               plot_aoa=true,
               aoa_ylims=(0.0, 15.0), 
               plot_heading=false, 
               plot_elevation=true,
               plot_azimuth=true, 
               plot_winch_force=true, 
               plot_set_values=false,
               plot_us=true,
               plot_tether_actual=false,
               plot_turn_radius=false,
               turn_radius_ylims=(0.0, 40.0),
               plot_cs=false,
               cs_ylims=(0.0, 0.02),
         )
      return fig
end

"""
    report_tether_direction_alignment(lg)

Print the key directions for tension derivation:
- Bridle (mid-LE → KCU/point 1)
- Segment 90 (point 1 → point 39)
"""
function report_tether_direction_alignment(lg)
    if !isempty(lg.syslog)
        idxs = unique([1, cld(length(lg.syslog), 2), length(lg.syslog)])
        for idx in idxs
            sl = lg.syslog[idx]
            p1 = [sl.X[1], sl.Y[1], sl.Z[1]]                # KCU/bridle hub
            ple12 = [sl.X[12], sl.Y[12], sl.Z[12]]
            ple14 = [sl.X[14], sl.Y[14], sl.Z[14]]
            p_le_mid = (ple12 .+ ple14) ./ 2
            bridle_dir = p1 .- p_le_mid
            midLE_to_KCU = norm(bridle_dir) > 0 ? bridle_dir / norm(bridle_dir) : bridle_dir

            # Segment 90 endpoints (KCU point 1 to point 39)
            p39 = [sl.X[39], sl.Y[39], sl.Z[39]]
            seg90_dir = p39 .- p1
            seg90_dir_unit = norm(seg90_dir) > 0 ? seg90_dir / norm(seg90_dir) : seg90_dir

            @info "Direction (sample $idx)" midLE_to_KCU seg90_dir_unit
        end
    end
end


"""
    compute_line_stretch(lg, sam; window_seconds=50.0)

Compute signed line stretch ratios (elongation positive, compression negative)
using the log data and segment rest lengths from `sam`. Returns a NamedTuple
with the time window used, per-category stretch matrices (rows = samples,
cols = segments), and per-pulley combined stretch vectors. Also logs mean/max
stretch for elongation and compression over the selected window, including the
sample/segment index where the maxima occur. Bridle category excludes segments
that belong to pulleys; those are grouped under a separate pulley category.

Optional `segment_l0_adjustments` lets you account for commanded line changes
(e.g. steering/power). Provide a Dict from segment idx to either a scalar offset
to add to the nominal `l0`, or a vector of offsets (one per sample) for
time-varying commands. Pulley pairs are evaluated using the combined length of
both segments (accounting for any `l0` adjustments).

Optional `tether_length` rescales tether segment rest lengths to match the
logged tether length before stretch ratios are computed; if omitted, the
logged winch length (`l_tether`) is used when available.
"""
function compute_line_stretch(lg, sam; window_seconds::Real=50.0, segment_l0_adjustments=nothing, tether_length=nothing)
   sl = hasproperty(lg, :syslog) ? lg.syslog : lg
   if isempty(sl)
      @warn "compute_line_stretch: empty log provided"
      return (window=(0.0, 0.0), ratio=Dict{Symbol, Matrix{Float64}}())
   end

   segments = sam.sys_struct.segments
   base_l0 = [Float64(seg.l0) for seg in segments]
   tether_seg_idxs = if !isempty(sam.sys_struct.tethers) && !isempty(sam.sys_struct.tethers[1].segment_idxs)
      Int.(sam.sys_struct.tethers[1].segment_idxs)
   else
      collect(90:95)
   end
   tether_seg_idxs = [idx for idx in tether_seg_idxs if 1 <= idx <= length(segments)]
   tether_seg_set = Set(tether_seg_idxs)
   nominal_total = sum(l0 for (i, l0) in enumerate(base_l0) if i in tether_seg_set && isfinite(l0) && l0 > 0)

   # Build a per-sample scaling factor for tether rest lengths.
   # Falls back to the logged l_tether (winch length) if no explicit value is provided.
   tether_len_series = Float64[]
   if tether_length === nothing
      # Prefer log values when available
      tether_len_series = fill(NaN, length(sl))
      if hasproperty(sl[1], :l_tether)
         for (idx, state) in enumerate(sl)
            lt_vec = getproperty(state, :l_tether)
            if !isempty(lt_vec)
               lt = float(lt_vec[1])
               tether_len_series[idx] = lt
            end
         end
      end
   elseif tether_length isa AbstractVector
      tether_len_series = [float(tether_length[min(i, length(tether_length))]) for i in 1:length(sl)]
   else
      tether_len_series = fill(float(tether_length), length(sl))
   end

   tether_scale = ones(Float64, length(sl))
   if nominal_total > 0 && !isempty(tether_seg_set)
      last_scale = 1.0
      for i in 1:length(sl)
         tl = tether_len_series[i]
         if isfinite(tl) && tl > 0
            last_scale = tl / nominal_total
         end
         tether_scale[i] = last_scale
      end
   elseif tether_length !== nothing && (nominal_total <= 0 || isempty(tether_seg_set))
      @warn "compute_line_stretch: no valid tether segments for scaling" nominal_total tether_seg_idxs
   end
   t_end = sl[end].time
   t_start = t_end - window_seconds
   start_idx = 1
   for i in 1:length(sl)
      if sl[i].time >= t_start
         start_idx = i
         break
      end
   end
   window = (sl[start_idx].time, t_end)
   window_span = window[2] - window[1]
   n_samples = length(sl) - start_idx + 1
   sample_times = [sl[i].time for i in start_idx:length(sl)]
   l0_adjustments = segment_l0_adjustments === nothing ? Dict{Int, Any}() : segment_l0_adjustments
   pulleys = sam.sys_struct.pulleys
   pulley_seg_set = Set{Int}()
   for p in pulleys
      push!(pulley_seg_set, Int(p.segment_idxs[1]))
      push!(pulley_seg_set, Int(p.segment_idxs[2]))
   end
   seg_len_for_pulleys = Dict(seg => fill(NaN, n_samples) for seg in pulley_seg_set)
   seg_l0_for_pulleys = Dict(seg => fill(NaN, n_samples) for seg in pulley_seg_set)
   bridle_seg_idxs = [i for i in 47:86 if !(i in pulley_seg_set)]  # Exclude control tapes (87-89)
   steering_tape_idxs = [87, 89]  # Steering tapes
   power_tape_idxs = [88]  # Power tape

   # Precompute lengths for pulley segments (kept separate from bridle category)
   if !isempty(pulley_seg_set)
      @inbounds for (sample_idx, log_idx) in enumerate(start_idx:length(sl))
         state = sl[log_idx]
         X, Y, Z = state.X, state.Y, state.Z
         for seg_idx in pulley_seg_set
            seg = segments[seg_idx]
            l0 = base_l0[seg_idx]
            if seg_idx in tether_seg_set
               l0 *= tether_scale[log_idx]
            end
            adj = get(l0_adjustments, seg_idx, nothing)
            if adj !== nothing
               if adj isa AbstractVector
                  idx = min(sample_idx, length(adj))
                  l0 += adj[idx]
               elseif adj isa Real
                  l0 += adj
               end
            end
            if !isfinite(l0) || l0 <= 0
               continue
            end
            p1, p2 = seg.point_idxs
            dx = X[p2] - X[p1]
            dy = Y[p2] - Y[p1]
            dz = Z[p2] - Z[p1]
            len = sqrt(dx * dx + dy * dy + dz * dz)
            if isfinite(len)
               seg_len_for_pulleys[seg_idx][sample_idx] = len
               seg_l0_for_pulleys[seg_idx][sample_idx] = l0
            end
         end
      end
   end

   categories = (
      (:tubular_frame, "Tubular frame", 1:19),
      (:te_wires_and_diagonals, "TE wires and diagonals", 20:46),
      (:bridles, "Bridles", bridle_seg_idxs),
      (:steering_tapes, "Steering tapes", steering_tape_idxs),
      (:power_tape, "Power tape", power_tape_idxs),
      (:tether, "Tether", 90:95),
   )

   ratio_by_category = Dict{Symbol, Matrix{Float64}}()

   for (key, label, seg_idxs) in categories
      ratios = fill(NaN, n_samples, length(seg_idxs))
      l0_used = fill(NaN, n_samples, length(seg_idxs))
      @inbounds for (sample_idx, log_idx) in enumerate(start_idx:length(sl))
         state = sl[log_idx]
         X, Y, Z = state.X, state.Y, state.Z
         for (col_idx, seg_idx) in enumerate(seg_idxs)
            seg = segments[seg_idx]
            l0 = base_l0[seg_idx]
            if seg_idx in tether_seg_set
               l0 *= tether_scale[log_idx]
            end
            adj = get(l0_adjustments, seg_idx, nothing)
            if adj !== nothing
               if adj isa AbstractVector
                  idx = min(sample_idx, length(adj))
                  l0 += adj[idx]
               elseif adj isa Real
                  l0 += adj
               end
            end
            if !isfinite(l0) || l0 <= 0
               continue
            end
            p1, p2 = seg.point_idxs
            dx = X[p2] - X[p1]
            dy = Y[p2] - Y[p1]
            dz = Z[p2] - Z[p1]
            len = sqrt(dx * dx + dy * dy + dz * dz)
            if isfinite(len)
               ratios[sample_idx, col_idx] = (len - l0) / l0
               l0_used[sample_idx, col_idx] = l0
               if seg_idx in pulley_seg_set
                  seg_len_for_pulleys[seg_idx][sample_idx] = len
                  seg_l0_for_pulleys[seg_idx][sample_idx] = l0
               end
            end
         end
      end

      finite_mask = isfinite.(ratios)
      if any(finite_mask)
         elong_mask = finite_mask .& (ratios .> 0)
         comp_mask = finite_mask .& (ratios .< 0)

         mean_elong = any(elong_mask) ? mean(ratios[elong_mask]) : NaN
         mean_comp = any(comp_mask) ? mean(ratios[comp_mask]) : NaN
         mean_elong_abs = any(elong_mask) ? mean(ratios[elong_mask] .* l0_used[elong_mask]) : NaN
         mean_comp_abs = any(comp_mask) ? mean(ratios[comp_mask] .* l0_used[comp_mask]) : NaN

         max_elong = NaN
         max_elong_info = missing
         max_elong_abs = NaN
         if any(elong_mask)
            masked = copy(ratios)
            masked[.!elong_mask] .= -Inf
            max_elong, linear_idx = findmax(masked)
            cart_idx = CartesianIndices(size(ratios))[linear_idx]
            row, col = cart_idx[1], cart_idx[2]
            max_elong_info = (; sample=row, segment=seg_idxs[col], time=sample_times[row])
            max_elong_abs = max_elong * l0_used[row, col]
         end

         max_comp = NaN
         max_comp_info = missing
         max_comp_abs = NaN
         if any(comp_mask)
            masked = copy(ratios)
            masked[.!comp_mask] .= Inf
            max_comp, linear_idx = findmin(masked)
            cart_idx = CartesianIndices(size(ratios))[linear_idx]
            row, col = cart_idx[1], cart_idx[2]
            max_comp_info = (; sample=row, segment=seg_idxs[col], time=sample_times[row])
            max_comp_abs = max_comp * l0_used[row, col]
         end

         elong_mean_str = any(elong_mask) ?
            "mean = $(round(mean_elong_abs, digits=4)) [m], $(round(mean_elong * 100, digits=4)) [%]" :
            "mean = n/a"
         elong_max_str = any(elong_mask) ?
            "max  = $(round(max_elong_abs, digits=4)) [m], $(round(max_elong * 100, digits=4)) [%] segment_idx = $(max_elong_info.segment), time = $(round(max_elong_info.time, digits=2)) [s]" :
            "max  = n/a"

         comp_mean_str = any(comp_mask) ?
            "mean = $(round(mean_comp_abs, digits=4)) [m], $(round(mean_comp * 100, digits=4)) [%]" :
            "mean = n/a"
         comp_max_str = any(comp_mask) ?
            "max  = $(round(max_comp_abs, digits=4)) [m], $(round(max_comp * 100, digits=4)) [%] segment_idx = $(max_comp_info.segment), time = $(round(max_comp_info.time, digits=2)) [s]" :
            "max  = n/a"

         msg = """
           Elongation
             $elong_mean_str
             $elong_max_str
           Compression
             $comp_mean_str
             $comp_max_str
         """
         @info "$label, last $(round(window_span, digits=2)) s\n$msg"
      else
         @warn "$label has no finite values in the selected window"
      end

      ratio_by_category[key] = ratios
   end

   pulley_ratio = Dict{Int, Vector{Float64}}()
   if !isempty(pulleys)
      np = length(pulleys)
      pulley_mat = fill(NaN, n_samples, np)
      pulley_l0 = fill(NaN, n_samples, np)
      for (col_idx, pulley) in enumerate(pulleys)
         s1, s2 = Int.(pulley.segment_idxs)
         len1 = get(seg_len_for_pulleys, s1, fill(NaN, n_samples))
         len2 = get(seg_len_for_pulleys, s2, fill(NaN, n_samples))
         l01 = get(seg_l0_for_pulleys, s1, fill(NaN, n_samples))
         l02 = get(seg_l0_for_pulleys, s2, fill(NaN, n_samples))
         total_l0 = l01 .+ l02
         total_len = len1 .+ len2
         ratios = (total_len .- total_l0) ./ total_l0
         pulley_mat[:, col_idx] .= ratios
         pulley_l0[:, col_idx] .= total_l0
         pulley_ratio[Int(pulley.idx)] = ratios
      end

      ratio_by_category[:pulleys] = pulley_mat

      finite_mask = isfinite.(pulley_mat)
      if any(finite_mask)
         elong_mask = finite_mask .& (pulley_mat .> 0)
         comp_mask = finite_mask .& (pulley_mat .< 0)

         mean_elong = any(elong_mask) ? mean(pulley_mat[elong_mask]) : NaN
         mean_comp = any(comp_mask) ? mean(pulley_mat[comp_mask]) : NaN
         mean_elong_abs = any(elong_mask) ? mean(pulley_mat[elong_mask] .* pulley_l0[elong_mask]) : NaN
         mean_comp_abs = any(comp_mask) ? mean(pulley_mat[comp_mask] .* pulley_l0[comp_mask]) : NaN

         max_elong = NaN
         max_elong_abs = NaN
         max_elong_info = missing
         if any(elong_mask)
            masked = copy(pulley_mat)
            masked[.!elong_mask] .= -Inf
            max_elong, linear_idx = findmax(masked)
            cart_idx = CartesianIndices(size(pulley_mat))[linear_idx]
            row, col = cart_idx[1], cart_idx[2]
            max_elong_abs = max_elong * pulley_l0[row, col]
            max_elong_info = (; sample=row, pulley=Int(pulleys[col].idx), segments=Tuple(Int.(pulleys[col].segment_idxs)), time=sample_times[row])
         end

         max_comp = NaN
         max_comp_abs = NaN
         max_comp_info = missing
         if any(comp_mask)
            masked = copy(pulley_mat)
            masked[.!comp_mask] .= Inf
            max_comp, linear_idx = findmin(masked)
            cart_idx = CartesianIndices(size(pulley_mat))[linear_idx]
            row, col = cart_idx[1], cart_idx[2]
            max_comp_abs = max_comp * pulley_l0[row, col]
            max_comp_info = (; sample=row, pulley=Int(pulleys[col].idx), segments=Tuple(Int.(pulleys[col].segment_idxs)), time=sample_times[row])
         end

         msg = """
           Elongation
             $(any(elong_mask) ? "mean = $(round(mean_elong_abs, digits=4)) [m], $(round(mean_elong * 100, digits=4)) [%]" : "mean = n/a")
             $(any(elong_mask) ? "max  = $(round(max_elong_abs, digits=4)) [m], $(round(max_elong * 100, digits=4)) [%] pulley_idx = $(max_elong_info.pulley), segments = $(max_elong_info.segments), time = $(round(max_elong_info.time, digits=2)) [s]" : "max  = n/a")
           Compression
             $(any(comp_mask) ? "mean = $(round(mean_comp_abs, digits=4)) [m], $(round(mean_comp * 100, digits=4)) [%]" : "mean = n/a")
             $(any(comp_mask) ? "max  = $(round(max_comp_abs, digits=4)) [m], $(round(max_comp * 100, digits=4)) [%] pulley_idx = $(max_comp_info.pulley), segments = $(max_comp_info.segments), time = $(round(max_comp_info.time, digits=2)) [s]" : "max  = n/a")
         """
         @info "Pulleys, last $(round(window_span, digits=2)) s\n$msg"
      else
         @warn "Pulleys has no finite values in the selected window"
      end
   end

   return (window=window, ratio=ratio_by_category, pulley_ratio=pulley_ratio)
end

# log_name = "batch_2026_01_07_10_58_14/circle__up_24_us_15_vw_8_lt_260_run_008_date_2026_01_07_11_26_11"
# log_name = "zenith_circle__up_40_us_20_vw_15_lt_275_date_2026_01_09_14_49"
# log_name = "circular_2019_batch_2026_01_10_12_01_04/circle__up_18_us_15_vw_9_lt_268_run_005_date_2026_01_10_12_11_49"
# log_name = "circular_2025_batch_2026_01_10_11_28_38/circle__up_42_us_20_vw_8_lt_262_run_007_date_2026_01_10_11_42_53"
log_name = "circular_up_varying_batch_2026_01_11_11_29_19/circle__up_38_us_8_vw_8_lt_262_run_062_date_2026_01_11_17_09_20"
# log_name = "zenith_circle__up_42_us_20_vw_8_lt_262_date_2026_01_10_09_54"
lg, sam, up, us, v_wind, lt = load_log_and_system(log_name=log_name)

# Log alignment info before plotting to decide tension source
# report_tether_direction_alignment(lg)

# Account for commanded steering/power when evaluating stretch of segments 87/88/89.
# Mapping matches v3_kite_circles.jl: 1400 mm steering span, 200 mm + 5 m * up for power.
up_fraction = up / 100
us_fraction = us / 100
seg87_nom = Float64(sam.sys_struct.segments[87].l0)
seg88_nom = Float64(sam.sys_struct.segments[88].l0)
seg89_nom = Float64(sam.sys_struct.segments[89].l0)
steering_tape_change = 1.4 * us_fraction
power_target_l0 = 0.2 + 5.0 * up_fraction
segment_l0_adjustments = Dict(
   87 => steering_tape_change,
   88 => power_target_l0 - seg88_nom,
   89 => -steering_tape_change,
)

# Compute line stretch ratios over the last 50 seconds
stretch_info = compute_line_stretch(
   lg, sam;
   window_seconds=50.0,
   segment_l0_adjustments=segment_l0_adjustments,
)

### plot time series
fig_time = plot_time_series(lg, sam)
## show 3D animation
scene = replay(lg, sam.sys_struct; autoplay=false, loop=true, show_panes=false)
### show 2D wing node plots
fig_wing = print_and_plot_wing(lg, sam,is_print=true)

scr1=display(fig_time)
wait(scr1)
scr2=display(scene)
wait(scr2)
scr3=display(fig_wing)
wait(scr3)

##TODO: record does not work
# record(scene, "v3_kite_circular_load_and_plot.mp4"; fps=30, duration=20)  # Adjust duration as needed
