using SymbolicAWEModels
using VortexStepMethod
using LinearAlgebra
using Statistics
using Dates
using StaticArrays
using KiteUtils

const WINDOW_SEC = 100.0

function parse_up_us_vw_lt(log_name::AbstractString)
    m = match(r"_up_([0-9]+)_us_([0-9._-]+)_vw_([0-9]+)_lt_([0-9]+)", log_name)
    m === nothing && return nothing
    up_raw = parse(Float64, m.captures[1])
    us_tokens = split(m.captures[2], "_")
    us_raw = parse(Float64, us_tokens[1])
    v_wind = parse(Int, m.captures[3])
    lt = parse(Int, m.captures[4])
    return up_raw / 100, us_raw / 100, v_wind, lt
end

function adjust_tether_length!(sys::SystemStructure, set::Settings, tether_length_raw; tether_point_idxs=39:44)
    tether_length = float(tether_length_raw)

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

function build_sys(; v_wind=10.0, tether_length=150.0)
    wing_type = SymbolicAWEModels.REFINE
    set_data_path("data/v3")
    set = Settings("system.yaml")
    set.v_wind = v_wind
    set.upwind_dir = -90.0
    if !isempty(set.l_tethers)
        set.l_tethers[1] = tether_length
    end

    model_name = "v3_refine"
    struc_yaml_path = joinpath("data", "v3", "CORRECT_struc_geometry.yaml")
    vsm_set_path = joinpath(get_data_path(), "CORRECT_vsm_settings.yaml")
    vsm_set = VortexStepMethod.VSMSettings(vsm_set_path; data_prefix=false)
    vsm_set.wings[1].n_panels = 36

    sys = load_sys_struct_from_yaml(struc_yaml_path;
        system_name=model_name, set, wing_type, vsm_set)
    adjust_tether_length!(sys, set, tether_length)
    return sys
end

function unwrap_phase!(vals::AbstractVector{<:Real}; period=2 * pi, thresh=pi)
    if isempty(vals)
        return vals
    end
    offset = 0.0
    prev = vals[1]
    for i in 2:length(vals)
        delta = vals[i] - prev
        if delta > thresh
            offset -= period
        elseif delta < -thresh
            offset += period
        end
        prev = vals[i]
        vals[i] += offset
    end
    return vals
end

function gradient_uniform(y::AbstractVector{<:Real}, ts::Real)
    n = length(y)
    grad = Vector{Float64}(undef, n)
    if n == 0
        return grad
    elseif n == 1
        grad[1] = 0.0
        return grad
    end
    grad[1] = (y[2] - y[1]) / ts
    for i in 2:(n - 1)
        grad[i] = (y[i + 1] - y[i - 1]) / (2 * ts)
    end
    grad[n] = (y[n] - y[n - 1]) / ts
    return grad
end

function moving_average_same(x::AbstractVector{<:Real}, window::Int)
    n = length(x)
    if window <= 1 || n == 0
        return Float64.(x)
    end
    left = div(window, 2)
    right = window - 1 - left
    padded = Vector{Float64}(undef, n + left + right)
    padded[1:left] .= 0.0
    padded[(left + 1):(left + n)] .= x
    padded[(left + n + 1):end] .= 0.0
    out = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        s = 0.0
        for k in 0:(window - 1)
            s += padded[i + k]
        end
        out[i] = s / window
    end
    return out
end

"""
    midle_to_kcu_dir(sl, k; eps=1e-12)

Compute the unit vector from the mid leading-edge (avg of points 12 & 14)
to the KCU/bridle hub (point 1) for sample `k` of a syslog entry.
Returns `nothing` if the required points are unavailable or degenerate.
"""
function midle_to_kcu_dir(sl, k; eps=1e-12)
    Xk = sl.X[k]; Yk = sl.Y[k]; Zk = sl.Z[k]
    if length(Xk) < 14 || length(Yk) < 14 || length(Zk) < 14
        return nothing
    end
    p1 = SVector{3, Float64}(Xk[1], Yk[1], Zk[1])
    ple12 = SVector{3, Float64}(Xk[12], Yk[12], Zk[12])
    ple14 = SVector{3, Float64}(Xk[14], Yk[14], Zk[14])
    p_le_mid = (ple12 + ple14) / 2
    dir = p1 - p_le_mid
    n = norm(dir)
    return n > eps ? dir / n : nothing
end

function calc_ref_area(sys::SystemStructure)
    isempty(sys.wings) && return NaN
    wing = sys.wings[1]
    hasproperty(wing, :vsm_aero) || return NaN
    panels = wing.vsm_aero.panels
    isempty(panels) && return NaN
    return sum(p.chord * p.width for p in panels)
end

function calculate_cs(sl, sys; rho=1.225, eps=1e-12)
    s_ref = calc_ref_area(sys)
    if !isfinite(s_ref) || s_ref <= eps
        return Float64[], Float64[]
    end

    n = length(sl.time)
    cs = Vector{Float64}(undef, n)

    @inbounds for k in 1:n
        v_kite = sl.vel_kite[k]
        v_wind = sl.v_wind_kite[k]
        v_a = v_kite - v_wind
        v_a_norm = norm(v_a)
        if v_a_norm <= eps
            cs[k] = NaN
            continue
        end
        drag_dir = -v_a / v_a_norm

        up_dir = midle_to_kcu_dir(sl, k; eps=eps)
        if up_dir === nothing
            cs[k] = NaN
            continue
        end
        up_dir = -up_dir

        side_raw = cross(drag_dir, up_dir)
        side_norm = norm(side_raw)
        if side_norm <= eps
            cs[k] = NaN
            continue
        end
        side_dir = side_raw / side_norm

        R_b_w = SymbolicAWEModels.quaternion_to_rotation_matrix(sl.orient[k])
        F_aero_b = sl.aero_force_b[k]
        F_aero_w = R_b_w * F_aero_b
        F_side = dot(F_aero_w, side_dir)
        cs[k] = F_side / (0.5 * rho * v_a_norm^2 * s_ref)
    end

    return cs, sl.time
end

function compute_turn_radius(sl_in, sys::SystemStructure; smooth_window=10, eps=1e-12)
    sl = hasproperty(sl_in, :syslog) ? sl_in.syslog : sl_in
    n = length(sl.time)
    if n < 2 || isempty(sl.vel_kite) || isempty(sl.orient)
        return nothing
    end
    if length(sl.vel_kite) < n || length(sl.orient) < n
        return nothing
    end

    ts = mean(diff(sl.time))
    ts = isfinite(ts) && ts > eps ? ts : eps

    v_x = Vector{Float64}(undef, n)
    v_y = Vector{Float64}(undef, n)
    v_z = Vector{Float64}(undef, n)
    @inbounds for k in 1:n
        v = sl.vel_kite[k]
        v_x[k] = v[1]
        v_y[k] = v[2]
        v_z[k] = v[3]
    end

    a_x = gradient_uniform(v_x, ts)
    a_y = gradient_uniform(v_y, ts)
    a_z = gradient_uniform(v_z, ts)

    if smooth_window > 1
        a_x = moving_average_same(a_x, smooth_window)
        a_y = moving_average_same(a_y, smooth_window)
        a_z = moving_average_same(a_z, smooth_window)
    end

    radius = Vector{Float64}(undef, n)
    @inbounds for k in 1:n
        v = SVector{3, Float64}(v_x[k], v_y[k], v_z[k])
        a = SVector{3, Float64}(a_x[k], a_y[k], a_z[k])
        v_norm = norm(v)
        if !isfinite(v_norm) || v_norm <= eps
            radius[k] = NaN
            continue
        end
        v_hat = v / v_norm
        a_t = dot(a, v_hat) * v_hat
        omega = cross(a - a_t, v) / (v_norm^2)
        omega_norm = norm(omega)
        if !isfinite(omega_norm) || omega_norm <= eps
            radius[k] = NaN
            continue
        end
        icr = cross(v, omega) / (omega_norm^2)
        R_b_w = SymbolicAWEModels.quaternion_to_rotation_matrix(sl.orient[k])
        e_x = SVector{3, Float64}(R_b_w[:, 1])
        det = e_x[1] * icr[2] - e_x[2] * icr[1]
        if !isfinite(det) || abs(det) <= eps
            radius[k] = NaN
        else
            radius[k] = -(det < 0 ? -1.0 : 1.0) * norm(icr)
        end
    end

    return radius, sl.time
end

function compute_ekf_yaw_and_rate(sl_in, sys::SystemStructure; eps=1e-12)
    sl = hasproperty(sl_in, :syslog) ? sl_in.syslog : sl_in
    n = length(sl.time)
    if n < 2 || isempty(sl.vel_kite)
        return nothing
    end
    if length(sys.wings) == 0 || length(sl.X) < n || length(sl.Y) < n || length(sl.Z) < n
        return nothing
    end

    kite_idx = sys.wings[1].origin_idx
    yaw = Vector{Float64}(undef, n)
    nan_count = 0

    @inbounds for k in 1:n
        pos = SVector{3, Float64}(sl.X[k][kite_idx], sl.Y[k][kite_idx], sl.Z[k][kite_idx])
        vel = SVector{3, Float64}(sl.vel_kite[k])

        npos = norm(pos)
        nvel = norm(vel)

        if npos > eps && nvel > eps
            radial = pos / npos
            tang_vel = vel - dot(vel, radial) * radial
            ntang = norm(tang_vel)

            if ntang > eps
                tang_vel_unit = tang_vel / ntang
                up_z = radial
                up_y_raw = SVector(-pos[2], pos[1], 0.0)
                nup_y = norm(up_y_raw)

                if nup_y > eps
                    up_y = up_y_raw / nup_y
                    up_x = cross(up_z, up_y)
                    nup_x = norm(up_x)

                    if nup_x > eps
                        up_x = up_x / nup_x
                        up_y = cross(up_z, up_x)
                        R_up = @SMatrix [up_x[1] up_y[1] up_z[1];
                                         up_x[2] up_y[2] up_z[2];
                                         up_x[3] up_y[3] up_z[3]]
                        heading_vec = R_up' * tang_vel_unit
                        yaw[k] = atan(heading_vec[2], heading_vec[1])
                        continue
                    end
                end
            end
        end

        yaw[k] = k > 1 ? yaw[k - 1] : NaN
        nan_count += 1
    end

    if nan_count > 0
        @info "compute_ekf_yaw_and_rate: $nan_count samples with degenerate geometry"
    end

    yaw_unwrapped = copy(yaw)
    unwrap_phase!(yaw_unwrapped)

    ts = mean(diff(sl.time))
    ts = isfinite(ts) && ts > eps ? ts : eps

    yaw_rate = gradient_uniform(yaw_unwrapped, ts)
    yaw_rate = moving_average_same(yaw_rate, 10)

    return yaw_unwrapped, rad2deg.(yaw_rate)
end

function unwrap_heading(heading::AbstractVector{<:Real})
    heading_unwrapped = copy(heading)
    for j in 2:length(heading_unwrapped)
        while heading_unwrapped[j] - heading_unwrapped[j - 1] > pi
            heading_unwrapped[j] -= 2 * pi
        end
        while heading_unwrapped[j] - heading_unwrapped[j - 1] < -pi
            heading_unwrapped[j] += 2 * pi
        end
    end
    return heading_unwrapped
end

function heading_rate(sl)
    heading_unwrapped = unwrap_heading(sl.heading)
    rates = diff(rad2deg.(heading_unwrapped)) ./ diff(sl.time)
    times = sl.time[1:end-1]
    return rates, times
end

function steering_command(sl, sys; steering_l0=nothing)
    seg_left = sys.segments[87]
    p_i, p_j = seg_left.point_idxs
    xs = sl.X
    ys = sl.Y
    zs = sl.Z
    n = length(sl.time)
    steering_len = zeros(Float64, n)
    @inbounds for k in 1:n
        p1 = SVector{3, Float64}(xs[k][p_i], ys[k][p_i], zs[k][p_i])
        p2 = SVector{3, Float64}(xs[k][p_j], ys[k][p_j], zs[k][p_j])
        steering_len[k] = norm(p2 - p1)
    end
    base_l0 = isnothing(steering_l0) ? steering_len[1] : steering_l0
    us_cmd = similar(steering_len)
    @inbounds for k in eachindex(us_cmd)
        delta = steering_len[k] - base_l0
        us_cmd[k] = abs(delta) > 1e-6 ? delta / 1.4 : 0.0
    end
    return us_cmd
end

function gk_series(sl, sys)
    heading_rate_deg, _ = heading_rate(sl)
    us_cmd = steering_command(sl, sys; steering_l0=1.6)
    v_app = sl.v_app[2:end]
    us_seg = us_cmd[2:end]
    gk = similar(heading_rate_deg)
    @inbounds for k in eachindex(gk)
        gk[k] = abs(us_seg[k]) > 1e-8 ? heading_rate_deg[k] / (v_app[k] * us_seg[k]) : NaN
    end
    times = sl.time[2:end]
    return gk, times
end

function gk_paper_series(sl, sys)
    n = length(sl.time)
    yaw = Vector{Float64}(undef, n)
    @inbounds for k in 1:n
        v = sl.vel_kite[k]
        w = sl.v_wind_kite[k]
        va_enu = w .- v
        va_ned = SVector{3, Float64}(va_enu[2], va_enu[1], -va_enu[3])
        yaw[k] = atan(va_ned[2], va_ned[1])
    end
    for k in 2:n
        dpsi = yaw[k] - yaw[k - 1]
        if dpsi > pi
            yaw[k] -= 2 * pi
        elseif dpsi < -pi
            yaw[k] += 2 * pi
        end
    end
    yaw_rate = diff(rad2deg.(yaw)) ./ diff(sl.time)

    us_cmd = steering_command(sl, sys)
    us_seg = us_cmd[2:end]
    v_app = sl.v_app[2:end]
    gk = similar(yaw_rate)
    @inbounds for k in eachindex(gk)
        gk[k] = abs(us_seg[k]) > 1e-8 ? yaw_rate[k] / (v_app[k] * us_seg[k]) : NaN
    end
    times = sl.time[2:end]
    return gk, times
end

function mean_last_window(values::AbstractVector{<:Real}, times::AbstractVector{<:Real};
                          window_sec::Real=WINDOW_SEC)
    @assert length(values) == length(times)
    t_end = times[end]
    mask = times .>= (t_end - window_sec)
    if !any(mask)
        mask = trues(length(times))
    end
    data = values[mask]
    data = data[isfinite.(data)]
    return isempty(data) ? NaN : mean(data)
end

function mean_at_time(values::AbstractVector{<:Real}, times::AbstractVector{<:Real},
                      target_time::Real; window_half::Real=0.5)
    @assert length(values) == length(times)
    mask = (times .>= (target_time - window_half)) .& (times .<= (target_time + window_half))
    if !any(mask)
        return NaN
    end
    data = values[mask]
    data = data[isfinite.(data)]
    return isempty(data) ? NaN : mean(data)
end

function analyze_log(lg, sys; window_sec::Real=WINDOW_SEC)
    sl = lg.syslog
    if length(sl.time) < 2
        return (
            aero_force=NaN, v_app=NaN, yaw_rate=NaN, yaw_rate_paper=NaN,
            gk=NaN, gk_paper=NaN, kite_vel=NaN, aoa=NaN, elevation=NaN, azimuth=NaN,
            cs=NaN
        )
    end

    aero_force_z = [sl.aero_force_b[i][3] for i in eachindex(sl.aero_force_b)]
    aero_force = mean_last_window(aero_force_z, sl.time; window_sec)

    v_app = mean_last_window(sl.v_app, sl.time; window_sec)

    yaw_rate_deg, yaw_rate_time = heading_rate(sl)
    yaw_rate = mean_last_window(yaw_rate_deg, yaw_rate_time; window_sec)

    ekf = compute_ekf_yaw_and_rate(lg, sys)
    if ekf === nothing
        yaw_rate_paper = yaw_rate
    else
        _, yaw_rate_ekf = ekf
        yaw_rate_paper = mean_last_window(yaw_rate_ekf, sl.time; window_sec)
    end

    gk_vals, gk_time = gk_series(sl, sys)
    gk = mean_last_window(gk_vals, gk_time; window_sec)

    gk_p_vals, gk_p_time = gk_paper_series(sl, sys)
    gk_paper = mean_last_window(gk_p_vals, gk_p_time; window_sec)

    v_kite_norm = [norm(v) for v in sl.vel_kite]
    kite_vel = mean_last_window(v_kite_norm, sl.time; window_sec)

    aoa_deg = rad2deg.(sl.AoA)
    aoa = mean_last_window(aoa_deg, sl.time; window_sec)

    elevation_deg = rad2deg.(sl.elevation)
    elevation = mean_last_window(elevation_deg, sl.time; window_sec)

    azimuth_deg = rad2deg.(sl.azimuth)
    azimuth = mean_last_window(azimuth_deg, sl.time; window_sec)

    cs_vals, cs_time = calculate_cs(sl, sys)
    cs = abs(mean_last_window(cs_vals, cs_time; window_sec))

    turn_radius_result = compute_turn_radius(sl, sys)
    if turn_radius_result === nothing
        turn_radius = NaN
    else
        turn_radius_vals, turn_radius_time = turn_radius_result
        turn_radius = abs(mean_last_window(turn_radius_vals, turn_radius_time; window_sec))
    end

    # Compute time-series snapshots at specific seconds
    us_cmd = steering_command(sl, sys)
    usva = us_cmd .* sl.v_app
    yaw_rate_deg, yaw_rate_time = heading_rate(sl)
    
    # Initialize time-series dictionaries
    usva_at = Dict{Int, Float64}()
    yaw_rate_at = Dict{Int, Float64}()
    va_at = Dict{Int, Float64}()
    aoa_at = Dict{Int, Float64}()
    
    for t_sec in 3:10
        usva_at[t_sec] = mean_at_time(usva, sl.time, Float64(t_sec))
        yaw_rate_at[t_sec] = mean_at_time(yaw_rate_deg, yaw_rate_time, Float64(t_sec))
        va_at[t_sec] = mean_at_time(sl.v_app, sl.time, Float64(t_sec))
        aoa_at[t_sec] = mean_at_time(aoa_deg, sl.time, Float64(t_sec))
    end

    return (
        aero_force=aero_force,
        v_app=v_app,
        yaw_rate=yaw_rate,
        yaw_rate_paper=yaw_rate_paper,
        gk=gk,
        gk_paper=gk_paper,
        kite_vel=kite_vel,
        aoa=aoa,
        elevation=elevation,
        azimuth=azimuth,
        cs=cs,
        turn_radius=turn_radius,
        usva_at=usva_at,
        yaw_rate_at=yaw_rate_at,
        va_at=va_at,
        aoa_at=aoa_at
    )
end

function find_log_names(batch_dir::AbstractString)
    isdir(batch_dir) || error("Batch folder not found: $batch_dir")
    files = readdir(batch_dir; join=true)
    names = String[]
    for file in files
        isfile(file) || continue
        endswith(file, ".txt") && continue
        name = splitext(basename(file))[1]
        parse_up_us_vw_lt(name) === nothing && continue
        push!(names, name)
    end
    return sort(unique(names))
end

function write_csv(path::AbstractString, rows)
    # Build header with time-series columns
    base_cols = "vw,up,us,lt,aero_force,v_app,yaw_rate,yaw_rate_paper,gk,gk_paper,kite_vel,aoa,elevation,azimuth,cs,turn_radius"
    time_cols = String[]
    for t in 3:10
        push!(time_cols, "usva_$t")
        push!(time_cols, "yaw_rate_$t")
        push!(time_cols, "va$(t)")
        push!(time_cols, "aoa$(t)")
    end
    header = base_cols * "," * join(time_cols, ",")
    
    open(path, "w") do io
        println(io, header)
        for r in rows
            base_vals = [
                r.vw, r.up, r.us, r.lt, r.aero_force, r.v_app, r.yaw_rate, r.yaw_rate_paper,
                r.gk, r.gk_paper, r.kite_vel, r.aoa, r.elevation, r.azimuth, r.cs, r.turn_radius
            ]
            time_vals = Float64[]
            for t in 3:10
                push!(time_vals, r.usva_at[t])
                push!(time_vals, r.yaw_rate_at[t])
                push!(time_vals, r.va_at[t])
                push!(time_vals, r.aoa_at[t])
            end
            println(io, join(vcat(base_vals, time_vals), ","))
        end
    end
end

function main()
    batch_name = isempty(ARGS) ? "" : strip(ARGS[1])
    batch_name = "circular_2019_batch_2026_01_10_12_01_04"
    if isempty(batch_name)
        print("Enter batch folder name (e.g. batch_2026_01_07_10_04_38): ")
        batch_name = strip(readline())
    end
    isempty(batch_name) && error("Batch folder name is required.")

    batch_dir = joinpath("processed_data", "v3_kite", batch_name)
    log_names = find_log_names(batch_dir)
    isempty(log_names) && error("No logs found in: $batch_dir")

    rows = NamedTuple[]
    sys_cache = Dict{Tuple{Int, Int}, SystemStructure}()

    for log_name in log_names
        tags = parse_up_us_vw_lt(log_name)
        tags === nothing && continue
        up, us, vw, lt = tags
        sys = get!(sys_cache, (vw, lt)) do
            build_sys(v_wind=float(vw), tether_length=float(lt))
        end
        lg = KiteUtils.load_log(log_name; path=batch_dir)
        metrics = analyze_log(lg, sys)
        push!(rows, (
            vw=vw,
            up=up,
            us=us,
            lt=lt,
            aero_force=metrics.aero_force,
            v_app=metrics.v_app,
            yaw_rate=metrics.yaw_rate,
            yaw_rate_paper=metrics.yaw_rate_paper,
            gk=metrics.gk,
            gk_paper=metrics.gk_paper,
            kite_vel=metrics.kite_vel,
            aoa=metrics.aoa,
            elevation=metrics.elevation,
            azimuth=metrics.azimuth,
            cs=metrics.cs,
            turn_radius=metrics.turn_radius,
            usva_at=metrics.usva_at,
            yaw_rate_at=metrics.yaw_rate_at,
            va_at=metrics.va_at,
            aoa_at=metrics.aoa_at
        ))
    end

    sort!(rows, by=r -> (r.vw, r.up, r.us, r.lt))

    out_path = joinpath(batch_dir, "circle_batch_analysis.csv")
    write_csv(out_path, rows)
    @info "Wrote batch analysis CSV" path=out_path rows=length(rows)
end

# if abspath(PROGRAM_FILE) == string(@__FILE__)
#     main()
# end

main()

nothing
