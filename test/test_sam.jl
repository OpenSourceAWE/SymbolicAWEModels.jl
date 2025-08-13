# SPDX-FileCopyrightText: 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

using Pkg
if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end
using Test, ControlSystemsBase, Printf
using SymbolicAWEModels, ControlPlots
using Statistics, LinearAlgebra, Serialization
using ModelingToolkit: @variables
using ModelingToolkit: t_nounits

tmpdir=mktempdir()
mkpath(joinpath(tmpdir, "data"))
old_data_path = get_data_path()
set_data_path(joinpath(tmpdir, "data"))
cp(old_data_path, get_data_path(); force=true)

set = Settings("system.yaml")
sam = SymbolicAWEModel(set, "ram")

tether_set = Settings("system.yaml")
tether_sam = SymbolicAWEModel(tether_set, "tether")
init!(tether_sam)

simple_set = Settings("system.yaml")
simple_sam = SymbolicAWEModel(simple_set, "simple_ram")
init!(simple_sam)

original_set = Settings("system.yaml")

function reset!(set::Settings)
    for field in fieldnames(Settings)
        setfield!(set, field, getfield(original_set, field))
    end
    return set
end

@testset verbose=true "SymbolicAWEModels Tests" begin

    function test_plot(sam)
        @testset "Plotting of SymbolicAWEModel" begin
            function plot_(zoom, front)
                plt.figure("Kite")
                lines, sc, txt = plot(sam, 0.0; zoom, front)
                plt.show(block=false)
                sleep(1)
                @test !isnothing(lines)
                @test length(lines) ≥ 1  # Should have at least one line
                @test !isnothing(sc)     # Should have scatter points
                @test !isnothing(txt)    # Should have time text
            end
            plot_(false, false)
            plot_(false, true)
            plot_(true, false)
            plot_(true, true)
        end
    end

    function init(elevation, azimuth, heading)
        set.elevation = elevation
        set.azimuth = azimuth
        set.heading = heading
        init!(sam)
        ss = SysState(sam)
        @test sam.sys_struct.wings[1].elevation ≈ deg2rad(set.elevation) atol=1e-2
        @test sam.sys_struct.wings[1].azimuth ≈ deg2rad(set.azimuth) atol=1e-2
        @test sam.sys_struct.wings[1].heading ≈ deg2rad(set.heading) atol=1e-2
        @test ss.elevation ≈ deg2rad(set.elevation) atol=1e-2
        @test ss.azimuth ≈ deg2rad(set.azimuth) atol=1e-2
        @test ss.heading ≈ deg2rad(set.heading) atol=1e-2
    end
    
    @testset "Initialization" begin
        @test sam isa SymbolicAWEModel
        @test sam isa AbstractKiteModel
        @test simple_sam isa SymbolicAWEModel
        @test simple_sam isa AbstractKiteModel
        @test tether_sam isa SymbolicAWEModel
        @test tether_sam isa AbstractKiteModel
        
        init_time = @elapsed init!(sam; prn=true)
        @show init_time
        @test init_time < 700
        init!(sam; prn=true)
        init_time = @elapsed init!(sam; prn=true)
        @test init_time < 0.3

        init(ones(3)...)
        init(zeros(3)...)
        test_plot(sam)
    end

    @testset "Tether properties" begin
        reset!(set)
        set.sample_freq = 600
        set.abs_tol = 1e-6
        set.rel_tol = 1e-6
        set.segments = 1
        one_seg_sam = SymbolicAWEModel(set, "ram")
        init!(one_seg_sam)
        one_seg_tether_sam = SymbolicAWEModel(set, "tether")
        init!(one_seg_tether_sam)

        # run twice to make sure state is reset properly
        axial_stiffness, axial_damping = 
            SymbolicAWEModels.calc_spring_props(one_seg_sam, one_seg_tether_sam)
        next_step!(one_seg_sam; dt=1.0)
        axial_stiffness, axial_damping = 
            SymbolicAWEModels.calc_spring_props(one_seg_sam, one_seg_tether_sam)
        segments = one_seg_sam.sys_struct.segments
        tethers = one_seg_sam.sys_struct.tethers
        segments = [segments[tether.segment_idxs[1]] for tether in tethers]
        real_axial_stiffness = [segment.axial_stiffness for segment in segments]
        real_axial_damping = [segment.axial_damping for segment in segments]
        @test isapprox(real_axial_stiffness, axial_stiffness; rtol=0.02)
        @test isapprox(real_axial_damping, axial_damping; rtol=0.2)

        println("\n--- Tether Spring Properties ---")
        # Print table headers
        @printf "%-8s | %-15s %-15s %-10s | %-15s %-15s %-10s\n" "Tether" "Calc. Stiffness" "Real Stiffness" "Error (%)" "Calc. Damping" "Real Damping" "Error (%)"
        # Print separator line
        println(repeat("-", 100))
        for i in 1:4
            # Calculate relative errors in percent
            stiffness_err = 100 * abs(axial_stiffness[i] - real_axial_stiffness[i]) / real_axial_stiffness[i]
            damping_err   = 100 * abs(axial_damping[i] - real_axial_damping[i]) / real_axial_damping[i]
            # Print data rows
            @printf "%-8d | %-15.2f %-15.2f %-10.2f | %-15.2f %-15.2f %-10.2f\n" i axial_stiffness[i] real_axial_stiffness[i] stiffness_err axial_damping[i] real_axial_damping[i] damping_err
        end
        println()
    end

    @testset "Oscillating simulation" begin
        function test_for_peak_at_steering_freq(sam, steering_freq)
            dt = 0.01
            sl, _ = sim_oscillate!(sam; total_time=5.0, steering_freq, dt)
            @test sl.syslog.elevation[begin] ≈ deg2rad(set.elevation) atol=1e-2
            @test sl.syslog.azimuth[begin] ≈ deg2rad(set.azimuth) atol=1e-2
            @test sl.syslog.heading[begin] ≈ deg2rad(set.heading) atol=1e-2
            @test isapprox(sl.syslog.time, collect(0.0:dt:5.0-dt))
            ControlPlots.plt.close_figs()
            plt = plot(sam.sys_struct, sl)
            display(plt)
            @test plt isa ControlPlots.PlotX
            savefig(joinpath(
                get_data_path(), "oscillate_$(sam.sys_struct.name)_$steering_freq.png"
            ))

            # --- Cross-Correlation Analysis with Linear Offset Removal (first/last only) ---
            heading_signal = sl.syslog.heading
            t = sl.syslog.time

            # Compute linear trend using endpoints
            trend = range(heading_signal[1], heading_signal[end], length=length(t))
            signal_detrended = heading_signal .- trend

            # Reference sine and cosine at steering frequency
            ref_sin = sin.(2 * π * steering_freq .* t)
            ref_cos = cos.(2 * π * steering_freq .* t)

            # Correlation at steering frequency
            corr_sin = dot(signal_detrended, ref_sin)
            corr_cos = dot(signal_detrended, ref_cos)
            magnitude_at_freq = sqrt(corr_sin^2 + corr_cos^2)

            # Compare to frequencies 50% lower and 50% higher
            freq_lower = steering_freq * 0.5
            ref_sin_lower = sin.(2 * π * freq_lower .* t)
            ref_cos_lower = cos.(2 * π * freq_lower .* t)
            mag_lower = sqrt(dot(signal_detrended, ref_sin_lower)^2 + dot(signal_detrended, ref_cos_lower)^2)

            freq_higher = steering_freq * 1.5
            ref_sin_higher = sin.(2 * π * freq_higher .* t)
            ref_cos_higher = cos.(2 * π * freq_higher .* t)
            mag_higher = sqrt(dot(signal_detrended, ref_sin_higher)^2 + dot(signal_detrended, ref_cos_higher)^2)

            @show magnitude_at_freq, mag_lower, mag_higher
            @test magnitude_at_freq > mag_lower && magnitude_at_freq > mag_higher
        end

        reset!(set)

        init!(sam)
        find_steady_state!(sam)
        test_for_peak_at_steering_freq(sam, 0.5)

        init!(sam)
        find_steady_state!(sam)
        SymbolicAWEModels.copy_to_simple!(sam, tether_sam, simple_sam)
        test_for_peak_at_steering_freq(simple_sam, 0.5)
    end

    @testset "Turning simulation" begin
        function unwrap!(v::AbstractVector, period::Real=2π)
            offset = 0.0
            for i in 2:length(v)
                diff = v[i] - v[i-1]
                if diff > period / 2
                    offset -= period
                elseif diff < -period / 2
                    offset += period
                end
                v[i] += offset
            end
            return v
        end
        function calc_heading(steering_time, steering_magnitude)
            reset!(set)
            init!(sam)
            find_steady_state!(sam)
            dt = 0.05
            sl, _ = sim_turn!(sam; total_time=10.0, steering_time, steering_magnitude, dt)
            unwrap!(sl.syslog.heading)
            @test sl.syslog.heading[begin] ≈ 0.0 atol=1e-1
            return sl.syslog.heading[end]
        end
        default_heading = calc_heading(1.0, 10.0)
        @test default_heading ≈ 1035 atol=10.0
        short_steer_heading = calc_heading(0.5, 10.0)
        soft_steer_heading = calc_heading(1.0, 5.0)
        # make sure less steering results in less final heading
        @test default_heading - short_steer_heading ≈ 93 atol=10.0
        @test default_heading - soft_steer_heading ≈ 150 atol=10.0
        @show default_heading, short_steer_heading, soft_steer_heading
    end

    @testset "Linearize" begin
        old_abs = set.abs_tol
        old_rel = set.rel_tol
        set.abs_tol = 1e-4
        set.rel_tol = 1e-4
        init!(sam)
        init!(simple_sam)

        (; A, B, C, D) = SymbolicAWEModels.linearize!(simple_sam)
        sys = ss(A,B,C,D)
        norm_A = norm(A)
        res = lsim(sys, repeat([-1.0 0.0 -1.0], 2)', [0.0, 0.5])
        println(res.y[:,2])
        @test isapprox(res.y[:,2], 
            [-0.0008037289321365251, 0.0004562826732837309, -0.020711457720341487, 
                        -0.0017333135190197818], rtol=0.1)

        find_steady_state!(sam; dt=3.0, t=10.0)
        (; A, B, C, D) = SymbolicAWEModels.simple_linearize!(sam; tstab=1.0)
        sys = ss(A,B,C,D)
        res = lsim(sys, repeat([-1.0 0.0 -1.0], 2)', [0.0, 0.5])
        println(res.y[:,2])
        @test isapprox(res.y[:,2],
            [0.015575316961016356, -0.0001989661253600774, -0.017933805715950355, 
                        6.679990358160092], rtol=0.1)

        # test that linearization is state-dependent
        next_step!(simple_sam; dt=1.0)
        (; A, B, C, D) = SymbolicAWEModels.linearize!(simple_sam)
        @test !isapprox(norm(A), norm_A; atol=1e-3)

        set.abs_tol = old_abs
        set.rel_tol = old_rel
    end

    @testset "Just a tether, without winch or kite" begin
        set.segments = 20
        dynamics_type = DYNAMIC

        points = Point[]
        segments = Segment[]

        points = push!(points, Point(1, zeros(3), STATIC; wing_idx=0))

        segment_idxs = Int[]
        for i in 1:set.segments
            point_idx = i+1
            pos = [0.0, 0.0, i * set.l_tether / set.segments]
            push!(points, Point(point_idx, pos, dynamics_type; wing_idx=0))
            segment_idx = i
            push!(segments, Segment(segment_idx, set, (point_idx-1, point_idx), BRIDLE))
            push!(segment_idxs, segment_idx)
        end

        transforms = [Transform(1, deg2rad(-80), 0.0, 0.0; 
            base_pos=[0.0, 0.0, 50.0], base_point_idx=points[1].idx, rot_point_idx=points[end].idx)]
        sys_struct = SymbolicAWEModels.SystemStructure("tether", set; points, segments, transforms)

        sam = SymbolicAWEModel(set, sys_struct)
        init!(sam; remake=false)
        sys = sam.prob.sys
        @test isapprox(sam.integrator[sys.pos[:, end]], [8.682408883346524, 0.0, 0.7596123493895988], atol=1e-2)
        for i in 1:100
            next_step!(sam)
        end
        @test sam.integrator[sys.pos[1, end]] > 0.8set.l_tether
        @test isapprox(sam.integrator[sys.pos[2, end]], 0.0, atol=1.0)
        test_plot(sam)
    end
    
    @testset "Serialization and Deserialization" begin
        points = [
            Point(1, zeros(3), DYNAMIC; wing_idx=0, transform_idx=1)
            Point(2, ones(3), DYNAMIC; wing_idx=0, transform_idx=1)
        ]
        segments = [
            Segment(1, set, (1,2), BRIDLE)
        ]
        transforms = [Transform(1, zeros(3)...; 
            base_pos=zeros(3), base_point_idx=1, rot_point_idx=1)]
        sys_struct = SystemStructure("one_point", set; points, segments, transforms)
        sam = SymbolicAWEModel(set, sys_struct)
        model_path = joinpath(get_data_path(), SymbolicAWEModels.get_model_name(set))

        function test_init_with_reset(create_prob, create_lin_prob, create_control_func)
            println("Create prob: $create_prob \t"*
                    "lin_prob: $create_lin_prob \t"*
                    "control_func: $create_control_func")
            rm(model_path; force=true)
            sam = SymbolicAWEModel(set, sys_struct)
            init!(sam; create_prob, create_lin_prob, create_control_func, prn=false)
            @test isnothing(sam.prob) == !create_prob
            @test isnothing(sam.lin_prob) == !create_lin_prob
            @test isnothing(sam.control_funcs) == !create_control_func
        end
        test_init_with_reset(false, false, false)
        test_init_with_reset(true, false, false)
        test_init_with_reset(false, true, false)
        test_init_with_reset(false, false, true)

        init!(sam; create_prob=true, create_lin_prob=true, create_control_func=true, prn=false)
        # check if not removed
        init!(sam; create_prob=false, create_lin_prob=false, create_control_func=false, prn=false)
        @test !isnothing(sam.prob)
        @test !isnothing(sam.lin_prob)
        @test !isnothing(sam.control_funcs)

        # same name, check if hash works
        push!(points, Point(3, [0, 0, 2], DYNAMIC; wing_idx=0, transform_idx=1))        
        sys_struct2 = SystemStructure("one_point", set; points, segments, transforms)
        sam.sys_struct = sys_struct2
        model_path2 = joinpath(get_data_path(), SymbolicAWEModels.get_model_name(set))
        @test model_path == model_path2
        # should create nothing because hashes are broken
        init!(sam; create_prob=false, create_lin_prob=false, create_control_func=false, prn=false)
        @test isnothing(sam.prob)
        @test isnothing(sam.lin_prob)
        @test isnothing(sam.control_funcs)
        init!(sam; create_prob=true, create_lin_prob=true, create_control_func=true, prn=false)
        @test !isnothing(sam.prob)
        @test !isnothing(sam.lin_prob)
        @test !isnothing(sam.control_funcs)

        # changing from 0 to 1 output
        old_ny = length(sam.outputs)
        outputs = [sam.prob.sys.pos[1,1]]
        init!(sam; outputs)
        @test old_ny == 0
        @test length(sam.outputs) == 1
        lin_model = linearize!(sam)
        @test size(lin_model.C)[1] == 1

        # changing from 1 to 2 outputs
        old_ny = length(sam.outputs)
        outputs = [sam.prob.sys.pos[1,1], sam.prob.sys.pos[2,1]]
        init!(sam; outputs)
        @test old_ny == 1
        @test length(sam.outputs) == 2
        lin_model = linearize!(sam)
        @test size(lin_model.C)[1] == 2
    end
end
set_data_path(old_data_path)
nothing

