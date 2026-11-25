# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Station twist modes: validate_station_modes, a rigid wing built with STATIC-twist
# stations (twist a prescribed control input), and the inertia a DYNAMIC station
# twists with.

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

using Test
using SymbolicAWEModels
using SymbolicAWEModels: VortexStepMethod, validate_station_modes, update_run_state!
using KiteUtils: init!, next_step!
using LinearAlgebra

pkg_root = dirname(@__DIR__)
set_data_path(joinpath(pkg_root, "data", "2plate_kite"))

"""The run of the smallest system with one station of `npoints` points and twist `type` on
a wing of `dynamics`, flown without aero."""
function station_check_model(npoints, type, dynamics)
    set = Settings("system.yaml")
    points = [Point(; name=Symbol(:p, i), pos_ENU=[0.0, Float64(i), 0.0])
              for i in 1:npoints]
    station = SpanwiseStation(; name=:g, wing=:wing,
                              points=[point.name for point in points],
                              model=StationModel(; type=string(type), moment_frac=0.25))
    wing = KiteGeometry.Wing(; name=:wing, pos_ENU=zeros(3), body=:wing,
        model=WingModel(; aero_mode="AeroNone", dynamics_type=string(dynamics),
                        stations=[:g],
                        origin_idx=dynamics == PARTICLE_DYNAMICS ? :p1 : nothing))
    body = Body(; name=:wing, pos_ENU=zeros(3), extra_mass=1.0,
                extra_inertia_KA=Matrix(Diagonal(ones(3))),
                model=BodyModel(; angular_damping=0.0))
    sys = system_definition("station_modes", set; points, spanwise_stations=[station],
                            wings=[wing], bodies=[body])
    return SymbolicAWEModel(set, sys)
end

"""Whether `validate_station_modes` accepts the station of
[`station_check_model`](@ref)."""
station_modes_valid(npoints, type, dynamics) =
    validate_station_modes(station_check_model(npoints, type, dynamics)) === nothing

@testset "validate_station_modes" begin
    # DYNAMIC on rigid + >=2 points -> ok
    @test station_modes_valid(2, DYNAMIC, RIGID_DYNAMICS)
    # STATIC on rigid, any point count -> ok
    @test station_modes_valid(1, STATIC, RIGID_DYNAMICS)
    @test station_modes_valid(3, STATIC, RIGID_DYNAMICS)
    # DYNAMIC on particle -> reject (needs rigid)
    @test_throws ErrorException station_modes_valid(2, DYNAMIC, PARTICLE_DYNAMICS)
    # DYNAMIC 1-point -> reject (needs bridle couple)
    @test_throws ErrorException station_modes_valid(1, DYNAMIC, RIGID_DYNAMICS)
    # STATIC on particle + multi-point -> ok (inert aero-section membership
    # marker; the free particles carry the deformation)
    @test station_modes_valid(2, STATIC, PARTICLE_DYNAMICS)
    # STATIC on particle + 1 point -> ok
    @test station_modes_valid(1, STATIC, PARTICLE_DYNAMICS)
end

@testset "STATIC twist on rigid VSM wing" begin
    # Work in a tmpdir copy so the temp STATIC-stations geometry never touches the
    # repo data dir.
    data_dir = mktempdir()
    cp(joinpath(pkg_root, "data", "2plate_kite"), data_dir; force=true)
    set_data_path(data_dir)

    set = Settings("system.yaml")
    set.g_earth = 0.0
    vsm_set = VortexStepMethod.VSMSettings(
        joinpath(get_data_path(), "vsm_settings.yaml"); data_prefix=false)

    # Build a temp structural geometry with STATIC stations.
    src_yaml = joinpath(get_data_path(), "rigid_structural_geometry.yaml")
    txt = read(src_yaml, String)
    txt = replace(txt,
        "[left, [le_left, te_left], DYNAMIC, 0.25, 100.0]" =>
            "[left, [le_left, te_left], STATIC, 0.25, 100.0]",
        "[center, [le_center, te_center], DYNAMIC, 0.25, 100.0]" =>
            "[center, [le_center, te_center], STATIC, 0.25, 100.0]",
        "[right, [le_right, te_right], DYNAMIC, 0.25, 100.0]" =>
            "[right, [le_right, te_right], STATIC, 0.25, 100.0]")
    fixed_yaml = joinpath(get_data_path(), "rigid_fixed_twist_geometry.yaml")
    write(fixed_yaml, txt)

    sys, transforms = load_authoring(fixed_yaml;
        name="2plate_fixed", set)
    sys.winches[:main_winch].brake = true
    @test all(station_type(g) == STATIC for g in sys.spanwise_stations)

    # Prescribe distinct twist angles per station.
    prescribed = Dict(:left => 0.05, :center => -0.10, :right => 0.08)
    for g in sys.spanwise_stations
        g.twist = prescribed[g.name]
    end

    sam = SymbolicAWEModel(set, sys; transforms, vsm_set)
    init!(sam; prn=false, remake=true, remake_vsm=false)
    next_step!(sam)

    # twist_angle must track the prescribed value (no dynamics, no clamp drift).
    for (station_idx, g) in enumerate(sys.spanwise_stations)
        @test isapprox(run_value(sam, :spanwise_stations, station_idx, :twist),
                       prescribed[g.name]; atol=1e-9)
    end
    rm(fixed_yaml; force=true)
end

"""The sum of point mass `quantity` over the points of station `station_idx` of `sam`."""
station_point_mass(sam, station_idx, quantity) =
    sum(run_value(sam, :points, i, quantity)
        for i in sam.definition.spanwise_stations[station_idx].points)

@testset "DYNAMIC station twists with the total mass of its points, $backend" for
        backend in (MonolithBackend(), KernelBackend())
    set_data_path(joinpath(pkg_root, "data", "2plate_kite"))
    set = Settings("system.yaml")
    set.g_earth = 0.0
    vsm_set = VortexStepMethod.VSMSettings(
        joinpath(get_data_path(), "vsm_settings.yaml"); data_prefix=false)
    sys, transforms = load_authoring(
        joinpath(get_data_path(), "rigid_structural_geometry.yaml");
        name="2plate_kite", set)
    sys.winches[:main_winch].brake = true
    sam = SymbolicAWEModel(set, sys; transforms, vsm_set, backend)
    init!(sam; prn=false, remake=false, remake_vsm=false)
    update_run_state!(sam.prob, sam.integrator, sam)

    stations = eachindex(sys.spanwise_stations)
    station(i, quantity) = run_value(sam, :spanwise_stations, i, quantity)
    @test all(station_point_mass(sam, i, :total_mass) >
              1.01 * station_point_mass(sam, i, :extra_mass)
              for i in stations)
    moment = [station(i, :aero_moment) + station(i, :tether_moment) for i in stations]
    inertia = [(station_point_mass(sam, i, :total_mass) + station(i, :body_mass)) *
               norm(station(i, :chord))^2 / 3 for i in stations]
    @test all(station(i, :twist_ω) == 0 for i in stations)

    # From rest the twist rate gains `moment / inertia · dt` to first order in `dt`.
    dt = 1e-6
    next_step!(sam; dt, vsm_interval=0)
    twist_acc = [station(i, :twist_ω) / dt for i in stations]
    @test twist_acc .* inertia ≈ moment rtol=2e-3
end
nothing
