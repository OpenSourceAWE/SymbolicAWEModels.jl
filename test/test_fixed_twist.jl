# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# FIXED twist mode: twist is a prescribed control input (no differential state,
# no algebraic equilibrium). Verifies validate_twist_surface_modes and a rigid
# wing built with FIXED-twist twist_surfaces.

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

using Test
using SymbolicAWEModels
using SymbolicAWEModels: VortexStepMethod, validate_twist_surface_modes,
    Wing, TwistSurface
using KiteUtils: init!, next_step!
using LinearAlgebra

pkg_root = dirname(@__DIR__)
set_data_path(joinpath(pkg_root, "data", "2plate_kite"))

@testset "validate_twist_surface_modes" begin
    rigid = Wing(:rigid, NameRef[], Matrix{Float64}(I, 3, 3),
        zeros(3), ones(3); dynamics_type=RIGID_DYNAMICS)
    rigid.idx = 1
    particle = Wing(:particle, NameRef[], Matrix{Float64}(I, 3, 3),
        zeros(3), ones(3); dynamics_type=PARTICLE_DYNAMICS)
    particle.idx = 1

    mktwist_surface(name, npoints, type) = begin
        g = TwistSurface(name, collect(1:npoints), type, 0.25)
        g.idx = 1
        g.point_idxs = collect(1:npoints)
        g
    end

    # DYNAMIC on rigid + >=2 points -> ok
    rigid.twist_surface_idxs = [1]
    @test validate_twist_surface_modes([mktwist_surface(:g, 2, DYNAMIC)], [rigid]) === nothing
    # QUASI_STATIC on rigid + >=2 points -> ok
    @test validate_twist_surface_modes([mktwist_surface(:g, 2, QUASI_STATIC)], [rigid]) === nothing
    # FIXED on rigid, any point count -> ok
    @test validate_twist_surface_modes([mktwist_surface(:g, 1, FIXED)], [rigid]) === nothing
    @test validate_twist_surface_modes([mktwist_surface(:g, 3, FIXED)], [rigid]) === nothing

    # DYNAMIC on particle -> reject (needs rigid)
    particle.twist_surface_idxs = [1]
    @test_throws ErrorException validate_twist_surface_modes(
        [mktwist_surface(:g, 2, DYNAMIC)], [particle])
    # DYNAMIC 1-point -> reject (needs bridle couple)
    @test_throws ErrorException validate_twist_surface_modes(
        [mktwist_surface(:g, 1, DYNAMIC)], [rigid])
    # QUASI_STATIC 1-point -> reject
    @test_throws ErrorException validate_twist_surface_modes(
        [mktwist_surface(:g, 1, QUASI_STATIC)], [rigid])
    # FIXED on particle + multi-point -> reject
    @test_throws ErrorException validate_twist_surface_modes(
        [mktwist_surface(:g, 2, FIXED)], [particle])
    # FIXED on particle + 1 point -> ok
    @test validate_twist_surface_modes([mktwist_surface(:g, 1, FIXED)], [particle]) === nothing
end

@testset "FIXED twist on rigid VSM wing" begin
    # Work in a tmpdir copy so the temp FIXED-twist_surfaces geometry never touches the
    # repo data dir.
    data_dir = mktempdir()
    cp(joinpath(pkg_root, "data", "2plate_kite"), data_dir; force=true)
    set_data_path(data_dir)

    set = Settings("system.yaml")
    set.g_earth = 0.0
    vsm_set = VortexStepMethod.VSMSettings(
        joinpath(get_data_path(), "vsm_settings.yaml"); data_prefix=false)

    # Build a temp structural geometry with FIXED twist_surfaces.
    src_yaml = joinpath(get_data_path(), "rigid_structural_geometry.yaml")
    txt = read(src_yaml, String)
    txt = replace(txt,
        "[left, [le_left, te_left], DYNAMIC, 0.25, 100.0]" =>
            "[left, [le_left, te_left], FIXED, 0.25, 100.0]",
        "[center, [le_center, te_center], DYNAMIC, 0.25, 100.0]" =>
            "[center, [le_center, te_center], FIXED, 0.25, 100.0]",
        "[right, [le_right, te_right], DYNAMIC, 0.25, 100.0]" =>
            "[right, [le_right, te_right], FIXED, 0.25, 100.0]")
    fixed_yaml = joinpath(get_data_path(), "rigid_fixed_twist_geometry.yaml")
    write(fixed_yaml, txt)

    sys = load_sys_struct_from_yaml(fixed_yaml;
        system_name="2plate_fixed", set, vsm_set)
    sys.winches[:main_winch].brake = true
    @test all(g.type == FIXED for g in sys.twist_surfaces)

    # Prescribe distinct twist angles per twist_surface.
    prescribed = Dict(:left => 0.05, :center => -0.10, :right => 0.08)
    for g in sys.twist_surfaces
        g.twist = prescribed[g.name]
    end

    sam = SymbolicAWEModel(set, sys)
    init!(sam; prn=false, remake=true, remake_vsm=false)
    next_step!(sam)

    # twist_angle must track the prescribed value (no dynamics, no clamp drift).
    for g in sam.sys_struct.twist_surfaces
        @test isapprox(g.twist, prescribed[g.name]; atol=1e-9)
    end
    rm(fixed_yaml; force=true)
end
nothing
