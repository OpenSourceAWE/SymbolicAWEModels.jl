# Copyright (c) 2025 Uwe Fechner, Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

using PrecompileTools: @setup_workload, @compile_workload

"""
    workload_fixture() -> String

Copy the 2-plate fixture into a scratch directory and return its path. The
workload builds real models and `init!` serialises one, so it must not write into
the package's own `data`.
"""
function workload_fixture()
    source = joinpath(@__DIR__, "..", "data", "2plate_kite")
    destination = mktempdir()
    for name in ("system.yaml", "settings.yaml", "vsm_settings.yaml",
                 "aero_geometry.yaml", "particle_structural_geometry.yaml",
                 "rigid_structural_geometry.yaml")
        cp(joinpath(source, name), joinpath(destination, name))
    end
    cp(joinpath(source, "polars"), joinpath(destination, "polars"))
    return destination
end

"""
    workload_model(fixture, geometry, system_name, backend) -> SymbolicAWEModel

Load the 2-plate structure from `geometry` and build it on `backend`, ready for
`init!`. Reloads the structure per call because building consumes it.
"""
function workload_model(fixture, geometry, system_name, backend)
    set = Settings("system.yaml")
    set.g_earth = 0.0
    vsm_set = VortexStepMethod.VSMSettings(joinpath(fixture, "vsm_settings.yaml");
                                           data_prefix=false)
    sys = load_sys_struct_from_yaml(joinpath(fixture, geometry);
                                    system_name, set, vsm_set, prn=false)
    haskey(sys.winches, :main_winch) && (sys.winches[:main_winch].brake = true)
    return SymbolicAWEModel(set, sys; backend)
end

"""
    run_workload(fixture)

Build, initialise and step the 2-plate models on both backends. Called from this
package's workload and again from the Makie extension's, where the same code is
compiled against the method table Makie leaves behind — loading Makie invalidates
around 19500 method instances, and a cold start is almost entirely inference.
"""
function run_workload(fixture)
    set_data_path(fixture)
    geometries = (("particle_structural_geometry.yaml", "workload_particle"),
                  ("rigid_structural_geometry.yaml", "workload_rigid"))
    for (geometry, system_name) in geometries, backend in (KernelBackend(),
                                                           MonolithBackend())
        sam = workload_model(fixture, geometry, system_name, backend)
        init!(sam; remake=false, remake_vsm=false, prn=false)
        next_step!(sam; dt=0.05)
    end
    return nothing
end

@setup_workload begin
    fixture = workload_fixture()
    @compile_workload begin
        run_workload(fixture)
    end
end
