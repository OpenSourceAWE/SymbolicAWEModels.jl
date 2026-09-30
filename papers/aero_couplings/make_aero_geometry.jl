# Copyright (c) 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# Copy V3Kite's data and the paper's particle projects into output/data, and write the
# NeuralFoil aero geometry both wings fly with, as V3Kite's
# examples/v3beam_aero_geometry.jl does for its own data directory.

using Pkg
Pkg.activate(@__DIR__)

using V3Kite
using SymbolicAWEModels
using VortexStepMethod
using VortexStepMethod.ObjAdapter: obj_to_yaml
using VortexStepMethod.AirfoilAero: ShrinkWrap, NeuralFoilSolver

include(joinpath(@__DIR__, "paper_setup.jl"))

REYNOLDS = 1.225 * 15.4 * 2.32 / 1.81e-5   # ρ v_app c_ref / μ at the flight condition

isdir(DATA_PATH) || cp(v3_data_path(), DATA_PATH)
chmod(DATA_PATH, 0o755; recursive=true)
for project in (PARTICLE, PARTICLE_LIVE)
    cp(joinpath(@__DIR__, project), joinpath(DATA_PATH, project); force=true)
end
obj_to_yaml(joinpath(pkgdir(VortexStepMethod), "data", "TUDELFT_V3_KITE", "V3_25.obj"),
    joinpath(DATA_PATH, "polars_neuralfoil");
    geometry_path=joinpath(DATA_PATH, "nf_aero_geometry.yaml"), n_sections=37,
    Re=REYNOLDS, alpha_range=-10:2:30, delta_range=nothing,
    aero_solver=NeuralFoilSolver(model_size="large"), wrap_method=ShrinkWrap(),
    wingtip_distance=0.05, crease_frac=V3BeamTopology().crease_frac,
    table_format=:arrow)
nothing
