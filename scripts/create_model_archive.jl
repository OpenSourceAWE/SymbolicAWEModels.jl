# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

using Timers
@info "Loading packages"
tic()
using SymbolicAWEModels 
toc()

prn = true
@info "Creating default models"
models = SymbolicAWEModels.create_default_models(; prn)
@info "Created all models"
toc()

version = VERSION.minor
SymbolicAWEModels.create_model_archive(get_data_path(), 
                                       joinpath(get_data_path(), "models_v1.$version.tar.gz"))

