<!--
SPDX-FileCopyrightText: 2025 Uwe Fechner

SPDX-License-Identifier: MIT
-->

# SymbolicAWEModels

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://OpenSourceAWE.github.io/SymbolicAWEModels.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://OpenSourceAWE.github.io/SymbolicAWEModels.jl/dev)
[![CI](https://github.com/OpenSourceAWE/SymbolicAWEModels.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/OpenSourceAWE/SymbolicAWEModels.jl/actions/workflows/CI.yml)
[![Coverage](https://codecov.io/gh/OpenSourceAWE/SymbolicAWEModels.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/OpenSourceAWE/SymbolicAWEModels.jl)
[![DOI](https://zenodo.org/badge/443855286.svg)](https://zenodo.org/doi/10.5281/zenodo.13310253)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

## Airborne Wind Energy (AWE) system models

This package provides modular symbolic models of Airborne Wind Energy (AWE) systems, 
which consist of a wing, one or more tethers, one or more winches and a bridle system with or without pulleys.
The kite is modeled as a deforming rigid body with orientation governed by quaternion dynamics. The aerodynamic forces and moments are computed using the Vortex Step Method. The tether is modeled as point masses connected by spring-damper elements, with aerodynamic drag modeled realistically. 
The winch is modeled as a motor/generator that can reel in or out the tethers.

The [`SymbolicAWEModel`](@ref) has the following subcomponents, implemented in separate packages:
- AtmosphericModel from [AtmosphericModels](https://github.com/aenarete/AtmosphericModels.jl)
- WinchModel from [WinchModels](https://github.com/aenarete/WinchModels.jl) 
- The aerodynamic forces and moments of some of the models are calculated using the package [VortexStepMethod](https://github.com/Albatross-Kite-Transport/VortexStepMethod.jl)

This package is part of [`SymbolicAWEModels`](https://github.com/OpenSourceAWE/SymbolicAWEModels.jl),
which in turn is part of the Julia Kite Power Tools, which consist of the following packages:
<p align="center"><img src="https://github.com/OpenSourceAWE/SymbolicAWEModels.jl/blob/main/docs/src/kite_power_tools.png" width="500" /></p>

## News
#### Work in progress

## Installation
If possible, install [Julia 1.11](https://OpenSourceAWE.github.io/2024/08/09/installing-julia-with-juliaup.html), if you haven't already. Julia 1.10 is still supported, but the performance is worse. On Linux, make sure that Python3 and Matplotlib are installed:
```
sudo apt install python3-matplotlib
```

Before installing this software it is suggested to create a new project, for example like this:
```bash
mkdir test
cd test
julia --project="."
```
Then add SymbolicAWEModels from  Julia's package manager, by typing:
```julia
using Pkg
pkg"add SymbolicAWEModels"
``` 
at the Julia prompt. You can run the unit tests with the command (careful, can take 60 min):
```julia
pkg"test SymbolicAWEModels"
```
You can copy the examples to your project with:
```julia
using SymbolicAWEModels
SymbolicAWEModels.install_examples()
```
This also adds the extra packages, needed for the examples to the project. Furthermore, it creates a folder `data`
with some example input files. You can now run the examples with the command:
```julia
include("examples/menu.jl")
```
You can also run the ram-air-kite example like this:
```julia
include("examples/ram_air_kite.jl")
```
This might take two minutes. To speed up the model initialization, you can create a system image:
```bash
cd bin
./create_sys_image
```
If you now launch Julia with `./bin/run_julia` and then run the above example again, it should be about three
times faster.

## Advanced installation
If you intend to modify or extend the code, it is suggested to fork the `SymbolicAWEModels.jl` repository and to check out your fork:
```bash
git clone https://github.com/USERNAME/SymbolicAWEModels.jl
```
where USERNAME is your github username.
Then compile a system image:
```bash
cd SymbolicAWEModels.jl/bin
./create_sys_image
```
If you now launch julia with:
```bash
cd ..
./bin/run_julia
```
You can run the examples with:
```julia
menu()
```
You can also run the ram-air-kite example like this:
```julia
include("examples/ram_air_kite.jl")
```

## Ram air kite model
This model represents the kite as a deforming rigid body, with orientation governed by quaternion dynamics. Aerodynamics are computed using the Vortex Step Method. The kite is controlled from the ground via four tethers.

## Tether
The tether is modeled as point masses, connected by spring-damper elements. Aerodynamic drag is modeled realistically. When reeling out or in the unstreched length of the spring-damper elements
is varied. This does not translate into physics directly, but it avoids adding point masses at run-time, which would be even worse because it would introduce discontinuities. When using
Dyneema or similar high-strength materials for the tether the resulting system is very stiff which is a challenge for the solver.

## Further reading

## Replaying log files
If you want to replay old flight log files in 2D and 3D to understand and explain better how kite power systems work, please have a look at [KiteViewer](https://github.com/ufechner7/KiteViewer) . How new log files can be created and replayed is explained in the documentation of [KiteSimulators](https://github.com/aenarete/KiteSimulators.jl) .

## Licence
This project is licensed under the MPL-2.0 License.

## Copyright notice
See the copyright notices in the source files, and the list of authors in [AUTHORS.md](AUTHORS.md).

## See also
- [Research Fechner](https://research.tudelft.nl/en/publications/?search=Fechner+wind&pageSize=50&ordering=rating&descending=true) for the scientic background of the winches and tethers.
- More kite models [KiteModels](https://github.com/ufechner7/KiteModels.jl)
- The meta-package [KiteSimulators](https://github.com/aenarete/KiteSimulators.jl)
- the package [KiteUtils](https://github.com/OpenSourceAWE/KiteUtils.jl)
- the packages [WinchModels](https://github.com/aenarete/WinchModels.jl) and [KitePodModels](https://github.com/aenarete/KitePodModels.jl) and [AtmosphericModels](https://github.com/aenarete/AtmosphericModels.jl)
- the packages [KiteControllers](https://github.com/aenarete/KiteControllers.jl) and [KiteViewers](https://github.com/aenarete/KiteViewers.jl)
- the [VortexStepMethod](https://github.com/Albatross-Kite-Transport/VortexStepMethod.jl)

Authors: Bart van de Lint (bart@vandelint.net), Uwe Fechner (uwe.fechner.msc@gmail.com)

**Documentation** [Stable Version](https://OpenSourceAWE.github.io/SymbolicAWEModels.jl/stable) --- [Development Version](https://OpenSourceAWE.github.io/SymbolicAWEModels.jl/dev)
