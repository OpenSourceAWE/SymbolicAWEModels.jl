```@meta
CurrentModule = SymbolicAWEModels
```

# Developer guide

This guide provides instructions and best practices for developers contributing
to `SymbolicAWEModels.jl`.

---

## Prerequisites

Before you begin, ensure you have the following software installed:

- **Julia**: Latest release version. Install on Linux using
  [juliaup](https://github.com/JuliaLang/juliaup):

  ```bash
  curl -fsSL https://install.julialang.org | sh
  juliaup add release
  juliaup default release
  ```

- **Git**: For version control.
- **Bash**: A Unix-like shell environment.
- **Code editor**: Your preferred code editor with Julia support.

For Windows or macOS, check [these](https://ufechner7.github.io/2024/08/09/installing-julia-with-juliaup.html) instructions.

---

## Getting started: development workflow

Follow these steps to set up your local development environment:

**Fork the Repository**\
[Fork](https://github.com/OpenSourceAWE/SymbolicAWEModels.jl/fork) the
`SymbolicAWEModels.jl` repository on GitHub to create your own copy.

**Clone Your Fork**\
Clone your forked repository to your local machine. Replace `<UserName>` with
your GitHub username.

```bash
git clone https://github.com/<UserName>/SymbolicAWEModels.jl
```

**Configure the Upstream Remote** Add the original `OpenSourceAWE` repository as
a remote named `upstream`. This allows you to pull in the latest changes from
the main project.

```bash
cd SymbolicAWEModels.jl
git remote add upstream https://github.com/OpenSourceAWE/SymbolicAWEModels.jl
```

**Install and precompile the packages**

```bash
cd bin
./install
```

If you have the time, also create a system image, which contains all packages but `SymbolicAWEModels.jl` itself. This
has the advantage of a much lower startup time and the disadvantage that you need to recreate the system image after
updating packages. On a laptop with an `AMD 7840U` CPU and 32 GB RAM on battery power this takes at least 15 minutes.


```bash
cd bin
./create_sys_image
```

This requires at least 48 GB memory. If you have 16GB RAM, create a swap file with 32 GB. Also 
close all other programs before creating the system image to avoid an out-of-memory error. 
On macOS this is handled automatically.

**Start Julia**

Always start Julia with

```bash
./bin/run_julia
```

or with
```bash
jl
```

The second form requires that the line:

```bash
alias jl='./bin/run_julia'
```

in your `.bashrc` file in your home directory (Linux and Windows). For Mac, add this line to the `.zshrc` file.

This has a few advantages:

- It will activate the current project
- it will set the required number of threads of the garbage collector
- it will use a system image if available
- it will provide the function `menu()` to launch any of the examples without the need to type the longish `include(...)` command.

## Contributing code: branches and pull requests

To contribute your changes, please follow this standard Git workflow:

**Sync with the Main Project** Before starting new work, fetch the latest
changes from the `upstream` repository and update your local `main` branch. This
helps prevent merge conflicts.

```bash
git fetch upstream
git checkout main
git rebase upstream/main
```
If rebase fails, you can also use the `git merge` command instead.

**Keep Your Feature Branch Up to Date** While working on your feature branch,
regularly rebase onto the latest changes from `main` to avoid conflicts later:

```bash
git fetch upstream
git checkout main
git rebase upstream/main
git checkout add_lei_model
git rebase main
```

This is especially important for long-running feature branches. Rebasing
frequently makes conflicts smaller and easier to resolve.

**Create a Feature Branch** Create a new branch from your up-to-date `main`
branch. Give it a short, descriptive name that summarizes your change.

```bash
# Create and switch to your new branch
git checkout -b add_lei_model
```

Good branch names include `add_lei_model`, `improve_plot_recipe`, or
`fix_winch_dynamics`.

**Make and Commit Your Changes** Work on your feature and commit your changes as
you go. Write clear and concise commit messages.

```bash
git add .
git commit -m "Add initial structure for LEI kite model"
```

**Push to Your Fork** Push your new branch to your forked repository on GitHub.

```bash
git push -u origin add_lei_model
```

**Create a Pull Request** Go to the GitHub page for your fork. You should see a
prompt to create a pull request from your new branch. Create a pull request that
targets the `main` branch of the original `OpenSourceAWE/SymbolicAWEModels.jl`
repository. Provide a clear title and a detailed description of your changes.

---

## Improving the development experience

### Use Revise.jl for faster workflow

We recommend adding
**[Revise.jl](https://timholy.github.io/Revise.jl/stable/)** to your global
Julia environment. It allows you to modify source code without restarting your
Julia session, which is essential for efficient development.

**Install Revise.jl globally:**

```julia
# Start Julia without a project
julia

# In the REPL
using Pkg
pkg"add Revise"
```

**Configure Revise to auto-load on startup:**

Create or edit `~/.julia/config/startup.jl` (on Linux/Mac) or
`%USERPROFILE%\.julia\config\startup.jl` (on Windows):

```julia
try
    @eval using Revise
catch e
    @warn "Error initializing Revise" exception=(e, catch_backtrace())
end
```

This will automatically load Revise every time you start Julia. The
`try`/`catch` block ensures Julia will still start even if Revise encounters an
issue.

**Verify it works:**

Start a new Julia session and you should see Revise load automatically. You can
verify by checking:

```julia
julia> @which Revise
```

Now any changes you make to package source code will be automatically reflected
in your Julia session!

### Running examples during development

When developing the package, you'll want to test your changes with the examples.
Here's how to set up the examples to use your local development version:

#### Launching Julia

1. **From the package root directory**:

   ```bash
   jl
   ```

#### Running examples

Now any changes you make to the source code will be immediately reflected when
you run the examples (thanks to Revise.jl):

```julia
include("examples/coupled_2plate_kite.jl")
include("examples/menu.jl")
```

The `examples/Project.toml` file already contains the necessary dependencies:

- `GLMakie` - for visualization
- `KiteUtils` - for utility functions
- `SymbolicAWEModels` - the package itself

The `examples` project gets automatically activated when you run one of the examples. You can also just type `menu()` to get a menu with the examples.

#### Managing package dependencies

**Understanding the Package Manager:**

Press `]` in the Julia REPL to enter package manager (Pkg) mode. The prompt
changes to show your current project:

```julia
julia> ]  # Press ] to enter Pkg mode
(examples) pkg>  # Prompt shows you're in the examples project
```

Press backspace to exit Pkg mode and return to the Julia REPL.

**Common Pkg commands:**

- `add PackageName` - Add a package to the current project
- `rm PackageName` - Remove a package
- `dev .` or `dev ..` - Use local source code instead of registered version
- `st` - Show status (list all packages and their versions)
- `up` - Update all packages
- `instantiate` - Install all packages from Project.toml
- `resolve` - Resolve possible conflicts. This can fail. If it fails, you have to disable the system image (delete it or rename it) and delete the `Manifest-v1.xx.toml` file of the active Julia version. When you now run instantiate or resolve, a new `Manifest.toml` will be created. Rename it manually to `Manifest-v1.xx.toml` with `xx` being your minor Julia version number.

**Adding packages to the examples:**

```bash
# Start Julia
jl
```
Use the package manager to activate the examples project and add your package:
```

]  # Enter Pkg mode - prompt shows (examples) pkg>
activate examples
add YourPackage
st  # Verify the package was added
```

**Adding packages to SymbolicAWEModels itself:**
```bash
# Start Julia
jl
```
Use the package manager to add your package:
```julia
]  # Enter Pkg mode - prompt shows (SymbolicAWEModels) pkg>
add YourPackage
st  # Verify the package was added
```

The prompt `(ProjectName) pkg>` always tells you which project you're modifying.

### Building documentation locally

To preview documentation changes as you work:

#### Using LiveServer (recommended)

1. **Start Julia**:

   ```bash
   jl
   ```

2. **Build the docs and show them with live reload**:

   ```julia
   include("scripts/build_docu.jl")
   ```

   This will:
   - Generate documentation figures, if needed
   - Build the documentation
   - Open it in your default browser
   - Watch for changes to documentation files
   - Automatically rebuild and refresh when you save changes

#### Manual build

Alternatively, you can build the documentation once without the live server:

```bash
jl
```

```julia
include("docs/make.jl")
```

Then open `docs/build/index.html` in your browser.

**Note**: If you make changes to the package source code (not just
documentation), you'll need to reload Julia or use Revise.jl for the changes to
be reflected in the built documentation.

---

## Testing

The test suite is designed around **component isolation**: each test file
builds a minimal model from constructors (no YAML, no full kite) and
verifies the physics of a single component against analytical solutions.
This proves that the underlying dynamics are physically correct — for
example, that angular momentum is conserved, that terminal velocity
matches the analytical prediction, and that spring-damper forces follow
the expected constitutive law.

### Running tests

```bash
# Run the full test suite
jl -e 'using Pkg; Pkg.test()'

# Run a single test file
jl test/test_point.jl
jl test/test_segment.jl
```

### Test files

| Test file | Component | What it verifies |
|-----------|-----------|------------------|
| `test_point` | [`Point`](@ref) | Gravity free-fall, damping, drag terminal velocity |
| `test_segment` | [`Segment`](@ref) | Spring-damper forces, stiffness, drag |
| `test_wing` | [`Wing`](@ref AbstractWing) | RIGID_DYNAMICS and PARTICLE_DYNAMICS construction, VSM coupling |
| `test_wing_dynamics` | [`Wing`](@ref AbstractWing) | Torque response, precession, angular momentum conservation |
| `test_tether_winch` | [`Tether`](@ref), [`Winch`](@ref) | Reel-out, Coulomb/viscous friction, terminal velocity |
| `test_pulley` | [`Pulley`](@ref) | Equal-tension constraints, multi-segment pulleys |
| `test_transform` | [`Transform`](@ref) | Spherical coordinate positioning |
| `test_quaternion_conversions` | — | Quaternion ↔ rotation matrix round-trips |
| `test_quaternion_auto_groups` | [`TwistSurface`](@ref) | Auto-generated twist DOFs |
| `test_principal_body_frame` | [`Wing`](@ref AbstractWing) | Principal vs body frame separation |
| `test_heading_calculation` | — | Kite heading from tether geometry |
| `test_section_alignment` | [`Wing`](@ref AbstractWing) | VSM section ↔ structural point mapping |
| `test_profile_law` | — | Atmospheric wind profile verification |
| `test_bench` | — | Performance regression tracking |

### Writing new tests

When adding a new component or equation, follow this pattern:

1. **Build a minimal model** using constructors — only include the
   components needed to test the behavior in question.
2. **Derive the expected result analytically** — free-fall distance,
   terminal velocity, oscillation frequency, etc.
3. **Simulate and compare** — run `next_step!` in a loop and check the
   result against the analytical solution with a tight tolerance.
4. **Keep tests independent** — each test file should build its own
   `SymbolicAWEModel` from scratch. Use `vsm_interval=0` and
   `AERO_NONE` when aerodynamics are not relevant.

---

## Coding style guidelines

Please adhere to the following style guidelines to maintain code quality and
readability:

- **Environment:** Add packages like `Revise` to your global Julia environment,
  not to the project's `Project.toml`.
- **No Magic Numbers:** Avoid hard-coded values (e.g., `9.81`). Define them as
  constants (e.g., `G_EARTH`) or read them from a configuration file.
- **Line Length:** Keep lines under 100 characters, including in documentation.
- **Operators:**
  - Use the tilde `~` for scalar equations in `ModelingToolkit` instead of the
    broadcasted `.~`.
  - Use the `\cdot` operator for the dot product (`⋅`) for improved readability.
  - Enclose binary operators (`+`, `*`, `=`) with single spaces (e.g.,
    `y = a * x + b`).
- **Spacing:** Use a space after a comma (e.g., `my_function(x, y)`).
- **Alignment:** Align assignment operators (`=`) in blocks of related
  assignments to improve readability:
  ```julia
  tether_rhs = [force_eqs[j, i].rhs for j in 1:3]
  kite_rhs   = [force_eqs[j, i+3].rhs for j in 1:3]
  f_xy       = dot(tether_rhs, e_z) * e_z
  ```
- **Settings:** Use the `Settings()` constructor to load the settings for the
  active project. You can specify a file with
  `set = Settings("my_settings.yaml")`. Use `set = Settings("")` to load the
  default settings file.

---

## Known issues and troubleshooting

### Segmentation fault when loading a cached `.bin` model

Cached `.bin` files contain serialized function pointers that are only valid for the
Julia and package version used to create them. Loading a stale `.bin` causes a segfault.

**Solution:** Remove the corrupt `.bin` file in the data directory:

```bash
rm data/2plate_kite/*.bin
```

---

## Source code organization

The source code is organized into modular directories:

- **`src/system_structure/`** — component types and assembly
  - `types.jl`: [`Point`](@ref), [`Segment`](@ref), [`Pulley`](@ref),
    [`Tether`](@ref), [`Winch`](@ref), [`TwistSurface`](@ref), [`Transform`](@ref)
  - `wing.jl`: [`AbstractWing`](@ref)/[`Wing`](@ref) types and the
    [`VSMWing`](@ref)/[`PlateWing`](@ref) constructors, aerodynamic setup
  - `system_structure_core.jl`: [`SystemStructure`](@ref) constructor, reference
    resolution
  - `named_collection.jl`: Symbol-based indexing ([`NamedCollection`](@ref))
  - `transforms.jl`: Spherical coordinate transforms
  - `utilities.jl`: Validation, tether creation, state management
- **`src/generate_system/`** — symbolic equation generation
  - `create_sys.jl`: Top-level orchestrator
  - `point_eqs.jl`, `segment_eqs.jl`, `wing_eqs.jl`, etc.: per-subsystem
    equations
- **`src/yaml_loader.jl`** — YAML configuration file parser
  ([`load_sys_struct_from_yaml`](@ref))
- **`src/linearize.jl`** — VSM linearization ([`linearize!`](@ref))
- **`src/simulate.jl`** — high-level simulation functions ([`sim!`](@ref))
