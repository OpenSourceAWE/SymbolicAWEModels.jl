# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# ==================== WEIGHTED REF POINTS ==================== #

"""Single point from a Symbol ref."""
WeightedRefPoints(ref::Symbol) =
    WeightedRefPoints(NameRef[ref], Int64[], [1.0])

"""Single point from a String ref (converted to Symbol)."""
WeightedRefPoints(ref::AbstractString) =
    WeightedRefPoints(Symbol(ref))

"""
    WeightedRefPoints(id::Integer)

Single point from a resolved index. Stores in `ids`
(not `refs`) so `resolve!` is a no-op.
"""
WeightedRefPoints(id::Integer) =
    WeightedRefPoints(NameRef[], Int64[Int64(id)], [1.0])

"""
    WeightedRefPoints(refs::AbstractVector)

Equal-weight average of multiple ref points, or weighted
if elements are `(name, weight)` tuples.

Supports:
- `[:le, :te]` → equal-weight average
- `[(:le, 0.7), (:te, 0.3)]` → weighted combination
"""
function WeightedRefPoints(refs::AbstractVector)
    isempty(refs) && error(
        "WeightedRefPoints requires at least one " *
        "reference point, got empty vector")
    if refs[1] isa Tuple
        names = NameRef[cell_ref(entry[1]) for entry in refs]
        weights = Float64[Float64(entry[2]) for entry in refs]
        validate_weights!(weights)
        return WeightedRefPoints(names, Int64[], weights)
    end
    names = NameRef[cell_ref(v) for v in refs]
    n = length(names)
    WeightedRefPoints(names, Int64[], fill(1.0 / n, n))
end

# ==================== VSM ENGINE ==================== #

"""
    count_aero_log_points(wings) -> Int

Total extra `SysState` log slots contributed by the wings' aero modes
([`n_aero_log_points`](@ref)): 4 per VSM panel for VSM modes, 4 per section
quad for flat-plate wings, 0 by default.
"""
function count_aero_log_points(wings)
    total = 0
    for wing in wings
        total += n_aero_log_points(wing.aero, wing)
    end
    return total
end

# ==================== CONSTRUCTORS ==================== #

"""Warn and normalize if weights don't sum to 1."""
function validate_weights!(weights::Vector{Float64})
    s = sum(weights)
    s > 0 || error(
        "Ref point weights sum to $s; " *
        "all weights must be positive")
    if !isapprox(s, 1.0; atol=1e-6)
        @warn "Ref point weights sum to $s, " *
            "normalizing to 1.0"
        weights ./= s
    end
    return weights
end

"""
    create_vsm_wing(set::Settings, vsm_set::VortexStepMethod.VSMSettings;
                    prn=true, sort_sections=true)

Create the aero-only `VortexStepMethod.Wing` from the settings: by `ObjWing` from the
`.obj` file `set.model` in the data path if there is one, at panel resolution with its
airfoils from the mesh and its `geometry.yaml` cached under `<model_dir>/obj_geometry`;
otherwise from each wing's `geometry_file`, or the project file's `aero_geometry:`.
"""
function create_vsm_wing(set::Settings, vsm_set::VortexStepMethod.VSMSettings;
                         prn=true, sort_sections=true)
    model_dir = get_data_path()
    obj_path = joinpath(model_dir, set.model)

    if isfile(obj_path)
        prn && @info "Generating wing aero geometry from .obj file"
        n_panels = isempty(vsm_set.wings) ? 56 : vsm_set.wings[1].n_panels
        return VortexStepMethod.ObjWing(obj_path;
            n_panels, crease_frac=set.crease_frac,
            output_dir=joinpath(model_dir, "obj_geometry"), verbose=prn)
    end

    prn && @info "Using provided VSMSettings for wing creation"
    aero_geometry = optional_project_file("aero_geometry")
    for wing_settings in vsm_set.wings
        geometry_file = wing_settings.geometry_file
        if isempty(geometry_file)
            geometry_file = something(aero_geometry, "")
        elseif !isabspath(geometry_file)
            geometry_file = joinpath(model_dir, basename(geometry_file))
        end
        wing_settings.geometry_file = geometry_file
    end
    return VortexStepMethod.Wing(vsm_set; sort_sections)
end

"""
    build_vsm_solver(vsm_aero, vsm_set, dynamics_type, n_stations)

Build the VortexStepMethod `Solver` with `vsm_set`'s solver settings, sized for the
unrefined sections the wing is solved with: one per station for a `PARTICLE_DYNAMICS`
wing with stations, which [`match_aero_sections_to_structure!`](@ref) re-sections onto
them, otherwise those of `vsm_aero`.
"""
function build_vsm_solver(vsm_aero::VortexStepMethod.BodyAerodynamics{P}, vsm_set,
                          dynamics_type, n_stations) where {P}
    n_sections = dynamics_type == PARTICLE_DYNAMICS && n_stations > 0 ?
        n_stations : VortexStepMethod.n_unrefined_sections(vsm_aero)
    return VortexStepMethod.Solver(P, n_sections;
        VortexStepMethod.solver_kwargs(vsm_set.solver_settings)...)
end

"""
    build_vsm_engine(set, vsm_set, dynamics_type; n_stations=0,
                     aero_scale_chord=0.0, aero_z_offset=0.0)

Build a [`VSMEngine`](@ref): the VortexStepMethod wing, aerodynamics and solver (sized by
`n_stations`, [`build_vsm_solver`](@ref)) and the linearization vectors, placeholders for
`RIGID_DYNAMICS` until [`setup_aero!`](@ref) sizes them. `aero_scale_chord` and
`aero_z_offset` adjust VSM forces and panels.
"""
function build_vsm_engine(set::Settings, vsm_set::VortexStepMethod.VSMSettings,
                          dynamics_type::WingType; n_stations=0,
                          aero_scale_chord=0.0, aero_z_offset=0.0,
                          unsteady=UnsteadyAero())
    vsm_wing = create_vsm_wing(set, vsm_set; prn=false, sort_sections=false)
    vsm_aero = VortexStepMethod.BodyAerodynamics([vsm_wing])
    vsm_solver = build_vsm_solver(vsm_aero, vsm_set, dynamics_type, n_stations)

    if dynamics_type == PARTICLE_DYNAMICS
        num_aero_outputs = 0
        num_aero_inputs = 0
    else
        n_stations_est = vsm_wing.n_unrefined_sections
        num_aero_outputs = 6 + n_stations_est
        num_aero_inputs = 5 + n_stations_est
    end

    return VSMEngine(vsm_aero, vsm_wing, vsm_solver,
        zeros(SimFloat, num_aero_inputs),
        zeros(SimFloat, num_aero_outputs),
        zeros(SimFloat, num_aero_outputs, num_aero_inputs),
        nothing, nothing, SimFloat(aero_scale_chord), SimFloat(aero_z_offset), unsteady)
end

# ==================== HELPER FUNCTIONS ==================== #

"""
    frame_sections(vsm_wing)

Every VSM section whose LE/TE are stored in the wing frame — `refined_sections`,
`non_deformed_sections` and `unrefined_sections` — as one iterator so a frame
transform reaches all three. The lists hold independent `Section` objects (VSM
copies on refine), so each is moved exactly once.
"""
frame_sections(vsm_wing) = Iterators.flatten((
    vsm_wing.refined_sections, vsm_wing.non_deformed_sections,
    vsm_wing.unrefined_sections))

"""
    adjust_vsm_panels_to_origin!(vsm_wing, origin_offset)

Subtract `origin_offset` (`[x, y, z]`) from the panel positions of `vsm_wing`, making them
relative to the body-frame origin after it has moved.
"""
function adjust_vsm_panels_to_origin!(vsm_wing, origin_offset)
    for section in frame_sections(vsm_wing)
        section.LE_point .-= origin_offset
        section.TE_point .-= origin_offset
    end
end

"""
    rotate_vsm_sections!(vsm_wing, R)

Rotate all VSM section LE/TE points by rotation matrix `R`.
"""
function rotate_vsm_sections!(vsm_wing, R)
    for section in frame_sections(vsm_wing)
        section.LE_point .= R * section.LE_point
        section.TE_point .= R * section.TE_point
    end
end

"""
    apply_aero_z_offset!(vsm_wing, aero_z_offset)

Shift the panels of `vsm_wing` by `aero_z_offset` [m] along body-frame `+z`, which moves a
`RIGID_DYNAMICS` wing's centre of pressure and so its moment arm. Apply after the COM
adjustment.
"""
function apply_aero_z_offset!(vsm_wing, aero_z_offset)
    abs(aero_z_offset) > 1e-10 || return nothing
    offset_vec = [0.0, 0.0, aero_z_offset]
    for section in frame_sections(vsm_wing)
        section.LE_point .+= offset_vec
        section.TE_point .+= offset_vec
    end
    return nothing
end

# ==================== BODY FRAMES ==================== #

"""The positions of `points`, each point's `field`, indexed as KiteGeometry reads them."""
struct PointPositions{P}
    points::P
    field::Symbol
end
Base.getindex(positions::PointPositions, i) =
    getproperty(positions.points[i], positions.field)

"""`ref_pt` as KiteGeometry's reference points: pairs of index and weight."""
ref_points(ref_pt::WeightedRefPoints) =
    [id => weight for (id, weight) in zip(ref_pt.ids, ref_pt.weights)]

"""
    get_ref_position_from_points(points, ref_pt; field=:pos_w)

Weighted position from structural points, through KiteGeometry's `ref_position`. `field`
selects which point coordinate to read (`:pos_w` or `:pos_ENU`).
"""
get_ref_position_from_points(points::AbstractVector, ref_pt::WeightedRefPoints;
                             field::Symbol=:pos_w) =
    KiteGeometry.ref_position(PointPositions(points, field),
                              zip(ref_pt.ids, ref_pt.weights))

"""
    calc_particle_dynamics_wing_frame(points, z_ref_points, y_ref_points, origin;
                                      field=:pos_w)

The `R_b_to_w` rotation matrix and origin position from structural point positions, each
point's `field`, by KiteGeometry's `wing_frame`.
"""
function calc_particle_dynamics_wing_frame(points::AbstractVector,
                                           z_ref_points::Tuple{WeightedRefPoints,
                                                               WeightedRefPoints},
                                           y_ref_points::Tuple{WeightedRefPoints,
                                                               WeightedRefPoints},
                                           origin::WeightedRefPoints;
                                           field::Symbol=:pos_w)
    frame = KiteGeometry.WingFrame(ref_points(origin), ref_points.(z_ref_points),
                                   ref_points.(y_ref_points))
    return KiteGeometry.wing_frame(PointPositions(points, field), frame)
end

