# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# ==================== CONSTRUCTORS ==================== #

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
        nothing, nothing, SimFloat(aero_scale_chord), SimFloat(aero_z_offset), unsteady,
        KVec3[], KVec3[], KVec3[], Vector{Int}[], Vector{Int}[], SimFloat[])
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

"""
    get_ref_position_from_points(points, refs; field=:pos_w)

Weighted position of the reference points `refs`, pairs of point index and weight, from
each point's `field` (`:pos_w`, `:vel_w` or `:pos_ENU`), by KiteGeometry's `ref_position`.
"""
get_ref_position_from_points(points::AbstractVector, refs; field::Symbol=:pos_w) =
    KiteGeometry.ref_position(PointPositions(points, field), refs)

"""
    calc_particle_dynamics_wing_frame(points, frame::KiteGeometry.WingFrame; field=:pos_w)

The `R_b_to_w` rotation matrix and origin position the reference points of `frame` give,
from each point's `field`, by KiteGeometry's `wing_frame`.
"""
calc_particle_dynamics_wing_frame(points::AbstractVector, frame::KiteGeometry.WingFrame;
                                  field::Symbol=:pos_w) =
    KiteGeometry.wing_frame(PointPositions(points, field), frame)

