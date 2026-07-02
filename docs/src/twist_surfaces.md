Yes, there is code in SymbolicAWEModels that assigns wing points to twist surfaces, through several mechanisms:

## 1. YAML / Programmatic Construction (User-Defined)

Points are assigned to twist surfaces explicitly at construction time via the `TwistSurface` constructor. Each `TwistSurface` stores `point_refs` (symbolic names) which are later resolved to numeric `point_idxs`:

```julia
TwistSurface(:ts1, [:p1, :p2], QUASI_STATIC, 0.25)
```

**Resolution** happens in `assign_indices_and_resolve!` in system_structure_core.jl (line 480):

```julia
for twist_surface in twist_surfaces
    twist_surface.point_idxs = Int64[resolve_ref(ref, point_names, "point") 
                                     for ref in twist_surface.point_refs]
end
```

For VSM wings, each twist surface typically has exactly 2 points: LE and TE (a wing section). For flat-plate (`AeroPlate`) wings, each twist surface can have just 1 point (`FIXED` type).

## 2. Auto-Creation for RIGID_DYNAMICS Wings

**`auto_create_twist_surfaces!`** in vsm_refine.jl (line ~547) automatically creates twist surfaces from wing segments (LE/TE pairs) when none are pre-defined:

```julia
for (le_idx, te_idx) in wing_segments
    twist_surface_idx = length(twist_surfaces) + 1
    new_twist_surface = TwistSurface(twist_surface_idx,
        [le_idx, te_idx], DYNAMIC, 0.0)
    new_twist_surface.point_idxs = [le_idx, te_idx]
    push!(twist_surfaces, new_twist_surface)
end
```

## 3. Geometry Check in Symbolic Equations

In point_eqs.jl (line 198), each point checks whether it belongs to a twist surface to determine in-group behavior and moment contributions:

```julia
for twist_surface_ in twist_surfaces
    if point.idx in twist_surface_.point_idxs
        twist_surface = twist_surface_
        found += 1
    end
end
```

## 4. Flat-Plate Aero Mode Lookup

In plate.jl (line 123), each structural wing point finds its twist surface by matching `point_idxs[1]`:

```julia
for gidx in wing.twist_surface_idxs
    if twist_surfaces[gidx].point_idxs[1] == point.idx
        ts_idx = gidx
        break
    end
end
```

## 5. VSM Section-to-TwistSurface Assignment (not per-point)

**`compute_spatial_twist_surface_mapping!`** in system_structure_core.jl (line 611) assigns entire VSM sections to twist surfaces by spatial proximity (nearest center), so a twist surface can own multiple adjacent sections driven by a single twist DOF.

---

**Summary:** The primary mechanism is that each `TwistSurface` object stores `point_idxs` — the point indices belonging to it. This can be set manually (YAML/programmatic), auto-created from wing segments, or looked up via geometry at equation-generation time. The `point_idxs` field is the canonical link between points and twist surfaces.