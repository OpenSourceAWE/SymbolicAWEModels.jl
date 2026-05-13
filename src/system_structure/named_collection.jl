# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

"""
NamedCollection - A wrapper that enables both numeric and symbolic indexing.

This file provides the `NamedCollection` type that wraps a vector of items
and provides lookup by symbolic name via a dictionary mapping.

# Examples
```julia
# Create collection with items that have names
points = [Point(1, ...; name=:kcu), Point(2, ...; name=:le_1)]
nc = NamedCollection(points)

# Access by index (unchanged behavior)
nc[1]  # -> Point 1

# Access by name (new capability)
nc[:kcu]  # -> Point 1
nc[:le_1]  # -> Point 2

# Check if name exists
haskey(nc, :kcu)  # -> true
```
"""

"""
    NamedCollection{T} <: AbstractVector{T}

A wrapper around a vector that enables both numeric and symbolic indexing.

By subtyping `AbstractVector`, this type works transparently with all existing
code that expects vectors, including `@unpack` macros and iteration.

Names are extracted from items that have a `name` field. Items without names
(or with `name=nothing`) are only accessible by numeric index.

$(TYPEDFIELDS)
"""
struct NamedCollection{T} <: AbstractVector{T}
    "The underlying vector of items"
    items::Vector{T}
    "Mapping from symbolic names to indices"
    name_to_idx::Dict{Symbol, Int64}
end

"""
    NamedCollection(items::Vector{T}) where T

Construct a NamedCollection from a vector of items.

Automatically builds the name→index mapping by reading the `name` field
from each item that has one. Items with `name=nothing` are skipped.
"""
function NamedCollection(items::Vector{T}) where T
    name_to_idx = Dict{Symbol, Int64}()

    for (i, item) in enumerate(items)
        # Use try-catch to handle both direct fields and delegated properties (e.g., VSMWing.name -> base.name)
        item_name = try
            item.name
        catch
            nothing
        end
        if !isnothing(item_name)
            name = item_name isa Symbol ? item_name : Symbol(item_name)
            if haskey(name_to_idx, name)
                error("Duplicate name '$name' found at indices $(name_to_idx[name]) and $i")
            end
            name_to_idx[name] = i
        end
    end

    return NamedCollection{T}(items, name_to_idx)
end

# ==================== INDEXING ==================== #

"""Access item by numeric index."""
Base.getindex(nc::NamedCollection, i::Integer) = nc.items[i]

"""Access item by symbolic name."""
function Base.getindex(nc::NamedCollection, name::Symbol)
    if !haskey(nc.name_to_idx, name)
        available = collect(keys(nc.name_to_idx))
        error("Name '$name' not found. Available names: $available")
    end
    return nc.items[nc.name_to_idx[name]]
end

"""Set item by numeric index."""
Base.setindex!(nc::NamedCollection, value, i::Integer) = (nc.items[i] = value)

"""Set item by symbolic name."""
function Base.setindex!(nc::NamedCollection, value, name::Symbol)
    if !haskey(nc.name_to_idx, name)
        error("Name '$name' not found in collection")
    end
    nc.items[nc.name_to_idx[name]] = value
end

# ==================== ITERATION ==================== #

Base.iterate(nc::NamedCollection) = iterate(nc.items)
Base.iterate(nc::NamedCollection, state) = iterate(nc.items, state)

# ==================== LENGTH & SIZING ==================== #

Base.length(nc::NamedCollection) = length(nc.items)
Base.size(nc::NamedCollection) = size(nc.items)
Base.eachindex(nc::NamedCollection) = eachindex(nc.items)
Base.firstindex(nc::NamedCollection) = firstindex(nc.items)
Base.lastindex(nc::NamedCollection) = lastindex(nc.items)
Base.axes(nc::NamedCollection) = axes(nc.items)
Base.isempty(nc::NamedCollection) = isempty(nc.items)

# ==================== DICTIONARY-LIKE INTERFACE ==================== #

"""Check if a symbolic name exists in the collection."""
Base.haskey(nc::NamedCollection, name::Symbol) = haskey(nc.name_to_idx, name)

"""Get all symbolic names in the collection."""
Base.keys(nc::NamedCollection) = keys(nc.name_to_idx)

"""Get all values (same as iterating)."""
Base.values(nc::NamedCollection) = nc.items

# ==================== ARRAY-LIKE INTERFACE ==================== #

Base.push!(nc::NamedCollection, item) = push!(nc.items, item)
Base.eltype(::Type{NamedCollection{T}}) where T = T
Base.eltype(::NamedCollection{T}) where T = T

# Allow findall and other array operations
Base.findall(f::Function, nc::NamedCollection) = findall(f, nc.items)
Base.findall(f::Base.Fix2{typeof(in)}, nc::NamedCollection) =
    findall(f, nc.items)
Base.filter(f::Function, nc::NamedCollection) = filter(f, nc.items)
Base.findfirst(f::Function, nc::NamedCollection) = findfirst(f, nc.items)

# ==================== SHOW ==================== #

function Base.show(io::IO, nc::NamedCollection{T}) where T
    n_named = length(nc.name_to_idx)
    n_total = length(nc.items)
    print(io, "NamedCollection{$T}($n_total items, $n_named named)")
end

function Base.show(io::IO, ::MIME"text/plain", nc::NamedCollection{T}) where T
    n_named = length(nc.name_to_idx)
    n_total = length(nc.items)
    println(io, "NamedCollection{$T} with $n_total items ($n_named named):")
    if n_named > 0
        println(io, "  Names: ", join(sort(collect(keys(nc.name_to_idx))), ", "))
    end
end

# ==================== HELPER FUNCTIONS ==================== #

"""
    get_idx(nc::NamedCollection, name::Symbol)

Get the numeric index for a given symbolic name.
"""
get_idx(nc::NamedCollection, name::Symbol) = nc.name_to_idx[name]

"""
    get_name(nc::NamedCollection, idx::Integer)

Get the symbolic name for a given numeric index, or nothing if unnamed.
"""
function get_name(nc::NamedCollection, idx::Integer)
    item = nc.items[idx]
    # Use try-catch to handle both direct fields and delegated properties
    return try
        item.name
    catch
        nothing
    end
end

"""
    names(nc::NamedCollection)

Return a vector of all symbolic names in order of their indices.
Names are `nothing` for unnamed items.
"""
function names(nc::NamedCollection)
    return [get_name(nc, i) for i in eachindex(nc.items)]
end
