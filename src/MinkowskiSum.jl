"""
    MinkowskiSum <: LazySet

Type that represents the Minkowski sum of two convex sets.

FIELDS:

- s1::LazySet : a convex set
- s2::LazySet : a convex set
"""
struct MinkowskiSum <: LazySet
    s1::LazySet
    s2::LazySet

    MinkowskiSum(s1, s2) = dim(s1) != dim(s2) ? throw(DimensionMismatch) : new(s1, s2)
end

import Base.+

function +(s1::LazySet, s2::LazySet)
    return MinkowskiSum(s1, s2)
end

# ambient dimension of the Minkowski sum of two sets
function dim(ms::MinkowskiSum)::Int64
    return dim(ms.s1)
end

# support vector of the Minkowski sum of two sets
function σ(d::Union{Vector{Float64}, SparseVector{Float64,Int64}}, ms::MinkowskiSum)::Vector{Float64}
    return σ(d, ms.s1) + σ(d, ms.s2)
end

export MinkowskiSum


# ========= Minkowski sum of sets (inspired by CartesianProductArray) ==========
"""
    MinkowskiSumArray <: LazySet

Type that represents the Minkowski sum of a finite number of sets.

FIELDS:

- ``sfarray`` -- array of sets

NOTE:

This type is optimized to be used on the left-hand side of additions only.

"""
mutable struct MinkowskiSumArray <: LazySet
    sfarray::Array{LazySet, 1}

    MinkowskiSumArray(sfarray) = new(sfarray)
end
MinkowskiSumArray() = MinkowskiSumArray(Array{LazySet, 1}(0))

function MinkowskiSumArray(n::Int64)::MinkowskiSumArray
    arr = Array{LazySet, 1}(0)
    sizehint!(arr, n)
    return MinkowskiSumArray(arr)
end

"""
    +(msa, sf)

Adds the support function to the array.
(This function is overridden for more specific types of `sf`).

INPUT:

- ``msa`` -- Minkowski sum array
- ``sf`` -- general support function
"""
function +(msa::MinkowskiSumArray, sf::LazySet)::MinkowskiSumArray
    push!(msa.sfarray, sf)
    return msa
end

"""
    +(msa1, msa2)

Appends the elements of the second array to the first array.

INPUT:

- ``msa1`` -- first Minkowski sum array
- ``msa2`` -- second Minkowski sum array
"""
function +(msa1::MinkowskiSumArray, msa2::MinkowskiSumArray)::MinkowskiSumArray
    append!(msa1.sfarray, msa2.sfarray)
    return msa1
end

"""
    +(msa, vs)

Returns the original array because addition with a void set is a no-op.

INPUT:

- ``msa`` -- Minkowski sum array
- ``vs`` -- void set
"""
function +(msa::MinkowskiSumArray, ::VoidSet)::MinkowskiSumArray
    return msa
end

"""
    dim(ms::MinkowskiSumArray)

Ambient dimension of the Minkowski sum of a finite number of sets.

INPUT:

- ``ms`` -- Minkowski sum array

NOTE:

We do not double-check that the dimensions always match.
"""
function dim(ms::MinkowskiSumArray)::Int64
    return length(ms.sfarray) == 0 ? 0 : dim(ms.sfarray[1])
end

"""
    σ(d::Vector{Float64}, ms::MinkowskiSumArray)

Support vector of the Minkowski sum of a finite number of sets.

INPUT:

- ``d`` -- direction

- ``ms`` -- Minkowski sum array
"""
function σ(d::Union{Vector{Float64}, SparseVector{Float64,Int64}}, ms::MinkowskiSumArray)::Vector{Float64}
    svec = zeros(length(d))
    for sj in ms.sfarray
        svec += σ(d, sj)
    end
    return svec
end

export MinkowskiSumArray
