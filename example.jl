using LazySets
using Plots
using Reexport

using ..LazySets: AbstractCentrallySymmetric
using LinearAlgebra: dot, I, checksquare, isposdef
using Random: AbstractRNG, GLOBAL_RNG
using ReachabilityBase.Arrays: inner
using ReachabilityBase.Distribution: reseed!

struct Ellipsoid2{N<:AbstractFloat,VN<:AbstractVector{N},
    MN<:AbstractMatrix{N}} <: AbstractCentrallySymmetric{N}

    center::VN 
    shape_matrix::MN

    function Ellipsoid2(c::VN, Q::MN; check_posdef::Bool=true) where
        {N<:AbstractFloat,VN<:AbstractVector{N},MN<:AbstractMatrix{N}}
        @assert size(Q, 1) == size(Q, 2) "Shape matrix Q must be square"
        @assert length(c) == size(Q, 1) "Center length must match Q's dimension"

        if check_posdef
            isposdef(Q) || throw(ArgumentError("err"))
        end
        return new{N,VN,MN}(c,Q)
    end
end

function LazySets.center(e::Ellipsoid2)
    return e.center
end

function σ(d::AbstractVector, E::Ellipsoid2)
    if iszero(norm(d, 2))
        return E.center
    end
    Qd = E.shape_matrix * d
    return E.center .+ Qd ./ sqrt(dot(d, Qd))
end


function Ellipsoid2(Q::AbstractMatrix{N}; check_posdef::Bool=true) where {N}
    return Ellipsoid2(zeros(N, size(Q, 1)), Q; check_posdef=check_posdef)
end

c = [1.0, 2.0]          # Center at (1, 2)
Q = [2.0 0.0; 0.0 1.0]  # Shape matrix (axis-aligned ellipsoid)
e = Ellipsoid2(c, Q)
elip = plot(e)
savefig(elip, "ellips.png")