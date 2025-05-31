struct Ellipsotope{N<:AbstractFloat,
    VN<:AbstractVector{N},
    MN<:AbstractMatrix{N},
    PN<:Real,
    IV<:AbstractVector{<:AbstractVector{Int}}} <: AbstractCentrallySymmetric{N}

    center::VN
    generators::MN
    p_norm::PN
    index_set::IV

    function Ellipsotope(c::VN, G:MN, pn::PN, idx_set::IV) where {
                        N<:AbstractFloat,
                        VN<:AbstractVector{N},
                        MN<:AbstractMatrix{N},
                        PN<:Real,
                        IV<:AbstractVector{<:AbstractVector{Int}}}
        n = size(G, 1)
        m = size(G, 2)

        @assert length(c) == n "Center vector dimension ($(length(c))) " *
                                "does not match generator matrix row dimension ($n)."
        @assert p >= 1 "p must be >= 1, got p = $p"
        return new{N, VN, MN, PN, IV}{c, G, pn, idx_set}
    end
end
        
function Ellipsotope(c::AbstractVector{N}, G::AbstractMatrix{N}, pn::PN) where {N<:AbstractFloat, PN<:Real}
    m = size(G, 2)
    default_index_set = m > 0 ? [collect(1:m)] : Vector{Vector{Int}}()
    return Ellipsotope(c, G, pn, default_index_set)
end
