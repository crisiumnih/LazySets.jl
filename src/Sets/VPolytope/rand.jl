"""
    rand(::Type{VPolytope}; [N]::Type{<:Real}=Float64, [dim]::Int=2,
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing,
         [num_vertices]::Int=-1)

Create a random polytope in vertex representation.

### Input

- `VPolytope`    -- type for dispatch
- `N`            -- (optional, default: `Float64`) numeric type
- `dim`          -- (optional, default: 2) dimension
- `rng`          -- (optional, default: `GLOBAL_RNG`) random number generator
- `seed`         -- (optional, default: `nothing`) seed for reseeding
- `num_vertices` -- (optional, default: `-1`) upper bound on the number of
                    vertices of the polytope (see comment below)

### Output

A random polytope in vertex representation.

### Algorithm

All numbers are normally distributed with mean 0 and standard deviation 1.

The number of vertices can be controlled with the argument `num_vertices`.
For a negative value we choose a random number in the range `dim:5*dim` (except
if `dim == 1`, in which case we choose in the range `1:2`).
Note that we do not guarantee that the vertices are not redundant.
"""
function rand(::Type{VPolytope};
              N::Type{<:Real}=Float64,
              dim::Int=2,
              rng::AbstractRNG=GLOBAL_RNG,
              seed::Union{Int,Nothing}=nothing,
              num_vertices::Int=-1)
    rng = reseed!(rng, seed)
    if num_vertices < 0
        num_vertices = (dim == 1) ? rand(rng, 1:2) : rand(rng, dim:(5 * dim))
    end
    vertices = [randn(rng, N, dim) for i in 1:num_vertices]
    return VPolytope(vertices)
end
