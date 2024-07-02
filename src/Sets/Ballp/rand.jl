"""
    rand(::Type{Ballp}; [N]::Type{<:Real}=Float64, [dim]::Int=2,
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing)

Create a random ball in the p-norm.

### Input

- `Ballp` -- type for dispatch
- `N`     -- (optional, default: `Float64`) numeric type
- `dim`   -- (optional, default: 2) dimension
- `rng`   -- (optional, default: `GLOBAL_RNG`) random number generator
- `seed`  -- (optional, default: `nothing`) seed for reseeding

### Output

A random ball in the p-norm.

### Algorithm

The center and radius are normally distributed with mean 0 and standard
deviation 1.
Additionally, the radius is nonnegative.
The p-norm is a normally distributed number ≥ 1 with mean 1 and standard
deviation 1.
"""
function rand(::Type{Ballp};
              N::Type{<:Real}=Float64,
              dim::Int=2,
              rng::AbstractRNG=GLOBAL_RNG,
              seed::Union{Int,Nothing}=nothing)
    rng = reseed!(rng, seed)
    p = one(N) + abs(randn(rng, N))
    center = randn(rng, N, dim)
    radius = abs(randn(rng, N))
    return Ballp(p, center, radius)
end
