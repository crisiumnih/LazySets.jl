"""
    rand(::Type{HPolytope}; [N]::Type{<:Real}=Float64, [dim]::Int=2,
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing)

Create a random polytope in constraint representation.

### Input

- `HPolytope`    -- type for dispatch
- `N`            -- (optional, default: `Float64`) numeric type
- `dim`          -- (optional, default: 2) dimension
- `rng`          -- (optional, default: `GLOBAL_RNG`) random number generator
- `seed`         -- (optional, default: `nothing`) seed for reseeding
- `num_vertices` -- (optional, default: `-1`) upper bound on the number of
                    vertices of the polytope (see comment below)

### Output

A random polytope in constraint representation.

### Algorithm

We create a random polytope in vertex representation and convert it to
constraint representation (hence the argument `num_vertices`).
See [`rand(::Type{VPolytope})`](@ref).
"""
function rand(::Type{HPolytope};
              N::Type{<:Real}=Float64,
              dim::Int=2,
              rng::AbstractRNG=GLOBAL_RNG,
              seed::Union{Int,Nothing}=nothing,
              num_vertices::Int=-1)
    require(@__MODULE__, :LazySets; fun_name="rand")

    if num_vertices == 1
        P = rand(Singleton; N=N, dim=dim, rng=rng, seed=seed)
    elseif dim == 1
        if num_vertices ∉ (-1, 2)
            throw(ArgumentError("creating a 1D random polytope is only supported for 2 vertices"))
        end
        P = rand(Interval; N=N, dim=dim, rng=rng, seed=seed)
    elseif dim == 2
        P = rand(VPolygon; N=N, dim=dim, rng=rng, seed=seed, num_vertices=num_vertices)
    else
        rng = reseed!(rng, seed)
        P = rand(VPolytope; N=N, dim=dim, rng=rng, seed=seed, num_vertices=num_vertices)
    end
    return convert(HPolytope, P)
end
