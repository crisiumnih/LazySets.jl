"""
    HParallelotope{N, VN<:AbstractVector{N}, MN<:AbstractMatrix{N}} <: AbstractZonotope{N}

Type that represents a parallelotope in constraint form.

### Fields

- `directions` -- square matrix where each row is the direction of two parallel
                  constraints
- `offset`     -- vector where each element is the offset of the corresponding
                  constraint

### Notes

Parallelotopes are centrally symmetric convex polytopes in ``ℝ^n``
having ``2n`` pairwise parallel constraints. Every parallelotope is a zonotope.
As such, parallelotopes can be represented in constraint form or in generator
form. The `HParallelotope` type represents parallelotopes in constraint form.

Let ``D ∈ ℝ^{n × n}`` be a matrix and let ``c ∈ ℝ^{2n}`` be
a vector. The parallelotope ``P ⊂ ℝ^n`` generated by the directions
matrix ``D`` and the offset vector ``c`` is given by the set of points
``x ∈ ℝ^n`` such that:

```math
    D_i ⋅ x ≤ c_{i},\\text{  and  } -D_i ⋅ x ≤ c_{n+i}
```
for ``i = 1, …, n``. Here ``D_i`` represents the ``i``-th row of ``D`` and
``c_i`` the ``i``-th component of ``c``.

Note that, although representing a zonotopic set, an `HParallelotope` can be
empty or unbounded if the constraints are unsuitably chosen. This may cause
problems with default methods because the library assumes that zonotopic sets
are non-empty and bounded. Thus such instances are considered illegal. The
default constructor thus checks these conditions, which can be deactivated by
passing the argument `check_consistency=false`.

For details as well as applications of parallelotopes in reachability analysis
we refer to [1] and [2]. For conversions between set representations we refer to
[3].

### References

[1] Tommaso Dreossi, Thao Dang, and Carla Piazza. *Reachability computation for
polynomial dynamical systems.* Formal Methods in System Design 50.1 (2017): 1-38.

[2] Tommaso Dreossi, Thao Dang, and Carla Piazza. *Parallelotope bundles for
polynomial reachability.* Proceedings of the 19th International Conference on
Hybrid Systems: Computation and Control. ACM, 2016.

[3] Matthias Althoff, Olaf Stursberg, and Martin Buss. *Computing reachable sets
of hybrid systems using a combination of zonotopes and polytopes.* Nonlinear
analysis: hybrid systems 4.2 (2010): 233-249.
"""
struct HParallelotope{N,VN<:AbstractVector{N},MN<:AbstractMatrix{N}} <: AbstractZonotope{N}
    directions::MN
    offset::VN

    # default constructor with dimension and consistency check
    function HParallelotope(D::MN, c::VN;
                            check_consistency::Bool=true) where {N,VN<:AbstractVector{N},
                                                                 MN<:AbstractMatrix{N}}
        @assert length(c) == 2 * checksquare(D) "the length of the offset " *
                                                "vector should be twice the size of the directions matrix, " *
                                                "but they have sizes $(length(c)) and $(size(D)) respectively"

        if check_consistency
            require(@__MODULE__, :LazySets; fun_name="HParallelotope")

            P = HPolyhedron(_constraints_list_hparallelotope(D, c, N, VN))
            if isempty(P)
                throw(ArgumentError("the constraints are contradictory"))
            elseif !isbounded(P)
                throw(ArgumentError("the constraints are not bounding"))
            end
        end

        return new{N,VN,MN}(D, c)
    end
end
