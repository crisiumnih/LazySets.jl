"""
# Extended help

    ==(X::LazySet, Y::LazySet)

### Algorithm

The default implementation recursively compares the fields of `X` and `Y` until
a mismatch is found.

### Examples

```jldoctest
julia> HalfSpace([1], 1) == HalfSpace([1.0], 1.0)
true

julia> HalfSpace([1], 1) == HalfSpace([1.00000001], 0.99999999)
false

julia> BallInf([0.0], 1.0) == Hyperrectangle([0.0], [1.0])
false
```
"""
function ==(X::LazySet, Y::LazySet)
    # if the common supertype of X and Y is abstract, they cannot be compared
    if isabstracttype(promote_type(typeof(X), typeof(Y)))
        return false
    end

    for f in fieldnames(typeof(X))
        if getfield(X, f) != getfield(Y, f)
            return false
        end
    end

    return true
end
