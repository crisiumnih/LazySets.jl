"""
# Extended help

    translate!(B::Ballp, v::AbstractVector)

### Algorithm

We add the vector to the center of the ball.
"""
function translate!(B::Ballp, v::AbstractVector)
    @assert length(v) == dim(B) "cannot translate a $(dim(B))-dimensional " *
                                "set by a $(length(v))-dimensional vector"
    c = center(B)
    c .+= v
    return B
end
